!>
!! @file m_jfnk.fpp
!! @brief Jacobian-free Newton-Krylov implicit time stepping

#:include 'case.fpp'
#:include 'macros.fpp'

!> Implicit time integration by Jacobian-free Newton-Krylov (Knoll & Keyes, JCP 193, 2004). Backward Euler is solved for the whole
!! discrete system at once,
!!
!!     R(u) = (u - u^n)/dt - RHS(u) = 0,
!!
!! by Newton's method, whose linear systems are solved with GMRES using only directional derivatives,
!!
!!     J(u) v ~ (R(u + eps*v) - R(u))/eps,
!!
!! so no Jacobian is ever formed. The point of doing it this way here is that the
!! acoustic stiffness is removed WITHOUT splitting the system: the semi-implicit
!! pressure projection had to compose a divergence, a gradient and a Laplacian
!! consistently, and on MFC's collocated grid it does not (see
!! examples/1D_contact_semiimplicit/README.md) -- a 1e-7 divergence error becomes
!! hundreds of Pa through a rho*c^2 gain of ~2.6e9. Here the split reappears only as
!! an optional preconditioner, where being approximately right costs Krylov
!! iterations instead of corrupting the answer.
!> @brief Jacobian-free Newton-Krylov implicit time stepping
module m_jfnk

    use m_derived_types
    use m_global_parameters
    use m_mpi_proxy
    use m_rhs, only: s_compute_rhs
    use m_boundary_common
    use m_helper
    use m_nvtx

    implicit none

    private; public :: s_initialize_jfnk_module, s_jfnk_step, s_finalize_jfnk_module

    !> Flat working vectors. The state is packed into a contiguous array so the Krylov algebra is plain vector arithmetic; only the
    !! residual evaluation needs the scalar_field layout back
    real(wp), allocatable, dimension(:)   :: u_n    !< state at t^n
    real(wp), allocatable, dimension(:)   :: u_k    !< current Newton iterate
    real(wp), allocatable, dimension(:)   :: res_k  !< R(u_k)
    real(wp), allocatable, dimension(:)   :: res_p  !< scratch inside the residual
    real(wp), allocatable, dimension(:)   :: pert   !< perturbed iterate u_k + eps*v
    real(wp), allocatable, dimension(:)   :: sol    !< Newton update
    real(wp), allocatable, dimension(:,:) :: kry    !< Krylov basis, (nloc, dim+1)
    real(wp), allocatable, dimension(:)   :: escal  !< per-element variable scale
    real(wp), allocatable, dimension(:)   :: uacc   !< accumulated known state for the current stage
    real(wp), allocatable, dimension(:,:) :: kstg   !< stage RHS values k_j = RHS(u_j)
    real(wp), allocatable, dimension(:)   :: pcd    !< preconditioner diagonal M^-1
    real(wp), allocatable, dimension(:)   :: zvec   !< M^-1 applied to a Krylov vector
    $:GPU_DECLARE(create='[u_n, u_k, res_k, res_p, pert, sol, kry, escal, uacc, kstg, pcd, zvec]')

    !> Dense GMRES workspace, small and kept on the host
    real(wp), allocatable, dimension(:,:) :: hess  !< Hessenberg matrix
    real(wp), allocatable, dimension(:)   :: gcos, gsin, gvec, yvec
    integer                               :: nloc  !< local packed length: sys_size*(m+1)*(n+1)*(p+1)

    !> ESDIRK tableau. Both higher-order options are one-step methods with an explicit first stage and a single repeated diagonal
    !! entry, so no solution history is stored and one preconditioner serves every stage
    real(wp), dimension(4, 4) :: esd_a     !< Butcher A
    real(wp), dimension(4)    :: esd_c     !< Butcher c (stage times)
    real(wp)                  :: esd_aii   !< repeated diagonal, used by the residual
    integer                   :: esd_ns    !< number of stages
    logical                   :: esd_expl  !< whether stage 1 is explicit

contains

    !> Allocate the flat vectors and the dense GMRES workspace
    impure subroutine s_initialize_jfnk_module

        integer  :: kd
        real(wp) :: gam

        nloc = sys_size*(m + 1)*(n + 1)*(p + 1)
        kd = jfnk_krylov_dim

        esd_a = 0._wp
        if (jfnk_order == 2) then
            ! ESDIRK2, three stages, L-stable and stiffly accurate. Also known as
            ! TR-BDF2 and as Kennedy & Carpenter ARK2(2)3L[2]SA -- despite the
            ! former name it is a one-step Runge-Kutta method, not a multistep one
            gam = 1._wp - 0.5_wp*sqrt(2._wp)
            esd_ns = 3; esd_expl = .true.; esd_aii = gam
            esd_c(1) = 0._wp; esd_c(2) = 2._wp*gam; esd_c(3) = 1._wp
            esd_a(2, 1) = gam; esd_a(2, 2) = gam
            esd_a(3, 1) = 0.25_wp*sqrt(2._wp); esd_a(3, 2) = 0.25_wp*sqrt(2._wp); esd_a(3, 3) = gam
        else if (jfnk_order == 3) then
            ! ESDIRK3, Kennedy & Carpenter ARK3(2)4L[2]SA, four stages, L-stable and
            ! stiffly accurate; gam is the root of x^3 - 3x^2 + 3x/2 - 1/6
            gam = 0.435866521508459_wp
            esd_ns = 4; esd_expl = .true.; esd_aii = gam
            esd_c(1) = 0._wp; esd_c(2) = 2._wp*gam; esd_c(3) = 0.6_wp; esd_c(4) = 1._wp
            esd_a(2, 1) = gam; esd_a(2, 2) = gam
            esd_a(3, 1) = 0.2576482460664272_wp; esd_a(3, 2) = -0.09351476757488625_wp; esd_a(3, 3) = gam
            esd_a(4, 1) = 0.1876410243467238_wp; esd_a(4, 2) = -0.5952974735769549_wp
            esd_a(4, 3) = 0.9717899277217721_wp; esd_a(4, 4) = gam
        else
            ! backward Euler as a one-stage, fully implicit tableau
            esd_ns = 1; esd_expl = .false.; esd_aii = 1._wp
            esd_c(1) = 1._wp
            esd_a(1, 1) = 1._wp
        end if

        @:ALLOCATE(u_n(1:nloc), u_k(1:nloc), res_k(1:nloc), res_p(1:nloc), pert(1:nloc), sol(1:nloc))
        @:ALLOCATE(kry(1:nloc, 1:kd + 1))
        @:ALLOCATE(escal(1:nloc))
        @:ALLOCATE(uacc(1:nloc))
        @:ALLOCATE(pcd(1:nloc))
        @:ALLOCATE(zvec(1:nloc))
        @:ALLOCATE(kstg(1:nloc, 1:4))

        allocate (hess(1:kd + 1,1:kd), gcos(1:kd + 1), gsin(1:kd + 1), gvec(1:kd + 1), yvec(1:kd))

    end subroutine s_initialize_jfnk_module

    !> Pack a conservative field array into the flat vector `vec`
    subroutine s_pack(q_vf, vec)

        type(scalar_field), dimension(sys_size), intent(in) :: q_vf
        real(wp), dimension(1:nloc), intent(out)            :: vec
        integer                                             :: i, j, k, l, idx

        $:GPU_PARALLEL_LOOP(collapse=4, private='[i, j, k, l, idx]')
        do i = 1, sys_size
            do l = 0, p
                do k = 0, n
                    do j = 0, m
                        idx = (((i - 1)*(p + 1) + l)*(n + 1) + k)*(m + 1) + j + 1
                        vec(idx) = real(q_vf(i)%sf(j, k, l), wp)
                    end do
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_pack

    !> Unpack the flat vector `vec` back into a conservative field array
    subroutine s_unpack(vec, q_vf)

        real(wp), dimension(1:nloc), intent(in)                :: vec
        type(scalar_field), dimension(sys_size), intent(inout) :: q_vf
        integer                                                :: i, j, k, l, idx

        $:GPU_PARALLEL_LOOP(collapse=4, private='[i, j, k, l, idx]')
        do i = 1, sys_size
            do l = 0, p
                do k = 0, n
                    do j = 0, m
                        idx = (((i - 1)*(p + 1) + l)*(n + 1) + k)*(m + 1) + j + 1
                        q_vf(i)%sf(j, k, l) = real(vec(idx), stp)
                    end do
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_unpack

    !> Global inner product of two flat vectors
    impure subroutine s_dot(a, b, res)

        real(wp), dimension(1:nloc), intent(in) :: a, b
        real(wp), intent(out)                   :: res
        real(wp)                                :: loc
        integer                                 :: idx

        loc = 0._wp
        $:GPU_PARALLEL_LOOP(private='[idx]', reduction='[[loc]]', reductionOp='[+]')
        do idx = 1, nloc
            loc = loc + a(idx)*b(idx)
        end do
        $:END_GPU_PARALLEL_LOOP()

        call s_mpi_allreduce_sum(loc, res)

    end subroutine s_dot

    !> Backward-Euler residual R(u) = (u - u^n)/dt - RHS(u), evaluated by unpacking the trial state, calling the existing RHS, and
    !! packing the result. This is the only place the solver touches the discretization, which is why nothing here has to know about
    !! fluxes, reconstruction or the equation set
    impure subroutine s_residual(uvec, rvec, q_cons_vf, q_T_sf, q_prim_vf, bc_type, rhs_vf, pb_in, rhs_pb, mv_in, rhs_mv, t_step)

        real(wp), dimension(1:nloc), intent(in)                                                    :: uvec
        real(wp), dimension(1:nloc), intent(out)                                                   :: rvec
        type(scalar_field), dimension(sys_size), intent(inout)                                     :: q_cons_vf, q_prim_vf, rhs_vf
        type(scalar_field), intent(inout)                                                          :: q_T_sf
        type(integer_field), dimension(1:num_dims,1:2), intent(in)                                 :: bc_type
        real(stp), dimension(idwbuff(1)%beg:,idwbuff(2)%beg:,idwbuff(3)%beg:,1:,1:), intent(inout) :: pb_in, mv_in
        real(wp), dimension(idwbuff(1)%beg:,idwbuff(2)%beg:,idwbuff(3)%beg:,1:,1:), intent(inout)  :: rhs_pb, rhs_mv
        integer, intent(in)                                                                        :: t_step
        integer                                                                                    :: idx

        ! uvec is in SCALED units: the conservative variables span alpha_rho ~ 1e3 to
        ! E ~ 8e8 for stiffened water, so an unscaled norm and finite-difference step
        ! are set entirely by the energy and the other equations are invisible to the
        ! Krylov space. Working with u/D makes every equation contribute comparably

        $:GPU_PARALLEL_LOOP(private='[idx]')
        do idx = 1, nloc
            res_p(idx) = uvec(idx)*escal(idx)
        end do
        $:END_GPU_PARALLEL_LOOP()
        call s_unpack(res_p, q_cons_vf)
        call s_compute_rhs(q_cons_vf, q_T_sf, q_prim_vf, bc_type, rhs_vf, pb_in, rhs_pb, mv_in, rhs_mv, t_step, 1)
        call s_pack(rhs_vf, rvec)

        $:GPU_PARALLEL_LOOP(private='[idx]')
        do idx = 1, nloc
            rvec(idx) = (uvec(idx) - uacc(idx))/(dt*esd_aii) - rvec(idx)/escal(idx)
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_residual

    !> Build the preconditioner M^-1, applied on the right inside GMRES.
    !!
    !! Currently the identity. A wave-speed diagonal, M^-1 = 1/(1/(dt*a_ii) +
    !! (|u| + c)/dx), was implemented and MEASURED WORSE: at ACFL 2 it took the
    !! Newton residual from 1e-5 down to only 1.09 while tripling the Krylov count
    !! (84-210 -> 558-570), and it turned the ACFL 5 case from converging into NaN.
    !! That is consistent with the theory -- a diagonal scaling does not touch the
    !! ACFL-dependent conditioning, which comes from the elliptic acoustic coupling
    !! in the off-diagonal blocks, and reshaping the spectrum without addressing it
    !! can make GMRES worse.
    !!
    !! The physics-based alternative -- one pressure-Helmholtz solve per application,
    !! the operator that inverts the acoustic coupling -- was then built and MEASURED
    !! WORSE STILL, in a way worth recording because the reasoning that motivates it
    !! is wrong. Writing z = tau*v and correcting it the way a projection does gives
    !!     pp - rho*c^2*tau^2*div(rho^-1 grad pp) = z_E/gamma - rho*c^2*tau*div(z_mom/rho)
    !! with z_mom -= tau*grad(pp) and z_E = gamma*pp. Damping that correction by a
    !! factor w and sweeping it gives sin(J M^-1 v, v) errors of 0.003, 11.3, 33.9,
    !! 67.8, 113.0 at w = 0, 0.1, 0.3, 0.6, 1 -- exactly linear, with no minimum. The
    !! correction is therefore not a mistuned approximation that better coefficients
    !! would rescue: J returns none of it. Adding the scheme's own upwind dissipation
    !! diagonal, which at an acoustic CFL of 2 is twice as large as 1/tau and so the
    !! biggest term in the momentum row, moves the error by under 10% across a
    !! coefficient sweep of 0 to 4, so that is not the missing piece either.
    !!
    !! The reason is structural, and it is the same one that defeated the projection
    !! method as a scheme (examples/1D_contact_semiimplicit/README.md). M^-1 makes the
    !! momentum block 265 times larger than its input -- correct linear algebra, since
    !! a pressure perturbation drives momentum -- so J has to annihilate that back down
    !! to the input to 0.4% for the preconditioner to help. Its gradient and M's
    !! gradient are different discrete operators, so it cannot. Any M built from the
    !! continuous acoustic equations fails here for the same reason at any CFL where
    !! the amplification is large. A useful M must be assembled from a discretization
    !! that shares J's stencils -- a first-order version of the same flux -- not
    !! derived from the PDEs; that is what the two-phase JFNK literature does, and it
    !! is why those codes precondition on staggered grids where the triple is exact.
    subroutine s_build_precond(q_prim_vf)

        type(scalar_field), dimension(sys_size), intent(in) :: q_prim_vf
        integer                                             :: idx

        $:GPU_PARALLEL_LOOP(private='[idx]')
        do idx = 1, nloc
            pcd(idx) = 1._wp
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_build_precond

    !> One implicit time step
    impure subroutine s_jfnk_step(q_cons_vf, q_T_sf, q_prim_vf, bc_type, rhs_vf, pb_in, rhs_pb, mv_in, rhs_mv, t_step)

        type(scalar_field), dimension(sys_size), intent(inout) :: q_cons_vf, q_prim_vf, rhs_vf
        type(scalar_field), intent(inout) :: q_T_sf
        type(integer_field), dimension(1:num_dims,1:2), intent(in) :: bc_type
        real(stp), dimension(idwbuff(1)%beg:,idwbuff(2)%beg:,idwbuff(3)%beg:,1:,1:), intent(inout) :: pb_in, mv_in
        real(wp), dimension(idwbuff(1)%beg:,idwbuff(2)%beg:,idwbuff(3)%beg:,1:,1:), intent(inout) :: rhs_pb, rhs_mv
        integer, intent(in) :: t_step
        real(wp) :: rnorm0, rnorm, unorm, vnorm, eps_fd, hij, beta, tmp, pcq
        integer :: newt, restart, jj, ii, idx, kd, nvar, stg
        integer :: nkry
        real(wp) :: acc, t_base
        real(wp) :: vloc, vglb

        kd = jfnk_krylov_dim

        call nvtxStartRange("TIMESTEP-JFNK")

        ! Per-variable scale D_i = max|u_i| over the domain. Without it the residual
        ! is dominated entirely by the energy equation -- E ~ 8e8 for stiffened water
        ! against alpha_rho ~ 1e3 -- so Newton and GMRES both converge on energy alone
        ! and the rest of the system is invisible to them. Working with u/D and R/D
        ! makes every equation contribute comparably
        call s_pack(q_cons_vf, u_n)
        nvar = (m + 1)*(n + 1)*(p + 1)
        do ii = 1, sys_size
            vloc = 0._wp
            $:GPU_PARALLEL_LOOP(private='[idx]', reduction='[[vloc]]', reductionOp='[max]')
            do idx = (ii - 1)*nvar + 1, ii*nvar
                vloc = max(vloc, abs(u_n(idx)))
            end do
            $:END_GPU_PARALLEL_LOOP()
            call s_mpi_allreduce_max(vloc, vglb)
            vglb = max(vglb, sgm_eps)
            $:GPU_PARALLEL_LOOP(private='[idx]')
            do idx = (ii - 1)*nvar + 1, ii*nvar
                escal(idx) = vglb
            end do
            $:END_GPU_PARALLEL_LOOP()
        end do

        $:GPU_PARALLEL_LOOP(private='[idx]')
        do idx = 1, nloc
            u_n(idx) = u_n(idx)/escal(idx)
            u_k(idx) = u_n(idx)
        end do
        $:END_GPU_PARALLEL_LOOP()

        ! Explicit first stage: k_1 = RHS(u^n). Evaluated through the residual with
        ! uacc = u^n, which returns -RHS since the difference term vanishes
        if (esd_expl) then
            mytime = t_base + esd_c(1)*dt
            $:GPU_UPDATE(device='[mytime]')
            $:GPU_PARALLEL_LOOP(private='[idx]')
            do idx = 1, nloc
                uacc(idx) = u_n(idx)
            end do
            $:END_GPU_PARALLEL_LOOP()
            call s_residual(u_n, res_k, q_cons_vf, q_T_sf, q_prim_vf, bc_type, rhs_vf, pb_in, rhs_pb, mv_in, rhs_mv, t_step)
            $:GPU_PARALLEL_LOOP(private='[idx]')
            do idx = 1, nloc
                kstg(idx, 1) = -res_k(idx)
            end do
            $:END_GPU_PARALLEL_LOOP()
        end if

        do stg = merge(2, 1, esd_expl), esd_ns
            ! uacc = u^n + dt*sum_{j<stg} a(stg,j)*k_j, the known part of this stage
            $:GPU_PARALLEL_LOOP(private='[idx]')
            do idx = 1, nloc
                uacc(idx) = u_n(idx)
            end do
            $:END_GPU_PARALLEL_LOOP()
            do ii = 1, stg - 1
                acc = dt*esd_a(stg, ii)
                $:GPU_PARALLEL_LOOP(private='[idx]')
                do idx = 1, nloc
                    uacc(idx) = uacc(idx) + acc*kstg(idx, ii)
                end do
                $:END_GPU_PARALLEL_LOOP()
            end do

            call s_residual(u_k, res_k, q_cons_vf, q_T_sf, q_prim_vf, bc_type, rhs_vf, pb_in, rhs_pb, mv_in, rhs_mv, t_step)
            call s_build_precond(q_prim_vf)
            call s_dot(res_k, res_k, tmp); rnorm0 = sqrt(tmp)

            nkry = 0; pcq = 0._wp
            do newt = 1, jfnk_max_newton
                call s_dot(res_k, res_k, tmp); rnorm = sqrt(tmp)
                if (rnorm <= jfnk_newton_tol*max(rnorm0, sgm_eps)) exit

                call s_dot(u_k, u_k, tmp); unorm = sqrt(tmp)

                ! GMRES on J*sol = -res_k, restarted
                $:GPU_PARALLEL_LOOP(private='[idx]')
                do idx = 1, nloc
                    sol(idx) = 0._wp
                end do
                $:END_GPU_PARALLEL_LOOP()

                do restart = 1, jfnk_max_restarts
                    ! initial Krylov direction is the current linear residual, which for a zero start is simply -R
                    $:GPU_PARALLEL_LOOP(private='[idx]')
                    do idx = 1, nloc
                        kry(idx, 1) = -res_k(idx)
                    end do
                    $:END_GPU_PARALLEL_LOOP()
                    call s_dot(kry(:,1), kry(:,1), tmp); beta = sqrt(tmp)
                    if (beta <= sgm_eps) exit
                    $:GPU_PARALLEL_LOOP(private='[idx]')
                    do idx = 1, nloc
                        kry(idx, 1) = kry(idx, 1)/beta
                    end do
                    $:END_GPU_PARALLEL_LOOP()

                    hess = 0._wp; gvec = 0._wp; gvec(1) = beta

                    do jj = 1, kd
                        ! Matrix-free directional derivative. The step is scaled by the
                        ! state norm so it stays well above round-off but small enough that
                        ! the difference quotient is still a derivative
                        ! Right preconditioning: the Arnoldi direction is M^-1 v, so
                        ! GMRES builds a Krylov space for J*M^-1 and the update is
                        ! recovered with one more application of M^-1. The outer Newton
                        ! iteration still sees the true discrete system, so M only has to
                        ! be approximately right -- a poor M costs iterations, not accuracy
                        $:GPU_PARALLEL_LOOP(private='[idx]')
                        do idx = 1, nloc
                            zvec(idx) = pcd(idx)*kry(idx, jj)
                        end do
                        $:END_GPU_PARALLEL_LOOP()
                        call s_dot(zvec, zvec, tmp); vnorm = sqrt(tmp)
                        eps_fd = sqrt((1._wp + unorm)*epsilon(1._wp))/max(vnorm, sgm_eps)

                        $:GPU_PARALLEL_LOOP(private='[idx]')
                        do idx = 1, nloc
                            pert(idx) = u_k(idx) + eps_fd*zvec(idx)
                        end do
                        $:END_GPU_PARALLEL_LOOP()
                        call s_residual(pert, kry(:,jj + 1), q_cons_vf, q_T_sf, q_prim_vf, bc_type, rhs_vf, pb_in, rhs_pb, mv_in, &
                                        & rhs_mv, t_step)
                        $:GPU_PARALLEL_LOOP(private='[idx]')
                        do idx = 1, nloc
                            kry(idx, jj + 1) = (kry(idx, jj + 1) - res_k(idx))/eps_fd
                        end do
                        $:END_GPU_PARALLEL_LOOP()

                        ! Preconditioner health, as the sine of the angle between
                        ! J*M^-1*v and v on the residual direction. A useful M^-1 leaves
                        ! the two nearly parallel; a value approaching one means M^-1
                        ! sends the residual somewhere J does not send it back from, so
                        ! GMRES is being actively hindered. Being an angle it is scale
                        ! free, so an M with the wrong units is not flattered, and it
                        ! costs two dot products in the step where a bad M would
                        ! otherwise only show up as a Krylov-count campaign
                        if (jj == 1 .and. run_time_info) then
                            call s_dot(kry(:,2), kry(:,2), tmp)
                            call s_dot(kry(:,2), kry(:,1), hij)
                            pcq = max(pcq, sqrt(max(1._wp - hij*hij/max(tmp, sgm_eps), 0._wp)))
                        end if

                        ! modified Gram-Schmidt
                        do ii = 1, jj
                            call s_dot(kry(:,jj + 1), kry(:,ii), hij)
                            hess(ii, jj) = hij
                            $:GPU_PARALLEL_LOOP(private='[idx]')
                            do idx = 1, nloc
                                kry(idx, jj + 1) = kry(idx, jj + 1) - hij*kry(idx, ii)
                            end do
                            $:END_GPU_PARALLEL_LOOP()
                        end do
                        call s_dot(kry(:,jj + 1), kry(:,jj + 1), tmp)
                        hess(jj + 1, jj) = sqrt(tmp)
                        if (hess(jj + 1, jj) > sgm_eps) then
                            $:GPU_PARALLEL_LOOP(private='[idx]')
                            do idx = 1, nloc
                                kry(idx, jj + 1) = kry(idx, jj + 1)/hess(jj + 1, jj)
                            end do
                            $:END_GPU_PARALLEL_LOOP()
                        end if

                        ! apply previous Givens rotations, then eliminate the subdiagonal
                        do ii = 1, jj - 1
                            tmp = gcos(ii)*hess(ii, jj) + gsin(ii)*hess(ii + 1, jj)
                            hess(ii + 1, jj) = -gsin(ii)*hess(ii, jj) + gcos(ii)*hess(ii + 1, jj)
                            hess(ii, jj) = tmp
                        end do
                        tmp = sqrt(hess(jj, jj)**2 + hess(jj + 1, jj)**2)
                        if (tmp <= sgm_eps) tmp = sgm_eps
                        gcos(jj) = hess(jj, jj)/tmp
                        gsin(jj) = hess(jj + 1, jj)/tmp
                        hess(jj, jj) = tmp
                        hess(jj + 1, jj) = 0._wp
                        gvec(jj + 1) = -gsin(jj)*gvec(jj)
                        gvec(jj) = gcos(jj)*gvec(jj)

                        nkry = nkry + 1
                        if (abs(gvec(jj + 1)) <= jfnk_krylov_tol*beta) exit
                    end do
                    jj = min(jj, kd)

                    ! back-substitute and form the update
                    do ii = jj, 1, -1
                        yvec(ii) = gvec(ii)
                        do idx = ii + 1, jj
                            yvec(ii) = yvec(ii) - hess(ii, idx)*yvec(idx)
                        end do
                        yvec(ii) = yvec(ii)/sign(max(abs(hess(ii, ii)), sgm_eps), hess(ii, ii))
                    end do
                    do ii = 1, jj
                        hij = yvec(ii)
                        $:GPU_PARALLEL_LOOP(private='[idx]')
                        do idx = 1, nloc
                            sol(idx) = sol(idx) + hij*pcd(idx)*kry(idx, ii)
                        end do
                        $:END_GPU_PARALLEL_LOOP()
                    end do

                    nkry = nkry + 1
                    if (abs(gvec(jj + 1)) <= jfnk_krylov_tol*beta) exit
                end do

                ! Newton update, then a fresh residual
                $:GPU_PARALLEL_LOOP(private='[idx]')
                do idx = 1, nloc
                    u_k(idx) = u_k(idx) + sol(idx)
                end do
                $:END_GPU_PARALLEL_LOOP()
                call s_residual(u_k, res_k, q_cons_vf, q_T_sf, q_prim_vf, bc_type, rhs_vf, pb_in, rhs_pb, mv_in, rhs_mv, t_step)
            end do

            if (proc_rank == 0 .and. run_time_info) then
                print '(A, ES10.3)', '   jfnk precond sin(J M^-1 v, v) ', pcq
                print '(A, I2, A, I3, A, I4, A, ES10.3, A, ES10.3)', '   jfnk stg ', stg, ' newton ', newt - 1, ' krylov ', nkry, &
                    & ' |R0| ', rnorm0, ' -> ', rnorm
            end if

            ! The converged stage satisfies (u_stg - uacc)/(dt*a_ii) = RHS(u_stg), so
            ! k_stg comes straight from the stage relation and costs no extra
            ! residual evaluation
            acc = 1._wp/(dt*esd_aii)
            $:GPU_PARALLEL_LOOP(private='[idx]')
            do idx = 1, nloc
                kstg(idx, stg) = (u_k(idx) - uacc(idx))*acc
            end do
            $:END_GPU_PARALLEL_LOOP()
        end do

        ! Both tableaux are stiffly accurate (b equals the last row of A), so the
        ! final stage IS the new state and no separate update is needed
        $:GPU_PARALLEL_LOOP(private='[idx]')
        do idx = 1, nloc
            pert(idx) = u_k(idx)*escal(idx)
        end do
        $:END_GPU_PARALLEL_LOOP()
        call s_unpack(pert, q_cons_vf)
        mytime = t_base + dt
        $:GPU_UPDATE(device='[mytime]')
        call nvtxEndRange

    end subroutine s_jfnk_step

    impure subroutine s_finalize_jfnk_module

        @:DEALLOCATE(u_n, u_k, res_k, res_p, pert, sol, kry, escal, uacc, kstg, pcd, zvec)
        deallocate (hess, gcos, gsin, gvec, yvec)

    end subroutine s_finalize_jfnk_module

end module m_jfnk
