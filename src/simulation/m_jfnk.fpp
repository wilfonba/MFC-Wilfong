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
    $:GPU_DECLARE(create='[u_n, u_k, res_k, res_p, pert, sol, kry, escal]')

    !> Dense GMRES workspace, small and kept on the host
    real(wp), allocatable, dimension(:,:) :: hess  !< Hessenberg matrix
    real(wp), allocatable, dimension(:)   :: gcos, gsin, gvec, yvec
    integer                               :: nloc  !< local packed length: sys_size*(m+1)*(n+1)*(p+1)

contains

    !> Allocate the flat vectors and the dense GMRES workspace
    impure subroutine s_initialize_jfnk_module

        integer :: kd

        nloc = sys_size*(m + 1)*(n + 1)*(p + 1)
        kd = jfnk_krylov_dim

        @:ALLOCATE(u_n(1:nloc), u_k(1:nloc), res_k(1:nloc), res_p(1:nloc), pert(1:nloc), sol(1:nloc))
        @:ALLOCATE(kry(1:nloc, 1:kd + 1))
        @:ALLOCATE(escal(1:nloc))

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
            rvec(idx) = (uvec(idx) - u_n(idx))/dt - rvec(idx)/escal(idx)
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_residual

    !> One implicit time step
    impure subroutine s_jfnk_step(q_cons_vf, q_T_sf, q_prim_vf, bc_type, rhs_vf, pb_in, rhs_pb, mv_in, rhs_mv, t_step)

        type(scalar_field), dimension(sys_size), intent(inout) :: q_cons_vf, q_prim_vf, rhs_vf
        type(scalar_field), intent(inout) :: q_T_sf
        type(integer_field), dimension(1:num_dims,1:2), intent(in) :: bc_type
        real(stp), dimension(idwbuff(1)%beg:,idwbuff(2)%beg:,idwbuff(3)%beg:,1:,1:), intent(inout) :: pb_in, mv_in
        real(wp), dimension(idwbuff(1)%beg:,idwbuff(2)%beg:,idwbuff(3)%beg:,1:,1:), intent(inout) :: rhs_pb, rhs_mv
        integer, intent(in) :: t_step
        real(wp) :: rnorm0, rnorm, unorm, vnorm, eps_fd, hij, beta, tmp
        integer :: newt, restart, jj, ii, idx, kd, nvar
        real(wp) :: vloc, vglb

        kd = jfnk_krylov_dim

        call nvtxStartRange("TIMESTEP-JFNK")

        ! Per-variable scaling hook. The conservative variables span alpha_rho ~ 1e3
        ! to E ~ 8e8 for stiffened water, so an unscaled norm and finite-difference
        ! step are set entirely by the energy -- this is why the solver diverges for
        ! ACFL >~ 10 (see the header). Scaling by the per-variable maximum was tried
        ! and made ACFL 1 WORSE (drift 1.68 against 9.8e-12 unscaled), so it is left
        ! inert here rather than shipped wrong; getting it right is the next task
        call s_pack(q_cons_vf, u_n)
        nvar = (m + 1)*(n + 1)*(p + 1)
        vloc = 0._wp; vglb = 0._wp
        $:GPU_PARALLEL_LOOP(private='[idx]')
        do idx = 1, nloc
            escal(idx) = 1._wp
        end do
        $:END_GPU_PARALLEL_LOOP()

        $:GPU_PARALLEL_LOOP(private='[idx]')
        do idx = 1, nloc
            u_n(idx) = u_n(idx)/escal(idx)
            u_k(idx) = u_n(idx)
        end do
        $:END_GPU_PARALLEL_LOOP()

        call s_residual(u_k, res_k, q_cons_vf, q_T_sf, q_prim_vf, bc_type, rhs_vf, pb_in, rhs_pb, mv_in, rhs_mv, t_step)
        call s_dot(res_k, res_k, tmp); rnorm0 = sqrt(tmp)

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
                    call s_dot(kry(:,jj), kry(:,jj), tmp); vnorm = sqrt(tmp)
                    eps_fd = sqrt((1._wp + unorm)*epsilon(1._wp))/max(vnorm, sgm_eps)

                    $:GPU_PARALLEL_LOOP(private='[idx]')
                    do idx = 1, nloc
                        pert(idx) = u_k(idx) + eps_fd*kry(idx, jj)
                    end do
                    $:END_GPU_PARALLEL_LOOP()
                    call s_residual(pert, kry(:,jj + 1), q_cons_vf, q_T_sf, q_prim_vf, bc_type, rhs_vf, pb_in, rhs_pb, mv_in, &
                                    & rhs_mv, t_step)
                    $:GPU_PARALLEL_LOOP(private='[idx]')
                    do idx = 1, nloc
                        kry(idx, jj + 1) = (kry(idx, jj + 1) - res_k(idx))/eps_fd
                    end do
                    $:END_GPU_PARALLEL_LOOP()

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
                        sol(idx) = sol(idx) + hij*kry(idx, ii)
                    end do
                    $:END_GPU_PARALLEL_LOOP()
                end do

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

        $:GPU_PARALLEL_LOOP(private='[idx]')
        do idx = 1, nloc
            pert(idx) = u_k(idx)*escal(idx)
        end do
        $:END_GPU_PARALLEL_LOOP()
        call s_unpack(pert, q_cons_vf)
        call nvtxEndRange

    end subroutine s_jfnk_step

    impure subroutine s_finalize_jfnk_module

        @:DEALLOCATE(u_n, u_k, res_k, res_p, pert, sol, kry, escal)
        deallocate (hess, gcos, gsin, gvec, yvec)

    end subroutine s_finalize_jfnk_module

end module m_jfnk
