!>
!! @file m_rhs_staggered.fpp
!! @brief Right-hand side assembly for the staggered discretization

#:include 'case.fpp'
#:include 'macros.fpp'

!> Right-hand side terms for the staggered scheme, built up one phase at a time. At present it carries the conservative transport of
!! a cell scalar by a prescribed face velocity -- the piece every conserved quantity in the five-equation model needs -- and the
!! test that verifies it.
!!
!! The flux is conservative by construction: cell j gains what face j-1 delivers and
!! loses what face j takes, with each face evaluated once, so the sum over cells
!! telescopes to the boundary regardless of what the reconstruction returns. Errors in
!! the reconstruction can therefore make the answer inaccurate but never unconservative,
!! which is why the test below reports accuracy and conservation separately.
!> @brief Right-hand side assembly for the staggered discretization
module m_rhs_staggered

    use m_derived_types
    use m_global_parameters
    use m_mpi_proxy
    use m_staggered
    use m_weno_staggered

    implicit none

    private; public :: s_initialize_rhs_staggered_module, s_scalar_transport_rhs, s_staggered_scalar_test, s_staggered_from_cons, &
        & s_staggered_to_cons, s_staggered_time_step, s_staggered_freestream_test, s_finalize_rhs_staggered_module

    !> Reconstruction and flux buffers, held for the life of the run. They must be device resident and must not be allocated inside
    !! the right-hand side, which is called several times per stage
    real(wp), allocatable, dimension(:,:,:) :: qfl, qfr, flx
    $:GPU_DECLARE(create='[qfl, qfr, flx]')

    !> The staggered solver owns its state for the life of a step, in working precision and on its own three-ghost extent. q_cons_vf
    !! is read at the start and written back for output, which keeps the storage-precision fields and the idwbuff extents out of
    !! every kernel below
    real(wp), allocatable, dimension(:,:,:,:) :: qs        !< cell state, one slot per equation
    real(wp), allocatable, dimension(:,:,:)   :: um        !< FACE momentum, the primary velocity unknown
    real(wp), allocatable, dimension(:,:,:)   :: pcl, rcl  !< cell pressure and mixture density
    real(wp), allocatable, dimension(:,:,:,:) :: kq, qsv   !< RK stage slope and stage state, cells
    real(wp), allocatable, dimension(:,:,:)   :: ku, umv   !< RK stage slope and stage state, faces
    real(wp), allocatable, dimension(:,:,:)   :: tmpa, tmpb
    real(wp), allocatable, dimension(:,:,:,:) :: fst       !< upwind face state of each cell variable
    real(wp), allocatable, dimension(:,:,:)   :: mfl       !< total mass flux per face, shared with the momentum control volume
    $:GPU_DECLARE(create='[qs, um, pcl, rcl, kq, qsv, ku, umv, tmpa, tmpb, mfl, fst]')

contains

    impure subroutine s_initialize_rhs_staggered_module

        @:ALLOCATE(qfl(-3:m + 3, -3:n + 3, -3:p + 3))
        @:ALLOCATE(qfr(-3:m + 3, -3:n + 3, -3:p + 3))
        @:ALLOCATE(flx(-3:m + 3, -3:n + 3, -3:p + 3))
        @:ALLOCATE(qs(-3:m + 3, -3:n + 3, -3:p + 3, 1:sys_size))
        @:ALLOCATE(qsv(-3:m + 3, -3:n + 3, -3:p + 3, 1:sys_size))
        @:ALLOCATE(kq(-3:m + 3, -3:n + 3, -3:p + 3, 1:sys_size))
        @:ALLOCATE(um(-3:m + 3, -3:n + 3, -3:p + 3))
        @:ALLOCATE(umv(-3:m + 3, -3:n + 3, -3:p + 3))
        @:ALLOCATE(ku(-3:m + 3, -3:n + 3, -3:p + 3))
        @:ALLOCATE(pcl(-3:m + 3, -3:n + 3, -3:p + 3))
        @:ALLOCATE(rcl(-3:m + 3, -3:n + 3, -3:p + 3))
        @:ALLOCATE(tmpa(-3:m + 3, -3:n + 3, -3:p + 3))
        @:ALLOCATE(tmpb(-3:m + 3, -3:n + 3, -3:p + 3))
        @:ALLOCATE(mfl(-3:m + 3, -3:n + 3, -3:p + 3))
        @:ALLOCATE(fst(-3:m + 3, -3:n + 3, -3:p + 3, 1:sys_size))

    end subroutine s_initialize_rhs_staggered_module

    !> Conservative transport of a cell scalar by a prescribed face velocity, rhs = -div(q*u). The face value is taken from
    !! whichever side the face velocity comes from, which is the whole of the upwinding for a passively advected scalar
    subroutine s_scalar_transport_rhs(qc, ufin, rhs)

        real(wp), dimension(-3:m + 3,-3:n + 3,-3:p + 3), intent(in)            :: qc
        real(wp), dimension(-3:m + 3,-3:n + 3,-3:p + 3,1:num_dims), intent(in) :: ufin
        real(wp), dimension(-3:m + 3,-3:n + 3,-3:p + 3), intent(out)           :: rhs
        integer                                                                :: j, k, l, d

        $:GPU_PARALLEL_LOOP(collapse=3, private='[j, k, l]')
        do l = 0, p
            do k = 0, n
                do j = 0, m
                    rhs(j, k, l) = 0._wp
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

        do d = 1, num_dims
            call s_weno_face_states(qc, qfl, qfr, d)

            ! one flux per face, evaluated once
            $:GPU_PARALLEL_LOOP(collapse=3, private='[j, k, l]')
            do l = -1, p + 1
                do k = -1, n + 1
                    do j = -1, m
                        if (ufin(j, k, l, d) >= 0._wp) then
                            flx(j, k, l) = ufin(j, k, l, d)*qfl(j, k, l)
                        else
                            flx(j, k, l) = ufin(j, k, l, d)*qfr(j, k, l)
                        end if
                    end do
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()

            #:for DIM, DXV, IDXV, IM1 in [(1, 'dx', 'j', 'j - 1, k, l'), &
                (2, 'dy', 'k', 'j, k - 1, l'), (3, 'dz', 'l', 'j, k, l - 1')]
                if (d == ${DIM}$) then
                    $:GPU_PARALLEL_LOOP(collapse=3, private='[j, k, l]')
                    do l = 0, p
                        do k = 0, n
                            do j = 0, m
                                rhs(j, k, l) = rhs(j, k, l) - (flx(j, k, l) - flx(${IM1}$))/${DXV}$(${IDXV}$)
                            end do
                        end do
                    end do
                    $:END_GPU_PARALLEL_LOOP()
                end if
            #:endfor
        end do

    end subroutine s_scalar_transport_rhs

    !> Acceptance gate for the staggered scalar transport.
    !!
    !! A smooth profile is advected one full period on a uniform, divergence-free face
    !! velocity, so the exact answer is the initial condition and any difference is error.
    !! Three things are reported:
    !!
    !!   1. Free-stream preservation. A constant field must stay exactly constant. This
    !!      is the cheapest check that the reconstruction weights sum to one and that the
    !!      flux differencing lines up with the face indexing; almost every stencil slip
    !!      breaks it, and it is independent of resolution.
    !!   2. Conservation. The summed quantity must not move at all. It is separate from
    !!      accuracy because the flux form guarantees it even when the reconstruction is
    !!      wrong, so agreement here is not evidence the scheme is right -- but a failure
    !!      means the flux is being evaluated inconsistently on the two sides of a face.
    !!   3. Accuracy after a full period, which at a small enough CFL is the spatial
    !!      order. Compare across resolutions to read it off.
    impure subroutine s_staggered_scalar_test

        real(wp), allocatable, dimension(:,:,:)   :: qc, q0, k1, qt
        real(wp), allocatable, dimension(:,:,:,:) :: vfp
        real(wp)                                  :: qx, uadv, dtl, tend, tnow, err, mass0, mass1, fs, vc, lx
        integer                                   :: j, k, l, d, nst

        allocate (qc(-3:m + 3,-3:n + 3,-3:p + 3), q0(-3:m + 3,-3:n + 3,-3:p + 3))
        allocate (k1(-3:m + 3,-3:n + 3,-3:p + 3), qt(-3:m + 3,-3:n + 3,-3:p + 3))
        allocate (vfp(-3:m + 3,-3:n + 3,-3:p + 3,1:num_dims))

        qx = 2._wp*pi/real(m + 1, wp)
        uadv = 1._wp
        lx = 0._wp
        do j = 0, m
            lx = lx + dx(j)
        end do

        do l = -3, p + 3
            do k = -3, n + 3
                do j = -3, m + 3
                    q0(j, k, l) = 1._wp + 0.5_wp*sin(qx*j)
                    do d = 1, num_dims
                        vfp(j, k, l, d) = uadv
                    end do
                end do
            end do
        end do

        ! 1. free stream: a constant field on this velocity must not move at all
        qc = 1._wp
        call s_scalar_transport_rhs(qc, vfp, k1)
        fs = 0._wp
        do l = 0, p
            do k = 0, n
                do j = 0, m
                    fs = max(fs, abs(k1(j, k, l)))
                end do
            end do
        end do

        ! one period at a CFL small enough that the third-order time error stays under the fifth-order spatial one
        qc = q0
        dtl = 0.1_wp*dx(0)/uadv
        tend = lx/uadv
        nst = max(int(tend/dtl) + 1, 1)
        dtl = tend/real(nst, wp)

        mass0 = 0._wp
        do l = 0, p
            do k = 0, n
                do j = 0, m
                    vc = dx(j)
                    mass0 = mass0 + vc*qc(j, k, l)
                end do
            end do
        end do

        ! SSP-RK3, with the periodic wrap applied straight to the ghost layers. The real
        ! boundary machinery arrives with the full right-hand side; here the point is the
        ! interior operator, and a wrap makes the exact solution exactly known
        do j = 1, nst
            call s_ssp_rk3_step(qc, vfp, k1, qt, dtl)
        end do

        mass1 = 0._wp; err = 0._wp
        do l = 0, p
            do k = 0, n
                do j = 0, m
                    vc = dx(j)
                    mass1 = mass1 + vc*qc(j, k, l)
                    err = max(err, abs(qc(j, k, l) - q0(j, k, l)))
                end do
            end do
        end do

        if (proc_rank == 0) then
            print '(A)', ' Staggered scalar transport test'
            print '(A, I6, A, I6, A)', '   cells ', m + 1, '   steps ', nst, '   one full period'
            print '(A, ES12.5)', '   free stream  max|rhs| on a constant field   = ', fs
            print '(A, ES12.5)', '   conservation |sum q dV - initial| / initial = ', abs(mass1 - mass0)/abs(mass0)
            print '(A, ES12.5)', '   accuracy     max|q - exact| after a period  = ', err
        end if

        deallocate (qc, q0, k1, qt, vfp)

    end subroutine s_staggered_scalar_test

    !> One SSP-RK3 step of the scalar transport, wrapping the ghost layers periodically
    subroutine s_ssp_rk3_step(qc, vfp, k1, qt, dtl)

        real(wp), dimension(-3:m + 3,-3:n + 3,-3:p + 3), intent(inout)         :: qc, k1, qt
        real(wp), dimension(-3:m + 3,-3:n + 3,-3:p + 3,1:num_dims), intent(in) :: vfp
        real(wp), intent(in)                                                   :: dtl
        integer                                                                :: j, k, l

        call s_wrap(qc)
        call s_scalar_transport_rhs(qc, vfp, k1)
        do l = 0, p
            do k = 0, n
                do j = 0, m
                    qt(j, k, l) = qc(j, k, l) + dtl*k1(j, k, l)
                end do
            end do
        end do

        call s_wrap(qt)
        call s_scalar_transport_rhs(qt, vfp, k1)
        do l = 0, p
            do k = 0, n
                do j = 0, m
                    qt(j, k, l) = 0.75_wp*qc(j, k, l) + 0.25_wp*(qt(j, k, l) + dtl*k1(j, k, l))
                end do
            end do
        end do

        call s_wrap(qt)
        call s_scalar_transport_rhs(qt, vfp, k1)
        do l = 0, p
            do k = 0, n
                do j = 0, m
                    qc(j, k, l) = (qc(j, k, l) + 2._wp*(qt(j, k, l) + dtl*k1(j, k, l)))/3._wp
                end do
            end do
        end do

    end subroutine s_ssp_rk3_step

    !> Periodic wrap of the three ghost layers, for the test only
    subroutine s_wrap(qc)

        real(wp), dimension(-3:m + 3,-3:n + 3,-3:p + 3), intent(inout) :: qc
        integer                                                        :: j, k, l, g

        do g = 1, 3
            do l = 0, p
                do k = 0, n
                    qc(-g, k, l) = qc(m + 1 - g, k, l)
                    qc(m + g, k, l) = qc(g - 1, k, l)
                end do
            end do
        end do

    end subroutine s_wrap

    !> Load the staggered state from the conservative fields. The face momentum is the average of the two cells it lies between,
    !! which is the transfer that keeps the momentum sum telescoping
    subroutine s_staggered_from_cons(q_cons_vf)

        type(scalar_field), dimension(sys_size), intent(in) :: q_cons_vf
        integer                                             :: i, j, k, l

        do i = 1, sys_size
            $:GPU_PARALLEL_LOOP(collapse=3, private='[j, k, l]')
            do l = 0, p
                do k = 0, n
                    do j = 0, m
                        qs(j, k, l, i) = real(q_cons_vf(i)%sf(j, k, l), wp)
                    end do
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()
        end do
        call s_staggered_bc_cells

        $:GPU_PARALLEL_LOOP(collapse=3, private='[j, k, l]')
        do l = -3, p + 2
            do k = -3, n + 2
                do j = -3, m + 2
                    um(j, k, l) = 0.5_wp*(qs(j, k, l, eqn_idx%mom%beg) + qs(j + 1, k, l, eqn_idx%mom%beg))
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()
        call s_staggered_bc_faces

    end subroutine s_staggered_from_cons

    !> Write the staggered state back, collapsing the face momentum onto cell centres so that output, the CFL check and every
    !! collocated consumer see the usual fields
    subroutine s_staggered_to_cons(q_cons_vf)

        type(scalar_field), dimension(sys_size), intent(inout) :: q_cons_vf
        integer                                                :: i, j, k, l

        do i = 1, sys_size
            $:GPU_PARALLEL_LOOP(collapse=3, private='[j, k, l]')
            do l = 0, p
                do k = 0, n
                    do j = 0, m
                        if (i == eqn_idx%mom%beg) then
                            q_cons_vf(i)%sf(j, k, l) = real(0.5_wp*(um(j - 1, k, l) + um(j, k, l)), stp)
                        else
                            q_cons_vf(i)%sf(j, k, l) = real(qs(j, k, l, i), stp)
                        end if
                    end do
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()
        end do

    end subroutine s_staggered_to_cons

    !> Ghost layers for the cell state: periodic wrap, otherwise zero gradient
    subroutine s_staggered_bc_cells

        integer :: i, j, k, l, g

        do i = 1, sys_size
            do g = 1, 3
                $:GPU_PARALLEL_LOOP(collapse=2, private='[k, l]')
                do l = 0, p
                    do k = 0, n
                        if (bc_x%beg == BC_PERIODIC) then
                            qs(-g, k, l, i) = qs(m + 1 - g, k, l, i)
                            qs(m + g, k, l, i) = qs(g - 1, k, l, i)
                        else
                            qs(-g, k, l, i) = qs(0, k, l, i)
                            qs(m + g, k, l, i) = qs(m, k, l, i)
                        end if
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()
            end do
        end do

    end subroutine s_staggered_bc_cells

    !> Ghost layers for the face momentum. Unlike a cell field a face at a wall is a primary unknown rather than something to
    !! extrapolate into, so a solid boundary sets it to zero rather than copying the interior
    subroutine s_staggered_bc_faces

        call s_staggered_bc_face(um)

    end subroutine s_staggered_bc_faces

    !> Ghost fill for any face-indexed quantity. Face -1 and face m are the same face on a periodic grid, which is what sets the
    !! offsets here
    subroutine s_staggered_bc_face(qf)

        real(wp), dimension(-3:m + 3,-3:n + 3,-3:p + 3), intent(inout) :: qf
        integer                                                        :: k, l, g

        do g = 1, 3
            $:GPU_PARALLEL_LOOP(collapse=2, private='[k, l]')
            do l = 0, p
                do k = 0, n
                    if (bc_x%beg == BC_PERIODIC) then
                        qf(-g, k, l) = qf(m + 1 - g, k, l)
                        qf(m + g, k, l) = qf(g - 1, k, l)
                    else
                        qf(-g, k, l) = qf(0, k, l)
                        qf(m + g, k, l) = qf(m, k, l)
                    end if
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()
        end do

    end subroutine s_staggered_bc_face

    !> Mixture density and pressure at cell centres. The kinetic energy is the piece that has to be defined compatibly with the
    !! staggered momentum: the two faces bounding a cell each carry half of it, so ke = rho/4 * (u_lo^2 + u_hi^2). Getting this
    !! wrong leaves the energy slowly drifting rather than failing outright
    subroutine s_staggered_primitives

        integer  :: i, j, k, l
        real(wp) :: gm, pf, qv, ke, ulo, uhi, rlo, rhi

        ! rho first, on its own. Writing it and reading its neighbour in one loop is a
        ! race, and the clamps that made that look safe were also hiding a read past the
        ! end of the array

        $:GPU_PARALLEL_LOOP(collapse=3, private='[i, j, k, l]')
        do l = -3, p + 3
            do k = -3, n + 3
                do j = -3, m + 3
                    rcl(j, k, l) = 0._wp
                    $:GPU_LOOP(parallelism='[seq]')
                    do i = 1, num_fluids
                        rcl(j, k, l) = rcl(j, k, l) + qs(j, k, l, i)
                    end do
                    rcl(j, k, l) = max(rcl(j, k, l), sgm_eps)
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

        $:GPU_PARALLEL_LOOP(collapse=3, private='[i, j, k, l, gm, pf, qv, ke, ulo, uhi, rlo, rhi]')
        do l = -2, p + 2
            do k = -2, n + 2
                do j = -2, m + 2
                    gm = 0._wp; pf = 0._wp; qv = 0._wp
                    $:GPU_LOOP(parallelism='[seq]')
                    do i = 1, num_fluids
                        gm = gm + qs(j, k, l, eqn_idx%adv%beg + i - 1)*gammas(i)
                        pf = pf + qs(j, k, l, eqn_idx%adv%beg + i - 1)*pi_infs(i)
                        qv = qv + qs(j, k, l, i)*qvs(i)
                    end do
                    rlo = 0.5_wp*(rcl(j, k, l) + rcl(j - 1, k, l))
                    rhi = 0.5_wp*(rcl(j, k, l) + rcl(j + 1, k, l))
                    ulo = um(j - 1, k, l)/max(rlo, sgm_eps)
                    uhi = um(j, k, l)/max(rhi, sgm_eps)
                    ke = 0.25_wp*rcl(j, k, l)*(ulo*ulo + uhi*uhi)
                    pcl(j, k, l) = (qs(j, k, l, eqn_idx%E) - ke - pf - qv)/max(gm, sgm_eps)
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

        ! pcl is one layer short of what the reconstruction reaches, because forming it
        ! needs rho at its own neighbours. Fill the last layer the way a cell field is
        ! filled rather than leaving the stencil to read uninitialised memory
        call s_staggered_bc_cell(pcl)

    end subroutine s_staggered_primitives

    !> Ghost fill for a single cell-indexed quantity
    subroutine s_staggered_bc_cell(qc)

        real(wp), dimension(-3:m + 3,-3:n + 3,-3:p + 3), intent(inout) :: qc
        integer                                                        :: k, l, g

        do g = 1, 3
            $:GPU_PARALLEL_LOOP(collapse=2, private='[k, l]')
            do l = 0, p
                do k = 0, n
                    if (bc_x%beg == BC_PERIODIC) then
                        qc(-g, k, l) = qc(m + 1 - g, k, l)
                        qc(m + g, k, l) = qc(g - 1, k, l)
                    else
                        qc(-g, k, l) = qc(0, k, l)
                        qc(m + g, k, l) = qc(m, k, l)
                    end if
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()
        end do

    end subroutine s_staggered_bc_cell

    !> Full right-hand side of the five-equation model on the staggered grid, in one dimension. Cell equations are fluxed by the
    !! face velocity; the face momentum is fluxed through the cell centres its own control volume is bounded by, and feels the
    !! pressure gradient over the same centre-to-centre distance the divergence is adjoint to -- which is the whole reason for
    !! staggering
    subroutine s_staggered_rhs

        integer  :: i, j, k, l
        real(wp) :: uc, dvg

        call s_staggered_primitives

        ! face velocity from the face momentum and the interpolated density
        $:GPU_PARALLEL_LOOP(collapse=3, private='[j, k, l]')
        do l = -3, p + 2
            do k = -3, n + 2
                do j = -3, m + 2
                    tmpa(j, k, l) = um(j, k, l)/max(0.5_wp*(rcl(j, k, l) + rcl(j + 1, k, l)), sgm_eps)
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

        ! div(u), needed by the volume-fraction equation
        $:GPU_PARALLEL_LOOP(collapse=3, private='[j, k, l]')
        do l = 0, p
            do k = 0, n
                do j = 0, m
                    tmpb(j, k, l) = (tmpa(j, k, l) - tmpa(j - 1, k, l))/dx(j)
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

        ! partial densities, accumulating the total mass flux each face carries
        $:GPU_PARALLEL_LOOP(collapse=3, private='[j, k, l]')
        do l = -3, p + 3
            do k = -3, n + 3
                do j = -3, m + 3
                    mfl(j, k, l) = 0._wp
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

        do i = eqn_idx%cont%beg, eqn_idx%cont%end
            call s_transport_1d(qs(:,:,:,i), tmpa, kq(:,:,:,i), fst(:,:,:,i))
            $:GPU_PARALLEL_LOOP(collapse=3, private='[j, k, l]')
            do l = 0, p
                do k = 0, n
                    do j = -1, m
                        mfl(j, k, l) = mfl(j, k, l) + flx(j, k, l)
                    end do
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()
        end do

        ! The momentum control volume reads the mass flux one face beyond where the
        ! reconstruction can produce it, so give it the same boundary treatment the face
        ! velocity gets
        call s_staggered_bc_face(mfl)

        ! volume fractions
        do i = eqn_idx%adv%beg, eqn_idx%adv%end
            call s_transport_1d(qs(:,:,:,i), tmpa, kq(:,:,:,i), fst(:,:,:,i))
        end do
        do i = eqn_idx%adv%beg, eqn_idx%adv%end
            $:GPU_PARALLEL_LOOP(collapse=3, private='[j, k, l]')
            do l = 0, p
                do k = 0, n
                    do j = 0, m
                        kq(j, k, l, i) = kq(j, k, l, i) + qs(j, k, l, i)*tmpb(j, k, l)
                    end do
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()
        end do

        ! Energy. The flux is (E + p)*u, but E must NOT be reconstructed as a
        ! conservative variable: WENO's weights depend on the field, so the weights it
        ! picks for E differ from those it picked for alpha, and the face energy then
        ! disagrees with the face volume fraction about where the interface is. Uniform
        ! pressure stops being uniform -- 31 Pa in a single step here, out of 1e5. Build
        ! the face energy instead from the face states the other equations already used,
        ! plus a reconstruction of the PRESSURE, which is exact when p is uniform because
        ! any reconstruction of a constant is that constant
        call s_weno_face_states(pcl, qfl, qfr, 1)
        $:GPU_PARALLEL_LOOP(collapse=3, private='[i, j, k, l, uc, dvg]')
        do l = 0, p
            do k = 0, n
                do j = -1, m
                    if (tmpa(j, k, l) >= 0._wp) then
                        uc = qfl(j, k, l)
                    else
                        uc = qfr(j, k, l)
                    end if
                    dvg = 0._wp
                    flx(j, k, l) = uc
                    $:GPU_LOOP(parallelism='[seq]')
                    do i = 1, num_fluids
                        flx(j, k, l) = flx(j, k, l) + fst(j, k, l, eqn_idx%adv%beg + i - 1)*(gammas(i)*uc + pi_infs(i)) + fst(j, &
                            & k, l, i)*qvs(i)
                        dvg = dvg + fst(j, k, l, i)
                    end do
                    flx(j, k, l) = tmpa(j, k, l)*(flx(j, k, l) + 0.5_wp*dvg*tmpa(j, k, l)*tmpa(j, k, l))
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

        $:GPU_PARALLEL_LOOP(collapse=3, private='[j, k, l]')
        do l = 0, p
            do k = 0, n
                do j = 0, m
                    kq(j, k, l, eqn_idx%E) = -(flx(j, k, l) - flx(j - 1, k, l))/dx(j)
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

        ! Face momentum. Its control volume is bounded by the cell centres, and its mass
        ! flux there MUST be the one the continuity equation used, averaged onto that
        ! boundary. Reconstructing rho independently instead leaves the momentum and the
        ! mass advecting on two different discrete operators, which at a density jump
        ! manufactures velocity out of nothing even though every equation is separately
        ! conservative and free-stream preserving. Reconstructing the face array puts the
        ! velocity states exactly at the cell centres, so index j-1 serves cell j
        call s_weno_face_states(tmpa, qfl, qfr, 1)
        $:GPU_PARALLEL_LOOP(collapse=3, private='[j, k, l, uc]')
        do l = 0, p
            do k = 0, n
                do j = 0, m + 1
                    uc = 0.5_wp*(mfl(j - 1, k, l) + mfl(j, k, l))
                    if (uc >= 0._wp) then
                        flx(j, k, l) = uc*qfl(j - 1, k, l)
                    else
                        flx(j, k, l) = uc*qfr(j - 1, k, l)
                    end if
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

        $:GPU_PARALLEL_LOOP(collapse=3, private='[j, k, l, dvg]')
        do l = 0, p
            do k = 0, n
                do j = 0, m
                    dvg = 0.5_wp*(dx(j) + dx(j + 1))
                    ku(j, k, l) = -(flx(j + 1, k, l) - flx(j, k, l))/dvg - (pcl(j + 1, k, l) - pcl(j, k, l))/dvg
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_staggered_rhs

    !> One-dimensional conservative transport of a cell field by a face velocity
    subroutine s_transport_1d(qc, ufc, rhs, fup)

        real(wp), dimension(-3:m + 3,-3:n + 3,-3:p + 3), intent(in)  :: qc, ufc
        real(wp), dimension(-3:m + 3,-3:n + 3,-3:p + 3), intent(out) :: rhs, fup
        integer                                                      :: j, k, l

        call s_weno_face_states(qc, qfl, qfr, 1)

        $:GPU_PARALLEL_LOOP(collapse=3, private='[j, k, l]')
        do l = 0, p
            do k = 0, n
                do j = -1, m
                    if (ufc(j, k, l) >= 0._wp) then
                        fup(j, k, l) = qfl(j, k, l)
                    else
                        fup(j, k, l) = qfr(j, k, l)
                    end if
                    flx(j, k, l) = ufc(j, k, l)*fup(j, k, l)
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

        $:GPU_PARALLEL_LOOP(collapse=3, private='[j, k, l]')
        do l = 0, p
            do k = 0, n
                do j = 0, m
                    rhs(j, k, l) = -(flx(j, k, l) - flx(j - 1, k, l))/dx(j)
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_transport_1d

    !> One SSP-RK3 step of the staggered scheme
    subroutine s_staggered_time_step(q_cons_vf, dtl)

        type(scalar_field), dimension(sys_size), intent(inout) :: q_cons_vf
        real(wp), intent(in)                                   :: dtl
        integer                                                :: i, j, k, l, st
        real(wp)                                               :: c1, c2

        call s_staggered_from_cons(q_cons_vf)

        do i = 1, sys_size
            $:GPU_PARALLEL_LOOP(collapse=3, private='[j, k, l]')
            do l = -3, p + 3
                do k = -3, n + 3
                    do j = -3, m + 3
                        qsv(j, k, l, i) = qs(j, k, l, i)
                    end do
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()
        end do
        $:GPU_PARALLEL_LOOP(collapse=3, private='[j, k, l]')
        do l = -3, p + 3
            do k = -3, n + 3
                do j = -3, m + 3
                    umv(j, k, l) = um(j, k, l)
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

        do st = 1, 3
            call s_staggered_rhs

            if (st == 1) then
                c1 = 0._wp; c2 = 1._wp
            else if (st == 2) then
                c1 = 0.75_wp; c2 = 0.25_wp
            else
                c1 = 1._wp/3._wp; c2 = 2._wp/3._wp
            end if

            do i = 1, sys_size
                if (i == eqn_idx%mom%beg) cycle
                $:GPU_PARALLEL_LOOP(collapse=3, private='[j, k, l]')
                do l = 0, p
                    do k = 0, n
                        do j = 0, m
                            qs(j, k, l, i) = c1*qsv(j, k, l, i) + c2*(qs(j, k, l, i) + dtl*kq(j, k, l, i))
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()
            end do
            $:GPU_PARALLEL_LOOP(collapse=3, private='[j, k, l]')
            do l = 0, p
                do k = 0, n
                    do j = 0, m
                        um(j, k, l) = c1*umv(j, k, l) + c2*(um(j, k, l) + dtl*ku(j, k, l))
                    end do
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()

            call s_staggered_bc_cells
            call s_staggered_bc_faces
        end do

        call s_staggered_to_cons(q_cons_vf)

    end subroutine s_staggered_time_step

    !> Free-stream test for the full system. A uniform state moving at a uniform speed is an exact steady solution of every
    !! equation, so each right-hand side must be exactly zero. Reported per equation, because which one is non-zero says which term
    !! is wrong
    impure subroutine s_staggered_freestream_test

        integer  :: i, j, k, l
        real(wp) :: gm, pf, qv, r0, u0, p0, mx, al

        ! Two phases either side of a sharp interface, at uniform pressure and uniform
        ! velocity. That is an exact solution translating at u0, so every right-hand side
        ! is again exactly zero -- but unlike the single-fluid case it is only zero if the
        ! energy flux and the volume-fraction equation are consistent with each other
        ! across the jump. This is the condition conservative multicomponent schemes are
        ! known to violate, and it is what a Riemann solver is doing for the collocated
        ! path, so a failure here is expected until that treatment is added

        u0 = 5._wp; p0 = 1.e5_wp

        ! Interior only, then let the boundary conditions fill the ghosts. Writing the
        ! profile straight into the ghost layers puts water where periodicity requires
        ! air, and the state is inconsistent at the seam before the solver has run
        do l = 0, p
            do k = 0, n
                do j = 0, m
                    if (j < (m + 1)/2) then
                        al = 1._wp
                    else
                        al = 0._wp
                    end if
                    qs(j, k, l, eqn_idx%adv%beg) = al
                    if (num_fluids > 1) qs(j, k, l, eqn_idx%adv%beg + 1) = 1._wp - al
                    qs(j, k, l, 1) = 1000._wp*al
                    if (num_fluids > 1) qs(j, k, l, 2) = 1.2_wp*(1._wp - al)
                    r0 = qs(j, k, l, 1)
                    gm = gammas(1)*al; pf = pi_infs(1)*al; qv = qvs(1)*qs(j, k, l, 1)
                    if (num_fluids > 1) then
                        r0 = r0 + qs(j, k, l, 2)
                        gm = gm + gammas(2)*(1._wp - al)
                        pf = pf + pi_infs(2)*(1._wp - al)
                        qv = qv + qvs(2)*qs(j, k, l, 2)
                    end if
                    qs(j, k, l, eqn_idx%E) = gm*p0 + pf + qv + 0.5_wp*r0*u0*u0
                end do
            end do
        end do
        call s_staggered_bc_cells

        ! face momentum from the face density, so the face velocity is exactly u0
        do l = 0, p
            do k = 0, n
                do j = -1, m
                    r0 = 0._wp
                    do i = 1, num_fluids
                        r0 = r0 + 0.5_wp*(qs(j, k, l, i) + qs(j + 1, k, l, i))
                    end do
                    um(j, k, l) = r0*u0
                end do
            end do
        end do
        call s_staggered_bc_face(um)

        ! A translating interface has a large right-hand side by construction, so the
        ! condition is not that the rates vanish. It is Abgrall's: a state at uniform
        ! pressure and uniform velocity must stay that way. One explicit step, then look
        ! at the spread of u and p
        call s_staggered_primitives
        mx = 0._wp; al = 0._wp
        do l = 0, p
            do k = 0, n
                do j = 1, m - 1
                    mx = max(mx, abs(pcl(j, k, l) - p0))
                    al = max(al, abs(um(j, k, l)/max(0.5_wp*(rcl(j, k, l) + rcl(j + 1, k, l)), sgm_eps) - u0))
                end do
            end do
        end do
        if (proc_rank == 0) print '(A, ES12.5, A, ES12.5)', '   initial state: max|p - p0| ', mx, '   max|u - u0| ', al

        call s_staggered_rhs

        do i = 1, sys_size
            do l = 0, p
                do k = 0, n
                    do j = 0, m
                        qs(j, k, l, i) = qs(j, k, l, i) + dt*kq(j, k, l, i)
                    end do
                end do
            end do
        end do
        do l = 0, p
            do k = 0, n
                do j = 0, m
                    um(j, k, l) = um(j, k, l) + dt*ku(j, k, l)
                end do
            end do
        end do
        call s_staggered_bc_cells
        call s_staggered_bc_face(um)
        call s_staggered_primitives

        mx = 0._wp; al = 0._wp
        do l = 0, p
            do k = 0, n
                do j = 1, m - 1
                    mx = max(mx, abs(pcl(j, k, l) - p0))
                    al = max(al, abs(um(j, k, l)/max(0.5_wp*(rcl(j, k, l) + rcl(j + 1, k, l)), sgm_eps) - u0))
                end do
            end do
        end do
        if (proc_rank == 0) then
            print '(A, ES12.5, A, ES12.5)', '   after one step: max|p - p0| ', mx, '   max|u - u0| ', al
        end if

    end subroutine s_staggered_freestream_test

    impure subroutine s_finalize_rhs_staggered_module

        @:DEALLOCATE(qfl, qfr, flx, qs, qsv, kq, um, umv, ku, pcl, rcl, tmpa, tmpb, mfl, fst)

    end subroutine s_finalize_rhs_staggered_module

end module m_rhs_staggered
