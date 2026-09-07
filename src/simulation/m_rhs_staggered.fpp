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

    private; public :: s_initialize_rhs_staggered_module, s_scalar_transport_rhs, s_staggered_scalar_test, &
        & s_finalize_rhs_staggered_module

    !> Reconstruction and flux buffers, held for the life of the run. They must be device resident and must not be allocated inside
    !! the right-hand side, which is called several times per stage
    real(wp), allocatable, dimension(:,:,:) :: qfl, qfr, flx
    $:GPU_DECLARE(create='[qfl, qfr, flx]')

contains

    impure subroutine s_initialize_rhs_staggered_module

        @:ALLOCATE(qfl(-3:m + 3, -3:n + 3, -3:p + 3))
        @:ALLOCATE(qfr(-3:m + 3, -3:n + 3, -3:p + 3))
        @:ALLOCATE(flx(-3:m + 3, -3:n + 3, -3:p + 3))

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
                    do j = -1, m + 1
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

    impure subroutine s_finalize_rhs_staggered_module

        @:DEALLOCATE(qfl, qfr, flx)

    end subroutine s_finalize_rhs_staggered_module

end module m_rhs_staggered
