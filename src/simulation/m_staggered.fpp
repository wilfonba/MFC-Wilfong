!>
!! @file m_staggered.fpp
!! @brief Staggered (MAC) grid state and difference operators

#:include 'case.fpp'
#:include 'macros.fpp'

!> Staggered, Harlow-Welch (MAC) placement of the unknowns: scalars stay at cell centres and each velocity component lives on the
!! faces normal to its own direction. The collocated grid MFC otherwise uses cannot compose a divergence, a gradient and a Laplacian
!! consistently -- its divergence and gradient are both 2*dx wide, so composing them gives the WIDE Laplacian, which is blind to a
!! checkerboard, while the natural Laplacian is compact. That single mismatch defeated the semi-implicit projection method as a
!! scheme (examples/1D_contact_semiimplicit/README.md) and every preconditioner tried for the implicit solver since.
!!
!! Staggering removes it by construction: the divergence of a face field and the gradient
!! of a cell field are exact negative adjoints, and composing them gives the compact
!! Laplacian with no null space. s_staggered_self_test measures exactly that, and is the
!! acceptance gate for this discretization -- it answers, before any physics is written,
!! the question the whole approach turns on.
!!
!! A second property matters just as much at a high density ratio. Face velocity is a
!! primary unknown here, so the light phase's momentum is never divided by the light
!! phase's density next to a heavy cell. On the collocated grid that quotient gives the
!! Jacobian an entry of order rho_heavy*c_heavy^2/rho_light ~ 1e9 across an 833:1
!! interface, which no preconditioner can match closely enough to help.
!!
!! Index convention: for direction d, index j of a face array means the face between
!! cells j and j+1, i.e. x_{j+1/2}. Cell j is therefore bounded by faces j-1 and j. The
!! face arrays carry all three spatial indices at full cell extent, which over-allocates
!! the two directions a given component is not staggered in; that is the standard MAC
!! layout choice and keeps every kernel a plain triple loop.
!> @brief Staggered (MAC) grid state and difference operators
module m_staggered

    use m_derived_types
    use m_global_parameters
    use m_mpi_proxy

    implicit none

    private; public :: s_initialize_staggered_module, s_cell_to_face, s_face_to_cell, s_face_divergence, s_cell_gradient, &
        & s_face_density, s_staggered_self_test, s_finalize_staggered_module

    !> Face-centred velocity and the harmonic face density that goes with it
    real(wp), allocatable, dimension(:,:,:,:) :: uf  !< face velocity, one component per direction
    real(wp), allocatable, dimension(:,:,:,:) :: rf  !< harmonic face density
    $:GPU_DECLARE(create='[uf, rf]')

contains

    !> Allocate the face state. The self-test is a separate entry point because it reads the grid spacings, which are not filled
    !! until the data files have been read
    impure subroutine s_initialize_staggered_module

        @:ALLOCATE(uf(-3:m + 3, -3:n + 3, -3:p + 3, 1:num_dims))
        @:ALLOCATE(rf(-3:m + 3, -3:n + 3, -3:p + 3, 1:num_dims))

    end subroutine s_initialize_staggered_module

    !> Arithmetic average of a cell field onto the faces of direction `d`. Used for the quantities that are conserved on the shifted
    !! control volume, where the average is what keeps the momentum sum telescoping
    subroutine s_cell_to_face(qc, qfd, d)

        real(wp), dimension(-3:m + 3,-3:n + 3,-3:p + 3), intent(in)  :: qc
        real(wp), dimension(-3:m + 3,-3:n + 3,-3:p + 3), intent(out) :: qfd
        integer, intent(in)                                          :: d
        integer                                                      :: j, k, l

        #:for DIM, IP1 in [(1, 'j + 1, k, l'), (2, 'j, k + 1, l'), (3, 'j, k, l + 1')]
            if (d == ${DIM}$) then
                $:GPU_PARALLEL_LOOP(collapse=3, private='[j, k, l]')
                do l = -2, p
                    do k = -2, n
                        do j = -2, m
                            qfd(j, k, l) = 0.5_wp*(qc(j, k, l) + qc(${IP1}$))
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()
            end if
        #:endfor

    end subroutine s_cell_to_face

    !> Average the two faces bounding a cell back to the cell centre
    subroutine s_face_to_cell(qfd, qc, d)

        real(wp), dimension(-3:m + 3,-3:n + 3,-3:p + 3), intent(in)  :: qfd
        real(wp), dimension(-3:m + 3,-3:n + 3,-3:p + 3), intent(out) :: qc
        integer, intent(in)                                          :: d
        integer                                                      :: j, k, l

        #:for DIM, IM1 in [(1, 'j - 1, k, l'), (2, 'j, k - 1, l'), (3, 'j, k, l - 1')]
            if (d == ${DIM}$) then
                $:GPU_PARALLEL_LOOP(collapse=3, private='[j, k, l]')
                do l = -1, p + 1
                    do k = -1, n + 1
                        do j = -1, m + 1
                            qc(j, k, l) = 0.5_wp*(qfd(j, k, l) + qfd(${IM1}$))
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()
            end if
        #:endfor

    end subroutine s_face_to_cell

    !> Divergence of a face field, cell by cell. Compact by construction: cell j sees only the two faces that bound it, so there is
    !! no odd-even decoupling to inherit
    subroutine s_face_divergence(ufin, dvg)

        real(wp), dimension(-3:m + 3,-3:n + 3,-3:p + 3,1:num_dims), intent(in) :: ufin
        real(wp), dimension(-3:m + 3,-3:n + 3,-3:p + 3), intent(out)           :: dvg
        integer                                                                :: j, k, l

        $:GPU_PARALLEL_LOOP(collapse=3, private='[j, k, l]')
        do l = 0, p
            do k = 0, n
                do j = 0, m
                    dvg(j, k, l) = 0._wp
                    #:for DIM, DXV, IDXV, IM1 in [(1, 'dx', 'j', 'j - 1, k, l'), &
                        (2, 'dy', 'k', 'j, k - 1, l'), (3, 'dz', 'l', 'j, k, l - 1')]
                        if (num_dims >= ${DIM}$) then
                            dvg(j, k, l) = dvg(j, k, l) + (ufin(j, k, l, ${DIM}$) - ufin(${IM1}$, ${DIM}$))/${DXV}$(${IDXV}$)
                        end if
                    #:endfor
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_face_divergence

    !> Gradient of a cell field onto the faces of direction `d`, over the centre-to-centre distance. This is the exact negative
    !! adjoint of s_face_divergence
    subroutine s_cell_gradient(qc, gfd, d)

        real(wp), dimension(-3:m + 3,-3:n + 3,-3:p + 3), intent(in)  :: qc
        real(wp), dimension(-3:m + 3,-3:n + 3,-3:p + 3), intent(out) :: gfd
        integer, intent(in)                                          :: d
        integer                                                      :: j, k, l

        #:for DIM, DXV, IDXV, IP1, IPX in [(1, 'dx', 'j', 'j + 1, k, l', 'j + 1'), &
            (2, 'dy', 'k', 'j, k + 1, l', 'k + 1'), (3, 'dz', 'l', 'j, k, l + 1', 'l + 1')]
            if (d == ${DIM}$) then
                $:GPU_PARALLEL_LOOP(collapse=3, private='[j, k, l]')
                do l = -2, p
                    do k = -2, n
                        do j = -2, m
                            gfd(j, k, l) = (qc(${IP1}$) - qc(j, k, l))/(0.5_wp*(${DXV}$(${IDXV}$) + ${DXV}$(${IPX}$)))
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()
            end if
        #:endfor

    end subroutine s_cell_gradient

    !> Harmonic face density. Across a large jump the harmonic mean is the flux-continuous choice -- the arithmetic one detonated
    !! the projection method's Helmholtz solve at high acoustic CFL where the harmonic one held to machine precision
    subroutine s_face_density(rhoc, rfd, d)

        real(wp), dimension(-3:m + 3,-3:n + 3,-3:p + 3), intent(in)  :: rhoc
        real(wp), dimension(-3:m + 3,-3:n + 3,-3:p + 3), intent(out) :: rfd
        integer, intent(in)                                          :: d
        integer                                                      :: j, k, l

        #:for DIM, IP1 in [(1, 'j + 1, k, l'), (2, 'j, k + 1, l'), (3, 'j, k, l + 1')]
            if (d == ${DIM}$) then
                $:GPU_PARALLEL_LOOP(collapse=3, private='[j, k, l]')
                do l = -2, p
                    do k = -2, n
                        do j = -2, m
                            rfd(j, k, l) = 2._wp*rhoc(j, k, l)*rhoc(${IP1}$)/max(rhoc(j, k, l) + rhoc(${IP1}$), sgm_eps)
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()
            end if
        #:endfor

    end subroutine s_face_density

    !> Acceptance gate for the staggered operators, run once at initialization.
    !!
    !! Three properties, in increasing order of what they buy:
    !!
    !!   1. Adjointness. Summation by parts on a periodic grid gives exactly
    !!        sum_c V_c p_c div(u)_c = - sum_f V_f u_f grad(p)_f,
    !!      with V_c the cell volume and V_f the face volume, for ANY p and u. This is
    !!      the discrete counterpart of integration by parts, and it is what makes the
    !!      pressure operator symmetric. It should hold to round-off.
    !!
    !!   2. Compactness. Composing the two gives the three-point Laplacian, not the wide
    !!      five-point one the collocated grid produces.
    !!
    !!   3. No checkerboard null space. The wide operator annihilates (-1)^(j+k+l)
    !!      entirely; the compact one returns 4/dx^2 per dimension on it. This is the
    !!      failure mode that has defeated every previous attempt, so it is checked
    !!      directly rather than inferred.
    impure subroutine s_staggered_self_test

        real(wp), allocatable, dimension(:,:,:)   :: pc, lap, chk
        real(wp), allocatable, dimension(:,:,:,:) :: gf, vf
        real(wp)                                  :: lhs, rhs, vc, vfa, adj, lapmax, chkmin, expect, sc, qx, qy, qz
        integer                                   :: j, k, l, d

        allocate (pc(-3:m + 3,-3:n + 3,-3:p + 3), lap(-3:m + 3,-3:n + 3,-3:p + 3))
        allocate (chk(-3:m + 3,-3:n + 3,-3:p + 3))
        allocate (gf(-3:m + 3,-3:n + 3,-3:p + 3,1:num_dims))
        allocate (vf(-3:m + 3,-3:n + 3,-3:p + 3,1:num_dims))

        ! Two unrelated fields with content down to the grid scale, so the identity is tested on something demanding rather than on
        ! a mode the operators happen to like. Both are exactly periodic in the index, which is what makes the summation by parts
        ! telescope with no boundary terms left over -- a test field that does not wrap fails this identity for a reason that has
        ! nothing to do with the operators
        qx = 2._wp*pi/real(m + 1, wp)
        qy = 2._wp*pi/real(n + 1, wp)
        qz = 2._wp*pi/real(p + 1, wp)
        do l = -3, p + 3
            do k = -3, n + 3
                do j = -3, m + 3
                    pc(j, k, l) = sin(3._wp*qx*j) + 0.3_wp*cos(7._wp*qx*j)
                    if (num_dims > 1) pc(j, k, l) = pc(j, k, l) + 0.7_wp*sin(2._wp*qy*k)
                    if (num_dims > 2) pc(j, k, l) = pc(j, k, l) + 0.5_wp*sin(qz*l)
                    chk(j, k, l) = real(1 - 2*modulo(j + k + l, 2), wp)
                    do d = 1, num_dims
                        vf(j, k, l, d) = cos(3._wp*qx*j + 0.5_wp*d) + 0.6_wp*sin(7._wp*qx*j)
                        if (num_dims > 1) vf(j, k, l, d) = vf(j, k, l, d) + 0.4_wp*cos(3._wp*qy*k)
                        if (num_dims > 2) vf(j, k, l, d) = vf(j, k, l, d) + 0.2_wp*cos(2._wp*qz*l)
                    end do
                end do
            end do
        end do

        ! 1. adjointness. Both sums are taken over the interior only, which is exact for the periodic wrap the test fields provide
        call s_face_divergence(vf, lap)
        lhs = 0._wp; sc = 0._wp
        do l = 0, p
            do k = 0, n
                do j = 0, m
                    vc = dx(j)
                    if (num_dims > 1) vc = vc*dy(k)
                    if (num_dims > 2) vc = vc*dz(l)
                    lhs = lhs + vc*pc(j, k, l)*lap(j, k, l)
                    sc = sc + abs(vc*pc(j, k, l)*lap(j, k, l))
                end do
            end do
        end do

        rhs = 0._wp
        do d = 1, num_dims
            call s_cell_gradient(pc, gf(:,:,:,d), d)
        end do
        do l = 0, p
            do k = 0, n
                do j = 0, m
                    #:for DIM, DXV, IDXV, IPX in [(1, 'dx', 'j', 'j + 1'), (2, 'dy', 'k', 'k + 1'), &
                        (3, 'dz', 'l', 'l + 1')]
                        if (num_dims >= ${DIM}$) then
                            vfa = 0.5_wp*(${DXV}$(${IDXV}$) + ${DXV}$(${IPX}$))
                            #:for OD, ODV, ODI in [(1, 'dx', 'j'), (2, 'dy', 'k'), (3, 'dz', 'l')]
                                if (num_dims >= ${OD}$ .and. ${OD}$ /= ${DIM}$) vfa = vfa*${ODV}$(${ODI}$)
                            #:endfor
                            rhs = rhs + vfa*vf(j, k, l, ${DIM}$)*gf(j, k, l, ${DIM}$)
                            sc = sc + abs(vfa*vf(j, k, l, ${DIM}$)*gf(j, k, l, ${DIM}$))
                        end if
                    #:endfor
                end do
            end do
        end do

        ! Divide by the accumulated size of the contributions, not by either sum. Each sum
        ! can vanish on its own when the test modes happen to be orthogonal, and
        ! normalising by a zero turns an exact result into an apparent failure
        adj = abs(lhs + rhs)/max(sc, sgm_eps)

        ! 2 and 3. Compose div(grad(.)) and look at what it does to a checkerboard. The
        ! wide operator the collocated grid produces returns exactly zero here
        do d = 1, num_dims
            call s_cell_gradient(chk, gf(:,:,:,d), d)
        end do
        call s_face_divergence(gf, lap)
        lapmax = 0._wp; chkmin = huge(1._wp)
        do l = 0, p
            do k = 0, n
                do j = 0, m
                    lapmax = max(lapmax, abs(lap(j, k, l)))
                    chkmin = min(chkmin, abs(lap(j, k, l)))
                end do
            end do
        end do

        expect = 0._wp
        #:for DIM, DXV, IDXV in [(1, 'dx', 'j'), (2, 'dy', 'k'), (3, 'dz', 'l')]
            if (num_dims >= ${DIM}$) expect = expect + 4._wp/(${DXV}$(0)*${DXV}$(0))
        #:endfor

        if (proc_rank == 0) then
            print '(A)', ' Staggered operator self-test'
            print '(A, ES12.5)', '   adjointness  |<p, div u> + <grad p, u>| / scale  = ', adj
            print '(A, ES12.5, A, ES12.5)', '   checkerboard div(grad) min ', chkmin, '  max ', lapmax
            print '(A, ES12.5)', '   checkerboard expected (compact, uniform grid)    = ', expect
        end if

        deallocate (pc, lap, chk, gf, vf)

    end subroutine s_staggered_self_test

    impure subroutine s_finalize_staggered_module

        @:DEALLOCATE(uf, rf)

    end subroutine s_finalize_staggered_module

end module m_staggered
