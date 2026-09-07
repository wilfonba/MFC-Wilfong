!>
!! @file m_weno_staggered.fpp
!! @brief WENO reconstruction of cell scalars onto staggered faces

#:include 'case.fpp'
#:include 'macros.fpp'

!> Fifth-order WENO reconstruction of a cell-centred scalar onto the faces of one direction, producing the left- and right-biased
!! states that an upwind flux selects between using the face velocity.
!!
!! This is deliberately separate from m_weno. That module reconstructs into the reshaped
!! Riemann buffers of the collocated pipeline, whose layout and setup are entangled with
!! the collocated RHS; the staggered path wants the same one-dimensional operation
!! applied directly to a cell array, with the face index meaning what it means everywhere
!! else here -- face j lies between cells j and j+1. Merging the two is worth revisiting
!! once the staggered scheme is doing real work, not before it is known to work at all.
!!
!! Stencil reach: the left state at face j needs cells j-2 .. j+2 and the right state
!! needs cells j-1 .. j+3, so faces -1 .. m require cells -3 .. m+3. That is why the
!! staggered arrays carry three ghost layers.
!> @brief WENO reconstruction of cell scalars onto staggered faces
module m_weno_staggered

    use m_derived_types
    use m_global_parameters

    implicit none

    private; public :: s_weno_face_states

    !> WENO5 optimal weights and the smoothness-indicator floor
    real(wp), parameter :: wd0 = 0.1_wp, wd1 = 0.6_wp, wd2 = 0.3_wp
    real(wp), parameter :: w13 = 13._wp/12._wp
    real(wp), parameter :: weps = 1.e-40_wp

contains

    !> Left- and right-biased face states in direction `d`. `qfl(j)` is the state on the low side of face j, reconstructed from
    !! cells j-2..j+2; `qfr(j)` is the state on the high side, reconstructed from cells j+3..j-1. The two are mirror images, so the
    !! same five-point formula serves both with its arguments reversed
    subroutine s_weno_face_states(qc, qfl, qfr, d)

        real(wp), dimension(-3:m + 3,-3:n + 3,-3:p + 3), intent(in)  :: qc
        real(wp), dimension(-3:m + 3,-3:n + 3,-3:p + 3), intent(out) :: qfl, qfr
        integer, intent(in)                                          :: d
        integer                                                      :: j, k, l

        #:for DIM, M2, M1, P1, P2, P3 in [ &
            (1, 'j - 2, k, l', 'j - 1, k, l', 'j + 1, k, l', 'j + 2, k, l', 'j + 3, k, l'), &
            (2, 'j, k - 2, l', 'j, k - 1, l', 'j, k + 1, l', 'j, k + 2, l', 'j, k + 3, l'), &
            (3, 'j, k, l - 2', 'j, k, l - 1', 'j, k, l + 1', 'j, k, l + 2', 'j, k, l + 3')]
            if (d == ${DIM}$) then
                $:GPU_PARALLEL_LOOP(collapse=3, private='[j, k, l]')
                do l = -1, p + 1
                    do k = -1, n + 1
                        do j = -1, m + 1
                            qfl(j, k, l) = f_weno5(qc(${M2}$), qc(${M1}$), qc(j, k, l), qc(${P1}$), qc(${P2}$))
                            qfr(j, k, l) = f_weno5(qc(${P3}$), qc(${P2}$), qc(${P1}$), qc(j, k, l), qc(${M1}$))
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()
            end if
        #:endfor

    end subroutine s_weno_face_states

    !> Fifth-order WENO-JS reconstruction at the right-hand face of the cell holding v0, from the five consecutive values centred on
    !! it
    pure elemental function f_weno5(vm2, vm1, v0, vp1, vp2) result(vf)

        real(wp), intent(in) :: vm2, vm1, v0, vp1, vp2
        real(wp)             :: vf
        real(wp)             :: p0, p1, p2, b0, b1, b2, a0, a1, a2, at

        ! candidate stencils

        p0 = (2._wp*vm2 - 7._wp*vm1 + 11._wp*v0)/6._wp
        p1 = (-vm1 + 5._wp*v0 + 2._wp*vp1)/6._wp
        p2 = (2._wp*v0 + 5._wp*vp1 - vp2)/6._wp

        ! smoothness indicators
        b0 = w13*(vm2 - 2._wp*vm1 + v0)**2 + 0.25_wp*(vm2 - 4._wp*vm1 + 3._wp*v0)**2
        b1 = w13*(vm1 - 2._wp*v0 + vp1)**2 + 0.25_wp*(vm1 - vp1)**2
        b2 = w13*(v0 - 2._wp*vp1 + vp2)**2 + 0.25_wp*(3._wp*v0 - 4._wp*vp1 + vp2)**2

        a0 = wd0/(weps + b0)**2
        a1 = wd1/(weps + b1)**2
        a2 = wd2/(weps + b2)**2
        at = a0 + a1 + a2

        vf = (a0*p0 + a1*p1 + a2*p2)/at

    end function f_weno5

end module m_weno_staggered
