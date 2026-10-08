!>
!! @file
!! @brief Contains module m_diffusion_sts

#:include 'case.fpp'
#:include 'macros.fpp'

!> @brief RKL2 super-time-stepping (Meyer, Balsara & Aslam, JCP 257, 2014) of heat conduction and species diffusion (diff_sts): s
!! stages of the explicit diffusion operator advance them stably over a step up to (s^2 + s - 2)/4 times its forward-Euler limit.
!! Run from the step's start state at fixed density and momentum, the result becomes rates (sts_rate) that every stage of the step
!! adds to its right-hand side, under the projection through the pressure equation as the unsplit diffusion is. Left as a
!! constant-volume pressure jump instead, the projection's relaxation of it loses energy.
module m_diffusion_sts

    use m_derived_types
    use m_global_parameters
    use m_rhs, only: s_compute_diffusion_rhs, sts_rate
    use m_ibm, only: ib_markers

    implicit none

    private; public :: s_initialize_diffusion_sts_module, s_diffusion_sts, f_sts_stages, f_sts_bound, &
        & s_finalize_diffusion_sts_module

    !> The stage start Y0, the stage before last, and the operator at Y0, over the diffused equations (energy and species)
    type(scalar_field), allocatable, dimension(:) :: y0, ym2, l0
    $:GPU_DECLARE(create='[y0, ym2, l0]')

    real(wp), parameter :: sts_safety = 0.8_wp  !< Fraction of the RKL2 stability bound a step may use

contains

    impure subroutine s_initialize_diffusion_sts_module

        integer :: i

        @:ALLOCATE(y0(1:sys_size), ym2(1:sys_size), l0(1:sys_size), sts_rate(1:sys_size))
        do i = 1, sys_size
            if (f_diffused(i)) then
                @:ALLOCATE(y0(i)%sf(0:m, 0:n, 0:p), ym2(i)%sf(0:m, 0:n, 0:p), l0(i)%sf(0:m, 0:n, 0:p), sts_rate(i)%sf(0:m, 0:n, &
                           & 0:p))
            else
                @:ALLOCATE(y0(i)%sf(0:0, 0:0, 0:0), ym2(i)%sf(0:0, 0:0, 0:0), l0(i)%sf(0:0, 0:0, 0:0), sts_rate(i)%sf(0:0, 0:0, &
                           & 0:0))
            end if
            @:ACC_SETUP_SFs(y0(i), ym2(i), l0(i), sts_rate(i))
        end do

    end subroutine s_initialize_diffusion_sts_module

    !> Equation i is diffused: the energy, and the species with chemistry
    logical function f_diffused(i)

        integer, intent(in) :: i

        f_diffused = i == eqn_idx%E .or. (chemistry .and. i >= eqn_idx%species%beg .and. i <= eqn_idx%species%end)

    end function f_diffused

    !> Largest thermal CFL number, D dt sum_d 1/dx_d^2 (forward-Euler limit 1/2), that s stages advance stably
    real(wp) function f_sts_bound(s)
        integer, intent(in) :: s

        f_sts_bound = 0.5_wp*sts_safety*real(s*s + s - 2, wp)/4._wp

    end function f_sts_bound

    !> Fewest stages (at least 2) whose bound, f_sts_bound, covers a step of thermal CFL number tcfl
    integer function f_sts_stages(tcfl)

        real(wp), intent(in) :: tcfl

        f_sts_stages = max(2, ceiling(0.5_wp*(sqrt(9._wp + 32._wp*tcfl/sts_safety) - 1._wp)))

    end function f_sts_stages

    !> RKL2 coefficient b_j
    pure real(wp) function f_b(j)
        integer, intent(in) :: j

        f_b = 1._wp/3._wp
        if (j > 2) f_b = real(j*j + j - 2, wp)/real(2*j*(j + 1), wp)

    end function f_b

    !> The diffusion rates of q_cons_vf over dt (sts_rate), from s RKL2 stages run in place and then undone (rhs_vf is scratch);
    !! fluid cells only, so the immersed-boundary ghost cells keep the boundary state they hold
    impure subroutine s_diffusion_sts(q_cons_vf, q_T_sf, bc_type, pb_in, mv_in, rhs_vf, dt_in, s)

        type(scalar_field), dimension(sys_size), intent(inout)                                     :: q_cons_vf
        type(scalar_field), intent(inout)                                                          :: q_T_sf
        type(integer_field), dimension(1:num_dims,1:2), intent(in)                                 :: bc_type
        real(stp), dimension(idwbuff(1)%beg:,idwbuff(2)%beg:,idwbuff(3)%beg:,1:,1:), intent(inout) :: pb_in, mv_in
        type(scalar_field), dimension(sys_size), intent(inout)                                     :: rhs_vf
        real(wp), intent(in)                                                                       :: dt_in
        integer, intent(in)                                                                        :: s
        real(wp)                                                                                   :: w1, mu, nu, mt, gt, td, y
        logical                                                                                    :: fluid
        integer                                                                                    :: i, j, k, l, st, i2

        w1 = 4._wp/real(s*s + s - 2, wp)
        td = dt_in
        i2 = merge(eqn_idx%species%end, eqn_idx%E, chemistry)

        ! Y0 and L(Y0); Y1 = Y0 + w1/3 dt L(Y0)
        call s_compute_diffusion_rhs(q_cons_vf, q_T_sf, bc_type, pb_in, mv_in, rhs_vf)
        mt = f_b(1)*w1
        $:GPU_PARALLEL_LOOP(collapse=3, private='[i, j, k, l, fluid]', firstprivate='[mt, td, i2]')
        do l = 0, p
            do k = 0, n
                do j = 0, m
                    fluid = .true.
                    if (ib) fluid = ib_markers%sf(j, k, l) == 0
                    $:GPU_LOOP(parallelism='[seq]')
                    do i = eqn_idx%E, i2
                        if (i == eqn_idx%E .or. i >= eqn_idx%species%beg) then
                            y0(i)%sf(j, k, l) = q_cons_vf(i)%sf(j, k, l)
                            ym2(i)%sf(j, k, l) = q_cons_vf(i)%sf(j, k, l)
                            l0(i)%sf(j, k, l) = rhs_vf(i)%sf(j, k, l)
                            if (fluid) q_cons_vf(i)%sf(j, k, l) = real(real(q_cons_vf(i)%sf(j, k, l), &
                                & wp) + mt*td*real(rhs_vf(i)%sf(j, k, l), wp), stp)
                        end if
                    end do
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

        ! Y_j = mu_j Y_j-1 + nu_j Y_j-2 + (1 - mu_j - nu_j) Y0 + mu~_j dt L(Y_j-1) + gamma~_j dt L(Y0)
        do st = 2, s
            mu = real(2*st - 1, wp)/real(st, wp)*f_b(st)/f_b(st - 1)
            nu = -real(st - 1, wp)/real(st, wp)*f_b(st)/f_b(st - 2)
            mt = mu*w1
            gt = -(1._wp - f_b(st - 1))*mt
            call s_compute_diffusion_rhs(q_cons_vf, q_T_sf, bc_type, pb_in, mv_in, rhs_vf)
            $:GPU_PARALLEL_LOOP(collapse=3, private='[i, j, k, l, y, fluid]', firstprivate='[mu, nu, mt, gt, td, i2]')
            do l = 0, p
                do k = 0, n
                    do j = 0, m
                        fluid = .true.
                        if (ib) fluid = ib_markers%sf(j, k, l) == 0
                        if (fluid) then
                            $:GPU_LOOP(parallelism='[seq]')
                            do i = eqn_idx%E, i2
                                if (i == eqn_idx%E .or. i >= eqn_idx%species%beg) then
                                    y = real(q_cons_vf(i)%sf(j, k, l), wp)
                                    q_cons_vf(i)%sf(j, k, l) = real(mu*y + nu*real(ym2(i)%sf(j, k, l), &
                                              & wp) + (1._wp - mu - nu)*real(y0(i)%sf(j, k, l), wp) + mt*td*real(rhs_vf(i)%sf(j, &
                                              & k, l), wp) + gt*td*real(l0(i)%sf(j, k, l), wp), stp)
                                    ym2(i)%sf(j, k, l) = real(y, stp)
                                end if
                            end do
                        end if
                    end do
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()
        end do

        $:GPU_PARALLEL_LOOP(collapse=3, private='[i, j, k, l]', firstprivate='[td, i2]')
        do l = 0, p
            do k = 0, n
                do j = 0, m
                    $:GPU_LOOP(parallelism='[seq]')
                    do i = eqn_idx%E, i2
                        if (i == eqn_idx%E .or. i >= eqn_idx%species%beg) then
                            sts_rate(i)%sf(j, k, l) = real((real(q_cons_vf(i)%sf(j, k, l), wp) - real(y0(i)%sf(j, k, l), wp))/td, &
                                     & stp)
                            q_cons_vf(i)%sf(j, k, l) = y0(i)%sf(j, k, l)
                        end if
                    end do
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_diffusion_sts

    impure subroutine s_finalize_diffusion_sts_module

        integer :: i

        do i = 1, sys_size
            @:DEALLOCATE(y0(i)%sf, ym2(i)%sf, l0(i)%sf, sts_rate(i)%sf)
        end do
        @:DEALLOCATE(y0, ym2, l0, sts_rate)

    end subroutine s_finalize_diffusion_sts_module

end module m_diffusion_sts
