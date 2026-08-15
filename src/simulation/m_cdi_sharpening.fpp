!>
!! @file
!! @brief Contains module m_cdi_sharpening

#:include 'macros.fpp'
#:include 'case.fpp'

!> @brief Conservative diffuse-interface (CDI) sharpening for the five-equation model (int_comp=3). Adds divergence-form
!! regularization fluxes to the RHS at every Runge-Kutta stage that drive material interfaces to a fixed equilibrium thickness
!! O(eps) while conserving phase mass, mixture momentum, and total energy. The volume-fraction flux for phase m is the N-phase
!! pairwise CDI flux a_m = Gamma*(eps*grad(alpha_m) - sum_{j/=m} alpha_m*alpha_j*nhat_mj), with consistency fluxes rho_m*a_m
!! (continuity), u*sum(rho_m*a_m) (momentum), and sum(a_m*(0.5*rho_m*|u|^2 + (rho*e)_m)) (energy). The energy flux carries phase
!! internal energy, not enthalpy, which preserves pressure/temperature/velocity equilibrium across interfaces. With surface tension,
!! the color function receives the same (two-phase) CDI flux so it stays co-located with the sharpened volume fraction; this term is
!! purely kinematic since c carries no mass and sigma does not enter the pressure inversion. References: S. R. Brill, B. J. Olson,
!! and G. T. Bokman, JCP 542 (2025) 114366 (Eqs. 38-40, 68); S. S. Jain et al., JCP 475 (2023) 111866 (divergence-form approach); S.
!! Mirjalili and A. Mani, JCP 498 (2024) 112657 (N-phase pairwise formulation).
module m_cdi_sharpening

    use m_derived_types
    use m_global_parameters
    use m_mpi_common, only: s_mpi_allreduce_max
    use m_variables_conversion, only: gammas, pi_infs, qvs
    use m_helper_basic, only: f_is_default

    implicit none

    private; public :: s_initialize_cdi_sharpening_module, s_compute_cdi_sharpening_rhs, s_finalize_cdi_sharpening_module

    !> Face-centered regularization fluxes for equations 1..eqn_idx%adv%end, reused per direction
    real(wp), allocatable, dimension(:,:,:,:) :: cdi_flux
    $:GPU_DECLARE(create='[cdi_flux]')

    !> Sharpening velocity scale Gamma: global max |u|, or ic_gamma when set by the user
    real(wp) :: cdi_gamma
    $:GPU_DECLARE(create='[cdi_gamma]')

contains

    !> @brief Pairwise volume fraction alpha_1/(alpha_1 + alpha_2), with clamped inputs
    pure function f_pair_frac(a1, a2) result(res)

        $:GPU_ROUTINE(parallelism='[seq]')
        real(wp), intent(in) :: a1, a2
        real(wp)             :: c1, c2, res

        c1 = min(max(a1, 0._wp), 1._wp)
        c2 = min(max(a2, 0._wp), 1._wp)
        res = c1/(c1 + c2 + sgm_eps)

    end function f_pair_frac

    !> @brief Allocate the CDI sharpening module arrays
    impure subroutine s_initialize_cdi_sharpening_module()

        integer :: flux_end

        flux_end = eqn_idx%adv%end
        if (surface_tension) flux_end = eqn_idx%c

        @:ALLOCATE(cdi_flux(idwbuff(1)%beg:idwbuff(1)%end, idwbuff(2)%beg:idwbuff(2)%end, idwbuff(3)%beg:idwbuff(3)%end, &
                   & 1:flux_end))

    end subroutine s_initialize_cdi_sharpening_module

    !> @brief Compute the CDI sharpening fluxes at cell faces and accumulate their divergence into the RHS
    impure subroutine s_compute_cdi_sharpening_rhs(q_prim_vf, rhs_vf)

        type(scalar_field), dimension(sys_size), intent(in)    :: q_prim_vf
        type(scalar_field), dimension(sys_size), intent(inout) :: rhs_vf

        #:if not MFC_CASE_OPTIMIZATION and USING_AMD
            real(wp), dimension(num_fluids_max) :: af_L, af_R, af_F, rho_F, sharp_t, a_reg
            real(wp), dimension(3)              :: vel_F
        #:else
            real(wp), dimension(num_fluids) :: af_L, af_R, af_F, rho_F, sharp_t, a_reg
            real(wp), dimension(num_vels)   :: vel_F
        #:endif
        real(wp) :: eps_face, gn, g1, g2, rmag, tpair, pres_F, velsq, flux_sum, vel_max_loc, cf_L, cf_R, cf_F
        integer  :: i, j, k, l, q1, q2, iq1, iq2

        ! Velocity scale Gamma = global max |u| unless the user set ic_gamma
        if (f_is_default(ic_gamma)) then
            vel_max_loc = 0._wp
            $:GPU_PARALLEL_LOOP(collapse=3, private='[i, j, k, l, velsq]', reduction='[[vel_max_loc]]', reductionOp='[max]')
            do l = 0, p
                do k = 0, n
                    do j = 0, m
                        velsq = 0._wp
                        $:GPU_LOOP(parallelism='[seq]')
                        do i = eqn_idx%mom%beg, eqn_idx%mom%end
                            velsq = velsq + q_prim_vf(i)%sf(j, k, l)*q_prim_vf(i)%sf(j, k, l)
                        end do
                        vel_max_loc = max(vel_max_loc, velsq)
                    end do
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()

            vel_max_loc = sqrt(vel_max_loc)

            if (num_procs == 1) then
                cdi_gamma = vel_max_loc
            else
                call s_mpi_allreduce_max(vel_max_loc, cdi_gamma)
            end if
        else
            cdi_gamma = ic_gamma
        end if
        $:GPU_UPDATE(device='[cdi_gamma]')

        #! Direction table: face (j, k, l) sits between a cell and its +1 neighbor along the normal. TPL builds index strings
        #! with offsets along the normal (n) and the two transverse directions (a, b). RTG guards inactive directions at
        #! runtime; T1G/T2G guard the transverse central differences (always true for directions the pass itself requires).
        #:for NORM_DIR, XYZ, TPL, RTG, NCC, NIX, DXI, T1CC, T1IX, T1G, T2CC, T2IX, T2G, JB, KB, LB in &
                [(1, 'x', 'j{n}, k{a}, l{b}', 'm > 0', 'x_cc', 'j', 'dx(j)', 'y_cc', 'k', 'n > 0', 'z_cc', 'l', 'p > 0', &
                  '-1, m', '0, n', '0, p'), &
                 (2, 'y', 'j{a}, k{n}, l{b}', 'n > 0', 'y_cc', 'k', 'dy(k)', 'x_cc', 'j', 'm > 0', 'z_cc', 'l', 'p > 0', &
                  '0, m', '-1, n', '0, p'), &
                 (3, 'z', 'j{a}, k{b}, l{n}', 'p > 0', 'z_cc', 'l', 'dz(l)', 'x_cc', 'j', 'm > 0', 'y_cc', 'k', 'n > 0', &
                  '0, m', '0, n', '-1, p')]
            #:set IX = lambda n='', a='', b='': TPL.format(n=n, a=a, b=b)
            #! Skip y/z passes entirely in case-optimized 1D/2D builds
            #:if NORM_DIR == 1 or not MFC_CASE_OPTIMIZATION or num_dims >= NORM_DIR
                if (${RTG}$) then
                    ! ${XYZ}$-direction face fluxes
                    $:GPU_PARALLEL_LOOP(collapse=3, private='[i, j, k, l, q1, q2, iq1, iq2, af_L, af_R, af_F, rho_F, sharp_t, &
                                        & a_reg, vel_F, eps_face, gn, g1, g2, rmag, tpair, pres_F, velsq, flux_sum, cf_L, cf_R, cf_F]')
                    do l = ${LB}$
                        do k = ${KB}$
                            do j = ${JB}$
                                eps_face = ${NCC}$(${NIX}$ + 1) - ${NCC}$(${NIX}$)

                                $:GPU_LOOP(parallelism='[seq]')
                                do i = 1, num_fluids
                                    af_L(i) = min(max(q_prim_vf(eqn_idx%adv%beg + i - 1)%sf(${IX()}$), 0._wp), 1._wp)
                                    af_R(i) = min(max(q_prim_vf(eqn_idx%adv%beg + i - 1)%sf(${IX(n=' + 1')}$), 0._wp), 1._wp)
                                    af_F(i) = 5e-1_wp*(af_L(i) + af_R(i))
                                    rho_F(i) = 5e-1_wp*(q_prim_vf(eqn_idx%cont%beg + i - 1)%sf(${IX()}$)/max(af_L(i), &
                                          & sgm_eps) + q_prim_vf(eqn_idx%cont%beg + i - 1)%sf(${IX(n=' + 1')}$)/max(af_R(i), &
                                          & sgm_eps))
                                    sharp_t(i) = 0._wp
                                end do

                                ! Pairwise sharpening: each unordered pair computed once and applied antisymmetrically so that
                                ! sum_m a_m = 0 holds bitwise (volume fraction compatibility)
                                $:GPU_LOOP(parallelism='[seq]')
                                do q1 = 1, num_fluids - 1
                                    $:GPU_LOOP(parallelism='[seq]')
                                    do q2 = q1 + 1, num_fluids
                                        iq1 = eqn_idx%adv%beg + q1 - 1
                                        iq2 = eqn_idx%adv%beg + q2 - 1

                                        gn = (f_pair_frac(q_prim_vf(iq1)%sf(${IX(n=' + 1')}$), &
                                              & q_prim_vf(iq2)%sf(${IX(n=' + 1')}$)) - f_pair_frac(q_prim_vf(iq1)%sf(${IX()}$), &
                                              & q_prim_vf(iq2)%sf(${IX()}$)))/eps_face

                                        g1 = 0._wp
                                        if (${T1G}$) then
                                            g1 = (f_pair_frac(q_prim_vf(iq1)%sf(${IX(a=' + 1')}$), &
                                                  & q_prim_vf(iq2)%sf(${IX(a=' + 1')}$)) &
                                                  & - f_pair_frac(q_prim_vf(iq1)%sf(${IX(a=' - 1')}$), &
                                                  & q_prim_vf(iq2)%sf(${IX(a=' - 1')}$)) &
                                                  & + f_pair_frac(q_prim_vf(iq1)%sf(${IX(n=' + 1', a=' + 1')}$), &
                                                  & q_prim_vf(iq2)%sf(${IX(n=' + 1', a=' + 1')}$)) &
                                                  & - f_pair_frac(q_prim_vf(iq1)%sf(${IX(n=' + 1', a=' - 1')}$), &
                                                  & q_prim_vf(iq2)%sf(${IX(n=' + 1', a=' - 1')}$)))/(2._wp*(${T1CC}$(${T1IX}$ + 1) &
                                                  & - ${T1CC}$(${T1IX}$ - 1)))
                                        end if

                                        g2 = 0._wp
                                        if (${T2G}$) then
                                            g2 = (f_pair_frac(q_prim_vf(iq1)%sf(${IX(b=' + 1')}$), &
                                                  & q_prim_vf(iq2)%sf(${IX(b=' + 1')}$)) &
                                                  & - f_pair_frac(q_prim_vf(iq1)%sf(${IX(b=' - 1')}$), &
                                                  & q_prim_vf(iq2)%sf(${IX(b=' - 1')}$)) &
                                                  & + f_pair_frac(q_prim_vf(iq1)%sf(${IX(n=' + 1', b=' + 1')}$), &
                                                  & q_prim_vf(iq2)%sf(${IX(n=' + 1', b=' + 1')}$)) &
                                                  & - f_pair_frac(q_prim_vf(iq1)%sf(${IX(n=' + 1', b=' - 1')}$), &
                                                  & q_prim_vf(iq2)%sf(${IX(n=' + 1', b=' - 1')}$)))/(2._wp*(${T2CC}$(${T2IX}$ + 1) &
                                                  & - ${T2CC}$(${T2IX}$ - 1)))
                                        end if

                                        rmag = sqrt(gn*gn + g1*g1 + g2*g2)

                                        if (rmag > verysmall) then
                                            tpair = af_F(q1)*af_F(q2)*gn/rmag
                                            sharp_t(q1) = sharp_t(q1) + tpair
                                            sharp_t(q2) = sharp_t(q2) - tpair
                                        end if
                                    end do
                                end do

                                ! Volume fraction and continuity fluxes; eps*d(alpha)/dn with eps = eps_face reduces to af_R - af_L
                                flux_sum = 0._wp
                                $:GPU_LOOP(parallelism='[seq]')
                                do i = 1, num_fluids
                                    a_reg(i) = cdi_gamma*((af_R(i) - af_L(i)) - sharp_t(i))
                                    cdi_flux(j, k, l, eqn_idx%adv%beg + i - 1) = a_reg(i)
                                    cdi_flux(j, k, l, eqn_idx%cont%beg + i - 1) = rho_F(i)*a_reg(i)
                                    flux_sum = flux_sum + rho_F(i)*a_reg(i)
                                end do

                                ! Momentum consistency flux
                                velsq = 0._wp
                                $:GPU_LOOP(parallelism='[seq]')
                                do i = 1, num_vels
                                    vel_F(i) = 5e-1_wp*(q_prim_vf(eqn_idx%mom%beg + i - 1)%sf(${IX()}$) &
                                          & + q_prim_vf(eqn_idx%mom%beg + i - 1)%sf(${IX(n=' + 1')}$))
                                    cdi_flux(j, k, l, eqn_idx%mom%beg + i - 1) = vel_F(i)*flux_sum
                                    velsq = velsq + vel_F(i)*vel_F(i)
                                end do

                                ! Energy consistency flux: kinetic + phase internal energy (not enthalpy), Brill et al. Eqs. 38-40
                                pres_F = 5e-1_wp*(q_prim_vf(eqn_idx%E)%sf(${IX()}$) + q_prim_vf(eqn_idx%E)%sf(${IX(n=' + 1')}$))
                                flux_sum = 0._wp
                                $:GPU_LOOP(parallelism='[seq]')
                                do i = 1, num_fluids
                                    flux_sum = flux_sum + a_reg(i)*(5e-1_wp*rho_F(i)*velsq + gammas(i)*pres_F + pi_infs(i) &
                                                                & + rho_F(i)*qvs(i))
                                end do
                                cdi_flux(j, k, l, eqn_idx%E) = flux_sum

                                ! Color function sharpening: same two-phase CDI flux, keeping c co-located with the sharpened
                                ! volume fraction. Kinematic only: c carries no mass and sigma does not enter the pressure
                                ! inversion, so no consistency terms are needed.
                                if (surface_tension) then
                                    cf_L = min(max(q_prim_vf(eqn_idx%c)%sf(${IX()}$), 0._wp), 1._wp)
                                    cf_R = min(max(q_prim_vf(eqn_idx%c)%sf(${IX(n=' + 1')}$), 0._wp), 1._wp)
                                    cf_F = 5e-1_wp*(cf_L + cf_R)

                                    gn = (cf_R - cf_L)/eps_face

                                    g1 = 0._wp
                                    if (${T1G}$) then
                                        g1 = (q_prim_vf(eqn_idx%c)%sf(${IX(a=' + 1')}$) &
                                              & - q_prim_vf(eqn_idx%c)%sf(${IX(a=' - 1')}$) &
                                              & + q_prim_vf(eqn_idx%c)%sf(${IX(n=' + 1', a=' + 1')}$) &
                                              & - q_prim_vf(eqn_idx%c)%sf(${IX(n=' + 1', a=' - 1')}$))/(2._wp*(${T1CC}$(${T1IX}$ &
                                              & + 1) - ${T1CC}$(${T1IX}$ - 1)))
                                    end if

                                    g2 = 0._wp
                                    if (${T2G}$) then
                                        g2 = (q_prim_vf(eqn_idx%c)%sf(${IX(b=' + 1')}$) &
                                              & - q_prim_vf(eqn_idx%c)%sf(${IX(b=' - 1')}$) &
                                              & + q_prim_vf(eqn_idx%c)%sf(${IX(n=' + 1', b=' + 1')}$) &
                                              & - q_prim_vf(eqn_idx%c)%sf(${IX(n=' + 1', b=' - 1')}$))/(2._wp*(${T2CC}$(${T2IX}$ &
                                              & + 1) - ${T2CC}$(${T2IX}$ - 1)))
                                    end if

                                    rmag = sqrt(gn*gn + g1*g1 + g2*g2)

                                    tpair = 0._wp
                                    if (rmag > verysmall) tpair = cf_F*(1._wp - cf_F)*gn/rmag

                                    cdi_flux(j, k, l, eqn_idx%c) = cdi_gamma*((cf_R - cf_L) - tpair)
                                end if
                            end do
                        end do
                    end do
                    $:END_GPU_PARALLEL_LOOP()

                    ! Accumulate +div(flux) into the RHS
                    $:GPU_PARALLEL_LOOP(collapse=3, private='[i, j, k, l]')
                    do l = 0, p
                        do k = 0, n
                            do j = 0, m
                                $:GPU_LOOP(parallelism='[seq]')
                                do i = 1, eqn_idx%adv%end
                                    rhs_vf(i)%sf(j, k, l) = rhs_vf(i)%sf(j, k, l) + (cdi_flux(j, k, l, &
                                           & i) - cdi_flux(${IX(n=' - 1')}$, i))/${DXI}$
                                end do
                                if (surface_tension) then
                                    rhs_vf(eqn_idx%c)%sf(j, k, l) = rhs_vf(eqn_idx%c)%sf(j, k, l) + (cdi_flux(j, k, l, &
                                           & eqn_idx%c) - cdi_flux(${IX(n=' - 1')}$, eqn_idx%c))/${DXI}$
                                end if
                            end do
                        end do
                    end do
                    $:END_GPU_PARALLEL_LOOP()
                end if
            #:endif
        #:endfor

    end subroutine s_compute_cdi_sharpening_rhs

    !> @brief Deallocate the CDI sharpening module arrays
    impure subroutine s_finalize_cdi_sharpening_module()

        @:DEALLOCATE(cdi_flux)

    end subroutine s_finalize_cdi_sharpening_module

end module m_cdi_sharpening
