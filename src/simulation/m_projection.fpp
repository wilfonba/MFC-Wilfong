!>
!! @file
!! @brief Contains module m_projection

#:include 'case.fpp'
#:include 'macros.fpp'

#! Conservative FV Helmholtz stencil row at cell (j, k, l): face coefficients
#! c_f = 1/(rho_face*d_face*d_cell) with face-averaged density, accumulated
#! into the off-diagonal term (reading neighbor pressures from `pfield`) and
#! the diagonal, scaled outside by coeff = rho*c^2*dt^2
#:def PROJECTION_STENCIL(pfield)
    coeff = real(rhoc2_cell(j, k, l), wp)*dt*dt
    offd = 0._wp
    diag = 0._wp
    rho_c = 0._wp
    $:GPU_LOOP(parallelism='[seq]')
    do i = 1, num_fluids
        rho_c = rho_c + real(q_cons_vf(i)%sf(j, k, l), wp)
    end do
    rho_nb = 0._wp
    $:GPU_LOOP(parallelism='[seq]')
    do i = 1, num_fluids
        rho_nb = rho_nb + real(q_cons_vf(i)%sf(j - 1, k, l), wp)
    end do
    c_f = 1._wp/(max(0.5_wp*(rho_c + rho_nb), sgm_eps)*0.5_wp*(dx(j - 1) + dx(j))*dx(j))
    offd = offd + c_f*real(${pfield}$ (j - 1, k, l), wp)
    diag = diag + c_f
    rho_nb = 0._wp
    $:GPU_LOOP(parallelism='[seq]')
    do i = 1, num_fluids
        rho_nb = rho_nb + real(q_cons_vf(i)%sf(j + 1, k, l), wp)
    end do
    c_f = 1._wp/(max(0.5_wp*(rho_c + rho_nb), sgm_eps)*0.5_wp*(dx(j) + dx(j + 1))*dx(j))
    offd = offd + c_f*real(${pfield}$ (j + 1, k, l), wp)
    diag = diag + c_f
    if (num_dims > 1) then
        rho_nb = 0._wp
        $:GPU_LOOP(parallelism='[seq]')
        do i = 1, num_fluids
            rho_nb = rho_nb + real(q_cons_vf(i)%sf(j, k - 1, l), wp)
        end do
        c_f = 1._wp/(max(0.5_wp*(rho_c + rho_nb), sgm_eps)*0.5_wp*(dy(k - 1) + dy(k))*dy(k))
        offd = offd + c_f*real(${pfield}$ (j, k - 1, l), wp)
        diag = diag + c_f
        rho_nb = 0._wp
        $:GPU_LOOP(parallelism='[seq]')
        do i = 1, num_fluids
            rho_nb = rho_nb + real(q_cons_vf(i)%sf(j, k + 1, l), wp)
        end do
        c_f = 1._wp/(max(0.5_wp*(rho_c + rho_nb), sgm_eps)*0.5_wp*(dy(k) + dy(k + 1))*dy(k))
        offd = offd + c_f*real(${pfield}$ (j, k + 1, l), wp)
        diag = diag + c_f
    end if
    if (num_dims > 2) then
        rho_nb = 0._wp
        $:GPU_LOOP(parallelism='[seq]')
        do i = 1, num_fluids
            rho_nb = rho_nb + real(q_cons_vf(i)%sf(j, k, l - 1), wp)
        end do
        c_f = 1._wp/(max(0.5_wp*(rho_c + rho_nb), sgm_eps)*0.5_wp*(dz(l - 1) + dz(l))*dz(l))
        offd = offd + c_f*real(${pfield}$ (j, k, l - 1), wp)
        diag = diag + c_f
        rho_nb = 0._wp
        $:GPU_LOOP(parallelism='[seq]')
        do i = 1, num_fluids
            rho_nb = rho_nb + real(q_cons_vf(i)%sf(j, k, l + 1), wp)
        end do
        c_f = 1._wp/(max(0.5_wp*(rho_c + rho_nb), sgm_eps)*0.5_wp*(dz(l) + dz(l + 1))*dz(l))
        offd = offd + c_f*real(${pfield}$ (j, k, l + 1), wp)
        diag = diag + c_f
    end if
#:enddef

!> @brief Semi-implicit pressure projection method (Kwatra et al.). Pressure is removed from the Riemann flux (advective-only wave
!! speeds), advected explicitly, then solved implicitly each RK stage from the Helmholtz equation p - rho*c^2*dt^2 div(grad(p)/rho)
!! = p_adv - rho*c^2*dt*div(u*), after which momentum is corrected and total energy is rebuilt from the EOS. This lifts the acoustic
!! CFL restriction; the time step is limited by the advective CFL only. Note: the divergence and pressure-gradient stencils are 2*dx
!! wide while the implicit Laplacian is compact, so the projection is not discretely exact; this matches the validated reference
!! discretization.
module m_projection

    use m_derived_types
    use m_global_parameters
    use m_mpi_common
    use m_mpi_proxy
    use m_boundary_common
    use m_body_forces, only: s_compute_acceleration
    use m_surface_tension, only: s_compute_capillary_source_flux
    use m_riemann_state, only: Re_avg_rsx_vf, vel_src_rsx_vf, Res_gs

    implicit none

    private; public :: s_initialize_projection_module, s_projection_directional_rhs, s_projection_add_flux_src, &
        & s_projection_apply, s_finalize_projection_module

    real(stp), allocatable, target, dimension(:,:,:) :: pres_proj      !< pressure iterate of the Helmholtz solve
    real(stp), allocatable, dimension(:,:,:)         :: pres_proj_old  !< previous iterate (Jacobi and Chebyshev)
    real(stp), allocatable, dimension(:,:,:)         :: d_cheb         !< Chebyshev search increment
    real(stp), allocatable, dimension(:,:,:)         :: pres_stage     !< pressure at the start of the current RK stage
    real(stp), allocatable, dimension(:,:,:)         :: pres_step0     !< pressure at the start of the time step (RK2/RK3 only)
    real(stp), allocatable, dimension(:,:,:)         :: rhs_p_adv      !< RHS of the explicit pressure advection equation
    real(stp), allocatable, dimension(:,:,:)         :: div_u_face     !< flux-form div(u) from Riemann face velocities
    real(stp), allocatable, dimension(:,:,:)         :: helm_rhs       !< RHS of the Helmholtz solve
    real(stp), allocatable, dimension(:,:,:)         :: rhoc2_cell     !< mixture rho*c^2 (Wood's sound speed)
    real(stp), allocatable, dimension(:,:,:)         :: flux_face_vel  !< Riemann face velocity (per direction sweep)
    real(stp), allocatable, dimension(:,:,:)         :: flux_pu        !< upwinded p*u face flux (per direction sweep)
    $:GPU_DECLARE(create='[pres_proj, pres_proj_old, d_cheb, pres_stage, pres_step0]')
    $:GPU_DECLARE(create='[rhs_p_adv, div_u_face, helm_rhs, rhoc2_cell]')
    $:GPU_DECLARE(create='[flux_face_vel, flux_pu]')

    type(scalar_field), dimension(1) :: pres_proj_sf
    $:GPU_DECLARE(create='[pres_proj_sf]')

    integer :: rb_offset  !< global parity offset of this rank for red-black coloring
    $:GPU_DECLARE(create='[rb_offset]')

    !> Zero face-velocity array for the capillary flux call: its velocity-work terms only feed the energy flux, which the projection
    !! EOS rebuild discards
    real(wp), allocatable, dimension(:,:,:,:) :: cap_vsrc
    type(int_bounds_info)                     :: cap_isx, cap_isy, cap_isz
    $:GPU_DECLARE(create='[cap_vsrc, cap_isx, cap_isy, cap_isz]')

contains

    !> Initialize the projection module
    impure subroutine s_initialize_projection_module()

        integer :: j, k, l

        @:ALLOCATE(pres_proj(idwbuff(1)%beg:idwbuff(1)%end, idwbuff(2)%beg:idwbuff(2)%end, idwbuff(3)%beg:idwbuff(3)%end))
        if (proj_iter_solver /= proj_iter_solver_gauss_seidel) then
            @:ALLOCATE(pres_proj_old(idwbuff(1)%beg:idwbuff(1)%end, idwbuff(2)%beg:idwbuff(2)%end, idwbuff(3)%beg:idwbuff(3)%end))
        end if
        if (proj_iter_solver == proj_iter_solver_chebyshev) then
            @:ALLOCATE(d_cheb(0:m, 0:n, 0:p))
            d_cheb = 0._stp
            $:GPU_UPDATE(device='[d_cheb]')
        end if

        @:ALLOCATE(pres_stage(0:m, 0:n, 0:p))
        if (time_stepper /= time_stepper_rk1) then
            @:ALLOCATE(pres_step0(0:m, 0:n, 0:p))
        end if
        @:ALLOCATE(rhs_p_adv(0:m, 0:n, 0:p))
        @:ALLOCATE(div_u_face(0:m, 0:n, 0:p))
        @:ALLOCATE(helm_rhs(0:m, 0:n, 0:p))
        @:ALLOCATE(rhoc2_cell(0:m, 0:n, 0:p))

        @:ALLOCATE(flux_face_vel(idwbuff(1)%beg:idwbuff(1)%end, idwbuff(2)%beg:idwbuff(2)%end, idwbuff(3)%beg:idwbuff(3)%end))
        @:ALLOCATE(flux_pu(idwbuff(1)%beg:idwbuff(1)%end, idwbuff(2)%beg:idwbuff(2)%end, idwbuff(3)%beg:idwbuff(3)%end))

        if (surface_tension) then
            @:ALLOCATE(cap_vsrc(-1:m + 1, -1:n + 1, -1:p + 1, 1:num_dims))
            cap_vsrc = 0._wp
            $:GPU_UPDATE(device='[cap_vsrc]')
        end if

        $:GPU_PARALLEL_LOOP(private='[j, k, l]', collapse=3)
        do l = idwbuff(3)%beg, idwbuff(3)%end
            do k = idwbuff(2)%beg, idwbuff(2)%end
                do j = idwbuff(1)%beg, idwbuff(1)%end
                    pres_proj(j, k, l) = 0._stp
                    if (proj_iter_solver /= proj_iter_solver_gauss_seidel) pres_proj_old(j, k, l) = 0._stp
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

        ! Global parity of this rank's first interior cell, so the red-black
        ! coloring is consistent across rank boundaries (start_idx is only set
        ! for parallel I/O, so it is recomputed from the block distribution)
        rb_offset = 0
        if (num_procs > 1) then
            rb_offset = f_global_offset(m_glb, num_procs_x, proc_coords(1))
            if (num_dims > 1) rb_offset = rb_offset + f_global_offset(n_glb, num_procs_y, proc_coords(2))
            if (num_dims > 2) rb_offset = rb_offset + f_global_offset(p_glb, num_procs_z, proc_coords(3))
            rb_offset = mod(rb_offset, 2)
        end if
        $:GPU_UPDATE(device='[rb_offset]')

        pres_proj_sf(1)%sf => pres_proj
        $:GPU_ENTER_DATA(copyin='[pres_proj_sf(1)%sf]')
        $:GPU_ENTER_DATA(attach='[pres_proj_sf(1)%sf]')

    end subroutine s_initialize_projection_module

    !> Global starting cell index of this rank along one direction under MFC's block distribution (same layout rule as the domain
    !! decomposition)
    pure function f_global_offset(glb_cells, nprocs_dir, coord) result(offset)

        integer, intent(in) :: glb_cells, nprocs_dir, coord
        integer             :: offset
        integer             :: base, rem

        base = (glb_cells + 1)/nprocs_dir
        rem = mod(glb_cells + 1, nprocs_dir)
        offset = base*coord + min(coord, rem)

    end function f_global_offset

    !> Compute the pressure-free convective RHS contribution of one direction sweep: an HLLC flux with advective-only wave speeds
    !! and no pressure terms, plus the face velocity (for div(u)) and the upwinded p*u flux (for the explicit pressure advection
    !! equation)
    subroutine s_projection_directional_rhs(id, q_faceL_rs_vf, q_faceR_rs_vf, q_prim_vf, flux_vf, flux_src_vf, rhs_vf)

        integer, intent(in)                                                                 :: id
        real(wp), dimension(idwbuff(1)%beg:,idwbuff(2)%beg:,idwbuff(3)%beg:,1:), intent(in) :: q_faceL_rs_vf, q_faceR_rs_vf
        type(scalar_field), dimension(sys_size), intent(in)                                 :: q_prim_vf
        type(scalar_field), dimension(sys_size), intent(inout)                              :: flux_vf
        type(scalar_field), dimension(sys_size), intent(inout)                              :: flux_src_vf
        type(scalar_field), dimension(sys_size), intent(inout)                              :: rhs_vf
        ! Plain scalar sweep bounds: loop-bound-only scalars are implicitly
        ! firstprivate on device, so no device residency is needed for them
        integer  :: isb1, ise1, isb2, ise2, isb3, ise3
        real(wp) :: rho_L, rho_R
        real(wp) :: u_L, u_R, pres_L, pres_R
        real(wp) :: s_L, s_R, s_star, rho_star
        real(wp) :: F_mass, F_mom, face_vel, pres_flux
        real(wp) :: nrm, vl, vr, alpha_f, a_flux, ar_c, a_c
        real(wp) :: inv_ds, f_m, f_p, divu_c
        real(wp) :: re_l, re_r
        integer  :: ibr
        integer  :: i, j, k, l, q

        if (id == 1) then
            isb1 = -1; ise1 = m; isb2 = 0; ise2 = n; isb3 = 0; ise3 = p
        else if (id == 2) then
            isb1 = -1; ise1 = n; isb2 = 0; ise2 = m; isb3 = 0; ise3 = p
        else
            isb1 = -1; ise1 = p; isb2 = 0; ise2 = n; isb3 = 0; ise3 = m
        end if

        if (id == 1 .and. bodyForces) call s_compute_acceleration(mytime)

        ! Pressure at the start of the stage, used by the p*div(u) source and
        ! the pressure-advection RK blend (captured before any state update)
        if (id == 1) then
            $:GPU_PARALLEL_LOOP(private='[j, k, l]', collapse=3)
            do l = 0, p
                do k = 0, n
                    do j = 0, m
                        pres_stage(j, k, l) = q_prim_vf(eqn_idx%E)%sf(j, k, l)
                    end do
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()
        end if

        #:for NORM_DIR, XYZ, SV, COORDS, JB, JE, KB, KE, LB, LE in &
            [(1, 'x', 'j', '{SI}, k, l', 'isb1', 'ise1', 'isb2', 'ise2', 'isb3', 'ise3'), &
             (2, 'y', 'k', 'j, {SI}, l', 'isb2', 'ise2', 'isb1', 'ise1', 'isb3', 'ise3'), &
             (3, 'z', 'l', 'j, k, {SI}', 'isb3', 'ise3', 'isb2', 'ise2', 'isb1', 'ise1')]
            #:set SF = lambda offs: COORDS.format(SI=SV + offs)
            if (id == ${NORM_DIR}$) then
                ! Face fluxes: left state at the face index, right state at face index + 1
                $:GPU_PARALLEL_LOOP(collapse=3, private='[i, j, k, l, q, rho_L, rho_R, u_L, u_R, pres_L, pres_R, s_L, s_R, &
                                    & s_star, rho_star, F_mass, F_mom, face_vel, pres_flux, nrm, vl, vr, alpha_f, a_flux, ar_c, &
                                    & a_c, re_l, re_r, ibr]')
                do l = ${LB}$, ${LE}$
                    do k = ${KB}$, ${KE}$
                        do j = ${JB}$, ${JE}$
                            rho_L = 0._wp; rho_R = 0._wp

                            $:GPU_LOOP(parallelism='[seq]')
                            do i = 1, num_fluids
                                rho_L = rho_L + q_faceL_rs_vf(${SF('')}$, i)
                                rho_R = rho_R + q_faceR_rs_vf(${SF(' + 1')}$, i)
                            end do

                            u_L = q_faceL_rs_vf(${SF('')}$, eqn_idx%cont%end + ${NORM_DIR}$)
                            u_R = q_faceR_rs_vf(${SF(' + 1')}$, eqn_idx%cont%end + ${NORM_DIR}$)
                            pres_L = q_faceL_rs_vf(${SF('')}$, eqn_idx%E)
                            pres_R = q_faceR_rs_vf(${SF(' + 1')}$, eqn_idx%E)

                            ! Advective-only wave speeds (no sound speed)
                            s_L = min(u_L, u_R)
                            s_R = max(u_L, u_R)
                            s_star = 0.5_wp*(u_L + u_R)

                            ! With advective wave speeds the star densities
                            ! rho_K*(s_K - u_K)/(s_K - s_star) reduce exactly to 0 (the
                            ! near velocity is the extremum, zero numerator) or 2*rho_K
                            ! (the far one is). The division form produces Inf when
                            ! u_L - u_R underflows while its half rounds to zero, so the
                            ! reduced form is used instead
                            rho_star = 0._wp
                            if (s_L >= 0._wp) then
                                ibr = 1
                                F_mass = rho_L*u_L
                                face_vel = u_L
                                pres_flux = pres_L*u_L
                            else if (s_R <= 0._wp) then
                                ibr = 2
                                F_mass = rho_R*u_R
                                face_vel = u_R
                                pres_flux = pres_R*u_R
                            else if (s_star >= 0._wp) then
                                ibr = 3
                                if (u_L > u_R) rho_star = 2._wp*rho_L
                                F_mass = rho_L*u_L + s_L*(rho_star - rho_L)
                                face_vel = s_star
                                pres_flux = pres_L*s_star
                            else
                                ibr = 4
                                if (u_R < u_L) rho_star = 2._wp*rho_R
                                F_mass = rho_R*u_R + s_R*(rho_star - rho_R)
                                face_vel = s_star
                                pres_flux = pres_R*s_star
                            end if

                            ! Momentum fluxes without the p*n term
                            $:GPU_LOOP(parallelism='[seq]')
                            do i = 1, num_dims
                                nrm = 0._wp
                                if (i == ${NORM_DIR}$) nrm = 1._wp
                                vl = q_faceL_rs_vf(${SF('')}$, eqn_idx%cont%end + i)
                                vr = q_faceR_rs_vf(${SF(' + 1')}$, eqn_idx%cont%end + i)
                                if (ibr == 1) then
                                    F_mom = rho_L*vl*u_L
                                else if (ibr == 2) then
                                    F_mom = rho_R*vr*u_R
                                else if (ibr == 3) then
                                    F_mom = rho_L*vl*u_L + s_L*(rho_star*(vl + (s_star - u_L)*nrm) - rho_L*vl)
                                else
                                    F_mom = rho_R*vr*u_R + s_R*(rho_star*(vr + (s_star - u_R)*nrm) - rho_R*vr)
                                end if
                                flux_vf(eqn_idx%cont%end + i)%sf(${SF('')}$) = real(F_mom, stp)
                            end do

                            ! Volume fraction and partial density fluxes, upwinded by the
                            ! contact/face velocity side (star region: by sign of s_star)
                            $:GPU_LOOP(parallelism='[seq]')
                            do i = 1, num_fluids
                                if (ibr == 1 .or. (ibr >= 3 .and. s_star >= 0._wp)) then
                                    alpha_f = q_faceL_rs_vf(${SF('')}$, eqn_idx%adv%beg + i - 1)
                                else
                                    alpha_f = q_faceR_rs_vf(${SF(' + 1')}$, eqn_idx%adv%beg + i - 1)
                                end if
                                a_flux = alpha_f*face_vel
                                flux_vf(eqn_idx%adv%beg + i - 1)%sf(${SF('')}$) = real(a_flux, stp)
                                if (num_fluids > 1) then
                                    ! Partial density flux from the cell upwind of the total mass flux
                                    if (F_mass >= 0._wp) then
                                        ar_c = real(q_prim_vf(i)%sf(${SF('')}$), wp)
                                        a_c = real(q_prim_vf(eqn_idx%adv%beg + i - 1)%sf(${SF('')}$), wp)
                                    else
                                        ar_c = real(q_prim_vf(i)%sf(${SF(' + 1')}$), wp)
                                        a_c = real(q_prim_vf(eqn_idx%adv%beg + i - 1)%sf(${SF(' + 1')}$), wp)
                                    end if
                                    flux_vf(i)%sf(${SF('')}$) = real(ar_c/max(a_c, sgm_eps)*a_flux, stp)
                                else
                                    flux_vf(i)%sf(${SF('')}$) = real(F_mass, stp)
                                end if
                            end do

                            ! Color function advects like the volume fractions
                            if (surface_tension) then
                                if (ibr == 1 .or. (ibr >= 3 .and. s_star >= 0._wp)) then
                                    alpha_f = q_faceL_rs_vf(${SF('')}$, eqn_idx%c)
                                else
                                    alpha_f = q_faceR_rs_vf(${SF(' + 1')}$, eqn_idx%c)
                                end if
                                flux_vf(eqn_idx%c)%sf(${SF('')}$) = real(alpha_f*face_vel, stp)
                            end if

                            ! Face-averaged Reynolds numbers for the viscous source flux
                            ! (harmonic mean, as in the Riemann solvers); the face velocity
                            ! array only feeds the discarded viscous energy work terms
                            if (viscous) then
                                $:GPU_LOOP(parallelism='[seq]')
                                do i = 1, num_vels
                                    vel_src_rsx_vf(${SF('')}$, i) = 0._wp
                                end do
                                $:GPU_LOOP(parallelism='[seq]')
                                do i = 1, 2
                                    re_l = dflt_real; re_r = dflt_real
                                    if (Re_size(i) > 0) then
                                        re_l = 0._wp; re_r = 0._wp
                                    end if
                                    $:GPU_LOOP(parallelism='[seq]')
                                    do q = 1, Re_size(i)
                                        re_l = re_l + q_faceL_rs_vf(${SF('')}$, eqn_idx%adv%beg + Re_idx(i, q) - 1)/Res_gs(i, q)
                                        re_r = re_r + q_faceR_rs_vf(${SF(' + 1')}$, eqn_idx%adv%beg + Re_idx(i, q) - 1)/Res_gs(i, q)
                                    end do
                                    re_l = 1._wp/max(re_l, sgm_eps)
                                    re_r = 1._wp/max(re_r, sgm_eps)
                                    Re_avg_rsx_vf(${SF('')}$, i) = 2._wp/(1._wp/re_l + 1._wp/re_r)
                                end do
                            end if

                            ! Energy flux is unused: total energy is rebuilt from the EOS
                            ! after every pressure solve, so its RHS never enters the state
                            flux_vf(eqn_idx%E)%sf(${SF('')}$) = 0._stp
                            flux_face_vel(${SF('')}$) = real(face_vel, stp)
                            flux_pu(${SF('')}$) = real(pres_flux, stp)
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()

                ! Flux differencing into the RHS, plus div(u) and the p*u flux divergence
                $:GPU_PARALLEL_LOOP(collapse=3, private='[i, j, k, l, inv_ds, f_m, f_p]')
                do l = 0, p
                    do k = 0, n
                        do j = 0, m
                            inv_ds = 1._wp/d${XYZ}$ (${SV}$)
                            $:GPU_LOOP(parallelism='[seq]')
                            do i = 1, sys_size
                                f_m = real(flux_vf(i)%sf(${SF(' - 1')}$), wp)
                                f_p = real(flux_vf(i)%sf(${SF('')}$), wp)
                                #:if NORM_DIR == 1
                                    rhs_vf(i)%sf(j, k, l) = real(inv_ds*(f_m - f_p), stp)
                                #:else
                                    rhs_vf(i)%sf(j, k, l) = rhs_vf(i)%sf(j, k, l) + real(inv_ds*(f_m - f_p), stp)
                                #:endif
                            end do
                            f_m = real(flux_face_vel(${SF(' - 1')}$), wp)
                            f_p = real(flux_face_vel(${SF('')}$), wp)
                            #:if NORM_DIR == 1
                                div_u_face(j, k, l) = real(inv_ds*(f_p - f_m), stp)
                            #:else
                                div_u_face(j, k, l) = div_u_face(j, k, l) + real(inv_ds*(f_p - f_m), stp)
                            #:endif
                            f_m = real(flux_pu(${SF(' - 1')}$), wp)
                            f_p = real(flux_pu(${SF('')}$), wp)
                            #:if NORM_DIR == 1
                                rhs_p_adv(j, k, l) = real(inv_ds*(f_m - f_p), stp)
                            #:else
                                rhs_p_adv(j, k, l) = rhs_p_adv(j, k, l) + real(inv_ds*(f_m - f_p), stp)
                            #:endif
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()

                ! Capillary momentum flux at the faces (explicit CSF stress), accumulated
                ! into flux_src_vf (zeroed by s_initialize_riemann_solver in the caller)
                ! and differenced by s_projection_add_flux_src after any viscous flux is
                ! added. The velocity-work energy terms vanish with the zero cap_vsrc;
                ! the energy RHS is dead under the projection anyway
                if (surface_tension) then
                    cap_isx%beg = ${JB}$; cap_isx%end = ${JE}$
                    cap_isy%beg = ${KB}$; cap_isy%end = ${KE}$
                    cap_isz%beg = ${LB}$; cap_isz%end = ${LE}$
                    $:GPU_UPDATE(device='[cap_isx, cap_isy, cap_isz]')

                    call s_compute_capillary_source_flux(cap_vsrc, flux_src_vf, ${NORM_DIR}$, cap_isx, cap_isy, cap_isz)
                end if
            end if
        #:endfor

        ! div(u) compatibility sources: alpha_k*div(u) for the volume fractions and p*div(u) for the pressure advection equation
        if (id == num_dims) then
            $:GPU_PARALLEL_LOOP(collapse=3, private='[i, j, k, l, divu_c]')
            do l = 0, p
                do k = 0, n
                    do j = 0, m
                        divu_c = real(div_u_face(j, k, l), wp)
                        $:GPU_LOOP(parallelism='[seq]')
                        do i = 1, num_fluids
                            rhs_vf(eqn_idx%adv%beg + i - 1)%sf(j, k, l) = rhs_vf(eqn_idx%adv%beg + i - 1)%sf(j, k, &
                                   & l) + real(real(q_prim_vf(eqn_idx%adv%beg + i - 1)%sf(j, k, l), wp)*divu_c, stp)
                        end do
                        if (surface_tension) then
                            rhs_vf(eqn_idx%c)%sf(j, k, l) = rhs_vf(eqn_idx%c)%sf(j, k, l) + real(real(q_prim_vf(eqn_idx%c)%sf(j, &
                                   & k, l), wp)*divu_c, stp)
                        end if
                        rhs_p_adv(j, k, l) = rhs_p_adv(j, k, l) + real(real(pres_stage(j, k, l), wp)*divu_c, stp)
                    end do
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()

            ! Body-force momentum source, inside the blended RHS so the star momentum
            ! (and div(u*) in the pressure solve) carries it: hydrostatic balance
            if (bodyForces) then
                $:GPU_PARALLEL_LOOP(collapse=3, private='[i, j, k, l, rho_L]')
                do l = 0, p
                    do k = 0, n
                        do j = 0, m
                            rho_L = 0._wp
                            $:GPU_LOOP(parallelism='[seq]')
                            do i = 1, num_fluids
                                rho_L = rho_L + real(q_prim_vf(i)%sf(j, k, l), wp)
                            end do
                            $:GPU_LOOP(parallelism='[seq]')
                            do i = 1, num_dims
                                rhs_vf(eqn_idx%mom%beg + i - 1)%sf(j, k, l) = rhs_vf(eqn_idx%mom%beg + i - 1)%sf(j, k, &
                                       & l) + real(rho_L*accel_bf(i), stp)
                            end do
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()
            end if
        end if

    end subroutine s_projection_directional_rhs

    !> Difference the accumulated source fluxes (viscous stress and capillary CSF) of one direction sweep into the RHS, matching the
    !! convective flux convention
    subroutine s_projection_add_flux_src(id, flux_src_vf, rhs_vf)

        integer, intent(in)                                    :: id
        type(scalar_field), dimension(sys_size), intent(in)    :: flux_src_vf
        type(scalar_field), dimension(sys_size), intent(inout) :: rhs_vf
        real(wp)                                               :: inv_ds, f_m, f_p
        integer                                                :: i, j, k, l

        #:for NORM_DIR, XYZ, SV, COORDS in &
            [(1, 'x', 'j', '{SI}, k, l'), &
             (2, 'y', 'k', 'j, {SI}, l'), &
             (3, 'z', 'l', 'j, k, {SI}')]
            #:set SF = lambda offs: COORDS.format(SI=SV + offs)
            if (id == ${NORM_DIR}$) then
                $:GPU_PARALLEL_LOOP(collapse=3, private='[i, j, k, l, inv_ds, f_m, f_p]')
                do l = 0, p
                    do k = 0, n
                        do j = 0, m
                            inv_ds = 1._wp/d${XYZ}$ (${SV}$)
                            $:GPU_LOOP(parallelism='[seq]')
                            do i = eqn_idx%mom%beg, eqn_idx%E
                                f_m = real(flux_src_vf(i)%sf(${SF(' - 1')}$), wp)
                                f_p = real(flux_src_vf(i)%sf(${SF('')}$), wp)
                                rhs_vf(i)%sf(j, k, l) = rhs_vf(i)%sf(j, k, l) + real(inv_ds*(f_m - f_p), stp)
                            end do
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()
            end if
        #:endfor

    end subroutine s_projection_add_flux_src

    !> Apply the implicit pressure solve and correction to the star state produced by the explicit RK blend: assemble the Helmholtz
    !! RHS from the blended advected pressure and div(u*), solve for the new pressure with Jacobi or red-black Gauss-Seidel, correct
    !! the momentum with the new pressure gradient, and rebuild total energy from the EOS
    impure subroutine s_projection_apply(q_cons_vf, bc_type, pb_in, mv_in, q_T_sf, rkc1, rkc2, rkc3, rkc4, stage, nstage)

        type(scalar_field), dimension(sys_size), intent(inout) :: q_cons_vf
        type(integer_field), dimension(1:num_dims,1:2), intent(in) :: bc_type
        real(stp), optional, dimension(idwbuff(1)%beg:,idwbuff(2)%beg:,idwbuff(3)%beg:,1:,1:), intent(inout) :: pb_in, mv_in
        type(scalar_field), optional, intent(inout) :: q_T_sf
        real(wp), intent(in) :: rkc1, rkc2, rkc3, rkc4
        integer, intent(in) :: stage, nstage
        real(wp) :: p_adv_blend, p0v, rhoc2_sum, blkmod, divs
        real(wp) :: rho_c, rho_nb, u_m, u_p
        real(wp) :: coeff, c_f, offd, diag, p_new, res_loc, res_glb
        real(wp) :: dpds, ke, gamma_mix, pi_inf_mix, qv_mix, mom_sq
        real(wp) :: zv, dv, mu_loc, mu_glb, sigma_ch, rho_ch, rho_prev, alpha_ch, beta_ch
        integer :: i, j, k, l, iter, color

        ! Star-state ghost cells (density and momentum feed the divergence and the Laplacian face densities)

        call s_populate_variables_buffers(bc_type, q_cons_vf, pb_in, mv_in, q_T_sf)

        ! Helmholtz RHS: blended advected pressure minus rho*c^2*dt*div(u*),
        ! with div(u*) as a face-averaged central difference of the star velocity
        $:GPU_PARALLEL_LOOP(collapse=3, private='[i, j, k, l, p_adv_blend, p0v, rhoc2_sum, blkmod, divs, rho_c, rho_nb, u_m, u_p]')
        do l = 0, p
            do k = 0, n
                do j = 0, m
                    if (nstage > 1) then
                        if (stage == 1) pres_step0(j, k, l) = pres_stage(j, k, l)
                        p0v = real(pres_step0(j, k, l), wp)
                    else
                        p0v = 0._wp
                    end if
                    p_adv_blend = (rkc1*real(pres_stage(j, k, l), wp) + rkc2*p0v + rkc3*dt*real(rhs_p_adv(j, k, l), wp))/rkc4

                    ! Wood's mixture sound speed: 1/(rho*c^2) = sum(alpha_k/(gamma_k*(p + pi_inf_k)))
                    rhoc2_sum = 0._wp
                    $:GPU_LOOP(parallelism='[seq]')
                    do i = 1, num_fluids
                        blkmod = ((gammas(i) + 1._wp)*real(pres_stage(j, k, l), wp) + pi_infs(i))/gammas(i)
                        rhoc2_sum = rhoc2_sum + real(q_cons_vf(eqn_idx%adv%beg + i - 1)%sf(j, k, l), wp)/max(blkmod, sgm_eps)
                    end do
                    rhoc2_cell(j, k, l) = real(1._wp/max(rhoc2_sum, sgm_eps), stp)

                    rho_c = 0._wp; rho_nb = 0._wp
                    $:GPU_LOOP(parallelism='[seq]')
                    do i = 1, num_fluids
                        rho_c = rho_c + real(q_cons_vf(i)%sf(j - 1, k, l), wp)
                        rho_nb = rho_nb + real(q_cons_vf(i)%sf(j + 1, k, l), wp)
                    end do
                    u_m = real(q_cons_vf(eqn_idx%mom%beg)%sf(j - 1, k, l), wp)/max(rho_c, sgm_eps)
                    u_p = real(q_cons_vf(eqn_idx%mom%beg)%sf(j + 1, k, l), wp)/max(rho_nb, sgm_eps)
                    divs = 0.5_wp*(u_p - u_m)/dx(j)

                    if (num_dims > 1) then
                        rho_c = 0._wp; rho_nb = 0._wp
                        $:GPU_LOOP(parallelism='[seq]')
                        do i = 1, num_fluids
                            rho_c = rho_c + real(q_cons_vf(i)%sf(j, k - 1, l), wp)
                            rho_nb = rho_nb + real(q_cons_vf(i)%sf(j, k + 1, l), wp)
                        end do
                        u_m = real(q_cons_vf(eqn_idx%mom%beg + 1)%sf(j, k - 1, l), wp)/max(rho_c, sgm_eps)
                        u_p = real(q_cons_vf(eqn_idx%mom%beg + 1)%sf(j, k + 1, l), wp)/max(rho_nb, sgm_eps)
                        divs = divs + 0.5_wp*(u_p - u_m)/dy(k)
                    end if

                    if (num_dims > 2) then
                        rho_c = 0._wp; rho_nb = 0._wp
                        $:GPU_LOOP(parallelism='[seq]')
                        do i = 1, num_fluids
                            rho_c = rho_c + real(q_cons_vf(i)%sf(j, k, l - 1), wp)
                            rho_nb = rho_nb + real(q_cons_vf(i)%sf(j, k, l + 1), wp)
                        end do
                        u_m = real(q_cons_vf(eqn_idx%mom%beg + 2)%sf(j, k, l - 1), wp)/max(rho_c, sgm_eps)
                        u_p = real(q_cons_vf(eqn_idx%mom%beg + 2)%sf(j, k, l + 1), wp)/max(rho_nb, sgm_eps)
                        divs = divs + 0.5_wp*(u_p - u_m)/dz(l)
                    end if

                    helm_rhs(j, k, l) = real(p_adv_blend - real(rhoc2_cell(j, k, l), wp)*dt*divs, stp)
                    pres_proj(j, k, l) = pres_stage(j, k, l)
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

        ! Ghost fill of the initial pressure iterate
        call s_populate_F_igr_buffers(bc_type, pres_proj_sf)

        if (proj_iter_solver /= proj_iter_solver_gauss_seidel) then
            $:GPU_PARALLEL_LOOP(private='[j, k, l]', collapse=3)
            do l = idwbuff(3)%beg, idwbuff(3)%end
                do k = idwbuff(2)%beg, idwbuff(2)%end
                    do j = idwbuff(1)%beg, idwbuff(1)%end
                        pres_proj_old(j, k, l) = pres_proj(j, k, l)
                    end do
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()
        end if

        sigma_ch = 0._wp; rho_ch = 0._wp; mu_glb = 0._wp
        if (proj_iter_solver == proj_iter_solver_chebyshev) then
            ! Gershgorin bound on the Jacobi iteration matrix: its rows are nonnegative
            ! with sum coeff*diag/(1 + coeff*diag) < 1, so the Jacobi-preconditioned
            ! operator spectrum is known in advance to lie in [1 - mu, 1 + mu], which
            ! is what lets Chebyshev run with no global reductions inside the loop
            mu_loc = 0._wp
            $:GPU_PARALLEL_LOOP(collapse=3, private='[i, j, k, l, coeff, c_f, offd, diag, rho_c, rho_nb]', &
                                & reduction='[[mu_loc]]', reductionOp='[max]')
            do l = 0, p
                do k = 0, n
                    do j = 0, m
                        @:PROJECTION_STENCIL(pres_proj)
                        mu_loc = max(mu_loc, coeff*diag/(1._wp + coeff*diag))
                    end do
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()
            call s_mpi_allreduce_max(mu_loc, mu_glb)
            mu_glb = min(max(mu_glb, sgm_eps), 1._wp - sgm_eps)
            sigma_ch = 1._wp/mu_glb
            rho_ch = mu_glb
        end if

        do iter = 1, proj_max_iters
            res_loc = 0._wp

            if (proj_iter_solver == proj_iter_solver_jacobi) then
                $:GPU_PARALLEL_LOOP(collapse=3, private='[i, j, k, l, coeff, c_f, offd, diag, p_new, rho_c, rho_nb]', &
                                    & reduction='[[res_loc]]', reductionOp='[max]')
                do l = 0, p
                    do k = 0, n
                        do j = 0, m
                            @:PROJECTION_STENCIL(pres_proj_old)
                            p_new = (real(helm_rhs(j, k, l), wp) + coeff*offd)/(1._wp + coeff*diag)
                            res_loc = max(res_loc, abs(p_new - real(pres_proj_old(j, k, l), wp)))
                            pres_proj(j, k, l) = real(p_new, stp)
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()

                call s_populate_F_igr_buffers(bc_type, pres_proj_sf)

                $:GPU_PARALLEL_LOOP(private='[j, k, l]', collapse=3)
                do l = idwbuff(3)%beg, idwbuff(3)%end
                    do k = idwbuff(2)%beg, idwbuff(2)%end
                        do j = idwbuff(1)%beg, idwbuff(1)%end
                            pres_proj_old(j, k, l) = pres_proj(j, k, l)
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()
            else if (proj_iter_solver == proj_iter_solver_chebyshev) then
                ! Chebyshev acceleration of the Jacobi iteration (Saad, Alg. 12.1) on
                ! [1 - mu, 1 + mu] (theta = 1, delta = mu): the Jacobi update supplies
                ! z = D^{-1} r, combined through the scalar rho recurrence
                if (iter == 1) then
                    alpha_ch = 1._wp
                    beta_ch = 0._wp
                else
                    rho_prev = rho_ch
                    rho_ch = 1._wp/(2._wp*sigma_ch - rho_prev)
                    alpha_ch = 2._wp*rho_ch/mu_glb
                    beta_ch = rho_ch*rho_prev
                end if

                $:GPU_PARALLEL_LOOP(collapse=3, private='[i, j, k, l, coeff, c_f, offd, diag, p_new, zv, dv, rho_c, rho_nb]', &
                                    & firstprivate='[alpha_ch, beta_ch]', reduction='[[res_loc]]', reductionOp='[max]')
                do l = 0, p
                    do k = 0, n
                        do j = 0, m
                            @:PROJECTION_STENCIL(pres_proj_old)
                            p_new = (real(helm_rhs(j, k, l), wp) + coeff*offd)/(1._wp + coeff*diag)
                            zv = p_new - real(pres_proj_old(j, k, l), wp)
                            dv = beta_ch*real(d_cheb(j, k, l), wp) + alpha_ch*zv
                            d_cheb(j, k, l) = real(dv, stp)
                            res_loc = max(res_loc, abs(dv))
                            pres_proj(j, k, l) = real(real(pres_proj_old(j, k, l), wp) + dv, stp)
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()

                call s_populate_F_igr_buffers(bc_type, pres_proj_sf)

                $:GPU_PARALLEL_LOOP(private='[j, k, l]', collapse=3)
                do l = idwbuff(3)%beg, idwbuff(3)%end
                    do k = idwbuff(2)%beg, idwbuff(2)%end
                        do j = idwbuff(1)%beg, idwbuff(1)%end
                            pres_proj_old(j, k, l) = pres_proj(j, k, l)
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()
            else  ! red-black Gauss-Seidel
                do color = 0, 1
                    $:GPU_PARALLEL_LOOP(collapse=3, private='[i, j, k, l, coeff, c_f, offd, diag, p_new, rho_c, rho_nb]', &
                                        & reduction='[[res_loc]]', reductionOp='[max]')
                    do l = 0, p
                        do k = 0, n
                            do j = 0, m
                                if (mod(j + k + l + rb_offset, 2) == color) then
                                    @:PROJECTION_STENCIL(pres_proj)
                                    p_new = (real(helm_rhs(j, k, l), wp) + coeff*offd)/(1._wp + coeff*diag)
                                    res_loc = max(res_loc, abs(p_new - real(pres_proj(j, k, l), wp)))
                                    pres_proj(j, k, l) = real(p_new, stp)
                                end if
                            end do
                        end do
                    end do
                    $:END_GPU_PARALLEL_LOOP()

                    call s_populate_F_igr_buffers(bc_type, pres_proj_sf)
                end do
            end if

            if (mod(iter, proj_check_iters) == 0 .or. iter == proj_max_iters) then
                call s_mpi_allreduce_max(res_loc, res_glb)
                if (res_glb < proj_tol) exit
            end if
        end do

        ! Momentum correction with the face-averaged new pressure gradient,
        ! then total energy rebuilt from the EOS with the new pressure
        $:GPU_PARALLEL_LOOP(collapse=3, private='[i, j, k, l, dpds, ke, gamma_mix, pi_inf_mix, qv_mix, mom_sq, rho_c]')
        do l = 0, p
            do k = 0, n
                do j = 0, m
                    dpds = 0.5_wp*(real(pres_proj(j + 1, k, l), wp) - real(pres_proj(j - 1, k, l), wp))/dx(j)
                    q_cons_vf(eqn_idx%mom%beg)%sf(j, k, l) = q_cons_vf(eqn_idx%mom%beg)%sf(j, k, l) - real(dt*dpds, stp)
                    if (num_dims > 1) then
                        dpds = 0.5_wp*(real(pres_proj(j, k + 1, l), wp) - real(pres_proj(j, k - 1, l), wp))/dy(k)
                        q_cons_vf(eqn_idx%mom%beg + 1)%sf(j, k, l) = q_cons_vf(eqn_idx%mom%beg + 1)%sf(j, k, l) - real(dt*dpds, stp)
                    end if
                    if (num_dims > 2) then
                        dpds = 0.5_wp*(real(pres_proj(j, k, l + 1), wp) - real(pres_proj(j, k, l - 1), wp))/dz(l)
                        q_cons_vf(eqn_idx%mom%beg + 2)%sf(j, k, l) = q_cons_vf(eqn_idx%mom%beg + 2)%sf(j, k, l) - real(dt*dpds, stp)
                    end if

                    rho_c = 0._wp; gamma_mix = 0._wp; pi_inf_mix = 0._wp; qv_mix = 0._wp
                    $:GPU_LOOP(parallelism='[seq]')
                    do i = 1, num_fluids
                        rho_c = rho_c + real(q_cons_vf(i)%sf(j, k, l), wp)
                        qv_mix = qv_mix + real(q_cons_vf(i)%sf(j, k, l), wp)*qvs(i)
                        gamma_mix = gamma_mix + real(q_cons_vf(eqn_idx%adv%beg + i - 1)%sf(j, k, l), wp)*gammas(i)
                        pi_inf_mix = pi_inf_mix + real(q_cons_vf(eqn_idx%adv%beg + i - 1)%sf(j, k, l), wp)*pi_infs(i)
                    end do

                    mom_sq = 0._wp
                    $:GPU_LOOP(parallelism='[seq]')
                    do i = 1, num_dims
                        mom_sq = mom_sq + real(q_cons_vf(eqn_idx%mom%beg + i - 1)%sf(j, k, l), wp)**2
                    end do
                    ke = 0.5_wp*mom_sq/max(rho_c, sgm_eps)

                    q_cons_vf(eqn_idx%E)%sf(j, k, l) = real(gamma_mix*real(pres_proj(j, k, l), wp) + pi_inf_mix + qv_mix + ke, stp)
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_projection_apply

    !> Finalize the projection module
    impure subroutine s_finalize_projection_module()

        $:GPU_EXIT_DATA(detach='[pres_proj_sf(1)%sf]')

        @:DEALLOCATE(pres_proj)
        if (proj_iter_solver /= proj_iter_solver_gauss_seidel) then
            @:DEALLOCATE(pres_proj_old)
        end if
        if (proj_iter_solver == proj_iter_solver_chebyshev) then
            @:DEALLOCATE(d_cheb)
        end if
        @:DEALLOCATE(pres_stage)
        if (time_stepper /= time_stepper_rk1) then
            @:DEALLOCATE(pres_step0)
        end if
        @:DEALLOCATE(rhs_p_adv, div_u_face, helm_rhs, rhoc2_cell)
        @:DEALLOCATE(flux_face_vel, flux_pu)
        if (surface_tension) then
            @:DEALLOCATE(cap_vsrc)
        end if

    end subroutine s_finalize_projection_module

end module m_projection
