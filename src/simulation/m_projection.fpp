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

#! Same Helmholtz stencil row on multigrid level `lv` at cell (j, k, l), reading
#! neighbor values of `pfield`, the level's cell densities/coefficients, and the
#! level's cell widths. Accumulates offd/diag; coeff is the level coefficient
#:def MG_STENCIL(lv, pfield)
    coeff = real(mg_coeff(${lv}$)%sf(j, k, l), wp)
    offd = 0._wp
    diag = 0._wp
    rho_c = real(mg_rho(${lv}$)%sf(j, k, l), wp)
    rho_nb = real(mg_rho(${lv}$)%sf(j - 1, k, l), wp)
    c_f = 1._wp/(max(0.5_wp*(rho_c + rho_nb), sgm_eps)*0.5_wp*(mg_dx(${lv}$, j - 1) + mg_dx(${lv}$, j))*mg_dx(${lv}$, j))
    offd = offd + c_f*real(${pfield}$ (j - 1, k, l), wp)
    diag = diag + c_f
    rho_nb = real(mg_rho(${lv}$)%sf(j + 1, k, l), wp)
    c_f = 1._wp/(max(0.5_wp*(rho_c + rho_nb), sgm_eps)*0.5_wp*(mg_dx(${lv}$, j) + mg_dx(${lv}$, j + 1))*mg_dx(${lv}$, j))
    offd = offd + c_f*real(${pfield}$ (j + 1, k, l), wp)
    diag = diag + c_f
    if (num_dims > 1) then
        rho_nb = real(mg_rho(${lv}$)%sf(j, k - 1, l), wp)
        c_f = 1._wp/(max(0.5_wp*(rho_c + rho_nb), sgm_eps)*0.5_wp*(mg_dy(${lv}$, k - 1) + mg_dy(${lv}$, k))*mg_dy(${lv}$, k))
        offd = offd + c_f*real(${pfield}$ (j, k - 1, l), wp)
        diag = diag + c_f
        rho_nb = real(mg_rho(${lv}$)%sf(j, k + 1, l), wp)
        c_f = 1._wp/(max(0.5_wp*(rho_c + rho_nb), sgm_eps)*0.5_wp*(mg_dy(${lv}$, k) + mg_dy(${lv}$, k + 1))*mg_dy(${lv}$, k))
        offd = offd + c_f*real(${pfield}$ (j, k + 1, l), wp)
        diag = diag + c_f
    end if
    if (num_dims > 2) then
        rho_nb = real(mg_rho(${lv}$)%sf(j, k, l - 1), wp)
        c_f = 1._wp/(max(0.5_wp*(rho_c + rho_nb), sgm_eps)*0.5_wp*(mg_dz(${lv}$, l - 1) + mg_dz(${lv}$, l))*mg_dz(${lv}$, l))
        offd = offd + c_f*real(${pfield}$ (j, k, l - 1), wp)
        diag = diag + c_f
        rho_nb = real(mg_rho(${lv}$)%sf(j, k, l + 1), wp)
        c_f = 1._wp/(max(0.5_wp*(rho_c + rho_nb), sgm_eps)*0.5_wp*(mg_dz(${lv}$, l) + mg_dz(${lv}$, l + 1))*mg_dz(${lv}$, l))
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
    use m_sim_helpers, only: proj_iters

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

    !> @name Geometric multigrid hierarchy (proj_iter_solver = 4). Each rank coarsens its local grid only: the Helmholtz identity
    !! term makes the coarsest system diagonally dominant once rho*c^2*dt^2/dX^2 < 1, so no rank agglomeration or global coarse
    !! solve is needed. Level 1 is the fine grid; mg_p(1) is unused (the fine solution lives in pres_proj with the standard halo
    !! machinery)
    !> @{
    integer, parameter                            :: mg_max_levels = 12
    integer, parameter                            :: mg_nu_pre = 2  !< pre-smoothing sweeps
    integer, parameter                            :: mg_nu_post = 2  !< post-smoothing sweeps
    integer, parameter                            :: mg_nu_coarse = 8  !< coarsest-level sweeps
    integer                                       :: nlev_mg = 0
    integer                                       :: nlev_eff = 0  !< per-solve V-cycle depth (screening truncates the rest)
    real(wp), dimension(mg_max_levels)            :: mg_dmin  !< per-level minimum cell width (host, set with the dx tables)
    integer, dimension(mg_max_levels)             :: mg_m, mg_n, mg_p_dim  !< local cells - 1 per level
    integer, dimension(mg_max_levels)             :: mg_rboff  !< red-black global parity per level
    type(scalar_field), allocatable, dimension(:) :: mg_p  !< level solution/correction (1-layer ghosts)
    type(scalar_field), allocatable, dimension(:) :: mg_rhs  !< level right-hand side
    type(scalar_field), allocatable, dimension(:) :: mg_res  !< level residual
    type(scalar_field), allocatable, dimension(:) :: mg_rho  !< level cell density (1-layer ghosts)
    type(scalar_field), allocatable, dimension(:) :: mg_coeff  !< level rho*c^2*dt^2
    real(wp), allocatable, dimension(:,:)         :: mg_dx, mg_dy, mg_dz  !< level cell widths (1 ghost each side)
    !> coarse-level halo slabs (wp), one column per direction
    real(wp), allocatable, dimension(:,:) :: mg_sbuf_b, mg_rbuf_b, mg_sbuf_e, mg_rbuf_e
    !> fine-level depth-1 halo slabs, one column per direction
    real(wp), allocatable, dimension(:,:) :: mg_fsbuf_b, mg_frbuf_b, mg_fsbuf_e, mg_frbuf_e
    logical                               :: mg_dx_built = .false.
    $:GPU_DECLARE(create='[mg_p, mg_rhs, mg_res, mg_rho, mg_coeff, mg_dx, mg_dy, mg_dz, mg_rboff]')
    $:GPU_DECLARE(create='[mg_sbuf_b, mg_rbuf_b, mg_sbuf_e, mg_rbuf_e]')
    $:GPU_DECLARE(create='[mg_fsbuf_b, mg_frbuf_b, mg_fsbuf_e, mg_frbuf_e]')
    !> @}

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
        if (proj_iter_solver == proj_iter_solver_multigrid) call s_initialize_projection_mg()

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
                                    ! Clamp so out-of-bounds cell alphas near sharpened interfaces cannot blow up the density ratio
                                    flux_vf(i)%sf(${SF('')}$) = real(max(ar_c, 0._wp)/min(max(a_c, sgm_eps), 1._wp)*a_flux, stp)
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
        real(wp) :: p_adv_blend, p0v, rhoc2_sum, blkmod, divs, a_cl, a_sum
        real(wp) :: rho_c, rho_nb, u_m, u_p
        real(wp) :: coeff, c_f, offd, diag, p_new, res_loc, res_glb
        real(wp) :: dpds, ke, gamma_mix, pi_inf_mix, qv_mix, mom_sq
        real(wp) :: zv, dv, mu_loc, mu_glb, sigma_ch, rho_ch, rho_prev, alpha_ch, beta_ch
        real(wp) :: tol_eff, pmax_loc, pmax_glb
        integer :: i, j, k, l, iter, color

        ! Post-blend clamp and renormalization of the star state (matches the
        ! reference SemiImplicitFV): floors the partial densities and volume
        ! fractions and renormalizes the fractions to sum to one before they
        ! feed the divergence, the Helmholtz coefficients, and the EOS rebuild

        if (proj_normalization) then
            call nvtxStartRange("TIMESTEP-PROJECTION-NORMALIZE")
            $:GPU_PARALLEL_LOOP(collapse=3, private='[i, j, k, l, a_sum]')
            do l = 0, p
                do k = 0, n
                    do j = 0, m
                        a_sum = 0._wp
                        $:GPU_LOOP(parallelism='[seq]')
                        do i = 1, num_fluids
                            q_cons_vf(i)%sf(j, k, l) = max(q_cons_vf(i)%sf(j, k, l), real(sgm_eps, stp))
                            q_cons_vf(eqn_idx%adv%beg + i - 1)%sf(j, k, l) = max(q_cons_vf(eqn_idx%adv%beg + i - 1)%sf(j, k, l), &
                                      & real(sgm_eps, stp))
                            a_sum = a_sum + real(q_cons_vf(eqn_idx%adv%beg + i - 1)%sf(j, k, l), wp)
                        end do
                        $:GPU_LOOP(parallelism='[seq]')
                        do i = 1, num_fluids
                            q_cons_vf(eqn_idx%adv%beg + i - 1)%sf(j, k, l) = real(real(q_cons_vf(eqn_idx%adv%beg + i - 1)%sf(j, &
                                      & k, l), wp)/max(a_sum, sgm_eps), stp)
                        end do
                    end do
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()
            call nvtxEndRange()
        end if

        ! Star-state ghost cells (density and momentum feed the divergence and the Laplacian face densities)
        call nvtxStartRange("TIMESTEP-PROJECTION-COMM")
        call s_populate_variables_buffers(bc_type, q_cons_vf, pb_in, mv_in, q_T_sf)
        call nvtxEndRange

        ! Helmholtz RHS: blended advected pressure minus rho*c^2*dt*div(u*),
        ! with div(u*) as a face-averaged central difference of the star velocity
        call nvtxStartRange("TIMESTEP-PROJECTION-HELM-RHS")
        $:GPU_PARALLEL_LOOP(collapse=3, private='[i, j, k, l, p_adv_blend, p0v, rhoc2_sum, blkmod, divs, a_cl, a_sum, rho_c, &
                            & rho_nb, u_m, u_p]')
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

                    ! Wood's mixture sound speed: 1/(rho*c^2) = sum(alpha_k/(gamma_k*(p + pi_inf_k))).
                    ! Star alphas can leave [0, 1] near sharpened interfaces (THINC,
                    ! mpp_lim), which would collapse the sum and blow the Helmholtz
                    ! coefficient up to 1/sgm_eps, so they are clamped and the sum
                    ! renormalized; a fully degenerate cell gets rhoc2 = 0 (identity row)
                    rhoc2_sum = 0._wp
                    a_sum = 0._wp
                    $:GPU_LOOP(parallelism='[seq]')
                    do i = 1, num_fluids
                        a_cl = min(max(real(q_cons_vf(eqn_idx%adv%beg + i - 1)%sf(j, k, l), wp), 0._wp), 1._wp)
                        blkmod = ((gammas(i) + 1._wp)*real(pres_stage(j, k, l), wp) + pi_infs(i))/gammas(i)
                        rhoc2_sum = rhoc2_sum + a_cl/max(blkmod, sgm_eps)
                        a_sum = a_sum + a_cl
                    end do
                    rhoc2_cell(j, k, l) = real(a_sum/max(rhoc2_sum, sgm_eps), stp)

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
        call nvtxEndRange

        ! Ghost fill of the initial pressure iterate
        call nvtxStartRange("TIMESTEP-PROJECTION-COMM")
        call s_populate_F_igr_buffers(bc_type, pres_proj_sf)
        call nvtxEndRange

        call nvtxStartRange("TIMESTEP-PROJECTION-ITER-SETUP")
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

        if (proj_iter_solver == proj_iter_solver_multigrid) call s_mg_setup_solve(q_cons_vf, bc_type)

        ! Convergence threshold: absolute (proj_tol), or relative to the maximum pressure magnitude when proj_tol_rel is specified
        tol_eff = proj_tol
        if (proj_tol_rel > 0._wp) then
            pmax_loc = 0._wp
            $:GPU_PARALLEL_LOOP(collapse=3, private='[j, k, l]', reduction='[[pmax_loc]]', reductionOp='[max]')
            do l = 0, p
                do k = 0, n
                    do j = 0, m
                        pmax_loc = max(pmax_loc, abs(real(pres_stage(j, k, l), wp)))
                    end do
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()
            call s_mpi_allreduce_max(pmax_loc, pmax_glb)
            tol_eff = proj_tol_rel*max(pmax_glb, sgm_eps)
        end if
        call nvtxEndRange

        do iter = 1, proj_max_iters
            res_loc = 0._wp
            call nvtxStartRange("TIMESTEP-PROJECTION-ITER")
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

                call nvtxStartRange("TIMESTEP-PROJECTION-COMM")
                call s_populate_F_igr_buffers(bc_type, pres_proj_sf)
                call nvtxEndRange

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

                call nvtxStartRange("TIMESTEP-PROJECTION-COMM")
                call s_populate_F_igr_buffers(bc_type, pres_proj_sf)
                call nvtxEndRange

                $:GPU_PARALLEL_LOOP(private='[j, k, l]', collapse=3)
                do l = idwbuff(3)%beg, idwbuff(3)%end
                    do k = idwbuff(2)%beg, idwbuff(2)%end
                        do j = idwbuff(1)%beg, idwbuff(1)%end
                            pres_proj_old(j, k, l) = pres_proj(j, k, l)
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()
            else if (proj_iter_solver == proj_iter_solver_multigrid) then
                ! One V-cycle per outer iteration; convergence measured as the cycle-to-cycle change of the fine iterate
                call s_mg_vcycle(bc_type)

                $:GPU_PARALLEL_LOOP(collapse=3, private='[j, k, l]', reduction='[[res_loc]]', reductionOp='[max]')
                do l = 0, p
                    do k = 0, n
                        do j = 0, m
                            res_loc = max(res_loc, abs(real(pres_proj(j, k, l), wp) - real(pres_proj_old(j, k, l), wp)))
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
                    call nvtxStartRange("TIMESTEP-PROJECTION-COMM")
                    call s_populate_F_igr_buffers(bc_type, pres_proj_sf)
                    call nvtxEndRange
                end do
            end if
            call nvtxEndRange
            ! Multigrid checks every cycle: a V-cycle costs far more than the
            ! reduction, so amortizing the check only overshoots converged solves
            if (proj_iter_solver == proj_iter_solver_multigrid .or. mod(iter, &
                & proj_check_iters) == 0 .or. iter == proj_max_iters) then
                call s_mpi_allreduce_max(res_loc, res_glb)
                if (res_glb < tol_eff) exit
            end if
        end do

        ! Step-line diagnostic: iterations summed over this step's RK stages
        if (stage == 1) then
            proj_iters = min(iter, proj_max_iters)
        else
            proj_iters = proj_iters + min(iter, proj_max_iters)
        end if

        ! Momentum correction with the face-averaged new pressure gradient,
        ! then total energy rebuilt from the EOS with the new pressure
        call nvtxStartRange("TIMESTEP-PROJECTION-CORRECT")
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
        call nvtxEndRange

    end subroutine s_projection_apply

    !> Finalize the projection module
    impure subroutine s_finalize_projection_module()

        integer :: i

        $:GPU_EXIT_DATA(detach='[pres_proj_sf(1)%sf]')

        @:DEALLOCATE(pres_proj)
        if (proj_iter_solver /= proj_iter_solver_gauss_seidel) then
            @:DEALLOCATE(pres_proj_old)
        end if
        if (proj_iter_solver == proj_iter_solver_chebyshev) then
            @:DEALLOCATE(d_cheb)
        end if
        if (proj_iter_solver == proj_iter_solver_multigrid) then
            $:GPU_EXIT_DATA(detach='[mg_p(1)%sf]')
            do i = 2, nlev_mg
                @:DEALLOCATE(mg_p(i)%sf)
            end do
            do i = 1, nlev_mg
                @:DEALLOCATE(mg_rho(i)%sf, mg_rhs(i)%sf, mg_res(i)%sf, mg_coeff(i)%sf)
            end do
            @:DEALLOCATE(mg_p, mg_rhs, mg_res, mg_rho, mg_coeff)
            @:DEALLOCATE(mg_dx, mg_dy, mg_dz)
            @:DEALLOCATE(mg_sbuf_b, mg_rbuf_b, mg_sbuf_e, mg_rbuf_e)
            @:DEALLOCATE(mg_fsbuf_b, mg_frbuf_b, mg_fsbuf_e, mg_frbuf_e)
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

    !> Build the multigrid hierarchy: each rank halves its local grid while every rank's local cell counts stay even (agreed by a
    !! global reduction), so coarse cells remain aligned with the global grid and the red-black parity is well defined on every
    !! level. Level 1 aliases the fine-grid solve arrays
    impure subroutine s_initialize_projection_mg()

        real(wp) :: nl_loc, nl_glb
        integer  :: lv, j, d, mx, cnt
        logical  :: can

        mg_m(1) = m; mg_n(1) = n; mg_p_dim(1) = p
        nl_loc = 1._wp
        do lv = 1, mg_max_levels - 1
            can = mod(mg_m(lv) + 1, 2) == 0 .and. (mg_m(lv) + 1)/2 >= 2
            if (n > 0) can = can .and. mod(mg_n(lv) + 1, 2) == 0 .and. (mg_n(lv) + 1)/2 >= 2
            if (p > 0) can = can .and. mod(mg_p_dim(lv) + 1, 2) == 0 .and. (mg_p_dim(lv) + 1)/2 >= 2
            if (.not. can) exit
            mg_m(lv + 1) = (mg_m(lv) + 1)/2 - 1
            mg_n(lv + 1) = mg_n(lv); if (n > 0) mg_n(lv + 1) = (mg_n(lv) + 1)/2 - 1
            mg_p_dim(lv + 1) = mg_p_dim(lv); if (p > 0) mg_p_dim(lv + 1) = (mg_p_dim(lv) + 1)/2 - 1
            nl_loc = real(lv + 1, wp)
        end do
        call s_mpi_allreduce_min(nl_loc, nl_glb)
        nlev_mg = int(nl_glb)

        ! Global red-black parity per level (offsets stay divisible because every
        ! rank's local counts are even on all levels above the agreed coarsest)
        do lv = 1, nlev_mg
            mg_rboff(lv) = 0
            if (num_procs > 1) then
                cnt = f_global_offset(m_glb, num_procs_x, proc_coords(1))/2**(lv - 1)
                if (num_dims > 1) cnt = cnt + f_global_offset(n_glb, num_procs_y, proc_coords(2))/2**(lv - 1)
                if (num_dims > 2) cnt = cnt + f_global_offset(p_glb, num_procs_z, proc_coords(3))/2**(lv - 1)
                mg_rboff(lv) = mod(cnt, 2)
            end if
        end do
        $:GPU_UPDATE(device='[mg_rboff]')

        @:ALLOCATE(mg_p(1:nlev_mg))
        @:ALLOCATE(mg_rhs(1:nlev_mg))
        @:ALLOCATE(mg_res(1:nlev_mg))
        @:ALLOCATE(mg_rho(1:nlev_mg))
        @:ALLOCATE(mg_coeff(1:nlev_mg))

        ! Level 1 solution is the fine iterate itself (standard halo machinery)
        mg_p(1)%sf => pres_proj
        $:GPU_ENTER_DATA(attach='[mg_p(1)%sf]')

        do lv = 1, nlev_mg
            if (lv > 1) then
                @:ALLOCATE(mg_p(lv)%sf(-1:mg_m(lv) + 1, -1:mg_n(lv) + 1, -1:mg_p_dim(lv) + 1))
                @:ACC_SETUP_SFs(mg_p(lv))
            end if
            @:ALLOCATE(mg_rho(lv)%sf(-1:mg_m(lv) + 1, -1:mg_n(lv) + 1, -1:mg_p_dim(lv) + 1))
            @:ALLOCATE(mg_rhs(lv)%sf(0:mg_m(lv), 0:mg_n(lv), 0:mg_p_dim(lv)))
            @:ALLOCATE(mg_res(lv)%sf(0:mg_m(lv), 0:mg_n(lv), 0:mg_p_dim(lv)))
            @:ALLOCATE(mg_coeff(lv)%sf(0:mg_m(lv), 0:mg_n(lv), 0:mg_p_dim(lv)))
            @:ACC_SETUP_SFs(mg_rho(lv), mg_rhs(lv), mg_res(lv), mg_coeff(lv))
        end do

        ! Level cell widths are built lazily at the first solve: the fine grid's
        ! ghost widths are only populated after the modules initialize
        mx = max(m, max(n, p))
        @:ALLOCATE(mg_dx(1:nlev_mg, -1:mx + 1))
        @:ALLOCATE(mg_dy(1:nlev_mg, -1:mx + 1))
        @:ALLOCATE(mg_dz(1:nlev_mg, -1:mx + 1))
        ! Coarse-level halo slabs (level 2 has the largest faces)
        cnt = max((mg_n(min(2, nlev_mg)) + 1)*(mg_p_dim(min(2, nlev_mg)) + 1), (mg_m(min(2, nlev_mg)) + 1)*(mg_p_dim(min(2, &
                  & nlev_mg)) + 1))
        cnt = max(cnt, (mg_m(min(2, nlev_mg)) + 1)*(mg_n(min(2, nlev_mg)) + 1))
        @:ALLOCATE(mg_sbuf_b(1:max(cnt, 1), 1:num_dims))
        @:ALLOCATE(mg_rbuf_b(1:max(cnt, 1), 1:num_dims))
        @:ALLOCATE(mg_sbuf_e(1:max(cnt, 1), 1:num_dims))
        @:ALLOCATE(mg_rbuf_e(1:max(cnt, 1), 1:num_dims))
        ! Fine-level depth-1 slabs for the overlapped smoother halo
        cnt = max((n + 1)*(p + 1), (m + 1)*(p + 1))
        cnt = max(cnt, (m + 1)*(n + 1))
        @:ALLOCATE(mg_fsbuf_b(1:max(cnt, 1), 1:num_dims))
        @:ALLOCATE(mg_frbuf_b(1:max(cnt, 1), 1:num_dims))
        @:ALLOCATE(mg_fsbuf_e(1:max(cnt, 1), 1:num_dims))
        @:ALLOCATE(mg_frbuf_e(1:max(cnt, 1), 1:num_dims))

    end subroutine s_initialize_projection_mg

    !> Build the per-level cell widths, pairwise-summed from the fine grid (supports grid stretching); ghost widths mirror the edge
    !! value except for single-rank periodic wrap. Rank-boundary coarse ghost widths are mirrored too: on stretched grids this only
    !! perturbs the coarse operators, which affects the convergence rate, never the converged fine-grid solution
    impure subroutine s_mg_build_dx()

        integer :: lv, j

        mg_dx = 1._wp; mg_dy = 1._wp; mg_dz = 1._wp
        mg_dx(1,-1:m + 1) = dx(-1:m + 1)
        if (n > 0) mg_dy(1,-1:n + 1) = dy(-1:n + 1)
        if (p > 0) mg_dz(1,-1:p + 1) = dz(-1:p + 1)
        do lv = 2, nlev_mg
            do j = 0, mg_m(lv)
                mg_dx(lv, j) = mg_dx(lv - 1, 2*j) + mg_dx(lv - 1, 2*j + 1)
            end do
            mg_dx(lv, -1) = mg_dx(lv, 0); mg_dx(lv, mg_m(lv) + 1) = mg_dx(lv, mg_m(lv))
            if (bc_x%beg == BC_PERIODIC) mg_dx(lv, -1) = mg_dx(lv, mg_m(lv))
            if (bc_x%end == BC_PERIODIC) mg_dx(lv, mg_m(lv) + 1) = mg_dx(lv, 0)
            if (n > 0) then
                do j = 0, mg_n(lv)
                    mg_dy(lv, j) = mg_dy(lv - 1, 2*j) + mg_dy(lv - 1, 2*j + 1)
                end do
                mg_dy(lv, -1) = mg_dy(lv, 0); mg_dy(lv, mg_n(lv) + 1) = mg_dy(lv, mg_n(lv))
                if (bc_y%beg == BC_PERIODIC) mg_dy(lv, -1) = mg_dy(lv, mg_n(lv))
                if (bc_y%end == BC_PERIODIC) mg_dy(lv, mg_n(lv) + 1) = mg_dy(lv, 0)
            end if
            if (p > 0) then
                do j = 0, mg_p_dim(lv)
                    mg_dz(lv, j) = mg_dz(lv - 1, 2*j) + mg_dz(lv - 1, 2*j + 1)
                end do
                mg_dz(lv, -1) = mg_dz(lv, 0); mg_dz(lv, mg_p_dim(lv) + 1) = mg_dz(lv, mg_p_dim(lv))
                if (bc_z%beg == BC_PERIODIC) mg_dz(lv, -1) = mg_dz(lv, mg_p_dim(lv))
                if (bc_z%end == BC_PERIODIC) mg_dz(lv, mg_p_dim(lv) + 1) = mg_dz(lv, 0)
            end if
        end do
        $:GPU_UPDATE(device='[mg_dx, mg_dy, mg_dz]')

        do lv = 1, nlev_mg
            mg_dmin(lv) = minval(mg_dx(lv,0:mg_m(lv)))
            if (n > 0) mg_dmin(lv) = min(mg_dmin(lv), minval(mg_dy(lv,0:mg_n(lv))))
            if (p > 0) mg_dmin(lv) = min(mg_dmin(lv), minval(mg_dz(lv,0:mg_p_dim(lv))))
        end do
        mg_dx_built = .true.

    end subroutine s_mg_build_dx

    !> Per-solve multigrid setup: level-1 density, coefficient, and right-hand side from the assembled fine system, then averaged
    !! down the hierarchy (rediscretized coarse operators)
    impure subroutine s_mg_setup_solve(q_cons_vf, bc_type)

        type(scalar_field), dimension(sys_size), intent(in)        :: q_cons_vf
        type(integer_field), dimension(1:num_dims,1:2), intent(in) :: bc_type
        real(wp)                                                   :: sm, cf_s
        integer                                                    :: lv, mml, nnl, ppl
        integer                                                    :: i, j, k, l, jf, kf, lf, nchild
        integer                                                    :: gy, gz

        call nvtxStartRange("TIMESTEP-PROJECTION-MG-SETUP")
        if (.not. mg_dx_built) call s_mg_build_dx()

        ! Ghost planes exist only in active dimensions: q_cons_vf has no
        ! y/z ghosts in 1D/2D, and the coarse stencil never reads them there
        gy = min(1, n); gz = min(1, p)
        $:GPU_PARALLEL_LOOP(collapse=3, private='[i, j, k, l, sm]')
        do l = -gz, p + gz
            do k = -gy, n + gy
                do j = -1, m + 1
                    sm = 0._wp
                    $:GPU_LOOP(parallelism='[seq]')
                    do i = 1, num_fluids
                        sm = sm + real(q_cons_vf(i)%sf(j, k, l), wp)
                    end do
                    mg_rho(1)%sf(j, k, l) = real(sm, stp)
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

        $:GPU_PARALLEL_LOOP(collapse=3, private='[j, k, l]')
        do l = 0, p
            do k = 0, n
                do j = 0, m
                    mg_coeff(1)%sf(j, k, l) = real(real(rhoc2_cell(j, k, l), wp)*dt*dt, stp)
                    mg_rhs(1)%sf(j, k, l) = helm_rhs(j, k, l)
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

        ! Effective depth: coarsen only while the screened operator is still
        ! stiff. On level lv the off-diagonal weight is ~coeff*2*num_dims/dx^2
        ! relative to the identity; once it drops below one the level is
        ! diagonally dominant and the coarsest sweeps finish the job, so
        ! deeper levels add halo latency without helping convergence
        cf_s = 0._wp
        $:GPU_PARALLEL_LOOP(collapse=3, private='[j, k, l]', reduction='[[cf_s]]', reductionOp='[max]')
        do l = 0, p
            do k = 0, n
                do j = 0, m
                    cf_s = max(cf_s, real(mg_coeff(1)%sf(j, k, l), wp))
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()
        call s_mpi_allreduce_max(cf_s, sm)
        nlev_eff = nlev_mg
        do lv = 1, nlev_mg
            if (sm*2._wp*real(num_dims, wp)/max(mg_dmin(lv)**2, sgm_eps) <= 1._wp) then
                nlev_eff = lv
                exit
            end if
        end do

        nchild = 2**num_dims
        do lv = 1, nlev_eff - 1
            mml = mg_m(lv + 1); nnl = mg_n(lv + 1); ppl = mg_p_dim(lv + 1)
            $:GPU_PARALLEL_LOOP(collapse=3, private='[j, k, l, jf, kf, lf, sm, cf_s]', firstprivate='[lv, nchild]')
            do l = 0, ppl
                do k = 0, nnl
                    do j = 0, mml
                        jf = 2*j; kf = k; lf = l
                        if (num_dims > 1) kf = 2*k
                        if (num_dims > 2) lf = 2*l
                        sm = real(mg_rho(lv)%sf(jf, kf, lf), wp) + real(mg_rho(lv)%sf(jf + 1, kf, lf), wp)
                        cf_s = real(mg_coeff(lv)%sf(jf, kf, lf), wp) + real(mg_coeff(lv)%sf(jf + 1, kf, lf), wp)
                        if (num_dims > 1) then
                            sm = sm + real(mg_rho(lv)%sf(jf, kf + 1, lf), wp) + real(mg_rho(lv)%sf(jf + 1, kf + 1, lf), wp)
                            cf_s = cf_s + real(mg_coeff(lv)%sf(jf, kf + 1, lf), wp) + real(mg_coeff(lv)%sf(jf + 1, kf + 1, lf), wp)
                        end if
                        if (num_dims > 2) then
                            sm = sm + real(mg_rho(lv)%sf(jf, kf, lf + 1), wp) + real(mg_rho(lv)%sf(jf + 1, kf, lf + 1), &
                                           & wp) + real(mg_rho(lv)%sf(jf, kf + 1, lf + 1), wp) + real(mg_rho(lv)%sf(jf + 1, &
                                           & kf + 1, lf + 1), wp)
                            cf_s = cf_s + real(mg_coeff(lv)%sf(jf, kf, lf + 1), wp) + real(mg_coeff(lv)%sf(jf + 1, kf, lf + 1), &
                                               & wp) + real(mg_coeff(lv)%sf(jf, kf + 1, lf + 1), &
                                               & wp) + real(mg_coeff(lv)%sf(jf + 1, kf + 1, lf + 1), wp)
                        end if
                        mg_rho(lv + 1)%sf(j, k, l) = real(sm/real(nchild, wp), stp)
                        mg_coeff(lv + 1)%sf(j, k, l) = real(cf_s/real(nchild, wp), stp)
                    end do
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()
            call s_mg_halo_coarse(lv + 1, mg_rho(lv + 1))
        end do

        call nvtxEndRange

    end subroutine s_mg_setup_solve

    !> One-layer halo exchange plus physical BCs for a coarse-level field: MPI neighbors from the domain-level bc codes, periodic
    !! wrap on a single rank, zero-order extrapolation otherwise. Face-varying bc patches reduce to the domain code here, which only
    !! perturbs the coarse operators. 5/7-point stencils read no corner ghosts, so faces alone suffice
    impure subroutine s_mg_halo_coarse(lv_in, f)

        integer, intent(in)               :: lv_in
        type(scalar_field), intent(inout) :: f
        integer                           :: mml, nnl, ppl, cnt, j, k, l, lv
        integer                           :: reqs(12), nreq

        call nvtxStartRange("TIMESTEP-PROJECTION-MG-COMM")
        lv = lv_in
        mml = mg_m(lv); nnl = mg_n(lv); ppl = mg_p_dim(lv)
        nreq = 0

        ! Pack the interior layers of every MPI side of every direction
        #:for DIR, BCV, T1, T1E, T2, T2E, GBEG, IBEG, GEND, IEND in &
            [(1, 'bc_x', 'k', 'nnl', 'l', 'ppl', '(-1, k, l)', '(0, k, l)', '(mml + 1, k, l)', '(mml, k, l)'), &
             (2, 'bc_y', 'j', 'mml', 'l', 'ppl', '(j, -1, l)', '(j, 0, l)', '(j, nnl + 1, l)', '(j, nnl, l)'), &
             (3, 'bc_z', 'j', 'mml', 'k', 'nnl', '(j, k, -1)', '(j, k, 0)', '(j, k, ppl + 1)', '(j, k, ppl)')]
            if (num_dims >= ${DIR}$) then
                if (${BCV}$%beg >= 0) then
                    $:GPU_PARALLEL_LOOP(collapse=2, private='[j, k, l]', firstprivate='[mml, nnl, ppl]')
                    do ${T2}$ = 0, ${T2E}$
                        do ${T1}$ = 0, ${T1E}$
                            mg_sbuf_b(1 + ${T1}$ + ${T2}$*(${T1E}$ + 1), ${DIR}$) = real(f%sf${IBEG}$, wp)
                        end do
                    end do
                    $:END_GPU_PARALLEL_LOOP()
                end if
                if (${BCV}$%end >= 0) then
                    $:GPU_PARALLEL_LOOP(collapse=2, private='[j, k, l]', firstprivate='[mml, nnl, ppl]')
                    do ${T2}$ = 0, ${T2E}$
                        do ${T1}$ = 0, ${T1E}$
                            mg_sbuf_e(1 + ${T1}$ + ${T2}$*(${T1E}$ + 1), ${DIR}$) = real(f%sf${IEND}$, wp)
                        end do
                    end do
                    $:END_GPU_PARALLEL_LOOP()
                end if
            end if
        #:endfor
        if (.not. rdma_mpi) then
            $:GPU_UPDATE(host='[mg_sbuf_b, mg_sbuf_e]')
        end if

        ! Post every direction's nonblocking exchange (even per-direction tag
        ! bases keep messages apart when partners repeat), then wait once
        #:for RDMA in [False, True]
            if (rdma_mpi .eqv. ${'.true.' if RDMA else '.false.'}$) then
                #:if RDMA
                    #:call GPU_HOST_DATA(use_device_addr='[mg_sbuf_b, mg_rbuf_b, mg_sbuf_e, mg_rbuf_e]')
                        #:for DIR, BCV, T1, T1E, T2, T2E, GBEG, IBEG, GEND, IEND in &
            [(1, 'bc_x', 'k', 'nnl', 'l', 'ppl', '(-1, k, l)', '(0, k, l)', '(mml + 1, k, l)', '(mml, k, l)'), &
             (2, 'bc_y', 'j', 'mml', 'l', 'ppl', '(j, -1, l)', '(j, 0, l)', '(j, nnl + 1, l)', '(j, nnl, l)'), &
             (3, 'bc_z', 'j', 'mml', 'k', 'nnl', '(j, k, -1)', '(j, k, 0)', '(j, k, ppl + 1)', '(j, k, ppl)')]
                            if (num_dims >= ${DIR}$) then
                                cnt = (${T1E}$ + 1)*(${T2E}$ + 1)
                                call s_mpi_iexchange_sides_wp(mg_sbuf_b(:,${DIR}$), mg_rbuf_b(:,${DIR}$), mg_sbuf_e(:,${DIR}$), &
                                                              & mg_rbuf_e(:,${DIR}$), cnt, ${BCV}$%beg, ${BCV}$%end, &
                                                              & 2*(${DIR}$ - 1), reqs, nreq)
                            end if
                        #:endfor
                        call s_mpi_wait_requests(reqs, nreq)
                    #:endcall GPU_HOST_DATA
                    $:GPU_WAIT()
                #:else
                    #:for DIR, BCV, T1, T1E, T2, T2E, GBEG, IBEG, GEND, IEND in &
            [(1, 'bc_x', 'k', 'nnl', 'l', 'ppl', '(-1, k, l)', '(0, k, l)', '(mml + 1, k, l)', '(mml, k, l)'), &
             (2, 'bc_y', 'j', 'mml', 'l', 'ppl', '(j, -1, l)', '(j, 0, l)', '(j, nnl + 1, l)', '(j, nnl, l)'), &
             (3, 'bc_z', 'j', 'mml', 'k', 'nnl', '(j, k, -1)', '(j, k, 0)', '(j, k, ppl + 1)', '(j, k, ppl)')]
                        if (num_dims >= ${DIR}$) then
                            cnt = (${T1E}$ + 1)*(${T2E}$ + 1)
                            call s_mpi_iexchange_sides_wp(mg_sbuf_b(:,${DIR}$), mg_rbuf_b(:,${DIR}$), mg_sbuf_e(:,${DIR}$), &
                                                          & mg_rbuf_e(:,${DIR}$), cnt, ${BCV}$%beg, ${BCV}$%end, 2*(${DIR}$ - 1), &
                                                          & reqs, nreq)
                        end if
                    #:endfor
                    call s_mpi_wait_requests(reqs, nreq)
                #:endif
            end if
        #:endfor
        if (.not. rdma_mpi) then
            $:GPU_UPDATE(device='[mg_rbuf_b, mg_rbuf_e]')
        end if

        ! Unpack MPI sides; physical sides wrap (single-rank periodic) or mirror
        #:for DIR, BCV, T1, T1E, T2, T2E, GBEG, IBEG, GEND, IEND in &
            [(1, 'bc_x', 'k', 'nnl', 'l', 'ppl', '(-1, k, l)', '(0, k, l)', '(mml + 1, k, l)', '(mml, k, l)'), &
             (2, 'bc_y', 'j', 'mml', 'l', 'ppl', '(j, -1, l)', '(j, 0, l)', '(j, nnl + 1, l)', '(j, nnl, l)'), &
             (3, 'bc_z', 'j', 'mml', 'k', 'nnl', '(j, k, -1)', '(j, k, 0)', '(j, k, ppl + 1)', '(j, k, ppl)')]
            if (num_dims >= ${DIR}$) then
                if (${BCV}$%beg >= 0) then
                    $:GPU_PARALLEL_LOOP(collapse=2, private='[j, k, l]', firstprivate='[mml, nnl, ppl]')
                    do ${T2}$ = 0, ${T2E}$
                        do ${T1}$ = 0, ${T1E}$
                            f%sf${GBEG}$ = real(mg_rbuf_b(1 + ${T1}$ + ${T2}$*(${T1E}$ + 1), ${DIR}$), stp)
                        end do
                    end do
                    $:END_GPU_PARALLEL_LOOP()
                else if (${BCV}$%beg == BC_PERIODIC) then
                    $:GPU_PARALLEL_LOOP(collapse=2, private='[j, k, l]', firstprivate='[mml, nnl, ppl]')
                    do ${T2}$ = 0, ${T2E}$
                        do ${T1}$ = 0, ${T1E}$
                            f%sf${GBEG}$ = f%sf${IEND}$
                        end do
                    end do
                    $:END_GPU_PARALLEL_LOOP()
                else
                    $:GPU_PARALLEL_LOOP(collapse=2, private='[j, k, l]', firstprivate='[mml, nnl, ppl]')
                    do ${T2}$ = 0, ${T2E}$
                        do ${T1}$ = 0, ${T1E}$
                            f%sf${GBEG}$ = f%sf${IBEG}$
                        end do
                    end do
                    $:END_GPU_PARALLEL_LOOP()
                end if
                if (${BCV}$%end >= 0) then
                    $:GPU_PARALLEL_LOOP(collapse=2, private='[j, k, l]', firstprivate='[mml, nnl, ppl]')
                    do ${T2}$ = 0, ${T2E}$
                        do ${T1}$ = 0, ${T1E}$
                            f%sf${GEND}$ = real(mg_rbuf_e(1 + ${T1}$ + ${T2}$*(${T1E}$ + 1), ${DIR}$), stp)
                        end do
                    end do
                    $:END_GPU_PARALLEL_LOOP()
                else if (${BCV}$%end == BC_PERIODIC) then
                    $:GPU_PARALLEL_LOOP(collapse=2, private='[j, k, l]', firstprivate='[mml, nnl, ppl]')
                    do ${T2}$ = 0, ${T2E}$
                        do ${T1}$ = 0, ${T1E}$
                            f%sf${GEND}$ = f%sf${IBEG}$
                        end do
                    end do
                    $:END_GPU_PARALLEL_LOOP()
                else
                    $:GPU_PARALLEL_LOOP(collapse=2, private='[j, k, l]', firstprivate='[mml, nnl, ppl]')
                    do ${T2}$ = 0, ${T2E}$
                        do ${T1}$ = 0, ${T1E}$
                            f%sf${GEND}$ = f%sf${IEND}$
                        end do
                    end do
                    $:END_GPU_PARALLEL_LOOP()
                end if
            end if
        #:endfor

        call nvtxEndRange

    end subroutine s_mg_halo_coarse

    !> Red-black Gauss-Seidel smoothing sweeps on one level (halo before each color, so freshly restricted or prolonged iterates
    !! enter consistently)
    !> Pack the depth-1 interior face layers of the fine pressure iterate and post every direction's nonblocking exchange. The
    !! caller overlaps host-side completion with device work and then calls s_mg_halo_fine_end
    impure subroutine s_mg_halo_fine_begin(reqs, nreq)

        integer, dimension(:), intent(inout) :: reqs
        integer, intent(inout)               :: nreq
        integer                              :: cnt, j, k, l

        #:for DIR, BCV, T1, T1E, T2, T2E, IBEG, IEND in &
            [(1, 'bc_x', 'k', 'n', 'l', 'p', '(0, k, l)', '(m, k, l)'), &
             (2, 'bc_y', 'j', 'm', 'l', 'p', '(j, 0, l)', '(j, n, l)'), &
             (3, 'bc_z', 'j', 'm', 'k', 'n', '(j, k, 0)', '(j, k, p)')]
            if (num_dims >= ${DIR}$) then
                if (${BCV}$%beg >= 0) then
                    $:GPU_PARALLEL_LOOP(collapse=2, private='[j, k, l]')
                    do ${T2}$ = 0, ${T2E}$
                        do ${T1}$ = 0, ${T1E}$
                            mg_fsbuf_b(1 + ${T1}$ + ${T2}$*(${T1E}$ + 1), ${DIR}$) = real(pres_proj${IBEG}$, wp)
                        end do
                    end do
                    $:END_GPU_PARALLEL_LOOP()
                end if
                if (${BCV}$%end >= 0) then
                    $:GPU_PARALLEL_LOOP(collapse=2, private='[j, k, l]')
                    do ${T2}$ = 0, ${T2E}$
                        do ${T1}$ = 0, ${T1E}$
                            mg_fsbuf_e(1 + ${T1}$ + ${T2}$*(${T1E}$ + 1), ${DIR}$) = real(pres_proj${IEND}$, wp)
                        end do
                    end do
                    $:END_GPU_PARALLEL_LOOP()
                end if
            end if
        #:endfor
        if (.not. rdma_mpi) then
            $:GPU_UPDATE(host='[mg_fsbuf_b, mg_fsbuf_e]')
        end if

        #:for RDMA in [False, True]
            if (rdma_mpi .eqv. ${'.true.' if RDMA else '.false.'}$) then
                #:if RDMA
                    #:call GPU_HOST_DATA(use_device_addr='[mg_fsbuf_b, mg_frbuf_b, mg_fsbuf_e, mg_frbuf_e]')
                        if (num_dims >= 1) then
                            cnt = (n + 1)*(p + 1)
                            call s_mpi_iexchange_sides_wp(mg_fsbuf_b(:,1), mg_frbuf_b(:,1), mg_fsbuf_e(:,1), mg_frbuf_e(:,1), &
                                                          & cnt, bc_x%beg, bc_x%end, 0, reqs, nreq)
                        end if
                        if (num_dims >= 2) then
                            cnt = (m + 1)*(p + 1)
                            call s_mpi_iexchange_sides_wp(mg_fsbuf_b(:,2), mg_frbuf_b(:,2), mg_fsbuf_e(:,2), mg_frbuf_e(:,2), &
                                                          & cnt, bc_y%beg, bc_y%end, 2, reqs, nreq)
                        end if
                        if (num_dims >= 3) then
                            cnt = (m + 1)*(n + 1)
                            call s_mpi_iexchange_sides_wp(mg_fsbuf_b(:,3), mg_frbuf_b(:,3), mg_fsbuf_e(:,3), mg_frbuf_e(:,3), &
                                                          & cnt, bc_z%beg, bc_z%end, 4, reqs, nreq)
                        end if
                    #:endcall GPU_HOST_DATA
                #:else
                    if (num_dims >= 1) then
                        cnt = (n + 1)*(p + 1)
                        call s_mpi_iexchange_sides_wp(mg_fsbuf_b(:,1), mg_frbuf_b(:,1), mg_fsbuf_e(:,1), mg_frbuf_e(:,1), cnt, &
                                                      & bc_x%beg, bc_x%end, 0, reqs, nreq)
                    end if
                    if (num_dims >= 2) then
                        cnt = (m + 1)*(p + 1)
                        call s_mpi_iexchange_sides_wp(mg_fsbuf_b(:,2), mg_frbuf_b(:,2), mg_fsbuf_e(:,2), mg_frbuf_e(:,2), cnt, &
                                                      & bc_y%beg, bc_y%end, 2, reqs, nreq)
                    end if
                    if (num_dims >= 3) then
                        cnt = (m + 1)*(n + 1)
                        call s_mpi_iexchange_sides_wp(mg_fsbuf_b(:,3), mg_frbuf_b(:,3), mg_fsbuf_e(:,3), mg_frbuf_e(:,3), cnt, &
                                                      & bc_z%beg, bc_z%end, 4, reqs, nreq)
                    end if
                #:endif
            end if
        #:endfor

    end subroutine s_mg_halo_fine_begin

    !> Unpack the completed fine-level exchange into the depth-1 ghost layer and apply the physical fills, reading the same
    !! per-point bc_type codes as the standard fine halo (periodic wrap; reflective and extrapolation coincide at depth one)
    impure subroutine s_mg_halo_fine_end(bc_type)

        type(integer_field), dimension(1:num_dims,1:2), intent(in) :: bc_type
        integer                                                    :: j, k, l, bcc

        if (.not. rdma_mpi) then
            $:GPU_UPDATE(device='[mg_frbuf_b, mg_frbuf_e]')
        end if

        #:for DIR, BCV, T1, T1E, T2, T2E, GBEG, IBEG, GEND, IEND, BCB, BCE in &
            [(1, 'bc_x', 'k', 'n', 'l', 'p', '(-1, k, l)', '(0, k, l)', '(m + 1, k, l)', '(m, k, l)', &
              & 'bc_type(1, 1)%sf(0, k, l)', 'bc_type(1, 2)%sf(0, k, l)'), &
             (2, 'bc_y', 'j', 'm', 'l', 'p', '(j, -1, l)', '(j, 0, l)', '(j, n + 1, l)', '(j, n, l)', &
              & 'bc_type(2, 1)%sf(j, 0, l)', 'bc_type(2, 2)%sf(j, 0, l)'), &
             (3, 'bc_z', 'j', 'm', 'k', 'n', '(j, k, -1)', '(j, k, 0)', '(j, k, p + 1)', '(j, k, p)', &
              & 'bc_type(3, 1)%sf(j, k, 0)', 'bc_type(3, 2)%sf(j, k, 0)')]
            if (num_dims >= ${DIR}$) then
                if (${BCV}$%beg >= 0) then
                    $:GPU_PARALLEL_LOOP(collapse=2, private='[j, k, l]')
                    do ${T2}$ = 0, ${T2E}$
                        do ${T1}$ = 0, ${T1E}$
                            pres_proj${GBEG}$ = real(mg_frbuf_b(1 + ${T1}$ + ${T2}$*(${T1E}$ + 1), ${DIR}$), stp)
                        end do
                    end do
                    $:END_GPU_PARALLEL_LOOP()
                else
                    $:GPU_PARALLEL_LOOP(collapse=2, private='[j, k, l, bcc]')
                    do ${T2}$ = 0, ${T2E}$
                        do ${T1}$ = 0, ${T1E}$
                            bcc = int(${BCB}$)
                            if (bcc == BC_PERIODIC) then
                                pres_proj${GBEG}$ = pres_proj${IEND}$
                            else
                                pres_proj${GBEG}$ = pres_proj${IBEG}$
                            end if
                        end do
                    end do
                    $:END_GPU_PARALLEL_LOOP()
                end if
                if (${BCV}$%end >= 0) then
                    $:GPU_PARALLEL_LOOP(collapse=2, private='[j, k, l]')
                    do ${T2}$ = 0, ${T2E}$
                        do ${T1}$ = 0, ${T1E}$
                            pres_proj${GEND}$ = real(mg_frbuf_e(1 + ${T1}$ + ${T2}$*(${T1E}$ + 1), ${DIR}$), stp)
                        end do
                    end do
                    $:END_GPU_PARALLEL_LOOP()
                else
                    $:GPU_PARALLEL_LOOP(collapse=2, private='[j, k, l, bcc]')
                    do ${T2}$ = 0, ${T2E}$
                        do ${T1}$ = 0, ${T1E}$
                            bcc = int(${BCE}$)
                            if (bcc == BC_PERIODIC) then
                                pres_proj${GEND}$ = pres_proj${IBEG}$
                            else
                                pres_proj${GEND}$ = pres_proj${IEND}$
                            end if
                        end do
                    end do
                    $:END_GPU_PARALLEL_LOOP()
                end if
            end if
        #:endfor

    end subroutine s_mg_halo_fine_end

    impure subroutine s_mg_smooth(lv_in, nsweeps, bc_type)

        integer, intent(in)                                        :: lv_in, nsweeps
        type(integer_field), dimension(1:num_dims,1:2), intent(in) :: bc_type
        real(wp)                                                   :: coeff, c_f, offd, diag, p_new, rho_c, rho_nb
        integer                                                    :: lv, mml, nnl, ppl, sweep, color, rboff
        integer                                                    :: i, j, k, l, jj
        integer                                                    :: kib, kie, lib, lie
        integer                                                    :: reqs(12), nreq

        call nvtxStartRange("TIMESTEP-PROJECTION-MG-SMOOTH")
        lv = lv_in
        mml = mg_m(lv); nnl = mg_n(lv); ppl = mg_p_dim(lv)
        rboff = mg_rboff(lv)

        if (lv == 1) then
            ! Fine level: the halo exchange is overlapped with the interior
            ! update. Pack + post, smooth the cells one layer in from every
            ! face while the host completes MPI, then unpack and finish the
            ! boundary shells. Same values in the same order as the plain
            ! sweep, so results are unchanged
            kib = min(1, n); kie = n - min(1, n)
            lib = min(1, p); lie = p - min(1, p)
            do sweep = 1, nsweeps
                do color = 0, 1
                    if (color == 0 .or. .not. proj_mg_single_halo) then
                        nreq = 0
                        call nvtxStartRange("TIMESTEP-PROJECTION-MG-COMM-POST")
                        call s_mg_halo_fine_begin(reqs, nreq)
                        call nvtxEndRange
                        $:GPU_PARALLEL_LOOP(collapse=3, private='[i, j, k, l, coeff, c_f, offd, diag, p_new, rho_c, rho_nb]', &
                                            & firstprivate='[lv, rboff, color, kib, kie, lib, lie]', extraAccArgs='async(1)', extraOmpArgs='nowait')
                        do l = lib, lie
                            do k = kib, kie
                                do j = 1, m - 1
                                    if (mod(j + k + l + rboff, 2) == color) then
                                        @:MG_STENCIL(lv, mg_p(lv)%sf)
                                        p_new = (real(mg_rhs(lv)%sf(j, k, l), wp) + coeff*offd)/(1._wp + coeff*diag)
                                        mg_p(lv)%sf(j, k, l) = real(p_new, stp)
                                    end if
                                end do
                            end do
                        end do
                        $:END_GPU_PARALLEL_LOOP()
                        call nvtxStartRange("TIMESTEP-PROJECTION-MG-COMM-WAIT")
                        call s_mpi_wait_requests(reqs, nreq)
                        $:GPU_WAIT()
                        call nvtxEndRange
                        call nvtxStartRange("TIMESTEP-PROJECTION-MG-COMM-UNPACK")
                        call s_mg_halo_fine_end(bc_type)
                        call nvtxEndRange
                        ! x-faces: j = 0 and j = m, full k/l extent
                        $:GPU_PARALLEL_LOOP(collapse=3, private='[i, j, k, l, jj, coeff, c_f, offd, diag, p_new, rho_c, rho_nb]', &
                                            & firstprivate='[lv, rboff, color]')
                        do l = 0, p
                            do k = 0, n
                                do jj = 0, 1
                                    j = jj*m
                                    if (mod(j + k + l + rboff, 2) == color) then
                                        @:MG_STENCIL(lv, mg_p(lv)%sf)
                                        p_new = (real(mg_rhs(lv)%sf(j, k, l), wp) + coeff*offd)/(1._wp + coeff*diag)
                                        mg_p(lv)%sf(j, k, l) = real(p_new, stp)
                                    end if
                                end do
                            end do
                        end do
                        $:END_GPU_PARALLEL_LOOP()
                        if (n > 0) then
                            ! y-faces: k = 0 and k = n, interior x extent
                            $:GPU_PARALLEL_LOOP(collapse=3, private='[i, j, k, l, jj, coeff, c_f, offd, diag, p_new, rho_c, &
                                                & rho_nb]', firstprivate='[lv, rboff, color]')
                            do l = 0, p
                                do jj = 0, 1
                                    do j = 1, m - 1
                                        k = jj*n
                                        if (mod(j + k + l + rboff, 2) == color) then
                                            @:MG_STENCIL(lv, mg_p(lv)%sf)
                                            p_new = (real(mg_rhs(lv)%sf(j, k, l), wp) + coeff*offd)/(1._wp + coeff*diag)
                                            mg_p(lv)%sf(j, k, l) = real(p_new, stp)
                                        end if
                                    end do
                                end do
                            end do
                            $:END_GPU_PARALLEL_LOOP()
                        end if
                        if (p > 0) then
                            ! z-faces: l = 0 and l = p, interior x/y extent
                            $:GPU_PARALLEL_LOOP(collapse=3, private='[i, j, k, l, jj, coeff, c_f, offd, diag, p_new, rho_c, &
                                                & rho_nb]', firstprivate='[lv, rboff, color]')
                            do jj = 0, 1
                                do k = 1, n - 1
                                    do j = 1, m - 1
                                        l = jj*p
                                        if (mod(j + k + l + rboff, 2) == color) then
                                            @:MG_STENCIL(lv, mg_p(lv)%sf)
                                            p_new = (real(mg_rhs(lv)%sf(j, k, l), wp) + coeff*offd)/(1._wp + coeff*diag)
                                            mg_p(lv)%sf(j, k, l) = real(p_new, stp)
                                        end if
                                    end do
                                end do
                            end do
                            $:END_GPU_PARALLEL_LOOP()
                        end if
                    else
                        ! proj_mg_single_halo second color: no exchange, full sweep
                        $:GPU_PARALLEL_LOOP(collapse=3, private='[i, j, k, l, coeff, c_f, offd, diag, p_new, rho_c, rho_nb]', &
                                            & firstprivate='[lv, rboff, color]')
                        do l = 0, p
                            do k = 0, n
                                do j = 0, m
                                    if (mod(j + k + l + rboff, 2) == color) then
                                        @:MG_STENCIL(lv, mg_p(lv)%sf)
                                        p_new = (real(mg_rhs(lv)%sf(j, k, l), wp) + coeff*offd)/(1._wp + coeff*diag)
                                        mg_p(lv)%sf(j, k, l) = real(p_new, stp)
                                    end if
                                end do
                            end do
                        end do
                        $:END_GPU_PARALLEL_LOOP()
                    end if
                end do
            end do
        else
            do sweep = 1, nsweeps
                do color = 0, 1
                    if (color == 0 .or. .not. proj_mg_single_halo) then
                        call s_mg_halo_coarse(lv, mg_p(lv))
                    end if
                    $:GPU_PARALLEL_LOOP(collapse=3, private='[i, j, k, l, coeff, c_f, offd, diag, p_new, rho_c, rho_nb]', &
                                        & firstprivate='[lv, mml, nnl, ppl, rboff, color]')
                    do l = 0, ppl
                        do k = 0, nnl
                            do j = 0, mml
                                if (mod(j + k + l + rboff, 2) == color) then
                                    @:MG_STENCIL(lv, mg_p(lv)%sf)
                                    p_new = (real(mg_rhs(lv)%sf(j, k, l), wp) + coeff*offd)/(1._wp + coeff*diag)
                                    mg_p(lv)%sf(j, k, l) = real(p_new, stp)
                                end if
                            end do
                        end do
                    end do
                    $:END_GPU_PARALLEL_LOOP()
                end do
            end do
        end if

        call nvtxEndRange

    end subroutine s_mg_smooth

    !> Residual r = rhs - A p on one level (with a fresh halo on p)
    impure subroutine s_mg_residual(lv_in, bc_type)

        integer, intent(in)                                        :: lv_in
        type(integer_field), dimension(1:num_dims,1:2), intent(in) :: bc_type
        real(wp)                                                   :: coeff, c_f, offd, diag, rho_c, rho_nb
        integer                                                    :: lv, mml, nnl, ppl
        integer                                                    :: i, j, k, l
        integer                                                    :: reqs(12), nreq

        call nvtxStartRange("TIMESTEP-PROJECTION-MG-RESIDUAL")
        lv = lv_in
        mml = mg_m(lv); nnl = mg_n(lv); ppl = mg_p_dim(lv)

        if (lv == 1) then
            ! The 7-point residual stencil reads one ghost layer, so exchange
            ! depth-1 slabs like the smoother rather than the buff_size-deep
            ! generic fine halo
            nreq = 0
            call s_mg_halo_fine_begin(reqs, nreq)
            call s_mpi_wait_requests(reqs, nreq)
            call s_mg_halo_fine_end(bc_type)
        else
            call s_mg_halo_coarse(lv, mg_p(lv))
        end if

        $:GPU_PARALLEL_LOOP(collapse=3, private='[i, j, k, l, coeff, c_f, offd, diag, rho_c, rho_nb]', firstprivate='[lv, mml, &
                            & nnl, ppl]')
        do l = 0, ppl
            do k = 0, nnl
                do j = 0, mml
                    @:MG_STENCIL(lv, mg_p(lv)%sf)
                    mg_res(lv)%sf(j, k, l) = real(real(mg_rhs(lv)%sf(j, k, l), wp) - (1._wp + coeff*diag)*real(mg_p(lv)%sf(j, k, &
                           & l), wp) + coeff*offd, stp)
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

        call nvtxEndRange

    end subroutine s_mg_residual

    !> Full-weighting restriction of the level residual into the next level's right-hand side; the coarse correction starts from
    !! zero
    impure subroutine s_mg_restrict(lv_in)

        integer, intent(in) :: lv_in
        real(wp)            :: sm
        integer             :: lv, mml, nnl, ppl, nchild
        integer             :: j, k, l, jf, kf, lf

        call nvtxStartRange("TIMESTEP-PROJECTION-MG-RESTRICT")
        lv = lv_in
        mml = mg_m(lv + 1); nnl = mg_n(lv + 1); ppl = mg_p_dim(lv + 1)
        nchild = 2**num_dims

        $:GPU_PARALLEL_LOOP(collapse=3, private='[j, k, l, jf, kf, lf, sm]', firstprivate='[lv, mml, nnl, ppl, nchild]')
        do l = -1, ppl + 1
            do k = -1, nnl + 1
                do j = -1, mml + 1
                    mg_p(lv + 1)%sf(j, k, l) = 0._stp
                    if (j >= 0 .and. j <= mml .and. k >= 0 .and. k <= nnl .and. l >= 0 .and. l <= ppl) then
                        jf = 2*j; kf = k; lf = l
                        if (num_dims > 1) kf = 2*k
                        if (num_dims > 2) lf = 2*l
                        sm = real(mg_res(lv)%sf(jf, kf, lf), wp) + real(mg_res(lv)%sf(jf + 1, kf, lf), wp)
                        if (num_dims > 1) then
                            sm = sm + real(mg_res(lv)%sf(jf, kf + 1, lf), wp) + real(mg_res(lv)%sf(jf + 1, kf + 1, lf), wp)
                        end if
                        if (num_dims > 2) then
                            sm = sm + real(mg_res(lv)%sf(jf, kf, lf + 1), wp) + real(mg_res(lv)%sf(jf + 1, kf, lf + 1), &
                                           & wp) + real(mg_res(lv)%sf(jf, kf + 1, lf + 1), wp) + real(mg_res(lv)%sf(jf + 1, &
                                           & kf + 1, lf + 1), wp)
                        end if
                        mg_rhs(lv + 1)%sf(j, k, l) = real(sm/real(nchild, wp), stp)
                    end if
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

        call nvtxEndRange

    end subroutine s_mg_restrict

    !> Piecewise-constant prolongation: add each coarse-cell correction to its child cells on the finer level
    impure subroutine s_mg_prolong(lv_in)

        integer, intent(in) :: lv_in
        integer             :: lv, mml, nnl, ppl
        integer             :: j, k, l, jc, kc, lc

        call nvtxStartRange("TIMESTEP-PROJECTION-MG-PROLONG")
        lv = lv_in
        mml = mg_m(lv); nnl = mg_n(lv); ppl = mg_p_dim(lv)

        $:GPU_PARALLEL_LOOP(collapse=3, private='[j, k, l, jc, kc, lc]', firstprivate='[lv, mml, nnl, ppl]')
        do l = 0, ppl
            do k = 0, nnl
                do j = 0, mml
                    jc = j/2; kc = k; lc = l
                    if (num_dims > 1) kc = k/2
                    if (num_dims > 2) lc = l/2
                    mg_p(lv)%sf(j, k, l) = real(real(mg_p(lv)%sf(j, k, l), wp) + real(mg_p(lv + 1)%sf(jc, kc, lc), wp), stp)
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

        call nvtxEndRange

    end subroutine s_mg_prolong

    !> One multigrid V-cycle on the assembled hierarchy
    impure subroutine s_mg_vcycle(bc_type)

        type(integer_field), dimension(1:num_dims,1:2), intent(in) :: bc_type
        integer                                                    :: lv

        do lv = 1, nlev_eff - 1
            call s_mg_smooth(lv, mg_nu_pre, bc_type)
            call s_mg_residual(lv, bc_type)
            call s_mg_restrict(lv)
        end do

        call s_mg_smooth(nlev_eff, mg_nu_coarse, bc_type)

        do lv = nlev_eff - 1, 1, -1
            call s_mg_prolong(lv)
            call s_mg_smooth(lv, mg_nu_post, bc_type)
        end do

    end subroutine s_mg_vcycle

end module m_projection
