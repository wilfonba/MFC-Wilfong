#:include 'macros.fpp'

!> @brief This module is used to compute source terms for hypoelastic model
module m_surface_tension

    ! Dependencies =============================================================

    use m_derived_types        !< Definitions of the derived types

    use m_global_parameters    !< Definitions of the global parameters

    use m_mpi_proxy            !< Message passing interface (MPI) module proxy

    use m_variables_conversion

    use m_weno

    use m_helper
    ! ==========================================================================

    implicit none

    private; public :: s_initialize_surface_tension_module, &
        s_compute_capilary_source_flux, &
        s_get_capilary, &
        s_finalize_surface_tension_module

    type(vector_field) :: c_divs
    !$acc declare create(c_divs)

    real(kind(0d0)), allocatable, dimension(:,:,:,:) :: gL_x, gR_x, gL_y, gR_y, gL_z, gR_z
    !$acc declare create(gL_x, gR_x, gL_y, gR_y, gL_z, gR_z)

    real(kind(0d0)), allocatable, dimension(:,:,:,:) :: cL_x, cR_x, cL_y, cR_y, cL_z, cR_z
    !$acc declare create(cL_x, cR_x, cL_y, cR_y, cL_z, cR_z)

    type(int_bounds_info) :: ix, iy, iz, is1, is2, is3, iv
    !$acc declare create(ix, iy, iz, is1, is2, is3, iv)

    real(kind(0d0)), allocatable, dimension(:) :: flux_src_H, flux_src_L
    !$acc declare create(flux_src_H, flux_src_L)

    integer :: j, k, l, i  

contains

    subroutine s_initialize_surface_tension_module()

        ! Configuring Coordinate Direction Indexes =========================
        ix%beg = -buff_size; iy%beg = 0; iz%beg = 0

        if (n > 0) iy%beg = -buff_size; if (p > 0) iz%beg = -buff_size

        ix%end = m - ix%beg; iy%end = n - iy%beg; iz%end = p - iz%beg
        ! ==================================================================

        @:ALLOCATE(c_divs%vf(1:num_dims + 1))

        do j = 1,num_dims + 1
            @:ALLOCATE(c_divs%vf(j)%sf(ix%beg:ix%end, iy%beg:iy%end, iz%beg:iz%end))
        end do

        @:ALLOCATE(gL_x(ix%beg:ix%end, iy%beg:iy%end, iz%beg:iz%end, num_dims + 1))
        @:ALLOCATE(gR_x(ix%beg:ix%end, iy%beg:iy%end, iz%beg:iz%end, num_dims + 1))
        @:ALLOCATE(cL_x(ix%beg:ix%end, iy%beg:iy%end, iz%beg:iz%end, c_idx:c_idx))
        @:ALLOCATE(cR_x(ix%beg:ix%end, iy%beg:iy%end, iz%beg:iz%end, c_idx:c_idx))

        @:ALLOCATE(gL_y(iy%beg:iy%end, ix%beg:ix%end, iz%beg:iz%end, num_dims + 1))
        @:ALLOCATE(gR_y(iy%beg:iy%end, ix%beg:ix%end, iz%beg:iz%end, num_dims + 1))
        @:ALLOCATE(cL_y(iy%beg:iy%end, ix%beg:ix%end, iz%beg:iz%end, c_idx:c_idx))
        @:ALLOCATE(cR_y(iy%beg:iy%end, ix%beg:ix%end, iz%beg:iz%end, c_idx:c_idx))
        
        if (p > 0) then
            @:ALLOCATE(gL_z(iz%beg:iz%end, ix%beg:ix%end, iy%beg:iy%end, num_dims + 1))
            @:ALLOCATE(gR_z(iz%beg:iz%end, ix%beg:ix%end, iy%beg:iy%end, num_dims + 1))
            @:ALLOCATE(cL_z(iz%beg:iz%end, ix%beg:ix%end, iy%beg:iy%end, c_idx:c_idx))
            @:ALLOCATE(cR_z(iz%beg:iz%end, ix%beg:ix%end, iy%beg:iy%end, c_idx:c_idx))
        end if

        @:ALLOCATE(flux_src_H(momxb:E_idx), flux_src_L(momxb:E_idx))

    end subroutine s_initialize_surface_tension_module

    subroutine s_compute_capilary_source_flux(q_prim_vf, &
                                              vSrc_rsx_vf, vSrc_rsy_vf, vSrc_rsz_vf, &
                                              flux_src_vf, &
                                              id, isx, isy, isz)

        type(int_bounds_info) :: isx, isy, isz
        type(scalar_field), dimension(sys_size) :: q_prim_vf
        real(kind(0d0)), dimension(-1:m,0:n,0:p,1:num_dims) :: vSrc_rsx_vf
        real(kind(0d0)), dimension(-1:n,0:m,0:p,1:num_dims) :: vSrc_rsy_vf
        real(kind(0d0)), dimension(-1:p,0:n,0:m,1:num_dims) :: vSrc_rsz_vf
        type(scalar_field), &
            dimension(sys_size), &
            intent(INOUT) :: flux_src_vf
        integer :: id

        real(kind(0d0)), dimension(num_dims, num_dims) :: Omega
        real(kind(0d0)) :: w1L, w1R, w2L, w2R, w3L, w3R, w1, w2, w3
        real(kind(0d0)) :: normWL, normWR, normW
        real(kind(0d0)) :: phi

        if (id == 1) then
            !$acc parallel loop collapse(3) gang vector default(present) private(Omega, &
            !$acc w1L, w2L, w3L, w1R, w2R, w3R, w1, w2, w3, normWL, normWR, normW, phi)
            do l = isz%beg, isz%end
                do k = isy%beg, isy%end
                    do j = isx%beg, isx%end

                        ! High order flux ======================================
                        w1L = gL_x(j, k, l, 1)
                        w2L = gL_x(j, k, l, 2)
                        w3L = 0d0
                        if (p > 0) w3L = gL_x(j, k, l, 3)

                        w1R = gR_x(j + 1, k, l, 1)
                        w2R = gR_x(j + 1, k, l, 2)
                        w3R = 0d0
                        if (p > 0) w3R = gR_x(j + 1, k, l, 3)

                        normWL = gL_x(j, k, l, num_dims + 1)
                        normWR = gR_x(j + 1, k, l, num_dims + 1)

                        w1 = (w1L + w1R)/2d0
                        w2 = (w2L + w2R)/2d0
                        w3 = (w3L + w3R)/2d0
                        normW = (normWL + normWR)/2d0

                        call s_compute_capilary_stress_tensor(w1, w2, w3, omega, normW)

                        do i = 1, num_dims

                            flux_src_H(momxb + i - 1) = flux_src_H(momxb + i - 1) + Omega(1,i)

                            flux_src_H(E_idx) = flux_src_H(E_idx) + &
                                Omega(1,i)*vSrc_rsx_vf(j,k,l,i)

                        end do

                        flux_src_H(E_idx) = flux_src_H(E_idx) + &
                            sigma*c_divs%vf(num_dims + 1)%sf(j,k,l)*vSrc_rsx_vf(j, k, l, 1)

                        call s_compute_flux_lim(j, k, l, phi, q_prim_vf, flux_src_vf, id)

                        ! Low order flux =======================================
                        w1L = c_divs%vf(1)%sf(j, k, l)
                        w2L = c_divs%vf(2)%sf(j, k, l)
                        w3L = 0d0
                        if (p > 0) w3L = c_divs%vf(3)%sf(j, k, l)

                        w1R = c_divs%vf(1)%sf(j + 1, k, l)
                        w2R = c_divs%vf(2)%sf(j + 1, k, l)
                        w3R = 0d0
                        if (p > 0) w3R = c_divs%vf(3)%sf(j + 1, k, l)

                        normWL = c_divs%vf(num_dims + 1)%sf(j, k, l)
                        normWR = c_divs%vf(num_dims + 1)%sf(j + 1, k, l)

                        w1 = (w1L + w1R)/2d0
                        w2 = (w2L + w2R)/2d0
                        w3 = (w3L + w3R)/2d0
                        normW = (normWL + normWR)/2d0

                        call s_compute_capilary_stress_tensor(w1, w2, w3, omega, normW)

                        do i = 1, num_dims

                            flux_src_L(momxb + i - 1) = flux_src_L(momxb + i - 1) + Omega(1,i)

                            flux_src_L(E_idx) = flux_src_L(E_idx) + &
                                Omega(1,i)*vSrc_rsx_vf(j,k,l,i)

                        end do

                        flux_src_L(E_idx) = flux_src_L(E_idx) + &
                            sigma*c_divs%vf(num_dims + 1)%sf(j,k,l)*vSrc_rsx_vf(j, k, l, 1)

                        ! Computing final flux =================================
                        call s_compute_flux_lim(j, k, l, phi, q_prim_vf, flux_src_vf, id)

                        do i = 1, num_dims
                            flux_src_vf(momxb + i - 1)%sf(j, k, l) = & 
                                flux_src_vf(momxb + i - 1)%sf(j, k, l) + &
                                flux_src_L(momxb + i - 1) + phi*(flux_src_H(momxb + i - 1) - &
                                flux_src_L(momxb + i - 1))
                        end do

                        flux_src_vf(E_idx)%sf(j, k, l) = flux_src_vf(E_idx)%sf(j, k, l) + &
                            flux_src_L(E_idx) + phi*(flux_src_H(E_idx) - flux_src_L(E_idx))

                    end do
                end do
            end do

        elseif (id == 2) then

            !$acc parallel loop collapse(3) gang vector default(present) private(Omega, &
            !$acc w1L, w2L, w3L, w1R, w2R, w3R, w1, w2, w3, normWL, normWR, normW, phi)
            do l = isz%beg, isz%end
                do k = isy%beg, isy%end
                    do j = isx%beg, isx%end

                        ! High order flux ======================================
                        w1L = gL_y(k, j, l, 1)
                        w2L = gL_y(k, j, l, 2)
                        w3L = 0d0
                        if (p > 0) w3L = gL_y(k, j, l, 3)

                        w1R = gR_y(k + 1, j, l, 1)
                        w2R = gR_y(k + 1, j, l, 2)
                        w3R = 0d0
                        if (p > 0) w3R = gR_y(k + 1, j, l, 3)

                        normWL = gL_y(k, j, l, num_dims + 1)
                        normWR = gR_y(k + 1, j, l, num_dims + 1)

                        w1 = (w1L + w1R)/2d0
                        w2 = (w2L + w2R)/2d0
                        w3 = (w3L + w3R)/2d0
                        normW = (normWL + normWR)/2d0

                        call s_compute_capilary_stress_tensor(w1, w2, w3, omega, normW)

                        do i = 1, num_dims

                            flux_src_H(momxb + i - 1) = flux_src_H(momxb + i - 1) + Omega(2,i)

                            flux_src_H(E_idx) = flux_src_H(E_idx) + &
                                Omega(2,i)*vSrc_rsy_vf(k, j, l, i)
                                
                        end do

                        flux_src_H(E_idx) = flux_src_H(E_idx) + &
                            sigma*c_divs%vf(num_dims + 1)%sf(j,k,l)*vSrc_rsy_vf(k, j, l, 2)

                        ! Low order flux =======================================
                        w1L = c_divs%vf(1)%sf(j, k, l)
                        w2L = c_divs%vf(2)%sf(j, k, l)
                        w3L = 0d0
                        if (p > 0) w3L = c_divs%vf(3)%sf(j, k, l)

                        w1R = c_divs%vf(1)%sf(j, k + 1, l)
                        w2R = c_divs%vf(2)%sf(j, k + 1, l)
                        w3R = 0d0
                        if (p > 0) w3R = c_divs%vf(3)%sf(j, k + 1, l)

                        normWL = c_divs%vf(num_dims + 1)%sf(j, k, l)
                        normWR = c_divs%vf(num_dims + 1)%sf(j, k + 1, l)

                        w1 = (w1L + w1R)/2d0
                        w2 = (w2L + w2R)/2d0
                        w3 = (w3L + w3R)/2d0
                        normW = (normWL + normWR)/2d0

                        call s_compute_capilary_stress_tensor(w1, w2, w3, omega, normW)

                        do i = 1, num_dims

                            flux_src_L(momxb + i - 1) = flux_src_L(momxb + i - 1) + Omega(2,i)

                            flux_src_L(E_idx) = flux_src_L(E_idx) + &
                                Omega(2,i)*vSrc_rsy_vf(k, j, l, i)
                                
                        end do

                        flux_src_L(E_idx) = flux_src_L(E_idx) + &
                            sigma*c_divs%vf(num_dims + 1)%sf(j,k,l)*vSrc_rsy_vf(k, j, l, 2)

                        ! Computing final flux =================================
                        call s_compute_flux_lim(j, k, l, phi, q_prim_vf, flux_src_vf, id)

                        do i = 1, num_dims
                            flux_src_vf(momxb + i - 1)%sf(j, k, l) = & 
                                flux_src_vf(momxb + i - 1)%sf(j, k, l) + &
                                flux_src_L(momxb + i - 1) + phi*(flux_src_H(momxb + i - 1) - &
                                flux_src_L(momxb + i - 1))
                        end do

                        flux_src_vf(E_idx)%sf(j, k, l) = flux_src_vf(E_idx)%sf(j, k, l) + &
                            flux_src_L(E_idx) + phi*(flux_src_H(E_idx) - flux_src_L(E_idx))

                    end do
                end do
            end do

        elseif (id == 3) then

            !$acc parallel loop collapse(3) gang vector default(present) private(Omega, &
            !$acc w1L, w2L, w3L, w1R, w2R, w3R, w1, w2, w3, normWL, normWR, normW)
            do l = isz%beg, isz%end
                do k = isy%beg, isy%end
                    do j = isx%beg, isx%end

                        ! High order flux ======================================
                        w1L = gL_z(l, j, k, 1)
                        w2L = gL_z(l, j, k, 2)
                        w3L = 0d0
                        if (p > 0) w3L = gL_z(l, j, k, 3)

                        w1R = gR_z(l + 1, j, k, 1)
                        w2R = gR_z(l + 1, j, k, 2)
                        w3R = 0d0
                        if (p > 0) w3R = gR_z(l + 1, j, k, 3)

                        normWL = gL_z(l, j, k, num_dims + 1)
                        normWR = gR_z(l + 1, j, k, num_dims + 1)

                        w1 = (w1L + w1R)/2d0
                        w2 = (w2L + w2R)/2d0
                        w3 = (w3L + w3R)/2d0
                        normW = (normWL + normWR)/2d0

                        call s_compute_capilary_stress_tensor(w1, w2, w3, omega, normW)

                        do i = 1, num_dims

                            flux_src_H(momxb + i - 1) = &
                                flux_src_H(momxb + i - 1) + Omega(3,i)

                            flux_src_H(E_idx) = flux_src_H(E_idx) + &
                                Omega(3,i)*vSrc_rsz_vf(l, j, k, i)
                                
                        end do

                        flux_src_H(E_idx) = flux_src_H(E_idx) + &
                            sigma*c_divs%vf(num_dims + 1)%sf(j,k,l)*vSrc_rsz_vf(l, j, k, 3)

                        call s_compute_flux_lim(j, k, l, phi, q_prim_vf, flux_src_vf, id)

                        ! Computing final flux =================================
                        call s_compute_flux_lim(j, k, l, phi, q_prim_vf, flux_src_vf, id)

                        do i = 1, num_dims
                            flux_src_L(momxb + i - 1) = flux_src_L(momxb + i - 1) + &
                                flux_src_L(momxb + i - 1) + phi*(flux_src_H(momxb + i - 1) - &
                                flux_src_L(momxb + i - 1))
                        end do

                        flux_src_L(E_idx) = flux_src_L(E_idx) + flux_src_L(E_idx) + &
                            phi*(flux_src_H(E_idx) - flux_src_L(E_idx))

                    end do
                end do
            end do

        end if

    end subroutine s_compute_capilary_source_flux

    subroutine s_compute_capilary_stress_tensor(w1, w2, w3, omega, normW)
        !$acc routine seq
        real(kind(0d0)) :: w1, w2, w3, normW
        real(kind(0d0)), dimension(num_dims, num_dims) :: omega

        if (normW > 1d-6) then
            Omega(1,1) = -sigma*(w2*w2 + w3*w3) / normW
        
            Omega(2,1) = sigma*w1*w2 / normW
            Omega(1,2) = Omega(2,1)
        
            Omega(2,2) = -sigma*(w1*w1 + w3*w3) / normW
            
            if (p > 0) then
        
                Omega(3,1) = sigma*w1*w3 / normW
                Omega(1,3) = Omega(3,1)
        
                Omega(3,2) = sigma*w2*w3 / normW
                Omega(2,3) = Omega(3,2)
        
                Omega(3,3) = -sigma*(w1*w1 + w2*w2) / normW
        
            end if
        else
            Omega(1,1) = 0d0
            Omega(1,2) = 0d0
            Omega(2,1) = 0d0
            Omega(2,2) = 0d0
            if (p > 0) then
                Omega(1,3) = 0d0
                Omega(3,1) = 0d0
                Omega(2,3) = 0d0
                Omega(3,2) = 0d0
                Omega(3,3) = 0d0
            end if
        end if

    end subroutine s_compute_capilary_stress_tensor

    subroutine s_compute_flux_lim(j, k, l, phi, q_prim_vf, flux_src_vf, id)
        !$acc routine seq
        real(kind(0d0)), intent(out) :: phi
        type(scalar_field), dimension(sys_size) :: q_prim_vf
        type(scalar_field), dimension(sys_size) :: flux_src_vf
        integer :: j, k, l, id, flux_lim
        real(kind(0d0)) :: slope, top, bottom

        if (id == 1) then
            if (flux_src_vf(advxb)%sf(j, k, l) <= 0d0) then
                top = q_prim_vf(advxb)%sf(j, k, l) - &
                    q_prim_vf(advxb)%sf(j - 1, k, l)
                bottom = q_prim_vf(advxb)%sf(j + 1, k, l) - &
                    q_prim_vf(advxb)%sf(j, k, l)
            else
                top = q_prim_vf(advxb)%sf(j + 2, k, l) - &
                    q_prim_vf(advxb)%sf(j + 1, k, l)
                bottom = q_prim_vf(advxb)%sf(j + 1, k, l) - &
                    q_prim_vf(advxb)%sf(j, k, l)    
            end if
        elseif(id == 2) then
            if (flux_src_vf(advxb)%sf(j, k, l) <= 0d0) then
                top = q_prim_vf(advxb)%sf(j, k, l) - &
                    q_prim_vf(advxb)%sf(j, k - 1, l)
                bottom = q_prim_vf(advxb)%sf(j, k + 1, l) - &
                    q_prim_vf(advxb)%sf(j, k, l)
            else
                top = q_prim_vf(advxb)%sf(j, k + 2, l) - &
                    q_prim_vf(advxb)%sf(j, k + 1, l)
                bottom = q_prim_vf(advxb)%sf(j, k + 1, l) - &
                    q_prim_vf(advxb)%sf(j, k, l)    
            end if
        else
            if (flux_src_vf(advxb)%sf(j, k, l) <= 0d0) then
                top = q_prim_vf(advxb)%sf(j, k, l) - &
                    q_prim_vf(advxb)%sf(j, k, l - 1)
                bottom = q_prim_vf(advxb)%sf(j, k, l + 1) - &
                    q_prim_vf(advxb)%sf(j, k, l)
            else
                top = q_prim_vf(advxb)%sf(j, k, l + 2) - &
                    q_prim_vf(advxb)%sf(j, k, l + 1)
                bottom = q_prim_vf(advxb)%sf(j, k, l + 1) - &
                    q_prim_vf(advxb)%sf(j, k, l)    
            end if
        end if

        if (abs(top) < 1d-8) top = 0d0
        if (abs(bottom) < 1d-8) bottom = 0d0

        if (top == bottom) then
            slope = 1d0
        else
            slope = (top*bottom)/max(bottom**2d0,sgm_eps)
        end if

        flux_lim = 4

        ! flux limiter function
        if (flux_lim == 1) then ! minmod (mm)
            phi = max(0d0,min(1d0,slope))
        elseif (flux_lim == 2) then ! muscl (mc)
            phi = max(0d0,min(2d0*slope,5d-1*(1d0+slope),2d0))
        elseif (flux_lim == 3) then ! ospre (op)
            phi = (15d-1*(slope**2d0+slope))/(slope**2d0+slope+1d0)
        elseif (flux_lim == 4) then ! superbee (sb)
            phi = max(0d0,min(1d0,2d0*slope),min(slope,2d0))
        elseif (flux_lim == 5) then ! sweby (sw) (beta = 1.5)
            phi = max(0d0,min(15d-1*slope,1d0),min(slope,15d-1))
        elseif (flux_lim == 6) then ! van albada (va)
            phi = (slope**2d0+slope)/(slope**2d0+1d0)
        elseif (flux_lim == 7) then ! van leer (vl)
            phi = (abs(slope) + slope)/(1d0 + abs(slope))
        end if
    end subroutine s_compute_flux_lim

    subroutine s_get_capilary(q_prim_vf)

        type(scalar_field), dimension(sys_size) :: q_prim_vf
        type(int_bounds_info) :: isx, isy, isz

        isx%beg = -1; isy%beg = 0; isz%beg = 0

        if (m > 0) isy%beg = -1; if (p > 0) isz%beg = -1

        isx%end = m; isy%end = n; isz%end = p

        iv%beg = c_idx; iv%end = c_idx
        do i = 1, num_dims
            call s_reconstruct_cell_boundary_values_capilary(q_prim_vf, cL_x, cL_y, cL_z, cR_x, cR_y, cR_z, i)
        end do

        ! compute gradient components
        !$acc parallel loop collapse(3) gang vector default(present)
        do l = 0, p
            do k = 0, n
                do j = 0, m
                    c_divs%vf(1)%sf(j, k, l) = 1d0/(x_cc(j+1) - x_cc(j-1)) * &
                        (q_prim_vf(c_idx)%sf(j + 1, k, l) - q_prim_vf(c_idx)%sf(j-1, k, l)) 
                    ! c_divs%vf(1)%sf(j, k, l) = &
                    !     1d0/((1d0+wa_flg)*dx(j)) &
                    !     *( wa_flg*cL_x(j + 1, k, l, c_idx) &
                    !     +        cR_x(j, k, l, c_idx) &
                    !     -        cL_x(j, k, l, c_idx) &
                    !     - wa_flg*cR_x(j - 1, k, l, c_idx))
                end do
            end do
        end do

        !$acc parallel loop collapse(3) gang vector default(present)
        do l = 0, p
            do k = 0, n
                do j = 0, m
                    c_divs%vf(2)%sf(j, k, l) = 1d0/(y_cc(j+1) -y_cc(j-1)) * &
                        (q_prim_vf(c_idx)%sf(j, k + 1,  l) - q_prim_vf(c_idx)%sf(j, k-1, l))
                    ! c_divs%vf(2)%sf(j, k, l) = &
                    !     1d0/((1d0+wa_flg)*dy(k)) &
                    !     *( wa_flg*cL_x(j, k + 1, l,c_idx) &
                    !     +        cR_x(j, k, l, c_idx) &
                    !     -        cL_x(j, k, l, c_idx) &
                    !     - wa_flg*cR_x(j, k - 1, l, c_idx) )
                end do
            end do
        end do

        if (p > 0) then
            !$acc parallel loop collapse(3) gang vector default(present)
            do l = 0, p
                do k = 0, n
                    do j = 0, m
                        ! c_divs%vf(3)%sf(j, k, l) = 1d0/(z_cc(l+1) - z_cc(l-1)) * &
                        !     (q_prim_vf(c_idx)%sf(j, k, l+1) - q_prim_vf(c_idx)%sf(j, k, l-1))
                    end do
                end do
            end do
        end if

        !$acc parallel loop collapse(3) gang vector default(present)
        do l = 0, p
            do k = 0, n
                do j = 0, m
                    c_divs%vf(num_dims + 1)%sf(j, k, l) = 0d0
                    !s$acc loop seq
                    do i = 1, num_dims
                        c_divs%vf(num_dims + 1)%sf(j, k, l) = &
                            c_divs%vf(num_dims + 1)%sf(j, k, l) + &
                            c_divs%vf(i)%sf(j, k, l) ** 2d0
                    end do
                    c_divs%vf(num_dims + 1)%sf(j, k, l) = &
                        sqrt(c_divs%vf(num_dims + 1)%sf(j, k, l))
                end do
            end do
        end do

        call s_populate_capillary_buffers()

        iv%beg = 1; iv%end = num_dims + 1

        ! reconstruct gradient components at cell boundaries
        do i = 1, num_dims
            call s_reconstruct_cell_boundary_values_capilary(c_divs%vf, gL_x, gL_y, gL_z, gR_x, gR_y, gR_z, i)
        end do

    end subroutine s_get_capilary

    subroutine s_populate_capillary_buffers()

        ! x - direction
        if (bc_x%beg <= -3) then !< ghost cell extrapolation   
            !$acc parallel loop collapse(4) gang vector default(present)
            do i = 1, num_dims + 1
                do l = 0, p
                    do k = 0, n
                        do j = 1, buff_size
                            c_divs%vf(i)%sf(-j, k, l) = &
                                c_divs%vf(i)%sf(0, k, l)
                        end do
                    end do
                end do
            end do
        elseif (bc_x%beg == -2) then !< Reflective
            !$acc parallel loop collapse(4) gang vector default(present)
            do i = 1, num_dims + 1
                do l = 0, p
                    do k = 0, n
                        do j = 1, buff_size
                            c_divs%vf(i)%sf(-j, k, l) = &
                                c_divs%vf(i)%sf(j - 1, k, l)
                        end do
                    end do
                end do
            end do
        elseif (bc_x%beg == -14) then
            !$acc parallel loop collapse(4) gang vector default(present)
            do i = 1, num_dims + 1
                do l = 0, p
                    do k = 0, n
                        do j = 1, buff_size
                            c_divs%vf(i)%sf(-j, k, l) = &
                                c_divs%vf(i)%sf(0, k, l)
                        end do
                    end do
                end do
            end do
        elseif (bc_x%beg == -1) then
            !$acc parallel loop collapse(4) gang vector default(present)
            do i = 1, num_dims + 1
                do l = 0, p
                    do k = 0, n
                        do j = 1, buff_size
                            c_divs%vf(i)%sf(-j, k, l) = &
                                c_divs%vf(i)%sf(m - (j - 1), k, l)
                        end do
                    end do
                end do
            end do
        else
            call s_mpi_sendrecv_capilary_variables_buffers(c_divs%vf, 1, -1)
        end if

        if (bc_x%end <= -3) then !< ghost-cell extrapolation
            !$acc parallel loop collapse(4) gang vector default(present)
            do i = 1, num_dims + 1
                do l = 0, p
                    do k = 0, n
                        do j = 1, buff_size
                            c_divs%vf(i)%sf(m + j, k, l) = &
                                c_divs%vf(i)%sf(m, k, l)
                        end do
                    end do
                end do
            end do
        elseif (bc_x%end == -2) then
            !$acc parallel loop collapse(4) default(present)
            do i = 1, num_dims + 1
                do l = 0, p
                    do k = 0, n
                        do j = 1, buff_size
                            c_divs%vf(i)%sf(m + j, k, l) = &
                                c_divs%vf(i)%sf(m - (j - 1), k, l)
                        end do
                    end do
                end do
            end do
        elseif (bc_x%end == -14) then
            !$acc parallel loop collapse(4) default(present)
            do i = 1, num_dims + 1
                do l = 0, p
                    do k = 0, n
                        do j = 1, buff_size
                            c_divs%vf(i)%sf(m + j, k, l) = &
                                c_divs%vf(i)%sf(m, k, l)
                        end do
                    end do
                end do
            end do
        else if (bc_x%end == -1) then
            !$acc parallel loop collapse(4) gang vector default(present)
            do i = 1, num_dims + 1
                do l = 0, p
                    do k = 0, n
                        do j = 1, buff_size
                            c_divs%vf(i)%sf(m + j, k, l) = &
                                c_divs%vf(i)%sf(j - 1, k, l)
                        end do
                    end do
                end do
            end do
        else
            call s_mpi_sendrecv_capilary_variables_buffers(c_divs%vf, 1, 1)
        end if

        if (n == 0) then
            return
        elseif (bc_y%beg <= -3) then !< ghost-cell extrapolation
            !$acc parallel loop collapse(4) gang vector default(present)
            do i = 1, num_dims + 1
                do k = 0, p
                    do j = 1, buff_size
                        do l = -buff_size, m + buff_size
                            c_divs%vf(i)%sf(l, -j, k) = &
                                c_divs%vf(i)%sf(l, 0, k)
                        end do
                    end do
                end do
            end do
        elseif (bc_y%beg == -2) then !< Reflective
            !$acc parallel loop collapse(4) gang vector default(present)
            do i = 1, num_dims + 1
                do k = 0, p
                    do j = 1, buff_size
                        do l = -buff_size, m + buff_size
                            c_divs%vf(i)%sf(l, -j, k) = &
                                c_divs%vf(i)%sf(l, j - 1, k)
                        end do
                    end do
                end do
            end do
        elseif (bc_y%beg == -14) then
            !$acc parallel loop collapse(4) gang vector default(present)
            do i = 1, num_dims + 1
                do k = 0, p
                    do j = 1, buff_size
                        do l = -buff_size, m + buff_size
                            c_divs%vf(i)%sf(l, -j, k) = &
                                c_divs%vf(i)%sf(l, 0, k)
                        end do
                    end do
                end do
            end do
        elseif (bc_y%beg == -1) then
            !$acc parallel loop collapse(4) gang vector default(present)
            do i = 1, num_dims + 1
                do k = 0, p
                    do j = 1, buff_size
                        do l = -buff_size, m + buff_size
                            c_divs%vf(i)%sf(l, -j, k) = &
                                c_divs%vf(i)%sf(l, n - (j - 1), k)
                        end do
                    end do
                end do
            end do
        else
            call s_mpi_sendrecv_capilary_variables_buffers(c_divs%vf, 2, -1)
        endif

        if (bc_y%end <= -3) then !< ghost-cell extrapolation
            !$acc parallel loop collapse(4) gang vector default(present)
            do i = 1, num_dims + 1
                do k = 0, p
                    do j = 1, buff_size
                        do l = -buff_size, m + buff_size
                            c_divs%vf(i)%sf(l, n + j, k) = &
                                c_divs%vf(i)%sf(l, n, k)
                        end do
                    end do
                end do
            end do
        elseif (bc_y%end == -2) then !< Reflective
            !$acc parallel loop collapse(4) gang vector default(present)
            do i = 1, num_dims + 1
                do k = 0, p
                    do j = 1, buff_size
                        do l = -buff_size, m + buff_size
                            c_divs%vf(i)%sf(l, n + j, k) = &
                                c_divs%vf(i)%sf(l, n - (j - 1), k)
                        end do
                    end do
                end do
            end do
        elseif (bc_y%end == -14) then
            !$acc parallel loop collapse(4) gang vector default(present)
            do i = 1, num_dims + 1
                do k = 0, p
                    do j = 1, buff_size
                        do l = -buff_size, m + buff_size
                            c_divs%vf(i)%sf(l, n + j, k) = &
                                c_divs%vf(i)%sf(l, n, k)
                        end do
                    end do
                end do
            end do
        elseif (bc_y%end == -1) then
            !$acc parallel loop collapse(4) gang vector default(present)
            do i = 1, num_dims + 1
                do k = 0, p
                    do j = 1, buff_size
                        do l = -buff_size, m + buff_size
                            c_divs%vf(i)%sf(l, n + j, k) = &
                                c_divs%vf(i)%sf(l, j - 1, k)
                        end do
                    end do
                end do
            end do
        else
            call s_mpi_sendrecv_capilary_variables_buffers(c_divs%vf, 2, 1)
        end if

        if (p == 0) then
            return
        elseif (bc_z%beg <= -3) then !< ghost-cell extrapolation
            !$acc parallel loop collapse(4) gang vector default(present)
            do i = 1, num_dims + 1
                do j = 1, buff_size
                    do l = -buff_size, n + buff_size
                        do k = -buff_size, m + buff_size
                            c_divs%vf(i)%sf(k, l, -j) = &
                                c_divs%vf(i)%sf(k, l, 0)
                        end do
                    end do
                end do
            end do
        elseif (bc_z%beg == -14) then !< slip wall
            !$acc parallel loop collapse(4) gang vector default(present)
            do i = 1, num_dims + 1
                do j = 1, buff_size
                    do l = -buff_size, n + buff_size
                        do k = -buff_size, m + buff_size
                            c_divs%vf(i)%sf(k, l, -j) = &
                                c_divs%vf(i)%sf(k, l, 0)
                        end do
                    end do
                end do
            end do
        elseif (bc_z%beg == -2) then !< Reflective
            !$acc parallel loop collapse(4) gang vector default(present)
            do i = 1, num_dims + 1 
                do j = 1, buff_size
                    do l = -buff_size, n + buff_size
                        do k = -buff_size, m + buff_size
                                c_divs%vf(i)%sf(k, l, -j) = &
                                    c_divs%vf(i)%sf(k, l, j - 1)
                        end do
                    end do
                end do
            end do
        elseif (bc_z%beg == -1) then
            !$acc parallel loop collapse(4) gang vector default(present)
            do i = 1, num_dims + 1
                do j = 1, buff_size
                    do l = -buff_size, n + buff_size
                        do k = -buff_size, m + buff_size
                            c_divs%vf(i)%sf(k, l, -j) = &
                                c_divs%vf(i)%sf(k, l, p - (j - 1))
                        end do
                    end do
                end do
            end do
        else
            call s_mpi_sendrecv_capilary_variables_buffers(c_divs%vf, 3, -1)
        end if
        
        if (bc_z%end <= -3) then !< ghost-cell extrapolation
            !$acc parallel loop collapse(4) gang vector default(present)
            do i = 1, num_dims + 1
                do j = 1, buff_size
                    do l = -buff_size, n + buff_size
                        do k = -buff_size, m + buff_size
                            c_divs%vf(i)%sf(k, l, p + j) = &
                                c_divs%vf(i)%sf(k, l, p)
                        end do
                    end do
                end do
            end do
        else if ( bc_Z%end == -14 ) then !< slip wall
            !$acc parallel loop collapse(4) gang vector default(present)
            do i = 1, num_dims + 1
                do j = 1, buff_size
                    do l = -buff_size, n + buff_size
                        do k = -buff_size, m + buff_size
                            c_divs%vf(i)%sf(k, l, p+j) = &
                                c_divs%vf(i)%sf(k, l, p)
                        end do
                    end do
                end do
            end do
        elseif (bc_z%end == -2) then !< Reflective
            !$acc parallel loop collapse(4) gang vector default(present)
            do i = 1, num_dims + 1
                do j = 1, buff_size
                    do l = -buff_size, n + buff_size
                        do k = -buff_size, m + buff_size
                            c_divs%vf(i)%sf(k, l, p + j) = &
                                c_divs%vf(i)%sf(k, l, p - (j - 1))
                        end do
                    end do
                end do
            end do
        elseif (bc_z%end == -1) then
            !$acc parallel loop collapse(4) gang vector default(present)
            do i = 1, num_dims + 1
                do j = 1, buff_size
                    do l = -buff_size, n + buff_size
                        do k = -buff_size, m + buff_size
                            c_divs%vf(i)%sf(k, l, p + j) = &
                                c_divs%vf(i)%sf(k, l, j - 1)
                        end do
                    end do
                end do
            end do
        else
            call s_mpi_sendrecv_capilary_variables_buffers(c_divs%vf, 3, 1)
        endif

    end subroutine s_populate_capillary_buffers

    subroutine s_reconstruct_cell_boundary_values_capilary(v_vf, vL_x, vL_y, vL_z, vR_x, vR_y, vR_z, & ! -
                                                             norm_dir)

        type(scalar_field), dimension(iv%beg:iv%end), intent(IN) :: v_vf

        real(kind(0d0)), dimension(startx:, starty:, startz:, iv%beg:), intent(INOUT) :: vL_x, vL_y, vL_z, vR_x, vR_y, vR_z 

        integer, intent(IN) :: norm_dir

        integer :: weno_dir !< Coordinate direction of the WENO reconstruction

        integer :: i, j, k, l
        ! Reconstruction in s1-direction ===================================

        if (norm_dir == 1) then
            is1 = ix; is2 = iy; is3 = iz
            weno_dir = 1; is1%beg = is1%beg + weno_polyn
            is1%end = is1%end - weno_polyn

        elseif (norm_dir == 2) then
            is1 = iy; is2 = ix; is3 = iz
            weno_dir = 2; is1%beg = is1%beg + weno_polyn
            is1%end = is1%end - weno_polyn

        else
            is1 = iz; is2 = iy; is3 = ix
            weno_dir = 3; is1%beg = is1%beg + weno_polyn
            is1%end = is1%end - weno_polyn

        end if

        !$acc update device(is1, is2, is3, iv)

        ! if (n > 0) then
        !     if (p > 0) then

        !         call s_weno(v_vf(iv%beg:iv%end), &
        !             vL_x(:, :, :, iv%beg:iv%end), vL_y(:, :, :, iv%beg:iv%end), vL_z(:, :, :, iv%beg:iv%end), vR_x(:, :, :, iv%beg:iv%end), vR_y(:, :, :, iv%beg:iv%end), vR_z(:, :, :, iv%beg:iv%end), &
        !                         norm_dir, weno_dir, &
        !                         is1, is2, is3)
        !     else
        !         call s_weno(v_vf(iv%beg:iv%end), &
        !             vL_x(:, :, :, iv%beg:iv%end), vL_y(:, :, :, iv%beg:iv%end), vL_z(:, :, :, :), vR_x(:, :, :, iv%beg:iv%end), vR_y(:, :, :, iv%beg:iv%end), vR_z(:, :, :, :), &
        !                         norm_dir, weno_dir, &
        !                         is1, is2, is3)
        !     end if
        ! else

        !     call s_weno(v_vf(iv%beg:iv%end), &
        !                 vL_x(:, :, :, iv%beg:iv%end), vL_y(:, :, :, :), vL_z(:, :, :, :), vR_x(:, :, :, iv%beg:iv%end), vR_y(:, :, :, :), vR_z(:, :, :, :), &
        !                     norm_dir, weno_dir, &
        !                     is1, is2, is3)
        ! end if

        if (weno_dir == 1) then
!$acc parallel loop collapse(4) default(present)
            do i = iv%beg, iv%end 
                    do l = is3%beg, is3%end
                        do k = is2%beg, is2%end
                            do j = is1%beg, is1%end
                                vL_x(j, k, l, i) = v_vf(i)%sf(j, k, l)
                                vR_x(j, k, l, i) = v_vf(i)%sf(j, k, l)
                            end do
                        end do
                    end do
                end do
                !$acc end parallel loop
        else if (weno_dir == 2) then
!$acc parallel loop collapse(4) default(present)
            do i = iv%beg, iv%end
                do l = is3%beg, is3%end
                    do k = is2%beg, is2%end
                        do j = is1%beg, is1%end
                            vL_y(j, k, l, i) = v_vf(i)%sf(k, j, l)
                            vR_y(j, k, l, i) = v_vf(i)%sf(k, j, l)
                        end do
                    end do
                end do
            end do
            !$acc end parallel loop
        else if (weno_dir == 3) then
!$acc parallel loop collapse(4) default(present)
            do i = iv%beg, iv%end
                do l = is3%beg, is3%end
                    do k = is2%beg, is2%end
                        do j = is1%beg, is1%end
                            vL_z(j, k, l, i) = v_vf(i)%sf(l, k, j)
                            vR_z(j, k, l, i) = v_vf(i)%sf(l, k, j)
                        end do
                    end do
                end do
            end do
            !$acc end parallel loop
        end if

    end subroutine s_reconstruct_cell_boundary_values_capilary

    subroutine s_finalize_surface_tension_module()

        do j = 1,num_dims
            @:DEALLOCATE(c_divs%vf(j)%sf)
        end do

        @:DEALLOCATE(c_divs%vf)

        @:DEALLOCATE(gL_x, gR_x, cL_x, cR_x)

        @:DEALLOCATE(gL_y, gR_y, cL_y, cR_y)
        
        if (p > 0) then
            @:DEALLOCATE(gL_z, gR_z, cL_z, cR_z)
        end if

        !@:DEALLOCATE(flux_src_H(momxb:E_idx), flux_src_L(momxb:E_idx))

    end subroutine s_finalize_surface_tension_module

end module m_surface_tension
