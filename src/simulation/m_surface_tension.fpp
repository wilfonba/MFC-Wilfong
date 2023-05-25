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

    real(kind(0d0)), allocatable, dimension(:, :) :: Omega
    !$acc declare create(c_divs, Omega)

    real(kind(0d0)), allocatable, dimension(:,:,:,:) :: gL_x, gR_x, gL_y, gR_y, gL_z, gR_z
    !$acc declare create(gL_x, gR_x, gL_y, gR_y, gL_z, gR_z)

    type(int_bounds_info) :: ix, iy, iz, is1, is2, is3, iv
    !$acc declare create(ix, iy, iz, is1, is2, is3)

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

        if (m > 0) then
            @:ALLOCATE(gL_y(iy%beg:iy%end, ix%beg:ix%end, iz%beg:iz%end, num_dims + 1))
            @:ALLOCATE(gR_y(iy%beg:iy%end, ix%beg:ix%end, iz%beg:iz%end, num_dims + 1))
            if (p > 0) then
                @:ALLOCATE(gL_z(iz%beg:iz%end, iy%beg:iy%end, ix%beg:ix%end, num_dims + 1))
                @:ALLOCATE(gR_z(iz%beg:iz%end, iy%beg:iy%end, ix%beg:ix%end, num_dims + 1))
            end if
        endif

        @:ALLOCATE(Omega(1:num_dims, 1:num_dims))

    end subroutine s_initialize_surface_tension_module

    subroutine s_compute_capilary_source_flux(q_prim_vf, &
                                              vSrc_rsx, vSrc_rsy, vSrc_rsz, &
                                              flux_src_vf, &
                                              id, isx, isy, isz)

        type(scalar_field), dimension(sys_size) :: q_prim_vf
        real(kind(0d0)), dimension(startx:,starty:,startz:,1:) :: vSrc_rsx, vSrc_rsy, vSrc_rsz
        type(scalar_field), &
            dimension(sys_size), &
            intent(INOUT) :: flux_src_vf
        type(int_bounds_info) :: isx, isy, isz
        integer :: id

        real(kind(0d0)) :: e0L, e0R

        if (id == 1) then
            dir_idx = (/1, 2, 3/)
        elseif (id == 2) then
            dir_idx = (/2, 1, 3/)
        else
            dir_idx = (/3, 1, 2/)
        end if

        if (id == 1) then

            !$acc parallel loop collapse(3) gang vector default(present) private(Omega)
            do l = isz%beg, isz%end
                do k = isy%beg, isy%end
                    do j = isx%beg, isx%end
                        call s_compute_capilary_stress_tensors(gL_x, j, k, l)

                        do i = 1, num_dims
                            flux_src_vf(momxb + i - 1)%sf(j, k, l) = &
                                flux_src_vf(momxb + i - 1)%sf(j, k, l) + Omega(1,i)

                            flux_src_vf(E_idx)%sf(j,k,l) = flux_src_vf(E_idx)%sf(j,k,l) + &
                                Omega(1,i)*vSrc_rsx(j,k,l,i)
                        end do

                        e0L = sigma*gL_x(j,k,l,num_dims+1)
                        flux_src_vf(E_idx)%sf(j,k,l) = flux_src_vf(E_idx)%sf(j,k,l) + &
                            e0L * vSrc_rsx(j,k,l,1)

                    end do
                end do
            end do

        elseif (id == 2) then

            !$acc parallel loop collapse(3) gang vector default(present) private(Omega)
            do l = isz%beg, isz%end
                do k = isy%beg, isy%end
                    do j = isx%beg, isx%end
                        call s_compute_capilary_stress_tensors(gL_y, j, k, l)

                        do i = 1, num_dims
                            flux_src_vf(momxb + i - 1)%sf(j, k, l) = &
                                flux_src_vf(momxb + i - 1)%sf(j, k, l) + Omega(2,i)

                            flux_src_vf(E_idx)%sf(j,k,l) = flux_src_vf(E_idx)%sf(j,k,l) + &
                                Omega(2,i)*vSrc_rsy(j, k, l, i)
                        end do

                        e0L = sigma*gL_y(j, k, l,num_dims+1)
                        flux_src_vf(E_idx)%sf(j,k,l) = flux_src_vf(E_idx)%sf(j,k,l) + &
                            e0L * vSrc_rsy(j, k, l, 2)

                    end do
                end do
            end do

        elseif (id == 3) then

            !$acc parallel loop collapse(3) gang vector default(present) private(Omega)
            do l = isz%beg, isz%end
                do k = isy%beg, isy%end
                    do j = isx%beg, isx%end
                        call s_compute_capilary_stress_tensors(gL_z, j, k, l)

                        do i = 1, num_dims
                            flux_src_vf(momxb + i - 1)%sf(j, k, l) = &
                                flux_src_vf(momxb + i - 1)%sf(j, k, l) + Omega(3,i)

                            flux_src_vf(E_idx)%sf(j,k,l) = flux_src_vf(E_idx)%sf(j,k,l) + &
                                Omega(3,i)*vSrc_rsz(j,k,l,i)
                        end do

                        e0L = sigma*gL_z(j,k,l,num_dims+1)
                        flux_src_vf(E_idx)%sf(j,k,l) = flux_src_vf(E_idx)%sf(j,k,l) + &
                            e0L * vSrc_rsz(j,k,l,3)

                    end do
                end do
            end do

        end if

    end subroutine s_compute_capilary_source_flux

    subroutine s_get_capilary(q_prim_vf, ix, iy, iz)

        type(scalar_field), dimension(sys_size) :: q_prim_vf
        type(int_bounds_info) :: ix, iy, iz

        ! compute gradient components
        !$acc parallel loop collapse(3) gang vector default(present)
        do l = iz%beg, iz%end
            do k = iy%beg, iy%end
                do j = ix%beg, ix%end
                    c_divs%vf(1)%sf(j, k, l) = &
                        (q_prim_vf(c_idx)%sf(j - 2, k, l) &
                        - 8d0*q_prim_vf(c_idx)%sf(j - 1, k, l) &
                        + 8d0*q_prim_vf(c_idx)%sf(j + 1, k, l) &
                        - q_prim_vf(c_idx)%sf(j + 2, k, l)) &
                        /(12d0*dx(j))
                end do
            end do
        end do

        ! 2D
        if (m > 0) then
            !$acc parallel loop collapse(3) gang vector default(present)
            do l = iz%beg, iz%end
                do k = iy%beg, iy%end
                    do j = ix%beg, ix%end
                        c_divs%vf(2)%sf(j, k, l) = &
                            (q_prim_vf(c_idx)%sf(j, k - 2, l) &
                            - 8d0*q_prim_vf(c_idx)%sf(j, k - 1, l) &
                            + 8d0*q_prim_vf(c_idx)%sf(j, k + 1, l) &
                            - q_prim_vf(c_idx)%sf(j, k + 2, l)) &
                            /(12d0*dy(k))
                    end do
                end do
            end do
            ! 3D
            if (p > 0) then
                !$acc parallel loop collapse(3) gang vector default(present)
                do l = iz%beg, iz%end
                    do k = iy%beg, iy%end
                        do j = ix%beg, ix%end
                            c_divs%vf(3)%sf(j, k, l) = &
                                (q_prim_vf(c_idx)%sf(j, k, l - 2) &
                                - 8d0*q_prim_vf(c_idx)%sf(j, k, l - 1) &
                                + 8d0*q_prim_vf(c_idx)%sf(j, k, l + 1) &
                                - q_prim_vf(c_idx)%sf(j, k, l + 2)) &
                                /(12d0*dz(l))
                        end do
                    end do
                end do
            end if
        end if

        iv%beg = 1; iv%end = num_dims

        ! reconstruct gradient components at cell boundaries
        do i = 1, num_dims
            call s_reconstruct_cell_boundary_values_capilary(c_divs%vf(1:num_dims), gL_x, gL_y, gL_z, gR_x, gR_y, gR_z, i)
        end do

        ! Compute gradient magnitude at cell boundaries
        !$acc parallel loop collapse(3) gang vector default(present)
        do l = iz%beg, iz%end
            do k = iy%beg, iy%end
                do j = ix%beg, ix%end
                    gL_x(j, k, l, num_dims+1) = 0d0
                    do i = 1, num_dims
                        gL_x(j, k, l, num_dims+1) = &
                            gL_x(j, k, l, num_dims+1) + &
                            gL_x(j, k, l, i)**2d0
                    end do
                    gL_x(j, k, l, num_dims+1) = sqrt(gL_x(j, k, l, num_dims+1))
                end do
            end do
        end do
        if (m > 0) then
            !$acc parallel loop collapse(3) gang vector default(present)
            do l = iz%beg, iz%end
                do k = iy%beg, iy%end
                    do j = ix%beg, ix%end
                        gL_y(j, k, l, num_dims+1) = 0d0
                        do i = 1, num_dims
                            gL_y(j, k, l, num_dims+1) = &
                                gL_y(j, k, l, num_dims+1) + &
                                gL_y(j, k, l, i)**2d0
                        end do
                        gL_y(j, k, l, num_dims+1) = sqrt(gL_y(j, k, l, num_dims+1))
                    end do
                end do
            end do
            if (p > 0) then
                !$acc parallel loop collapse(3) gang vector default(present)
                do l = iz%beg, iz%end
                    do k = iy%beg, iy%end
                        do j = ix%beg, ix%end
                            gL_z(l, k, j, num_dims+1) = 0d0
                            do i = 1, num_dims
                                gL_z(l, k, j, num_dims+1) = &
                                    gL_z(l, k, j, num_dims+1) + &
                                    gL_z(l, k, j, i)**2d0
                            end do
                            gL_z(l, k, j, num_dims+1) = sqrt(gL_z(l, k, j, num_dims+1))
                        end do
                    end do
                end do
            end if
        end if

    end subroutine s_get_capilary

    subroutine s_compute_capilary_stress_tensors(gL, j, k, l)

        integer :: j, k, l
        real(kind(0d0)), dimension(startx:,starty:,startz:,:) :: gL

        !$acc routine seq
        Omega(1,1) = -sigma*(gL(j,k,l,2)**2) / &
                        max(gL(j,k,l,num_dims+1),sgm_eps)
        if (m > 0) then
            Omega(1,2) = sigma*gL(j, k, l, 1) * gL(j, k, l, 2) / &
                                max(gL(j, k, l, num_dims + 1),sgm_eps)
            Omega(2,1) = Omega(1,2)

            Omega(2,2) = -sigma*(gL(j,k,l,1)**2) / &
                            max(gL(j,k,l,num_dims+1),sgm_eps)
            if (p > 0) then
                Omega(1,3) = sigma*(gL(j, k, l, 1) * gL(j, k, l, 3) / &
                                    max(gL(j, k, l, num_dims + 1),sgm_eps))
                Omega(3,1) = Omega(1,3)

                Omega(2,3) = sigma*gL(j, k, l, 2) * gL(j, k, l, 3) / &
                                    max(gL(j, k, l, num_dims + 1),sgm_eps)
                Omega(3,2) = Omega(2,3)

                Omega(3,3) = -sigma*(gL(j, k, l, 1)**2 + gL(j, k, l, 2)**2) / &
                                    max(gL(j, k, l, num_dims+1),sgm_eps)
                Omega(1,1) = -sigma*(gL(j, k, l, 2)**2 + gL(j, k, l, 3)**2) / &
                                    max(gL(j, k, l, num_dims+1), sgm_eps)
                Omega(2,2) = -sigma*(gL(j, k, l, 1)**2 + gL(j, k, l, 3)**2) / &
                                    max(gL(j, k, l, num_dims+1),sgm_eps) 
            end if
        end if

    end subroutine s_compute_capilary_stress_tensors

        !>  The purpose of this subroutine is to WENO-reconstruct the
        !!      left and the right cell-boundary values, including values
        !!      at the Gaussian quadrature points, from the cell-averaged
        !!      variables.
        !!  @param v_vf Cell-average variables
        !!  @param vL_qp Left WENO-reconstructed, cell-boundary values including
        !!          the values at the quadrature points, of the cell-average variables
        !!  @param vR_qp Right WENO-reconstructed, cell-boundary values including
        !!          the values at the quadrature points, of the cell-average variables
        !!  @param norm_dir Splitting coordinate direction
    subroutine s_reconstruct_cell_boundary_values_capilary(v_vf, vL_x, vL_y, vL_z, vR_x, vR_y, vR_z, & ! -
                                                      norm_dir)

        type(scalar_field), dimension(iv%beg:iv%end), intent(IN) :: v_vf

        real(kind(0d0)), dimension(startx:, starty:, startz:, 1:), intent(INOUT) :: vL_x, vL_y, vL_z, vR_x, vR_y, vR_z

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

        if (n > 0) then
            if (p > 0) then

                call s_weno(v_vf(iv%beg:iv%end), &
                    vL_x(:, :, :, iv%beg:iv%end), vL_y(:, :, :, iv%beg:iv%end), vL_z(:, :, :, iv%beg:iv%end), vR_x(:, :, :, iv%beg:iv%end), vR_y(:, :, :, iv%beg:iv%end), vR_z(:, :, :, iv%beg:iv%end), &
                                norm_dir, weno_dir, &
                                is1, is2, is3)
            else
                call s_weno(v_vf(iv%beg:iv%end), &
                    vL_x(:, :, :, iv%beg:iv%end), vL_y(:, :, :, iv%beg:iv%end), vL_z(:, :, :, :), vR_x(:, :, :, iv%beg:iv%end), vR_y(:, :, :, iv%beg:iv%end), vR_z(:, :, :, :), &
                                norm_dir, weno_dir, &
                                is1, is2, is3)
            end if
        else

            call s_weno(v_vf(iv%beg:iv%end), &
                        vL_x(:, :, :, iv%beg:iv%end), vL_y(:, :, :, :), vL_z(:, :, :, :), vR_x(:, :, :, iv%beg:iv%end), vR_y(:, :, :, :), vR_z(:, :, :, :), &
                            norm_dir, weno_dir, &
                            is1, is2, is3)
        end if

    end subroutine s_reconstruct_cell_boundary_values_capilary  

    subroutine s_finalize_surface_tension_module()

        do j = 1,num_dims
            @:DEALLOCATE(c_divs%vf(j)%sf)
        end do

        @:DEALLOCATE(c_divs%vf)

        @:DEALLOCATE(Omega)
        @:DEALLOCATE(gL_x, gR_x)
        if (n > 0) then
            @:DEALLOCATE(gL_y, gR_y)
            if (p > 0) then
                @:DEALLOCATE(gL_z, gR_z)
            end if
        end if

    end subroutine s_finalize_surface_tension_module

end module m_surface_tension
