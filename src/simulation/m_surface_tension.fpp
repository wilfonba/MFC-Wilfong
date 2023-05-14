#:include 'macros.fpp'

!> @brief This module is used to compute source terms for hypoelastic model
module m_surface_tension

    ! Dependencies =============================================================

    use m_derived_types        !< Definitions of the derived types

    use m_global_parameters    !< Definitions of the global parameters

    use m_mpi_proxy            !< Message passing interface (MPI) module proxy

    use m_rhs
    ! ==========================================================================

    implicit none

    private; public :: s_initialize_surface_tension_module, &
        s_get_surface_tension, &
        s_finalize_surface_tension_module

    type(vector_field) :: c_divs
    type(scalar_field) :: gm

    real(kind(0d0)), allocatable, dimension(:, :) :: OmegaL, OmegaR
    !$acc declare create(c_divs, gm, omegaL, omegaR)

    real(kind(0d0)), allocatable, dimension(:,:,:,:) :: gL_x, gR_x, gL_y, gR_y, gL_z, gR_z
    !$acc declare create(gL_x, gR_x, gL_y, gR_y, gL_z, gR_z)

    real(kind(0d0)) :: gmL, gmR, w1L, w1R, w2L, w2R, w3L, w3R
    !$acc declare create(gmL, gmR, w1L, w1R, w2L, w2R, w3L, w3R)

    type(int_bounds_info) :: ix, iy, iz
    !$acc declare create(ix, iy, iz)

    integer :: j, k, l, q  

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

        @:ALLOCATE(gm%sf(ix%beg:ix%end, iy%beg:iy%end, iz%beg:iz%end))

        @:ALLOCATE(gL_x(ix%beg:ix%end, iy%beg:iy%end, iz%beg:iz%end, num_dims + 1))
        @:ALLOCATE(gR_x(ix%beg:ix%end, iy%beg:iy%end, iz%beg:iz%end, num_dims + 1))

        if (m > 0) then
            @:ALLOCATE(gL_y(ix%beg:ix%end, iy%beg:iy%end, iz%beg:iz%end, num_dims + 1))
            @:ALLOCATE(gR_y(ix%beg:ix%end, iy%beg:iy%end, iz%beg:iz%end, num_dims + 1))
        elseif (p > 0) then
            @:ALLOCATE(gL_z(ix%beg:ix%end, iy%beg:iy%end, iz%beg:iz%end, num_dims + 1))
            @:ALLOCATE(gR_z(ix%beg:ix%end, iy%beg:iy%end, iz%beg:iz%end, num_dims + 1))
        endif

        @:ALLOCATE(OmegaL(1:num_dims, 1:num_dims), OmegaR(1:num_dims, 1:num_dims))

    end subroutine s_initialize_surface_tension_module

    subroutine s_get_surface_tension(q_cons_vf, q_prim_vf)

        type(scalar_field), dimension(sys_size), intent(IN) :: q_prim_vf
        type(scalar_field), dimension(sys_size), intent(IN) :: q_cons_vf
        real(kind(0d0)) :: gm !< grad(c) magnitude

        ! Compute derivatives of the shape function
        !$acc parallel loop collapse(3) gang vector default(present)
        do l = 0, p
            do k = 0, n
                do j = 0, m
                    c_divs%vf(1)%sf(j, k, l) = &
                        (q_prim_vf(c_idx)%sf(j - 2, k, l) &
                         - 8d0*q_prim_vf(c_idx)%sf(j - 1, k, l) &
                         + 8d0*q_prim_vf(c_idx)%sf(j + 1, k, l) &
                         - q_prim_vf(c_idx)%sf(j + 2, k, l)) &
                        /(12d0*dx(j))
                end do
            end do
        end do

        if (m > 0) then
            !$acc parallel loop collapse(3) gang vector default(present)
            do l = 0, p
                do k = 0, n
                    do j = 0, m
                        c_divs%vf(2)%sf(j, k, l) = &
                            (q_prim_vf(c_idx)%sf(j, k - 2, l) &
                            - 8d0*q_prim_vf(c_idx)%sf(j, k - 1, l) &
                            + 8d0*q_prim_vf(c_idx)%sf(j, k + 1, l) &
                            - q_prim_vf(c_idx)%sf(j, k + 2, l)) &
                            /(12d0*dx(j))
                    end do
                end do
            end do

            if (p > 0) then
                !$acc parallel loop collapse(3) gang vector default(present)
                do l = 0, p
                    do k = 0, n
                        do j = 0, m
                            c_divs%vf(3)%sf(j, k, l) = &
                                (q_prim_vf(c_idx)%sf(j, k, l - 2) &
                                - 8d0*q_prim_vf(c_idx)%sf(j, k, l - 1) &
                                + 8d0*q_prim_vf(c_idx)%sf(j, k, l + 1) &
                                - q_prim_vf(c_idx)%sf(j, k, l + 2)) &
                                /(12d0*dx(j))
                        end do
                    end do
                end do
            end if
        end if

        ! Cell centered gradient magnitude
        !$acc parallel loop collapse(3) gang vector default(present)
        do l = 0, p
            do k = 0, n
                do j = 0, m
                    c_divs%vf(num_dims + 1)%sf(j,k,l) = 0d0
                    do q = 1, num_dims
                        c_divs%vf(num_dims + 1)%sf(j,k,l) = &
                            c_divs%vf(num_dims + 1)%sf(j,k,l) + &
                            c_divs%vf(q)%sf(j,k,l)**2d0
                    end do
                    c_divs%vf(num_dims + 1)%sf(j,k,l) = &
                        sqrt(c_divs%vf(num_dims + 1)%sf(j,k,l))
                end do
            end do
        end do

        ! Reconstruct gradient at cell boundaries
        do q = 1, num_dims
            call s_reconstruct_cell_boundary_values(c_divs%vf(1:num_dims), gL_x, gL_y, gL_z, gR_x, gR_y, gR_z, q)
        end do

        ! < x-direction
        !$acc parallel loop collapse(3) gang vector default(present)
        do l = 0, p
            do k = 0, n
                do j = 0, m
                    call s_compute_capillary_stress_tensors(gL_x, gR_x, j, k, l)
                    do q = 1,num_dims
                        q_cons_vf(momxb + q - 1)%sf(j,k,l) = &
                            q_prim_vf(momxb + q - 1)%sf(j,k,l) - dt/dx(j) * &
                            (OmegaR(1,q) - OmegaL(1,q))
                    end do
                    q_cons_vf(E_idx)%sf(j,k,l) = q_prim_vf(E_idx)%sf(j,k,l) - &
                        dt/dx(j) * (sigma*gR_x(j,k,l,num_dims+1)*q_prim_vf(momxb)%sf(j,k,l) + &
                        OmegaR(1,1)*q_prim_vf(momxb)%sf(j,k,l) + &
                        OmegaR(1,2)*q_prim_vf(momxb+1)%sf(j,k,l)) + &
                        !OmegaR(1,3)*q_prim_vf(momxe)%sf(j,k,l)) + &
                        dt/dx(j) * (sigma*gL_x(j,k,l,num_dims+1)*q_prim_vf(momxb)%sf(j,k,l) + &
                        OmegaL(1,1)*q_prim_vf(momxb)%sf(j,k,l) + &
                        OmegaL(1,2)*q_prim_vf(momxb+1)%sf(j,k,l)) ! + &
                        !OmegaL(1,3)*q_prim_vf(momxe)%sf(j,k,l))
                end do
            end do
        end do    
        
        if (m > 0) then
            !$acc parallel loop collapse(3) gang vector default(present)
            do l = 0, p
                do k = 0, n
                    do j = 0, m
                        call s_compute_capillary_stress_tensors(gL_y, gR_y, j, k, l)
                        do q = 1,num_dims
                            q_cons_vf(momxb + q - 1)%sf(j,k,l) = &
                                q_prim_vf(momxb + q - 1)%sf(j,k,l) - dt/dy(j) * &
                                (OmegaR(2,q) - OmegaL(2,q))
                        end do
                        q_cons_vf(E_idx)%sf(j,k,l) = q_prim_vf(E_idx)%sf(j,k,l) - &
                            dt/dy(k) * (sigma*gR_y(j,k,l,num_dims+1)*q_prim_vf(momxb + 1)%sf(j,k,l) + &
                            OmegaR(2,1)*q_prim_vf(momxb)%sf(j,k,l) + &
                            OmegaR(2,2)*q_prim_vf(momxb+1)%sf(j,k,l)) + &
                            !OmegaR(2,3)*q_prim_vf(momxe)%sf(j,k,l)) + &
                            dt/dy(k) * (sigma*gL_y(j,k,l,num_dims+1)*q_prim_vf(momxb + 1)%sf(j,k,l) + &
                            OmegaL(2,1)*q_prim_vf(momxb)%sf(j,k,l) + &
                            OmegaL(2,2)*q_prim_vf(momxb+1)%sf(j,k,l)) ! + &
                            !OmegaL(2,3)*q_prim_vf(momxe)%sf(j,k,l))
                    end do
                end do
            end do
        end if

        ! if (p > 0) then
        !     do l = 0, p
        !         do k = 0, n
        !             do j = 0, m
        !                 call s_compute_capillary_stress_tensors(gL_z, gR_z, j, k, l)
        !                 do q = 1,num_dims
        !                     q_cons_vf(momxb + q - 1)%sf(j,k,l) = &
        !                         q_prim_vf(momxb + q - 1)%sf(j,k,l) - dt/dz(l) * &
        !                         (OmegaR(3,q) - OmegaL(3,q))
        !                 end do
        !                 q_cons_vf(E_idx)%sf(j,k,l) = q_prim_vf(E_idx)%sf(j,k,l) - &
        !                     dt/dz(l) * (sigma*gR_z(j,k,l,num_dims+1)*q_prim_vf(momxe)%sf(j,k,l) + &
        !                     OmegaR(3,1)*q_prim_vf(momxb)%sf(j,k,l) + &
        !                     OmegaR(3,2)*q_prim_vf(momxb+1)%sf(j,k,l) + &
        !                     OmegaR(3,3)*q_prim_vf(momxe)%sf(j,k,l)) + &
        !                     dt/dz(l) * (sigma*gL_z(j,k,l,num_dims+1)*q_prim_vf(momxe)%sf(j,k,l) + &
        !                     OmegaL(3,1)*q_prim_vf(momxb)%sf(j,k,l) + &
        !                     OmegaL(3,2)*q_prim_vf(momxb+1)%sf(j,k,l) + &
        !                     OmegaL(3,3)*q_prim_vf(momxe)%sf(j,k,l))
        !             end do
        !         end do
        !     end do 
        ! end if

    end subroutine s_get_surface_tension

    subroutine s_compute_capillary_stress_tensors(gL, gR, j, k, l)

        integer :: j, k, l
        real(kind(0d0)), dimension(ix%beg:,iy%beg:,iz%beg:,:) :: gL, gR

        !$acc routine seq
        #:for S in ['L','R']
        Omega${S}$(1,1) = -sigma*(g${S}$(j,k,l,num_dims+1) - &
                            g${S}$(j, k, l, 1) * g${S}$(j, k, l, 1) / &
                            max(g${S}$(j, k, l, num_dims+1),1e-16))
        if (m > 0) then
            Omega${S}$(1,2) = -sigma*(g${S}$(j, k, l, 1) * g${S}$(j, k, l, 2) / &
                                max(g${S}$(j, k, l, num_dims + 1),1e-16))
            Omega${S}$(2,2) = -sigma*(g${S}$(j,k,l,num_dims+1) - &
                                g${S}$(j, k, l, 2) * g${S}$(j, k, l, 2) / &
                                max(g${S}$(j, k, l, num_dims+1),1e-16))
            if (p > 0) then
                Omega${S}$(1,3) = -sigma*(g${S}$(j, k, l, 1) * g${S}$(j, k, l, 3) / &
                                    max(g${S}$(j, k, l, num_dims + 1),1e-16))
                Omega${S}$(2,3) = -sigma*(g${S}$(j, k, l, 2) * g${S}$(j, k, l, 3) / &
                                    max(g${S}$(j, k, l, num_dims + 1),1e-16))
                Omega${S}$(3,3) = -sigma*(g${S}$(j,k,l,num_dims+1) - &
                                    g${S}$(j, k, l, 3) * g${S}$(j, k, l, 3) / &
                                    max(g${S}$(j, k, l, num_dims+1),1e-16))
            end if
        end if
        #:endfor

    end subroutine s_compute_capillary_stress_tensors

    subroutine s_finalize_surface_tension_module()

        do j = 1,num_dims
            @:DEALLOCATE(c_divs%vf(j)%sf)
        end do

        @:DEALLOCATE(c_divs%vf)

    end subroutine s_finalize_surface_tension_module

end module m_surface_tension