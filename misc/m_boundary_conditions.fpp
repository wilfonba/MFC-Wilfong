!>
!! @file m_perturbation.fpp
!! @brief Contains module m_perturbation

!> @brief This module contains
module m_boundary_conditions

    use m_derived_types
    use m_global_parameters

    implicit none

    type(scalar_field), dimension(:,:), allocatable :: buff_vals
    type(integer_field), dimension(:,:), allocatable :: bc_type

    real(wp) :: x_centroid, y_centroid, z_centroid
    real(wp) :: length_x, length_y, length_z
    type(bounds_info) :: x_boundary, y_boundary, z_boundary  !<

    private; public :: s_initialize_boundary_conditions_module

contains

    subroutine s_initialize_boundary_conditions_module

        allocate(buff_vals(1:num_dims,-1:1))
        allocate(bc_type(1:num_dims,-1:1))

        if (p == 0) then
            allocate(buff_vals(1,-1)%sf(1:sys_size,0:n,0))
            allocate(buff_vals(1,1)%sf(1:sys_size,0:n,0))
            allocate(buff_vals(2,-1)%sf(0:m,1:sys_size,0))
            allocate(buff_vals(2,1)%sf(0:m,1:sys_size,0))

            allocate(bc_type(1,-1)%sf(0,0:n,0))
            allocate(bc_type(1,1)%sf(0,0:n,0))
            allocate(bc_type(2,-1)%sf(0:m,0,0))
            allocate(bc_type(2,1)%sf(0:m,0,0))
        else
            allocate(buff_vals(1,-1)%sf(1:sys_size,0:n,0:p))
            allocate(buff_vals(1,1)%sf(1:sys_size,0:n,0:p))
            allocate(buff_vals(2,-1)%sf(0:m,1:sys_size,0:p))
            allocate(buff_vals(2,1)%sf(0:m,1:sys_size,0:p))
            allocate(buff_vals(3,-1)%sf(0:m,0:n,1:sys_size))
            allocate(buff_vals(3,1)%sf(0:m,0:n,1:sys_size))

            allocate(bc_type(1,-1)%sf(0,0:n,0:p))
            allocate(bc_type(1,1)%sf(0,0:n,0:p))
            allocate(bc_type(2,-1)%sf(0:m,0,0:p))
            allocate(bc_type(2,1)%sf(0:m,0,0:p))
            allocate(bc_type(3,-1)%sf(0:m,0:n,0))
            allocate(bc_type(3,1)%sf(0:m,0:n,0))
        end if

    end subroutine s_initialize_boundary_conditions_module

    subroutine s_line_segment_bc(patch_id)

        integer, intent(in) :: patch_id

        integer :: i, j, k !< Generic loop operators

        ! x_beg and x_end
        if (patch_bc(patch_id)%dir == 1) then
            y_centroid = patch_bc(patch_id)%centroid(2)
            length_y = patch_bc(patch_id)%length(2)

            y_boundary%beg = y_centroid - 0.5_wp*length_y
            y_boundary%end = y_centroid + 0.5_wp*length_y

            #:for LOC in [-1, 1]
            if (patch_bc(patch_id)%loc == ${LOC}$) then
                do i = 0, n
                    if (y_cc(i) > y_boundary%beg .and. y_cc(i) < y_boundary%end) then
                        bc_type(1,${LOC}$)%sf(0,i,0) = patch_bc(patch_id)%type

                        ! Velocities
                        do j = 1, num_dims
                            buff_vals(1,${LOC}$)%sf(momxb+j-1,i,0) = patch_bc%vel(j)
                        end do

                        ! Density and volume fraction
                        do j = 1, num_fluids
                            buff_vals(1,${LOC}$)%sf(contxb+j-1,i,0) = patch_bc%alpha_rho(j)
                            buff_vals(1,${LOC}$)%sf(advxb+j-1,i,0) = patch_bc%alpha(j)
                        end do

                        ! Pressure
                        buff_vals(1,${LOC}$)%sf(E_idx,i,0) = patch_bc%pres

                    end if
                end do
            end if
            #:endfor

        ! y_beg and y_end
        else if (patch_bc(patch_id)%dir == 2) then
            x_centroid = patch_bc(patch_id)%centroid(1)
            length_x = patch_bc(patch_id)%length(1)

            x_boundary%beg = x_centroid - 0.5_wp*length_x
            x_boundary%end = x_centroid + 0.5_wp*length_x

            #:for LOC in [-1, 1]
            if (patch_bc(patch_id)%loc == ${LOC}$) then
                do i = 0, n
                    if (x_cc(i) > x_boundary%beg .and. x_cc(i) < x_boundary%end) then
                        bc_type(2,${LOC}$)%sf(i,0,0) = patch_bc(patch_id)%type

                        ! Velocities
                        do j = 1, num_dims
                            buff_vals(2,${LOC}$)%sf(i,momxb+j-1,0) = patch_bc%vel(j)
                        end do

                        ! Density and volume fraction
                        do j = 1, num_fluids
                            buff_vals(2,${LOC}$)%sf(i,contxb+j-1,0) = patch_bc%alpha_rho(j)
                            buff_vals(2,${LOC}$)%sf(i,advxb+j-1,0) = patch_bc%alpha(j)
                        end do

                        ! Pressure
                        buff_vals(2,${LOC}$)%sf(i,E_idx,0) = patch_bc%pres

                    end if
                end do
            end if
            #:endfor
        end if

    end subroutine s_line_segment_bc

end module m_boundary_conditions
