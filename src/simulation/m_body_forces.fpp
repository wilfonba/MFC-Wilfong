#:include 'macros.fpp'

module m_body_forces

    ! Dependencies =============================================================

    use m_derived_types        !< Definitions of the derived types

    use m_global_parameters    !< Definitions of the global parameters

    use m_variables_conversion

    use m_nvtx

    use openacc
    ! ==========================================================================

    implicit none

    private; public :: s_compute_acceleration, &
        s_compute_mixture_density, &
        s_compute_body_forces_rhs, &
        s_initialize_body_forces_module, &
        s_finalize_body_forces_module

    real(kind(0d0)), allocatable, dimension(:, :, :) :: rhoM
    real(kind(0d0)) :: rhoF

    !$acc declare create(rhoM, rhoF)

contains

    subroutine s_initialize_body_forces_module()
        
        ! Simulation is at least 2D
        if (n > 0) then
            ! Simulation is 3D
            if (p > 0) then
                @:ALLOCATE (rhoM(-buff_size:buff_size + m, &
                            -buff_size:buff_size + n, &
                            -buff_size:buff_size + p))
            ! Simulation is 2D
            else
                @:ALLOCATE (rhoM(-buff_size:buff_size + m, &
                            -buff_size:buff_size + n, &
                            0:0))
            end if
        ! Simulation is 1D
        else
            @:ALLOCATE (rhoM(-buff_size:buff_size + m, &
                        0:0, &
                        0:0))
        end if

    end subroutine s_initialize_body_forces_module

    subroutine s_compute_acceleration(t)
        
        real(kind(0d0)) :: t

        if (m > 0) then
            accel_bf(1) = 0
            if (bf_x == 1) then
                accel_bf(1) = k_x*sin(w_x*t - p_x)
            elseif (bf_x == 2) then !< analytic
                accel_bf(1) = 1
            endif
            if (n > 0) then
                accel_bf(2) = 0
                if (bf_y == 1) then
                    accel_bf(2) = k_y*sin(w_y*t - p_y)
                elseif (bf_y == 2) then
                    accel_bf(2) = 1
                end if
                if (p > 0) then
                    accel_bf(3) = 0
                    if (bf_z == 1) then
                        accel_bf(3) = k_z*sin(w_z*t - p_z)
                    elseif (bf_z == 2) then
                        accel_bf(3) = 1
                    end if
                end if
            end if
        end if

        !$acc update device(accel_bf)

    end subroutine s_compute_acceleration

    subroutine s_compute_mixture_density(q_prim_vf)

        type(scalar_field), dimension(sys_size), intent(IN) :: q_prim_vf
        integer :: i, j, k, l !< standard iterators

        !$acc parallel loop collapse(3) gang vector default(present)
        do l = 0, p
            do k = 0, n
                do j = 0, m
                    rhoM(j,k,l) = 0d0
                end do
            end do
        end do

        !$acc parallel loop collapse(4) gang vector default(present)  
        do l = 0, p
            do k = 0, n
                do j = 0, m
                    do i = 1, num_fluids
                        rhoM(j,k,l) = rhoM(j,k,l) + &
                            q_prim_vf(contxb + i - 1)%sf(j,k,l) 
                    end do
                end do
            end do
        end do

    end subroutine s_compute_mixture_density

    subroutine s_compute_body_forces_rhs(idir, q_prim_vf, rhs_vf)

        type(scalar_field), dimension(sys_size), intent(IN) :: q_prim_vf
        type(scalar_field), dimension(sys_size), intent(INOUT) :: rhs_vf
        integer, intent(IN) :: idir

        integer :: i, j, k, l, q !< Loop variables

        if (idir == 1 .and. bf_x .ne. dflt_int) then

            !$acc parallel loop collapse(3) gang vector default(present)
            do l = 0,p
                do k = 0,n
                    do j = 0,m
                        rhs_vf(momxb)%sf(j,k,l) = rhs_vf(momxb)%sf(j,k,l) + &
                            rhoM(j,k,l)*accel_bf(1)
                        rhs_vf(E_idx)%sf(j,k,l) = rhs_vf(E_idx)%sf(j,k,l) + &
                            rhoM(j,k,l)*h_prim_vf(momxb)%sf(j,k,l)*accel_bf(1)
                    end do
                end do
            end do

            ! Six equation model
            if (model_eqns == 3) then
                !$acc parallel loop collapse(4) gang vector default(present) private(rhoF)
                do l = 0,p
                    do k = 0,n
                        do j = 0,m 
                            do q = 1,num_fluids
                                rhoF = q_prim_vf(contxb + q - 1)%sf(j,k,l)/&
                                    q_prim_vf(advxb + q - 1)%sf(j,k,l)
                                rhs_vf(intxb + q - 1)%sf(j,k,l) = &
                                    rhs_vf(intxb + q - 1)%sf(j,k,l) + &
                                    rhoF*q_prim_vf(momxb)%sf(j,k,l)*accel_bf(1)
                            end do
                        end do
                    end do
                end do
            end if

        elseif (idir == 2 .and. bf_y .ne. dflt_int) then

            !$acc parallel loop collapse(3) gang vector default(present)
            do l = 0,p
                do k = 0,n
                    do j = 0,m
                        rhs_vf(momxb+1)%sf(j,k,l) = rhs_vf(momxb+1)%sf(j,k,l) + &
                            (1d3 - rhoM(j,k,l))*accel_bf(2)
                        rhs_vf(E_idx)%sf(j,k,l) = rhs_vf(E_idx)%sf(j,k,l) + &
                            rhoM(j,k,l)*q_prim_vf(momxb+1)%sf(j,k,l)*accel_bf(2)
                    end do
                end do
            end do

            ! Six equation model
            if (model_eqns == 3) then
                !$acc parallel loop collapse(4) gang vector default(present) private(rhoF)
                do l = 0,p
                    do k = 0,n
                        do j = 0,m          
                            do q = 1,num_fluids
                                rhoF = q_prim_vf(contxb + q - 1)%sf(j,k,l)/&
                                    q_prim_vf(advxb + q - 1)%sf(j,k,l)
                                rhs_vf(intxb + q - 1)%sf(j,k,l) = &
                                    rhs_vf(intxb + q - 1)%sf(j,k,l) + &
                                    rhoF*q_prim_vf(momxb+1)%sf(j,k,l)*accel_bf(2)
                            end do
                        end do
                    end do
                end do
            end if

        elseif (idir == 3 .and. bf_z .ne. dflt_int) then

            !$acc parallel loop collapse(3) gang vector default(present)
            do l = 0,p
                do k = 0,n
                    do j = 0,m
                        rhs_vf(momxe)%sf(j,k,l) = rhs_vf(momxe)%sf(j,k,l) + &
                            rhoM(j,k,l)*accel_bf(3)
                        rhs_vf(E_idx)%sf(j,k,l) = rhs_vf(E_idx)%sf(j,k,l) + &
                            rhoM(j,k,l)*q_prim_vf(momxe)%sf(j,k,l)*accel_bf(3)
                    end do
                end do
            end do
            
            ! Six equation model
            if (model_eqns == 3) then
                !$acc parallel loop collapse(4) gang vector default(present) private(rhoF)
                do l = 0,p
                    do k = 0,n
                        do j = 0,m
                            do q = 1,num_fluids
                                rhoF = q_prim_vf(contxb + q - 1)%sf(j,k,l)/&
                                    q_prim_vf(advxb + q - 1)%sf(j,k,l)
                                rhs_vf(intxb + q - 1)%sf(j,k,l) = &
                                    rhs_vf(intxb + q - 1)%sf(j,k,l) + &
                                    rhoF*q_prim_vf(momxb)%sf(j,k,l)*accel_bf(3)
                            end do
                        end do
                    end do
                end do
            end if
        end if

    end subroutine s_compute_body_forces_rhs

    subroutine s_finalize_body_forces_module()

        @:DEALLOCATE(rhoM)

    end subroutine s_finalize_body_forces_module

end module m_body_forces
