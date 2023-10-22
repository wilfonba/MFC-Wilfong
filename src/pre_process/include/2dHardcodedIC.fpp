#:def Hardcoded2DVariables()

    ! HCID = 200
    real(kind(0d0)), dimension(0:999) :: interfaceX200, interfaceY200
    real(kind(0d0)) :: dx200 = 0.2/999d0
    integer :: indexX200
    real(kind(0d0)) :: minValX200
    integer :: q200

    !HCID = 201
    real(kind(0d0)), dimension(0:999) :: interfaceX201, interfaceY201
    real(kind(0d0)) :: dx201 = 0.2/999d0
    integer :: indexX201
    real(kind(0d0)) :: minValX201
    integer :: q201

#:enddef  

#:def Hardcoded2D()
    
    select case(patch_icpp(patch_id)%hcid)
        case(200) ! Rayleigh-Taylor Case 1

            ! Pressure
            if (y_cc(j) > 0.7) then
                q_prim_vf(E_idx)%sf(i, j, 0) = 1d5 + 3d0*9.81*(1.2 - y_cc(j))
            else
                q_prim_vf(E_idx)%sf(i, j, 0) = 1d5 + 3d0*9.81*(0.7 - y_cc(j)) + &
                    1d0*9.81*(0.7 - y_cc(j))
            end if
            
            ! Everything Else
            if (y_cc(j)  > interfaceY200(indexX200)) then
                q_prim_vf(advxb)%sf(i, j, 0) = 1d0 - 1d-8
                q_prim_vf(advxe)%sf(i, j, 0) = 1d-8
                q_prim_vf(contxb)%sf(i, j, 0) = (1d0 - 1d-8)*3d0
                q_prim_vf(contxe)%sf(i, j, 0) = 1d-8*1d0
            else
                q_prim_vf(advxb)%sf(i, j, 0) = 1d-8
                q_prim_vf(advxe)%sf(i, j, 0) = 1d0 - 1d-8
                q_prim_vf(contxb)%sf(i, j, 0) = 1d-8*3d0
                q_prim_vf(contxe)%sf(i, j, 0) = (1d0 - 1d-8)*1
            end if

        case(201) ! Rayleigh-Taylor Case 2

            ! Pressure
            if (y_cc(j) > 0.7) then
                q_prim_vf(E_idx)%sf(i, j, 0) = 1d5 + 19d0*9.81*(1.2 - y_cc(j))
            else
                q_prim_vf(E_idx)%sf(i, j, 0) = 1d5 + 19d0*9.81*(0.7 - y_cc(j)) + &
                    1d0*9.81*(0.7 - y_cc(j))
            end if
            
            ! Everything Else
            ! if (y_cc(j)  > interfaceY201(indexX201)) then
            if (y_cc(j) > (5d-2/(2*pi/0.2))*sin(10*pi*x_cc(i) + pi/2d0) + 0.7) then 
                q_prim_vf(advxb)%sf(i, j, 0) = 1d0 - 1d-8
                q_prim_vf(advxe)%sf(i, j, 0) = 1d-8
                q_prim_vf(contxb)%sf(i, j, 0) = (1d0 - 1d-8)*19d0
                q_prim_vf(contxe)%sf(i, j, 0) = 1d-8*1d0
            else
                q_prim_vf(advxb)%sf(i, j, 0) = 1d-8
                q_prim_vf(advxe)%sf(i, j, 0) = 1d0 - 1d-8
                q_prim_vf(contxb)%sf(i, j, 0) = 1d-8*19d0
                q_prim_vf(contxe)%sf(i, j, 0) = (1d0 - 1d-8)*1d0
            end if 

        case(202)  ! 2D_hardcoded_ic example case

            if (y_cc(j) <= (-x_cc(i)**3 + 1)**(1d0/3d0)) then
                ! Volume Fractions
                q_prim_vf(advxb)%sf(i, j, 0) = 1d-8
                q_prim_vf(advxe)%sf(i, j, 0) = 1d0-1d-8
                ! Denssities
                q_prim_vf(contxb)%sf(i, j, 0) = 1d-8*1000d0
                q_prim_vf(contxe)%sf(i, j, 0) = (1d0-1d-8)*1d0
                ! Pressure
                q_prim_vf(E_idx)%sf(i, j, 0) = 1000d0
            end if

        case(203) ! 2D Interface

            !print*, y_cc(j), 50d-3*(r(i) - 0.5)
            if (y_cc(j) <  1e-2*(sin(10*pi*x_cc(i))**2d0 + sin(25*pi*x_cc(i)))) then
                ! Volume Fractions
                q_prim_vf(advxb)%sf(i, j, 0) = 1d0 - 1e-8
                q_prim_vf(advxe)%sf(i, j, 0) = 1d-8
                q_prim_vf(contxb)%sf(i, j, 0) = (1d0-1e-8)*1d3
                q_prim_vf(contxe)%sf(i, j, 0) = 1d-8*1d0
            else
                q_prim_vf(advxb)%sf(i, j, 0) = 1d-8
                q_prim_vf(advxe)%sf(i, j, 0) = 1d0 - 1d-8
                q_prim_vf(contxb)%sf(i, j, 0) = 1d-8*1d3
                q_prim_vf(contxe)%sf(i, j, 0) = (1d0 - 1d-8)*1d0
            endif

        case default
            
            if (proc_rank == 0) then
                call s_int_to_str(patch_id, iStr)
                call s_mpi_abort("Invalid hcid specified for patch " // trim(iStr))
            end if

    end select

#:enddef
