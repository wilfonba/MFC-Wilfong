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

    !HCID = 203
    real(kind(0d0)) :: num

#:enddef  

#:def Hardcoded2D()
    
    select case(patch_icpp(patch_id)%hcid)
        case(200) ! Rayleigh-Taylor Case 1

            ! do q200 = 0,999
            !     interfaceX200(q200) = q200*dx200
            !     interfaceY200(q200) = (5d-2/(2*pi/0.2))*sin(2.*pi*interfaceX200(q200)/0.2 + pi/2) + 0.7
            ! end do

            ! indexX200 = 0
            ! minValX200 = 1d0
            ! do q200 = 0, 999
            !     if (abs(x_cc(i) - interfaceX200(q200)) < minValX200) then
            !         minValX200 = abs(x_cc(i) - interfaceX200(q200))
            !         indexX200 = q200
            !     end if
            ! end do

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

            ! do q201 = 0,999
            !     interfaceX201(q201) = q201*dx201
            !     interfaceY201(q201) = (5d-2/(2*pi/0.2))*sin(2.*pi*interfaceX201(q201)/0.2 + pi/2) + 0.7
            ! end do

            ! indexX201 = 0
            ! minValX201 = 1d0
            ! do q201 = 0, 999
            !     if (abs(x_cc(i) - interfaceX201(q201)) < minValX201) then
            !         minValX201 = abs(x_cc(i) - interfaceX201(q201))
            !         indexX201 = q201
            !     end if
            ! end do

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
            
            num = 0.002*(sin(2*pi*x_cc(i)) + abs(0.4*cos(10*pi*x_cc(i))) - abs(0.4*sin(10*pi*x_cc(i) - pi)))

            if (y_cc(j) > num) then
                ! Volume Fractions
                q_prim_vf(advxb)%sf(i, j, 0) = 1d0 - 1d-8
                q_prim_vf(advxe)%sf(i, j, 0) = 1d-8
                q_prim_vf(contxb)%sf(i, j, 0) = (1d0 - 1d-8)*19d0
                q_prim_vf(contxe)%sf(i, j, 0) = 1d-8*1d0
            else
                q_prim_vf(advxb)%sf(i, j, 0) = 1d-8
                q_prim_vf(advxe)%sf(i, j, 0) = 1d0 - 1d-8
                q_prim_vf(contxb)%sf(i, j, 0) = 1d-8*19d0
                q_prim_vf(contxe)%sf(i, j, 0) = (1d0 - 1d-8)*1d0
            endif

        case default
            
            if (proc_rank == 0) then
                call s_int_to_str(patch_id, iStr)
                call s_mpi_abort("Invalid hcid specified for patch " // trim(iStr))
            end if

    end select

#:enddef
