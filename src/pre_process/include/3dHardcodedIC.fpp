#:def Hardcoded3DVariables()
    ! Place any declaration of intermediate variables here
    

#:enddef

#:def Hardcoded3D()

    select case(patch_icpp(patch_id)%hcid)
        case(300) ! Rayleigh-Taylor Case 1

            ! Pressure
            if (y_cc(j) > 0.7) then
                q_prim_vf(E_idx)%sf(i, j, 0) = 1d5 + 3d0*9.81*(1.2 - y_cc(j))
            else
                q_prim_vf(E_idx)%sf(i, j, 0) = 1d5 + 3d0*9.81*(0.7 - y_cc(j)) + &
                    1d0*9.81*(0.7 - y_cc(j))
            end if

            ! Everything Else
            if (y_cc(j)  > (5d-2/(4*pi/0.2))**(sin(5*pi*x_cc(i) + pi/2d0) + sin(5*pi*z_cc(k) + pi/2d0) + 1.4)) then

                q_prim_vf(advxb)%sf(i, j, 0) = 1d0 - 1d-9
                q_prim_vf(advxe)%sf(i, j, 0) = 1d-9
                q_prim_vf(contxb)%sf(i, j, 0) = (1d0 - 1d-9)*3d0
                q_prim_vf(contxe)%sf(i, j, 0) = 1d-9*1d0
            else
                q_prim_vf(advxb)%sf(i, j, 0) = 1d-9
                q_prim_vf(advxe)%sf(i, j, 0) = 1d0 - 1d-9
                q_prim_vf(contxb)%sf(i, j, 0) = 1d-9*3d0
                q_prim_vf(contxe)%sf(i, j, 0) = (1d0 - 1d-9)*1d0
            end if

        case(301) ! Rayleigh-Taylor Case 2

            ! Pressure
            if (y_cc(j) > 0.7) then
                q_prim_vf(E_idx)%sf(i, j, 0) = 1d5 + 19d0*9.81*(1.2 - y_cc(j))
            else
                q_prim_vf(E_idx)%sf(i, j, 0) = 1d5 + 19d0*9.81*(0.7 - y_cc(j)) + &
                    1d0*9.81*(0.7 - y_cc(j))
            end if

            ! Everything Else
            if (y_cc(j)  > (5d-2/(4*pi/0.2))*(sin(5*pi*x_cc(i) + pi/2d0) + sin(5*pi*z_cc(k) + pi/2d0) + 1.4)) then
                q_prim_vf(advxb)%sf(i, j, 0) = 1d0 - 1d-9
                q_prim_vf(advxe)%sf(i, j, 0) = 1d-9
                q_prim_vf(contxb)%sf(i, j, 0) = (1d0 - 1d-9)*19d0
                q_prim_vf(contxe)%sf(i, j, 0) = 1d-9*1d0
            else
                q_prim_vf(advxb)%sf(i, j, 0) = 1d-9
                q_prim_vf(advxe)%sf(i, j, 0) = 1d0 - 1d-9
                q_prim_vf(contxb)%sf(i, j, 0) = 1d-9*19d0
                q_prim_vf(contxe)%sf(i, j, 0) = (1d0 - 1d-9)*1d0
            end if 

  
        case default
            call s_int_to_str(patch_id, iStr)
            call s_mpi_abort("Invalid hcid specified for patch " // trim(iStr))
    end select

#:enddef
