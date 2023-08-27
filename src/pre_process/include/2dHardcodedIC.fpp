#:def Hardcoded2DVariables()

    select case (patch_icpp(patch_id)%hcid)
        case(200) ! Rayleigh Taylor case 1
            !< Everything else
            if (y_cc(j) > (5d-2/(2*pi/0.2))*sin(10*pi*x_cc(i) + pi/2d0) + 0.7) then
                q_prim_vf(advxb)%sf(i, j, 0) = 1d0 - 1d-9
                q_prim_vf(advxe)%sf(i, j, 0) = 1d-9
                q_prim_vf(contxb)%sf(i, j, 0) = (1e0 - 1d-9)*3
                q_prim_vf(contxe)%sf(i, j, 0) = (1d-9)*1
            else
                q_prim_vf(advxb)%sf(i, j, 0) = 1d-9
                q_prim_vf(advxe)%sf(i, j, 0) = 1d0 - 1d-9
                q_prim_vf(contxb)%sf(i, j, 0) = (1d-9)*3
                q_prim_vf(contxe)%sf(i, j, 0) = (1d0 - 1d-9)*1
            end if

            !< Pressure
            if (y_cc(j) > 0.7) then
                q_prim_vf(E_idx)%sf(i, j, 0) = 1d5 + 3d0*9.81*(1.2 - y_cc(j))
            else
                q_prim_vf(E_idx)%sf(i, j, 0) = 1d5 + 3d0*9.81*(1.2 - 0.7) + 1d0*9.81*(0.7 - y_cc(j))
            end if
        case(201)
            !< Everything else
            if (y_cc(j) > (5d-2/(2*pi/0.2))*sin(10*pi*x_cc(i) + pi/2d0) + 0.7) then
                q_prim_vf(advxb)%sf(i, j, 0) = 1d0 - 1d-9
                q_prim_vf(advxe)%sf(i, j, 0) = 1d-9
                q_prim_vf(contxb)%sf(i, j, 0) = (1e0 - 1d-9)*19
                q_prim_vf(contxe)%sf(i, j, 0) = (1d-9)*1
            else
                q_prim_vf(advxb)%sf(i, j, 0) = 1d-9
                q_prim_vf(advxe)%sf(i, j, 0) = 1d0 - 1d-9
                q_prim_vf(contxb)%sf(i, j, 0) = (1d-9)*19
                q_prim_vf(contxe)%sf(i, j, 0) = (1d0 - 1d-9)*1
            end if

            !< Pressure
            if (y_cc(j) > 0.7) then
                q_prim_vf(E_idx)%sf(i, j, 0) = 1d5 + 19d0*9.81*(1.2 - y_cc(j))
            else
                q_prim_vf(E_idx)%sf(i, j, 0) = 1d5 + 19d0*9.81*(1.2 - 0.7) + 1d0*9.81*(0.7 - y_cc(j))
            end if 
        case default
            call s_mpi_abort("Invalid hcid")
    end select

#:enddef
