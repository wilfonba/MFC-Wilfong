#:def Hardcoded3DVariables()
    ! Place any declaration of intermediate variables here
    
    ! HCID = 300
    real(kind(0d0)), dimension(0:999, 0:999) :: interfaceX300, interfaceY300, interfaceZ300
    real(kind(0d0)) :: dx300 = 0.2/999d0
    real(kind(0d0)) :: dz300 = 0.2/999d0
    integer :: indexX300, indexZ300
    real(kind(0d0)) :: minValX300, minValZ300
    integer :: q300, p300

#:enddef

#:def Hardcoded3D()

    select case(patch_icpp(patch_id)%hcid)
        case(300) ! Rayleigh-Taylor Case 1

            do q300 = 0,999
                do p300 = 0,999
                    interfaceX300(q300, p300) = q300*dx300
                    interfaceZ300(q300, p300) = p300*dz300
                    interfaceY300(q300, p300) = (5d-2/(4*pi/0.2)) * &
                        (sin(2.*pi*interfaceX300(q300, p300)/0.2 + pi/2) + &
                        sin(2.*pi*interfaceZ300(q300, p300)/0.2 + pi/2)) + 0.7
                end do
            end do

            indexX300 = 0; indexZ300 = 0
            minValX300 = 1d0; minValZ300 = 1d0
            do q300 = 0, 999
                do p300 = 0, 999
                    if (abs(x_cc(i) - interfaceX300(q300, p300)) < minValX300) then
                        minValX300 = abs(x_cc(i) - interfaceX300(q300, p300))
                        indexX300 = q300
                    end if
                    if (abs(z_cc(k) - interfaceZ300(q300, p300)) < minValZ300) then
                        minValZ300 = abs(z_cc(k) - interfaceZ300(q300, p300))
                        indexZ300 = p300
                    end if
                end do
            end do

            ! Pressure
            if (y_cc(j) > 0.7) then
                q_prim_vf(E_idx)%sf(i, j, k) = 1d5 + 3d0*9.81*(1.2 - y_cc(j))
            else
                q_prim_vf(E_idx)%sf(i, j, k) = 1d5 + 3d0*9.81*0.5 + &
                    1d0*9.81*(0.7 - y_cc(j))
            end if
 
            ! Everything Else
            if (y_cc(j)  > interfaceY300(indexX300, indexZ300)) then

                q_prim_vf(advxb)%sf(i, j, k) = 1d0 - 1d-9
                q_prim_vf(advxe)%sf(i, j, k) = 1d-9
                q_prim_vf(contxb)%sf(i, j, k) = (1d0 - 1d-9)*3d0
                q_prim_vf(contxe)%sf(i, j, k) = 1d-9*1d0
            else
                q_prim_vf(advxb)%sf(i, j, k) = 1d-9
                q_prim_vf(advxe)%sf(i, j, k) = 1d0 - 1d-9
                q_prim_vf(contxb)%sf(i, j, k) = 1d-9*3d0
                q_prim_vf(contxe)%sf(i, j, k) = (1d0 - 1d-9)*1d0
            end if

        case(301) ! Rayleigh-Taylor Case 2

            ! Pressure
            if (y_cc(j) > 0.7) then
                q_prim_vf(E_idx)%sf(i, j, k) = 1d5 + 19d0*9.81*(1.2 - y_cc(j))
            else
                q_prim_vf(E_idx)%sf(i, j, k) = 1d5 + 19d0*9.81*0.5 + &
                    1d0*9.81*(0.7 - y_cc(j))
            end if

            ! Everything Else
            if (y_cc(j)  > (5d-1/(4*pi/0.2))*(sin(10*pi*x_cc(i) + pi/2d0) + sin(10*pi*z_cc(k) + pi/2d0)) + 0.7) then
                q_prim_vf(advxb)%sf(i, j, k) = 1d0 - 1d-9
                q_prim_vf(advxe)%sf(i, j, k) = 1d-9
                q_prim_vf(contxb)%sf(i, j, k) = (1d0 - 1d-9)*19d0
                q_prim_vf(contxe)%sf(i, j, k) = 1d-9*1d0
            else
                q_prim_vf(advxb)%sf(i, j, k) = 1d-9
                q_prim_vf(advxe)%sf(i, j, k) = 1d0 - 1d-9
                q_prim_vf(contxb)%sf(i, j, k) = 1d-9*19d0
                q_prim_vf(contxe)%sf(i, j, k) = (1d0 - 1d-9)*1d0
            end if 

  
        case default
            call s_int_to_str(patch_id, iStr)
            call s_mpi_abort("Invalid hcid specified for patch " // trim(iStr))
    end select

#:enddef
