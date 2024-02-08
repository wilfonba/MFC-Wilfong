#:def Hardcoded2DVariables()

    ! HCID = 200
    real(kind(0d0)), dimension(0:999) :: interfaceX200, interfaceY200
    real(kind(0d0)) :: dx200 = 0.2/999d0
    integer :: indexX200
    real(kind(0d0)) :: minValX200
    integer :: q200
    real(kind(0d0)) :: pInterface

    !HCID = 201
    real(kind(0d0)), dimension(0:999) :: interfaceX201, interfaceY201
    real(kind(0d0)) :: dx201 = 0.2/999d0
    integer :: indexX201
    real(kind(0d0)) :: minValX201
    integer :: q201

    real(kind(0d0)) :: ih, alph

#:enddef  

#:def Hardcoded2D()
    
    select case(patch_icpp(patch_id)%hcid)
        case(200) ! Rayleigh-Taylor Case 1

            ih = 3 - 0.1*sin(2*pi*x_cc(i) - pi/2)

            ! Everything Else
            if (y_cc(j)  > ih) then
                q_prim_vf(advxb)%sf(i, j, 0) = 1d0 - 1d-8
                q_prim_vf(advxe)%sf(i, j, 0) = 1d-8
                q_prim_vf(contxb)%sf(i, j, 0) = (1d0 - 1d-8)*1000d0
                q_prim_vf(contxe)%sf(i, j, 0) = 1d-8*1d0
                q_prim_vf(E_idx)%sf(i, j, 0) = 1d5 + 1000*30*9.81*(4 - y_cc(j))
            else
                q_prim_vf(advxb)%sf(i, j, 0) = 1d-8
                q_prim_vf(advxe)%sf(i, j, 0) = 1d0 - 1d-8
                q_prim_vf(contxb)%sf(i, j, 0) = 1d-8*1000d0
                q_prim_vf(contxe)%sf(i, j, 0) = (1d0 - 1d-8)*1d0
                pInterface = 1d5 + 1000*30*9.81*(4 - ih)
                q_prim_vf(E_idx)%sf(i, j, 0) = pInterface + 1d0*30*9.81*(ih - y_cc(j))
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

        case(2300) ! 2D Interface

            ih = 0.016 - 0.00186/20*(sin((2*pi/1.86d-3)*x_cc(i)+pi/2))
            alph = 5d-1*(1 + tanh((y_cc(j) - ih)/1d-16))

            if (alph < 1e-6) alph = 1e-6
            if (alph > 1 - 1e-6) alph = 1 - 1e-6

            if (sigma .ne. dflt_real) q_prim_vf(c_idx)%sf(i, j, 0) = alph
            q_prim_vf(advxb)%sf(i, j, 0) = 1 - alph
            q_prim_vf(advxe)%sf(i, j, 0) = alph
            q_prim_vf(contxb)%sf(i, j, 0) = (1 - alph)*950d0
            q_prim_vf(contxe)%sf(i, j, 0) = alph*1d0

            if (y_cc(j)  > 0.016) then
                pInterface = 1d5 + 950*37.2*9.81*0.0016
                q_prim_vf(E_idx)%sf(i, j, 0) = pInterface + 1*37.2*9.81*(y_cc(j) - 0.016)
            else
                q_prim_vf(E_idx)%sf(i, j, 0) = 1d5 + 950d0*37.2*9.81*(y_cc(j))
            end if
        
        case(23002) ! 2D Interface

            ih = 0.0436 - 0.00186/20*(sin((2*pi/1.86d-3)*x_cc(i)+pi/2))
            alph = 5d-1*(1 + tanh((y_cc(j) - ih)/1d-16))

            if (alph < 1e-6) alph = 1e-6
            if (alph > 1 - 1e-6) alph = 1 - 1e-6

            if (sigma .ne. dflt_real) q_prim_vf(c_idx)%sf(i, j, 0) = alph
            q_prim_vf(advxb)%sf(i, j, 0) = 1 - alph
            q_prim_vf(advxe)%sf(i, j, 0) = alph
            q_prim_vf(contxb)%sf(i, j, 0) = (1 - alph)*950d0
            q_prim_vf(contxe)%sf(i, j, 0) = alph*1d0

            if (y_cc(j)  > 0.0436) then
                pInterface = 1d5 + 950*37.2*9.81*0.0439
                q_prim_vf(E_idx)%sf(i, j, 0) = pInterface + 1*37.2*9.81*(y_cc(j) - 0.0436)
            else
                q_prim_vf(E_idx)%sf(i, j, 0) = 1d5 + 950d0*37.2*9.81*(y_cc(j))
            end if

        case(23003) ! 2D Interface

            ih = 0.02065 - 0.00186/20*(sin((2*pi/1.86d-3)*x_cc(i)+pi/2))
            alph = 5d-1*(1 + tanh((y_cc(j) - ih)/1d-16))

            if (alph < 1e-6) alph = 1e-6
            if (alph > 1 - 1e-6) alph = 1 - 1e-6

            if (sigma .ne. dflt_real) q_prim_vf(c_idx)%sf(i, j, 0) = alph
            q_prim_vf(advxb)%sf(i, j, 0) = 1 - alph
            q_prim_vf(advxe)%sf(i, j, 0) = alph
            q_prim_vf(contxb)%sf(i, j, 0) = (1 - alph)*950d0
            q_prim_vf(contxe)%sf(i, j, 0) = alph*1d0

            if (y_cc(j)  > 0.02065) then
                pInterface = 1d5 + 950*37.2*9.81*(0.02065 + 0.13274)
                q_prim_vf(E_idx)%sf(i, j, 0) = pInterface + 1*37.2*9.81*(y_cc(j) - 0.02065)
            else
                q_prim_vf(E_idx)%sf(i, j, 0) = 1d5 + 950d0*37.2*9.81*(y_cc(j) + 0.13274)
            end if

        case(23004) ! 2D Interface

            ih = 0.0436 - 0.00186/20*(sin((2*pi/1.86d-3)*x_cc(i)+pi/2))
            alph = 5d-1*(1 + tanh((y_cc(j) - ih)/1d-16))

            if (alph < 1e-6) alph = 1e-6
            if (alph > 1 - 1e-6) alph = 1 - 1e-6

            if (sigma .ne. dflt_real) q_prim_vf(c_idx)%sf(i, j, 0) = alph
            q_prim_vf(advxb)%sf(i, j, 0) = 1 - alph
            q_prim_vf(advxe)%sf(i, j, 0) = alph
            q_prim_vf(contxb)%sf(i, j, 0) = (1 - alph)*950d0
            q_prim_vf(contxe)%sf(i, j, 0) = alph*1d0

            if (y_cc(j)  > 0.0436) then
                pInterface = 1d5 + 950*9.81*0.0439
                q_prim_vf(E_idx)%sf(i, j, 0) = pInterface + 1d0*9.81*(y_cc(j) - 0.0436)
            else
                q_prim_vf(E_idx)%sf(i, j, 0) = 1d5 + 950d0*9.81*(y_cc(j))
            end if


        case default
            
            if (proc_rank == 0) then
                call s_int_to_str(patch_id, iStr)
                call s_mpi_abort("Invalid hcid specified for patch " // trim(iStr))
            end if

    end select

#:enddef
