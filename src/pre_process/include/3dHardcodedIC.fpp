#:def Hardcoded3DVariables()
    ! Place any declaration of intermediate variables here
    
    ! HCID = 300
    real(kind(0d0)), dimension(0:999, 0:999) :: interfaceX300, interfaceY300, interfaceZ300
    real(kind(0d0)) :: dx300 = 0.2/999d0
    real(kind(0d0)) :: dz300 = 0.2/999d0
    integer :: indexX300, indexZ300
    real(kind(0d0)) :: minValX300, minValZ300
    integer :: q300, p300

    real(kind(0d0)) :: ih3, alph3, ih, alph, pInterface

    integer :: i1, i2
    real(kind(0d0)), dimension(0:199,0:199) :: ihPerlin

    open(32,file="/storage/scratch1/6/bwilfong3/software/MFC-Wilfong3/src/pre_process/include/perlinNoise.csv")

    do i1 = 0,199
        read(32,*) (ihPerlin(i1,i2),i2=0,199)
    end do

    print*, ihPerlin(0:2,0:2)

#:enddef

#:def Hardcoded3D()

    select case(patch_icpp(patch_id)%hcid)
        case(300) ! Rayleigh-Taylor Case 1

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

            ih3 =  (5d-2/(4*pi/0.2))*(sin(10*pi*x_cc(i) + pi/2d0) + sin(10*pi*z_cc(k) + pi/2d0)) + 0.7
            alph3 = 5d-1*(1d0 + tanh((y_cc(j) - ih3)/0.0025))

     
            ! Pressure
            if (y_cc(j) > 0.7) then
                q_prim_vf(E_idx)%sf(i, j, k) = 1d5 + 19d0*490.5*(1.2 - y_cc(j))
            else
                q_prim_vf(E_idx)%sf(i, j, k) = 1d5 + 19d0*490.5*0.5 + &
                    1d0*490.5*(0.7 - y_cc(j))
            end if

                q_prim_vf(advxb)%sf(i, j, k) = alph3
                q_prim_vf(advxe)%sf(i, j, k) = (1d0-alph3)
                q_prim_vf(contxb)%sf(i, j, k) = alph3*19d0
                q_prim_vf(contxe)%sf(i, j, k) = (1d0-alph3)*1d0
        case(302) ! 3D Shaking Interface

            ih3 = 0.5 + 0.007*(sin(pi*x_cc(i)) + sin(pi*z_cc(k)) + sin(2*pi*x_cc(i)) + sin(2*pi*z_cc(k)) + &
                sin(4*pi*x_cc(i)) + sin(4*pi*z_cc(k)) + sin(8*pi*x_cc(i)) + sin(8*pi*z_cc(k)) + &
                sin(16*pi*x_cc(i)) + sin(16*pi*z_cc(k)) + sin(32*pi*x_cc(i)) + sin(32*pi*z_cc(k)) + &
                sin(64*pi*x_cc(i)) + sin(64*pi*z_cc(k)))
            alph3 = 5d-1*(1 + tanh((y_cc(j) - ih3)/0.005))

            q_prim_vf(advxb)%sf(i, j, k) = alph3
            q_prim_vf(advxe)%sf(i, j, k) = 1d0 - alph3
            q_prim_vf(contxb)%sf(i, j, k) = alph3*2d0
            q_prim_vf(contxe)%sf(i, j, k) = (1d0-alph3)*1d0
  
        case(303) ! Perlin noise

            ih3 = 0.025*ihPerlin(i,k) + 1.5d0
            alph3 = 5d-1*(1 + tanh((y_cc(j) - ih3)/0.01))

            if (alph3 < 1e-6) alph3 = 1e-6
            if (alph3 > 1 - 1e-6) alph3 = 1 - 1e-6

            q_prim_vf(advxb)%sf(i, j, k) = 1d0 - alph3
            q_prim_vf(advxe)%sf(i, j, k) = alph3
            q_prim_vf(contxb)%sf(i, j, k) = (1d0 - alph3)*1000
            q_prim_vf(contxe)%sf(i, j, k) = alph3*1d0

            ! ! Pressure
            ! if (y_cc(j) > ih3) then
            !     q_prim_vf(E_idx)%sf(i, j, k) = 1d5 + 1d0*9.8*(ih3 - y_cc(j))
            ! else
            !     q_prim_vf(E_idx)%sf(i, j, k) = 1d5 + 1d0*9.8*(2 - ih3) + &
            !         1000d0*9.8*(ih3 - y_cc(j))
            ! end if

        case(3300) ! 3D Interface

            ih = 0.00186 - 0.00186/40*(sin(2*pi/1.86d-3*z_cc(k) + pi/2) + sin((2*pi/1.86d-3)*x_cc(i)+pi/2))
            alph = 5d-1*(1 + tanh((y_cc(j) - ih)/1d-16))

            if (alph < 1e-6) alph = 1e-6
            if (alph > 1 - 1e-6) alph = 1 - 1e-6

            if (sigma .ne. dflt_real) q_prim_vf(c_idx)%sf(i, j, 0) = alph
            q_prim_vf(advxb)%sf(i, j, k) = 1 - alph
            q_prim_vf(advxe)%sf(i, j, k) = alph
            q_prim_vf(contxb)%sf(i, j, k) = (1 - alph)*950d0
            q_prim_vf(contxe)%sf(i, j, k) = alph*1d0

            if (y_cc(j)  > 0.00186) then
                pInterface = 1d5 + 950*9.81*0.00186
                q_prim_vf(E_idx)%sf(i, j, k) = pInterface + 1d0*9.81*(y_cc(j) - 0.00186)
            else
                q_prim_vf(E_idx)%sf(i, j, k) = 1d5 + 950d0*9.81*(y_cc(j))
            end if


        case default
            call s_int_to_str(patch_id, iStr)
            call s_mpi_abort("Invalid hcid specified for patch " // trim(iStr))
    end select

#:enddef
