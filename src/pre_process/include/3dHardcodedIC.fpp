#:def Hardcoded3DVariables()
    ! Place any declaration of intermediate variables here

    real(wp) :: eps

    ! Case 302 - Single M10 Jet
    real(wp) :: r, ux_th, ux_am, p_th, p_am, rho_th, rho_am, y_th, z_th, r_th, eps_smooth, x0, SS, VV

    ! Case 303 - 7 Jet
    real(wp), dimension(:), allocatable :: y_th_arr, z_th_arr

    ! Case 303 - Additionfor three jet
    real(wp) :: z_th2, r2, z_th3, r3, rcut, xcut, rn
    real(wp), dimension(0:n, 0:p) :: rcut_arr, x0_arr
    integer :: l, q, s

    integer :: pos, start, end

    character(len=1000) :: line
    character(len=25) :: value

    integer :: nn, NJet
    integer,allocatable :: seed(:)

    call random_seed(size=nn)
    allocate(seed(nn))
    seed = proc_rank + 1   ! putting arbitrary seed to all elements
    call random_seed(put=seed)
    deallocate(seed)

    ! Basic Parameters
    NJet = 19

    eps = 1e-4_wp

    SS = 6._wp/10._wp
    VV = dsin(60._wp*pi/180._wp)*SS

    ux_th = 9 ! 10, !11.0*sqrt(1.4)
    ux_am = 0.0*7.5
    p_th = 7.2_wp
    p_am = 0.4_wp
    rho_th = 2._wp
    rho_am = 1._wp
    x0 = 0.0_wp
    r_th = 0.0_wp
    eps_smooth = 0.3_wp

    open(unit=10, file="njet.txt", status="old", action="read")
    read(10,*) NJet
    close(10)

    allocate(y_th_arr(0:NJet - 1))
    allocate(z_th_arr(0:NJet - 1))

    open(unit=10, file="jets.csv", status="old", action="read")
    do q = 0, NJet - 1
        read(10, '(A)') line  ! Read a full line as a string
        start = 1

        do l = 0, 1
            end = index(line(start:), ',')  ! Find the next comma
            if (end == 0) then
                value = trim(adjustl(line(start:)))  ! Last value in the line
            else
                value = trim(adjustl(line(start:start+end-2)))  ! Extract substring
                start = start + end  ! Move to next value
            end if
            if (l == 0) then
                read(value, *) y_th_arr(q)  ! Convert string to numeric value
            else
                read(value, *) z_th_arr(q)
            end if
        end do
        y_th_arr(q) = SS*y_th_arr(q)
        z_th_arr(q) = SS*z_th_arr(q)
    end do
    close(10)

    do q = 0, p
        do l = 0, n
            rcut = 0._wp
            do s = 0, NJet - 1
                r = sqrt((y_cc(l) - y_th_arr(s))**2._wp + (z_cc(q) - z_th_arr(s))**2._wp)
                rcut = rcut + f_cut_on(r - r_th,eps_smooth)
            end do
            rcut_arr(l,q) = rcut
            call random_number(rn)
            x0_arr(l,q) = x0 !+ 0.05*(r_th + eps_smooth)*(rn - 0.5)
        end do
    end do

#:enddef

#:def Hardcoded3D()

    select case (patch_icpp(patch_id)%hcid)

    case (303) ! Multi Jet 2 Fluid

        rcut = rcut_arr(j,k)
        xcut = f_cut_on(x_cc(i) - x0_arr(j,k),eps_smooth)

        q_prim_vf(momxb)%sf(i,j,k) = (ux_th - ux_am) * rcut * xcut + ux_am
        q_prim_vf(momxb+1)%sf(i,j,k) = 0._wp
        q_prim_vf(momxe)%sf(i,j,k) = 0._wp

        q_prim_vf(advxb)%sf(i,j,k) = (1._wp - 2._wp*eps) * rcut * xcut + eps
        q_prim_vf(advxe)%sf(i,j,k) = 1._wp - q_prim_vf(advxb)%sf(i,j,k)

        q_prim_vf(contxb)%sf(i,j,k) = q_prim_vf(advxb)%sf(i,j,k)*rho_th
        q_prim_vf(contxe)%sf(i,j,k) = q_prim_vf(advxe)%sf(i,j,k)*rho_am

        q_prim_vf(E_idx)%sf(i,j,k) = (p_th - p_am) * rcut * xcut + p_am

    case (306) ! Multi Jet 1 Fluid

        rcut = rcut_arr(j,k)
        xcut = f_cut_on(x_cc(i) - x0_arr(j,k),eps_smooth)

        q_prim_vf(momxb)%sf(i,j,k) = (ux_th - ux_am) * rcut * xcut + ux_am
        q_prim_vf(momxb+1)%sf(i,j,k) = 0._wp
        q_prim_vf(momxe)%sf(i,j,k) = 0._wp

        q_prim_vf(advxb)%sf(i,j,k) = 1._wp

        q_prim_vf(contxb)%sf(i,j,k) = (rho_th - rho_am) * rcut * xcut + rho_am

        q_prim_vf(E_idx)%sf(i,j,k) = (p_th - p_am) * rcut * xcut + p_am

    case default
        call s_int_to_str(patch_id, iStr)
        call s_mpi_abort("Invalid hcid specified for patch "//trim(iStr))
    end select

#:enddef

