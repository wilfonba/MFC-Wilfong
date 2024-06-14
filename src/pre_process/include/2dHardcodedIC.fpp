#:def Hardcoded2DVariables()

    real(kind(0d0)) :: eps
    real(kind(0d0)) :: r, rmax, gam, umax, p0

    real(kind(0d0)) :: rhoH, rhoL, pRef, pInt, h, lam, wl, amp, intH, alph

    eps = 1e-9


    real(kind(0d0)), dimension(0:m_glb) :: ihCsv

    integer :: i1, ios
    character(len=50) :: line
    real(kind(0d0)) :: r1

    eps = 1e-9

    if (perlinNoise) then

        open(unit=10, file='perlin.txt', status='old', action='read')

        do i = 0,(m_glb+1)

            ! Read a line from the file into the 'line' variable
            read(10, '(A)', iostat=ios) line
            if (ios /= 0) exit  ! Exit the loop if end of file or error

            ! Read two integers followed by a float from the line
            read(line, *) i1, r1

            ihCsv(i1) = r1
        end do

    end if

#:enddef

#:def Hardcoded2D()

    select case (patch_icpp(patch_id)%hcid) ! 2D_hardcoded_ic example case
    case (200)
        if (y_cc(j) <= (-x_cc(i)**3 + 1)**(1d0/3d0)) then
            ! Volume Fractions
            q_prim_vf(advxb)%sf(i, j, 0) = eps
            q_prim_vf(advxe)%sf(i, j, 0) = 1d0 - eps
            ! Denssities
            q_prim_vf(contxb)%sf(i, j, 0) = eps*1000d0
            q_prim_vf(contxe)%sf(i, j, 0) = (1d0 - eps)*1d0
            ! Pressure
            q_prim_vf(E_idx)%sf(i, j, 0) = 1000d0
        end if
    case (202) ! Gresho vortex (Gouasmi et al 2022 JCP)
        r = ((x_cc(i) - 0.5d0)**2 + (y_cc(j) - 0.5d0)**2)**0.5d0
        rmax = 0.2

        gam = 1d0 + 1d0/fluid_pp(1)%gamma
        umax = 2*pi*rmax*patch_icpp(patch_id)%vel(2)
        p0 = umax**2*(1d0/(gam*patch_icpp(patch_id)%vel(2)**2) - 0.5d0)

        if (r < rmax) then
            q_prim_vf(momxb)%sf(i, j, 0) = -(y_cc(j) - 0.5d0)*umax/rmax
            q_prim_vf(momxe)%sf(i, j, 0) = (x_cc(i) - 0.5d0)*umax/rmax
            q_prim_vf(E_idx)%sf(i, j, 0) = p0 + umax**2*((r/rmax)**2/2d0)
        else if (r < 2*rmax) then
            q_prim_vf(momxb)%sf(i, j, 0) = -((y_cc(j) - 0.5d0)/r)*umax*(2d0 - r/rmax)
            q_prim_vf(momxe)%sf(i, j, 0) = ((x_cc(i) - 0.5d0)/r)*umax*(2d0 - r/rmax)
            q_prim_vf(E_idx)%sf(i, j, 0) = p0 + umax**2*((r/rmax)**2/2d0 + 4*(1 - (r/rmax) + log(r/rmax)))
        else
            q_prim_vf(momxb)%sf(i, j, 0) = 0d0
            q_prim_vf(momxe)%sf(i, j, 0) = 0d0
            q_prim_vf(E_idx)%sf(i, j, 0) = p0 + umax**2*(-2 + 4*log(2.))
        end if
    case (203) ! Gresho vortex (Gouasmi et al 2022 JCP) with density correction
        r = ((x_cc(i) - 0.5d0)**2 + (y_cc(j) - 0.5d0)**2)**0.5d0
        rmax = 0.2

        gam = 1d0 + 1d0/fluid_pp(1)%gamma
        umax = 2*pi*rmax*patch_icpp(patch_id)%vel(2)
        p0 = umax**2*(1d0/(gam*patch_icpp(patch_id)%vel(2)**2) - 0.5d0)

        if (r < rmax) then
            q_prim_vf(momxb)%sf(i, j, 0) = -(y_cc(j) - 0.5d0)*umax/rmax
            q_prim_vf(momxe)%sf(i, j, 0) = (x_cc(i) - 0.5d0)*umax/rmax
            q_prim_vf(E_idx)%sf(i, j, 0) = p0 + umax**2*((r/rmax)**2/2d0)
        else if (r < 2*rmax) then
            q_prim_vf(momxb)%sf(i, j, 0) = -((y_cc(j) - 0.5d0)/r)*umax*(2d0 - r/rmax)
            q_prim_vf(momxe)%sf(i, j, 0) = ((x_cc(i) - 0.5d0)/r)*umax*(2d0 - r/rmax)
            q_prim_vf(E_idx)%sf(i, j, 0) = p0 + umax**2*((r/rmax)**2/2d0 + 4*(1 - (r/rmax) + log(r/rmax)))
        else
            q_prim_vf(momxb)%sf(i, j, 0) = 0d0
            q_prim_vf(momxe)%sf(i, j, 0) = 0d0
            q_prim_vf(E_idx)%sf(i, j, 0) = p0 + umax**2*(-2 + 4*log(2.))
        end if

        q_prim_vf(contxb)%sf(i, j, 0) = q_prim_vf(E_idx)%sf(i, j, 0)**(1d0/gam)

    case (204) ! Rayleigh-Taylor instability
        rhoH = 3
        rhoL = 1
        pRef = 1e5
        pInt = pRef
        h = 0.7
        lam = 0.2
        wl = 2*pi/lam
        amp = 0.05/wl

        intH = amp*sin(2*pi*x_cc(i)/lam - pi/2) + h

        alph = 5d-1*(1 + tanh((y_cc(j) - intH)/2.5e-3))

        if (alph < eps) alph = eps
        if (alph > 1 - eps) alph = 1 - eps

        if (y_cc(j) > intH) then
            q_prim_vf(advxb)%sf(i, j, 0) = alph
            q_prim_vf(advxe)%sf(i, j, 0) = 1 - alph
            q_prim_vf(contxb)%sf(i, j, 0) = alph*rhoH
            q_prim_vf(contxe)%sf(i, j, 0) = (1 - alph)*rhoL
            q_prim_vf(E_idx)%sf(i, j, 0) = pref + rhoH*9.81*(1.2 - y_cc(j))
        else
            q_prim_vf(advxb)%sf(i, j, 0) = alph
            q_prim_vf(advxe)%sf(i, j, 0) = 1 - alph
            q_prim_vf(contxb)%sf(i, j, 0) = alph*rhoH
            q_prim_vf(contxe)%sf(i, j, 0) = (1 - alph)*rhoL
            pInt = pref + rhoH*9.81*(1.2 - intH)
            q_prim_vf(E_idx)%sf(i, j, 0) = pInt + rhoL*9.81*(intH - y_cc(j))
        end if

    case (21)

        lam = 0.030
        acc = 9.81

        ih = 5*lam + ihCsv(start_idx(1) + i

        alph = 5d-1*(1 + tanh((y_cc(j) - ih)/1d-3))

        if (alph < 1e-6) alph = 1e-6
        if (alph > 1 - 1e-6) alph = 1 - 1e-6

        if (sigma .ne. dflt_real) q_prim_vf(c_idx)%sf(i, j, k) = alph
        q_prim_vf(advxb)%sf(i, j, 0) = 1 - alph
        q_prim_vf(advxe)%sf(i, j, 0) = alph
        q_prim_vf(contxb)%sf(i, j, 0) = (1 - alph)*950d0
        q_prim_vf(contxe)%sf(i, j, 9) = alph*1d0

        if (y_cc(j)  > ih) then
            q_prim_vf(E_idx)%sf(i, j, 0) = 1d5 - 1d0*acc*(y_cc(j) - ih)
        else
            q_prim_vf(E_idx)%sf(i, j, 0) = 1d5 + 950d0*acc*(ih - y_cc(j))
        end if

    case (22)

    case (23)

    case (24)

    case (241)

    case (242)

    case (243)

    case (244)

    case (11)

    case (12)

    case (13)

    case (14)

    case (31)

    case (32)

    case (33)

    case (34)

    case default
        if (proc_rank == 0) then
            call s_int_to_str(patch_id, iStr)
            call s_mpi_abort("Invalid hcid specified for patch "//trim(iStr))
        end if
    end select

#:enddef
