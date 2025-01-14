#:def Hardcoded3DVariables()
    ! Place any declaration of intermediate variables here

    real(wp) :: rhoH, rhoL, pRef, pInt, h, lam, wl, amp, intH, alph

    real(wp) :: eps

    ! 3D vortex collision
    real(wp) :: r1, r2, off1, off2, v1, vr, rc, xl, yl, zl, theta, rp, vrp, phi, num, den

    eps = 1e-9_wp
#:enddef

#:def Hardcoded3D()

    select case (patch_icpp(patch_id)%hcid)
    case (300) ! Rayleigh-Taylor instability
        rhoH = 3._wp
        rhoL = 1._wp
        pRef = 1.e5_wp
        pInt = pRef
        h = 0.7_wp
        lam = 0.2_wp
        wl = 2._wp*pi/lam
        amp = 0.025_wp/wl

        intH = amp*(sin(2._wp*pi*x_cc(i)/lam - pi/2._wp) + sin(2._wp*pi*z_cc(k)/lam - pi/2._wp)) + h

        alph = 5e-1_wp*(1._wp + tanh((y_cc(j) - intH)/2.5e-3_wp))

        if (alph < eps) alph = eps
        if (alph > 1._wp - eps) alph = 1._wp - eps

        if (y_cc(j) > intH) then
            q_prim_vf(advxb)%sf(i, j, k) = alph
            q_prim_vf(advxe)%sf(i, j, k) = 1._wp - alph
            q_prim_vf(contxb)%sf(i, j, k) = alph*rhoH
            q_prim_vf(contxe)%sf(i, j, k) = (1._wp - alph)*rhoL
            q_prim_vf(E_idx)%sf(i, j, k) = pref + rhoH*9.81_wp*(1.2_wp - y_cc(j))
        else
            q_prim_vf(advxb)%sf(i, j, k) = alph
            q_prim_vf(advxe)%sf(i, j, k) = 1._wp - alph
            q_prim_vf(contxb)%sf(i, j, k) = alph*rhoH
            q_prim_vf(contxe)%sf(i, j, k) = (1._wp - alph)*rhoL
            pInt = pref + rhoH*9.81_wp*(1.2_wp - intH)
            q_prim_vf(E_idx)%sf(i, j, k) = pInt + rhoL*9.81_wp*(intH - y_cc(j))
        end if

    case (301) ! (3D lung geometry in X direction, |sin(*)+sin(*)|)
        h = 0.0_wp
        lam = 1.0_wp
        amp = patch_icpp(patch_id)%a(2)
        intH = amp*abs((sin(2*pi*y_cc(j)/lam - pi/2) + sin(2*pi*z_cc(k)/lam - pi/2)) + h)
        if (x_cc(i) > intH) then
            q_prim_vf(contxb)%sf(i, j, k) = patch_icpp(1)%alpha_rho(1)
            q_prim_vf(contxe)%sf(i, j, k) = patch_icpp(1)%alpha_rho(2)
            q_prim_vf(E_idx)%sf(i, j, k) = patch_icpp(1)%pres
            q_prim_vf(advxb)%sf(i, j, k) = patch_icpp(1)%alpha(1)
            q_prim_vf(advxe)%sf(i, j, k) = patch_icpp(1)%alpha(2)
        end if

    case (302) ! 3D vertex collision
        r1 = 1._wp
        r2 = 0.25_wp
        v1 = 1
        vr = 1
        off1 = -1.1*r2
        off2 = 1.1*r2

        if ((sqrt(z_cc(k)**2._wp + y_cc(j)**2._wp) - r1)**2._wp + (x_cc(i) + off1)**2._wp < r2**2._wp) then
            xl = x_cc(i) + off1
            yl = y_cc(j)
            zl = z_cc(k)

            theta = atan2(yl,zl)
            rp = sqrt(xl**2 + (yl - r1*sin(theta))**2 + (zl - r1*cos(theta))**2)
            vrp = -(rp/r2)*vr
            if (sqrt(yl**2 + zl**2) > r1) then
                num = sqrt((yl - r1*sin(theta))**2 + (zl - r1*cos(theta))**2)
            else
                num = -sqrt((yl - r1*sin(theta))**2 + (zl - r1*cos(theta))**2)
            end if
            den = xl
            phi = atan2(den, num) - pi/2

            q_prim_vf(momxb)%sf(i,j,k) = sin(phi)*vrp - v1
            q_prim_vf(momxb+1)%sf(i,j,k) = -cos(phi)*sin(theta)*vrp
            q_prim_vf(momxe)%sf(i,j,k) = cos(phi)*cos(theta)*vrp
            q_prim_vf(e_idx)%sf(i,j,k) = 1
            q_prim_vf(advxb)%sf(i,j,k) = 1e-9
            q_prim_vf(contxb)%sf(i,j,k) = 1e-9
            q_prim_vf(advxb+1)%sf(i,j,k) = 1-eps
            q_prim_vf(contxb+1)%sf(i,j,k) = 1-eps

        elseif ((sqrt(z_cc(k)**2._wp + y_cc(j)**2._wp) - r1)**2._wp + (x_cc(i) + off2)**2._wp < r2**2._wp) then
            xl = x_cc(i) + off2
            yl = y_cc(j)
            zl = z_cc(k)

            theta = atan2(yl,zl)
            rp = sqrt(xl**2 + (yl - r1*sin(theta))**2 + (zl - r1*cos(theta))**2)
            vrp = (rp/r2)*vr
            if (sqrt(yl**2 + zl**2) > r1) then
                num = sqrt((yl - r1*sin(theta))**2 + (zl - r1*cos(theta))**2)
            else
                num = -sqrt((yl - r1*sin(theta))**2 + (zl - r1*cos(theta))**2)
            end if
            den = xl
            phi = atan2(den, num) - pi/2

            q_prim_vf(momxb)%sf(i,j,k) = sin(phi)*vrp + v1
            q_prim_vf(momxb+1)%sf(i,j,k) = -cos(phi)*sin(theta)*vrp
            q_prim_vf(momxe)%sf(i,j,k) = cos(phi)*cos(theta)*vrp
            q_prim_vf(e_idx)%sf(i,j,k) = 1
            q_prim_vf(advxb)%sf(i,j,k) = 1e-9
            q_prim_vf(contxb)%sf(i,j,k) = 1e-9
            q_prim_vf(advxb+2)%sf(i,j,k) = 1-eps
            q_prim_vf(contxb+2)%sf(i,j,k) = 1-eps
        end if


    ! Put your variable assignments here
    case default
        call s_int_to_str(patch_id, iStr)
        call s_mpi_abort("Invalid hcid specified for patch "//trim(iStr))
    end select

#:enddef
