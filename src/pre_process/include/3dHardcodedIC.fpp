#:def Hardcoded3DVariables()
! Place any declaration of intermediate variables here

real(kind(0d0)) :: rhoH, rhoL, pRef, pInt, h, lam, wl, amp, intH, alph

real(kind(0d0)) :: eps

real(kind(0d0)) :: acc, ih, pInterface, r, seed(10)

eps = 1e-9

seed(1:10) = (/123, 135, 531, 324, 654, 123456, 125236, 6643, 1298, 12451234 /)
#:enddef

#:def Hardcoded3D()

select case (patch_icpp(patch_id)%hcid)
case (300) ! Rayleigh-Taylor instability
    rhoH = 3
    rhoL = 1
    pRef = 1e5
    pInt = pRef
    h = 0.7
    lam = 0.2
    wl = 2*pi/lam
    amp = 0.025/wl

    intH = amp*(sin(2*pi*x_cc(i)/lam - pi/2) + sin(2*pi*z_cc(k)/lam - pi/2)) + h

    alph = 5d-1*(1 + tanh((y_cc(j) - intH)/2.5e-3))

    if (alph < eps) alph = eps
    if (alph > 1 - eps) alph = 1 - eps

    if (y_cc(j) > intH) then
        q_prim_vf(advxb)%sf(i, j, k) = alph
        q_prim_vf(advxe)%sf(i, j, k) = 1 - alph
        q_prim_vf(contxb)%sf(i, j, k) = alph*rhoH
        q_prim_vf(contxe)%sf(i, j, k) = (1 - alph)*rhoL
        q_prim_vf(E_idx)%sf(i, j, k) = pref + rhoH*9.81*(1.2 - y_cc(j))
    else
        q_prim_vf(advxb)%sf(i, j, k) = alph
        q_prim_vf(advxe)%sf(i, j, k) = 1 - alph
        q_prim_vf(contxb)%sf(i, j, k) = alph*rhoH
        q_prim_vf(contxe)%sf(i, j, k) = (1 - alph)*rhoL
        pInt = pref + rhoH*9.81*(1.2 - intH)
        q_prim_vf(E_idx)%sf(i, j, k) = pInt + rhoL*9.81*(intH - y_cc(j))
    end if
case (301) ! 3D interface shake

    lam = 0.030
    acc = 9.81

    ih = 5*lam - lam/20*(sin(2*pi/lam*z_cc(k) + pi/2) + sin((2*pi/lam)*x_cc(i)+pi/2)) + &
        0.0005*(sin(pi/3*x_cc(i)/lam) + sin(pi/5*z_cc(k)/lam) + &
                sin(pi/11*x_cc(i)/lam + sin(pi/7*z_cc(k)/lam)))
    alph = 5d-1*(1 + tanh((y_cc(j) - ih)/1d-3))

    if (alph < 1e-6) alph = 1e-6
    if (alph > 1 - 1e-6) alph = 1 - 1e-6

    if (sigma .ne. dflt_real) q_prim_vf(c_idx)%sf(i, j, k) = alph
    q_prim_vf(advxb)%sf(i, j, k) = 1 - alph
    q_prim_vf(advxe)%sf(i, j, k) = alph
    q_prim_vf(contxb)%sf(i, j, k) = (1 - alph)*5d0
    q_prim_vf(contxe)%sf(i, j, k) = alph*1d0

    if (y_cc(j)  > ih) then
        q_prim_vf(E_idx)%sf(i, j, k) = 1d5 - 1d0*acc*(y_cc(j) - ih)
    else
        q_prim_vf(E_idx)%sf(i, j, k) = 1d5 + 5d0*acc*(ih - y_cc(j))
    end if

    if ((y_cc(j) < ih + lam/20) .and. (y_cc(j) > ih - lam/20)) then
       call random_number(r)
       q_prim_vf(momxb)%sf(i, j, k) = (r - 0.5)/4
       call random_number(r)
       q_prim_vf(momxe)%sf(i, j, k) = (r - 0.5)/4
    end if

case (302) ! 3D interface wtih exponential perturbation

    lam = 0.030
    acc = 9.81

    ih = 5*lam + 0.00075*exp(sin(2*pi*x_cc(i)/lam - pi/2) + sin(2*pi*z_cc(k)/lam - pi/2)) + &
        0.000125*(sin(13*pi*x_cc(i)/lam) + sin(11*pi*z_cc(k)/lam) + sin(5*pi*x_cc(i)/lam) + sin(7*pi*y_cc(k)/lam))
    alph = 5d-1*(1 + tanh((y_cc(j) - ih)/1d-3))

    if (alph < 1e-6) alph = 1e-6
    if (alph > 1 - 1e-6) alph = 1 - 1e-6

    if (sigma .ne. dflt_real) q_prim_vf(c_idx)%sf(i, j, k) = alph
    q_prim_vf(advxb)%sf(i, j, k) = 1 - alph
    q_prim_vf(advxe)%sf(i, j, k) = alph
    q_prim_vf(contxb)%sf(i, j, k) = (1 - alph)*5d0
    q_prim_vf(contxe)%sf(i, j, k) = alph*1d0

    if (y_cc(j)  > ih) then
        q_prim_vf(E_idx)%sf(i, j, k) = 1d5 - 1d0*acc*(y_cc(j) - ih)
    else
        q_prim_vf(momxb)%sf(i, j, k) = 0.01
        q_prim_vf(E_idx)%sf(i, j, k) = 1d5 + 5d0*acc*(ih - y_cc(j))
    end if

    ! Put your variable assignments here
    case default
        call s_int_to_str(patch_id, iStr)
        call s_mpi_abort("Invalid hcid specified for patch "//trim(iStr))
    end select

#:enddef
