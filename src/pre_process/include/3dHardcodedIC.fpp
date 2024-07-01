#:def Hardcoded3DVariables()
! Place any declaration of intermediate variables here

real(kind(0d0)) :: rhoH, rhoL, pRef, pInt, h, lam, wl, amp, intH, alph

real(kind(0d0)) :: eps

real(kind(0d0)) :: acc, ih, pInterface, r, seed(10)

real(kind(0d0)), dimension(0:m_glb, 0:p_glb) :: ihCsv

real(kind(0d0)) :: w

integer :: i1, i2, ios
character(len=50) :: line
real(kind(0d0)) :: r1

eps = 1e-9

open(unit=10, file='perlin.txt', status='old', action='read')

do i = 0,(m_glb+1)*(p_glb + 1)

    ! Read a line from the file into the 'line' variable
    read(10, '(A)', iostat=ios) line
    if (ios /= 0) exit  ! Exit the loop if end of file or error

    ! Read two integers followed by a float from the line
    read(line, *) i1, i2, r1

    ihCsv(i1,i2) = 2e-2*r1
end do

close (10)

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

    ih = 3*lam + ihCsv(start_idx(1) + i, start_idx(3) + k)

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

    ! Put your variable assignments here
    case default
        call s_int_to_str(patch_id, iStr)
        call s_mpi_abort("Invalid hcid specified for patch "//trim(iStr))
    end select

#:enddef
