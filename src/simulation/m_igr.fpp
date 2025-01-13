#:include 'macros.fpp'

module m_igr


    use m_derived_types        !< Definitions of the derived types

    use m_global_parameters

    use m_variables_conversion

    use m_mpi_proxy

    implicit none

    private; public :: s_initialize_igr_module, &
        s_reconstruct_deriv, &
        s_igr_jacobi_iteration, &
        s_igr_riemann_solver, &
        s_initialize_igr, &
        s_reconstruct_prim_vars_igr, &
        s_get_viscous_igr, &
        s_igr_flux_add

    real(wp), allocatable, dimension(:, :, :) :: rho_igr, dux, duy, dvx, dvy, fd_coeff, duz, dvz, dwz, dwx, dwy
    real(wp), allocatable, dimension(:, :, :) :: jac, jac_old, rhs_igr, jac_rhs
    real(wp), allocatable, dimension(:, :, :) :: duLx, duLy, dvLx, dvLy, duLz, dvLz, dwLz, dwLx, dwLy, FL
    real(wp), allocatable, dimension(:, :, :) :: duRx, duRy, dvRx, dvRy, duRz, dvRz, dwRz, dwRx, dwRy, FR
    !$acc declare create(rho_igr, dux, duy, dvx, dvy, fd_coeff,jac, jac_old, rhs_igr, jac_rhs, duz, dvz, dwz, dwx, dwy)
    !$acc declare create(duLx, duLy, dvLx, dvLy, duLz, dvLz, dwLz, dwLx, dwLy, FL)
    !$acc declare create(duRx, duRy, dvRx, dvRy, duRz, dvRz, dwRz, dwRx, dwRy, FR)

    real(kind(0d0)), allocatable, dimension(:, :, :, :) :: qL_rs_vf, qR_rs_vf


    real(kind(0d0)) :: alf_igr, omega, mu, bcxb, bcxe, bcyb, bcye, bczb, bcze
    !$acc declare create(alf_igr, omega, mu, bcxb, bcxe, bcyb, bcye, bczb, bcze)

    type(int_bounds_info) :: ix, iy, iz
    !$acc declare create(ix, iy, iz)

    integer :: i, j, k, l, q

contains

    subroutine s_initialize_igr_module()

        bcxb = bc_x%beg; bcxe = bc_x%end; bcyb = bc_y%beg; bcye = bc_y%end; bczb = bc_z%beg; bcze = bc_z%end
        !bcxb = -1; bcxe = -1; bcyb = -1; bcye = -1; bczb = -1; bcze = -1
        !$acc update device(bcxb, bcxe, bcyb, bcye, bczb, bcze)


        #:for VAR in ['rho_igr', 'dux', 'duy', 'duz', 'dvx', 'dvy', 'dvz', &
                    & 'dwx', 'dwy', 'dwz', 'jac', 'jac_rhs', 'jac_old',     &
                    & 'fd_coeff']
                  @:ALLOCATE(${VAR}$(idwbuff(1)%beg:idwbuff(1)%end, &
                                     idwbuff(2)%beg:idwbuff(2)%end, &
                                     idwbuff(3)%beg:idwbuff(3)%end))
        #:endfor

        @:ALLOCATE(qL_rs_vf(idwbuff(1)%beg:idwbuff(1)%end, &
                            idwbuff(2)%beg:idwbuff(2)%end, &
                            idwbuff(3)%beg:idwbuff(3)%end, &
                            1:sys_size))
        @:ALLOCATE(qR_rs_vf(idwbuff(1)%beg:idwbuff(1)%end, &
                            idwbuff(2)%beg:idwbuff(2)%end, &
                            idwbuff(3)%beg:idwbuff(3)%end, &
                            1:sys_size))


        if(any(Re_size > 0)) then
            #:for VAR in ['duLx', 'duLy', 'dvLx', 'dvLy', 'duLz', 'dvLz',       &
                        & 'dwLz', 'dwLx', 'dwLy', 'FL', 'duRx', 'duRy', 'dvRx', &
                        & 'dvRy', 'duRz', 'dvRz', 'dwRz', 'dwRx', 'dwRy', 'FR']
                @:ALLOCATE(${VAR}$(idwbuff(1)%beg:idwbuff(1)%end, &
                                      idwbuff(2)%beg:idwbuff(2)%end, &
                                      idwbuff(3)%beg:idwbuff(3)%end))
            #:endfor
        else
            @:ALLOCATE(FR(idwbuff(1)%beg:idwbuff(1)%end, &
                          idwbuff(2)%beg:idwbuff(2)%end, &
                          idwbuff(3)%beg:idwbuff(3)%end))
            @:ALLOCATE(FL(idwbuff(1)%beg:idwbuff(1)%end, &
                          idwbuff(2)%beg:idwbuff(2)%end, &
                          idwbuff(3)%beg:idwbuff(3)%end))
        end if

        !$acc parallel loop collapse(3) gang vector default(present)
        do l = idwbuff(3)%beg, idwbuff(3)%end
            do k = idwbuff(2)%beg, idwbuff(2)%end
                do j = idwbuff(1)%beg, idwbuff(1)%end
                    jac(j, k, l) = 0d0
                    jac_old(j, k, l) = 0d0
                    jac_rhs(j, k, l) = 0d0
               end do
            end do
        end do

        mu = 1d0 / fluid_pp(1)%Re(1)
        !$acc update device(mu)

    end subroutine s_initialize_igr_module

    ! WENO like reconstructions using only the optimial weights as given in
    ! Balsara & Shu (1999)
    subroutine s_reconstruct_deriv(qL, qR, q_prim, idir)

        real(kind(0d0)), &
            dimension(startx:, starty:, startz:), &
            intent(INOUT) :: qL, qR, q_prim
        integer, intent(IN) :: idir

        if(idir == 1) then
            !$acc parallel loop collapse(3) gang vector default(present)
            do l = idwbuff(3)%beg + 2, idwbuff(3)%end - 2
                do k = idwbuff(2)%beg + 2, idwbuff(2)%end - 2
                    do j = idwbuff(1)%beg + 2, idwbuff(1)%end - 2
                        qL(j, k, l) = (1d0/60d0) * (-3d0 * q_prim(j-2, k, l) + &
                                                    27d0 * q_prim(j-1, k, l) + &
                                                    47d0 * q_prim(j, k, l) -   &
                                                    13d0 * q_prim(j+1, k, l) + &
                                                    2d0 * q_prim(j+2, k, l))
                        qR(j, k, l) = (1d0/60d0) * (-3d0 * q_prim(j+2, k, l) + &
                                                    27d0 * q_prim(j+1, k, l) + &
                                                    47d0 * q_prim(j, k, l) -   &
                                                    13d0 * q_prim(j-1, k, l) + &
                                                    2d0 * q_prim(j-2, k, l))
                    end do
                end do
            end do
        else if(idir == 2) then
            !$acc parallel loop collapse(3) gang vector default(present)
            do l = idwbuff(3)%beg + 2, idwbuff(3)%end - 2
                do k = idwbuff(2)%beg + 2, idwbuff(2)%end - 2
                    do j = idwbuff(1)%beg + 2, idwbuff(1)%end - 2
                        qL(j, k, l) = (1d0/60d0) * (-3d0 * q_prim(j, k-2, l) + &
                                                    27d0 * q_prim(j, k-1, l) + &
                                                    47d0 * q_prim(j, k, l) -   &
                                                    13d0 * q_prim(j, k+1, l) + &
                                                    2d0 * q_prim(j, k+2, l))
                        qR(j, k, l) = (1d0/60d0) * (-3d0 * q_prim(j, k+2, l) + &
                                                    27d0 * q_prim(j, k+1, l) + &
                                                    47d0 * q_prim(j, k, l) -   &
                                                    13d0 * q_prim(j, k-1, l) + &
                                                    2d0 * q_prim(j, k-2, l))
                    end do
                end do
            end do
        else
            !$acc parallel loop collapse(3) gang vector default(present)
            do l = idwbuff(3)%beg + 2, idwbuff(3)%end - 2
                do k = idwbuff(2)%beg + 2, idwbuff(2)%end - 2
                    do j = idwbuff(1)%beg + 2, idwbuff(1)%end - 2
                        qL(j, k, l) = (1d0/60d0) * (-3d0 * q_prim(j, k, l-2) + &
                                                    27d0 * q_prim(j, k, l-1) + &
                                                    47d0 * q_prim(j, k, l) -   &
                                                    13d0 * q_prim(j, k, l+1) + &
                                                    2d0 * q_prim(j, k, l+2))
                        qR(j, k, l) = (1d0/60d0) * (-3d0 * q_prim(j, k, l+2) + &
                                                    27d0 * q_prim(j, k, l+1) + &
                                                    47d0 * q_prim(j, k, l) -   &
                                                    13d0 * q_prim(j, k, l-1) + &
                                                    2d0 * q_prim(j, k, l-2))
                    end do
                end do
            end do
        end if

    end subroutine s_reconstruct_deriv

    subroutine s_igr_jacobi_iteration()

        real(kind(0d0)) :: rho_rx, rho_ry, rho_rz, rho_lx, rho_ly, rho_lz

        do q = 1, num_igr_iters

            !$acc parallel loop collapse(3) gang vector default(present) private(rho_lx, rho_rx, rho_ly, rho_ry, rho_lz, rho_rz)
            do l = 0, p
                do k = 0, n
                    do j = 0, m
                        rho_lx = 0.5d0 *(1d0 / rho_igr(j,k,l) + 1d0 / rho_igr(j-1,k,l))
                        rho_rx = 0.5d0 *(1d0 / rho_igr(j,k,l) + 1d0 / rho_igr(j+1,k,l))
                        rho_ly = 0.5d0 *(1d0 / rho_igr(j,k,l) + 1d0 / rho_igr(j,k-1,l))
                        rho_ry = 0.5d0 *(1d0 / rho_igr(j,k,l) + 1d0 / rho_igr(j,k+1,l))

                        if(p > 0) then
                            rho_lz = 0.5d0 *(1d0 / rho_igr(j,k,l) + 1d0 / rho_igr(j,k,l-1))
                            rho_rz = 0.5d0 *(1d0 / rho_igr(j,k,l) + 1d0 / rho_igr(j,k,l+1))
                        end if

                        jac(j, k, l) = jac_rhs(j, k, l)
                        jac(j, k, l) = jac(j, k, l) + alf_igr * (1d0 / dx(j)**2d0) * (rho_lx* jac_old(j-1,k,l) + rho_rx*jac_old(j+1,k,l))
                        jac(j, k, l) = jac(j, k, l) + alf_igr * (1d0 / dy(k)**2d0) * (rho_ly* jac_old(j,k-1,l) + rho_ry*jac_old(j,k+1,l))
                        if(p > 0) then
                            jac(j, k, l) = jac(j, k, l) + alf_igr * (1d0 / dz(l)**2d0) * (rho_lz* jac_old(j,k,l-1) + rho_rz*jac_old(j,k,l+1))
                        end if
                        jac(j, k, l) = omega * (1 / fd_coeff(j,k,l))*jac(j,k,l) + (1 - omega)*jac_old(j, k, l)
                    end do
                end do
            end do

            if(bcxb >= -1) then
                if(bcxb >= 0) then
                    call s_mpi_sendrecv_F_igr(jac, 1, -1)
                else
                    !$acc parallel loop gang vector collapse(3) default(present)
                    do l = 0, p
                        do k = 0, n
                            do j = 1, buff_size
                                jac(-j, k, l) = jac(m-j+1,k,l)
                            end do
                        end do
                    end do
                end if
            end if

            if(bcxe >= -1) then
                if(bcxe >= 0) then
                    call s_mpi_sendrecv_F_igr(jac, 1, 1)
                else
                    !$acc parallel loop collapse(3) gang vector default(present)
                    do l = 0, p
                        do k = 0, n
                            do j = 1, buff_size
                                jac(m+j, k, l) = jac(j-1,k,l)
                            end do
                        end do
                    end do
                end if
            end if

            if(bcyb >= -1) then
                if(bcyb >= 0) then
                    call s_mpi_sendrecv_F_igr(jac, 2, -1)
                else
                    !$acc parallel loop collapse(3) gang vector default(present)
                    do l = 0, p
                        do k = 1, buff_size
                            do j = idwbuff(1)%beg, idwbuff(1)%end
                                jac(j,-k,l) = jac(j,n-k+1,l)
                            end do
                        end do
                    end do
                end if
            end if

            if(bcye >= -1) then
                if(bcye >= 0) then
                    call s_mpi_sendrecv_F_igr(jac, 2, 1)
                else
                    !$acc parallel loop collapse(3) gang vector default(present)
                    do l = 0, p
                        do k = 1, buff_size
                            do j = idwbuff(1)%beg, idwbuff(1)%end
                                jac(j,n+k,l) = jac(j,k-1,l)
                            end do
                        end do
                    end do
                end if
            end if

            if(p > 0) then
                if(bczb >= -1) then
                    if(bczb >= 0) then
                        call s_mpi_sendrecv_F_igr(jac, 3, -1)
                    else
                        !$acc parallel loop collapse(3) gang vector default(present)
                        do l = 1, buff_size
                            do k = idwbuff(2)%beg, idwbuff(2)%end
                                do j = idwbuff(1)%beg, idwbuff(1)%end
                                    jac(j,k,-l) = jac(j,k,p-l+1)
                                end do
                            end do
                        end do
                    end if
                end if

                if(bcze >= -1) then
                    if(bcze >= 0) then
                        call s_mpi_sendrecv_F_igr(jac, 3, 1)
                    else
                    !$acc parallel loop gang vector collapse(3) default(present)
                        do l = 1, buff_size
                            do k = idwbuff(2)%beg, idwbuff(2)%end
                                do j = idwbuff(1)%beg, idwbuff(1)%end
                                    jac(j,k,p+l) = jac(j,k,l-1)
                                end do
                            end do
                        end do
                    end if
                end if
            end if

            !$acc parallel loop collapse(3) gang vector default(present)
            do l = idwbuff(3)%beg, idwbuff(3)%end
                do k = idwbuff(2)%beg, idwbuff(2)%end
                    do j = idwbuff(1)%beg, idwbuff(1)%end
                        jac_old(j, k, l) = jac(j, k, l)
                    end do
                end do
            end do

        end do

    end subroutine s_igr_jacobi_iteration

    subroutine s_igr_riemann_solver(flux_vf, idir)

        type(scalar_field), &
            dimension(sys_size), &
            intent(inout) :: flux_vf

        integer, intent(in) :: idir

        if (idir == 1) then

            if(p == 0) then
                !$acc parallel loop collapse(3) gang vector default(present)
                do l = 0, p
                    do k = -1, n+1
                        do j = -1, m+1

                            flux_vf(contxb)%sf(j,k,l) = &
                                0.5d0 * (qL_rs_vf(j+1,k,l, contxb) * &
                                qL_rs_vf(j+1,k,l, momxb)) + &
                                0.5d0 * (qR_rs_vf(j,k,l, contxb) * &
                                qR_rs_vf(j,k,l, momxb)) + &
                                 250d0 * (qR_rs_vf(j, k, l, contxb) - qL_rs_vf(j+1, k, l, contxb))

                            ! Momentum -> rho*u^2 + p + [[F_igr]]
                            flux_vf(momxb)%sf(j,k,l) = &
                                 0.5d0* ( qL_rs_vf(j+1,k,l,contxb) * &
                                (qL_rs_vf(j+1,k,l,momxb)**2.0) + &
                                 qL_rs_vf(j+1,k,l,E_idx) + &
                                FL(j+1, k, l) ) + &
                                 0.5d0* ( qR_rs_vf(j,k,l,contxb) * &
                                (qR_rs_vf(j,k,l,momxb)**2.0) + &
                                 qR_rs_vf(j,k,l,E_idx) + &
                                FR(j, k, l) ) + &
                                 250d0 * (qR_rs_vf(j, k, l, momxb) - qL_rs_vf(j+1, k, l, momxb))

                            flux_vf(momxb+1)%sf(j, k, l) = 0.5d0* qL_rs_vf(j+1,k,l,contxb) * &
                                qL_rs_vf(j+1,k,l,momxb)*qL_rs_vf(j+1,k,l,momxb+1) + &
                                0.5d0*  qR_rs_vf(j,k,l,contxb) * &
                                qR_rs_vf(j,k,l,momxb)*qR_rs_vf(j,k,l,momxb+1) + &
                                250d0 * (qR_rs_vf(j, k, l, momxb+1) - qL_rs_vf(j+1, k, l, momxb+1))

                             flux_vf(E_idx)%sf(j, k, l) = &
                             0.5d0 * ( qL_rs_vf(j+1,k,l,momxb) * (qL_rs_vf(j+1,k,l,E_idx)*gammas(1) + pi_infs(1) + &
                                0.5d0*qL_rs_vf(j+1, k, l,contxb) * (qL_rs_vf(j+1, k, l,momxb)**2d0 + qL_rs_vf(j+1, k, l,momxb+1)**2d0 ) + &
                               qL_rs_vf(j+1,k,l,E_idx) + FL(j+1, k, l)) ) + &
                            0.5d0 * ( qR_rs_vf(j,k,l,momxb) * (qR_rs_vf(j,k,l,E_idx)*gammas(1) + pi_infs(1) + &
                                0.5d0*qR_rs_vf(j, k, l,contxb) * (qR_rs_vf(j, k, l,momxb)**2d0 + qR_rs_vf(j, k, l,momxb+1)**2d0 ) + &
                               qR_rs_vf(j,k,l,E_idx) + FR(j, k, l)) ) + &
                               250d0 * (qR_rs_vf(j, k, l, E_idx) - qL_rs_vf(j+1, k, l, E_idx))

                            if(any(Re_size>0)) then
                                flux_vf(momxb)%sf(j, k, l) = flux_vf(momxb)%sf(j, k, l) - &
                                            0.5d0 * mu*((4d0/3d0)*duLx(j+1, k, l) - (2d0/3d0)*dvLy(j+1, k, l)) - &
                                            0.5d0 * mu*((4d0/3d0)*duRx(j, k, l) - (2d0/3d0)*dvRy(j, k, l))

                                flux_vf(momxb+1)%sf(j, k, l) = flux_vf(momxb+1)%sf(j, k, l) - &
                                           0.5d0*mu*(duLy(j+1, k, l) + dvLx(j+1, k, l))  -  &
                                           0.5d0*mu*(duRy(j, k, l) + dvRx(j, k, l))

                                flux_vf(E_idx)%sf(j, k, l) = flux_vf(E_idx)%sf(j, k, l) - &
                                    0.5d0*mu*qL_rs_vf(j+1, k, l, momxb)*((4d0/3d0)*duLx(j+1, k, l) - (2d0/3d0)*dvLy(j+1, k, l)) - &
                                    0.5d0*mu*qL_rs_vf(j+1, k, l, momxb)*(duLy(j+1, k, l) + dvLx(j+1, k, l))   - &
                                    0.5d0*mu*qR_rs_vf(j, k, l, momxb)*((4d0/3d0)*duRx(j, k, l) - (2d0/3d0)*dvRy(j, k, l)) - &
                                    0.5d0*mu*qR_rs_vf(j, k, l, momxb)*(duRy(j, k, l) + dvRx(j, k, l))
                            end if

                        end do
                    end do
                end do
            else
                !$acc parallel loop collapse(3) gang vector default(present)
                do l = -1, p+1
                    do k = -1, n+1
                        do j = -1, m+1

                            flux_vf(contxb)%sf(j,k,l) = &
                                0.5d0 * (qL_rs_vf(j+1,k,l, contxb) * &
                                qL_rs_vf(j+1,k,l, momxb)) + &
                                0.5d0 * (qR_rs_vf(j,k,l, contxb) * &
                                qR_rs_vf(j,k,l, momxb)) + &
                                 250d0 * (qR_rs_vf(j, k, l, contxb) - qL_rs_vf(j+1, k, l, contxb))

                            ! Momentum -> rho*u^2 + p + [[F_igr]]
                            flux_vf(momxb)%sf(j,k,l) = &
                                 0.5d0* ( qL_rs_vf(j+1,k,l,contxb) * &
                                (qL_rs_vf(j+1,k,l,momxb)**2.0) + &
                                 qL_rs_vf(j+1,k,l,E_idx) + &
                                FL(j+1, k, l) ) + &
                                 0.5d0* ( qR_rs_vf(j,k,l,contxb) * &
                                (qR_rs_vf(j,k,l,momxb)**2.0) + &
                                 qR_rs_vf(j,k,l,E_idx) + &
                                FR(j, k, l) ) + &
                                 250d0 * (qR_rs_vf(j, k, l, momxb) - qL_rs_vf(j+1, k, l, momxb))

                            flux_vf(momxb+1)%sf(j, k, l) = 0.5d0* qL_rs_vf(j+1,k,l,contxb) * &
                                (qL_rs_vf(j+1,k,l,momxb)*qL_rs_vf(j+1,k,l,momxb+1)) + &
                                 0.5d0* qR_rs_vf(j,k,l,contxb) * &
                                (qR_rs_vf(j,k,l,momxb)*qR_rs_vf(j,k,l,momxb+1)) + &
                                 250d0 * (qR_rs_vf(j, k, l, momxb+1) - qL_rs_vf(j+1, k, l, momxb+1))

                            flux_vf(momxb+2)%sf(j, k, l) = 0.5d0* qL_rs_vf(j+1,k,l,contxb) * &
                                (qL_rs_vf(j+1,k,l,momxb)*qL_rs_vf(j+1,k,l,momxb+2)) + &
                                 0.5d0*  qR_rs_vf(j,k,l,contxb) * &
                                (qR_rs_vf(j,k,l,momxb)*qR_rs_vf(j,k,l,momxb+2)) + &
                                 250d0 * (qR_rs_vf(j, k, l, momxb+2) - qL_rs_vf(j+1, k, l, momxb+2))

                            flux_vf(E_idx)%sf(j, k, l) = &
                             0.5d0 * ( qL_rs_vf(j+1,k,l,momxb) * (qL_rs_vf(j+1,k,l,E_idx)*gammas(1) + pi_infs(1) + &
                                0.5d0*qL_rs_vf(j+1, k, l,contxb) * (qL_rs_vf(j+1, k, l,momxb)**2d0 + qL_rs_vf(j+1, k, l,momxb+1)**2d0 + qL_rs_vf(j+1, k, l,momxb+2)**2d0) + &
                               qL_rs_vf(j+1,k,l,E_idx) + FL(j+1, k, l)) ) + &
                            0.5d0 * ( qR_rs_vf(j,k,l,momxb) * (qR_rs_vf(j,k,l,E_idx)*gammas(1) + pi_infs(1) + &
                                0.5d0*qR_rs_vf(j, k, l,contxb) * (qR_rs_vf(j, k, l,momxb)**2d0 + qR_rs_vf(j, k, l,momxb+1)**2d0 + qR_rs_vf(j, k, l,momxb+2)**2d0) + &
                               qR_rs_vf(j,k,l,E_idx) + FR(j, k, l)) ) + &
                               250d0 * (qR_rs_vf(j, k, l, E_idx) - qL_rs_vf(j+1, k, l, E_idx))

                            if(any(Re_size>0)) then
                                flux_vf(momxb+2)%sf(j, k, l) = flux_vf(momxb+2)%sf(j, k, l) - &
                                           0.5d0*mu*(duLz(j+1, k, l) + dwLx(j+1, k, l))  -  &
                                           0.5d0*mu*(duRz(j, k, l) + dwRx(j, k, l))

                                flux_vf(momxb+1)%sf(j, k, l) = flux_vf(momxb+1)%sf(j, k, l) - &
                                           0.5d0*mu*(duLy(j+1, k, l) + dvLx(j+1, k, l))  -  &
                                           0.5d0*mu*(duRy(j, k, l) + dvRx(j, k, l))

                                flux_vf(momxb)%sf(j, k, l) = flux_vf(momxb)%sf(j, k, l) - &
                                            0.5d0 * mu*((4d0/3d0)*duLx(j+1, k, l) - (2d0/3d0)*dvLy(j+1, k, l) - (2d0/3d0) * dwLz(j+1, k, l)) - &
                                            0.5d0 * mu*((4d0/3d0)*duRx(j, k, l) - (2d0/3d0)*dvRy(j, k, l) - (2d0/3d0) * dwRz(j, k, l))

                                flux_vf(E_idx)%sf(j, k, l) = flux_vf(E_idx)%sf(j, k, l) - &
                                    0.5d0*mu*qL_rs_vf(j+1, k, l, momxb)*((4d0/3d0)*duLx(j+1, k, l) - (2d0/3d0)*dvLy(j+1, k, l) - (2d0/3d0)*dwLz(j+1, k, l)) - &
                                    0.5d0*mu*qL_rs_vf(j+1, k, l, momxb)*(duLy(j+1, k, l) + dvLx(j+1, k, l))   - &
                                    0.5d0*mu*qL_rs_vf(j+1, k, l, momxb)*(duLz(j+1, k, l) + dwLx(j+1, k, l))   - &
                                    0.5d0*mu*qR_rs_vf(j, k, l, momxb)*((4d0/3d0)*duRx(j, k, l) - (2d0/3d0)*dvRy(j, k, l) - (2d0/3d0)*dwRz(j, k, l)) - &
                                    0.5d0*mu*qR_rs_vf(j, k, l, momxb)*(duRy(j, k, l) + dvRx(j, k, l)) - &
                                    0.5d0*mu*qR_rs_vf(j, k, l, momxb)*(duRz(j, k, l) + dwRx(j, k, l))
                            end if
                        end do
                    end do
                end do
            end if
        else if (idir == 2) then
            if(p == 0) then
                !$acc parallel loop collapse(3) gang vector default(present)
                do l = 0, p
                    do k = -1, m+1
                        do j = -1, n+1

                            flux_vf(contxb)%sf(j, k, l) = &
                                0.5d0 * (qL_rs_vf(j+1,k,l, contxb) * &
                                qL_rs_vf(j+1,k,l, momxb+1)) + &
                                0.5d0 * (qR_rs_vf(j,k,l, contxb) * &
                                qR_rs_vf(j,k,l, momxb+1)) + &
                                 250d0 * (qR_rs_vf(j, k, l, contxb) - qL_rs_vf(j+1, k, l, contxb))

                            flux_vf(momxb+1)%sf(j, k, l) = &
                                 0.5d0* ( qL_rs_vf(j+1,k,l,contxb) * &
                                (qL_rs_vf(j+1,k,l,momxb+1)**2.0) + &
                                 qL_rs_vf(j+1,k,l,E_idx) + &
                                FL(j+1, k, l) ) + &
                                 0.5d0* ( qR_rs_vf(j,k,l,contxb) * &
                                (qR_rs_vf(j,k,l,momxb+1)**2.0) + &
                                 qR_rs_vf(j,k,l,E_idx) + &
                                FR(j, k, l) ) + &
                                 250d0 * (qR_rs_vf(j, k, l, momxb+1) - qL_rs_vf(j+1, k, l, momxb+1))

                            flux_vf(momxb)%sf(j, k, l) = 0.5d0* qL_rs_vf(j+1,k,l,contxb) * &
                                (qL_rs_vf(j+1,k,l,momxb)*qL_rs_vf(j+1,k,l,momxb+1)) + &
                                 0.5d0* qR_rs_vf(j,k,l,contxb) * &
                                (qR_rs_vf(j,k,l,momxb)*qR_rs_vf(j,k,l,momxb+1)) + &
                                 250d0 * (qR_rs_vf(j, k, l, momxb) - qL_rs_vf(j+1, k, l, momxb))

                            flux_vf(E_idx)%sf(j, k, l) = &
                             0.5d0 * ( qL_rs_vf(j+1,k,l,momxb+1) * (qL_rs_vf(j+1,k,l,E_idx)*gammas(1) + pi_infs(1) + &
                                0.5d0*qL_rs_vf(j+1, k, l,contxb) * (qL_rs_vf(j+1, k, l,momxb)**2d0 + qL_rs_vf(j+1, k, l,momxb+1)**2d0 ) + &
                               qL_rs_vf(j+1,k,l,E_idx) + FL(j+1, k, l)) ) + &
                            0.5d0 * ( qR_rs_vf(j,k,l,momxb+1) * (qR_rs_vf(j,k,l,E_idx)*gammas(1) + pi_infs(1) + &
                                0.5d0*qR_rs_vf(j, k, l,contxb) * (qR_rs_vf(j, k, l,momxb)**2d0 + qR_rs_vf(j, k, l,momxb+1)**2d0 ) + &
                               qR_rs_vf(j,k,l,E_idx) + FR(j, k, l)) ) + &
                               250d0 * (qR_rs_vf(j, k, l, E_idx) - qL_rs_vf(j+1, k, l, E_idx))

                            if(any(Re_size>0)) then
                                flux_vf(momxb+1)%sf(j, k, l) = flux_vf(momxb+1)%sf(j, k, l) - &
                                            0.5d0 * mu*((4d0/3d0)*dvLy(j+1, k, l) - (2d0/3d0)*duLx(j+1, k, l)) - &
                                            0.5d0 * mu*((4d0/3d0)*dvRy(j, k, l) - (2d0/3d0)*duRx(j, k, l) )

                                flux_vf(momxb)%sf(j, k, l) = flux_vf(momxb)%sf(j, k, l) - &
                                           0.5d0*mu*(duLy(j+1, k, l) + dvLx(j+1, k, l))  -  &
                                           0.5d0*mu*(duRy(j, k, l) + dvRx(j, k, l))

                                flux_vf(E_idx)%sf(j, k, l) = flux_vf(E_idx)%sf(j, k, l) - &
                                    0.5d0*mu*qL_rs_vf(j+1, k, l, momxb+1)*((4d0/3d0)*dvLy(j+1, k, l) - (2d0/3d0)*duLx(j+1, k, l)) - &
                                    0.5d0*mu*qL_rs_vf(j+1, k, l, momxb+1)*(duLy(j+1, k, l) + dvLx(j+1, k, l))   - &
                                    0.5d0*mu*qR_rs_vf(j, k, l, momxb+1)*((4d0/3d0)*dvRy(j, k, l) - (2d0/3d0)*duRx(j, k, l)) - &
                                    0.5d0*mu*qR_rs_vf(j, k, l, momxb+1)*(duRy(j, k, l) + dvRx(j, k, l))
                            end if

                        end do
                    end do
                end do
            else
                !$acc parallel loop collapse(3) gang vector default(present)
                do l = -1, p+1
                    do k = -1, n+1
                         do j = -1, m+1

                            flux_vf(contxb)%sf(j, k, l) = &
                                0.5d0 * (qL_rs_vf(j,k+1,l, contxb) * &
                                qL_rs_vf(j,k+1,l, momxb+1)) + &
                                0.5d0 * (qR_rs_vf(j,k,l, contxb) * &
                                qR_rs_vf(j,k,l, momxb+1)) + &
                                 250d0 * (qR_rs_vf(j, k, l, contxb) - qL_rs_vf(j, k+1, l, contxb))

                            flux_vf(momxb+1)%sf(j, k, l) = &
                                 0.5d0* ( qL_rs_vf(j,k+1,l,contxb) * &
                                (qL_rs_vf(j,k+1,l,momxb+1)**2.0) + &
                                 qL_rs_vf(j,k+1,l,E_idx) + &
                                FL(j, k+1, l) ) + &
                                 0.5d0* ( qR_rs_vf(j,k,l,contxb) * &
                                (qR_rs_vf(j,k,l,momxb+1)**2.0) + &
                                 qR_rs_vf(j,k,l,E_idx) + &
                                FR(j, k, l) ) + &
                                 250d0 * (qR_rs_vf(j, k, l, momxb+1) - qL_rs_vf(j, k+1, l, momxb+1))

                            flux_vf(momxb)%sf(j, k, l) = 0.5d0* qL_rs_vf(j,k+1,l,contxb) * &
                                (qL_rs_vf(j,k+1,l,momxb)*qL_rs_vf(j,k+1,l,momxb+1)) + &
                                 0.5d0* qR_rs_vf(j,k,l,contxb) * &
                                (qR_rs_vf(j,k,l,momxb)*qR_rs_vf(j,k,l,momxb+1)) + &
                                 250d0 * (qR_rs_vf(j, k, l, momxb) - qL_rs_vf(j, k+1, l, momxb))

                            flux_vf(momxb+2)%sf(j, k, l) = 0.5d0* qL_rs_vf(j,k+1,l,contxb) * &
                                (qL_rs_vf(j,k+1,l,momxb+1)*qL_rs_vf(j,k+1,l,momxb+2)) + &
                                 0.5d0* qR_rs_vf(j,k,l,contxb) * &
                                (qR_rs_vf(j,k,l,momxb+1)*qR_rs_vf(j,k,l,momxb+2)) + &
                                 250d0 * (qR_rs_vf(j, k, l, momxb+2) - qL_rs_vf(j, k+1, l, momxb+2))

                            flux_vf(E_idx)%sf(j, k, l) = &
                             0.5d0 * ( qL_rs_vf(j,k+1,l,momxb+1) * (qL_rs_vf(j,k+1,l,E_idx)*gammas(1) + pi_infs(1) + &
                                0.5d0*qL_rs_vf(j, k+1, l,contxb) * (qL_rs_vf(j, k+1, l,momxb)**2d0 + qL_rs_vf(j, k+1, l,momxb+1)**2d0 + qL_rs_vf(j,k+1,l,momxb+2)**2d0) + &
                               qL_rs_vf(j,k+1,l,E_idx) + FL(j, k+1, l)) ) + &
                            0.5d0 * ( qR_rs_vf(j,k,l,momxb+1) * (qR_rs_vf(j,k,l,E_idx)*gammas(1) + pi_infs(1) + &
                                0.5d0*qR_rs_vf(j, k, l,contxb) * (qR_rs_vf(j, k, l,momxb)**2d0 + qR_rs_vf(j, k, l,momxb+1)**2d0 + qR_rs_vf(j,k,l,momxb+2)**2d0 ) + &
                               qR_rs_vf(j,k,l,E_idx) + FR(j, k, l)) ) + &
                               250d0 * (qR_rs_vf(j, k, l, E_idx) - qL_rs_vf(j, k+1, l, E_idx))

                            if(any(Re_size>0)) then
                                flux_vf(momxb+1)%sf(j, k, l) = flux_vf(momxb+1)%sf(j, k, l) - &
                                            0.5d0 * mu*((4d0/3d0)*dvLy(j, k+1, l) - (2d0/3d0)*duLx(j, k+1, l) - (2d0/3d0)*dwLz(j,k+1,l)) - &
                                            0.5d0 * mu*((4d0/3d0)*dvRy(j, k, l) - (2d0/3d0)*duRx(j, k, l) - (2d0/3d0)*dwRz(j,k,l) )

                                flux_vf(momxb)%sf(j, k, l) = flux_vf(momxb)%sf(j, k, l) - &
                                           0.5d0*mu*(duLy(j, k+1, l) + dvLx(j, k+1, l))  -  &
                                           0.5d0*mu*(duRy(j, k, l) + dvRx(j, k, l))

                                flux_vf(momxb+2)%sf(j, k, l) = flux_vf(momxb+2)%sf(j, k, l) - &
                                           0.5d0*mu*(dvLz(j, k+1, l) + dwLy(j, k+1, l))  -  &
                                           0.5d0*mu*(dvRz(j, k, l) + dwRy(j, k, l))

                                flux_vf(E_idx)%sf(j, k, l) = flux_vf(E_idx)%sf(j, k, l) - &
                                    0.5d0*mu*qL_rs_vf(j, k+1, l, momxb+1)*((4d0/3d0)*dvLy(j, k+1, l) - (2d0/3d0)*duLx(j, k+1, l) - (2d0/3d0)*dwLz(j,k+1 ,l)) - &
                                    0.5d0*mu*qL_rs_vf(j, k+1, l, momxb+1)*(duLy(j, k+1, l) + dvLx(j, k+1, l))   - &
                                    0.5d0*mu*qL_rs_vf(j, k+1, l, momxb+1)*(dwLy(j, k+1, l) + dvLz(j, k+1, l))   - &
                                    0.5d0*mu*qR_rs_vf(j, k, l, momxb+1)*((4d0/3d0)*dvRy(j, k, l) - (2d0/3d0)*duRx(j, k, l) - (2d0/3d0)*dwRz(j ,k ,l)) - &
                                    0.5d0*mu*qR_rs_vf(j, k, l, momxb+1)*(duRy(j, k, l) + dvRx(j, k, l))   - &
                                    0.5d0*mu*qR_rs_vf(j, k, l, momxb+1)*(dwRy(j, k, l) + dvRz(j, k, l))
                            end if
                        end do
                    end do
                end do
            end if
        elseif (idir == 3) then
            !$acc parallel loop collapse(3) gang vector default(present)
            do l = -1, p+1
                do k = -1, n+1
                    do j = -1, m+1

                        flux_vf(contxb)%sf(j, k, l) = &
                            0.5d0 * (qL_rs_vf(j,k,l+1, contxb) * &
                            qL_rs_vf(j,k,l+1, momxb+2)) + &
                            0.5d0 * (qR_rs_vf(j,k,l, contxb) * &
                            qR_rs_vf(j,k,l, momxb+2)) + &
                             250d0 * (qR_rs_vf(j, k, l, contxb) - qL_rs_vf(j, k, l+1, contxb))

                        flux_vf(momxb+2)%sf(j, k, l) = &
                             0.5d0* ( qL_rs_vf(j,k,l+1,contxb) * &
                            (qL_rs_vf(j,k,l+1,momxb+2)**2.0) + &
                             qL_rs_vf(j,k,l+1,E_idx) + &
                            FL(j, k, l+1) ) + &
                             0.5d0* ( qR_rs_vf(j,k,l,contxb) * &
                            (qR_rs_vf(j,k,l,momxb+2)**2.0) + &
                             qR_rs_vf(j,k,l,E_idx) + &
                            FR(j, k, l) ) + &
                             250d0 * (qR_rs_vf(j, k, l, momxb+2) - qL_rs_vf(j, k, l+1, momxb+2))

                        flux_vf(momxb)%sf(j, k, l) = 0.5d0* qL_rs_vf(j+1,k,l,contxb) * &
                            (qL_rs_vf(j,k,l+1,momxb)*qL_rs_vf(j,k,l+1,momxb+2)) + &
                             0.5d0* qR_rs_vf(j,k,l,contxb) * &
                            (qR_rs_vf(j,k,l,momxb)*qR_rs_vf(j,k,l,momxb+2)) + &
                             250d0 * (qR_rs_vf(j, k, l, momxb) - qL_rs_vf(j, k, l+1, momxb))

                        flux_vf(momxb+1)%sf(j, k, l) = 0.5d0*  qL_rs_vf(j,k,l+1,contxb) * &
                            (qL_rs_vf(j,k,l+1,momxb+1)*qL_rs_vf(j,k,l+1,momxb+2)) + &
                             0.5d0* qR_rs_vf(j,k,l,contxb) * &
                            (qR_rs_vf(j,k,l,momxb+1)*qR_rs_vf(j,k,l,momxb+2)) + &
                             250d0 * (qR_rs_vf(j, k, l, momxb+1) - qL_rs_vf(j, k, l+1, momxb+1))

                         flux_vf(E_idx)%sf(j, k, l) = &
                         0.5d0 * ( qL_rs_vf(j,k,l+1,momxb+2) * (qL_rs_vf(j,k,l+1,E_idx)*gammas(1) + pi_infs(1) + &
                            0.5d0*qL_rs_vf(j, k, l+1,contxb) * (qL_rs_vf(j, k, l+1,momxb)**2d0 + qL_rs_vf(j, k, l+1,momxb+1)**2d0 + qL_rs_vf(j,k,l+1,momxb+2)**2d0) + &
                           qL_rs_vf(j,k,l+1,E_idx) + FL(j, k, l+1)) ) + &
                        0.5d0 * ( qR_rs_vf(j,k,l,momxb+2) * (qR_rs_vf(j,k,l,E_idx)*gammas(1) + pi_infs(1) + &
                            0.5d0*qR_rs_vf(j, k, l,contxb) * (qR_rs_vf(j, k, l,momxb)**2d0 + qR_rs_vf(j, k, l,momxb+1)**2d0 + qR_rs_vf(j,k,l,momxb+2)**2d0 ) + &
                           qR_rs_vf(j,k,l,E_idx) + FR(j, k, l)) ) + &
                           250d0 * (qR_rs_vf(j, k, l, E_idx) - qL_rs_vf(j, k, l+1, E_idx))

                        if(any(Re_size>0)) then
                            flux_vf(momxb+2)%sf(j, k, l) = flux_vf(momxb+2)%sf(j, k, l) - &
                                        0.5d0 * mu*((4d0/3d0)*dwLz(j, k, l+1) - (2d0/3d0)*duLx(j, k, l+1) - (2d0/3d0)*dvLy(j,k,l+1)) - &
                                        0.5d0 * mu*((4d0/3d0)*dwRz(j, k, l) - (2d0/3d0)*duRx(j, k, l) - (2d0/3d0)*dvRy(j,k,l) )

                            flux_vf(momxb)%sf(j, k, l) = flux_vf(momxb)%sf(j, k, l) - &
                                       0.5d0*mu*(duLz(j, k, l+1) + dwLx(j, k, l+1))  -  &
                                       0.5d0*mu*(duRz(j, k, l) + dwRx(j, k, l))

                            flux_vf(momxb+1)%sf(j, k, l) = flux_vf(momxb+1)%sf(j, k, l) - &
                                       0.5d0*mu*(dvLz(j, k, l+1) + dwLy(j, k, l+1))  -  &
                                       0.5d0*mu*(dvRz(j, k, l) + dwRy(j, k, l))

                            flux_vf(E_idx)%sf(j, k, l) = flux_vf(E_idx)%sf(j, k, l) - &
                            0.5d0*mu*qL_rs_vf(j, k, l+1, momxb+2)*((4d0/3d0)*dwLz(j, k, l+1) - (2d0/3d0)*duLx(j, k, l+1) - (2d0/3d0)*dvLy(j ,k ,l+1)) - &
                            0.5d0*mu*qL_rs_vf(j, k, l+1, momxb+2)*(duLz(j, k, l+1) + dwLx(j, k, l+1))   - &
                            0.5d0*mu*qL_rs_vf(j, k, l+1, momxb+2)*(dvLz(j, k, l+1) + dwLy(j, k, l+1))   - &
                            0.5d0*mu*qR_rs_vf(j, k, l, momxb+2)*((4d0/3d0)*dwRz(j, k, l) - (2d0/3d0)*duRx(j, k, l) - (2d0/3d0)*dvRy(j ,k ,l)) - &
                            0.5d0*mu*qR_rs_vf(j, k, l, momxb+2)*(duRz(j, k, l) + dwRx(j, k, l))   - &
                            0.5d0*mu*qR_rs_vf(j, k, l, momxb+2)*(dvRz(j, k, l) + dwRy(j, k, l))
                        end if
                    end do
                end do
            end do
        end if

    end subroutine s_igr_riemann_solver

    subroutine s_initialize_igr(q_prim_vf)

        type(scalar_field), dimension(sys_size), intent(in) :: q_prim_vf

        real(kind(0d0)) :: rho_rx, rho_ry, rho_rz, rho_lx, rho_ly, rho_lz

        alf_igr = 0d0*(dx(1)**2)
        !$acc update device(alf_igr)

        omega = 1d0
        !$acc update device(omega)

        !$acc parallel loop collapse(3) gang vector default(present)
        do l = idwbuff(3)%beg, idwbuff(3)%end
            do k = idwbuff(2)%beg, idwbuff(2)%end
                do j = idwbuff(1)%beg, idwbuff(1)%end
                    rho_igr(j,k,l) = q_prim_vf(contxb)%sf(j,k,l)
                end do
            end do
        end do

        if(p == 0) then
            !$acc parallel loop collapse(3) gang vector default(present)
            do l = 0, p
                do k = idwbuff(2)%beg + 2, idwbuff(2)%end - 2
                    do j = idwbuff(1)%beg+2, idwbuff(1)%end - 2
                        dux(j,k,l) = (1/(12d0*dx(j))) * ( &
                            8d0*q_prim_vf(momxb)%sf(j+1,k,l) - &
                            8d0*q_prim_vf(momxb)%sf(j-1,k,l) + &
                            q_prim_vf(momxb)%sf(j-2,k,l) - &
                            q_prim_vf(momxb)%sf(j+2,k,l) )

                        duy(j,k,l) = (1/(12d0*dy(k))) * ( &
                            8d0*q_prim_vf(momxb)%sf(j,k+1,l) - &
                            8d0*q_prim_vf(momxb)%sf(j,k-1,l) + &
                            q_prim_vf(momxb)%sf(j,k-2,l) - &
                            q_prim_vf(momxb)%sf(j,k+2,l) )

                    end do
                end do
            end do

            !$acc parallel loop collapse(3) gang vector default(present)
            do l = 0, p
                do k = idwbuff(2)%beg + 2, idwbuff(2)%end - 2
                    do j = idwbuff(1)%beg+2, idwbuff(1)%end - 2
                        dvx(j,k,l) = (1/(12d0*dx(j))) * ( &
                            8d0*q_prim_vf(momxb+1)%sf(j+1,k,l) - &
                            8d0*q_prim_vf(momxb+1)%sf(j-1,k,l) + &
                            q_prim_vf(momxb+1)%sf(j-2,k,l) - &
                            q_prim_vf(momxb+1)%sf(j+2,k,l) )

                        dvy(j,k,l) = (1/(12d0*dy(k))) * ( &
                            8d0*q_prim_vf(momxb+1)%sf(j,k+1,l) - &
                            8d0*q_prim_vf(momxb+1)%sf(j,k-1,l) + &
                            q_prim_vf(momxb+1)%sf(j,k-2,l) - &
                            q_prim_vf(momxb+1)%sf(j,k+2,l) )

                    end do
                end do
            end do

            !$acc parallel loop collapse(3) gang vector default(present)
            do l = 0, p
                do k = idwbuff(2)%beg + 1, idwbuff(2)%end - 1
                    do j = idwbuff(1)%beg + 1, idwbuff(1)%end - 1
                        jac_rhs(j, k, l) = alf_igr*( &
                                              dux(j,k,l)*dux(j,k,l) + dvx(j,k,l) * duy(j,k,l) + &
                                              duy(j,k,l)*dvx(j,k,l) + dvy(j,k,l) * dvy(j,k,l) +  &
                                              (dux(j,k,l) + dvy(j,k,l))**2d0)
                   end do
                end do
            end do

            !$acc parallel loop collapse(3) gang vector default(present) private(rho_lx, rho_rx, rho_ly, rho_ry)
            do l = 0, p
                do k = idwbuff(2)%beg + 1, idwbuff(2)%end - 1
                    do j = idwbuff(1)%beg + 1, idwbuff(1)%end - 1
                        rho_lx = 0.5d0 *(1d0 / rho_igr(j,k,l) + 1d0 / rho_igr(j-1,k,l))
                        rho_rx = 0.5d0 *(1d0 / rho_igr(j,k,l) + 1d0 / rho_igr(j+1,k,l))
                        rho_ly = 0.5d0 *(1d0 / rho_igr(j,k,l) + 1d0 / rho_igr(j,k-1,l))
                        rho_ry = 0.5d0 *(1d0 / rho_igr(j,k,l) + 1d0 / rho_igr(j,k+1,l))

                        fd_coeff(j,k,l) = 1d0 / rho_igr(j, k, l)

                        fd_coeff(j,k,l) = fd_coeff(j,k,l) + alf_igr * ( (1d0 / dx(j)**2d0) * (rho_lx + rho_rx) +  (1d0 / dy(k)**2d0) *(rho_ly + rho_ry) )

                    end do
                end do
            end do

        !!! 3D
        else
            !$acc parallel loop collapse(3) gang vector default(present)
            do l = idwbuff(3)%beg + 2, idwbuff(3)%end - 2
                do k = idwbuff(2)%beg + 2, idwbuff(2)%end - 2
                    do j = idwbuff(1)%beg+2, idwbuff(1)%end - 2
                        dux(j,k,l) = (1/(12d0*dx(j))) * ( &
                            8d0*q_prim_vf(momxb)%sf(j+1,k,l) - &
                            8d0*q_prim_vf(momxb)%sf(j-1,k,l) + &
                            q_prim_vf(momxb)%sf(j-2,k,l) - &
                            q_prim_vf(momxb)%sf(j+2,k,l) )

                        duy(j,k,l) = (1/(12d0*dy(k))) * ( &
                            8d0*q_prim_vf(momxb)%sf(j,k+1,l) - &
                            8d0*q_prim_vf(momxb)%sf(j,k-1,l) + &
                            q_prim_vf(momxb)%sf(j,k-2,l) - &
                            q_prim_vf(momxb)%sf(j,k+2,l) )

                        duz(j,k,l) = (1/(12d0*dz(l))) * ( &
                            8d0*q_prim_vf(momxb)%sf(j,k,l+1) - &
                            8d0*q_prim_vf(momxb)%sf(j,k,l-1) + &
                            q_prim_vf(momxb)%sf(j,k,l-2) - &
                            q_prim_vf(momxb)%sf(j,k,l+2) )

                    end do
                end do
            end do

            !$acc parallel loop collapse(3) gang vector default(present)
            do l = idwbuff(3)%beg + 2, idwbuff(3)%end - 2
                do k = idwbuff(2)%beg + 2, idwbuff(2)%end - 2
                    do j = idwbuff(1)%beg+2, idwbuff(1)%end - 2
                        dvx(j,k,l) = (1/(12d0*dx(j))) * ( &
                            8d0*q_prim_vf(momxb+1)%sf(j+1,k,l) - &
                            8d0*q_prim_vf(momxb+1)%sf(j-1,k,l) + &
                            q_prim_vf(momxb+1)%sf(j-2,k,l) - &
                            q_prim_vf(momxb+1)%sf(j+2,k,l) )

                        dvy(j,k,l) = (1/(12d0*dy(k))) * ( &
                            8d0*q_prim_vf(momxb+1)%sf(j,k+1,l) - &
                            8d0*q_prim_vf(momxb+1)%sf(j,k-1,l) + &
                            q_prim_vf(momxb+1)%sf(j,k-2,l) - &
                            q_prim_vf(momxb+1)%sf(j,k+2,l) )

                        dvz(j,k,l) = (1/(12d0*dz(l))) * ( &
                            8d0*q_prim_vf(momxb+1)%sf(j,k,l+1) - &
                            8d0*q_prim_vf(momxb+1)%sf(j,k,l-1) + &
                            q_prim_vf(momxb+1)%sf(j,k,l-2) - &
                            q_prim_vf(momxb+1)%sf(j,k,l+2) )

                    end do
                end do
            end do

            !$acc parallel loop collapse(3) gang vector default(present)
            do l = idwbuff(3)%beg + 2, idwbuff(3)%end - 2
                do k = idwbuff(2)%beg + 2, idwbuff(2)%end - 2
                    do j = idwbuff(1)%beg+2, idwbuff(1)%end - 2
                        dwx(j,k,l) = (1/(12d0*dx(j))) * ( &
                            8d0*q_prim_vf(momxb+2)%sf(j+1,k,l) - &
                            8d0*q_prim_vf(momxb+2)%sf(j-1,k,l) + &
                            q_prim_vf(momxb+2)%sf(j-2,k,l) - &
                            q_prim_vf(momxb+2)%sf(j+2,k,l) )

                        dwy(j,k,l) = (1/(12d0*dy(k))) * ( &
                            8d0*q_prim_vf(momxb+2)%sf(j,k+1,l) - &
                            8d0*q_prim_vf(momxb+2)%sf(j,k-1,l) + &
                            q_prim_vf(momxb+2)%sf(j,k-2,l) - &
                            q_prim_vf(momxb+2)%sf(j,k+2,l) )

                        dwz(j,k,l) = (1/(12d0*dz(l))) * ( &
                            8d0*q_prim_vf(momxb+2)%sf(j,k,l+1) - &
                            8d0*q_prim_vf(momxb+2)%sf(j,k,l-1) + &
                            q_prim_vf(momxb+2)%sf(j,k,l-2) - &
                            q_prim_vf(momxb+2)%sf(j,k,l+2) )

                    end do
                end do
            end do

            !$acc parallel loop collapse(3) gang vector default(present)
            do l = idwbuff(3)%beg + 1, idwbuff(3)%end - 1
                do k = idwbuff(2)%beg + 1, idwbuff(2)%end - 1
                    do j = idwbuff(1)%beg + 1, idwbuff(1)%end - 1
                        jac_rhs(j, k, l) = alf_igr*((dux(j,k,l)*dux(j,k,l) + dvx(j,k,l) * duy(j,k,l) + dwx(j,k,l) * duz(j,k,l)) + &
                                              (duy(j,k,l)*dvx(j,k,l) + dvy(j,k,l) * dvy(j,k,l) + dwy(j,k,l) * dvz(j,k,l)) + &
                                              (duz(j,k,l)*dwx(j,k,l) + dvz(j,k,l) * dwy(j,k,l) + dwz(j,k,l) * dwz(j,k,l)) + &
                                              (dux(j,k,l) + dvy(j,k,l) + dwz(j,k,l))**2d0 )
                   end do
                end do
            end do

            !$acc parallel loop collapse(3) gang vector default(present) private(rho_lx, rho_rx, rho_ly, rho_ry, rho_lz, rho_rz)
            do l = idwbuff(3)%beg + 1, idwbuff(3)%end - 1
                do k = idwbuff(2)%beg + 1, idwbuff(2)%end - 1
                    do j = idwbuff(1)%beg + 1, idwbuff(1)%end - 1
                        rho_lx = 0.5d0 *(1d0 / rho_igr(j,k,l) + 1d0 / rho_igr(j-1,k,l))
                        rho_rx = 0.5d0 *(1d0 / rho_igr(j,k,l) + 1d0 / rho_igr(j+1,k,l))
                        rho_ly = 0.5d0 *(1d0 / rho_igr(j,k,l) + 1d0 / rho_igr(j,k-1,l))
                        rho_ry = 0.5d0 *(1d0 / rho_igr(j,k,l) + 1d0 / rho_igr(j,k+1,l))
                        rho_lz = 0.5d0 *(1d0 / rho_igr(j,k,l) + 1d0 / rho_igr(j,k,l-1))
                        rho_rz = 0.5d0 *(1d0 / rho_igr(j,k,l) + 1d0 / rho_igr(j,k,l+1))

                        fd_coeff(j,k,l) = 1d0 / rho_igr(j, k, l)

                        fd_coeff(j,k,l) = fd_coeff(j,k,l) + alf_igr * ( (1d0 / dx(j)**2d0) * (rho_lx + rho_rx) +  (1d0 / dy(k)**2d0) *(rho_ly + rho_ry) + (1d0 / dz(l)**2d0) * (rho_lz + rho_rz) )

                    end do
                end do
            end do
        end if

    end subroutine s_initialize_igr

    subroutine s_reconstruct_prim_vars_igr(q_prim_vf, idir)

        type(scalar_field), intent(in), dimension(sys_size) :: q_prim_vf
        integer, intent(in) :: idir

        if (idir == 1) then
            if(p == 0) then
                !$acc parallel loop collapse(4) gang vector default(present)
                do i = 1, sys_size - 1
                    do l = 0, p
                        do k = idwbuff(2)%beg + 2, idwbuff(2)%end - 2
                            do j = idwbuff(1)%beg + 2, idwbuff(1)%end - 2
                                qL_rs_vf(j, k, l, i) = (1d0/60d0) * (-3d0 * q_prim_vf(i)%sf(j-2, k, l) + 27d0 * q_prim_vf(i)%sf(j-1, k, l) +  47d0 * q_prim_vf(i)%sf(j, k, l) -13d0 * q_prim_vf(i)%sf(j+1, k, l) + 2d0 * q_prim_vf(i)%sf(j+2, k, l))
                                qR_rs_vf(j, k, l, i) = (1d0/60d0) * (-3d0 * q_prim_vf(i)%sf(j+2, k, l) + 27d0 * q_prim_vf(i)%sf(j+1, k, l) +  47d0 * q_prim_vf(i)%sf(j, k, l) -13d0 * q_prim_vf(i)%sf(j-1, k, l) + 2d0 * q_prim_vf(i)%sf(j-2, k, l))
                            end do
                        end do
                    end do
                end do
            else
                !$acc parallel loop collapse(4) gang vector default(present)
                do i = 1, sys_size - 1
                    do l = idwbuff(3)%beg + 2, idwbuff(3)%end - 2
                        do k = idwbuff(2)%beg + 2, idwbuff(2)%end - 2
                            do j = idwbuff(1)%beg + 2, idwbuff(1)%end - 2
                                qL_rs_vf(j, k, l, i) = (1d0/60d0) * (-3d0 * q_prim_vf(i)%sf(j-2, k, l) + 27d0 * q_prim_vf(i)%sf(j-1, k, l) +  47d0 * q_prim_vf(i)%sf(j, k, l) -13d0 * q_prim_vf(i)%sf(j+1, k, l) + 2d0 * q_prim_vf(i)%sf(j+2, k, l))
                                qR_rs_vf(j, k, l, i) = (1d0/60d0) * (-3d0 * q_prim_vf(i)%sf(j+2, k, l) + 27d0 * q_prim_vf(i)%sf(j+1, k, l) +  47d0 * q_prim_vf(i)%sf(j, k, l) -13d0 * q_prim_vf(i)%sf(j-1, k, l) + 2d0 * q_prim_vf(i)%sf(j-2, k, l))
                            end do
                        end do
                    end do
                end do
            end if
        else if (idir == 2) then
            if(p == 0) then
                !$acc parallel loop collapse(4) gang vector default(present)
                do i = 1, sys_size - 1
                    do l = 0, p
                        do k = idwbuff(2)%beg + 2, idwbuff(2)%end - 2
                            do j = idwbuff(1)%beg + 2, idwbuff(1)%end - 2
                                qL_rs_vf(j, k, l, i) = (1d0/60d0) * (-3d0 * q_prim_vf(i)%sf(j, k-2, l) + 27d0 * q_prim_vf(i)%sf(j, k-1, l) +  47d0 * q_prim_vf(i)%sf(j, k, l) -13d0 * q_prim_vf(i)%sf(j, k+1, l) + 2d0 * q_prim_vf(i)%sf(j, k+2, l))
                                qR_rs_vf(j, k, l, i) = (1d0/60d0) * (-3d0 * q_prim_vf(i)%sf(j, k+2, l) + 27d0 * q_prim_vf(i)%sf(j, k+1, l) +  47d0 * q_prim_vf(i)%sf(j, k, l) -13d0 * q_prim_vf(i)%sf(j, k-1, l) + 2d0 * q_prim_vf(i)%sf(j, k-2, l))
                            end do
                        end do
                    end do
                end do
            else
                !$acc parallel loop collapse(4) gang vector default(present)
                do i = 1, sys_size - 1
                    do l = idwbuff(3)%beg + 2, idwbuff(3)%end - 2
                        do k = idwbuff(2)%beg + 2, idwbuff(2)%end - 2
                            do j = idwbuff(1)%beg + 2, idwbuff(1)%end - 2
                                qL_rs_vf(j, k, l, i) = (1d0/60d0) * (-3d0 * q_prim_vf(i)%sf(j, k-2, l) + 27d0 * q_prim_vf(i)%sf(j, k-1, l) +  47d0 * q_prim_vf(i)%sf(j, k, l) -13d0 * q_prim_vf(i)%sf(j, k+1, l) + 2d0 * q_prim_vf(i)%sf(j, k+2, l))
                                qR_rs_vf(j, k, l, i) = (1d0/60d0) * (-3d0 * q_prim_vf(i)%sf(j, k+2, l) + 27d0 * q_prim_vf(i)%sf(j, k+1, l) +  47d0 * q_prim_vf(i)%sf(j, k, l) -13d0 * q_prim_vf(i)%sf(j, k-1, l) + 2d0 * q_prim_vf(i)%sf(j, k-2, l))
                            end do
                        end do
                    end do
                end do
            end if
        else if (idir == 3) then
            !$acc parallel loop collapse(4) gang vector default(present)
            do i = 1, sys_size - 1
                do l = idwbuff(3)%beg + 2, idwbuff(3)%end - 2
                    do k = idwbuff(2)%beg + 2, idwbuff(2)%end - 2
                        do j = idwbuff(1)%beg + 2, idwbuff(1)%end - 2
                            qL_rs_vf(j, k, l, i) = (1d0/60d0) * (-3d0 * q_prim_vf(i)%sf(j, k, l-2) + 27d0 * q_prim_vf(i)%sf(j, k, l-1) +  47d0 * q_prim_vf(i)%sf(j, k, l) -13d0 * q_prim_vf(i)%sf(j, k, l+1) + 2d0 * q_prim_vf(i)%sf(j, k, l+2))
                            qR_rs_vf(j, k, l, i) = (1d0/60d0) * (-3d0 * q_prim_vf(i)%sf(j, k, l+2) + 27d0 * q_prim_vf(i)%sf(j, k, l+1) +  47d0 * q_prim_vf(i)%sf(j, k, l) -13d0 * q_prim_vf(i)%sf(j, k, l-1) + 2d0 * q_prim_vf(i)%sf(j, k, l-2))
                        end do
                    end do
                end do
            end do
        end if

        call s_reconstruct_deriv(FL, FR, jac, idir)

    end subroutine s_reconstruct_prim_vars_igr

    subroutine s_get_viscous_igr(idir)

        integer, intent(in) :: idir

        if (idir == 1) then
            call s_reconstruct_deriv(duLx, duRx, dux, 1)
            call s_reconstruct_deriv(dvLx, dvRx, dvx, 1)
            call s_reconstruct_deriv(dwLx, dwRx, dwx, 1)

            call s_reconstruct_deriv(duLy, duRy, duy, 1)
            call s_reconstruct_deriv(dvLy, dvRy, dvy, 1)
            call s_reconstruct_deriv(dwLy, dwRy, dwy, 1)

            call s_reconstruct_deriv(duLz, duRz, duz, 1)
            call s_reconstruct_deriv(dvLz, dvRz, dvz, 1)
            call s_reconstruct_deriv(dwLz, dwRz, dwz, 1)
        elseif (idir == 2) then
            call s_reconstruct_deriv(duLx, duRx, dux, 2)
            call s_reconstruct_deriv(dvLx, dvRx, dvx, 2)
            call s_reconstruct_deriv(dwLx, dwRx, dwx, 2)

            call s_reconstruct_deriv(duLy, duRy, duy, 2)
            call s_reconstruct_deriv(dvLy, dvRy, dvy, 2)
            call s_reconstruct_deriv(dwLy, dwRy, dwy, 2)

            call s_reconstruct_deriv(duLz, duRz, duz, 2)
            call s_reconstruct_deriv(dvLz, dvRz, dvz, 2)
            call s_reconstruct_deriv(dwLz, dwRz, dwz, 2)
        elseif (idir == 3) then
            call s_reconstruct_deriv(duLx, duRx, dux, 3)
            call s_reconstruct_deriv(dvLx, dvRx, dvx, 3)
            call s_reconstruct_deriv(dwLx, dwRx, dwx, 3)

            call s_reconstruct_deriv(duLy, duRy, duy, 3)
            call s_reconstruct_deriv(dvLy, dvRy, dvy, 3)
            call s_reconstruct_deriv(dwLy, dwRy, dwy, 3)

            call s_reconstruct_deriv(duLz, duRz, duz, 3)
            call s_reconstruct_deriv(dvLz, dvRz, dvz, 3)
            call s_reconstruct_deriv(dwLz, dwRz, dwz, 3)
        end if

    end subroutine s_get_viscous_igr

    subroutine s_igr_flux_add(rhs_vf, flux_vf, idir)

        type(scalar_field), &
            dimension(sys_size), &
            intent(inout) :: flux_vf, rhs_vf

        integer, intent(in) :: idir

        if (idir == 1) then
            !$acc parallel loop collapse(4) gang vector default(present)
            do i = 1, sys_size - 1
                do l = 0, p
                    do k = 0, n
                        do j = 0, m
                            rhs_vf(i)%sf(j, k, l) = 1d0/dx(j)* &
                                (flux_vf(i)%sf(j - 1 , k, l) &
                                 - flux_vf(i)%sf(j , k, l))
                        end do
                    end do
                end do
            end do
        elseif (idir == 2) then
            !$acc parallel loop collapse(4) gang vector default(present)
            do i = 1, sys_size - 1
                do l = 0, p
                    do k = 0, n
                        do j = 0, m
                            rhs_vf(i)%sf(j, k, l) = &
                                rhs_vf(i)%sf(j, k, l) + 1d0/dy(k)* &
                                (flux_vf(i)%sf(k - 1, j  , l) &
                                 - flux_vf(i)%sf(k, j , l))
                        end do
                    end do
                end do
            end do
        elseif (idir == 3) then
           !$acc parallel loop collapse(4) gang vector default(present)
            do i = 1, sys_size - 1
                do l = 0, p
                    do k = 0, n
                        do j = 0, m
                            rhs_vf(i)%sf(j, k, l) = &
                                rhs_vf(i)%sf(j, k, l) + 1d0/dz(l)* &
                                (flux_vf(i)%sf(j, k, l-1) &
                                 - flux_vf(i)%sf(j, k, l ))
                        end do
                    end do
                end do
            end do
        end if

    end subroutine s_igr_flux_add

end module m_igr
