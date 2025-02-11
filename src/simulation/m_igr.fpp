#:include 'macros.fpp'

module m_igr

    use m_derived_types        !< Definitions of the derived types

    use m_global_parameters

    use m_variables_conversion

    use m_mpi_proxy

    use m_helper

    implicit none

    private; public :: s_initialize_igr_module, &
        s_reconstruct_igr, &
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

    real(wp), allocatable, dimension(:, :, :, :) :: qL_rs_vf, qR_rs_vf

    real(wp) :: alf_igr, omega, mu, bcxb, bcxe, bcyb, bcye, bczb, bcze
    !$acc declare create(alf_igr, omega, mu, bcxb, bcxe, bcyb, bcye, bczb, bcze)

    type(int_bounds_info) :: ix, iy, iz
    !$acc declare create(ix, iy, iz)

    real(wp), allocatable, dimension(:, :) :: Res
    !$acc declare create(Res)

    integer :: i, j, k, l, q

contains

    subroutine s_initialize_igr_module()

        bcxb = bc_x%beg; bcxe = bc_x%end; bcyb = bc_y%beg; bcye = bc_y%end; bczb = bc_z%beg; bcze = bc_z%end
        !bcxb = -1; bcxe = -1; bcyb = -1; bcye = -1; bczb = -1; bcze = -1
        !$acc update device(bcxb, bcxe, bcyb, bcye, bczb, bcze)


        #:for VAR in ['rho_igr', 'dux', 'duy', 'dvx', 'dvy', 'jac', 'jac_rhs', &
                    & 'jac_old', 'fd_coeff']
            @:ALLOCATE(${VAR}$(idwbuff(1)%beg:idwbuff(1)%end, &
                         idwbuff(2)%beg:idwbuff(2)%end, &
                         idwbuff(3)%beg:idwbuff(3)%end))
        #:endfor

        if (p > 0) then
            #:for VAR in ['duz', 'dvz', 'dwx', 'dwy', 'dwz']
                @:ALLOCATE(${VAR}$(idwbuff(1)%beg:idwbuff(1)%end, &
                         idwbuff(2)%beg:idwbuff(2)%end, &
                         idwbuff(3)%beg:idwbuff(3)%end))

            #:endfor
        end if

        @:ALLOCATE(qL_rs_vf(idwbuff(1)%beg:idwbuff(1)%end, &
                            idwbuff(2)%beg:idwbuff(2)%end, &
                            idwbuff(3)%beg:idwbuff(3)%end, &
                            1:sys_size))
        @:ALLOCATE(qR_rs_vf(idwbuff(1)%beg:idwbuff(1)%end, &
                            idwbuff(2)%beg:idwbuff(2)%end, &
                            idwbuff(3)%beg:idwbuff(3)%end, &
                            1:sys_size))

        if(viscous) then
            #:for VAR in ['duLx', 'duLy', 'dvLx', 'dvLy', 'duRx', 'duRy', &
                        & 'dvRx', 'dvRy']
                @:ALLOCATE(${VAR}$(idwbuff(1)%beg:idwbuff(1)%end, &
                                  idwbuff(2)%beg:idwbuff(2)%end, &
                                  idwbuff(3)%beg:idwbuff(3)%end))
            #:endfor

            if (p > 0) then
                #:for VAR in ['duLz', 'dvLz', 'dwLz', 'dwLx', 'dwLy', &
                            & 'duRz', 'dvRz', 'dwRz', 'dwRx', 'dwRy']
                    @:ALLOCATE(${VAR}$(idwbuff(1)%beg:idwbuff(1)%end, &
                                      idwbuff(2)%beg:idwbuff(2)%end, &
                                      idwbuff(3)%beg:idwbuff(3)%end))
                #:endfor
            end if
        end if

        @:ALLOCATE(FR(idwbuff(1)%beg:idwbuff(1)%end, &
                      idwbuff(2)%beg:idwbuff(2)%end, &
                      idwbuff(3)%beg:idwbuff(3)%end))
        @:ALLOCATE(FL(idwbuff(1)%beg:idwbuff(1)%end, &
                      idwbuff(2)%beg:idwbuff(2)%end, &
                      idwbuff(3)%beg:idwbuff(3)%end))

        !$acc parallel loop collapse(3) gang vector default(present)
        do l = idwbuff(3)%beg, idwbuff(3)%end
            do k = idwbuff(2)%beg, idwbuff(2)%end
                do j = idwbuff(1)%beg, idwbuff(1)%end
                    jac(j, k, l) = 0._wp
                    jac_old(j, k, l) = 0._wp
                    jac_rhs(j, k, l) = 0._wp
               end do
            end do
        end do

        if (viscous) then
            @:ALLOCATE(Res(1:2, 1:maxval(Re_size)))
            do i = 1, 2
                do j = 1, Re_size(i)
                    Res(i, j) = fluid_pp(Re_idx(i, j))%Re(i)
                end do
            end do
            !$acc update device(Res, Re_idx, Re_size)
        end if

    end subroutine s_initialize_igr_module

    ! WENO like reconstructions using only the optimial weights as given in
    ! Balsara & Shu (1999)
    subroutine s_reconstruct_igr(qL, qR, q_prim, idir)

        real(wp), &
            dimension(idwbuff(1)%beg:,idwbuff(2)%beg:,idwbuff(3)%beg:), &
            intent(INOUT) :: qL, qR, q_prim
        integer, intent(IN) :: idir

        if (p == 0) then
            if(idir == 1) then
                !$acc parallel loop collapse(3) gang vector default(present)
                do l = 0, 0
                    do k = idwbuff(2)%beg + 2, idwbuff(2)%end + 2
                        do j = idwbuff(1)%beg + 2, idwbuff(1)%end + 2
                            qL(j, k, l) = (1._wp/60._wp) * (-3._wp * q_prim(j-2, k, l) + &
                                                        27._wp * q_prim(j-1, k, l) + &
                                                        47._wp * q_prim(j, k, l) -   &
                                                        13._wp * q_prim(j+1, k, l) + &
                                                        2._wp * q_prim(j+2, k, l))
                            qR(j, k, l) = (1._wp/60._wp) * (-3._wp * q_prim(j+2, k, l) + &
                                                        27._wp * q_prim(j+1, k, l) + &
                                                        47._wp * q_prim(j, k, l) -   &
                                                        13._wp * q_prim(j-1, k, l) + &
                                                        2._wp * q_prim(j-2, k, l))
                        end do
                    end do
                end do
            else if(idir == 2) then
                !$acc parallel loop collapse(3) gang vector default(present)
                do l = 0, 0
                    do k = idwbuff(2)%beg + 2, idwbuff(2)%end + 2
                        do j = idwbuff(1)%beg + 2, idwbuff(1)%end + 2
                            qL(j, k, l) = (1._wp/60._wp) * (-3._wp * q_prim(j, k-2, l) + &
                                                        27._wp * q_prim(j, k-1, l) + &
                                                        47._wp * q_prim(j, k, l) -   &
                                                        13._wp * q_prim(j, k+1, l) + &
                                                        2._wp * q_prim(j, k+2, l))
                            qR(j, k, l) = (1._wp/60._wp) * (-3._wp * q_prim(j, k+2, l) + &
                                                        27._wp * q_prim(j, k+1, l) + &
                                                        47._wp * q_prim(j, k, l) -   &
                                                        13._wp * q_prim(j, k-1, l) + &
                                                        2._wp * q_prim(j, k-2, l))
                        end do
                    end do
                end do
            end if
        else
            if(idir == 1) then
                !$acc parallel loop collapse(3) gang vector default(present)
                do l = idwbuff(3)%beg + 2, idwbuff(3)%end + 2
                    do k = idwbuff(2)%beg + 2, idwbuff(2)%end + 2
                        do j = idwbuff(1)%beg + 2, idwbuff(1)%end + 2
                            qL(j, k, l) = (1._wp/60._wp) * (-3._wp * q_prim(j-2, k, l) + &
                                                        27._wp * q_prim(j-1, k, l) + &
                                                        47._wp * q_prim(j, k, l) -   &
                                                        13._wp * q_prim(j+1, k, l) + &
                                                        2._wp * q_prim(j+2, k, l))
                            qR(j, k, l) = (1._wp/60._wp) * (-3._wp * q_prim(j+2, k, l) + &
                                                        27._wp * q_prim(j+1, k, l) + &
                                                        47._wp * q_prim(j, k, l) -   &
                                                        13._wp * q_prim(j-1, k, l) + &
                                                        2._wp * q_prim(j-2, k, l))
                        end do
                    end do
                end do
            else if(idir == 2) then
                !$acc parallel loop collapse(3) gang vector default(present)
                do l = idwbuff(3)%beg + 2, idwbuff(3)%end + 2
                    do k = idwbuff(2)%beg + 2, idwbuff(2)%end + 2
                        do j = idwbuff(1)%beg + 2, idwbuff(1)%end + 2
                            qL(j, k, l) = (1._wp/60._wp) * (-3._wp * q_prim(j, k-2, l) + &
                                                        27._wp * q_prim(j, k-1, l) + &
                                                        47._wp * q_prim(j, k, l) -   &
                                                        13._wp * q_prim(j, k+1, l) + &
                                                        2._wp * q_prim(j, k+2, l))
                            qR(j, k, l) = (1._wp/60._wp) * (-3._wp * q_prim(j, k+2, l) + &
                                                        27._wp * q_prim(j, k+1, l) + &
                                                        47._wp * q_prim(j, k, l) -   &
                                                        13._wp * q_prim(j, k-1, l) + &
                                                        2._wp * q_prim(j, k-2, l))
                        end do
                    end do
                end do
            else
                !$acc parallel loop collapse(3) gang vector default(present)
                do l = idwbuff(3)%beg + 2, idwbuff(3)%end + 2
                    do k = idwbuff(2)%beg + 2, idwbuff(2)%end + 2
                        do j = idwbuff(1)%beg + 2, idwbuff(1)%end + 2
                            qL(j, k, l) = (1._wp/60._wp) * (-3._wp * q_prim(j, k, l-2) + &
                                                        27._wp * q_prim(j, k, l-1) + &
                                                        47._wp * q_prim(j, k, l) -   &
                                                        13._wp * q_prim(j, k, l+1) + &
                                                        2._wp * q_prim(j, k, l+2))
                            qR(j, k, l) = (1._wp/60._wp) * (-3._wp * q_prim(j, k, l+2) + &
                                                        27._wp * q_prim(j, k, l+1) + &
                                                        47._wp * q_prim(j, k, l) -   &
                                                        13._wp * q_prim(j, k, l-1) + &
                                                        2._wp * q_prim(j, k, l-2))
                        end do
                    end do
                end do
            end if
        end if

    end subroutine s_reconstruct_igr

    subroutine s_igr_jacobi_iteration()

        real(wp) :: rho_rx, rho_ry, rho_rz, rho_lx, rho_ly, rho_lz

        do q = 1, num_igr_iters

            !$acc parallel loop collapse(3) gang vector default(present) private(rho_lx, rho_rx, rho_ly, rho_ry, rho_lz, rho_rz)
            do l = 0, p
                do k = 0, n
                    do j = 0, m
                        rho_lx = 0.5_wp *(1._wp / rho_igr(j,k,l) + 1._wp / rho_igr(j-1,k,l))
                        rho_rx = 0.5_wp *(1._wp / rho_igr(j,k,l) + 1._wp / rho_igr(j+1,k,l))
                        rho_ly = 0.5_wp *(1._wp / rho_igr(j,k,l) + 1._wp / rho_igr(j,k-1,l))
                        rho_ry = 0.5_wp *(1._wp / rho_igr(j,k,l) + 1._wp / rho_igr(j,k+1,l))

                        if(p > 0) then
                            rho_lz = 0.5_wp *(1._wp / rho_igr(j,k,l) + 1._wp / rho_igr(j,k,l-1))
                            rho_rz = 0.5_wp *(1._wp / rho_igr(j,k,l) + 1._wp / rho_igr(j,k,l+1))
                        end if

                        jac(j, k, l) = jac_rhs(j, k, l)
                        jac(j, k, l) = jac(j, k, l) + alf_igr * (1._wp / dx(j)**2._wp) * (rho_lx* jac_old(j-1,k,l) + rho_rx*jac_old(j+1,k,l))
                        jac(j, k, l) = jac(j, k, l) + alf_igr * (1._wp / dy(k)**2._wp) * (rho_ly* jac_old(j,k-1,l) + rho_ry*jac_old(j,k+1,l))
                        if(p > 0) then
                            jac(j, k, l) = jac(j, k, l) + alf_igr * (1._wp / dz(l)**2._wp) * (rho_lz* jac_old(j,k,l-1) + rho_rz*jac_old(j,k,l+1))
                        end if

                        jac(j, k, l) = omega * (1._wp / fd_coeff(j,k,l))*jac(j,k,l) + (1._wp - omega)*jac_old(j, k, l)
                    end do
                end do
            end do

            if(bcxb >= -12) then
                if(bcxb >= 0) then
                    call s_mpi_sendrecv_F_igr(jac, 1, -1)
                else if (bcxb == -1) then
                    !$acc parallel loop gang vector collapse(3) default(present)
                    do l = 0, p
                        do k = 0, n
                            do j = 1, buff_size
                                jac(-j, k, l) = jac(m-j+1,k,l)
                            end do
                        end do
                    end do
                else if (bcxb == -2) then
                    !$acc parallel loop gang vector collapse(3) default(present)
                    do l = 0, p
                        do k = 0, n
                            do j = 1, buff_size
                                jac(-j, k, l) = jac(j - 1,k,l)
                            end do
                        end do
                    end do
                else
                    !$acc parallel loop gang vector collapse(3) default(present)
                    do l = 0, p
                        do k = 0, n
                            do j = 1, buff_size
                                jac(-j, k, l) = jac(0,k,l)
                            end do
                        end do
                    end do

                end if
            end if

            if(bcxe >= -12) then
                if(bcxe >= 0) then
                    call s_mpi_sendrecv_F_igr(jac, 1, 1)
                else if (bcxe == -1) then
                    !$acc parallel loop collapse(3) gang vector default(present)
                    do l = 0, p
                        do k = 0, n
                            do j = 1, buff_size
                                jac(m+j, k, l) = jac(j-1,k,l)
                            end do
                        end do
                    end do
                else if (bcxe == -2) then
                    !$acc parallel loop collapse(3) gang vector default(present)
                    do l = 0, p
                        do k = 0, n
                            do j = 1, buff_size
                                jac(m+j, k, l) = jac(m - (j - 1),k,l)
                            end do
                        end do
                    end do
                else
                    !$acc parallel loop collapse(3) gang vector default(present)
                    do l = 0, p
                        do k = 0, n
                            do j = 1, buff_size
                                jac(m+j, k, l) = jac(m,k,l)
                            end do
                        end do
                    end do
                end if
            end if

            if(bcyb >= -12) then
                if(bcyb >= 0) then
                    call s_mpi_sendrecv_F_igr(jac, 2, -1)
                else if (bcyb == -1) then
                    !$acc parallel loop collapse(3) gang vector default(present)
                    do l = 0, p
                        do k = 1, buff_size
                            do j = idwbuff(1)%beg, idwbuff(1)%end
                                jac(j,-k,l) = jac(j,n-k+1,l)
                            end do
                        end do
                    end do
                else if (bcyb == -2) then
                    !$acc parallel loop collapse(3) gang vector default(present)
                    do l = 0, p
                        do k = 1, buff_size
                            do j = idwbuff(1)%beg, idwbuff(1)%end
                                jac(j,-k,l) = jac(j,k-1,l)
                            end do
                        end do
                    end do
                else
                    !$acc parallel loop collapse(3) gang vector default(present)
                    do l = 0, p
                        do k = 1, buff_size
                            do j = idwbuff(1)%beg, idwbuff(1)%end
                                jac(j,-k,l) = jac(j,0,l)
                            end do
                        end do
                    end do
                end if
            end if

            if(bcye >= -12) then
                if(bcye >= 0) then
                    call s_mpi_sendrecv_F_igr(jac, 2, 1)
                else if (bcye == -1) then
                    !$acc parallel loop collapse(3) gang vector default(present)
                    do l = 0, p
                        do k = 1, buff_size
                            do j = idwbuff(1)%beg, idwbuff(1)%end
                                jac(j,n+k,l) = jac(j,k-1,l)
                            end do
                        end do
                    end do
                else if (bcye == -2) then
                    !$acc parallel loop collapse(3) gang vector default(present)
                    do l = 0, p
                        do k = 1, buff_size
                            do j = idwbuff(1)%beg, idwbuff(1)%end
                                jac(j,n+k,l) = jac(j,n - (k-1),l)
                            end do
                        end do
                    end do
                else
                    !$acc parallel loop collapse(3) gang vector default(present)
                    do l = 0, p
                        do k = 1, buff_size
                            do j = idwbuff(1)%beg, idwbuff(1)%end
                                jac(j,n+k,l) = jac(j,n,l)
                            end do
                        end do
                    end do

                end if
            end if

            if(p > 0) then
                if(bczb >= -12) then
                    if(bczb >= 0) then
                        call s_mpi_sendrecv_F_igr(jac, 3, -1)
                    else if (bczb == -1) then
                        !$acc parallel loop collapse(3) gang vector default(present)
                        do l = 1, buff_size
                            do k = idwbuff(2)%beg, idwbuff(2)%end
                                do j = idwbuff(1)%beg, idwbuff(1)%end
                                    jac(j,k,-l) = jac(j,k,p-l+1)
                                end do
                            end do
                        end do
                    else if (bczb == -2) then
                        !$acc parallel loop collapse(3) gang vector default(present)
                        do l = 1, buff_size
                            do k = idwbuff(2)%beg, idwbuff(2)%end
                                do j = idwbuff(1)%beg, idwbuff(1)%end
                                    jac(j,k,-l) = jac(j,k,l-1)
                                end do
                            end do
                        end do
                    else
                        !$acc parallel loop collapse(3) gang vector default(present)
                        do l = 1, buff_size
                            do k = idwbuff(2)%beg, idwbuff(2)%end
                                do j = idwbuff(1)%beg, idwbuff(1)%end
                                    jac(j,k,-l) = jac(j,k,0)
                                end do
                            end do
                        end do

                    end if
                end if

                if(bcze >= -12) then
                    if(bcze >= 0) then
                        call s_mpi_sendrecv_F_igr(jac, 3, 1)
                    else if (bcze == -1) then
                        !$acc parallel loop gang vector collapse(3) default(present)
                        do l = 1, buff_size
                            do k = idwbuff(2)%beg, idwbuff(2)%end
                                do j = idwbuff(1)%beg, idwbuff(1)%end
                                    jac(j,k,p+l) = jac(j,k,l-1)
                                end do
                            end do
                        end do
                    else if (bcze == -2) then
                        !$acc parallel loop gang vector collapse(3) default(present)
                        do l = 1, buff_size
                            do k = idwbuff(2)%beg, idwbuff(2)%end
                                do j = idwbuff(1)%beg, idwbuff(1)%end
                                    jac(j,k,p+l) = jac(j,k,p - (l-1))
                                end do
                            end do
                        end do
                    else
                        !$acc parallel loop gang vector collapse(3) default(present)
                        do l = 1, buff_size
                            do k = idwbuff(2)%beg, idwbuff(2)%end
                                do j = idwbuff(1)%beg, idwbuff(1)%end
                                    jac(j,k,p+l) = jac(j,k,p)
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

    subroutine s_compute_mixture(rho_L, rho_R, gamma_L, gamma_R, pi_inf_L, pi_inf_R, &
                                    Re_L, Re_R, cfl, j, k, l, idir, alpha_L, alpha_rho_L, alpha_R, alpha_rho_R, vel_L, vel_R, pres_L, pres_R)
#ifdef _CRAYFTN
        !DIR$ INLINEALWAYS s_compute_mixture
#else
        !$acc routine seq
#endif
        real(wp) :: rho_L, gamma_L, pi_inf_L, Re_L
        real(wp) :: rho_R, gamma_R, pi_inf_R, Re_R
        real(wp), dimension(num_fluids) :: alpha_L, alpha_rho_L, alpha_R, alpha_rho_R
        real(wp), dimension(num_dims) :: vel_L, vel_R
        real(wp) :: pres_L, pres_R
        real(wp) :: a_L, a_R, cfl
        real(wp) :: u_L, u_R, v_L, v_R, w_L, w_R
        integer :: i, j, k, l, q, idir

        rho_L = 0._wp; gamma_L = 0._wp; pi_inf_L = 0._wp

        do i = 1, num_fluids
            rho_L = rho_L + alpha_rho_L(i)
            gamma_L = gamma_L + alpha_L(i)*gammas(i)
            pi_inf_L = pi_inf_L + alpha_L(i)*pi_infs(i)
        end do

        a_L = sqrt((pres_L*(1._wp/gamma_L + 1._wp) + pi_inf_L / gamma_L) / rho_L)
        u_L = vel_L(1)
        v_L = vel_L(2)
        if (p > 0) w_L = vel_L(3)

        if (viscous) then
            Re_L = dflt_real
            if (Re_size(1) > 0) Re_L = 0._wp
            !$acc loop seq
            do q = 1, Re_size(1)
                Re_L =  alpha_L(Re_idx(1, q)) / Res(1, q) &
                          + Re_L
            end do
            Re_L = 1._wp/max(Re_L, sgm_eps)
        end if

        rho_R = 0._wp; gamma_R = 0._wp; pi_inf_R = 0._wp

        do i = 1, num_fluids
            rho_R = rho_R + alpha_rho_R(i)
            gamma_R = gamma_R + alpha_R(i)*gammas(i)
            pi_inf_R = pi_inf_R + alpha_R(i)*pi_infs(i)
        end do

        a_R = sqrt((pres_R*(1._wp/gamma_R + 1._wp) + pi_inf_R / gamma_R) / rho_R)
        u_R = vel_R(1)
        v_R = vel_R(2)
        if (p > 0) w_R = vel_R(3)

        if (viscous) then
            Re_R = dflt_real
            if (Re_size(1) > 0) Re_R = 0._wp
            !$acc loop seq
            do q = 1, Re_size(1)
                Re_R = alpha_R(Re_idx(1, q))/Res(1, q) &
                          + Re_R
            end do
            Re_R = 1._wp/max(Re_R, sgm_eps)
        end if

        if (p == 0) then
            cfl = max(sqrt(u_L**2._wp + v_L**2._wp), sqrt(u_R**2._wp + v_R**2._wp)) + max(a_L, a_R)
        else
            cfl = max(sqrt(u_L**2._wp + v_L**2._wp + w_L**2._wp), sqrt(u_R**2._wp + v_R**2._wp + w_R**2._wp)) + max(a_L, a_R)
        end if

    end subroutine s_compute_mixture

    subroutine s_igr_riemann_solver(flux_vf, idir)

        type(scalar_field), &
            dimension(sys_size), &
            intent(inout) :: flux_vf
        integer, intent(in) :: idir

        real(wp) :: rho_L, gamma_L, pi_inf_L, Re_L, mu_L
        real(wp) :: rho_R, gamma_R, pi_inf_R, Re_R, mu_R
        real(wp) :: cfl
        real(wp), dimension(num_fluids) :: alpha_L, alpha_rho_L, alpha_R, alpha_rho_R
        real(wp), dimension(num_dims) :: vel_L, vel_R
        real(wp) :: pres_L, pres_R

        if (idir == 1) then

            if(p == 0) then
                !$acc parallel loop collapse(3) gang vector default(present) private(alpha_L, alpha_rho_L, alpha_R, alpha_rho_R, vel_L, vel_R, pres_L, pres_R)
                do l = 0, p
                    do k = -1, n+1
                        do j = -1, m+1

                            !$acc loop seq
                            do i = 1, num_fluids
                                alpha_rho_L(i) = qL_rs_vf(j+1, k, l, i)
                                alpha_L(i) = qL_rs_vf(j+1, k, l, E_idx + i)

                                alpha_rho_R(i) = qR_rs_vf(j, k, l, i)
                                alpha_R(i) = qR_rs_vf(j, k, l, E_idx + i)
                            end do

                            !$acc loop seq
                            do i = 1, num_dims
                                vel_L(i) = qL_rs_vf(j+1, k, l, contxe + i)
                                vel_R(i) = qR_rs_vf(j, k, l, contxe + i)
                            end do

                            pres_L = qL_rs_vf(j+1, k, l, E_idx)
                            pres_R = qR_rs_vf(j, k, l, E_idx)

                            call s_compute_mixture(rho_L, rho_R, gamma_L, gamma_R, &
                                                     pi_inf_L, pi_inf_R, Re_L, Re_r, cfl, j, k, l, idir, alpha_L, alpha_rho_L, alpha_R, alpha_rho_R, vel_L, vel_R, pres_L, pres_R)

                            do i = 1, num_fluids
                                flux_vf(advxb+i-1)%sf(j,k,l) = &
                                    0.5_wp * (qL_rs_vf(j+1,k,l,advxb+i-1) * &
                                    qL_rs_vf(j+1,k,l,momxb)) + &
                                    0.5_wp * (qR_rs_vf(j,k,l,advxb+i-1) * &
                                    qR_rs_vf(j,k,l,momxb)) + &
                                    05_wp*cfl * (qR_rs_vf(j, k, l, advxb+i-1) - qL_rs_vf(j+1, k, l, advxb+i-1))
                                flux_vf(i)%sf(j,k,l) = &
                                    0.5_wp * (qL_rs_vf(j+1,k,l,i) * &
                                    qL_rs_vf(j+1,k,l, momxb)) + &
                                    0.5_wp * (qR_rs_vf(j,k,l,i) * &
                                    qR_rs_vf(j,k,l, momxb)) + &
                                    0.5_wp*cfl * (qR_rs_vf(j, k, l, i) - qL_rs_vf(j+1, k, l, i))
                            end do

                            ! Momentum -> rho*u^2 + p + [[F_igr]]
                            flux_vf(momxb)%sf(j,k,l) = &
                                 0.5_wp* (rho_L * (qL_rs_vf(j+1,k,l,momxb)**2.0) + &
                                 qL_rs_vf(j+1,k,l,E_idx) + FL(j+1, k, l) ) + &
                                 0.5_wp* (rho_R * (qR_rs_vf(j,k,l,momxb)**2.0) + &
                                 qR_rs_vf(j,k,l,E_idx) + FR(j, k, l) ) + &
                                 0.5_wp*cfl * (qR_rs_vf(j, k, l, momxb) - qL_rs_vf(j+1, k, l, momxb))

                            flux_vf(momxb+1)%sf(j, k, l) =  &
                                0.5_wp * rho_L * qL_rs_vf(j+1,k,l,momxb)*qL_rs_vf(j+1,k,l,momxb+1) + &
                                0.5_wp * rho_R * qR_rs_vf(j,k,l,momxb)*qR_rs_vf(j,k,l,momxb+1) + &
                                0.5_wp*cfl * (qR_rs_vf(j, k, l, momxb+1) - qL_rs_vf(j+1, k, l, momxb+1))

                             flux_vf(E_idx)%sf(j, k, l) = &
                                0.5_wp * (qL_rs_vf(j+1,k,l,momxb) * (qL_rs_vf(j+1,k,l,E_idx)*gamma_L + pi_inf_L + &
                                0.5_wp * rho_L * (qL_rs_vf(j+1, k, l,momxb)**2._wp + qL_rs_vf(j+1, k, l,momxb+1)**2._wp ) + &
                                qL_rs_vf(j+1,k,l,E_idx) + FL(j+1, k, l)) ) + &
                                0.5_wp * (qR_rs_vf(j,k,l,momxb) * (qR_rs_vf(j,k,l,E_idx)*gamma_R + pi_inf_R + &
                                0.5_wp * rho_R * (qR_rs_vf(j, k, l,momxb)**2._wp + qR_rs_vf(j, k, l,momxb+1)**2._wp ) + &
                                qR_rs_vf(j,k,l,E_idx) + FR(j, k, l)) ) + &
                                0.5_wp*cfl * (qR_rs_vf(j, k, l, E_idx) - qL_rs_vf(j+1, k, l, E_idx))

                            if(viscous) then
                                mu_L = 1/Re_L; mu_R = 1/Re_R

                                flux_vf(momxb)%sf(j, k, l) = flux_vf(momxb)%sf(j, k, l) - &
                                            0.5_wp * mu_L * ((4._wp/3._wp)*duLx(j+1, k, l) - (2._wp/3._wp)*dvLy(j+1, k, l)) - &
                                            0.5_wp * mu_R * ((4._wp/3._wp)*duRx(j, k, l) - (2._wp/3._wp)*dvRy(j, k, l))

                                flux_vf(momxb+1)%sf(j, k, l) = flux_vf(momxb+1)%sf(j, k, l) - &
                                           0.5_wp * mu_L * (duLy(j+1, k, l) + dvLx(j+1, k, l))  -  &
                                           0.5_wp * mu_R * (duRy(j, k, l) + dvRx(j, k, l))

                                flux_vf(E_idx)%sf(j, k, l) = flux_vf(E_idx)%sf(j, k, l) - &
                                    0.5_wp*mu_L*qL_rs_vf(j+1, k, l, momxb)*((4._wp/3._wp)*duLx(j+1, k, l) - (2._wp/3._wp)*dvLy(j+1, k, l)) - &
                                    0.5_wp*mu_L*qL_rs_vf(j+1, k, l, momxb)*(duLy(j+1, k, l) + dvLx(j+1, k, l))   - &
                                    0.5_wp*mu_R*qR_rs_vf(j, k, l, momxb)*((4._wp/3._wp)*duRx(j, k, l) - (2._wp/3._wp)*dvRy(j, k, l)) - &
                                    0.5_wp*mu_R*qR_rs_vf(j, k, l, momxb)*(duRy(j, k, l) + dvRx(j, k, l))
                            end if

                        end do
                    end do
                end do
            else
                !$acc parallel loop collapse(3) gang vector default(present) private(alpha_L, alpha_rho_L, alpha_R, alpha_rho_R, vel_L, vel_R, pres_L, pres_R)
                do l = -1, p+1
                    do k = -1, n+1
                        do j = -1, m+1

                            !$acc loop seq
                            do i = 1, num_fluids
                                alpha_rho_L(i) = qL_rs_vf(j+1, k, l, i)
                                alpha_L(i) = qL_rs_vf(j+1, k, l, E_idx + i)

                                alpha_rho_R(i) = qR_rs_vf(j, k, l, i)
                                alpha_R(i) = qR_rs_vf(j, k, l, E_idx + i)
                            end do

                            !$acc loop seq
                            do i = 1, num_dims
                                vel_L(i) = qL_rs_vf(j+1, k, l, contxe + i)
                                vel_R(i) = qR_rs_vf(j, k, l, contxe + i)
                            end do

                            pres_L = qL_rs_vf(j+1, k, l, E_idx)
                            pres_R = qR_rs_vf(j, k, l, E_idx)

                            call s_compute_mixture(rho_L, rho_R, gamma_L, gamma_R, &
                                                     pi_inf_L, pi_inf_R, Re_L, Re_r, cfl, j, k, l, idir, alpha_L, alpha_rho_L, alpha_R, alpha_rho_R, vel_L, vel_R, pres_L, pres_R)

                            !$acc loop seq
                            do i = 1, num_fluids
                                flux_vf(advxb+i-1)%sf(j,k,l) = &
                                    0.5_wp * (qL_rs_vf(j+1,k,l,advxb+i-1) * &
                                    qL_rs_vf(j+1,k,l,momxb)) + &
                                    0.5_wp * (qR_rs_vf(j,k,l,advxb+i-1) * &
                                    qR_rs_vf(j,k,l,momxb)) + &
                                    0.5_wp*cfl * (qR_rs_vf(j, k, l, advxb+i-1) - qL_rs_vf(j+1, k, l, advxb+i-1))
                                flux_vf(i)%sf(j,k,l) = &
                                    0.5_wp * (qL_rs_vf(j+1,k,l,i) * &
                                    qL_rs_vf(j+1,k,l, momxb)) + &
                                    0.5_wp * (qR_rs_vf(j,k,l,i) * &
                                    qR_rs_vf(j,k,l, momxb)) + &
                                    0.5_wp*cfl * (qR_rs_vf(j, k, l, i) - qL_rs_vf(j+1, k, l, i))
                            end do

                            ! Momentum -> rho*u^2 + p + [[F_igr]]
                            flux_vf(momxb)%sf(j,k,l) = &
                                0.5_wp* (rho_L * (qL_rs_vf(j+1,k,l,momxb)**2.0) + &
                                qL_rs_vf(j+1,k,l,E_idx) + FL(j+1, k, l) ) + &
                                0.5_wp* (rho_R * (qR_rs_vf(j,k,l,momxb)**2.0) + &
                                qR_rs_vf(j,k,l,E_idx) + FR(j, k, l) ) + &
                                0.5_wp*cfl * (qR_rs_vf(j, k, l, momxb) - qL_rs_vf(j+1, k, l, momxb))

                            flux_vf(momxb+1)%sf(j, k, l) = &
                                0.5_wp * rho_L * (qL_rs_vf(j+1,k,l,momxb)*qL_rs_vf(j+1,k,l,momxb+1)) + &
                                0.5_wp * rho_R * (qR_rs_vf(j,k,l,momxb)*qR_rs_vf(j,k,l,momxb+1)) + &
                                0.5_wp*cfl * (qR_rs_vf(j, k, l, momxb+1) - qL_rs_vf(j+1, k, l, momxb+1))

                            flux_vf(momxe)%sf(j, k, l) = &
                                0.5_wp * rho_L * (qL_rs_vf(j+1,k,l,momxb)*qL_rs_vf(j+1,k,l,momxe)) + &
                                0.5_wp*  rho_R * (qR_rs_vf(j,k,l,momxb)*qR_rs_vf(j,k,l,momxe)) + &
                                0.5_wp*cfl * (qR_rs_vf(j, k, l, momxe) - qL_rs_vf(j+1, k, l, momxe))

                            flux_vf(E_idx)%sf(j, k, l) = &
                                0.5_wp * ( qL_rs_vf(j+1,k,l,momxb) * (qL_rs_vf(j+1,k,l,E_idx)*gamma_L + pi_inf_L + &
                                0.5_wp * rho_L * (qL_rs_vf(j+1, k, l,momxb)**2._wp + qL_rs_vf(j+1, k, l,momxb+1)**2._wp + qL_rs_vf(j+1, k, l,momxe)**2._wp) + &
                                qL_rs_vf(j+1,k,l,E_idx) + FL(j+1, k, l)) ) + &
                                0.5_wp * ( qR_rs_vf(j,k,l,momxb) * (qR_rs_vf(j,k,l,E_idx)*gamma_R + pi_inf_R + &
                                0.5_wp * rho_R * (qR_rs_vf(j, k, l,momxb)**2._wp + qR_rs_vf(j, k, l,momxb+1)**2._wp + qR_rs_vf(j, k, l,momxe)**2._wp) + &
                                qR_rs_vf(j,k,l,E_idx) + FR(j, k, l)) ) + &
                                0.5_wp*cfl * (qR_rs_vf(j, k, l, E_idx) - qL_rs_vf(j+1, k, l, E_idx))

                            if(viscous) then
                                mu_L = 1/Re_L; mu_R = 1/Re_R

                                flux_vf(momxe)%sf(j, k, l) = flux_vf(momxe)%sf(j, k, l) - &
                                           0.5_wp*mu_L*(duLz(j+1, k, l) + dwLx(j+1, k, l))  -  &
                                           0.5_wp*mu_R*(duRz(j, k, l) + dwRx(j, k, l))

                                flux_vf(momxb+1)%sf(j, k, l) = flux_vf(momxb+1)%sf(j, k, l) - &
                                           0.5_wp*mu_L*(duLy(j+1, k, l) + dvLx(j+1, k, l))  -  &
                                           0.5_wp*mu_R*(duRy(j, k, l) + dvRx(j, k, l))

                                flux_vf(momxb)%sf(j, k, l) = flux_vf(momxb)%sf(j, k, l) - &
                                            0.5_wp * mu_L * ((4._wp/3._wp)*duLx(j+1, k, l) - (2._wp/3._wp)*dvLy(j+1, k, l) - (2._wp/3._wp) * dwLz(j+1, k, l)) - &
                                            0.5_wp * mu_R * ((4._wp/3._wp)*duRx(j, k, l) - (2._wp/3._wp)*dvRy(j, k, l) - (2._wp/3._wp) * dwRz(j, k, l))

                                flux_vf(E_idx)%sf(j, k, l) = flux_vf(E_idx)%sf(j, k, l) - &
                                    0.5_wp*mu_L*qL_rs_vf(j+1, k, l, momxb)*((4._wp/3._wp)*duLx(j+1, k, l) - (2._wp/3._wp)*dvLy(j+1, k, l) - (2._wp/3._wp)*dwLz(j+1, k, l)) - &
                                    0.5_wp*mu_L*qL_rs_vf(j+1, k, l, momxb)*(duLy(j+1, k, l) + dvLx(j+1, k, l))   - &
                                    0.5_wp*mu_L*qL_rs_vf(j+1, k, l, momxb)*(duLz(j+1, k, l) + dwLx(j+1, k, l))   - &
                                    0.5_wp*mu_R*qR_rs_vf(j, k, l, momxb)*((4._wp/3._wp)*duRx(j, k, l) - (2._wp/3._wp)*dvRy(j, k, l) - (2._wp/3._wp)*dwRz(j, k, l)) - &
                                    0.5_wp*mu_R*qR_rs_vf(j, k, l, momxb)*(duRy(j, k, l) + dvRx(j, k, l)) - &
                                    0.5_wp*mu_R*qR_rs_vf(j, k, l, momxb)*(duRz(j, k, l) + dwRx(j, k, l))
                            end if
                        end do
                    end do
                end do
            end if
        else if (idir == 2) then
            if(p == 0) then
                !$acc parallel loop collapse(3) gang vector default(present) private(alpha_L, alpha_rho_L, alpha_R, alpha_rho_R, vel_L, vel_R, pres_L, pres_R)
                do l = 0, p
                    do k = -1, m+1
                        do j = -1, n+1

                            !$acc loop seq
                            do i = 1, num_fluids
                                alpha_rho_L(i) = qL_rs_vf(j, k+1, l, i)
                                alpha_L(i) = qL_rs_vf(j, k+1, l, E_idx + i)

                                alpha_rho_R(i) = qR_rs_vf(j, k, l, i)
                                alpha_R(i) = qR_rs_vf(j, k, l, E_idx + i)
                            end do


                            !$acc loop seq
                            do i = 1, num_dims
                                vel_L(i) = qL_rs_vf(j, k+1, l, contxe + i)
                                vel_R(i) = qR_rs_vf(j, k, l, contxe + i)
                            end do

                            pres_L = qL_rs_vf(j, k+1, l, E_idx)
                            pres_R = qR_rs_vf(j, k, l, E_idx)

                            call s_compute_mixture(rho_L, rho_R, gamma_L, gamma_R, &
                                                     pi_inf_L, pi_inf_R, Re_L, Re_r, cfl, j, k, l, idir, alpha_L, alpha_rho_L, alpha_R, alpha_rho_R, vel_L, vel_R, pres_L, pres_R)

                            do i = 1, num_fluids
                                flux_vf(advxb+i-1)%sf(j,k,l) = &
                                    0.5_wp * (qL_rs_vf(j,k+1,l,advxb+i-1) * &
                                    qL_rs_vf(j,k+1,l,momxb+1)) + &
                                    0.5_wp * (qR_rs_vf(j,k,l,advxb+i-1) * &
                                    qR_rs_vf(j,k,l,momxb+1)) + &
                                    0.5_wp*cfl * (qR_rs_vf(j, k, l, advxb+i-1) - qL_rs_vf(j, k+1, l, advxb+i-1))
                                flux_vf(i)%sf(j,k,l) = &
                                    0.5_wp * (qL_rs_vf(j,k+1,l,i) * &
                                    qL_rs_vf(j,k+1,l, momxb+1)) + &
                                    0.5_wp * (qR_rs_vf(j,k,l,i) * &
                                    qR_rs_vf(j,k,l, momxb+1)) + &
                                    0.5_wp*cfl * (qR_rs_vf(j, k, l, i) - qL_rs_vf(j, k+1, l, i))
                            end do

                            flux_vf(momxb+1)%sf(j, k, l) = &
                                 0.5_wp * (rho_L * (qL_rs_vf(j,k+1,l,momxb+1)**2.0) + &
                                 qL_rs_vf(j,k+1,l,E_idx) + FL(j, k+1, l) ) + &
                                 0.5_wp * (rho_R * (qR_rs_vf(j,k,l,momxb+1)**2.0) + &
                                 qR_rs_vf(j,k,l,E_idx) + FR(j, k, l) ) + &
                                 0.5_wp*cfl * (qR_rs_vf(j, k, l, momxb+1) - qL_rs_vf(j, k+1, l, momxb+1))

                            flux_vf(momxb)%sf(j, k, l) = &
                                0.5_wp * rho_L * (qL_rs_vf(j,k+1,l,momxb)*qL_rs_vf(j,k+1,l,momxb+1)) + &
                                0.5_wp * rho_R * (qR_rs_vf(j,k,l,momxb)*qR_rs_vf(j,k,l,momxb+1)) + &
                                0.5_wp*cfl * (qR_rs_vf(j, k, l, momxb) - qL_rs_vf(j, k+1, l, momxb))

                            flux_vf(E_idx)%sf(j, k, l) = &
                                0.5_wp * ( qL_rs_vf(j,k+1,l,momxb+1) * (qL_rs_vf(j,k+1,l,E_idx)*gamma_L + pi_inf_L + &
                                0.5_wp * rho_L * (qL_rs_vf(j, k+1, l,momxb)**2._wp + qL_rs_vf(j, k+1, l,momxb+1)**2._wp ) + &
                                qL_rs_vf(j,k+1,l,E_idx) + FL(j, k+1, l)) ) + &
                                0.5_wp * ( qR_rs_vf(j,k,l,momxb+1) * (qR_rs_vf(j,k,l,E_idx)*gamma_R + pi_inf_R + &
                                0.5_wp * rho_R * (qR_rs_vf(j, k, l,momxb)**2._wp + qR_rs_vf(j, k, l,momxb+1)**2._wp ) + &
                                qR_rs_vf(j,k,l,E_idx) + FR(j, k, l)) ) + &
                                0.5_wp*cfl * (qR_rs_vf(j, k, l, E_idx) - qL_rs_vf(j, k+1, l, E_idx))

                            if(viscous) then
                                mu_L = 1/Re_L; mu_R = 1/Re_R

                                flux_vf(momxb+1)%sf(j, k, l) = flux_vf(momxb+1)%sf(j, k, l) - &
                                            0.5_wp * mu_L * ((4._wp/3._wp)*dvLy(j, k+1, l) - (2._wp/3._wp)*duLx(j, k+1, l)) - &
                                            0.5_wp * mu_R * ((4._wp/3._wp)*dvRy(j, k, l) - (2._wp/3._wp)*duRx(j, k, l) )

                                flux_vf(momxb)%sf(j, k, l) = flux_vf(momxb)%sf(j, k, l) - &
                                           0.5_wp * mu_L * (duLy(j, k+1, l) + dvLx(j, k+1, l))  -  &
                                           0.5_wp * mu_R * (duRy(j, k, l) + dvRx(j, k, l))

                                flux_vf(E_idx)%sf(j, k, l) = flux_vf(E_idx)%sf(j, k, l) - &
                                    0.5_wp*mu_L*qL_rs_vf(j, k+1, l, momxb+1)*((4._wp/3._wp)*dvLy(j, k+1, l) - (2._wp/3._wp)*duLx(j,k+1, l)) - &
                                    0.5_wp*mu_L*qL_rs_vf(j, k+1, l, momxb+1)*(duLy(j, k+1, l) + dvLx(j, k+1, l))   - &
                                    0.5_wp*mu_R*qR_rs_vf(j, k, l, momxb+1)*((4._wp/3._wp)*dvRy(j, k, l) - (2._wp/3._wp)*duRx(j, k, l)) - &
                                    0.5_wp*mu_R*qR_rs_vf(j, k, l, momxb+1)*(duRy(j, k, l) + dvRx(j, k, l))
                            end if
                        end do
                    end do
                end do
            else
                !$acc parallel loop collapse(3) gang vector default(present) private(alpha_L, alpha_rho_L, alpha_R, alpha_rho_R, vel_L, vel_R, pres_L, pres_R)
                do l = -1, p+1
                    do k = -1, n+1
                         do j = -1, m+1

                            !$acc loop seq
                            do i = 1, num_fluids
                                alpha_rho_L(i) = qL_rs_vf(j, k+1, l, i)
                                alpha_L(i) = qL_rs_vf(j, k+1, l, E_idx + i)

                                alpha_rho_R(i) = qR_rs_vf(j, k, l, i)
                                alpha_R(i) = qR_rs_vf(j, k, l, E_idx + i)
                            end do

                            !$acc loop seq
                            do i = 1, num_dims
                                vel_L(i) = qL_rs_vf(j, k+1, l, contxe + i)
                                vel_R(i) = qR_rs_vf(j, k, l, contxe + i)
                            end do

                            pres_L = qL_rs_vf(j, k+1, l, E_idx)
                            pres_R = qR_rs_vf(j, k, l, E_idx)

                            call s_compute_mixture(rho_L, rho_R, gamma_L, gamma_R, &
                                                     pi_inf_L, pi_inf_R, Re_L, Re_r, cfl, j, k, l, idir, alpha_L, alpha_rho_L, alpha_R, alpha_rho_R, vel_L, vel_R, pres_L, pres_R)

                            do i = 1, num_fluids
                                flux_vf(advxb+i-1)%sf(j,k,l) = &
                                    0.5_wp * (qL_rs_vf(j,k+1,l,advxb+i-1) * &
                                    qL_rs_vf(j,k+1,l,momxb+1)) + &
                                    0.5_wp * (qR_rs_vf(j,k,l,advxb+i-1) * &
                                    qR_rs_vf(j,k,l,momxb+1)) + &
                                    0.5_wp*cfl * (qR_rs_vf(j, k, l, advxb+i-1) - qL_rs_vf(j, k+1, l, advxb+i-1))
                                flux_vf(i)%sf(j,k,l) = &
                                    0.5_wp * (qL_rs_vf(j,k+1,l,i) * &
                                    qL_rs_vf(j,k+1,l, momxb+1)) + &
                                    0.5_wp * (qR_rs_vf(j,k,l,i) * &
                                    qR_rs_vf(j,k,l, momxb+1)) + &
                                    0.5_wp*cfl* (qR_rs_vf(j, k, l, i) - qL_rs_vf(j, k+1, l, i))
                            end do

                            flux_vf(momxb)%sf(j, k, l) = &
                                0.5_wp * rho_L * (qL_rs_vf(j,k+1,l,momxb)*qL_rs_vf(j,k+1,l,momxb+1)) + &
                                0.5_wp * rho_R * (qR_rs_vf(j,k,l,momxb)*qR_rs_vf(j,k,l,momxb+1)) + &
                                0.5_wp*cfl * (qR_rs_vf(j, k, l, momxb) - qL_rs_vf(j, k+1, l, momxb))

                            flux_vf(momxb+1)%sf(j, k, l) = &
                                0.5_wp* (rho_L * (qL_rs_vf(j,k+1,l,momxb+1)**2.0) + &
                                qL_rs_vf(j,k+1,l,E_idx) + FL(j, k+1, l) ) + &
                                0.5_wp* (rho_R * (qR_rs_vf(j,k,l,momxb+1)**2.0) + &
                                qR_rs_vf(j,k,l,E_idx) + FR(j, k, l) ) + &
                                0.5_wp*cfl * (qR_rs_vf(j, k, l, momxb+1) - qL_rs_vf(j, k+1, l, momxb+1))

                            flux_vf(momxe)%sf(j, k, l) = &
                                0.5_wp* rho_L * (qL_rs_vf(j,k+1,l,momxb+1)*qL_rs_vf(j,k+1,l,momxe)) + &
                                0.5_wp* rho_R * (qR_rs_vf(j,k,l,momxb+1)*qR_rs_vf(j,k,l,momxe)) + &
                                0.5_wp*cfl * (qR_rs_vf(j, k, l, momxe) - qL_rs_vf(j, k+1, l, momxe))

                            flux_vf(E_idx)%sf(j, k, l) = &
                                0.5_wp * ( qL_rs_vf(j,k+1,l,momxb+1) * (qL_rs_vf(j,k+1,l,E_idx)*gamma_L + pi_inf_L + &
                                0.5_wp * rho_L * (qL_rs_vf(j, k+1, l,momxb)**2._wp + qL_rs_vf(j, k+1, l,momxb+1)**2._wp + qL_rs_vf(j,k+1,l,momxe)**2._wp) + &
                                qL_rs_vf(j,k+1,l,E_idx) + FL(j, k+1, l)) ) + &
                                0.5_wp * ( qR_rs_vf(j,k,l,momxb+1) * (qR_rs_vf(j,k,l,E_idx)*gamma_R + pi_inf_R + &
                                0.5_wp * rho_R * (qR_rs_vf(j, k, l,momxb)**2._wp + qR_rs_vf(j, k, l,momxb+1)**2._wp + qR_rs_vf(j,k,l,momxe)**2._wp ) + &
                                qR_rs_vf(j,k,l,E_idx) + FR(j, k, l)) ) + &
                                0.5_wp*cfl * (qR_rs_vf(j, k, l, E_idx) - qL_rs_vf(j, k+1, l, E_idx))

                            if(viscous) then
                                mu_L = 1/Re_L; mu_R = 1/Re_R

                                flux_vf(momxb+1)%sf(j, k, l) = flux_vf(momxb+1)%sf(j, k, l) - &
                                            0.5_wp * mu_L*((4._wp/3._wp)*dvLy(j, k+1, l) - (2._wp/3._wp)*duLx(j, k+1, l) - (2._wp/3._wp)*dwLz(j,k+1,l)) - &
                                            0.5_wp * mu_R*((4._wp/3._wp)*dvRy(j, k, l) - (2._wp/3._wp)*duRx(j, k, l) - (2._wp/3._wp)*dwRz(j,k,l) )

                                flux_vf(momxb)%sf(j, k, l) = flux_vf(momxb)%sf(j, k, l) - &
                                           0.5_wp*mu_L*(duLy(j, k+1, l) + dvLx(j, k+1, l))  -  &
                                           0.5_wp*mu_R*(duRy(j, k, l) + dvRx(j, k, l))

                                flux_vf(momxe)%sf(j, k, l) = flux_vf(momxe)%sf(j, k, l) - &
                                           0.5_wp*mu_L*(dvLz(j, k+1, l) + dwLy(j, k+1, l))  -  &
                                           0.5_wp*mu_R*(dvRz(j, k, l) + dwRy(j, k, l))

                                flux_vf(E_idx)%sf(j, k, l) = flux_vf(E_idx)%sf(j, k, l) - &
                                    0.5_wp*mu_L*qL_rs_vf(j, k+1, l, momxb+1)*((4._wp/3._wp)*dvLy(j, k+1, l) - (2._wp/3._wp)*duLx(j, k+1, l) - (2._wp/3._wp)*dwLz(j,k+1 ,l)) - &
                                    0.5_wp*mu_L*qL_rs_vf(j, k+1, l, momxb+1)*(duLy(j, k+1, l) + dvLx(j, k+1, l))   - &
                                    0.5_wp*mu_L*qL_rs_vf(j, k+1, l, momxb+1)*(dwLy(j, k+1, l) + dvLz(j, k+1, l))   - &
                                    0.5_wp*mu_R*qR_rs_vf(j, k, l, momxb+1)*((4._wp/3._wp)*dvRy(j, k, l) - (2._wp/3._wp)*duRx(j, k, l) - (2._wp/3._wp)*dwRz(j ,k ,l)) - &
                                    0.5_wp*mu_R*qR_rs_vf(j, k, l, momxb+1)*(duRy(j, k, l) + dvRx(j, k, l))   - &
                                    0.5_wp*mu_R*qR_rs_vf(j, k, l, momxb+1)*(dwRy(j, k, l) + dvRz(j, k, l))
                            end if
                        end do
                    end do
                end do
            end if
        elseif (idir == 3) then
            !$acc parallel loop collapse(3) gang vector default(present) private(alpha_L, alpha_rho_L, alpha_R, alpha_rho_R, vel_L, vel_R, pres_L, pres_R)
            do l = -1, p+1
                do k = -1, n+1
                    do j = -1, m+1

                        !$acc loop seq
                        do i = 1, num_fluids
                            alpha_rho_L(i) = qL_rs_vf(j, k, l+1, i)
                            alpha_L(i) = qL_rs_vf(j, k, l+1, E_idx + i)

                            alpha_rho_R(i) = qR_rs_vf(j, k, l, i)
                            alpha_R(i) = qR_rs_vf(j, k, l, E_idx + i)
                        end do

                        !$acc loop seq
                        do i = 1, num_dims
                            vel_L(i) = qL_rs_vf(j, k, l+1, contxe + i)
                            vel_R(i) = qR_rs_vf(j, k, l, contxe + i)
                        end do

                        pres_L = qL_rs_vf(j, k, l+1, E_idx)
                        pres_R = qR_rs_vf(j, k, l, E_idx)

                        call s_compute_mixture(rho_L, rho_R, gamma_L, gamma_R, &
                                             pi_inf_L, pi_inf_R, Re_L, Re_r, cfl, j, k, l, idir, alpha_L, alpha_rho_L, alpha_R, alpha_rho_R, vel_L, vel_R, pres_L, pres_R)


                        do i = 1, num_fluids
                            flux_vf(advxb+i-1)%sf(j,k,l) = &
                                0.5_wp * (qL_rs_vf(j,k,l+1,advxb+i-1) * &
                                qL_rs_vf(j,k,l+1,momxe)) + &
                                0.5_wp * (qR_rs_vf(j,k,l,advxb+i-1) * &
                                qR_rs_vf(j,k,l,momxe)) + &
                                0.5_wp*cfl * (qR_rs_vf(j, k, l, advxb+i-1) - qL_rs_vf(j, k, l+1, advxb+i-1))
                            flux_vf(i)%sf(j,k,l) = &
                                0.5_wp * (qL_rs_vf(j,k,l+1,i) * &
                                qL_rs_vf(j,k,l+1, momxe)) + &
                                0.5_wp * (qR_rs_vf(j,k,l,i) * &
                                qR_rs_vf(j,k,l, momxe)) + &
                                0.5_wp * cfl * (qR_rs_vf(j, k, l, i) - qL_rs_vf(j, k, l+1, i))
                        end do

                        flux_vf(momxb)%sf(j, k, l) = &
                            0.5_wp * rho_L * (qL_rs_vf(j,k,l+1,momxb)*qL_rs_vf(j,k,l+1,momxe)) + &
                            0.5_wp * rho_R * (qR_rs_vf(j,k,l,momxb)*qR_rs_vf(j,k,l,momxe)) + &
                            0.5_wp * cfl * (qR_rs_vf(j, k, l, momxb) - qL_rs_vf(j, k, l+1, momxb))

                        flux_vf(momxb+1)%sf(j, k, l) = &
                            0.5_wp * rho_L * (qL_rs_vf(j,k,l+1,momxb+1)*qL_rs_vf(j,k,l+1,momxe)) + &
                            0.5_wp * rho_R * (qR_rs_vf(j,k,l,momxb+1)*qR_rs_vf(j,k,l,momxe)) + &
                            0.5_wp * cfl * (qR_rs_vf(j, k, l, momxb+1) - qL_rs_vf(j, k, l+1, momxb+1))

                        flux_vf(momxe)%sf(j, k, l) = &
                             0.5_wp* (rho_L * (qL_rs_vf(j,k,l+1,momxe)**2.0) + &
                             qL_rs_vf(j,k,l+1,E_idx) + FL(j, k, l+1) ) + &
                             0.5_wp* (rho_R * (qR_rs_vf(j,k,l,momxe)**2.0) + &
                             qR_rs_vf(j,k,l,E_idx) + FR(j, k, l) ) + &
                             0.5_wp * cfl * (qR_rs_vf(j, k, l, momxe) - qL_rs_vf(j, k, l+1, momxe))

                         flux_vf(E_idx)%sf(j, k, l) = &
                            0.5_wp * (qL_rs_vf(j,k,l+1,momxe) * (qL_rs_vf(j,k,l+1,E_idx)*gamma_L + pi_inf_L + &
                            0.5_wp * rho_L * (qL_rs_vf(j, k, l+1,momxb)**2._wp + qL_rs_vf(j, k, l+1,momxb+1)**2._wp + qL_rs_vf(j,k,l+1,momxe)**2._wp) + &
                            qL_rs_vf(j,k,l+1,E_idx) + FL(j, k, l+1)) ) + &
                            0.5_wp * ( qR_rs_vf(j,k,l,momxe) * (qR_rs_vf(j,k,l,E_idx)*gamma_R + pi_inf_R + &
                            0.5_wp * rho_R * (qR_rs_vf(j, k, l,momxb)**2._wp + qR_rs_vf(j, k, l,momxb+1)**2._wp + qR_rs_vf(j,k,l,momxe)**2._wp ) + &
                            qR_rs_vf(j,k,l,E_idx) + FR(j, k, l)) ) + &
                            0.5_wp * cfl * (qR_rs_vf(j, k, l, E_idx) - qL_rs_vf(j, k, l+1, E_idx))

                        if(viscous) then
                            mu_L = 1/Re_L; mu_R = 1/Re_R

                            flux_vf(momxe)%sf(j, k, l) = flux_vf(momxe)%sf(j, k, l) - &
                                        0.5_wp * mu_L*((4._wp/3._wp)*dwLz(j, k, l+1) - (2._wp/3._wp)*duLx(j, k, l+1) - (2._wp/3._wp)*dvLy(j,k,l+1)) - &
                                        0.5_wp * mu_R*((4._wp/3._wp)*dwRz(j, k, l) - (2._wp/3._wp)*duRx(j, k, l) - (2._wp/3._wp)*dvRy(j,k,l) )

                            flux_vf(momxb)%sf(j, k, l) = flux_vf(momxb)%sf(j, k, l) - &
                                       0.5_wp*mu_L*(duLz(j, k, l+1) + dwLx(j, k, l+1))  -  &
                                       0.5_wp*mu_R*(duRz(j, k, l) + dwRx(j, k, l))

                            flux_vf(momxb+1)%sf(j, k, l) = flux_vf(momxb+1)%sf(j, k, l) - &
                                       0.5_wp*mu_L*(dvLz(j, k, l+1) + dwLy(j, k, l+1))  -  &
                                       0.5_wp*mu_R*(dvRz(j, k, l) + dwRy(j, k, l))

                            flux_vf(E_idx)%sf(j, k, l) = flux_vf(E_idx)%sf(j, k, l) - &
                            0.5_wp*mu_L*qL_rs_vf(j, k, l+1, momxe)*((4._wp/3._wp)*dwLz(j, k, l+1) - (2._wp/3._wp)*duLx(j, k, l+1) - (2._wp/3._wp)*dvLy(j ,k ,l+1)) - &
                            0.5_wp*mu_L*qL_rs_vf(j, k, l+1, momxe)*(duLz(j, k, l+1) + dwLx(j, k, l+1))   - &
                            0.5_wp*mu_L*qL_rs_vf(j, k, l+1, momxe)*(dvLz(j, k, l+1) + dwLy(j, k, l+1))   - &
                            0.5_wp*mu_R*qR_rs_vf(j, k, l, momxe)*((4._wp/3._wp)*dwRz(j, k, l) - (2._wp/3._wp)*duRx(j, k, l) - (2._wp/3._wp)*dvRy(j ,k ,l)) - &
                            0.5_wp*mu_R*qR_rs_vf(j, k, l, momxe)*(duRz(j, k, l) + dwRx(j, k, l))   - &
                            0.5_wp*mu_R*qR_rs_vf(j, k, l, momxe)*(dvRz(j, k, l) + dwRy(j, k, l))
                        end if
                    end do
                end do
            end do
        end if

    end subroutine s_igr_riemann_solver

    subroutine s_initialize_igr(q_prim_vf)

        type(scalar_field), dimension(sys_size), intent(in) :: q_prim_vf

        real(wp) :: rho_rx, rho_ry, rho_rz, rho_lx, rho_ly, rho_lz

        if (p == 0) then
            alf_igr = alf_factor*max(dx(1), dy(1))**2._wp
        else
            alf_igr = alf_factor*max(dx(1), dy(1), dz(1))**2._wp
        end if
        !$acc update device(alf_igr)

        omega = 1.0_wp
        !$acc update device(omega)

        !$acc parallel loop collapse(3) gang vector default(present)
        do l = idwbuff(3)%beg, idwbuff(3)%end
            do k = idwbuff(2)%beg, idwbuff(2)%end
                do j = idwbuff(1)%beg, idwbuff(1)%end
                    do i = 1, num_fluids
                        rho_igr(j,k,l) = q_prim_vf(contxb + i - 1)%sf(j,k,l)
                    end do
                end do
            end do
        end do

        if(p == 0) then
            !$acc parallel loop collapse(3) gang vector default(present)
            do l = 0, p
                do k = idwbuff(2)%beg + 2, idwbuff(2)%end - 2
                    do j = idwbuff(1)%beg+2, idwbuff(1)%end - 2
                        dux(j,k,l) = (1/(12._wp*dx(j))) * ( &
                            8._wp*q_prim_vf(momxb)%sf(j+1,k,l) - &
                            8._wp*q_prim_vf(momxb)%sf(j-1,k,l) + &
                            q_prim_vf(momxb)%sf(j-2,k,l) - &
                            q_prim_vf(momxb)%sf(j+2,k,l) )

                        duy(j,k,l) = (1/(12._wp*dy(k))) * ( &
                            8._wp*q_prim_vf(momxb)%sf(j,k+1,l) - &
                            8._wp*q_prim_vf(momxb)%sf(j,k-1,l) + &
                            q_prim_vf(momxb)%sf(j,k-2,l) - &
                            q_prim_vf(momxb)%sf(j,k+2,l) )

                    end do
                end do
            end do

            !$acc parallel loop collapse(3) gang vector default(present)
            do l = 0, p
                do k = idwbuff(2)%beg + 2, idwbuff(2)%end - 2
                    do j = idwbuff(1)%beg+2, idwbuff(1)%end - 2
                        dvx(j,k,l) = (1/(12._wp*dx(j))) * ( &
                            8._wp*q_prim_vf(momxb+1)%sf(j+1,k,l) - &
                            8._wp*q_prim_vf(momxb+1)%sf(j-1,k,l) + &
                            q_prim_vf(momxb+1)%sf(j-2,k,l) - &
                            q_prim_vf(momxb+1)%sf(j+2,k,l) )

                        dvy(j,k,l) = (1/(12._wp*dy(k))) * ( &
                            8._wp*q_prim_vf(momxb+1)%sf(j,k+1,l) - &
                            8._wp*q_prim_vf(momxb+1)%sf(j,k-1,l) + &
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
                                              (dux(j,k,l) + dvy(j,k,l))**2._wp)
                   end do
                end do
            end do

            !$acc parallel loop collapse(3) gang vector default(present) private(rho_lx, rho_rx, rho_ly, rho_ry)
            do l = 0, p
                do k = idwbuff(2)%beg + 1, idwbuff(2)%end - 1
                    do j = idwbuff(1)%beg + 1, idwbuff(1)%end - 1
                        rho_lx = 0.5_wp *(1._wp / rho_igr(j,k,l) + 1._wp / rho_igr(j-1,k,l))
                        rho_rx = 0.5_wp *(1._wp / rho_igr(j,k,l) + 1._wp / rho_igr(j+1,k,l))
                        rho_ly = 0.5_wp *(1._wp / rho_igr(j,k,l) + 1._wp / rho_igr(j,k-1,l))
                        rho_ry = 0.5_wp *(1._wp / rho_igr(j,k,l) + 1._wp / rho_igr(j,k+1,l))

                        fd_coeff(j,k,l) = 1._wp / rho_igr(j, k, l)

                        fd_coeff(j,k,l) = fd_coeff(j,k,l) + alf_igr * ( (1._wp / dx(j)**2._wp) * (rho_lx + rho_rx) +  (1._wp / dy(k)**2._wp) *(rho_ly + rho_ry) )
                    end do
                end do
            end do

        !!! 3D
        else
            !$acc parallel loop collapse(3) gang vector default(present)
            do l = idwbuff(3)%beg + 2, idwbuff(3)%end - 2
                do k = idwbuff(2)%beg + 2, idwbuff(2)%end - 2
                    do j = idwbuff(1)%beg+2, idwbuff(1)%end - 2
                        dux(j,k,l) = (1/(12._wp*dx(j))) * ( &
                            8._wp*q_prim_vf(momxb)%sf(j+1,k,l) - &
                            8._wp*q_prim_vf(momxb)%sf(j-1,k,l) + &
                            q_prim_vf(momxb)%sf(j-2,k,l) - &
                            q_prim_vf(momxb)%sf(j+2,k,l) )

                        duy(j,k,l) = (1/(12._wp*dy(k))) * ( &
                            8._wp*q_prim_vf(momxb)%sf(j,k+1,l) - &
                            8._wp*q_prim_vf(momxb)%sf(j,k-1,l) + &
                            q_prim_vf(momxb)%sf(j,k-2,l) - &
                            q_prim_vf(momxb)%sf(j,k+2,l) )

                        duz(j,k,l) = (1/(12._wp*dz(l))) * ( &
                            8._wp*q_prim_vf(momxb)%sf(j,k,l+1) - &
                            8._wp*q_prim_vf(momxb)%sf(j,k,l-1) + &
                            q_prim_vf(momxb)%sf(j,k,l-2) - &
                            q_prim_vf(momxb)%sf(j,k,l+2) )
                    end do
                end do
            end do

            !$acc parallel loop collapse(3) gang vector default(present)
            do l = idwbuff(3)%beg + 2, idwbuff(3)%end - 2
                do k = idwbuff(2)%beg + 2, idwbuff(2)%end - 2
                    do j = idwbuff(1)%beg+2, idwbuff(1)%end - 2
                        dvx(j,k,l) = (1/(12._wp*dx(j))) * ( &
                            8._wp*q_prim_vf(momxb+1)%sf(j+1,k,l) - &
                            8._wp*q_prim_vf(momxb+1)%sf(j-1,k,l) + &
                            q_prim_vf(momxb+1)%sf(j-2,k,l) - &
                            q_prim_vf(momxb+1)%sf(j+2,k,l) )

                        dvy(j,k,l) = (1/(12._wp*dy(k))) * ( &
                            8._wp*q_prim_vf(momxb+1)%sf(j,k+1,l) - &
                            8._wp*q_prim_vf(momxb+1)%sf(j,k-1,l) + &
                            q_prim_vf(momxb+1)%sf(j,k-2,l) - &
                            q_prim_vf(momxb+1)%sf(j,k+2,l) )

                        dvz(j,k,l) = (1/(12._wp*dz(l))) * ( &
                            8._wp*q_prim_vf(momxb+1)%sf(j,k,l+1) - &
                            8._wp*q_prim_vf(momxb+1)%sf(j,k,l-1) + &
                            q_prim_vf(momxb+1)%sf(j,k,l-2) - &
                            q_prim_vf(momxb+1)%sf(j,k,l+2) )
                    end do
                end do
            end do

            !$acc parallel loop collapse(3) gang vector default(present)
            do l = idwbuff(3)%beg + 2, idwbuff(3)%end - 2
                do k = idwbuff(2)%beg + 2, idwbuff(2)%end - 2
                    do j = idwbuff(1)%beg+2, idwbuff(1)%end - 2
                        dwx(j,k,l) = (1/(12._wp*dx(j))) * ( &
                            8._wp*q_prim_vf(momxe)%sf(j+1,k,l) - &
                            8._wp*q_prim_vf(momxe)%sf(j-1,k,l) + &
                            q_prim_vf(momxe)%sf(j-2,k,l) - &
                            q_prim_vf(momxe)%sf(j+2,k,l) )

                        dwy(j,k,l) = (1/(12._wp*dy(k))) * ( &
                            8._wp*q_prim_vf(momxe)%sf(j,k+1,l) - &
                            8._wp*q_prim_vf(momxe)%sf(j,k-1,l) + &
                            q_prim_vf(momxe)%sf(j,k-2,l) - &
                            q_prim_vf(momxe)%sf(j,k+2,l) )

                        dwz(j,k,l) = (1/(12._wp*dz(l))) * ( &
                            8._wp*q_prim_vf(momxe)%sf(j,k,l+1) - &
                            8._wp*q_prim_vf(momxe)%sf(j,k,l-1) + &
                            q_prim_vf(momxe)%sf(j,k,l-2) - &
                            q_prim_vf(momxe)%sf(j,k,l+2) )
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
                                              (dux(j,k,l) + dvy(j,k,l) + dwz(j,k,l))**2._wp )
                   end do
                end do
            end do

            !$acc parallel loop collapse(3) gang vector default(present) private(rho_lx, rho_rx, rho_ly, rho_ry, rho_lz, rho_rz)
            do l = idwbuff(3)%beg + 1, idwbuff(3)%end - 1
                do k = idwbuff(2)%beg + 1, idwbuff(2)%end - 1
                    do j = idwbuff(1)%beg + 1, idwbuff(1)%end - 1
                        rho_lx = 0.5_wp *(1._wp / rho_igr(j,k,l) + 1._wp / rho_igr(j-1,k,l))
                        rho_rx = 0.5_wp *(1._wp / rho_igr(j,k,l) + 1._wp / rho_igr(j+1,k,l))
                        rho_ly = 0.5_wp *(1._wp / rho_igr(j,k,l) + 1._wp / rho_igr(j,k-1,l))
                        rho_ry = 0.5_wp *(1._wp / rho_igr(j,k,l) + 1._wp / rho_igr(j,k+1,l))
                        rho_lz = 0.5_wp *(1._wp / rho_igr(j,k,l) + 1._wp / rho_igr(j,k,l-1))
                        rho_rz = 0.5_wp *(1._wp / rho_igr(j,k,l) + 1._wp / rho_igr(j,k,l+1))

                        fd_coeff(j,k,l) = 1._wp / rho_igr(j, k, l)

                        fd_coeff(j,k,l) = fd_coeff(j,k,l) + alf_igr * ( (1._wp / dx(j)**2._wp) * (rho_lx + rho_rx) +  (1._wp / dy(k)**2._wp) *(rho_ly + rho_ry) + (1._wp / dz(l)**2._wp) * (rho_lz + rho_rz) )

                    end do
                end do
            end do
        end if

    end subroutine s_initialize_igr

    subroutine s_reconstruct_prim_vars_igr(q_prim_vf, idir)

        type(scalar_field), intent(in), dimension(sys_size) :: q_prim_vf
        integer, intent(in) :: idir

        do i = 1,sys_size
            call s_reconstruct_igr(qL_rs_vf(:,:,:,i), qR_rs_vf(:,:,:,i), q_prim_vf(i)%sf, idir)
        end do

        call s_reconstruct_igr(FL, FR, jac, idir)

    end subroutine s_reconstruct_prim_vars_igr

    subroutine s_get_viscous_igr(idir)

        integer, intent(in) :: idir

        if (idir == 1) then
            call s_reconstruct_igr(duLx, duRx, dux, 1)
            call s_reconstruct_igr(dvLx, dvRx, dvx, 1)

            call s_reconstruct_igr(duLy, duRy, duy, 1)
            call s_reconstruct_igr(dvLy, dvRy, dvy, 1)

            if (p > 0) then
                call s_reconstruct_igr(dwLx, dwRx, dwx, 1)
                call s_reconstruct_igr(dwLy, dwRy, dwy, 1)

                call s_reconstruct_igr(duLz, duRz, duz, 1)
                call s_reconstruct_igr(dvLz, dvRz, dvz, 1)
                call s_reconstruct_igr(dwLz, dwRz, dwz, 1)
            end if
        elseif (idir == 2) then
            call s_reconstruct_igr(duLx, duRx, dux, 2)
            call s_reconstruct_igr(dvLx, dvRx, dvx, 2)

            call s_reconstruct_igr(duLy, duRy, duy, 2)
            call s_reconstruct_igr(dvLy, dvRy, dvy, 2)

            if (p > 0) then
                call s_reconstruct_igr(dwLx, dwRx, dwx, 2)
                call s_reconstruct_igr(dwLy, dwRy, dwy, 2)

                call s_reconstruct_igr(duLz, duRz, duz, 2)
                call s_reconstruct_igr(dvLz, dvRz, dvz, 2)
                call s_reconstruct_igr(dwLz, dwRz, dwz, 2)
            end if
        elseif (idir == 3) then
            call s_reconstruct_igr(duLx, duRx, dux, 3)
            call s_reconstruct_igr(dvLx, dvRx, dvx, 3)
            call s_reconstruct_igr(dwLx, dwRx, dwx, 3)

            call s_reconstruct_igr(duLy, duRy, duy, 3)
            call s_reconstruct_igr(dvLy, dvRy, dvy, 3)
            call s_reconstruct_igr(dwLy, dwRy, dwy, 3)

            call s_reconstruct_igr(duLz, duRz, duz, 3)
            call s_reconstruct_igr(dvLz, dvRz, dvz, 3)
            call s_reconstruct_igr(dwLz, dwRz, dwz, 3)
        end if

    end subroutine s_get_viscous_igr

    subroutine s_igr_flux_add(q_prim_vf, rhs_vf, flux_vf, idir)

        type(scalar_field), &
            dimension(sys_size), &
            intent(inout) :: q_prim_vf, flux_vf, rhs_vf

        integer, intent(in) :: idir

        if (idir == 1) then
            !$acc parallel loop collapse(4) gang vector default(present)
            do i = 1, sys_size
                do l = 0, p
                    do k = 0, n
                        do j = 0, m
                            rhs_vf(i)%sf(j, k, l) = 1._wp/dx(j)* &
                                (flux_vf(i)%sf(j - 1 , k, l) &
                                 - flux_vf(i)%sf(j , k, l))
                        end do
                    end do
                end do
            end do

            !$acc parallel loop collapse(4) gang vector default(present)
            do i = advxb, advxe
                do l = 0, p
                    do k = 0, n
                        do j = 0, m
                            !rhs_vf(i)%sf(j, k, l) = rhs_vf(i)%sf(j, k, l) + &
                                !q_prim_vf(i)%sf(j,k,l) * dux(j,k,l)
                            rhs_vf(i)%sf(j, k, l) = rhs_vf(i)%sf(j, k, l) + &
                                q_prim_vf(i)%sf(j,k,l) * 1/dx(j)*(qR_rs_vf(j,k,l,momxb) - qL_rs_vf(j,k,l,momxb))
                            !rhs_vf(i)%sf(j,k,l) = &
                                !q_prim_vf(momxb)%sf(j,k,l) * &
                                !(1/(2*dx(j))) * (qL_rs_vf(j,k,l,i) + qR_rs_vf(j-1,k,l,i) - &
                                !qR_rs_vf(j,k,l,i) - qL_rs_vf(j+1,k,l,i))
                        end do
                    end do
                end do
            end do
        elseif (idir == 2) then
            !$acc parallel loop collapse(4) gang vector default(present)
            do i = 1, sys_size
                do l = 0, p
                    do k = 0, n
                        do j = 0, m
                            rhs_vf(i)%sf(j, k, l) = &
                                rhs_vf(i)%sf(j, k, l) + 1._wp/dy(k)* &
                                (flux_vf(i)%sf(j, k - 1, l) &
                                 - flux_vf(i)%sf(j, k , l))
                        end do
                    end do
                end do
            end do

            !$acc parallel loop collapse(4) gang vector default(present)
            do i = advxb, advxe
                do l = 0, p
                    do k = 0, n
                        do j = 0, m
                            !rhs_vf(i)%sf(j, k, l) = rhs_vf(i)%sf(j, k, l) + &
                                !q_prim_vf(i)%sf(j,k,l) * dvy(j,k,l)
                            rhs_vf(i)%sf(j, k, l) = rhs_vf(i)%sf(j, k, l) + &
                                q_prim_vf(i)%sf(j,k,l) * (1/dy(k))*(qR_rs_vf(j,k,l,momxb+1) - qL_rs_vf(j,k,l,momxb+1))
                            !rhs_vf(i)%sf(j,k,l) = rhs_vf(i)%sf(j,k,l) + &
                                !q_prim_vf(momxb+1)%sf(j,k,l) * &
                                !(1/(2*dy(k))) * (qL_rs_vf(j,k,l,i) + qR_rs_vf(j,k-1,l,i) - &
                                !qR_rs_vf(j,k,l,i) - qL_rs_vf(j,k+1,l,i))
                        end do
                    end do
                end do
            end do
        elseif (idir == 3) then
            !$acc parallel loop collapse(4) gang vector default(present)
            do i = 1, E_idx
                do l = 0, p
                    do k = 0, n
                        do j = 0, m
                            rhs_vf(i)%sf(j, k, l) = &
                                rhs_vf(i)%sf(j, k, l) + 1._wp/dz(l)* &
                                (flux_vf(i)%sf(j, k, l-1) &
                                 - flux_vf(i)%sf(j, k, l ))
                        end do
                    end do
                end do
            end do

            !$acc parallel loop collapse(4) gang vector default(present)
            do i = advxb, advxe
                do l = 0, p
                    do k = 0, n
                        do j = 0, m
                            !rhs_vf(i)%sf(j, k, l) = rhs_vf(i)%sf(j, k, l) + &
                                !q_prim_vf(i)%sf(j,k,l) * dwz(j,k,l)
                            rhs_vf(i)%sf(j, k, l) = rhs_vf(i)%sf(j, k, l) + &
                                q_prim_vf(i)%sf(j,k,l) * (1/dz(l))*(qR_rs_vf(j,k,l,momxe) - qL_rs_vf(j,k,l,momxe))
                            !rhs_vf(i)%sf(j,k,l) = rhs_vf(i)%sf(j,k,l) - &
                                !q_prim_vf(momxe)%sf(j,k,l) * &
                                !(1/(12._wp*dz(l))) * ( &
                                !8._wp*q_prim_vf(i)%sf(j,k,l+1) - &
                                !8._wp*q_prim_vf(i)%sf(j,k,l-1) + &
                                !q_prim_vf(i)%sf(j,k,l-2) - &
                                !q_prim_vf(i)%sf(j,k,l+1) )
                        end do
                    end do
                end do
            end do
        end if

    end subroutine s_igr_flux_add

end module m_igr
