#:include 'macros.fpp'

module m_igr

    use m_derived_types        !< Definitions of the derived types

    use m_global_parameters

    use m_variables_conversion

    use m_mpi_proxy

    use m_helper

    implicit none

    private; public :: s_initialize_igr_module, &
        s_igr_jacobi_iteration, &
        s_igr_riemann_solver, &
        s_igr_sigma, &
        s_initialize_igr, &
        s_igr_flux_add

    real(wp), allocatable, dimension(:, :, :) :: fd_coeff
    real(wp), allocatable, dimension(:, :, :) :: jac,jac_rhs,rho_igr
    !$acc declare create(fd_coeff,jac, jac_rhs,rho_igr)

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

        if (viscous) then
            @:ALLOCATE(Res(1:2, 1:maxval(Re_size)))
            do i = 1, 2
                do j = 1, Re_size(i)
                    Res(i, j) = fluid_pp(Re_idx(i, j))%Re(i)
                end do
            end do
            !$acc update device(Res, Re_idx, Re_size)
        end if

        if(igr) then 
            if(num_fluids > 1) then 
                 @:ALLOCATE(rho_igr(idwbuff(1)%beg:idwbuff(1)%end, &
                    idwbuff(2)%beg:idwbuff(2)%end, &
                    idwbuff(3)%beg:idwbuff(3)%end))
            end if
            #:for VAR in [ 'jac','jac_rhs','fd_coeff']
                @:ALLOCATE(${VAR}$(idwbuff(1)%beg:idwbuff(1)%end, &
                             idwbuff(2)%beg:idwbuff(2)%end, &
                             idwbuff(3)%beg:idwbuff(3)%end))
            #:endfor

            !$acc parallel loop collapse(3) gang vector default(present)
            do l = idwbuff(3)%beg, idwbuff(3)%end
                do k = idwbuff(2)%beg, idwbuff(2)%end
                    do j = idwbuff(1)%beg, idwbuff(1)%end
                        jac(j, k, l) = 0._wp
                   end do
                end do
            end do
        end if

    end subroutine s_initialize_igr_module

    subroutine s_igr_jacobi_iteration(q_prim_vf, t_step)

        type(scalar_field), &
            dimension(sys_size), &
            intent(inout) :: q_prim_vf
        real(wp) :: rho_rx, rho_ry, rho_rz, rho_lx, rho_ly, rho_lz
        real(wp) :: resid
        integer :: num_iters, t_step

        num_iters = num_igr_iters

        if(num_fluids > 1) then 
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

                        fd_coeff(j,k,l) = (1._wp / rho_igr(j,k,l))

                        fd_coeff(j,k,l) = fd_coeff(j,k,l) + alf_igr * ( (1._wp / dx(j)**2._wp) * (rho_lx + rho_rx) +  (1._wp / dy(k)**2._wp) *(rho_ly + rho_ry) )

                        if(p > 0) then 
                            fd_coeff(j,k,l) = fd_coeff(j,k,l) + alf_igr* (1._wp / dz(l)**2._wp) * (rho_lz + rho_rz)
                        end if
                    end do 
                end do 
            end do
        else 
            !$acc parallel loop collapse(3) gang vector default(present) private(rho_lx, rho_rx, rho_ly, rho_ry, rho_lz, rho_rz)
            do l = 0, p
                do k = 0, n
                    do j = 0, m
                        rho_lx = 0.5_wp *(1._wp / q_prim_vf(1)%sf(j,k,l) + 1._wp / q_prim_vf(1)%sf(j-1,k,l))
                        rho_rx = 0.5_wp *(1._wp / q_prim_vf(1)%sf(j,k,l) + 1._wp / q_prim_vf(1)%sf(j+1,k,l))
                        rho_ly = 0.5_wp *(1._wp / q_prim_vf(1)%sf(j,k,l) + 1._wp / q_prim_vf(1)%sf(j,k-1,l))
                        rho_ry = 0.5_wp *(1._wp / q_prim_vf(1)%sf(j,k,l) + 1._wp / q_prim_vf(1)%sf(j,k+1,l))

                        if(p > 0) then
                            rho_lz = 0.5_wp *(1._wp / q_prim_vf(1)%sf(j,k,l) + 1._wp / q_prim_vf(1)%sf(j,k,l-1))
                            rho_rz = 0.5_wp *(1._wp / q_prim_vf(1)%sf(j,k,l) + 1._wp / q_prim_vf(1)%sf(j,k,l+1))
                        end if

                        fd_coeff(j,k,l) = (1._wp / q_prim_vf(1)%sf(j,k,l))

                        fd_coeff(j,k,l) = fd_coeff(j,k,l) + alf_igr * ( (1._wp / dx(j)**2._wp) * (rho_lx + rho_rx) +  (1._wp / dy(k)**2._wp) *(rho_ly + rho_ry) )

                        if(p > 0) then 
                            fd_coeff(j,k,l) = fd_coeff(j,k,l) + alf_igr* (1._wp / dz(l)**2._wp) * (rho_lz + rho_rz)
                        end if
                    end do 
                end do 
            end do
        end if

        do q = 1, num_iters
            if(num_fluids > 1) then 
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

                            if(p > 0) then
                                jac(j, k, l) = (alf_igr / fd_coeff(j,k,l)) * ( (1._wp / dx(j)**2._wp) * (rho_lx* jac(j-1,k,l) + rho_rx*jac(j+1,k,l)) + &
                                                        (1._wp / dy(k)**2._wp) * (rho_ly* jac(j,k-1,l) + rho_ry*jac(j,k+1,l)) + &
                                                        (1._wp / dz(l)**2._wp) * (rho_lz* jac(j,k,l-1) + rho_rz*jac(j,k,l+1)) ) + &
                                                        jac_rhs(j,k,l) / fd_coeff(j, k, l)
                           else 
                                jac(j, k, l) = (alf_igr / fd_coeff(j,k,l)) * ( (1._wp / dx(j)**2._wp) * (rho_lx* jac(j-1,k,l) + rho_rx*jac(j+1,k,l)) + &
                                                        (1._wp / dy(k)**2._wp) * (rho_ly* jac(j,k-1,l) + rho_ry*jac(j,k+1,l)) ) + &
                                                        jac_rhs(j,k,l) / fd_coeff(j, k, l)                              
                           end if
                        end do
                    end do
                end do
            else 
                !$acc parallel loop collapse(3) gang vector default(present) private(rho_lx, rho_rx, rho_ly, rho_ry, rho_lz, rho_rz)
                do l = 0, p
                    do k = 0, n
                        do j = 0, m
                            rho_lx = 0.5_wp *(1._wp / q_prim_vf(1)%sf(j,k,l) + 1._wp / q_prim_vf(1)%sf(j-1,k,l))
                            rho_rx = 0.5_wp *(1._wp / q_prim_vf(1)%sf(j,k,l) + 1._wp / q_prim_vf(1)%sf(j+1,k,l))
                            rho_ly = 0.5_wp *(1._wp / q_prim_vf(1)%sf(j,k,l) + 1._wp / q_prim_vf(1)%sf(j,k-1,l))
                            rho_ry = 0.5_wp *(1._wp / q_prim_vf(1)%sf(j,k,l) + 1._wp / q_prim_vf(1)%sf(j,k+1,l))

                            if(p > 0) then
                                rho_lz = 0.5_wp *(1._wp / q_prim_vf(1)%sf(j,k,l) + 1._wp / q_prim_vf(1)%sf(j,k,l-1))
                                rho_rz = 0.5_wp *(1._wp / q_prim_vf(1)%sf(j,k,l) + 1._wp / q_prim_vf(1)%sf(j,k,l+1))
                            end if

                            if(p > 0) then
                                jac(j, k, l) = (alf_igr / fd_coeff(j,k,l)) * ( (1._wp / dx(j)**2._wp) * (rho_lx* jac(j-1,k,l) + rho_rx*jac(j+1,k,l)) + &
                                                        (1._wp / dy(k)**2._wp) * (rho_ly* jac(j,k-1,l) + rho_ry*jac(j,k+1,l)) + &
                                                        (1._wp / dz(l)**2._wp) * (rho_lz* jac(j,k,l-1) + rho_rz*jac(j,k,l+1)) ) + &
                                                        jac_rhs(j,k,l) / fd_coeff(j, k, l) 
                           else 
                                jac(j, k, l) = (alf_igr / fd_coeff(j,k,l)) * ( (1._wp / dx(j)**2._wp) * (rho_lx* jac(j-1,k,l) + rho_rx*jac(j+1,k,l)) + &
                                                        (1._wp / dy(k)**2._wp) * (rho_ly* jac(j,k-1,l) + rho_ry*jac(j,k+1,l)) ) + &
                                                        jac_rhs(j,k,l) / fd_coeff(j, k, l)                               
                           end if
                        end do
                    end do
                end do
            end if

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
        end do
    end subroutine s_igr_jacobi_iteration

    subroutine s_igr_sigma(q_prim_vf, rhs_vf, idir)

        type(scalar_field), &
            dimension(sys_size), &
            intent(inout) :: rhs_vf
        type(scalar_field), &
            dimension(sys_size), &
            intent(inout) :: q_prim_vf
        integer, intent(in) :: idir

        real(wp) :: F_L, vel_L,rho_L
        real(wp), dimension(num_fluids) :: alpha_rho_L

        if(idir == 1) then 
            if(p == 0) then 
                !$acc parallel loop collapse(3) gang vector default(present) private(F_L,vel_L,alpha_rho_L)
                do l = 0, p
                    do k = 0, n
                        do j = -buff_size+3, m+buff_size-3

                            F_L = (1._wp/60._wp) * (-3._wp * jac(j-1, k, l) + &
                                                27._wp * jac(j, k, l) + &
                                                47._wp * jac(j+1, k, l) -   &
                                                13._wp * jac(j+2, k, l) + &
                                                2._wp * jac(j+3, k, l))

                            !$acc loop seq 
                            do i = 1, num_fluids
                                alpha_rho_L(i) = (1._wp/60._wp) * (-3._wp * q_prim_vf(i)%sf(j-1, k, l) + &
                                                27._wp * q_prim_vf(i)%sf(j, k, l) + &
                                                47._wp * q_prim_vf(i)%sf(j+1, k, l) -   &
                                                13._wp * q_prim_vf(i)%sf(j+2, k, l) + &
                                                2._wp * q_prim_vf(i)%sf(j+3, k, l))
                            end do

                            vel_L = (1._wp/60._wp) * (-3._wp * q_prim_vf(momxb)%sf(j-1,k,l) + &
                                                27._wp * q_prim_vf(momxb)%sf(j,k,l) + &
                                                47._wp * q_prim_vf(momxb)%sf(j+1,k,l) -   &
                                                13._wp * q_prim_vf(momxb)%sf(j+2,k,l) + &
                                                2._wp * q_prim_vf(momxb)%sf(j+3,k,l)) / sum(alpha_rho_L)

                            !$acc atomic
                            rhs_vf(momxb)%sf(j+1,k,l) = rhs_vf(momxb)%sf(j+1,k,l) + &
                                                      0.5_wp * F_L * (1._wp/dx(j+1))

                            !$acc atomic
                            rhs_vf(E_idx)%sf(j+1,k,l) = rhs_vf(E_idx)%sf(j+1,k,l) + &
                                                      0.5_wp * vel_L * F_L * (1._wp/dx(j+1))


                            !$acc atomic
                            rhs_vf(momxb)%sf(j,k,l) = rhs_vf(momxb)%sf(j,k,l) - &
                                                      0.5_wp * F_L * (1._wp/dx(j))

                            !$acc atomic
                            rhs_vf(E_idx)%sf(j,k,l) = rhs_vf(E_idx)%sf(j,k,l) - &
                                                      0.5_wp * vel_L * F_L * (1._wp/dx(j))


                            !$acc loop seq 
                            do i = 1, num_fluids
                                alpha_rho_L(i) = (1._wp/60._wp) * (-3._wp * q_prim_vf(i)%sf(j+2, k, l) + &
                                                27._wp * q_prim_vf(i)%sf(j+1, k, l) + &
                                                47._wp * q_prim_vf(i)%sf(j, k, l) -   &
                                                13._wp * q_prim_vf(i)%sf(j-1, k, l) + &
                                                2._wp * q_prim_vf(i)%sf(j-2, k, l))
                            end do

                            vel_L = (1._wp/60._wp) * (-3._wp * q_prim_vf(momxb)%sf(j+2,k,l) + &
                                                27._wp * q_prim_vf(momxb)%sf(j+1,k,l) + &
                                                47._wp * q_prim_vf(momxb)%sf(j,k,l) -   &
                                                13._wp * q_prim_vf(momxb)%sf(j-1,k,l) + &
                                                2._wp * q_prim_vf(momxb)%sf(j-2,k,l)) / sum(alpha_rho_L)

                            !$acc atomic
                            rhs_vf(momxb)%sf(j+1,k,l) = rhs_vf(momxb)%sf(j+1,k,l) + &
                                                      0.5_wp * F_L * (1._wp/dx(j+1))

                            !$acc atomic
                            rhs_vf(E_idx)%sf(j+1,k,l) = rhs_vf(E_idx)%sf(j+1,k,l) + &
                                                      0.5_wp * vel_L * F_L * (1._wp/dx(j+1))


                            !$acc atomic
                            rhs_vf(momxb)%sf(j,k,l) = rhs_vf(momxb)%sf(j,k,l) - &
                                                      0.5_wp * F_L * (1._wp/dx(j))

                            !$acc atomic
                            rhs_vf(E_idx)%sf(j,k,l) = rhs_vf(E_idx)%sf(j,k,l) - &
                                                      0.5_wp * vel_L * F_L * (1._wp/dx(j))                   

                        end do 
                    end do 
                end do
            else 
                !$acc parallel loop collapse(3) gang vector default(present) private(F_L, vel_L,alpha_rho_L)
                do l = 0, p
                    do k = 0, n
                        do j = -buff_size+3, m+buff_size-3

                           F_L = (1._wp/60._wp) * (-3._wp * jac(j-1, k, l) + &
                                                27._wp * jac(j, k, l) + &
                                                47._wp * jac(j+1, k, l) -   &
                                                13._wp * jac(j+2, k, l) + &
                                                2._wp * jac(j+3, k, l))

                            !$acc loop seq 
                            do i = 1, num_fluids
                                alpha_rho_L(i) = (1._wp/60._wp) * (-3._wp * q_prim_vf(i)%sf(j-1, k, l) + &
                                                27._wp * q_prim_vf(i)%sf(j, k, l) + &
                                                47._wp * q_prim_vf(i)%sf(j+1, k, l) -   &
                                                13._wp * q_prim_vf(i)%sf(j+2, k, l) + &
                                                2._wp * q_prim_vf(i)%sf(j+3, k, l))
                            end do

                            vel_L = (1._wp/60._wp) * (-3._wp * q_prim_vf(momxb)%sf(j-1,k,l) + &
                                                27._wp * q_prim_vf(momxb)%sf(j,k,l) + &
                                                47._wp * q_prim_vf(momxb)%sf(j+1,k,l) -   &
                                                13._wp * q_prim_vf(momxb)%sf(j+2,k,l) + &
                                                2._wp * q_prim_vf(momxb)%sf(j+3,k,l)) / sum(alpha_rho_L)

                            !$acc atomic
                            rhs_vf(momxb)%sf(j+1,k,l) = rhs_vf(momxb)%sf(j+1,k,l) + &
                                                      0.5_wp * F_L * (1._wp/dx(j+1))

                            !$acc atomic
                            rhs_vf(E_idx)%sf(j+1,k,l) = rhs_vf(E_idx)%sf(j+1,k,l) + &
                                                      0.5_wp * vel_L * F_L * (1._wp/dx(j+1))


                            !$acc atomic
                            rhs_vf(momxb)%sf(j,k,l) = rhs_vf(momxb)%sf(j,k,l) - &
                                                      0.5_wp * F_L * (1._wp/dx(j))

                            !$acc atomic
                            rhs_vf(E_idx)%sf(j,k,l) = rhs_vf(E_idx)%sf(j,k,l) - &
                                                      0.5_wp * vel_L * F_L * (1._wp/dx(j))


                            !$acc loop seq 
                            do i = 1, num_fluids
                                alpha_rho_L(i) = (1._wp/60._wp) * (-3._wp * q_prim_vf(i)%sf(j+2, k, l) + &
                                                27._wp * q_prim_vf(i)%sf(j+1, k, l) + &
                                                47._wp * q_prim_vf(i)%sf(j, k, l) -   &
                                                13._wp * q_prim_vf(i)%sf(j-1, k, l) + &
                                                2._wp * q_prim_vf(i)%sf(j-2, k, l))
                            end do

                            vel_L = (1._wp/60._wp) * (-3._wp * q_prim_vf(momxb)%sf(j+2,k,l) + &
                                                27._wp * q_prim_vf(momxb)%sf(j+1,k,l) + &
                                                47._wp * q_prim_vf(momxb)%sf(j,k,l) -   &
                                                13._wp * q_prim_vf(momxb)%sf(j-1,k,l) + &
                                                2._wp * q_prim_vf(momxb)%sf(j-2,k,l)) / sum(alpha_rho_L)

                            !$acc atomic
                            rhs_vf(momxb)%sf(j+1,k,l) = rhs_vf(momxb)%sf(j+1,k,l) + &
                                                      0.5_wp * F_L * (1._wp/dx(j+1))

                            !$acc atomic
                            rhs_vf(E_idx)%sf(j+1,k,l) = rhs_vf(E_idx)%sf(j+1,k,l) + &
                                                      0.5_wp * vel_L * F_L * (1._wp/dx(j+1))


                            !$acc atomic
                            rhs_vf(momxb)%sf(j,k,l) = rhs_vf(momxb)%sf(j,k,l) - &
                                                      0.5_wp * F_L * (1._wp/dx(j))

                            !$acc atomic
                            rhs_vf(E_idx)%sf(j,k,l) = rhs_vf(E_idx)%sf(j,k,l) - &
                                                      0.5_wp * vel_L * F_L * (1._wp/dx(j)) 

                        end do 
                    end do 
                end do                
            end if
        end if

    end subroutine s_igr_sigma

    subroutine s_igr_riemann_solver(q_prim_vf, rhs_vf, idir)

        type(scalar_field), &
            dimension(sys_size), &
            intent(inout) :: rhs_vf
        type(scalar_field), &
            dimension(sys_size), &
            intent(inout) :: q_prim_vf
        integer, intent(in) :: idir

        real(wp) :: rho_L, gamma_L, pi_inf_L, a_L, cfl, pres_L, dvel1, dvel2, F_L, E_L, pres_R, a_R, rho_R, gamma_R, pi_inf_R, E_R
        real(wp), dimension(num_fluids) :: alpha_rho_L, alpha_L, alpha_R, alpha_rho_R
        real(wp), dimension(num_dims,-3:2) :: vel_L, vel_R
        real(wp), dimension(-2:2,num_dims) :: rho_sf
        real(wp), dimension(-3:2) :: mu_L, mu_R

        if (idir == 1) then
            if(p == 0) then
                !$acc parallel loop collapse(3) gang vector default(present) private(rho_L,gamma_L,pi_inf_L,mu_L,a_L,vel_L,vel_R, pres_L,alpha_L,alpha_R,alpha_rho_L,cfl,dvel1,dvel2,F_L,E_L,mu_R,rho_sf,alpha_rho_R)
                do l = 0, p
                    do k = 0, n
                        do j = -buff_size+5, m+buff_size-5

                            !$acc loop seq 
                            do i = 1, num_fluids
                                alpha_rho_L(i) = (1._wp/60._wp) * (-3._wp * q_prim_vf(i)%sf(j-1, k, l) + &
                                                27._wp * q_prim_vf(i)%sf(j, k, l) + &
                                                47._wp * q_prim_vf(i)%sf(j+1, k, l) -   &
                                                13._wp * q_prim_vf(i)%sf(j+2, k, l) + &
                                                2._wp * q_prim_vf(i)%sf(j+3, k, l))

                                alpha_rho_R(i) = (1._wp/60._wp) * (-3._wp * q_prim_vf(i)%sf(j+2, k, l) + &
                                                27._wp * q_prim_vf(i)%sf(j+1, k, l) + &
                                                47._wp * q_prim_vf(i)%sf(j, k, l) -   &
                                                13._wp * q_prim_vf(i)%sf(j-1, k, l) + &
                                                2._wp * q_prim_vf(i)%sf(j-2, k, l))

                                alpha_L(i) =  (1._wp/60._wp) * (-3._wp * q_prim_vf(E_idx+i)%sf(j-1, k, l) + &
                                                27._wp * q_prim_vf(E_idx+i)%sf(j, k, l) + &
                                                47._wp * q_prim_vf(E_idx+i)%sf(j+1, k, l) -   &
                                                13._wp * q_prim_vf(E_idx+i)%sf(j+2, k, l) + &
                                                2._wp * q_prim_vf(E_idx+i)%sf(j+3, k, l))
                            end do 

                            if(num_fluids > 1) then  
                                !$acc loop seq 
                                do i = 1,num_fluids                           
                                    alpha_R(i) = (1._wp/60._wp) * (-3._wp * q_prim_vf(E_idx+i)%sf(j+2, k, l) + &
                                                        27._wp * q_prim_vf(E_idx+i)%sf(j+1, k, l) + &
                                                        47._wp * q_prim_vf(E_idx+i)%sf(j, k, l) -   &
                                                        13._wp * q_prim_vf(E_idx+i)%sf(j-1, k, l) + &
                                                        2._wp * q_prim_vf(E_idx+i)%sf(j-2, k, l))
                                end do 
                            else 
                                alpha_R(1) = 1._wp
                            end if

                            !$acc loop seq 
                            do i = -2, 2
                                rho_L = 0._wp
                                !$acc loop seq 
                                do q = 1, num_fluids
                                    rho_L = rho_L + q_prim_vf(q)%sf(j+i,k,l)
                                end do
                                rho_sf(i,1) = rho_L
                            end do

                            !$acc loop seq 
                            do i = -2, 2
                                rho_L = 0._wp
                                !$acc loop seq 
                                do q = 1, num_fluids
                                    rho_L = rho_L + q_prim_vf(q)%sf(j,k+i,l)
                                end do
                                rho_sf(i,2) = rho_L
                            end do

                            rho_L = sum(alpha_rho_L)
                            gamma_L = sum(alpha_L*gammas)
                            pi_inf_L = sum(alpha_L*pi_infs)

                            rho_R = sum(alpha_rho_R)
                            gamma_R = sum(alpha_R*gammas)
                            pi_inf_R = sum(alpha_R*pi_infs)

                            !$acc loop seq 
                            do q = -3, 2
                                !$acc loop seq 
                                do i = 1, num_dims
                                    vel_L(i,q) = (1._wp/60._wp) * (-3._wp * q_prim_vf(momxb+i-1)%sf(j-1+q,k,l) + &
                                                    27._wp * q_prim_vf(momxb+i-1)%sf(j+q,k,l) + &
                                                    47._wp * q_prim_vf(momxb+i-1)%sf(j+1+q,k,l) -   &
                                                    13._wp * q_prim_vf(momxb+i-1)%sf(j+2+q,k,l) + &
                                                    2._wp * q_prim_vf(momxb+i-1)%sf(j+3+q,k,l)) / rho_L
                                end do 

                                !$acc loop seq 
                                do i = 1, num_dims
                                    vel_R(i,q) = (1._wp/60._wp) * (-3._wp * q_prim_vf(momxb+i-1)%sf(j+2+q,k,l) + &
                                                    27._wp * q_prim_vf(momxb+i-1)%sf(j+1+q,k,l) + &
                                                    47._wp * q_prim_vf(momxb+i-1)%sf(j+q,k,l) -   &
                                                    13._wp * q_prim_vf(momxb+i-1)%sf(j-1+q,k,l) + &
                                                    2._wp * q_prim_vf(momxb+i-1)%sf(j-2+q,k,l)) / rho_R
                                end do
                            end do

                            E_L = (1._wp/60._wp) * (-3._wp * q_prim_vf(E_idx)%sf(j-1,k,l) + &
                                                27._wp * q_prim_vf(E_idx)%sf(j,k,l) + &
                                                47._wp * q_prim_vf(E_idx)%sf(j+1,k,l) -   &
                                                13._wp * q_prim_vf(E_idx)%sf(j+2,k,l) + &
                                                2._wp * q_prim_vf(E_idx)%sf(j+3,k,l))

                            E_R = (1._wp/60._wp) * (-3._wp * q_prim_vf(E_idx)%sf(j+2,k,l) + &
                                                27._wp * q_prim_vf(E_idx)%sf(j+1,k,l) + &
                                                47._wp * q_prim_vf(E_idx)%sf(j,k,l) -   &
                                                13._wp * q_prim_vf(E_idx)%sf(j-1,k,l) + &
                                                2._wp * q_prim_vf(E_idx)%sf(j-2,k,l))       

                            pres_L = (E_L - pi_inf_L - 0.5_wp*rho_L*(vel_L(1,0)**2._wp + vel_L(2,0)**2._wp))/gamma_L                    

                            pres_R = (E_R - pi_inf_R - 0.5_wp*rho_R*(vel_R(1,0)**2._wp + vel_R(2,0)**2._wp))/gamma_R 

                            a_L = sqrt((pres_L*(1._wp/gamma_L + 1._wp) + pi_inf_L / gamma_L) / rho_L)
                            a_R = sqrt((pres_R*(1._wp/gamma_R + 1._wp) + pi_inf_R / gamma_R) / rho_R)

                            cfl = (max(sqrt(vel_L(1,0)**2._wp + vel_L(2,0)**2._wp),sqrt(vel_R(1,0)**2._wp + vel_R(2,0)**2._wp)) + max(a_L,a_R))

                            !$acc loop seq
                            do i = 1, num_fluids
                                !$acc atomic
                                rhs_vf(i)%sf(j+1,k,l) = rhs_vf(i)%sf(j+1,k,l) + &
                                    (0.5_wp * (alpha_rho_L(i) * &
                                    vel_L(1,0))*(1._wp/dx(j+1)) - &
                                    0.5_wp*cfl *(alpha_rho_L(i))*(1._wp/dx(j+1)))

                                !$acc atomic
                                rhs_vf(i)%sf(j,k,l) = rhs_vf(i)%sf(j,k,l) - &
                                    (0.5_wp * (alpha_rho_L(i) * &
                                    vel_L(1,0))*(1._wp/dx(j)) - &
                                    0.5_wp*cfl * (alpha_rho_L(i))*(1._wp/dx(j))) 
                            end do

                            if(num_fluids > 1) then 
                                !$acc loop seq
                                do i = 1, num_fluids
                                    !$acc atomic
                                    rhs_vf(advxb+i-1)%sf(j+1,k,l) = rhs_vf(advxb+i-1)%sf(j+1,k,l) + &
                                        (0.5_wp * (alpha_L(i) * &
                                        vel_L(1,0))*(1._wp/dx(j+1)) - &
                                        0.5_wp*cfl*(alpha_L(i))*(1._wp/dx(j+1)))

                                    !$acc atomic
                                    rhs_vf(advxb+i-1)%sf(j+1,k,l) = rhs_vf(advxb+i-1)%sf(j+1,k,l) &
                                    - (0.5_wp * q_prim_vf(advxb+i-1)%sf(j+1,k,l) * vel_L(1,0)*(1._wp/dx(j+1)))

                                    !$acc atomic
                                    rhs_vf(advxb+i-1)%sf(j,k,l) = rhs_vf(advxb+i-1)%sf(j,k,l) - &
                                        (0.5_wp * (alpha_L(i) * &
                                        vel_L(1,0))*(1._wp/dx(j)) - &
                                        0.5_wp*cfl*(alpha_L(i))*(1._wp/dx(j)))

                                    !$acc atomic
                                    rhs_vf(advxb+i-1)%sf(j,k,l) = rhs_vf(advxb+i-1)%sf(j,k,l) &
                                    + (0.5_wp * q_prim_vf(advxb+i-1)%sf(j,k,l)  * vel_L(1,0)*(1._wp/dx(j)))
                                end do
                            end if

                            !$acc atomic
                            rhs_vf(momxb)%sf(j+1,k,l) = rhs_vf(momxb)%sf(j+1,k,l) + &
                                 (0.5_wp* (rho_L * (vel_L(1,0))**2.0 + &
                                 pres_L)*(1._wp/dx(j+1)) - &
                                 0.5_wp*cfl * (rho_L*vel_L(1,0))*(1._wp/dx(j+1)))

                            !$acc atomic
                            rhs_vf(momxb+1)%sf(j+1, k, l) =  rhs_vf(momxb+1)%sf(j+1, k, l) + &
                                (0.5_wp * rho_L * vel_L(1,0)*vel_L(2,0)*(1._wp/dx(j+1)) - &
                                0.5_wp*cfl * (rho_L*vel_L(2,0))*(1._wp/dx(j+1)))

                            !$acc atomic
                            rhs_vf(E_idx)%sf(j+1, k, l) = rhs_vf(E_idx)%sf(j+1, k, l) +  &
                                (0.5_wp * (vel_L(1,0) * (E_L + &
                                pres_L) )*(1._wp/dx(j+1)) - &
                                0.5_wp*cfl * (E_L)*(1._wp/dx(j+1)))

                            !$acc atomic
                            rhs_vf(momxb)%sf(j,k,l) = rhs_vf(momxb)%sf(j,k,l) - &
                                 (0.5_wp* (rho_L * (vel_L(1,0))**2.0 + &
                                 pres_L)*(1._wp/dx(j)) - &
                                 0.5_wp*cfl * (rho_L*vel_L(1,0))*(1._wp/dx(j)))

                            !$acc atomic
                            rhs_vf(momxb+1)%sf(j, k, l) =  rhs_vf(momxb+1)%sf(j, k, l) - &
                                (0.5_wp * rho_L * vel_L(1,0)*vel_L(2,0)*(1._wp/dx(j)) - &
                                0.5_wp*cfl * (rho_L*vel_L(2,0))*(1._wp/dx(j)))

                            !$acc atomic
                             rhs_vf(E_idx)%sf(j, k, l) = rhs_vf(E_idx)%sf(j, k, l) -  &
                                (0.5_wp * (vel_L(1,0) * (E_L + &
                                pres_L) )*(1._wp/dx(j)) - &
                                0.5_wp*cfl * (E_L)*(1._wp/dx(j)))

                            !! duy & dvx
                            dvel1 = (1/(12._wp*dy(k))) * ( &
                            8._wp*q_prim_vf(momxb)%sf(j,k+1,l)/rho_sf(1,2) - &
                            8._wp*q_prim_vf(momxb)%sf(j,k-1,l)/rho_sf(-1,2) + &
                            q_prim_vf(momxb)%sf(j,k-2,l)/rho_sf(-2,2) - &
                            q_prim_vf(momxb)%sf(j,k+2,l)/rho_sf(2,2) )

                            dvel2 = (1/(12._wp*dx(j))) * ( &
                            8._wp*q_prim_vf(momxb+1)%sf(j+1,k,l)/rho_sf(1,1) - &
                            8._wp*q_prim_vf(momxb+1)%sf(j-1,k,l)/rho_sf(-1,1) + &
                            q_prim_vf(momxb+1)%sf(j-2,k,l)/rho_sf(-2,1) - &
                            q_prim_vf(momxb+1)%sf(j+2,k,l)/rho_sf(2,1) )

                            jac_rhs(j, k, l) = alf_igr* (2._wp*dvel1*dvel2)

                            if(num_fluids > 1) then 
                                !$acc loop seq
                                do q = -3, 2 
                                    mu_L(q) = 0._wp
                                    mu_R(q) = 0._wp
                                    !$acc loop seq 
                                    do i = 1, num_fluids
                                        mu_L(q) = (1._wp/60._wp) * (-3._wp * q_prim_vf(E_idx+i)%sf(j-1+q, k, l) + &
                                                27._wp * q_prim_vf(E_idx+i)%sf(j+q, k, l) + &
                                                47._wp * q_prim_vf(E_idx+i)%sf(j+1+q, k, l) -   &
                                                13._wp * q_prim_vf(E_idx+i)%sf(j+2+q, k, l) + &
                                                2._wp * q_prim_vf(E_idx+i)%sf(j+3+q, k, l)) / Res(1, i) + mu_L(q)

                                        mu_R(q) = (1._wp/60._wp) * (-3._wp * q_prim_vf(E_idx+i)%sf(j+2+q, k, l) + &
                                                27._wp * q_prim_vf(E_idx+i)%sf(j+1+q, k, l) + &
                                                47._wp * q_prim_vf(E_idx+i)%sf(j+q, k, l) -   &
                                                13._wp * q_prim_vf(E_idx+i)%sf(j-1+q, k, l) + &
                                                2._wp * q_prim_vf(E_idx+i)%sf(j-2+q, k, l)) / Res(1, i) + mu_R(q)
                                    end do
                                end do
                            else
                                !$acc loop seq
                                do q = -3, 2
                                    mu_L(q) = 1._wp / Res(1, 1) 
                                    mu_R(q) = 1._wp / Res(1, 1) 
                                end do
                            end if

                            if(viscous) then

                                !$acc atomic
                                rhs_vf(momxb+1)%sf(j+2,k,l) = rhs_vf(momxb+1)%sf(j+2,k,l) - 0.5_wp*mu_L(1)*((-3._wp/60._wp)*(dvel1+dvel2))*(1._wp/dx(j+2))
                                !$acc atomic
                                rhs_vf(momxb+1)%sf(j+1,k,l) = rhs_vf(momxb+1)%sf(j+1,k,l) - 0.5_wp*mu_L(0)*((27._wp/60._wp)*(dvel1+dvel2))*(1._wp/dx(j+1))
                                !$acc atomic
                                rhs_vf(momxb+1)%sf(j,k,l) = rhs_vf(momxb+1)%sf(j,k,l) - 0.5_wp*mu_L(-1)*((47._wp/60._wp)*(dvel1+dvel2))*(1._wp/dx(j))
                                !$acc atomic
                                rhs_vf(momxb+1)%sf(j-1,k,l) = rhs_vf(momxb+1)%sf(j-1,k,l) - 0.5_wp*mu_L(-2)*((-13._wp/60._wp)*(dvel1+dvel2))*(1._wp/dx(j-1))
                                !$acc atomic
                                rhs_vf(momxb+1)%sf(j-2,k,l) = rhs_vf(momxb+1)%sf(j-2,k,l) - 0.5_wp*mu_L(-3)*((2._wp/60._wp)*(dvel1+dvel2))*(1._wp/dx(j-2))

                                !$acc atomic
                                rhs_vf(E_idx)%sf(j+2,k,l) = rhs_vf(E_idx)%sf(j+2,k,l) - 0.5_wp*mu_L(1)*((-3._wp/60._wp)*(dvel1+dvel2))*vel_L(2,1)*(1._wp/dx(j+2))
                                !$acc atomic
                                rhs_vf(E_idx)%sf(j+1,k,l) = rhs_vf(E_idx)%sf(j+1,k,l) - 0.5_wp*mu_L(0)*((27._wp/60._wp)*(dvel1+dvel2))*vel_L(2,0)*(1._wp/dx(j+1))
                                !$acc atomic
                                rhs_vf(E_idx)%sf(j,k,l) = rhs_vf(E_idx)%sf(j,k,l) - 0.5_wp*mu_L(-1)*((47._wp/60._wp)*(dvel1+dvel2))*vel_L(2,-1)*(1._wp/dx(j))
                                !$acc atomic
                                rhs_vf(E_idx)%sf(j-1,k,l) = rhs_vf(E_idx)%sf(j-1,k,l) - 0.5_wp*mu_L(-2)*((-13._wp/60._wp)*(dvel1+dvel2))*vel_L(2,-2)*(1._wp/dx(j-1))
                                !$acc atomic
                                rhs_vf(E_idx)%sf(j-2,k,l) = rhs_vf(E_idx)%sf(j-2,k,l) - 0.5_wp*mu_L(-3)*((2._wp/60._wp)*(dvel1+dvel2))*vel_L(2,-3)*(1._wp/dx(j-2))

                                !$acc atomic
                                rhs_vf(momxb+1)%sf(j+1,k,l) = rhs_vf(momxb+1)%sf(j+1,k,l) + 0.5_wp*mu_L(1)*((-3._wp/60._wp)*(dvel1+dvel2))*(1._wp/dx(j+1))
                                !$acc atomic
                                rhs_vf(momxb+1)%sf(j,k,l) = rhs_vf(momxb+1)%sf(j,k,l) + 0.5_wp*mu_L(0)*((27._wp/60._wp)*(dvel1+dvel2))*(1._wp/dx(j))
                                !$acc atomic
                                rhs_vf(momxb+1)%sf(j-1,k,l) = rhs_vf(momxb+1)%sf(j-1,k,l) + 0.5_wp*mu_L(-1)*((47._wp/60._wp)*(dvel1+dvel2))*(1._wp/dx(j-1))
                                !$acc atomic
                                rhs_vf(momxb+1)%sf(j-2,k,l) = rhs_vf(momxb+1)%sf(j-2,k,l) + 0.5_wp*mu_L(-2)*((-13._wp/60._wp)*(dvel1+dvel2))*(1._wp/dx(j-2))
                                !$acc atomic
                                rhs_vf(momxb+1)%sf(j-3,k,l) = rhs_vf(momxb+1)%sf(j-3,k,l) + 0.5_wp*mu_L(-3)*((2._wp/60._wp)*(dvel1+dvel2))*(1._wp/dx(j-3))

                                !$acc atomic
                                rhs_vf(E_idx)%sf(j+1,k,l) = rhs_vf(E_idx)%sf(j+1,k,l) + 0.5_wp*mu_L(1)*((-3._wp/60._wp)*(dvel1+dvel2))*vel_L(2,1)*(1._wp/dx(j+1))
                                !$acc atomic
                                rhs_vf(E_idx)%sf(j,k,l) = rhs_vf(E_idx)%sf(j,k,l) + 0.5_wp*mu_L(0)*((27._wp/60._wp)*(dvel1+dvel2))*vel_L(2,0)*(1._wp/dx(j))
                                !$acc atomic
                                rhs_vf(E_idx)%sf(j-1,k,l) = rhs_vf(E_idx)%sf(j-1,k,l) + 0.5_wp*mu_L(-1)*((47._wp/60._wp)*(dvel1+dvel2))*vel_L(2,-1)*(1._wp/dx(j-1))
                                !$acc atomic
                                rhs_vf(E_idx)%sf(j-2,k,l) = rhs_vf(E_idx)%sf(j-2,k,l) + 0.5_wp*mu_L(-2)*((-13._wp/60._wp)*(dvel1+dvel2))*vel_L(2,-2)*(1._wp/dx(j-2))
                                !$acc atomic
                                rhs_vf(E_idx)%sf(j-3,k,l) = rhs_vf(E_idx)%sf(j-3,k,l) + 0.5_wp*mu_L(-3)*((2._wp/60._wp)*(dvel1+dvel2))*vel_L(2,-3)*(1._wp/dx(j-3))
                               
                                !$acc atomic
                                rhs_vf(momxb+1)%sf(j-1,k,l) = rhs_vf(momxb+1)%sf(j-1,k,l) - 0.5_wp*mu_R(-2)*((-3._wp/60._wp)*(dvel1+dvel2))*(1._wp/dx(j-1))
                                !$acc atomic
                                rhs_vf(momxb+1)%sf(j,k,l) = rhs_vf(momxb+1)%sf(j,k,l) - 0.5_wp*mu_R(-1)*((27._wp/60._wp)*(dvel1+dvel2))*(1._wp/dx(j))
                                !$acc atomic
                                rhs_vf(momxb+1)%sf(j+1,k,l) = rhs_vf(momxb+1)%sf(j+1,k,l) - 0.5_wp*mu_R(0)*((47._wp/60._wp)*(dvel1+dvel2))*(1._wp/dx(j+1))
                                !$acc atomic
                                rhs_vf(momxb+1)%sf(j+2,k,l) = rhs_vf(momxb+1)%sf(j+2,k,l) - 0.5_wp*mu_R(1)*((-13._wp/60._wp)*(dvel1+dvel2))*(1._wp/dx(j+2))
                                !$acc atomic
                                rhs_vf(momxb+1)%sf(j+3,k,l) = rhs_vf(momxb+1)%sf(j+3,k,l) - 0.5_wp*mu_R(2)*((2._wp/60._wp)*(dvel1+dvel2))*(1._wp/dx(j+3))

                                !$acc atomic
                                rhs_vf(E_idx)%sf(j-1,k,l) = rhs_vf(E_idx)%sf(j-1,k,l) - 0.5_wp*mu_R(-2)*((-3._wp/60._wp)*(dvel1+dvel2))*vel_R(2,-2)*(1._wp/dx(j-1))
                                !$acc atomic
                                rhs_vf(E_idx)%sf(j,k,l) = rhs_vf(E_idx)%sf(j,k,l) - 0.5_wp*mu_R(-1)*((27._wp/60._wp)*(dvel1+dvel2))*vel_R(2,-1)*(1._wp/dx(j))
                                !$acc atomic
                                rhs_vf(E_idx)%sf(j+1,k,l) = rhs_vf(E_idx)%sf(j+1,k,l) - 0.5_wp*mu_R(0)*((47._wp/60._wp)*(dvel1+dvel2))*vel_R(2,0)*(1._wp/dx(j+1))
                                !$acc atomic
                                rhs_vf(E_idx)%sf(j+2,k,l) = rhs_vf(E_idx)%sf(j+2,k,l) - 0.5_wp*mu_R(1)*((-13._wp/60._wp)*(dvel1+dvel2))*vel_R(2,1)*(1._wp/dx(j+2))
                                !$acc atomic
                                rhs_vf(E_idx)%sf(j+3,k,l) = rhs_vf(E_idx)%sf(j+3,k,l) - 0.5_wp*mu_R(2)*((2._wp/60._wp)*(dvel1+dvel2))*vel_R(2,2)*(1._wp/dx(j+3)) 

                                !$acc atomic
                                rhs_vf(momxb+1)%sf(j-2,k,l) = rhs_vf(momxb+1)%sf(j-2,k,l) + 0.5_wp*mu_R(-2)*((-3._wp/60._wp)*(dvel1+dvel2))*(1._wp/dx(j-2))
                                !$acc atomic
                                rhs_vf(momxb+1)%sf(j-1,k,l) = rhs_vf(momxb+1)%sf(j-1,k,l) + 0.5_wp*mu_R(-1)*((27._wp/60._wp)*(dvel1+dvel2))*(1._wp/dx(j-1))
                                !$acc atomic
                                rhs_vf(momxb+1)%sf(j,k,l) = rhs_vf(momxb+1)%sf(j,k,l) + 0.5_wp*mu_R(0)*((47._wp/60._wp)*(dvel1+dvel2))*(1._wp/dx(j))
                                !$acc atomic
                                rhs_vf(momxb+1)%sf(j+1,k,l) = rhs_vf(momxb+1)%sf(j+1,k,l) + 0.5_wp*mu_R(1)*((-13._wp/60._wp)*(dvel1+dvel2))*(1._wp/dx(j+1))
                                !$acc atomic
                                rhs_vf(momxb+1)%sf(j+2,k,l) = rhs_vf(momxb+1)%sf(j+2,k,l) + 0.5_wp*mu_R(2)*((2._wp/60._wp)*(dvel1+dvel2))*(1._wp/dx(j+2))

                                !$acc atomic
                                rhs_vf(E_idx)%sf(j-2,k,l) = rhs_vf(E_idx)%sf(j-2,k,l) + 0.5_wp*mu_R(-2)*((-3._wp/60._wp)*(dvel1+dvel2))*vel_R(2,-2)*(1._wp/dx(j-2))
                                !$acc atomic
                                rhs_vf(E_idx)%sf(j-1,k,l) = rhs_vf(E_idx)%sf(j-1,k,l) + 0.5_wp*mu_R(-1)*((27._wp/60._wp)*(dvel1+dvel2))*vel_R(2,-1)*(1._wp/dx(j-1))
                                !$acc atomic
                                rhs_vf(E_idx)%sf(j,k,l) = rhs_vf(E_idx)%sf(j,k,l) + 0.5_wp*mu_R(0)*((47._wp/60._wp)*(dvel1+dvel2))*vel_R(2,0)*(1._wp/dx(j))
                                !$acc atomic
                                rhs_vf(E_idx)%sf(j+1,k,l) = rhs_vf(E_idx)%sf(j+1,k,l) + 0.5_wp*mu_R(1)*((-13._wp/60._wp)*(dvel1+dvel2))*vel_R(2,1)*(1._wp/dx(j+1))
                                !$acc atomic
                                rhs_vf(E_idx)%sf(j+2,k,l) = rhs_vf(E_idx)%sf(j+2,k,l) + 0.5_wp*mu_R(2)*((2._wp/60._wp)*(dvel1+dvel2))*vel_R(2,2)*(1._wp/dx(j+2))

                            end if

                            !! dux & dvy
                            dvel1 = (1/(12._wp*dx(j))) * ( &
                            8._wp*q_prim_vf(momxb)%sf(j+1,k,l)/rho_sf(1,1) - &
                            8._wp*q_prim_vf(momxb)%sf(j-1,k,l)/rho_sf(-1,1) + &
                            q_prim_vf(momxb)%sf(j-2,k,l)/rho_sf(-2,1) - &
                            q_prim_vf(momxb)%sf(j+2,k,l)/rho_sf(2,1) )


                            dvel2 = (1/(12._wp*dy(k))) * ( &
                            8._wp*q_prim_vf(momxb+1)%sf(j,k+1,l)/rho_sf(1,2) - &
                            8._wp*q_prim_vf(momxb+1)%sf(j,k-1,l)/rho_sf(-1,2) + &
                            q_prim_vf(momxb+1)%sf(j,k-2,l)/rho_sf(-2,2) - &
                            q_prim_vf(momxb+1)%sf(j,k+2,l)/rho_sf(2,2) )

                            jac_rhs(j, k, l) = jac_rhs(j,k,l) + alf_igr* (dvel1**2_wp + dvel2**2_wp + (dvel1 + dvel2)**2_wp)

                            if(viscous) then 

                                !$acc atomic
                                rhs_vf(momxb)%sf(j+2,k,l) = rhs_vf(momxb)%sf(j+2,k,l) - 0.5_wp*mu_L(1)*((-3._wp/60._wp)*(4._wp*dvel1-2._wp*dvel2)/3._wp)*(1._wp/dx(j+2))
                                !$acc atomic
                                rhs_vf(momxb)%sf(j+1,k,l) = rhs_vf(momxb)%sf(j+1,k,l) - 0.5_wp*mu_L(0)*((27._wp/60._wp)*(4._wp*dvel1-2._wp*dvel2)/3._wp)*(1._wp/dx(j+1))
                                !$acc atomic
                                rhs_vf(momxb)%sf(j,k,l) = rhs_vf(momxb)%sf(j,k,l) - 0.5_wp*mu_L(-1)*((47._wp/60._wp)*(4._wp*dvel1-2._wp*dvel2)/3._wp)*(1._wp/dx(j))
                                !$acc atomic
                                rhs_vf(momxb)%sf(j-1,k,l) = rhs_vf(momxb)%sf(j-1,k,l) - 0.5_wp*mu_L(-2)*((-13._wp/60._wp)*(4._wp*dvel1-2._wp*dvel2)/3._wp)*(1._wp/dx(j-1))
                                !$acc atomic
                                rhs_vf(momxb)%sf(j-2,k,l) = rhs_vf(momxb)%sf(j-2,k,l) - 0.5_wp*mu_L(-3)*((2._wp/60._wp)*(4._wp*dvel1-2._wp*dvel2)/3._wp)*(1._wp/dx(j-2))

                                !$acc atomic
                                rhs_vf(E_idx)%sf(j+2,k,l) = rhs_vf(E_idx)%sf(j+2,k,l) - 0.5_wp*mu_L(1)*((-3._wp/60._wp)*(4._wp*dvel1-2._wp*dvel2)/3._wp)*vel_L(1,1)*(1._wp/dx(j+2))
                                !$acc atomic
                                rhs_vf(E_idx)%sf(j+1,k,l) = rhs_vf(E_idx)%sf(j+1,k,l) - 0.5_wp*mu_L(0)*((27._wp/60._wp)*(4._wp*dvel1-2._wp*dvel2)/3._wp)*vel_L(1,0)*(1._wp/dx(j+1))
                                !$acc atomic
                                rhs_vf(E_idx)%sf(j,k,l) = rhs_vf(E_idx)%sf(j,k,l) - 0.5_wp*mu_L(-1)*((47._wp/60._wp)*(4._wp*dvel1-2._wp*dvel2)/3._wp)*vel_L(1,-1)*(1._wp/dx(j))
                                !$acc atomic
                                rhs_vf(E_idx)%sf(j-1,k,l) = rhs_vf(E_idx)%sf(j-1,k,l) - 0.5_wp*mu_L(-2)*((-13._wp/60._wp)*(4._wp*dvel1-2._wp*dvel2)/3._wp)*vel_L(1,-2)*(1._wp/dx(j-1))
                                !$acc atomic
                                rhs_vf(E_idx)%sf(j-2,k,l) = rhs_vf(E_idx)%sf(j-2,k,l) - 0.5_wp*mu_L(-3)*((2._wp/60._wp)*(4._wp*dvel1-2._wp*dvel2)/3._wp)*vel_L(1,-3)*(1._wp/dx(j-2))

                                !$acc atomic
                                rhs_vf(momxb)%sf(j+1,k,l) = rhs_vf(momxb)%sf(j+1,k,l) + 0.5_wp*mu_L(1)*((-3._wp/60._wp)*(4._wp*dvel1-2._wp*dvel2)/3._wp)*(1._wp/dx(j+1))
                                !$acc atomic
                                rhs_vf(momxb)%sf(j,k,l) = rhs_vf(momxb)%sf(j,k,l) + 0.5_wp*mu_L(0)*((27._wp/60._wp)*(4._wp*dvel1-2._wp*dvel2)/3._wp)*(1._wp/dx(j))
                                !$acc atomic
                                rhs_vf(momxb)%sf(j-1,k,l) = rhs_vf(momxb)%sf(j-1,k,l) + 0.5_wp*mu_L(-1)*((47._wp/60._wp)*(4._wp*dvel1-2._wp*dvel2)/3._wp)*(1._wp/dx(j-1))
                                !$acc atomic
                                rhs_vf(momxb)%sf(j-2,k,l) = rhs_vf(momxb)%sf(j-2,k,l) + 0.5_wp*mu_L(-2)*((-13._wp/60._wp)*(4._wp*dvel1-2._wp*dvel2)/3._wp)*(1._wp/dx(j-2))
                                !$acc atomic
                                rhs_vf(momxb)%sf(j-3,k,l) = rhs_vf(momxb)%sf(j-3,k,l) + 0.5_wp*mu_L(-3)*((2._wp/60._wp)*(4._wp*dvel1-2._wp*dvel2)/3._wp)*(1._wp/dx(j-3))

                                !$acc atomic
                                rhs_vf(E_idx)%sf(j+1,k,l) = rhs_vf(E_idx)%sf(j+1,k,l) + 0.5_wp*mu_L(1)*((-3._wp/60._wp)*(4._wp*dvel1-2._wp*dvel2)/3._wp)*vel_L(1,1)*(1._wp/dx(j+1))
                                !$acc atomic
                                rhs_vf(E_idx)%sf(j,k,l) = rhs_vf(E_idx)%sf(j,k,l) + 0.5_wp*mu_L(0)*((27._wp/60._wp)*(4._wp*dvel1-2._wp*dvel2)/3._wp)*vel_L(1,0)*(1._wp/dx(j))
                                !$acc atomic
                                rhs_vf(E_idx)%sf(j-1,k,l) = rhs_vf(E_idx)%sf(j-1,k,l) + 0.5_wp*mu_L(-1)*((47._wp/60._wp)*(4._wp*dvel1-2._wp*dvel2)/3._wp)*vel_L(1,-1)*(1._wp/dx(j-1))
                                !$acc atomic
                                rhs_vf(E_idx)%sf(j-2,k,l) = rhs_vf(E_idx)%sf(j-2,k,l) + 0.5_wp*mu_L(-2)*((-13._wp/60._wp)*(4._wp*dvel1-2._wp*dvel2)/3._wp)*vel_L(1,-2)*(1._wp/dx(j-2))
                                !$acc atomic
                                rhs_vf(E_idx)%sf(j-3,k,l) = rhs_vf(E_idx)%sf(j-3,k,l) + 0.5_wp*mu_L(-3)*((2._wp/60._wp)*(4._wp*dvel1-2._wp*dvel2)/3._wp)*vel_L(1,-3)*(1._wp/dx(j-3))

                               
                                !$acc atomic
                                rhs_vf(momxb)%sf(j-1,k,l) = rhs_vf(momxb)%sf(j-1,k,l) - 0.5_wp*mu_R(-2)*((-3._wp/60._wp)*(4._wp*dvel1-2._wp*dvel2)/3._wp)*(1._wp/dx(j-1))
                                !$acc atomic
                                rhs_vf(momxb)%sf(j,k,l) = rhs_vf(momxb)%sf(j,k,l) - 0.5_wp*mu_R(-1)*((27._wp/60._wp)*(4._wp*dvel1-2._wp*dvel2)/3._wp)*(1._wp/dx(j))
                                !$acc atomic
                                rhs_vf(momxb)%sf(j+1,k,l) = rhs_vf(momxb)%sf(j+1,k,l) - 0.5_wp*mu_R(0)*((47._wp/60._wp)*(4._wp*dvel1-2._wp*dvel2)/3._wp)*(1._wp/dx(j+1))
                                !$acc atomic
                                rhs_vf(momxb)%sf(j+2,k,l) = rhs_vf(momxb)%sf(j+2,k,l) - 0.5_wp*mu_R(1)*((-13._wp/60._wp)*(4._wp*dvel1-2._wp*dvel2)/3._wp)*(1._wp/dx(j+2))
                                !$acc atomic
                                rhs_vf(momxb)%sf(j+3,k,l) = rhs_vf(momxb)%sf(j+3,k,l) - 0.5_wp*mu_R(2)*((2._wp/60._wp)*(4._wp*dvel1-2._wp*dvel2)/3._wp)*(1._wp/dx(j+3))

                                !$acc atomic
                                rhs_vf(E_idx)%sf(j-1,k,l) = rhs_vf(E_idx)%sf(j-1,k,l) - 0.5_wp*mu_R(-2)*((-3._wp/60._wp)*(4._wp*dvel1-2._wp*dvel2)/3._wp)*vel_R(1,-2)*(1._wp/dx(j-1))
                                !$acc atomic
                                rhs_vf(E_idx)%sf(j,k,l) = rhs_vf(E_idx)%sf(j,k,l) - 0.5_wp*mu_R(-1)*((27._wp/60._wp)*(4._wp*dvel1-2._wp*dvel2)/3._wp)*vel_R(1,-1)*(1._wp/dx(j))
                                !$acc atomic
                                rhs_vf(E_idx)%sf(j+1,k,l) = rhs_vf(E_idx)%sf(j+1,k,l) - 0.5_wp*mu_R(0)*((47._wp/60._wp)*(4._wp*dvel1-2._wp*dvel2)/3._wp)*vel_R(1,0)*(1._wp/dx(j+1))
                                !$acc atomic
                                rhs_vf(E_idx)%sf(j+2,k,l) = rhs_vf(E_idx)%sf(j+2,k,l) - 0.5_wp*mu_R(1)*((-13._wp/60._wp)*(4._wp*dvel1-2._wp*dvel2)/3._wp)*vel_R(1,1)*(1._wp/dx(j+2))
                                !$acc atomic
                                rhs_vf(E_idx)%sf(j+3,k,l) = rhs_vf(E_idx)%sf(j+3,k,l) - 0.5_wp*mu_R(2)*((2._wp/60._wp)*(4._wp*dvel1-2._wp*dvel2)/3._wp)*vel_R(1,2)*(1._wp/dx(j+3)) 

                                !$acc atomic
                                rhs_vf(momxb)%sf(j-2,k,l) = rhs_vf(momxb)%sf(j-2,k,l) + 0.5_wp*mu_R(-2)*((-3._wp/60._wp)*(4._wp*dvel1-2._wp*dvel2)/3._wp)*(1._wp/dx(j-2))
                                !$acc atomic
                                rhs_vf(momxb)%sf(j-1,k,l) = rhs_vf(momxb)%sf(j-1,k,l) + 0.5_wp*mu_R(-1)*((27._wp/60._wp)*(4._wp*dvel1-2._wp*dvel2)/3._wp)*(1._wp/dx(j-1))
                                !$acc atomic
                                rhs_vf(momxb)%sf(j,k,l) = rhs_vf(momxb)%sf(j,k,l) + 0.5_wp*mu_R(0)*((47._wp/60._wp)*(4._wp*dvel1-2._wp*dvel2)/3._wp)*(1._wp/dx(j))
                                !$acc atomic
                                rhs_vf(momxb)%sf(j+1,k,l) = rhs_vf(momxb)%sf(j+1,k,l) + 0.5_wp*mu_R(1)*((-13._wp/60._wp)*(4._wp*dvel1-2._wp*dvel2)/3._wp)*(1._wp/dx(j+1))
                                !$acc atomic
                                rhs_vf(momxb)%sf(j+2,k,l) = rhs_vf(momxb)%sf(j+2,k,l) + 0.5_wp*mu_R(2)*((2._wp/60._wp)*(4._wp*dvel1-2._wp*dvel2)/3._wp)*(1._wp/dx(j+2))

                                !$acc atomic
                                rhs_vf(E_idx)%sf(j-2,k,l) = rhs_vf(E_idx)%sf(j-2,k,l) + 0.5_wp*mu_R(-2)*((-3._wp/60._wp)*(4._wp*dvel1-2._wp*dvel2)/3._wp)*vel_R(1,-2)*(1._wp/dx(j-2))
                                !$acc atomic
                                rhs_vf(E_idx)%sf(j-1,k,l) = rhs_vf(E_idx)%sf(j-1,k,l) + 0.5_wp*mu_R(-1)*((27._wp/60._wp)*(4._wp*dvel1-2._wp*dvel2)/3._wp)*vel_R(1,-1)*(1._wp/dx(j-1))
                                !$acc atomic
                                rhs_vf(E_idx)%sf(j,k,l) = rhs_vf(E_idx)%sf(j,k,l) + 0.5_wp*mu_R(0)*((47._wp/60._wp)*(4._wp*dvel1-2._wp*dvel2)/3._wp)*vel_R(1,0)*(1._wp/dx(j))
                                !$acc atomic
                                rhs_vf(E_idx)%sf(j+1,k,l) = rhs_vf(E_idx)%sf(j+1,k,l) + 0.5_wp*mu_R(1)*((-13._wp/60._wp)*(4._wp*dvel1-2._wp*dvel2)/3._wp)*vel_R(1,1)*(1._wp/dx(j+1))
                                !$acc atomic
                                rhs_vf(E_idx)%sf(j+2,k,l) = rhs_vf(E_idx)%sf(j+2,k,l) + 0.5_wp*mu_R(2)*((2._wp/60._wp)*(4._wp*dvel1-2._wp*dvel2)/3._wp)*vel_R(1,2)*(1._wp/dx(j+2))

                            end if

                            !$acc loop seq
                            do i = 1, num_fluids
                                !$acc atomic
                                rhs_vf(i)%sf(j+1,k,l) = rhs_vf(i)%sf(j+1,k,l) + &
                                    (0.5_wp * (alpha_rho_R(i) * &
                                    vel_R(1,0))*(1._wp/dx(j+1)) + &
                                    0.5_wp*cfl * (alpha_rho_R(i))*(1._wp/dx(j+1)))

                                !$acc atomic
                                rhs_vf(i)%sf(j,k,l) = rhs_vf(i)%sf(j,k,l) - &
                                    (0.5_wp * (alpha_rho_R(i) * &
                                    vel_R(1,0))*(1._wp/dx(j)) + &
                                    0.5_wp*cfl * (alpha_rho_R(i))*(1._wp/dx(j)))
                            end do

                            if(num_fluids > 1) then 
                                !$acc loop seq
                                do i = 1, num_fluids
                                    !$acc atomic
                                    rhs_vf(advxb+i-1)%sf(j+1,k,l) = rhs_vf(advxb+i-1)%sf(j+1,k,l) + &
                                        (0.5_wp * (alpha_R(i) * &
                                        vel_R(1,0))*(1._wp/dx(j+1)) + &
                                        0.5_wp*cfl*(alpha_R(i))*(1._wp/dx(j+1)))

                                    !$acc atomic
                                    rhs_vf(advxb+i-1)%sf(j+1,k,l) = rhs_vf(advxb+i-1)%sf(j+1,k,l) &
                                    - (0.5_wp * q_prim_vf(advxb+i-1)%sf(j+1,k,l)  * vel_R(1,0)*(1._wp/dx(j+1)))

                                    !$acc atomic
                                    rhs_vf(advxb+i-1)%sf(j,k,l) = rhs_vf(advxb+i-1)%sf(j,k,l) - &
                                        (0.5_wp * (alpha_R(i) * &
                                        vel_R(1,0))*(1._wp/dx(j)) + &
                                        0.5_wp*cfl*(alpha_R(i))*(1._wp/dx(j)))

                                    !$acc atomic
                                    rhs_vf(advxb+i-1)%sf(j,k,l) = rhs_vf(advxb+i-1)%sf(j,k,l) &
                                    + (0.5_wp * q_prim_vf(advxb+i-1)%sf(j,k,l)  * vel_R(1,0)*(1._wp/dx(j)))
                                end do
                            end if

                            !$acc atomic
                            rhs_vf(momxb)%sf(j+1,k,l) = rhs_vf(momxb)%sf(j+1,k,l) + &
                                 (0.5_wp* (rho_R * (vel_R(1,0))**2.0 + &
                                 pres_R)*(1._wp/dx(j+1)) + &
                                 0.5_wp*cfl * (rho_R*vel_R(1,0))*(1._wp/dx(j+1)))

                            !$acc atomic
                            rhs_vf(momxb+1)%sf(j+1, k, l) =  rhs_vf(momxb+1)%sf(j+1, k, l) + &
                                (0.5_wp * rho_R * vel_R(1,0)*vel_R(2,0)*(1._wp/dx(j+1)) + &
                                0.5_wp*cfl * (rho_R*vel_R(2,0))*(1._wp/dx(j+1)))

                            !$acc atomic
                             rhs_vf(E_idx)%sf(j+1, k, l) = rhs_vf(E_idx)%sf(j+1, k, l) +  &
                                (0.5_wp * (vel_R(1,0) * (E_R + &
                                pres_R) )*(1._wp/dx(j+1)) + &
                                0.5_wp*cfl * (E_R)*(1._wp/dx(j+1)))

                            !$acc atomic
                            rhs_vf(momxb)%sf(j,k,l) = rhs_vf(momxb)%sf(j,k,l) - &
                                 (0.5_wp* (rho_R * (vel_R(1,0))**2.0 + &
                                 pres_R)*(1._wp/dx(j)) + &
                                 0.5_wp*cfl * (rho_R*vel_R(1,0))*(1._wp/dx(j)))

                            !$acc atomic
                            rhs_vf(momxb+1)%sf(j, k, l) =  rhs_vf(momxb+1)%sf(j, k, l) - &
                                (0.5_wp * rho_R * vel_R(1,0)*vel_R(2,0)*(1._wp/dx(j)) + &
                                0.5_wp*cfl * (rho_R*vel_R(2,0))*(1._wp/dx(j)))

                            !$acc atomic
                             rhs_vf(E_idx)%sf(j, k, l) = rhs_vf(E_idx)%sf(j, k, l) -  &
                                (0.5_wp * (vel_R(1,0) * (E_R + &
                                pres_R) )*(1._wp/dx(j)) + &
                                0.5_wp*cfl * (E_R)*(1._wp/dx(j)))
                        end do
                    end do
                end do
            else
                !$acc parallel loop collapse(3) gang vector default(present) private(rho_L,gamma_L,pi_inf_L,mu_L,a_L,vel_L,vel_R, pres_L,alpha_L,alpha_R,alpha_rho_L,cfl,dvel1,dvel2,F_L,E_L,mu_R,rho_sf,alpha_rho_R)
                do l = 0, p
                    do k = 0, n
                        do j = -buff_size+5, m+buff_size-5

                            !$acc loop seq 
                            do i = 1, num_fluids
                                alpha_rho_L(i) = (1._wp/60._wp) * (-3._wp * q_prim_vf(i)%sf(j-1, k, l) + &
                                                27._wp * q_prim_vf(i)%sf(j, k, l) + &
                                                47._wp * q_prim_vf(i)%sf(j+1, k, l) -   &
                                                13._wp * q_prim_vf(i)%sf(j+2, k, l) + &
                                                2._wp * q_prim_vf(i)%sf(j+3, k, l))

                                alpha_rho_R(i) = (1._wp/60._wp) * (-3._wp * q_prim_vf(i)%sf(j+2, k, l) + &
                                                27._wp * q_prim_vf(i)%sf(j+1, k, l) + &
                                                47._wp * q_prim_vf(i)%sf(j, k, l) -   &
                                                13._wp * q_prim_vf(i)%sf(j-1, k, l) + &
                                                2._wp * q_prim_vf(i)%sf(j-2, k, l))

                                alpha_L(i) =  (1._wp/60._wp) * (-3._wp * q_prim_vf(E_idx+i)%sf(j-1, k, l) + &
                                                27._wp * q_prim_vf(E_idx+i)%sf(j, k, l) + &
                                                47._wp * q_prim_vf(E_idx+i)%sf(j+1, k, l) -   &
                                                13._wp * q_prim_vf(E_idx+i)%sf(j+2, k, l) + &
                                                2._wp * q_prim_vf(E_idx+i)%sf(j+3, k, l))
                            end do 
                            
                            if(num_fluids > 1) then  
                                !$acc loop seq 
                                do i = 1,num_fluids                           
                                    alpha_R(i) = (1._wp/60._wp) * (-3._wp * q_prim_vf(E_idx+i)%sf(j+2, k, l) + &
                                                        27._wp * q_prim_vf(E_idx+i)%sf(j+1, k, l) + &
                                                        47._wp * q_prim_vf(E_idx+i)%sf(j, k, l) -   &
                                                        13._wp * q_prim_vf(E_idx+i)%sf(j-1, k, l) + &
                                                        2._wp * q_prim_vf(E_idx+i)%sf(j-2, k, l))
                                end do 
                            else 
                                alpha_R(1) = 1._wp
                            end if

                            !$acc loop seq 
                            do i = -2, 2
                                rho_L = 0._wp
                                !$acc loop seq 
                                do q = 1, num_fluids
                                    rho_L = rho_L + q_prim_vf(q)%sf(j+i,k,l)
                                end do
                                rho_sf(i,1) = rho_L
                            end do

                            !$acc loop seq 
                            do i = -2, 2
                                rho_L = 0._wp
                                !$acc loop seq 
                                do q = 1, num_fluids
                                    rho_L = rho_L + q_prim_vf(q)%sf(j,k+i,l)
                                end do
                                rho_sf(i,2) = rho_L
                            end do

                            !$acc loop seq 
                            do i = -2, 2
                                rho_L = 0._wp
                                !$acc loop seq 
                                do q = 1, num_fluids
                                    rho_L = rho_L + q_prim_vf(q)%sf(j,k,l+i)
                                end do
                                rho_sf(i,3) = rho_L
                            end do

                            rho_L = sum(alpha_rho_L)
                            gamma_L = sum(alpha_L*gammas)
                            pi_inf_L = sum(alpha_L*pi_infs)

                            rho_R = sum(alpha_rho_R)
                            gamma_R = sum(alpha_R*gammas)
                            pi_inf_R = sum(alpha_R*pi_infs)

                            !$acc loop seq 
                            do q = -3, 2
                                !$acc loop seq 
                                do i = 1, num_dims
                                    vel_L(i,q) = (1._wp/60._wp) * (-3._wp * q_prim_vf(momxb+i-1)%sf(j-1+q,k,l) + &
                                                    27._wp * q_prim_vf(momxb+i-1)%sf(j+q,k,l) + &
                                                    47._wp * q_prim_vf(momxb+i-1)%sf(j+1+q,k,l) -   &
                                                    13._wp * q_prim_vf(momxb+i-1)%sf(j+2+q,k,l) + &
                                                    2._wp * q_prim_vf(momxb+i-1)%sf(j+3+q,k,l)) / rho_L
                                end do 

                                !$acc loop seq 
                                do i = 1, num_dims
                                    vel_R(i,q) = (1._wp/60._wp) * (-3._wp * q_prim_vf(momxb+i-1)%sf(j+2+q,k,l) + &
                                                    27._wp * q_prim_vf(momxb+i-1)%sf(j+1+q,k,l) + &
                                                    47._wp * q_prim_vf(momxb+i-1)%sf(j+q,k,l) -   &
                                                    13._wp * q_prim_vf(momxb+i-1)%sf(j-1+q,k,l) + &
                                                    2._wp * q_prim_vf(momxb+i-1)%sf(j-2+q,k,l)) / rho_R
                                end do
                            end do

                            E_L = (1._wp/60._wp) * (-3._wp * q_prim_vf(E_idx)%sf(j-1,k,l) + &
                                                27._wp * q_prim_vf(E_idx)%sf(j,k,l) + &
                                                47._wp * q_prim_vf(E_idx)%sf(j+1,k,l) -   &
                                                13._wp * q_prim_vf(E_idx)%sf(j+2,k,l) + &
                                                2._wp * q_prim_vf(E_idx)%sf(j+3,k,l))

                            E_R = (1._wp/60._wp) * (-3._wp * q_prim_vf(E_idx)%sf(j+2,k,l) + &
                                                27._wp * q_prim_vf(E_idx)%sf(j+1,k,l) + &
                                                47._wp * q_prim_vf(E_idx)%sf(j,k,l) -   &
                                                13._wp * q_prim_vf(E_idx)%sf(j-1,k,l) + &
                                                2._wp * q_prim_vf(E_idx)%sf(j-2,k,l))       

                            pres_L = (E_L - pi_inf_L - 0.5_wp*rho_L*(vel_L(1,0)**2._wp + vel_L(2,0)**2._wp + vel_L(3,0)**2._wp))/gamma_L                    

                            pres_R = (E_R - pi_inf_R - 0.5_wp*rho_R*(vel_R(1,0)**2._wp + vel_R(2,0)**2._wp + vel_R(3,0)**2._wp))/gamma_R 

                            a_L = sqrt((pres_L*(1._wp/gamma_L + 1._wp) + pi_inf_L / gamma_L) / rho_L)
                            a_R = sqrt((pres_R*(1._wp/gamma_R + 1._wp) + pi_inf_R / gamma_R) / rho_R)

                            cfl = (max(sqrt(vel_L(1,0)**2._wp + vel_L(2,0)**2._wp + vel_L(3,0)**2._wp),sqrt(vel_R(1,0)**2._wp + vel_R(2,0)**2._wp + vel_R(3,0)**2._wp)) + max(a_L,a_R))

                            !$acc loop seq
                            do i = 1, num_fluids
                                !$acc atomic
                                rhs_vf(i)%sf(j+1,k,l) = rhs_vf(i)%sf(j+1,k,l) + &
                                    (0.5_wp * (alpha_rho_L(i) * &
                                    vel_L(1,0))*(1._wp/dx(j+1)) - &
                                    0.5_wp*cfl * (alpha_rho_L(i))*(1._wp/dx(j+1)))

                                !$acc atomic
                                rhs_vf(i)%sf(j,k,l) = rhs_vf(i)%sf(j,k,l) - &
                                    (0.5_wp * (alpha_rho_L(i) * &
                                    vel_L(1,0))*(1._wp/dx(j)) - &
                                    0.5_wp*cfl * (alpha_rho_L(i))*(1._wp/dx(j)))
                            end do

                            if(num_fluids > 1) then 
                                !$acc loop seq
                                do i = 1, num_fluids
                                    !$acc atomic
                                    rhs_vf(advxb+i-1)%sf(j+1,k,l) = rhs_vf(advxb+i-1)%sf(j+1,k,l) + &
                                        (0.5_wp * (alpha_L(i) * &
                                        vel_L(1,0))*(1._wp/dx(j+1)) - &
                                        0.5_wp*cfl*(alpha_L(i))*(1._wp/dx(j+1)))

                                    !$acc atomic
                                    rhs_vf(advxb+i-1)%sf(j+1,k,l) = rhs_vf(advxb+i-1)%sf(j+1,k,l) &
                                    - (0.5_wp * q_prim_vf(advxb+i-1)%sf(j+1,k,l)  * vel_L(1,0)*(1._wp/dx(j+1)))

                                    !$acc atomic
                                    rhs_vf(advxb+i-1)%sf(j,k,l) = rhs_vf(advxb+i-1)%sf(j,k,l) - &
                                        (0.5_wp * (alpha_L(i) * &
                                        vel_L(1,0))*(1._wp/dx(j)) - &
                                        0.5_wp*cfl*(alpha_L(i))*(1._wp/dx(j)))

                                    !$acc atomic
                                    rhs_vf(advxb+i-1)%sf(j,k,l) = rhs_vf(advxb+i-1)%sf(j,k,l) &
                                    + (0.5_wp * q_prim_vf(advxb+i-1)%sf(j,k,l)  * vel_L(1,0)*(1._wp/dx(j)))
                                end do
                            end if

                            !$acc atomic
                            rhs_vf(momxb)%sf(j+1,k,l) = rhs_vf(momxb)%sf(j+1,k,l) + &
                                 (0.5_wp* (rho_L * (vel_L(1,0))**2.0 + &
                                 pres_L)*(1._wp/dx(j+1)) - &
                                 0.5_wp*cfl * (rho_L*vel_L(1,0))*(1._wp/dx(j+1)))

                            !$acc atomic
                            rhs_vf(momxb+1)%sf(j+1, k, l) =  rhs_vf(momxb+1)%sf(j+1, k, l) + &
                                (0.5_wp * rho_L * vel_L(1,0)*vel_L(2,0)*(1._wp/dx(j+1)) - &
                                0.5_wp*cfl * (rho_L*vel_L(2,0))*(1._wp/dx(j+1)))

                            !$acc atomic
                            rhs_vf(momxb+2)%sf(j+1, k, l) =  rhs_vf(momxb+2)%sf(j+1, k, l) + &
                                (0.5_wp * rho_L * vel_L(1,0)*vel_L(3,0)*(1._wp/dx(j+1)) - &
                                0.5_wp*cfl * (rho_L*vel_L(3,0))*(1._wp/dx(j+1)))

                            !$acc atomic
                             rhs_vf(E_idx)%sf(j+1, k, l) = rhs_vf(E_idx)%sf(j+1, k, l) +  &
                                (0.5_wp * (vel_L(1,0) * (E_L + &
                                pres_L) )*(1._wp/dx(j+1)) - &
                                0.5_wp*cfl * (E_L)*(1._wp/dx(j+1)))

                            !$acc atomic
                            rhs_vf(momxb)%sf(j,k,l) = rhs_vf(momxb)%sf(j,k,l) - &
                                 (0.5_wp* (rho_L * (vel_L(1,0))**2.0 + &
                                 pres_L)*(1._wp/dx(j)) - &
                                 0.5_wp*cfl * (rho_L*vel_L(1,0))*(1._wp/dx(j)))

                            !$acc atomic
                            rhs_vf(momxb+1)%sf(j, k, l) =  rhs_vf(momxb+1)%sf(j, k, l) - &
                                (0.5_wp * rho_L * vel_L(1,0)*vel_L(2,0)*(1._wp/dx(j)) - &
                                0.5_wp*cfl * (rho_L*vel_L(2,0))*(1._wp/dx(j)))

                            !$acc atomic
                            rhs_vf(momxb+2)%sf(j, k, l) =  rhs_vf(momxb+2)%sf(j, k, l) - &
                                (0.5_wp * rho_L * vel_L(1,0)*vel_L(3,0)*(1._wp/dx(j)) - &
                                0.5_wp*cfl * (rho_L*vel_L(3,0))*(1._wp/dx(j)))

                            !$acc atomic
                             rhs_vf(E_idx)%sf(j, k, l) = rhs_vf(E_idx)%sf(j, k, l) -  &
                                (0.5_wp * (vel_L(1,0) * (E_L + &
                                pres_L) )*(1._wp/dx(j)) - &
                                0.5_wp*cfl * (E_L)*(1._wp/dx(j)))

                            !! duy & dvx
                            dvel1 = (1/(12._wp*dy(k))) * ( &
                            8._wp*q_prim_vf(momxb)%sf(j,k+1,l)/rho_sf(1,2) - &
                            8._wp*q_prim_vf(momxb)%sf(j,k-1,l)/rho_sf(-1,2) + &
                            q_prim_vf(momxb)%sf(j,k-2,l)/rho_sf(-2,2) - &
                            q_prim_vf(momxb)%sf(j,k+2,l)/rho_sf(2,2) )

                            dvel2 = (1/(12._wp*dx(j))) * ( &
                            8._wp*q_prim_vf(momxb+1)%sf(j+1,k,l)/rho_sf(1,1) - &
                            8._wp*q_prim_vf(momxb+1)%sf(j-1,k,l)/rho_sf(-1,1) + &
                            q_prim_vf(momxb+1)%sf(j-2,k,l)/rho_sf(-2,1) - &
                            q_prim_vf(momxb+1)%sf(j+2,k,l)/rho_sf(2,1) )

                            jac_rhs(j, k, l) = alf_igr* (2._wp*dvel1*dvel2)

                            if(num_fluids > 1) then 
                                !$acc loop seq
                                do q = -3, 2 
                                    mu_L(q) = 0._wp
                                    mu_R(q) = 0._wp
                                    !$acc loop seq 
                                    do i = 1, num_fluids
                                        mu_L(q) = (1._wp/60._wp) * (-3._wp * q_prim_vf(E_idx+i)%sf(j-1+q, k, l) + &
                                                27._wp * q_prim_vf(E_idx+i)%sf(j+q, k, l) + &
                                                47._wp * q_prim_vf(E_idx+i)%sf(j+1+q, k, l) -   &
                                                13._wp * q_prim_vf(E_idx+i)%sf(j+2+q, k, l) + &
                                                2._wp * q_prim_vf(E_idx+i)%sf(j+3+q, k, l)) / Res(1, i) + mu_L(q)
                                        
                                        mu_R(q) = (1._wp/60._wp) * (-3._wp * q_prim_vf(E_idx+i)%sf(j+2+q, k, l) + &
                                                27._wp * q_prim_vf(E_idx+i)%sf(j+1+q, k, l) + &
                                                47._wp * q_prim_vf(E_idx+i)%sf(j+q, k, l) -   &
                                                13._wp * q_prim_vf(E_idx+i)%sf(j-1+q, k, l) + &
                                                2._wp * q_prim_vf(E_idx+i)%sf(j-2+q, k, l)) / Res(1, i) + mu_R(q)
                                    end do
                                end do
                            else
                                !$acc loop seq
                                do q = -3, 2
                                    mu_L(q) = 1._wp / Res(1, 1) 
                                    mu_R(q) = 1._wp / Res(1, 1) 
                                end do
                            end if

                            if(viscous) then

                                !$acc atomic
                                rhs_vf(momxb+1)%sf(j+2,k,l) = rhs_vf(momxb+1)%sf(j+2,k,l) - 0.5_wp*mu_L(1)*((-3._wp/60._wp)*(dvel1+dvel2))*(1._wp/dx(j+2))
                                !$acc atomic
                                rhs_vf(momxb+1)%sf(j+1,k,l) = rhs_vf(momxb+1)%sf(j+1,k,l) - 0.5_wp*mu_L(0)*((27._wp/60._wp)*(dvel1+dvel2))*(1._wp/dx(j+1))
                                !$acc atomic
                                rhs_vf(momxb+1)%sf(j,k,l) = rhs_vf(momxb+1)%sf(j,k,l) - 0.5_wp*mu_L(-1)*((47._wp/60._wp)*(dvel1+dvel2))*(1._wp/dx(j))
                                !$acc atomic
                                rhs_vf(momxb+1)%sf(j-1,k,l) = rhs_vf(momxb+1)%sf(j-1,k,l) - 0.5_wp*mu_L(-2)*((-13._wp/60._wp)*(dvel1+dvel2))*(1._wp/dx(j-1))
                                !$acc atomic
                                rhs_vf(momxb+1)%sf(j-2,k,l) = rhs_vf(momxb+1)%sf(j-2,k,l) - 0.5_wp*mu_L(-3)*((2._wp/60._wp)*(dvel1+dvel2))*(1._wp/dx(j-2))

                                !$acc atomic
                                rhs_vf(E_idx)%sf(j+2,k,l) = rhs_vf(E_idx)%sf(j+2,k,l) - 0.5_wp*mu_L(1)*((-3._wp/60._wp)*(dvel1+dvel2))*vel_L(2,1)*(1._wp/dx(j+2))
                                !$acc atomic
                                rhs_vf(E_idx)%sf(j+1,k,l) = rhs_vf(E_idx)%sf(j+1,k,l) - 0.5_wp*mu_L(0)*((27._wp/60._wp)*(dvel1+dvel2))*vel_L(2,0)*(1._wp/dx(j+1))
                                !$acc atomic
                                rhs_vf(E_idx)%sf(j,k,l) = rhs_vf(E_idx)%sf(j,k,l) - 0.5_wp*mu_L(-1)*((47._wp/60._wp)*(dvel1+dvel2))*vel_L(2,-1)*(1._wp/dx(j))
                                !$acc atomic
                                rhs_vf(E_idx)%sf(j-1,k,l) = rhs_vf(E_idx)%sf(j-1,k,l) - 0.5_wp*mu_L(-2)*((-13._wp/60._wp)*(dvel1+dvel2))*vel_L(2,-2)*(1._wp/dx(j-1))
                                !$acc atomic
                                rhs_vf(E_idx)%sf(j-2,k,l) = rhs_vf(E_idx)%sf(j-2,k,l) - 0.5_wp*mu_L(-3)*((2._wp/60._wp)*(dvel1+dvel2))*vel_L(2,-3)*(1._wp/dx(j-2))

                                !$acc atomic
                                rhs_vf(momxb+1)%sf(j+1,k,l) = rhs_vf(momxb+1)%sf(j+1,k,l) + 0.5_wp*mu_L(1)*((-3._wp/60._wp)*(dvel1+dvel2))*(1._wp/dx(j+1))
                                !$acc atomic
                                rhs_vf(momxb+1)%sf(j,k,l) = rhs_vf(momxb+1)%sf(j,k,l) + 0.5_wp*mu_L(0)*((27._wp/60._wp)*(dvel1+dvel2))*(1._wp/dx(j))
                                !$acc atomic
                                rhs_vf(momxb+1)%sf(j-1,k,l) = rhs_vf(momxb+1)%sf(j-1,k,l) + 0.5_wp*mu_L(-1)*((47._wp/60._wp)*(dvel1+dvel2))*(1._wp/dx(j-1))
                                !$acc atomic
                                rhs_vf(momxb+1)%sf(j-2,k,l) = rhs_vf(momxb+1)%sf(j-2,k,l) + 0.5_wp*mu_L(-2)*((-13._wp/60._wp)*(dvel1+dvel2))*(1._wp/dx(j-2))
                                !$acc atomic
                                rhs_vf(momxb+1)%sf(j-3,k,l) = rhs_vf(momxb+1)%sf(j-3,k,l) + 0.5_wp*mu_L(-3)*((2._wp/60._wp)*(dvel1+dvel2))*(1._wp/dx(j-3))

                                !$acc atomic
                                rhs_vf(E_idx)%sf(j+1,k,l) = rhs_vf(E_idx)%sf(j+1,k,l) + 0.5_wp*mu_L(1)*((-3._wp/60._wp)*(dvel1+dvel2))*vel_L(2,1)*(1._wp/dx(j+1))
                                !$acc atomic
                                rhs_vf(E_idx)%sf(j,k,l) = rhs_vf(E_idx)%sf(j,k,l) + 0.5_wp*mu_L(0)*((27._wp/60._wp)*(dvel1+dvel2))*vel_L(2,0)*(1._wp/dx(j))
                                !$acc atomic
                                rhs_vf(E_idx)%sf(j-1,k,l) = rhs_vf(E_idx)%sf(j-1,k,l) + 0.5_wp*mu_L(-1)*((47._wp/60._wp)*(dvel1+dvel2))*vel_L(2,-1)*(1._wp/dx(j-1))
                                !$acc atomic
                                rhs_vf(E_idx)%sf(j-2,k,l) = rhs_vf(E_idx)%sf(j-2,k,l) + 0.5_wp*mu_L(-2)*((-13._wp/60._wp)*(dvel1+dvel2))*vel_L(2,-2)*(1._wp/dx(j-2))
                                !$acc atomic
                                rhs_vf(E_idx)%sf(j-3,k,l) = rhs_vf(E_idx)%sf(j-3,k,l) + 0.5_wp*mu_L(-3)*((2._wp/60._wp)*(dvel1+dvel2))*vel_L(2,-3)*(1._wp/dx(j-3))

                                !$acc atomic
                                rhs_vf(momxb+1)%sf(j-1,k,l) = rhs_vf(momxb+1)%sf(j-1,k,l) - 0.5_wp*mu_R(-2)*((-3._wp/60._wp)*(dvel1+dvel2))*(1._wp/dx(j-1))
                                !$acc atomic
                                rhs_vf(momxb+1)%sf(j,k,l) = rhs_vf(momxb+1)%sf(j,k,l) - 0.5_wp*mu_R(-1)*((27._wp/60._wp)*(dvel1+dvel2))*(1._wp/dx(j))
                                !$acc atomic
                                rhs_vf(momxb+1)%sf(j+1,k,l) = rhs_vf(momxb+1)%sf(j+1,k,l) - 0.5_wp*mu_R(0)*((47._wp/60._wp)*(dvel1+dvel2))*(1._wp/dx(j+1))
                                !$acc atomic
                                rhs_vf(momxb+1)%sf(j+2,k,l) = rhs_vf(momxb+1)%sf(j+2,k,l) - 0.5_wp*mu_R(1)*((-13._wp/60._wp)*(dvel1+dvel2))*(1._wp/dx(j+2))
                                !$acc atomic
                                rhs_vf(momxb+1)%sf(j+3,k,l) = rhs_vf(momxb+1)%sf(j+3,k,l) - 0.5_wp*mu_R(2)*((2._wp/60._wp)*(dvel1+dvel2))*(1._wp/dx(j+3))

                                !$acc atomic
                                rhs_vf(E_idx)%sf(j-1,k,l) = rhs_vf(E_idx)%sf(j-1,k,l) - 0.5_wp*mu_R(-2)*((-3._wp/60._wp)*(dvel1+dvel2))*vel_R(2,-2)*(1._wp/dx(j-1))
                                !$acc atomic
                                rhs_vf(E_idx)%sf(j,k,l) = rhs_vf(E_idx)%sf(j,k,l) - 0.5_wp*mu_R(-1)*((27._wp/60._wp)*(dvel1+dvel2))*vel_R(2,-1)*(1._wp/dx(j))
                                !$acc atomic
                                rhs_vf(E_idx)%sf(j+1,k,l) = rhs_vf(E_idx)%sf(j+1,k,l) - 0.5_wp*mu_R(0)*((47._wp/60._wp)*(dvel1+dvel2))*vel_R(2,0)*(1._wp/dx(j+1))
                                !$acc atomic
                                rhs_vf(E_idx)%sf(j+2,k,l) = rhs_vf(E_idx)%sf(j+2,k,l) - 0.5_wp*mu_R(1)*((-13._wp/60._wp)*(dvel1+dvel2))*vel_R(2,1)*(1._wp/dx(j+2))
                                !$acc atomic
                                rhs_vf(E_idx)%sf(j+3,k,l) = rhs_vf(E_idx)%sf(j+3,k,l) - 0.5_wp*mu_R(2)*((2._wp/60._wp)*(dvel1+dvel2))*vel_R(2,2)*(1._wp/dx(j+3)) 

                                !$acc atomic
                                rhs_vf(momxb+1)%sf(j-2,k,l) = rhs_vf(momxb+1)%sf(j-2,k,l) + 0.5_wp*mu_R(-2)*((-3._wp/60._wp)*(dvel1+dvel2))*(1._wp/dx(j-2))
                                !$acc atomic
                                rhs_vf(momxb+1)%sf(j-1,k,l) = rhs_vf(momxb+1)%sf(j-1,k,l) + 0.5_wp*mu_R(-1)*((27._wp/60._wp)*(dvel1+dvel2))*(1._wp/dx(j-1))
                                !$acc atomic
                                rhs_vf(momxb+1)%sf(j,k,l) = rhs_vf(momxb+1)%sf(j,k,l) + 0.5_wp*mu_R(0)*((47._wp/60._wp)*(dvel1+dvel2))*(1._wp/dx(j))
                                !$acc atomic
                                rhs_vf(momxb+1)%sf(j+1,k,l) = rhs_vf(momxb+1)%sf(j+1,k,l) + 0.5_wp*mu_R(1)*((-13._wp/60._wp)*(dvel1+dvel2))*(1._wp/dx(j+1))
                                !$acc atomic
                                rhs_vf(momxb+1)%sf(j+2,k,l) = rhs_vf(momxb+1)%sf(j+2,k,l) + 0.5_wp*mu_R(2)*((2._wp/60._wp)*(dvel1+dvel2))*(1._wp/dx(j+2))

                                !$acc atomic
                                rhs_vf(E_idx)%sf(j-2,k,l) = rhs_vf(E_idx)%sf(j-2,k,l) + 0.5_wp*mu_R(-2)*((-3._wp/60._wp)*(dvel1+dvel2))*vel_R(2,-2)*(1._wp/dx(j-2))
                                !$acc atomic
                                rhs_vf(E_idx)%sf(j-1,k,l) = rhs_vf(E_idx)%sf(j-1,k,l) + 0.5_wp*mu_R(-1)*((27._wp/60._wp)*(dvel1+dvel2))*vel_R(2,-1)*(1._wp/dx(j-1))
                                !$acc atomic
                                rhs_vf(E_idx)%sf(j,k,l) = rhs_vf(E_idx)%sf(j,k,l) + 0.5_wp*mu_R(0)*((47._wp/60._wp)*(dvel1+dvel2))*vel_R(2,0)*(1._wp/dx(j))
                                !$acc atomic
                                rhs_vf(E_idx)%sf(j+1,k,l) = rhs_vf(E_idx)%sf(j+1,k,l) + 0.5_wp*mu_R(1)*((-13._wp/60._wp)*(dvel1+dvel2))*vel_R(2,1)*(1._wp/dx(j+1))
                                !$acc atomic
                                rhs_vf(E_idx)%sf(j+2,k,l) = rhs_vf(E_idx)%sf(j+2,k,l) + 0.5_wp*mu_R(2)*((2._wp/60._wp)*(dvel1+dvel2))*vel_R(2,2)*(1._wp/dx(j+2))

                            end if

                            !! duz & dwx
                            dvel1 = (1/(12._wp*dz(l))) * ( &
                            8._wp*q_prim_vf(momxb)%sf(j,k,l+1)/rho_sf(1,3) - &
                            8._wp*q_prim_vf(momxb)%sf(j,k,l-1)/rho_sf(-1,3) + &
                            q_prim_vf(momxb)%sf(j,k,l-2)/rho_sf(-2,3) - &
                            q_prim_vf(momxb)%sf(j,k,l+2)/rho_sf(2,3) )

                            dvel2 = (1/(12._wp*dx(j))) * ( &
                            8._wp*q_prim_vf(momxb+2)%sf(j+1,k,l)/rho_sf(1,1) - &
                            8._wp*q_prim_vf(momxb+2)%sf(j-1,k,l)/rho_sf(-1,1) + &
                            q_prim_vf(momxb+2)%sf(j-2,k,l)/rho_sf(-2,1) - &
                            q_prim_vf(momxb+2)%sf(j+2,k,l)/rho_sf(2,1) )

                            jac_rhs(j, k, l) = jac_rhs(j,k,l) + alf_igr* (2_wp * dvel1* dvel2)

                            if(viscous) then

                                !$acc atomic
                                rhs_vf(momxb+2)%sf(j+2,k,l) = rhs_vf(momxb+2)%sf(j+2,k,l) - 0.5_wp*mu_L(1)*((-3._wp/60._wp)*(dvel1+dvel2))*(1._wp/dx(j+2))
                                !$acc atomic
                                rhs_vf(momxb+2)%sf(j+1,k,l) = rhs_vf(momxb+2)%sf(j+1,k,l) - 0.5_wp*mu_L(0)*((27._wp/60._wp)*(dvel1+dvel2))*(1._wp/dx(j+1))
                                !$acc atomic
                                rhs_vf(momxb+2)%sf(j,k,l) = rhs_vf(momxb+2)%sf(j,k,l) - 0.5_wp*mu_L(-1)*((47._wp/60._wp)*(dvel1+dvel2))*(1._wp/dx(j))
                                !$acc atomic
                                rhs_vf(momxb+2)%sf(j-1,k,l) = rhs_vf(momxb+2)%sf(j-1,k,l) - 0.5_wp*mu_L(-2)*((-13._wp/60._wp)*(dvel1+dvel2))*(1._wp/dx(j-1))
                                !$acc atomic
                                rhs_vf(momxb+2)%sf(j-2,k,l) = rhs_vf(momxb+2)%sf(j-2,k,l) - 0.5_wp*mu_L(-3)*((2._wp/60._wp)*(dvel1+dvel2))*(1._wp/dx(j-2))

                                !$acc atomic
                                rhs_vf(E_idx)%sf(j+2,k,l) = rhs_vf(E_idx)%sf(j+2,k,l) - 0.5_wp*mu_L(1)*((-3._wp/60._wp)*(dvel1+dvel2))*vel_L(3,1)*(1._wp/dx(j+2))
                                !$acc atomic
                                rhs_vf(E_idx)%sf(j+1,k,l) = rhs_vf(E_idx)%sf(j+1,k,l) - 0.5_wp*mu_L(0)*((27._wp/60._wp)*(dvel1+dvel2))*vel_L(3,0)*(1._wp/dx(j+1))
                                !$acc atomic
                                rhs_vf(E_idx)%sf(j,k,l) = rhs_vf(E_idx)%sf(j,k,l) - 0.5_wp*mu_L(-1)*((47._wp/60._wp)*(dvel1+dvel2))*vel_L(3,-1)*(1._wp/dx(j))
                                !$acc atomic
                                rhs_vf(E_idx)%sf(j-1,k,l) = rhs_vf(E_idx)%sf(j-1,k,l) - 0.5_wp*mu_L(-2)*((-13._wp/60._wp)*(dvel1+dvel2))*vel_L(3,-2)*(1._wp/dx(j-1))
                                !$acc atomic
                                rhs_vf(E_idx)%sf(j-2,k,l) = rhs_vf(E_idx)%sf(j-2,k,l) - 0.5_wp*mu_L(-3)*((2._wp/60._wp)*(dvel1+dvel2))*vel_L(3,-3)*(1._wp/dx(j-2))

                                !$acc atomic
                                rhs_vf(momxb+2)%sf(j+1,k,l) = rhs_vf(momxb+2)%sf(j+1,k,l) + 0.5_wp*mu_L(1)*((-3._wp/60._wp)*(dvel1+dvel2))*(1._wp/dx(j+1))
                                !$acc atomic
                                rhs_vf(momxb+2)%sf(j,k,l) = rhs_vf(momxb+2)%sf(j,k,l) + 0.5_wp*mu_L(0)*((27._wp/60._wp)*(dvel1+dvel2))*(1._wp/dx(j))
                                !$acc atomic
                                rhs_vf(momxb+2)%sf(j-1,k,l) = rhs_vf(momxb+2)%sf(j-1,k,l) + 0.5_wp*mu_L(-1)*((47._wp/60._wp)*(dvel1+dvel2))*(1._wp/dx(j-1))
                                !$acc atomic
                                rhs_vf(momxb+2)%sf(j-2,k,l) = rhs_vf(momxb+2)%sf(j-2,k,l) + 0.5_wp*mu_L(-2)*((-13._wp/60._wp)*(dvel1+dvel2))*(1._wp/dx(j-2))
                                !$acc atomic
                                rhs_vf(momxb+2)%sf(j-3,k,l) = rhs_vf(momxb+2)%sf(j-3,k,l) + 0.5_wp*mu_L(-3)*((2._wp/60._wp)*(dvel1+dvel2))*(1._wp/dx(j-3))

                                !$acc atomic
                                rhs_vf(E_idx)%sf(j+1,k,l) = rhs_vf(E_idx)%sf(j+1,k,l) + 0.5_wp*mu_L(1)*((-3._wp/60._wp)*(dvel1+dvel2))*vel_L(3,1)*(1._wp/dx(j+1))
                                !$acc atomic
                                rhs_vf(E_idx)%sf(j,k,l) = rhs_vf(E_idx)%sf(j,k,l) + 0.5_wp*mu_L(0)*((27._wp/60._wp)*(dvel1+dvel2))*vel_L(3,0)*(1._wp/dx(j))
                                !$acc atomic
                                rhs_vf(E_idx)%sf(j-1,k,l) = rhs_vf(E_idx)%sf(j-1,k,l) + 0.5_wp*mu_L(-1)*((47._wp/60._wp)*(dvel1+dvel2))*vel_L(3,-1)*(1._wp/dx(j-1))
                                !$acc atomic
                                rhs_vf(E_idx)%sf(j-2,k,l) = rhs_vf(E_idx)%sf(j-2,k,l) + 0.5_wp*mu_L(-2)*((-13._wp/60._wp)*(dvel1+dvel2))*vel_L(3,-2)*(1._wp/dx(j-2))
                                !$acc atomic
                                rhs_vf(E_idx)%sf(j-3,k,l) = rhs_vf(E_idx)%sf(j-3,k,l) + 0.5_wp*mu_L(-3)*((2._wp/60._wp)*(dvel1+dvel2))*vel_L(3,-3)*(1._wp/dx(j-3))
                               
                                !$acc atomic
                                rhs_vf(momxb+2)%sf(j-1,k,l) = rhs_vf(momxb+2)%sf(j-1,k,l) - 0.5_wp*mu_R(-2)*((-3._wp/60._wp)*(dvel1+dvel2))*(1._wp/dx(j-1))
                                !$acc atomic
                                rhs_vf(momxb+2)%sf(j,k,l) = rhs_vf(momxb+2)%sf(j,k,l) - 0.5_wp*mu_R(-1)*((27._wp/60._wp)*(dvel1+dvel2))*(1._wp/dx(j))
                                !$acc atomic
                                rhs_vf(momxb+2)%sf(j+1,k,l) = rhs_vf(momxb+2)%sf(j+1,k,l) - 0.5_wp*mu_R(0)*((47._wp/60._wp)*(dvel1+dvel2))*(1._wp/dx(j+1))
                                !$acc atomic
                                rhs_vf(momxb+2)%sf(j+2,k,l) = rhs_vf(momxb+2)%sf(j+2,k,l) - 0.5_wp*mu_R(1)*((-13._wp/60._wp)*(dvel1+dvel2))*(1._wp/dx(j+2))
                                !$acc atomic
                                rhs_vf(momxb+2)%sf(j+3,k,l) = rhs_vf(momxb+2)%sf(j+3,k,l) - 0.5_wp*mu_R(2)*((2._wp/60._wp)*(dvel1+dvel2))*(1._wp/dx(j+3))

                                !$acc atomic
                                rhs_vf(E_idx)%sf(j-1,k,l) = rhs_vf(E_idx)%sf(j-1,k,l) - 0.5_wp*mu_R(-2)*((-3._wp/60._wp)*(dvel1+dvel2))*vel_R(3,-2)*(1._wp/dx(j-1))
                                !$acc atomic
                                rhs_vf(E_idx)%sf(j,k,l) = rhs_vf(E_idx)%sf(j,k,l) - 0.5_wp*mu_R(-1)*((27._wp/60._wp)*(dvel1+dvel2))*vel_R(3,-1)*(1._wp/dx(j))
                                !$acc atomic
                                rhs_vf(E_idx)%sf(j+1,k,l) = rhs_vf(E_idx)%sf(j+1,k,l) - 0.5_wp*mu_R(0)*((47._wp/60._wp)*(dvel1+dvel2))*vel_R(3,0)*(1._wp/dx(j+1))
                                !$acc atomic
                                rhs_vf(E_idx)%sf(j+2,k,l) = rhs_vf(E_idx)%sf(j+2,k,l) - 0.5_wp*mu_R(1)*((-13._wp/60._wp)*(dvel1+dvel2))*vel_R(3,1)*(1._wp/dx(j+2))
                                !$acc atomic
                                rhs_vf(E_idx)%sf(j+3,k,l) = rhs_vf(E_idx)%sf(j+3,k,l) - 0.5_wp*mu_R(2)*((2._wp/60._wp)*(dvel1+dvel2))*vel_R(3,2)*(1._wp/dx(j+3)) 

                                !$acc atomic
                                rhs_vf(momxb+2)%sf(j-2,k,l) = rhs_vf(momxb+2)%sf(j-2,k,l) + 0.5_wp*mu_R(-2)*((-3._wp/60._wp)*(dvel1+dvel2))*(1._wp/dx(j-2))
                                !$acc atomic
                                rhs_vf(momxb+2)%sf(j-1,k,l) = rhs_vf(momxb+2)%sf(j-1,k,l) + 0.5_wp*mu_R(-1)*((27._wp/60._wp)*(dvel1+dvel2))*(1._wp/dx(j-1))
                                !$acc atomic
                                rhs_vf(momxb+2)%sf(j,k,l) = rhs_vf(momxb+2)%sf(j,k,l) + 0.5_wp*mu_R(0)*((47._wp/60._wp)*(dvel1+dvel2))*(1._wp/dx(j))
                                !$acc atomic
                                rhs_vf(momxb+2)%sf(j+1,k,l) = rhs_vf(momxb+2)%sf(j+1,k,l) + 0.5_wp*mu_R(1)*((-13._wp/60._wp)*(dvel1+dvel2))*(1._wp/dx(j+1))
                                !$acc atomic
                                rhs_vf(momxb+2)%sf(j+2,k,l) = rhs_vf(momxb+2)%sf(j+2,k,l) + 0.5_wp*mu_R(2)*((2._wp/60._wp)*(dvel1+dvel2))*(1._wp/dx(j+2))

                                !$acc atomic
                                rhs_vf(E_idx)%sf(j-2,k,l) = rhs_vf(E_idx)%sf(j-2,k,l) + 0.5_wp*mu_R(-2)*((-3._wp/60._wp)*(dvel1+dvel2))*vel_R(3,-2)*(1._wp/dx(j-2))
                                !$acc atomic
                                rhs_vf(E_idx)%sf(j-1,k,l) = rhs_vf(E_idx)%sf(j-1,k,l) + 0.5_wp*mu_R(-1)*((27._wp/60._wp)*(dvel1+dvel2))*vel_R(3,-1)*(1._wp/dx(j-1))
                                !$acc atomic
                                rhs_vf(E_idx)%sf(j,k,l) = rhs_vf(E_idx)%sf(j,k,l) + 0.5_wp*mu_R(0)*((47._wp/60._wp)*(dvel1+dvel2))*vel_R(3,0)*(1._wp/dx(j))
                                !$acc atomic
                                rhs_vf(E_idx)%sf(j+1,k,l) = rhs_vf(E_idx)%sf(j+1,k,l) + 0.5_wp*mu_R(1)*((-13._wp/60._wp)*(dvel1+dvel2))*vel_R(3,1)*(1._wp/dx(j+1))
                                !$acc atomic
                                rhs_vf(E_idx)%sf(j+2,k,l) = rhs_vf(E_idx)%sf(j+2,k,l) + 0.5_wp*mu_R(2)*((2._wp/60._wp)*(dvel1+dvel2))*vel_R(3,2)*(1._wp/dx(j+2))

                            end if

                            !! dvz & dwy
                            dvel1 = (1/(12._wp*dz(l))) * ( &
                            8._wp*q_prim_vf(momxb+1)%sf(j,k,l+1)/rho_sf(1,3) - &
                            8._wp*q_prim_vf(momxb+1)%sf(j,k,l-1)/rho_sf(-1,3) + &
                            q_prim_vf(momxb+1)%sf(j,k,l-2)/rho_sf(-2,3) - &
                            q_prim_vf(momxb+1)%sf(j,k,l+2)/rho_sf(2,3) )

                            dvel2 = (1/(12._wp*dy(k))) * ( &
                            8._wp*q_prim_vf(momxb+2)%sf(j,k+1,l)/rho_sf(1,2) - &
                            8._wp*q_prim_vf(momxb+2)%sf(j,k-1,l)/rho_sf(-1,2) + &
                            q_prim_vf(momxb+2)%sf(j,k-2,l)/rho_sf(-2,2) - &
                            q_prim_vf(momxb+2)%sf(j,k+2,l)/rho_sf(2,2))

                            jac_rhs(j,k,l) = jac_rhs(j,k,l) + alf_igr*(2_wp*dvel1*dvel2)

                            !! dux & dvy
                            dvel1 = (1/(12._wp*dx(j))) * ( &
                            8._wp*q_prim_vf(momxb)%sf(j+1,k,l)/rho_sf(1,1) - &
                            8._wp*q_prim_vf(momxb)%sf(j-1,k,l)/rho_sf(-1,1) + &
                            q_prim_vf(momxb)%sf(j-2,k,l)/rho_sf(-2,1) - &
                            q_prim_vf(momxb)%sf(j+2,k,l)/rho_sf(2,1) )


                            dvel2 = (1/(12._wp*dy(k))) * ( &
                            8._wp*q_prim_vf(momxb+1)%sf(j,k+1,l)/rho_sf(1,2) - &
                            8._wp*q_prim_vf(momxb+1)%sf(j,k-1,l)/rho_sf(-1,2) + &
                            q_prim_vf(momxb+1)%sf(j,k-2,l)/rho_sf(-2,2) - &
                            q_prim_vf(momxb+1)%sf(j,k+2,l)/rho_sf(2,2) )

                            jac_rhs(j, k, l) = jac_rhs(j,k,l) + alf_igr* (dvel1**2_wp + dvel2**2_wp)

                            if(viscous) then 

                                !$acc atomic
                                rhs_vf(momxb)%sf(j+2,k,l) = rhs_vf(momxb)%sf(j+2,k,l) - 0.5_wp*mu_L(1)*((-3._wp/60._wp)*(4._wp*dvel1-2._wp*dvel2)/3._wp)*(1._wp/dx(j+2))
                                !$acc atomic
                                rhs_vf(momxb)%sf(j+1,k,l) = rhs_vf(momxb)%sf(j+1,k,l) - 0.5_wp*mu_L(0)*((27._wp/60._wp)*(4._wp*dvel1-2._wp*dvel2)/3._wp)*(1._wp/dx(j+1))
                                !$acc atomic
                                rhs_vf(momxb)%sf(j,k,l) = rhs_vf(momxb)%sf(j,k,l) - 0.5_wp*mu_L(-1)*((47._wp/60._wp)*(4._wp*dvel1-2._wp*dvel2)/3._wp)*(1._wp/dx(j))
                                !$acc atomic
                                rhs_vf(momxb)%sf(j-1,k,l) = rhs_vf(momxb)%sf(j-1,k,l) - 0.5_wp*mu_L(-2)*((-13._wp/60._wp)*(4._wp*dvel1-2._wp*dvel2)/3._wp)*(1._wp/dx(j-1))
                                !$acc atomic
                                rhs_vf(momxb)%sf(j-2,k,l) = rhs_vf(momxb)%sf(j-2,k,l) - 0.5_wp*mu_L(-3)*((2._wp/60._wp)*(4._wp*dvel1-2._wp*dvel2)/3._wp)*(1._wp/dx(j-2))

                                !$acc atomic
                                rhs_vf(E_idx)%sf(j+2,k,l) = rhs_vf(E_idx)%sf(j+2,k,l) - 0.5_wp*mu_L(1)*((-3._wp/60._wp)*(4._wp*dvel1-2._wp*dvel2)/3._wp)*vel_L(1,1)*(1._wp/dx(j+2))
                                !$acc atomic
                                rhs_vf(E_idx)%sf(j+1,k,l) = rhs_vf(E_idx)%sf(j+1,k,l) - 0.5_wp*mu_L(0)*((27._wp/60._wp)*(4._wp*dvel1-2._wp*dvel2)/3._wp)*vel_L(1,0)*(1._wp/dx(j+1))
                                !$acc atomic
                                rhs_vf(E_idx)%sf(j,k,l) = rhs_vf(E_idx)%sf(j,k,l) - 0.5_wp*mu_L(-1)*((47._wp/60._wp)*(4._wp*dvel1-2._wp*dvel2)/3._wp)*vel_L(1,-1)*(1._wp/dx(j))
                                !$acc atomic
                                rhs_vf(E_idx)%sf(j-1,k,l) = rhs_vf(E_idx)%sf(j-1,k,l) - 0.5_wp*mu_L(-2)*((-13._wp/60._wp)*(4._wp*dvel1-2._wp*dvel2)/3._wp)*vel_L(1,-2)*(1._wp/dx(j-1))
                                !$acc atomic
                                rhs_vf(E_idx)%sf(j-2,k,l) = rhs_vf(E_idx)%sf(j-2,k,l) - 0.5_wp*mu_L(-3)*((2._wp/60._wp)*(4._wp*dvel1-2._wp*dvel2)/3._wp)*vel_L(1,-3)*(1._wp/dx(j-2))

                                !$acc atomic
                                rhs_vf(momxb)%sf(j+1,k,l) = rhs_vf(momxb)%sf(j+1,k,l) + 0.5_wp*mu_L(1)*((-3._wp/60._wp)*(4._wp*dvel1-2._wp*dvel2)/3._wp)*(1._wp/dx(j+1))
                                !$acc atomic
                                rhs_vf(momxb)%sf(j,k,l) = rhs_vf(momxb)%sf(j,k,l) + 0.5_wp*mu_L(0)*((27._wp/60._wp)*(4._wp*dvel1-2._wp*dvel2)/3._wp)*(1._wp/dx(j))
                                !$acc atomic
                                rhs_vf(momxb)%sf(j-1,k,l) = rhs_vf(momxb)%sf(j-1,k,l) + 0.5_wp*mu_L(-1)*((47._wp/60._wp)*(4._wp*dvel1-2._wp*dvel2)/3._wp)*(1._wp/dx(j-1))
                                !$acc atomic
                                rhs_vf(momxb)%sf(j-2,k,l) = rhs_vf(momxb)%sf(j-2,k,l) + 0.5_wp*mu_L(-2)*((-13._wp/60._wp)*(4._wp*dvel1-2._wp*dvel2)/3._wp)*(1._wp/dx(j-2))
                                !$acc atomic
                                rhs_vf(momxb)%sf(j-3,k,l) = rhs_vf(momxb)%sf(j-3,k,l) + 0.5_wp*mu_L(-3)*((2._wp/60._wp)*(4._wp*dvel1-2._wp*dvel2)/3._wp)*(1._wp/dx(j-3))

                                !$acc atomic
                                rhs_vf(E_idx)%sf(j+1,k,l) = rhs_vf(E_idx)%sf(j+1,k,l) + 0.5_wp*mu_L(1)*((-3._wp/60._wp)*(4._wp*dvel1-2._wp*dvel2)/3._wp)*vel_L(1,1)*(1._wp/dx(j+1))
                                !$acc atomic
                                rhs_vf(E_idx)%sf(j,k,l) = rhs_vf(E_idx)%sf(j,k,l) + 0.5_wp*mu_L(0)*((27._wp/60._wp)*(4._wp*dvel1-2._wp*dvel2)/3._wp)*vel_L(1,0)*(1._wp/dx(j))
                                !$acc atomic
                                rhs_vf(E_idx)%sf(j-1,k,l) = rhs_vf(E_idx)%sf(j-1,k,l) + 0.5_wp*mu_L(-1)*((47._wp/60._wp)*(4._wp*dvel1-2._wp*dvel2)/3._wp)*vel_L(1,-1)*(1._wp/dx(j-1))
                                !$acc atomic
                                rhs_vf(E_idx)%sf(j-2,k,l) = rhs_vf(E_idx)%sf(j-2,k,l) + 0.5_wp*mu_L(-2)*((-13._wp/60._wp)*(4._wp*dvel1-2._wp*dvel2)/3._wp)*vel_L(1,-2)*(1._wp/dx(j-2))
                                !$acc atomic
                                rhs_vf(E_idx)%sf(j-3,k,l) = rhs_vf(E_idx)%sf(j-3,k,l) + 0.5_wp*mu_L(-3)*((2._wp/60._wp)*(4._wp*dvel1-2._wp*dvel2)/3._wp)*vel_L(1,-3)*(1._wp/dx(j-3))
                               
                                !$acc atomic
                                rhs_vf(momxb)%sf(j-1,k,l) = rhs_vf(momxb)%sf(j-1,k,l) - 0.5_wp*mu_R(-2)*((-3._wp/60._wp)*(4._wp*dvel1-2._wp*dvel2)/3._wp)*(1._wp/dx(j-1))
                                !$acc atomic
                                rhs_vf(momxb)%sf(j,k,l) = rhs_vf(momxb)%sf(j,k,l) - 0.5_wp*mu_R(-1)*((27._wp/60._wp)*(4._wp*dvel1-2._wp*dvel2)/3._wp)*(1._wp/dx(j))
                                !$acc atomic
                                rhs_vf(momxb)%sf(j+1,k,l) = rhs_vf(momxb)%sf(j+1,k,l) - 0.5_wp*mu_R(0)*((47._wp/60._wp)*(4._wp*dvel1-2._wp*dvel2)/3._wp)*(1._wp/dx(j+1))
                                !$acc atomic
                                rhs_vf(momxb)%sf(j+2,k,l) = rhs_vf(momxb)%sf(j+2,k,l) - 0.5_wp*mu_R(1)*((-13._wp/60._wp)*(4._wp*dvel1-2._wp*dvel2)/3._wp)*(1._wp/dx(j+2))
                                !$acc atomic
                                rhs_vf(momxb)%sf(j+3,k,l) = rhs_vf(momxb)%sf(j+3,k,l) - 0.5_wp*mu_R(2)*((2._wp/60._wp)*(4._wp*dvel1-2._wp*dvel2)/3._wp)*(1._wp/dx(j+3))

                                !$acc atomic
                                rhs_vf(E_idx)%sf(j-1,k,l) = rhs_vf(E_idx)%sf(j-1,k,l) - 0.5_wp*mu_R(-2)*((-3._wp/60._wp)*(4._wp*dvel1-2._wp*dvel2)/3._wp)*vel_R(1,-2)*(1._wp/dx(j-1))
                                !$acc atomic
                                rhs_vf(E_idx)%sf(j,k,l) = rhs_vf(E_idx)%sf(j,k,l) - 0.5_wp*mu_R(-1)*((27._wp/60._wp)*(4._wp*dvel1-2._wp*dvel2)/3._wp)*vel_R(1,-1)*(1._wp/dx(j))
                                !$acc atomic
                                rhs_vf(E_idx)%sf(j+1,k,l) = rhs_vf(E_idx)%sf(j+1,k,l) - 0.5_wp*mu_R(0)*((47._wp/60._wp)*(4._wp*dvel1-2._wp*dvel2)/3._wp)*vel_R(1,0)*(1._wp/dx(j+1))
                                !$acc atomic
                                rhs_vf(E_idx)%sf(j+2,k,l) = rhs_vf(E_idx)%sf(j+2,k,l) - 0.5_wp*mu_R(1)*((-13._wp/60._wp)*(4._wp*dvel1-2._wp*dvel2)/3._wp)*vel_R(1,1)*(1._wp/dx(j+2))
                                !$acc atomic
                                rhs_vf(E_idx)%sf(j+3,k,l) = rhs_vf(E_idx)%sf(j+3,k,l) - 0.5_wp*mu_R(2)*((2._wp/60._wp)*(4._wp*dvel1-2._wp*dvel2)/3._wp)*vel_R(1,2)*(1._wp/dx(j+3)) 

                                !$acc atomic
                                rhs_vf(momxb)%sf(j-2,k,l) = rhs_vf(momxb)%sf(j-2,k,l) + 0.5_wp*mu_R(-2)*((-3._wp/60._wp)*(4._wp*dvel1-2._wp*dvel2)/3._wp)*(1._wp/dx(j-2))
                                !$acc atomic
                                rhs_vf(momxb)%sf(j-1,k,l) = rhs_vf(momxb)%sf(j-1,k,l) + 0.5_wp*mu_R(-1)*((27._wp/60._wp)*(4._wp*dvel1-2._wp*dvel2)/3._wp)*(1._wp/dx(j-1))
                                !$acc atomic
                                rhs_vf(momxb)%sf(j,k,l) = rhs_vf(momxb)%sf(j,k,l) + 0.5_wp*mu_R(0)*((47._wp/60._wp)*(4._wp*dvel1-2._wp*dvel2)/3._wp)*(1._wp/dx(j))
                                !$acc atomic
                                rhs_vf(momxb)%sf(j+1,k,l) = rhs_vf(momxb)%sf(j+1,k,l) + 0.5_wp*mu_R(1)*((-13._wp/60._wp)*(4._wp*dvel1-2._wp*dvel2)/3._wp)*(1._wp/dx(j+1))
                                !$acc atomic
                                rhs_vf(momxb)%sf(j+2,k,l) = rhs_vf(momxb)%sf(j+2,k,l) + 0.5_wp*mu_R(2)*((2._wp/60._wp)*(4._wp*dvel1-2._wp*dvel2)/3._wp)*(1._wp/dx(j+2))

                                !$acc atomic
                                rhs_vf(E_idx)%sf(j-2,k,l) = rhs_vf(E_idx)%sf(j-2,k,l) + 0.5_wp*mu_R(-2)*((-3._wp/60._wp)*(4._wp*dvel1-2._wp*dvel2)/3._wp)*vel_R(1,-2)*(1._wp/dx(j-2))
                                !$acc atomic
                                rhs_vf(E_idx)%sf(j-1,k,l) = rhs_vf(E_idx)%sf(j-1,k,l) + 0.5_wp*mu_R(-1)*((27._wp/60._wp)*(4._wp*dvel1-2._wp*dvel2)/3._wp)*vel_R(1,-1)*(1._wp/dx(j-1))
                                !$acc atomic
                                rhs_vf(E_idx)%sf(j,k,l) = rhs_vf(E_idx)%sf(j,k,l) + 0.5_wp*mu_R(0)*((47._wp/60._wp)*(4._wp*dvel1-2._wp*dvel2)/3._wp)*vel_R(1,0)*(1._wp/dx(j))
                                !$acc atomic
                                rhs_vf(E_idx)%sf(j+1,k,l) = rhs_vf(E_idx)%sf(j+1,k,l) + 0.5_wp*mu_R(1)*((-13._wp/60._wp)*(4._wp*dvel1-2._wp*dvel2)/3._wp)*vel_R(1,1)*(1._wp/dx(j+1))
                                !$acc atomic
                                rhs_vf(E_idx)%sf(j+2,k,l) = rhs_vf(E_idx)%sf(j+2,k,l) + 0.5_wp*mu_R(2)*((2._wp/60._wp)*(4._wp*dvel1-2._wp*dvel2)/3._wp)*vel_R(1,2)*(1._wp/dx(j+2))

                            end if

                            !! dux + dvy, dwz

                            dvel1 = dvel1 + dvel2

                            dvel2 = (1/(12._wp*dz(l))) * ( &
                            8._wp*q_prim_vf(momxb+2)%sf(j,k,l+1)/rho_sf(1,3) - &
                            8._wp*q_prim_vf(momxb+2)%sf(j,k,l-1)/rho_sf(-1,3) + &
                            q_prim_vf(momxb+2)%sf(j,k,l-2)/rho_sf(-2,3) - &
                            q_prim_vf(momxb+2)%sf(j,k,l+2)/rho_sf(2,3) )

                            jac_rhs(j,k,l) = jac_rhs(j,k,l) + alf_igr * (dvel2 **2_wp + (dvel1 + dvel2)**2_wp)

                            if(viscous) then 

                                !$acc atomic
                                rhs_vf(momxb)%sf(j+2,k,l) = rhs_vf(momxb)%sf(j+2,k,l) - 0.5_wp*mu_L(1)*((-3._wp/60._wp)*(-2._wp*dvel2)/3._wp)*(1._wp/dx(j+2))
                                !$acc atomic
                                rhs_vf(momxb)%sf(j+1,k,l) = rhs_vf(momxb)%sf(j+1,k,l) - 0.5_wp*mu_L(0)*((27._wp/60._wp)*(-2._wp*dvel2)/3._wp)*(1._wp/dx(j+1))
                                !$acc atomic
                                rhs_vf(momxb)%sf(j,k,l) = rhs_vf(momxb)%sf(j,k,l) - 0.5_wp*mu_L(-1)*((47._wp/60._wp)*(-2._wp*dvel2)/3._wp)*(1._wp/dx(j))
                                !$acc atomic
                                rhs_vf(momxb)%sf(j-1,k,l) = rhs_vf(momxb)%sf(j-1,k,l) - 0.5_wp*mu_L(-2)*((-13._wp/60._wp)*(-2._wp*dvel2)/3._wp)*(1._wp/dx(j-1))
                                !$acc atomic
                                rhs_vf(momxb)%sf(j-2,k,l) = rhs_vf(momxb)%sf(j-2,k,l) - 0.5_wp*mu_L(-3)*((2._wp/60._wp)*(-2._wp*dvel2)/3._wp)*(1._wp/dx(j-2))

                                !$acc atomic
                                rhs_vf(E_idx)%sf(j+2,k,l) = rhs_vf(E_idx)%sf(j+2,k,l) - 0.5_wp*mu_L(1)*((-3._wp/60._wp)*(-2._wp*dvel2)/3._wp)*vel_L(1,1)*(1._wp/dx(j+2))
                                !$acc atomic
                                rhs_vf(E_idx)%sf(j+1,k,l) = rhs_vf(E_idx)%sf(j+1,k,l) - 0.5_wp*mu_L(0)*((27._wp/60._wp)*(-2._wp*dvel2)/3._wp)*vel_L(1,0)*(1._wp/dx(j+1))
                                !$acc atomic
                                rhs_vf(E_idx)%sf(j,k,l) = rhs_vf(E_idx)%sf(j,k,l) - 0.5_wp*mu_L(-1)*((47._wp/60._wp)*(-2._wp*dvel2)/3._wp)*vel_L(1,-1)*(1._wp/dx(j))
                                !$acc atomic
                                rhs_vf(E_idx)%sf(j-1,k,l) = rhs_vf(E_idx)%sf(j-1,k,l) - 0.5_wp*mu_L(-2)*((-13._wp/60._wp)*(-2._wp*dvel2)/3._wp)*vel_L(1,-2)*(1._wp/dx(j-1))
                                !$acc atomic
                                rhs_vf(E_idx)%sf(j-2,k,l) = rhs_vf(E_idx)%sf(j-2,k,l) - 0.5_wp*mu_L(-3)*((2._wp/60._wp)*(-2._wp*dvel2)/3._wp)*vel_L(1,-3)*(1._wp/dx(j-2))

                                !$acc atomic
                                rhs_vf(momxb)%sf(j+1,k,l) = rhs_vf(momxb)%sf(j+1,k,l) + 0.5_wp*mu_L(1)*((-3._wp/60._wp)*(-2._wp*dvel2)/3._wp)*(1._wp/dx(j+1))
                                !$acc atomic
                                rhs_vf(momxb)%sf(j,k,l) = rhs_vf(momxb)%sf(j,k,l) + 0.5_wp*mu_L(0)*((27._wp/60._wp)*(-2._wp*dvel2)/3._wp)*(1._wp/dx(j))
                                !$acc atomic
                                rhs_vf(momxb)%sf(j-1,k,l) = rhs_vf(momxb)%sf(j-1,k,l) + 0.5_wp*mu_L(-1)*((47._wp/60._wp)*(-2._wp*dvel2)/3._wp)*(1._wp/dx(j-1))
                                !$acc atomic
                                rhs_vf(momxb)%sf(j-2,k,l) = rhs_vf(momxb)%sf(j-2,k,l) + 0.5_wp*mu_L(-2)*((-13._wp/60._wp)*(-2._wp*dvel2)/3._wp)*(1._wp/dx(j-2))
                                !$acc atomic
                                rhs_vf(momxb)%sf(j-3,k,l) = rhs_vf(momxb)%sf(j-3,k,l) + 0.5_wp*mu_L(-3)*((2._wp/60._wp)*(-2._wp*dvel2)/3._wp)*(1._wp/dx(j-3))

                                !$acc atomic
                                rhs_vf(E_idx)%sf(j+1,k,l) = rhs_vf(E_idx)%sf(j+1,k,l) + 0.5_wp*mu_L(1)*((-3._wp/60._wp)*(-2._wp*dvel2)/3._wp)*vel_L(1,1)*(1._wp/dx(j+1))
                                !$acc atomic
                                rhs_vf(E_idx)%sf(j,k,l) = rhs_vf(E_idx)%sf(j,k,l) + 0.5_wp*mu_L(0)*((27._wp/60._wp)*(-2._wp*dvel2)/3._wp)*vel_L(1,0)*(1._wp/dx(j))
                                !$acc atomic
                                rhs_vf(E_idx)%sf(j-1,k,l) = rhs_vf(E_idx)%sf(j-1,k,l) + 0.5_wp*mu_L(-1)*((47._wp/60._wp)*(-2._wp*dvel2)/3._wp)*vel_L(1,-1)*(1._wp/dx(j-1))
                                !$acc atomic
                                rhs_vf(E_idx)%sf(j-2,k,l) = rhs_vf(E_idx)%sf(j-2,k,l) + 0.5_wp*mu_L(-2)*((-13._wp/60._wp)*(-2._wp*dvel2)/3._wp)*vel_L(1,-2)*(1._wp/dx(j-2))
                                !$acc atomic
                                rhs_vf(E_idx)%sf(j-3,k,l) = rhs_vf(E_idx)%sf(j-3,k,l) + 0.5_wp*mu_L(-3)*((2._wp/60._wp)*(-2._wp*dvel2)/3._wp)*vel_L(1,-3)*(1._wp/dx(j-3))
                               
                                !$acc atomic
                                rhs_vf(momxb)%sf(j-1,k,l) = rhs_vf(momxb)%sf(j-1,k,l) - 0.5_wp*mu_R(-2)*((-3._wp/60._wp)*(-2._wp*dvel2)/3._wp)*(1._wp/dx(j-1))
                                !$acc atomic
                                rhs_vf(momxb)%sf(j,k,l) = rhs_vf(momxb)%sf(j,k,l) - 0.5_wp*mu_R(-1)*((27._wp/60._wp)*(-2._wp*dvel2)/3._wp)*(1._wp/dx(j))
                                !$acc atomic
                                rhs_vf(momxb)%sf(j+1,k,l) = rhs_vf(momxb)%sf(j+1,k,l) - 0.5_wp*mu_R(0)*((47._wp/60._wp)*(-2._wp*dvel2)/3._wp)*(1._wp/dx(j+1))
                                !$acc atomic
                                rhs_vf(momxb)%sf(j+2,k,l) = rhs_vf(momxb)%sf(j+2,k,l) - 0.5_wp*mu_R(1)*((-13._wp/60._wp)*(-2._wp*dvel2)/3._wp)*(1._wp/dx(j+2))
                                !$acc atomic
                                rhs_vf(momxb)%sf(j+3,k,l) = rhs_vf(momxb)%sf(j+3,k,l) - 0.5_wp*mu_R(2)*((2._wp/60._wp)*(-2._wp*dvel2)/3._wp)*(1._wp/dx(j+3))

                                !$acc atomic
                                rhs_vf(E_idx)%sf(j-1,k,l) = rhs_vf(E_idx)%sf(j-1,k,l) - 0.5_wp*mu_R(-2)*((-3._wp/60._wp)*(-2._wp*dvel2)/3._wp)*vel_R(1,-2)*(1._wp/dx(j-1))
                                !$acc atomic
                                rhs_vf(E_idx)%sf(j,k,l) = rhs_vf(E_idx)%sf(j,k,l) - 0.5_wp*mu_R(-1)*((27._wp/60._wp)*(-2._wp*dvel2)/3._wp)*vel_R(1,-1)*(1._wp/dx(j))
                                !$acc atomic
                                rhs_vf(E_idx)%sf(j+1,k,l) = rhs_vf(E_idx)%sf(j+1,k,l) - 0.5_wp*mu_R(0)*((47._wp/60._wp)*(-2._wp*dvel2)/3._wp)*vel_R(1,0)*(1._wp/dx(j+1))
                                !$acc atomic
                                rhs_vf(E_idx)%sf(j+2,k,l) = rhs_vf(E_idx)%sf(j+2,k,l) - 0.5_wp*mu_R(1)*((-13._wp/60._wp)*(-2._wp*dvel2)/3._wp)*vel_R(1,1)*(1._wp/dx(j+2))
                                !$acc atomic
                                rhs_vf(E_idx)%sf(j+3,k,l) = rhs_vf(E_idx)%sf(j+3,k,l) - 0.5_wp*mu_R(2)*((2._wp/60._wp)*(-2._wp*dvel2)/3._wp)*vel_R(1,2)*(1._wp/dx(j+3)) 

                                !$acc atomic
                                rhs_vf(momxb)%sf(j-2,k,l) = rhs_vf(momxb)%sf(j-2,k,l) + 0.5_wp*mu_R(-2)*((-3._wp/60._wp)*(-2._wp*dvel2)/3._wp)*(1._wp/dx(j-2))
                                !$acc atomic
                                rhs_vf(momxb)%sf(j-1,k,l) = rhs_vf(momxb)%sf(j-1,k,l) + 0.5_wp*mu_R(-1)*((27._wp/60._wp)*(-2._wp*dvel2)/3._wp)*(1._wp/dx(j-1))
                                !$acc atomic
                                rhs_vf(momxb)%sf(j,k,l) = rhs_vf(momxb)%sf(j,k,l) + 0.5_wp*mu_R(0)*((47._wp/60._wp)*(-2._wp*dvel2)/3._wp)*(1._wp/dx(j))
                                !$acc atomic
                                rhs_vf(momxb)%sf(j+1,k,l) = rhs_vf(momxb)%sf(j+1,k,l) + 0.5_wp*mu_R(1)*((-13._wp/60._wp)*(-2._wp*dvel2)/3._wp)*(1._wp/dx(j+1))
                                !$acc atomic
                                rhs_vf(momxb)%sf(j+2,k,l) = rhs_vf(momxb)%sf(j+2,k,l) + 0.5_wp*mu_R(2)*((2._wp/60._wp)*(-2._wp*dvel2)/3._wp)*(1._wp/dx(j+2))

                                !$acc atomic
                                rhs_vf(E_idx)%sf(j-2,k,l) = rhs_vf(E_idx)%sf(j-2,k,l) + 0.5_wp*mu_R(-2)*((-3._wp/60._wp)*(-2._wp*dvel2)/3._wp)*vel_R(1,-2)*(1._wp/dx(j-2))
                                !$acc atomic
                                rhs_vf(E_idx)%sf(j-1,k,l) = rhs_vf(E_idx)%sf(j-1,k,l) + 0.5_wp*mu_R(-1)*((27._wp/60._wp)*(-2._wp*dvel2)/3._wp)*vel_R(1,-1)*(1._wp/dx(j-1))
                                !$acc atomic
                                rhs_vf(E_idx)%sf(j,k,l) = rhs_vf(E_idx)%sf(j,k,l) + 0.5_wp*mu_R(0)*((47._wp/60._wp)*(-2._wp*dvel2)/3._wp)*vel_R(1,0)*(1._wp/dx(j))
                                !$acc atomic
                                rhs_vf(E_idx)%sf(j+1,k,l) = rhs_vf(E_idx)%sf(j+1,k,l) + 0.5_wp*mu_R(1)*((-13._wp/60._wp)*(-2._wp*dvel2)/3._wp)*vel_R(1,1)*(1._wp/dx(j+1))
                                !$acc atomic
                                rhs_vf(E_idx)%sf(j+2,k,l) = rhs_vf(E_idx)%sf(j+2,k,l) + 0.5_wp*mu_R(2)*((2._wp/60._wp)*(-2._wp*dvel2)/3._wp)*vel_R(1,2)*(1._wp/dx(j+2))

                            end if

                            !$acc loop seq
                            do i = 1, num_fluids
                                !$acc atomic
                                rhs_vf(i)%sf(j+1,k,l) = rhs_vf(i)%sf(j+1,k,l) + &
                                    (0.5_wp * (alpha_rho_R(i) * &
                                    vel_R(1,0))*(1._wp/dx(j+1)) + &
                                    0.5_wp*cfl * (alpha_rho_R(i))*(1._wp/dx(j+1)))

                                !$acc atomic
                                rhs_vf(i)%sf(j,k,l) = rhs_vf(i)%sf(j,k,l) - &
                                    (0.5_wp * (alpha_rho_R(i) * &
                                    vel_R(1,0))*(1._wp/dx(j)) + &
                                    0.5_wp*cfl * (alpha_rho_R(i))*(1._wp/dx(j))) 
                            end do

                            if(num_fluids > 1) then 
                                !$acc loop seq
                                do i = 1, num_fluids
                                    !$acc atomic
                                    rhs_vf(advxb+i-1)%sf(j+1,k,l) = rhs_vf(advxb+i-1)%sf(j+1,k,l) + &
                                        (0.5_wp * (alpha_R(i) * &
                                        vel_R(1,0))*(1._wp/dx(j+1)) + &
                                        0.5_wp*cfl*(alpha_R(i))*(1._wp/dx(j+1)))

                                    !$acc atomic
                                    rhs_vf(advxb+i-1)%sf(j+1,k,l) = rhs_vf(advxb+i-1)%sf(j+1,k,l) &
                                    - (0.5_wp * q_prim_vf(advxb+i-1)%sf(j+1,k,l)  * vel_R(1,0)*(1._wp/dx(j+1)))

                                    !$acc atomic
                                    rhs_vf(advxb+i-1)%sf(j,k,l) = rhs_vf(advxb+i-1)%sf(j,k,l) - &
                                        (0.5_wp * (alpha_R(i) * &
                                        vel_R(1,0))*(1._wp/dx(j)) + &
                                        0.5_wp*cfl*(alpha_R(i))*(1._wp/dx(j)))

                                    !$acc atomic
                                    rhs_vf(advxb+i-1)%sf(j,k,l) = rhs_vf(advxb+i-1)%sf(j,k,l) &
                                    + (0.5_wp * q_prim_vf(advxb+i-1)%sf(j,k,l) * vel_R(1,0)*(1._wp/dx(j)))
                                end do
                            end if

                            !$acc atomic
                            rhs_vf(momxb)%sf(j+1,k,l) = rhs_vf(momxb)%sf(j+1,k,l) + &
                                 (0.5_wp* (rho_R * (vel_R(1,0))**2.0 + &
                                 pres_R)*(1._wp/dx(j+1)) + &
                                 0.5_wp*cfl * (rho_R*vel_R(1,0))*(1._wp/dx(j+1)))

                            !$acc atomic
                            rhs_vf(momxb+1)%sf(j+1, k, l) =  rhs_vf(momxb+1)%sf(j+1, k, l) + &
                                (0.5_wp * rho_R * vel_R(1,0)*vel_R(2,0)*(1._wp/dx(j+1)) + &
                                0.5_wp*cfl * (rho_R*vel_R(2,0))*(1._wp/dx(j+1)))

                            !$acc atomic
                            rhs_vf(momxb+2)%sf(j+1, k, l) =  rhs_vf(momxb+2)%sf(j+1, k, l) + &
                                (0.5_wp * rho_R * vel_R(1,0)*vel_R(3,0)*(1._wp/dx(j+1)) + &
                                0.5_wp*cfl * (rho_R*vel_R(3,0))*(1._wp/dx(j+1)))

                            !$acc atomic
                             rhs_vf(E_idx)%sf(j+1, k, l) = rhs_vf(E_idx)%sf(j+1, k, l) +  &
                                (0.5_wp * (vel_R(1,0) * (E_R + &
                                pres_R) )*(1._wp/dx(j+1)) + &
                                0.5_wp*cfl * (E_R)*(1._wp/dx(j+1)))

                            !$acc atomic
                            rhs_vf(momxb)%sf(j,k,l) = rhs_vf(momxb)%sf(j,k,l) - &
                                 (0.5_wp* (rho_R * (vel_R(1,0))**2.0 + &
                                 pres_R)*(1._wp/dx(j)) + &
                                 0.5_wp*cfl * (rho_R*vel_R(1,0))*(1._wp/dx(j)))

                            !$acc atomic
                            rhs_vf(momxb+1)%sf(j, k, l) =  rhs_vf(momxb+1)%sf(j, k, l) - &
                                (0.5_wp * rho_R * vel_R(1,0)*vel_R(2,0)*(1._wp/dx(j)) + &
                                0.5_wp*cfl * (rho_R*vel_R(2,0))*(1._wp/dx(j)))

                            !$acc atomic
                            rhs_vf(momxb+2)%sf(j, k, l) =  rhs_vf(momxb+2)%sf(j, k, l) - &
                                (0.5_wp * rho_R * vel_R(1,0)*vel_R(3,0)*(1._wp/dx(j)) + &
                                0.5_wp*cfl * (rho_R*vel_R(3,0))*(1._wp/dx(j)))

                            !$acc atomic
                             rhs_vf(E_idx)%sf(j, k, l) = rhs_vf(E_idx)%sf(j, k, l) -  &
                                (0.5_wp * (vel_R(1,0) * (E_R + &
                                pres_R) )*(1._wp/dx(j)) + &
                                0.5_wp*cfl * (E_R)*(1._wp/dx(j)))
                        end do
                    end do
                end do
            end if
        else if (idir == 2) then
            if(p == 0) then
                !$acc parallel loop collapse(3) gang vector default(present) private(rho_L,gamma_L,pi_inf_L,mu_L,a_L,vel_L,vel_R, pres_L,alpha_L,alpha_R,alpha_rho_L,cfl,dvel1,dvel2,F_L,E_L,mu_R,rho_sf,alpha_rho_R)
                do l = 0, p
                    do k = -buff_size+5, n+buff_size-5
                        do j = 0, m

                            !$acc loop seq 
                            do i = 1, num_fluids
                                alpha_rho_L(i) = (1._wp/60._wp) * (-3._wp * q_prim_vf(i)%sf(j, k-1, l) + &
                                                27._wp * q_prim_vf(i)%sf(j, k, l) + &
                                                47._wp * q_prim_vf(i)%sf(j, k+1, l) -   &
                                                13._wp * q_prim_vf(i)%sf(j, k+2, l) + &
                                                2._wp * q_prim_vf(i)%sf(j, k+3, l))

                                alpha_rho_R(i) = (1._wp/60._wp) * (-3._wp * q_prim_vf(i)%sf(j, k+2, l) + &
                                                27._wp * q_prim_vf(i)%sf(j, k+1, l) + &
                                                47._wp * q_prim_vf(i)%sf(j, k, l) -   &
                                                13._wp * q_prim_vf(i)%sf(j, k-1, l) + &
                                                2._wp * q_prim_vf(i)%sf(j, k-2, l))

                                alpha_L(i) =  (1._wp/60._wp) * (-3._wp * q_prim_vf(E_idx+i)%sf(j, k-1, l) + &
                                                27._wp * q_prim_vf(E_idx+i)%sf(j, k, l) + &
                                                47._wp * q_prim_vf(E_idx+i)%sf(j, k+1, l) -   &
                                                13._wp * q_prim_vf(E_idx+i)%sf(j, k+2, l) + &
                                                2._wp * q_prim_vf(E_idx+i)%sf(j, k+3, l))
                            end do 

                            if(num_fluids > 1) then  
                                !$acc loop seq 
                                do i = 1,num_fluids                           
                                    alpha_R(i) = (1._wp/60._wp) * (-3._wp * q_prim_vf(E_idx+i)%sf(j, k+2, l) + &
                                                        27._wp * q_prim_vf(E_idx+i)%sf(j, k+1, l) + &
                                                        47._wp * q_prim_vf(E_idx+i)%sf(j, k, l) -   &
                                                        13._wp * q_prim_vf(E_idx+i)%sf(j, k-1, l) + &
                                                        2._wp * q_prim_vf(E_idx+i)%sf(j, k-2, l))
                                end do
                            else 
                                alpha_R(1) = 1._wp 
                            end if

                            F_L = (1._wp/60._wp) * (-3._wp * jac(j, k-1, l) + &
                                                27._wp * jac(j, k, l) + &
                                                47._wp * jac(j, k+1, l) -   &
                                                13._wp * jac(j, k+2, l) + &
                                                2._wp * jac(j, k+3, l))

                            !$acc loop seq 
                            do i = -2, 2
                                rho_L = 0._wp
                                !$acc loop seq 
                                do q = 1, num_fluids
                                    rho_L = rho_L + q_prim_vf(q)%sf(j+i,k,l)
                                end do
                                rho_sf(i,1) = rho_L
                            end do

                            !$acc loop seq 
                            do i = -2, 2
                                rho_L = 0._wp
                                !$acc loop seq 
                                do q = 1, num_fluids
                                    rho_L = rho_L + q_prim_vf(q)%sf(j,k+i,l)
                                end do
                                rho_sf(i,2) = rho_L
                            end do

                            rho_L = sum(alpha_rho_L)
                            gamma_L = sum(alpha_L*gammas)
                            pi_inf_L = sum(alpha_L*pi_infs)

                            rho_R = sum(alpha_rho_R)
                            gamma_R = sum(alpha_R*gammas)
                            pi_inf_R = sum(alpha_R*pi_infs)

                            !$acc loop seq 
                            do q = -3, 2
                                !$acc loop seq 
                                do i = 1, num_dims
                                    vel_L(i,q) = (1._wp/60._wp) * (-3._wp * q_prim_vf(momxb+i-1)%sf(j,k-1+q,l) + &
                                                    27._wp * q_prim_vf(momxb+i-1)%sf(j,k+q,l) + &
                                                    47._wp * q_prim_vf(momxb+i-1)%sf(j,k+1+q,l) -   &
                                                    13._wp * q_prim_vf(momxb+i-1)%sf(j,k+2+q,l) + &
                                                    2._wp * q_prim_vf(momxb+i-1)%sf(j,k+3+q,l)) / rho_L
                                end do 

                                !$acc loop seq 
                                do i = 1, num_dims
                                    vel_R(i,q) = (1._wp/60._wp) * (-3._wp * q_prim_vf(momxb+i-1)%sf(j,k+2+q,l) + &
                                                    27._wp * q_prim_vf(momxb+i-1)%sf(j,k+1+q,l) + &
                                                    47._wp * q_prim_vf(momxb+i-1)%sf(j,k+q,l) -   &
                                                    13._wp * q_prim_vf(momxb+i-1)%sf(j,k-1+q,l) + &
                                                    2._wp * q_prim_vf(momxb+i-1)%sf(j,k-2+q,l)) / rho_R
                                end do
                            end do

                            E_L = (1._wp/60._wp) * (-3._wp * q_prim_vf(E_idx)%sf(j,k-1,l) + &
                                                27._wp * q_prim_vf(E_idx)%sf(j,k,l) + &
                                                47._wp * q_prim_vf(E_idx)%sf(j,k+1,l) -   &
                                                13._wp * q_prim_vf(E_idx)%sf(j,k+2,l) + &
                                                2._wp * q_prim_vf(E_idx)%sf(j,k+3,l))

                            E_R = (1._wp/60._wp) * (-3._wp * q_prim_vf(E_idx)%sf(j,k+2,l) + &
                                                27._wp * q_prim_vf(E_idx)%sf(j,k+1,l) + &
                                                47._wp * q_prim_vf(E_idx)%sf(j,k,l) -   &
                                                13._wp * q_prim_vf(E_idx)%sf(j,k-1,l) + &
                                                2._wp * q_prim_vf(E_idx)%sf(j,k-2,l))       

                            pres_L = (E_L - pi_inf_L - 0.5_wp*rho_L*(vel_L(1,0)**2._wp + vel_L(2,0)**2._wp ))/gamma_L                    

                            pres_R = (E_R - pi_inf_R - 0.5_wp*rho_R*(vel_R(1,0)**2._wp + vel_R(2,0)**2._wp ))/gamma_R 

                            a_L = sqrt((pres_L*(1._wp/gamma_L + 1._wp) + pi_inf_L / gamma_L) / rho_L)
                            a_R = sqrt((pres_R*(1._wp/gamma_R + 1._wp) + pi_inf_R / gamma_R) / rho_R)

                            cfl = (max(sqrt(vel_L(1,0)**2._wp + vel_L(2,0)**2._wp),sqrt(vel_R(1,0)**2._wp + vel_R(2,0)**2._wp)) + max(a_L,a_R))


                            !$acc loop seq
                            do i = 1, num_fluids
                                !$acc atomic
                                rhs_vf(i)%sf(j,k+1,l) = rhs_vf(i)%sf(j,k+1,l) + &
                                    (0.5_wp * (alpha_rho_L(i) * &
                                    vel_L(2,0))*(1._wp/dy(k+1)) - &
                                    0.5_wp*cfl * (alpha_rho_L(i))*(1._wp/dy(k+1)))

                                !$acc atomic
                                rhs_vf(i)%sf(j,k,l) = rhs_vf(i)%sf(j,k,l) - &
                                    (0.5_wp * (alpha_rho_L(i) * &
                                    vel_L(2,0))*(1._wp/dy(k)) - &
                                    0.5_wp*cfl * (alpha_rho_L(i))*(1._wp/dy(k)) )
                            end do

                            if(num_fluids > 1) then 
                                !$acc loop seq
                                do i = 1, num_fluids
                                    !$acc atomic
                                    rhs_vf(advxb+i-1)%sf(j,k+1,l) = rhs_vf(advxb+i-1)%sf(j,k+1,l) + &
                                        (0.5_wp * (alpha_L(i) * &
                                        vel_L(2,0))*(1._wp/dy(k+1)) - &
                                        0.5_wp*cfl*(alpha_L(i))*(1._wp/dy(k+1)))

                                    !$acc atomic
                                    rhs_vf(advxb+i-1)%sf(j,k+1,l) = rhs_vf(advxb+i-1)%sf(j,k+1,l) &
                                    - (0.5_wp * q_prim_vf(advxb+i-1)%sf(j,k+1,l)  * vel_L(2,0)*(1._wp/dy(k+1)))

                                    !$acc atomic
                                    rhs_vf(advxb+i-1)%sf(j,k,l) = rhs_vf(advxb+i-1)%sf(j,k,l) - &
                                        (0.5_wp * (alpha_L(i) * &
                                        vel_L(2,0))*(1._wp/dy(k)) - &
                                        0.5_wp*cfl*(alpha_L(i))*(1._wp/dy(k)))

                                    !$acc atomic
                                    rhs_vf(advxb+i-1)%sf(j,k,l) = rhs_vf(advxb+i-1)%sf(j,k,l) &
                                    + (0.5_wp * q_prim_vf(advxb+i-1)%sf(j,k,l)  * vel_L(2,0)*(1._wp/dy(k)))
                                end do
                            end if

                            !$acc atomic
                            rhs_vf(momxb+1)%sf(j,k+1,l) = rhs_vf(momxb+1)%sf(j,k+1,l) + &
                                 (0.5_wp* (rho_L * (vel_L(2,0))**2.0 + &
                                 pres_L + F_L)*(1._wp/dy(k+1)) - &
                                 0.5_wp*cfl * (rho_L*vel_L(2,0))*(1._wp/dy(k+1)))

                            !$acc atomic
                            rhs_vf(momxb)%sf(j, k+1, l) =  rhs_vf(momxb)%sf(j, k+1, l) + &
                                (0.5_wp * rho_L * vel_L(1,0)*vel_L(2,0)*(1._wp/dy(k+1)) - &
                                0.5_wp*cfl * (rho_L*vel_L(1,0))*(1._wp/dy(k+1)))

                            !$acc atomic
                            rhs_vf(E_idx)%sf(j, k+1, l) = rhs_vf(E_idx)%sf(j, k+1, l) +  &
                                (0.5_wp * (vel_L(2,0) * (E_L + &
                                pres_L + F_L) )*(1._wp/dy(k+1)) - &
                                0.5_wp*cfl * (E_L)*(1._wp/dy(k+1)))

                            !$acc atomic
                            rhs_vf(momxb+1)%sf(j,k,l) = rhs_vf(momxb+1)%sf(j,k,l) - &
                                 (0.5_wp* (rho_L * (vel_L(2,0))**2.0 + &
                                 pres_L + F_L)*(1._wp/dy(k)) - &
                                 0.5_wp*cfl * (rho_L*vel_L(2,0))*(1._wp/dy(k)))

                            !$acc atomic
                            rhs_vf(momxb)%sf(j, k, l) =  rhs_vf(momxb)%sf(j, k, l) - &
                                (0.5_wp * rho_L * vel_L(1,0)*vel_L(2,0)*(1._wp/dy(k)) - &
                                0.5_wp*cfl * (rho_L*vel_L(1,0))*(1._wp/dy(k)))

                            !$acc atomic
                            rhs_vf(E_idx)%sf(j, k, l) = rhs_vf(E_idx)%sf(j, k, l) -  &
                                (0.5_wp * (vel_L(2,0) * (E_L + &
                                pres_L + F_L) )*(1._wp/dy(k)) - &
                                0.5_wp*cfl * (E_L)*(1._wp/dy(k)))

                            !! duy & dvx
                            dvel1 = (1/(12._wp*dy(k))) * ( &
                            8._wp*q_prim_vf(momxb)%sf(j,k+1,l)/rho_sf(1,2) - &
                            8._wp*q_prim_vf(momxb)%sf(j,k-1,l)/rho_sf(-1,2) + &
                            q_prim_vf(momxb)%sf(j,k-2,l)/rho_sf(-2,2) - &
                            q_prim_vf(momxb)%sf(j,k+2,l)/rho_sf(2,2) )


                            dvel2 = (1/(12._wp*dx(j))) * ( &
                            8._wp*q_prim_vf(momxb+1)%sf(j+1,k,l)/rho_sf(1,1) - &
                            8._wp*q_prim_vf(momxb+1)%sf(j-1,k,l)/rho_sf(-1,1) + &
                            q_prim_vf(momxb+1)%sf(j-2,k,l)/rho_sf(-2,1) - &
                            q_prim_vf(momxb+1)%sf(j+2,k,l)/rho_sf(2,1) )

                            if(num_fluids > 1) then 
                                !$acc loop seq
                                do q = -3, 2 
                                    mu_L(q) = 0._wp
                                    mu_R(q) = 0._wp
                                    !$acc loop seq 
                                    do i = 1, num_fluids
                                        mu_L(q) = (1._wp/60._wp) * (-3._wp * q_prim_vf(E_idx+i)%sf(j, k-1+q, l) + &
                                                27._wp * q_prim_vf(E_idx+i)%sf(j, k+q, l) + &
                                                47._wp * q_prim_vf(E_idx+i)%sf(j, k+1+q, l) -   &
                                                13._wp * q_prim_vf(E_idx+i)%sf(j, k+2+q, l) + &
                                                2._wp * q_prim_vf(E_idx+i)%sf(j, k+3+q, l)) / Res(1, i) + mu_L(q)
                                        
                                        mu_R(q) = (1._wp/60._wp) * (-3._wp * q_prim_vf(E_idx+i)%sf(j, k+2+q, l) + &
                                                27._wp * q_prim_vf(E_idx+i)%sf(j, k+1+q, l) + &
                                                47._wp * q_prim_vf(E_idx+i)%sf(j, k+q, l) -   &
                                                13._wp * q_prim_vf(E_idx+i)%sf(j, k-1+q, l) + &
                                                2._wp * q_prim_vf(E_idx+i)%sf(j, k-2+q, l)) / Res(1, i) + mu_R(q)
                                    end do
                                end do
                            else
                                !$acc loop seq
                                do q = -3, 2
                                    mu_L(q) = 1._wp / Res(1, 1) 
                                    mu_R(q) = 1._wp / Res(1, 1) 
                                end do
                            end if

                            if(viscous) then

                                !$acc atomic
                                rhs_vf(momxb)%sf(j,k+2,l) = rhs_vf(momxb)%sf(j,k+2,l) - 0.5_wp*mu_L(1)*((-3._wp/60._wp)*(dvel1+dvel2))*(1._wp/dy(k+2))
                                !$acc atomic
                                rhs_vf(momxb)%sf(j,k+1,l) = rhs_vf(momxb)%sf(j,k+1,l) - 0.5_wp*mu_L(0)*((27._wp/60._wp)*(dvel1+dvel2))*(1._wp/dy(k+1))
                                !$acc atomic
                                rhs_vf(momxb)%sf(j,k,l) = rhs_vf(momxb)%sf(j,k,l) - 0.5_wp*mu_L(-1)*((47._wp/60._wp)*(dvel1+dvel2))*(1._wp/dy(k))
                                !$acc atomic
                                rhs_vf(momxb)%sf(j,k-1,l) = rhs_vf(momxb)%sf(j,k-1,l) - 0.5_wp*mu_L(-2)*((-13._wp/60._wp)*(dvel1+dvel2))*(1._wp/dy(k-1))
                                !$acc atomic
                                rhs_vf(momxb)%sf(j,k-2,l) = rhs_vf(momxb)%sf(j,k-2,l) - 0.5_wp*mu_L(-3)*((2._wp/60._wp)*(dvel1+dvel2))*(1._wp/dy(k-2))

                                !$acc atomic
                                rhs_vf(E_idx)%sf(j,k+2,l) = rhs_vf(E_idx)%sf(j,k+2,l) - 0.5_wp*mu_L(1)*((-3._wp/60._wp)*(dvel1+dvel2))*vel_L(1,1)*(1._wp/dy(k+2))
                                !$acc atomic
                                rhs_vf(E_idx)%sf(j,k+1,l) = rhs_vf(E_idx)%sf(j,k+1,l) - 0.5_wp*mu_L(0)*((27._wp/60._wp)*(dvel1+dvel2))*vel_L(1,0)*(1._wp/dy(k+1))
                                !$acc atomic
                                rhs_vf(E_idx)%sf(j,k,l) = rhs_vf(E_idx)%sf(j,k,l) - 0.5_wp*mu_L(-1)*((47._wp/60._wp)*(dvel1+dvel2))*vel_L(1,-1)*(1._wp/dy(k))
                                !$acc atomic
                                rhs_vf(E_idx)%sf(j,k-1,l) = rhs_vf(E_idx)%sf(j,k-1,l) - 0.5_wp*mu_L(-2)*((-13._wp/60._wp)*(dvel1+dvel2))*vel_L(1,-2)*(1._wp/dy(k-1))
                                !$acc atomic
                                rhs_vf(E_idx)%sf(j,k-2,l) = rhs_vf(E_idx)%sf(j,k-2,l) - 0.5_wp*mu_L(-3)*((2._wp/60._wp)*(dvel1+dvel2))*vel_L(1,-3)*(1._wp/dy(k-2))

                                !$acc atomic
                                rhs_vf(momxb)%sf(j,k+1,l) = rhs_vf(momxb)%sf(j,k+1,l) + 0.5_wp*mu_L(1)*((-3._wp/60._wp)*(dvel1+dvel2))*(1._wp/dy(k+1))
                                !$acc atomic
                                rhs_vf(momxb)%sf(j,k,l) = rhs_vf(momxb)%sf(j,k,l) + 0.5_wp*mu_L(0)*((27._wp/60._wp)*(dvel1+dvel2))*(1._wp/dy(k))
                                !$acc atomic
                                rhs_vf(momxb)%sf(j,k-1,l) = rhs_vf(momxb)%sf(j,k-1,l) + 0.5_wp*mu_L(-1)*((47._wp/60._wp)*(dvel1+dvel2))*(1._wp/dy(k-1))
                                !$acc atomic
                                rhs_vf(momxb)%sf(j,k-2,l) = rhs_vf(momxb)%sf(j,k-2,l) + 0.5_wp*mu_L(-2)*((-13._wp/60._wp)*(dvel1+dvel2))*(1._wp/dy(k-2))
                                !$acc atomic
                                rhs_vf(momxb)%sf(j,k-3,l) = rhs_vf(momxb)%sf(j,k-3,l) + 0.5_wp*mu_L(-3)*((2._wp/60._wp)*(dvel1+dvel2))*(1._wp/dy(k-3))

                                !$acc atomic
                                rhs_vf(E_idx)%sf(j,k+1,l) = rhs_vf(E_idx)%sf(j,k+1,l) + 0.5_wp*mu_L(1)*((-3._wp/60._wp)*(dvel1+dvel2))*vel_L(1,1)*(1._wp/dy(k+1))
                                !$acc atomic
                                rhs_vf(E_idx)%sf(j,k,l) = rhs_vf(E_idx)%sf(j,k,l) + 0.5_wp*mu_L(0)*((27._wp/60._wp)*(dvel1+dvel2))*vel_L(1,0)*(1._wp/dy(k))
                                !$acc atomic
                                rhs_vf(E_idx)%sf(j,k-1,l) = rhs_vf(E_idx)%sf(j,k-1,l) + 0.5_wp*mu_L(-1)*((47._wp/60._wp)*(dvel1+dvel2))*vel_L(1,-1)*(1._wp/dy(k-1))
                                !$acc atomic
                                rhs_vf(E_idx)%sf(j,k-2,l) = rhs_vf(E_idx)%sf(j,k-2,l) + 0.5_wp*mu_L(-2)*((-13._wp/60._wp)*(dvel1+dvel2))*vel_L(1,-2)*(1._wp/dy(k-2))
                                !$acc atomic
                                rhs_vf(E_idx)%sf(j,k-3,l) = rhs_vf(E_idx)%sf(j,k-3,l) + 0.5_wp*mu_L(-3)*((2._wp/60._wp)*(dvel1+dvel2))*vel_L(1,-3)*(1._wp/dy(k-3))
                               
                                !$acc atomic
                                rhs_vf(momxb)%sf(j,k-1,l) = rhs_vf(momxb)%sf(j,k-1,l) - 0.5_wp*mu_R(-2)*((-3._wp/60._wp)*(dvel1+dvel2))*(1._wp/dy(k-1))
                                !$acc atomic
                                rhs_vf(momxb)%sf(j,k,l) = rhs_vf(momxb)%sf(j,k,l) - 0.5_wp*mu_R(-1)*((27._wp/60._wp)*(dvel1+dvel2))*(1._wp/dy(k))
                                !$acc atomic
                                rhs_vf(momxb)%sf(j,k+1,l) = rhs_vf(momxb)%sf(j,k+1,l) - 0.5_wp*mu_R(0)*((47._wp/60._wp)*(dvel1+dvel2))*(1._wp/dy(k+1))
                                !$acc atomic
                                rhs_vf(momxb)%sf(j,k+2,l) = rhs_vf(momxb)%sf(j,k+2,l) - 0.5_wp*mu_R(1)*((-13._wp/60._wp)*(dvel1+dvel2))*(1._wp/dy(k+2))
                                !$acc atomic
                                rhs_vf(momxb)%sf(j,k+3,l) = rhs_vf(momxb)%sf(j,k+3,l) - 0.5_wp*mu_R(2)*((2._wp/60._wp)*(dvel1+dvel2))*(1._wp/dy(k+3))

                                !$acc atomic
                                rhs_vf(E_idx)%sf(j,k-1,l) = rhs_vf(E_idx)%sf(j,k-1,l) - 0.5_wp*mu_R(-2)*((-3._wp/60._wp)*(dvel1+dvel2))*vel_R(1,-2)*(1._wp/dy(k-1))
                                !$acc atomic
                                rhs_vf(E_idx)%sf(j,k,l) = rhs_vf(E_idx)%sf(j,k,l) - 0.5_wp*mu_R(-1)*((27._wp/60._wp)*(dvel1+dvel2))*vel_R(1,-1)*(1._wp/dy(k))
                                !$acc atomic
                                rhs_vf(E_idx)%sf(j,k+1,l) = rhs_vf(E_idx)%sf(j,k+1,l) - 0.5_wp*mu_R(0)*((47._wp/60._wp)*(dvel1+dvel2))*vel_R(1,0)*(1._wp/dy(k+1))
                                !$acc atomic
                                rhs_vf(E_idx)%sf(j,k+2,l) = rhs_vf(E_idx)%sf(j,k+2,l) - 0.5_wp*mu_R(1)*((-13._wp/60._wp)*(dvel1+dvel2))*vel_R(1,1)*(1._wp/dy(k+2))
                                !$acc atomic
                                rhs_vf(E_idx)%sf(j,k+3,l) = rhs_vf(E_idx)%sf(j,k+3,l) - 0.5_wp*mu_R(2)*((2._wp/60._wp)*(dvel1+dvel2))*vel_R(1,2)*(1._wp/dy(k+3)) 

                                !$acc atomic
                                rhs_vf(momxb)%sf(j,k-2,l) = rhs_vf(momxb)%sf(j,k-2,l) + 0.5_wp*mu_R(-2)*((-3._wp/60._wp)*(dvel1+dvel2))*(1._wp/dy(k-2))
                                !$acc atomic
                                rhs_vf(momxb)%sf(j,k-1,l) = rhs_vf(momxb)%sf(j,k-1,l) + 0.5_wp*mu_R(-1)*((27._wp/60._wp)*(dvel1+dvel2))*(1._wp/dy(k-1))
                                !$acc atomic
                                rhs_vf(momxb)%sf(j,k,l) = rhs_vf(momxb)%sf(j,k,l) + 0.5_wp*mu_R(0)*((47._wp/60._wp)*(dvel1+dvel2))*(1._wp/dy(k))
                                !$acc atomic
                                rhs_vf(momxb)%sf(j,k+1,l) = rhs_vf(momxb)%sf(j,k+1,l) + 0.5_wp*mu_R(1)*((-13._wp/60._wp)*(dvel1+dvel2))*(1._wp/dy(k+1))
                                !$acc atomic
                                rhs_vf(momxb)%sf(j,k+2,l) = rhs_vf(momxb)%sf(j,k+2,l) + 0.5_wp*mu_R(2)*((2._wp/60._wp)*(dvel1+dvel2))*(1._wp/dy(k+2))

                                !$acc atomic
                                rhs_vf(E_idx)%sf(j,k-2,l) = rhs_vf(E_idx)%sf(j,k-2,l) + 0.5_wp*mu_R(-2)*((-3._wp/60._wp)*(dvel1+dvel2))*vel_R(1,-2)*(1._wp/dy(k-2))
                                !$acc atomic
                                rhs_vf(E_idx)%sf(j,k-1,l) = rhs_vf(E_idx)%sf(j,k-1,l) + 0.5_wp*mu_R(-1)*((27._wp/60._wp)*(dvel1+dvel2))*vel_R(1,-1)*(1._wp/dy(k-1))
                                !$acc atomic
                                rhs_vf(E_idx)%sf(j,k,l) = rhs_vf(E_idx)%sf(j,k,l) + 0.5_wp*mu_R(0)*((47._wp/60._wp)*(dvel1+dvel2))*vel_R(1,0)*(1._wp/dy(k))
                                !$acc atomic
                                rhs_vf(E_idx)%sf(j,k+1,l) = rhs_vf(E_idx)%sf(j,k+1,l) + 0.5_wp*mu_R(1)*((-13._wp/60._wp)*(dvel1+dvel2))*vel_R(1,1)*(1._wp/dy(k+1))
                                !$acc atomic
                                rhs_vf(E_idx)%sf(j,k+2,l) = rhs_vf(E_idx)%sf(j,k+2,l) + 0.5_wp*mu_R(2)*((2._wp/60._wp)*(dvel1+dvel2))*vel_R(1,2)*(1._wp/dy(k+2))

                            end if

                            !! dvy & dux
                            dvel1 = (1/(12._wp*dy(k))) * ( &
                            8._wp*q_prim_vf(momxb+1)%sf(j,k+1,l)/rho_sf(1,2) - &
                            8._wp*q_prim_vf(momxb+1)%sf(j,k-1,l)/rho_sf(-1,2) + &
                            q_prim_vf(momxb+1)%sf(j,k-2,l)/rho_sf(-2,2) - &
                            q_prim_vf(momxb+1)%sf(j,k+2,l)/rho_sf(2,2) )

                            dvel2 = (1/(12._wp*dx(j))) * ( &
                            8._wp*q_prim_vf(momxb)%sf(j+1,k,l)/rho_sf(1,1) - &
                            8._wp*q_prim_vf(momxb)%sf(j-1,k,l)/rho_sf(-1,1) + &
                            q_prim_vf(momxb)%sf(j-2,k,l)/rho_sf(-2,1) - &
                            q_prim_vf(momxb)%sf(j+2,k,l)/rho_sf(2,1) )

                            if(viscous) then

                                !$acc atomic
                                rhs_vf(momxb+1)%sf(j,k+2,l) = rhs_vf(momxb+1)%sf(j,k+2,l) - 0.5_wp*mu_L(1)*((-3._wp/60._wp)*(4._wp*dvel1-2._wp*dvel2)/3._wp)*(1._wp/dy(k+2))
                                !$acc atomic
                                rhs_vf(momxb+1)%sf(j,k+1,l) = rhs_vf(momxb+1)%sf(j,k+1,l) - 0.5_wp*mu_L(0)*((27._wp/60._wp)*(4._wp*dvel1-2._wp*dvel2)/3._wp)*(1._wp/dy(k+1))
                                !$acc atomic
                                rhs_vf(momxb+1)%sf(j,k,l) = rhs_vf(momxb+1)%sf(j,k,l) - 0.5_wp*mu_L(-1)*((47._wp/60._wp)*(4._wp*dvel1-2._wp*dvel2)/3._wp)*(1._wp/dy(k))
                                !$acc atomic
                                rhs_vf(momxb+1)%sf(j,k-1,l) = rhs_vf(momxb+1)%sf(j,k-1,l) - 0.5_wp*mu_L(-2)*((-13._wp/60._wp)*(4._wp*dvel1-2._wp*dvel2)/3._wp)*(1._wp/dy(k-1))
                                !$acc atomic
                                rhs_vf(momxb+1)%sf(j,k-2,l) = rhs_vf(momxb+1)%sf(j,k-2,l) - 0.5_wp*mu_L(-3)*((2._wp/60._wp)*(4._wp*dvel1-2._wp*dvel2)/3._wp)*(1._wp/dy(k-2))

                                !$acc atomic
                                rhs_vf(E_idx)%sf(j,k+2,l) = rhs_vf(E_idx)%sf(j,k+2,l) - 0.5_wp*mu_L(1)*((-3._wp/60._wp)*(4._wp*dvel1-2._wp*dvel2)/3._wp)*vel_L(2,1)*(1._wp/dy(k+2))
                                !$acc atomic
                                rhs_vf(E_idx)%sf(j,k+1,l) = rhs_vf(E_idx)%sf(j,k+1,l) - 0.5_wp*mu_L(0)*((27._wp/60._wp)*(4._wp*dvel1-2._wp*dvel2)/3._wp)*vel_L(2,0)*(1._wp/dy(k+1))
                                !$acc atomic
                                rhs_vf(E_idx)%sf(j,k,l) = rhs_vf(E_idx)%sf(j,k,l) - 0.5_wp*mu_L(-1)*((47._wp/60._wp)*(4._wp*dvel1-2._wp*dvel2)/3._wp)*vel_L(2,-1)*(1._wp/dy(k))
                                !$acc atomic
                                rhs_vf(E_idx)%sf(j,k-1,l) = rhs_vf(E_idx)%sf(j,k-1,l) - 0.5_wp*mu_L(-2)*((-13._wp/60._wp)*(4._wp*dvel1-2._wp*dvel2)/3._wp)*vel_L(2,-2)*(1._wp/dy(k-1))
                                !$acc atomic
                                rhs_vf(E_idx)%sf(j,k-2,l) = rhs_vf(E_idx)%sf(j,k-2,l) - 0.5_wp*mu_L(-3)*((2._wp/60._wp)*(4._wp*dvel1-2._wp*dvel2)/3._wp)*vel_L(2,-3)*(1._wp/dy(k-2))

                                !$acc atomic
                                rhs_vf(momxb+1)%sf(j,k+1,l) = rhs_vf(momxb+1)%sf(j,k+1,l) + 0.5_wp*mu_L(1)*((-3._wp/60._wp)*(4._wp*dvel1-2._wp*dvel2)/3._wp)*(1._wp/dy(k+1))
                                !$acc atomic
                                rhs_vf(momxb+1)%sf(j,k,l) = rhs_vf(momxb+1)%sf(j,k,l) + 0.5_wp*mu_L(0)*((27._wp/60._wp)*(4._wp*dvel1-2._wp*dvel2)/3._wp)*(1._wp/dy(k))
                                !$acc atomic
                                rhs_vf(momxb+1)%sf(j,k-1,l) = rhs_vf(momxb+1)%sf(j,k-1,l) + 0.5_wp*mu_L(-1)*((47._wp/60._wp)*(4._wp*dvel1-2._wp*dvel2)/3._wp)*(1._wp/dy(k-1))
                                !$acc atomic
                                rhs_vf(momxb+1)%sf(j,k-2,l) = rhs_vf(momxb+1)%sf(j,k-2,l) + 0.5_wp*mu_L(-2)*((-13._wp/60._wp)*(4._wp*dvel1-2._wp*dvel2)/3._wp)*(1._wp/dy(k-2))
                                !$acc atomic
                                rhs_vf(momxb+1)%sf(j,k-3,l) = rhs_vf(momxb+1)%sf(j,k-3,l) + 0.5_wp*mu_L(-3)*((2._wp/60._wp)*(4._wp*dvel1-2._wp*dvel2)/3._wp)*(1._wp/dy(k-3))

                                !$acc atomic
                                rhs_vf(E_idx)%sf(j,k+1,l) = rhs_vf(E_idx)%sf(j,k+1,l) + 0.5_wp*mu_L(1)*((-3._wp/60._wp)*(4._wp*dvel1-2._wp*dvel2)/3._wp)*vel_L(2,1)*(1._wp/dy(k+1))
                                !$acc atomic
                                rhs_vf(E_idx)%sf(j,k,l) = rhs_vf(E_idx)%sf(j,k,l) + 0.5_wp*mu_L(0)*((27._wp/60._wp)*(4._wp*dvel1-2._wp*dvel2)/3._wp)*vel_L(2,0)*(1._wp/dy(k))
                                !$acc atomic
                                rhs_vf(E_idx)%sf(j,k-1,l) = rhs_vf(E_idx)%sf(j,k-1,l) + 0.5_wp*mu_L(-1)*((47._wp/60._wp)*(4._wp*dvel1-2._wp*dvel2)/3._wp)*vel_L(2,-1)*(1._wp/dy(k-1))
                                !$acc atomic
                                rhs_vf(E_idx)%sf(j,k-2,l) = rhs_vf(E_idx)%sf(j,k-2,l) + 0.5_wp*mu_L(-2)*((-13._wp/60._wp)*(4._wp*dvel1-2._wp*dvel2)/3._wp)*vel_L(2,-2)*(1._wp/dy(k-2))
                                !$acc atomic
                                rhs_vf(E_idx)%sf(j,k-3,l) = rhs_vf(E_idx)%sf(j,k-3,l) + 0.5_wp*mu_L(-3)*((2._wp/60._wp)*(4._wp*dvel1-2._wp*dvel2)/3._wp)*vel_L(2,-3)*(1._wp/dy(k-3))

                                !$acc atomic
                                rhs_vf(momxb+1)%sf(j,k-1,l) = rhs_vf(momxb+1)%sf(j,k-1,l) - 0.5_wp*mu_R(-2)*((-3._wp/60._wp)*(4._wp*dvel1-2._wp*dvel2)/3._wp)*(1._wp/dy(k-1))
                                !$acc atomic
                                rhs_vf(momxb+1)%sf(j,k,l) = rhs_vf(momxb+1)%sf(j,k,l) - 0.5_wp*mu_R(-1)*((27._wp/60._wp)*(4._wp*dvel1-2._wp*dvel2)/3._wp)*(1._wp/dy(k))
                                !$acc atomic
                                rhs_vf(momxb+1)%sf(j,k+1,l) = rhs_vf(momxb+1)%sf(j,k+1,l) - 0.5_wp*mu_R(0)*((47._wp/60._wp)*(4._wp*dvel1-2._wp*dvel2)/3._wp)*(1._wp/dy(k+1))
                                !$acc atomic
                                rhs_vf(momxb+1)%sf(j,k+2,l) = rhs_vf(momxb+1)%sf(j,k+2,l) - 0.5_wp*mu_R(1)*((-13._wp/60._wp)*(4._wp*dvel1-2._wp*dvel2)/3._wp)*(1._wp/dy(k+2))
                                !$acc atomic
                                rhs_vf(momxb+1)%sf(j,k+3,l) = rhs_vf(momxb+1)%sf(j,k+3,l) - 0.5_wp*mu_R(2)*((2._wp/60._wp)*(4._wp*dvel1-2._wp*dvel2)/3._wp)*(1._wp/dy(k+3))

                                !$acc atomic
                                rhs_vf(E_idx)%sf(j,k-1,l) = rhs_vf(E_idx)%sf(j,k-1,l) - 0.5_wp*mu_R(-2)*((-3._wp/60._wp)*(4._wp*dvel1-2._wp*dvel2)/3._wp)*vel_R(2,-2)*(1._wp/dy(k-1))
                                !$acc atomic
                                rhs_vf(E_idx)%sf(j,k,l) = rhs_vf(E_idx)%sf(j,k,l) - 0.5_wp*mu_R(-1)*((27._wp/60._wp)*(4._wp*dvel1-2._wp*dvel2)/3._wp)*vel_R(2,-1)*(1._wp/dy(k))
                                !$acc atomic
                                rhs_vf(E_idx)%sf(j,k+1,l) = rhs_vf(E_idx)%sf(j,k+1,l) - 0.5_wp*mu_R(0)*((47._wp/60._wp)*(4._wp*dvel1-2._wp*dvel2)/3._wp)*vel_R(2,0)*(1._wp/dy(k+1))
                                !$acc atomic
                                rhs_vf(E_idx)%sf(j,k+2,l) = rhs_vf(E_idx)%sf(j,k+2,l) - 0.5_wp*mu_R(1)*((-13._wp/60._wp)*(4._wp*dvel1-2._wp*dvel2)/3._wp)*vel_R(2,1)*(1._wp/dy(k+2))
                                !$acc atomic
                                rhs_vf(E_idx)%sf(j,k+3,l) = rhs_vf(E_idx)%sf(j,k+3,l) - 0.5_wp*mu_R(2)*((2._wp/60._wp)*(4._wp*dvel1-2._wp*dvel2)/3._wp)*vel_R(2,2)*(1._wp/dy(k+3)) 

                                !$acc atomic
                                rhs_vf(momxb+1)%sf(j,k-2,l) = rhs_vf(momxb+1)%sf(j,k-2,l) + 0.5_wp*mu_R(-2)*((-3._wp/60._wp)*(4._wp*dvel1-2._wp*dvel2)/3._wp)*(1._wp/dy(k-2))
                                !$acc atomic
                                rhs_vf(momxb+1)%sf(j,k-1,l) = rhs_vf(momxb+1)%sf(j,k-1,l) + 0.5_wp*mu_R(-1)*((27._wp/60._wp)*(4._wp*dvel1-2._wp*dvel2)/3._wp)*(1._wp/dy(k-1))
                                !$acc atomic
                                rhs_vf(momxb+1)%sf(j,k,l) = rhs_vf(momxb+1)%sf(j,k,l) + 0.5_wp*mu_R(0)*((47._wp/60._wp)*(4._wp*dvel1-2._wp*dvel2)/3._wp)*(1._wp/dy(k))
                                !$acc atomic
                                rhs_vf(momxb+1)%sf(j,k+1,l) = rhs_vf(momxb+1)%sf(j,k+1,l) + 0.5_wp*mu_R(1)*((-13._wp/60._wp)*(4._wp*dvel1-2._wp*dvel2)/3._wp)*(1._wp/dy(k+1))
                                !$acc atomic
                                rhs_vf(momxb+1)%sf(j,k+2,l) = rhs_vf(momxb+1)%sf(j,k+2,l) + 0.5_wp*mu_R(2)*((2._wp/60._wp)*(4._wp*dvel1-2._wp*dvel2)/3._wp)*(1._wp/dy(k+2))

                                !$acc atomic
                                rhs_vf(E_idx)%sf(j,k-2,l) = rhs_vf(E_idx)%sf(j,k-2,l) + 0.5_wp*mu_R(-2)*((-3._wp/60._wp)*(4._wp*dvel1-2._wp*dvel2)/3._wp)*vel_R(2,-2)*(1._wp/dy(k-2))
                                !$acc atomic
                                rhs_vf(E_idx)%sf(j,k-1,l) = rhs_vf(E_idx)%sf(j,k-1,l) + 0.5_wp*mu_R(-1)*((27._wp/60._wp)*(4._wp*dvel1-2._wp*dvel2)/3._wp)*vel_R(2,-1)*(1._wp/dy(k-1))
                                !$acc atomic
                                rhs_vf(E_idx)%sf(j,k,l) = rhs_vf(E_idx)%sf(j,k,l) + 0.5_wp*mu_R(0)*((47._wp/60._wp)*(4._wp*dvel1-2._wp*dvel2)/3._wp)*vel_R(2,0)*(1._wp/dy(k))
                                !$acc atomic
                                rhs_vf(E_idx)%sf(j,k+1,l) = rhs_vf(E_idx)%sf(j,k+1,l) + 0.5_wp*mu_R(1)*((-13._wp/60._wp)*(4._wp*dvel1-2._wp*dvel2)/3._wp)*vel_R(2,1)*(1._wp/dy(k+1))
                                !$acc atomic
                                rhs_vf(E_idx)%sf(j,k+2,l) = rhs_vf(E_idx)%sf(j,k+2,l) + 0.5_wp*mu_R(2)*((2._wp/60._wp)*(4._wp*dvel1-2._wp*dvel2)/3._wp)*vel_R(2,2)*(1._wp/dy(k+2))

                            end if

                            F_L = (1._wp/60._wp) * (-3._wp * jac(j, k+2, l) + &
                                                27._wp * jac(j, k+1, l) + &
                                                47._wp * jac(j, k, l) -   &
                                                13._wp * jac(j, k-1, l) + &
                                                2._wp * jac(j, k-2, l))

                            !$acc loop seq
                            do i = 1, num_fluids
                                !$acc atomic
                                rhs_vf(i)%sf(j,k+1,l) = rhs_vf(i)%sf(j,k+1,l) + &
                                    (0.5_wp * (alpha_rho_R(i) * &
                                    vel_R(2,0))*(1._wp/dy(k+1)) + &
                                    0.5_wp*cfl * (alpha_rho_R(i))*(1._wp/dy(k+1)))

                                !$acc atomic
                                rhs_vf(i)%sf(j,k,l) = rhs_vf(i)%sf(j,k,l) - &
                                    (0.5_wp * (alpha_rho_R(i) * &
                                    vel_R(2,0))*(1._wp/dy(k)) + &
                                    0.5_wp*cfl * (alpha_rho_R(i))*(1._wp/dy(k))) 
                            end do

                            if(num_fluids > 1) then 
                                !$acc loop seq
                                do i = 1, num_fluids
                                    !$acc atomic
                                    rhs_vf(advxb+i-1)%sf(j,k+1,l) = rhs_vf(advxb+i-1)%sf(j,k+1,l) + &
                                        (0.5_wp * (alpha_R(i) * &
                                        vel_R(2,0))*(1._wp/dy(k+1)) + &
                                        0.5_wp*cfl*(alpha_R(i))*(1._wp/dy(k+1)))

                                    !$acc atomic
                                    rhs_vf(advxb+i-1)%sf(j,k+1,l) = rhs_vf(advxb+i-1)%sf(j,k+1,l) &
                                    - (0.5_wp * q_prim_vf(advxb+i-1)%sf(j,k+1,l)  * vel_R(2,0)*(1._wp/dy(k+1)))

                                    !$acc atomic
                                    rhs_vf(advxb+i-1)%sf(j,k,l) = rhs_vf(advxb+i-1)%sf(j,k,l) - &
                                        (0.5_wp * (alpha_R(i) * &
                                        vel_R(2,0))*(1._wp/dy(k)) + &
                                        0.5_wp*cfl*(alpha_R(i))*(1._wp/dy(k)))

                                    !$acc atomic
                                    rhs_vf(advxb+i-1)%sf(j,k,l) = rhs_vf(advxb+i-1)%sf(j,k,l) &
                                    + (0.5_wp * q_prim_vf(advxb+i-1)%sf(j,k,l)  * vel_R(2,0)*(1._wp/dy(k)))
                                end do
                            end if

                            !$acc atomic
                            rhs_vf(momxb+1)%sf(j,k+1,l) = rhs_vf(momxb+1)%sf(j,k+1,l) + &
                                 (0.5_wp* (rho_R * (vel_R(2,0))**2.0 + &
                                 pres_R+F_L)*(1._wp/dy(k+1)) + &
                                 0.5_wp*cfl * (rho_R*vel_R(2,0))*(1._wp/dy(k+1)))

                            !$acc atomic
                            rhs_vf(momxb)%sf(j, k+1, l) =  rhs_vf(momxb)%sf(j, k+1, l) + &
                                (0.5_wp * rho_R * vel_R(2,0)*vel_R(1,0)*(1._wp/dy(k+1)) + &
                                0.5_wp*cfl * (rho_R*vel_R(1,0))*(1._wp/dy(k+1)))

                            !$acc atomic
                             rhs_vf(E_idx)%sf(j, k+1, l) = rhs_vf(E_idx)%sf(j, k+1, l) +  &
                                (0.5_wp * (vel_R(2,0) * (E_R + &
                                pres_R+F_L ) )*(1._wp/dy(k+1)) + &
                                0.5_wp*cfl * (E_R)*(1._wp/dy(k+1)))

                            !$acc atomic
                            rhs_vf(momxb+1)%sf(j,k,l) = rhs_vf(momxb+1)%sf(j,k,l) - &
                                 (0.5_wp* (rho_R * (vel_R(2,0))**2.0 + &
                                 pres_R+F_L)*(1._wp/dy(k)) + &
                                 0.5_wp*cfl * (rho_R*vel_R(2,0))*(1._wp/dy(k)))

                            !$acc atomic
                            rhs_vf(momxb)%sf(j, k, l) =  rhs_vf(momxb)%sf(j, k, l) - &
                                (0.5_wp * rho_R * vel_R(2,0)*vel_R(1,0)*(1._wp/dy(k)) + &
                                0.5_wp*cfl * (rho_R*vel_R(1,0))*(1._wp/dy(k)))

                            !$acc atomic
                             rhs_vf(E_idx)%sf(j, k, l) = rhs_vf(E_idx)%sf(j, k, l) -  &
                                (0.5_wp * (vel_R(2,0) * (E_R + &
                                pres_R+F_L) )*(1._wp/dy(k)) + &
                                0.5_wp*cfl * (E_R)*(1._wp/dy(k)))

                        end do
                    end do
                end do
            else
                !$acc parallel loop collapse(3) gang vector default(present) private(rho_L,gamma_L,pi_inf_L,mu_L,a_L,vel_L,vel_R, pres_L,alpha_L,alpha_R,alpha_rho_L,cfl,dvel1,dvel2,F_L,E_L,mu_R,rho_sf,alpha_rho_R)
                do l = 0, p
                    do k = -buff_size+5, n+buff_size-5
                        do j = 0, m

                            !$acc loop seq 
                            do i = 1, num_fluids
                                alpha_rho_L(i) = (1._wp/60._wp) * (-3._wp * q_prim_vf(i)%sf(j, k-1, l) + &
                                                27._wp * q_prim_vf(i)%sf(j, k, l) + &
                                                47._wp * q_prim_vf(i)%sf(j, k+1, l) -   &
                                                13._wp * q_prim_vf(i)%sf(j, k+2, l) + &
                                                2._wp * q_prim_vf(i)%sf(j, k+3, l))

                                alpha_rho_R(i) = (1._wp/60._wp) * (-3._wp * q_prim_vf(i)%sf(j, k+2, l) + &
                                                27._wp * q_prim_vf(i)%sf(j, k+1, l) + &
                                                47._wp * q_prim_vf(i)%sf(j, k, l) -   &
                                                13._wp * q_prim_vf(i)%sf(j, k-1, l) + &
                                                2._wp * q_prim_vf(i)%sf(j, k-2, l))

                                alpha_L(i) =  (1._wp/60._wp) * (-3._wp * q_prim_vf(E_idx+i)%sf(j, k-1, l) + &
                                                27._wp * q_prim_vf(E_idx+i)%sf(j, k, l) + &
                                                47._wp * q_prim_vf(E_idx+i)%sf(j, k+1, l) -   &
                                                13._wp * q_prim_vf(E_idx+i)%sf(j, k+2, l) + &
                                                2._wp * q_prim_vf(E_idx+i)%sf(j, k+3, l))
                            end do 

                            if(num_fluids > 1) then  
                                !$acc loop seq 
                                do i = 1,num_fluids                           
                                    alpha_R(i) = (1._wp/60._wp) * (-3._wp * q_prim_vf(E_idx+i)%sf(j, k+2, l) + &
                                                        27._wp * q_prim_vf(E_idx+i)%sf(j, k+1, l) + &
                                                        47._wp * q_prim_vf(E_idx+i)%sf(j, k, l) -   &
                                                        13._wp * q_prim_vf(E_idx+i)%sf(j, k-1, l) + &
                                                        2._wp * q_prim_vf(E_idx+i)%sf(j, k-2, l))
                                end do 
                            else 
                                alpha_R(1) = 1._wp
                            end if

                            F_L = (1._wp/60._wp) * (-3._wp * jac(j, k-1, l) + &
                                                27._wp * jac(j, k, l) + &
                                                47._wp * jac(j, k+1, l) -   &
                                                13._wp * jac(j, k+2, l) + &
                                                2._wp * jac(j, k+3, l))
                            !$acc loop seq 
                            do i = -2, 2
                                rho_L = 0._wp
                                !$acc loop seq 
                                do q = 1, num_fluids
                                    rho_L = rho_L + q_prim_vf(q)%sf(j+i,k,l)
                                end do
                                rho_sf(i,1) = rho_L
                            end do

                            !$acc loop seq 
                            do i = -2, 2
                                rho_L = 0._wp
                                !$acc loop seq 
                                do q = 1, num_fluids
                                    rho_L = rho_L + q_prim_vf(q)%sf(j,k+i,l)
                                end do
                                rho_sf(i,2) = rho_L
                            end do

                            !$acc loop seq 
                            do i = -2, 2
                                rho_L = 0._wp
                                !$acc loop seq 
                                do q = 1, num_fluids
                                    rho_L = rho_L + q_prim_vf(q)%sf(j,k,l+i)
                                end do
                                rho_sf(i,3) = rho_L
                            end do

                            rho_L = sum(alpha_rho_L)
                            gamma_L = sum(alpha_L*gammas)
                            pi_inf_L = sum(alpha_L*pi_infs)

                            rho_R = sum(alpha_rho_R)
                            gamma_R = sum(alpha_R*gammas)
                            pi_inf_R = sum(alpha_R*pi_infs)

                            !$acc loop seq 
                            do q = -3, 2
                                !$acc loop seq 
                                do i = 1, num_dims
                                    vel_L(i,q) = (1._wp/60._wp) * (-3._wp * q_prim_vf(momxb+i-1)%sf(j,k-1+q,l) + &
                                                    27._wp * q_prim_vf(momxb+i-1)%sf(j,k+q,l) + &
                                                    47._wp * q_prim_vf(momxb+i-1)%sf(j,k+1+q,l) -   &
                                                    13._wp * q_prim_vf(momxb+i-1)%sf(j,k+2+q,l) + &
                                                    2._wp * q_prim_vf(momxb+i-1)%sf(j,k+3+q,l)) / rho_L
                                end do 

                                !$acc loop seq 
                                do i = 1, num_dims
                                    vel_R(i,q) = (1._wp/60._wp) * (-3._wp * q_prim_vf(momxb+i-1)%sf(j,k+2+q,l) + &
                                                    27._wp * q_prim_vf(momxb+i-1)%sf(j,k+1+q,l) + &
                                                    47._wp * q_prim_vf(momxb+i-1)%sf(j,k+q,l) -   &
                                                    13._wp * q_prim_vf(momxb+i-1)%sf(j,k-1+q,l) + &
                                                    2._wp * q_prim_vf(momxb+i-1)%sf(j,k-2+q,l)) / rho_R
                                end do
                            end do

                            E_L = (1._wp/60._wp) * (-3._wp * q_prim_vf(E_idx)%sf(j,k-1,l) + &
                                                27._wp * q_prim_vf(E_idx)%sf(j,k,l) + &
                                                47._wp * q_prim_vf(E_idx)%sf(j,k+1,l) -   &
                                                13._wp * q_prim_vf(E_idx)%sf(j,k+2,l) + &
                                                2._wp * q_prim_vf(E_idx)%sf(j,k+3,l))

                            E_R = (1._wp/60._wp) * (-3._wp * q_prim_vf(E_idx)%sf(j,k+2,l) + &
                                                27._wp * q_prim_vf(E_idx)%sf(j,k+1,l) + &
                                                47._wp * q_prim_vf(E_idx)%sf(j,k,l) -   &
                                                13._wp * q_prim_vf(E_idx)%sf(j,k-1,l) + &
                                                2._wp * q_prim_vf(E_idx)%sf(j,k-2,l))       

                            pres_L = (E_L - pi_inf_L - 0.5_wp*rho_L*(vel_L(1,0)**2._wp + vel_L(2,0)**2._wp + vel_L(3,0)**2._wp ))/gamma_L                    

                            pres_R = (E_R - pi_inf_R - 0.5_wp*rho_R*(vel_R(1,0)**2._wp + vel_R(2,0)**2._wp + vel_R(3,0)**2._wp ))/gamma_R 

                            a_L = sqrt((pres_L*(1._wp/gamma_L + 1._wp) + pi_inf_L / gamma_L) / rho_L)
                            a_R = sqrt((pres_R*(1._wp/gamma_R + 1._wp) + pi_inf_R / gamma_R) / rho_R)

                            cfl = (max(sqrt(vel_L(1,0)**2._wp + vel_L(2,0)**2._wp + vel_L(3,0)**2._wp),sqrt(vel_R(1,0)**2._wp + vel_R(2,0)**2._wp + vel_R(3,0)**2._wp)) + max(a_L,a_R))

                            !$acc loop seq
                            do i = 1, num_fluids
                                !$acc atomic
                                rhs_vf(i)%sf(j,k+1,l) = rhs_vf(i)%sf(j,k+1,l) + &
                                    (0.5_wp * (alpha_rho_L(i) * &
                                    vel_L(2,0))*(1._wp/dy(k+1)) - &
                                    0.5_wp*cfl * (alpha_rho_L(i))*(1._wp/dy(k+1)))

                                !$acc atomic
                                rhs_vf(i)%sf(j,k,l) = rhs_vf(i)%sf(j,k,l) - &
                                    (0.5_wp * (alpha_rho_L(i) * &
                                    vel_L(2,0))*(1._wp/dy(k)) - &
                                    0.5_wp*cfl * (alpha_rho_L(i))*(1._wp/dy(k))) 
                            end do

                            if(num_fluids > 1) then 
                                !$acc loop seq
                                do i = 1, num_fluids
                                    !$acc atomic
                                    rhs_vf(advxb+i-1)%sf(j,k+1,l) = rhs_vf(advxb+i-1)%sf(j,k+1,l) + &
                                        (0.5_wp * (alpha_L(i) * &
                                        vel_L(2,0))*(1._wp/dy(k+1)) - &
                                        0.5_wp*cfl*(alpha_L(i))*(1._wp/dy(k+1)))

                                    !$acc atomic
                                    rhs_vf(advxb+i-1)%sf(j,k+1,l) = rhs_vf(advxb+i-1)%sf(j,k+1,l) &
                                    - (0.5_wp * q_prim_vf(advxb+i-1)%sf(j,k+1,l)  * vel_L(2,0)*(1._wp/dy(k+1)))

                                    !$acc atomic
                                    rhs_vf(advxb+i-1)%sf(j,k,l) = rhs_vf(advxb+i-1)%sf(j,k,l) - &
                                        (0.5_wp * (alpha_L(i) * &
                                        vel_L(2,0))*(1._wp/dy(k)) - &
                                        0.5_wp*cfl*(alpha_L(i))*(1._wp/dy(k)))

                                    !$acc atomic
                                    rhs_vf(advxb+i-1)%sf(j,k,l) = rhs_vf(advxb+i-1)%sf(j,k,l) &
                                    + (0.5_wp * q_prim_vf(advxb+i-1)%sf(j,k,l)  * vel_L(2,0)*(1._wp/dy(k)))
                                end do
                            end if

                            !$acc atomic
                            rhs_vf(momxb+1)%sf(j,k+1,l) = rhs_vf(momxb+1)%sf(j,k+1,l) + &
                                 (0.5_wp* (rho_L * (vel_L(2,0))**2.0 + &
                                 pres_L+F_L)*(1._wp/dy(k+1)) - &
                                 0.5_wp*cfl * (rho_L*vel_L(2,0))*(1._wp/dy(k+1)))

                            !$acc atomic
                            rhs_vf(momxb)%sf(j, k+1, l) =  rhs_vf(momxb)%sf(j, k+1, l) + &
                                (0.5_wp * rho_L * vel_L(1,0)*vel_L(2,0)*(1._wp/dy(k+1)) - &
                                0.5_wp*cfl * (rho_L*vel_L(1,0))*(1._wp/dy(k+1)))

                            !$acc atomic
                            rhs_vf(momxb+2)%sf(j, k+1, l) =  rhs_vf(momxb+2)%sf(j, k+1, l) + &
                                (0.5_wp * rho_L * vel_L(3,0)*vel_L(2,0)*(1._wp/dy(k+1)) - &
                                0.5_wp*cfl * (rho_L*vel_L(3,0))*(1._wp/dy(k+1)))

                            !$acc atomic
                            rhs_vf(E_idx)%sf(j, k+1, l) = rhs_vf(E_idx)%sf(j, k+1, l) +  &
                                (0.5_wp * (vel_L(2,0) * (E_L + &
                                pres_L+F_L) )*(1._wp/dy(k+1)) - &
                                0.5_wp*cfl * (E_L)*(1._wp/dy(k+1)))

                            !$acc atomic
                            rhs_vf(momxb+1)%sf(j,k,l) = rhs_vf(momxb+1)%sf(j,k,l) - &
                                 (0.5_wp* (rho_L * (vel_L(2,0))**2.0 + &
                                 pres_L+F_L)*(1._wp/dy(k)) - &
                                 0.5_wp*cfl * (rho_L*vel_L(2,0))*(1._wp/dy(k)))

                            !$acc atomic
                            rhs_vf(momxb)%sf(j, k, l) =  rhs_vf(momxb)%sf(j, k, l) - &
                                (0.5_wp * rho_L * vel_L(1,0)*vel_L(2,0)*(1._wp/dy(k)) - &
                                0.5_wp*cfl * (rho_L*vel_L(1,0))*(1._wp/dy(k)))

                            !$acc atomic
                            rhs_vf(momxb+2)%sf(j, k, l) =  rhs_vf(momxb+2)%sf(j, k, l) - &
                                (0.5_wp * rho_L * vel_L(3,0)*vel_L(2,0)*(1._wp/dy(k)) - &
                                0.5_wp*cfl * (rho_L*vel_L(3,0))*(1._wp/dy(k)))

                            !$acc atomic
                            rhs_vf(E_idx)%sf(j, k, l) = rhs_vf(E_idx)%sf(j, k, l) -  &
                                (0.5_wp * (vel_L(2,0) * (E_L + &
                                pres_L+F_L) )*(1._wp/dy(k)) - &
                                0.5_wp*cfl * (E_L)*(1._wp/dy(k)))

                            !! duy & dvx
                            dvel1 = (1/(12._wp*dy(k))) * ( &
                            8._wp*q_prim_vf(momxb)%sf(j,k+1,l)/rho_sf(1,2) - &
                            8._wp*q_prim_vf(momxb)%sf(j,k-1,l)/rho_sf(-1,2) + &
                            q_prim_vf(momxb)%sf(j,k-2,l)/rho_sf(-2,2) - &
                            q_prim_vf(momxb)%sf(j,k+2,l)/rho_sf(2,2) )

                            dvel2 = (1/(12._wp*dx(j))) * ( &
                            8._wp*q_prim_vf(momxb+1)%sf(j+1,k,l)/rho_sf(1,1) - &
                            8._wp*q_prim_vf(momxb+1)%sf(j-1,k,l)/rho_sf(-1,1) + &
                            q_prim_vf(momxb+1)%sf(j-2,k,l)/rho_sf(-2,1) - &
                            q_prim_vf(momxb+1)%sf(j+2,k,l)/rho_sf(2,1) )

                            if(num_fluids > 1) then 
                                !$acc loop seq
                                do q = -3, 2 
                                    mu_L(q) = 0._wp
                                    mu_R(q) = 0._wp
                                    !$acc loop seq 
                                    do i = 1, num_fluids
                                        mu_L(q) = (1._wp/60._wp) * (-3._wp * q_prim_vf(E_idx+i)%sf(j, k-1+q, l) + &
                                                27._wp * q_prim_vf(E_idx+i)%sf(j, k+q, l) + &
                                                47._wp * q_prim_vf(E_idx+i)%sf(j, k+1+q, l) -   &
                                                13._wp * q_prim_vf(E_idx+i)%sf(j, k+2+q, l) + &
                                                2._wp * q_prim_vf(E_idx+i)%sf(j, k+3+q, l)) / Res(1, i) + mu_L(q)
                                        
                                        mu_R(q) = (1._wp/60._wp) * (-3._wp * q_prim_vf(E_idx+i)%sf(j, k+2+q, l) + &
                                                27._wp * q_prim_vf(E_idx+i)%sf(j, k+1+q, l) + &
                                                47._wp * q_prim_vf(E_idx+i)%sf(j, k+q, l) -   &
                                                13._wp * q_prim_vf(E_idx+i)%sf(j, k-1+q, l) + &
                                                2._wp * q_prim_vf(E_idx+i)%sf(j, k-2+q, l)) / Res(1, i) + mu_R(q)
                                    end do
                                end do
                            else
                                !$acc loop seq
                                do q = -3, 2
                                    mu_L(q) = 1._wp / Res(1, 1) 
                                    mu_R(q) = 1._wp / Res(1, 1) 
                                end do
                            end if
                            
                            if(viscous) then

                                !$acc atomic
                                rhs_vf(momxb)%sf(j,k+2,l) = rhs_vf(momxb)%sf(j,k+2,l) - 0.5_wp*mu_L(1)*((-3._wp/60._wp)*(dvel1+dvel2))*(1._wp/dy(k+2))
                                !$acc atomic
                                rhs_vf(momxb)%sf(j,k+1,l) = rhs_vf(momxb)%sf(j,k+1,l) - 0.5_wp*mu_L(0)*((27._wp/60._wp)*(dvel1+dvel2))*(1._wp/dy(k+1))
                                !$acc atomic
                                rhs_vf(momxb)%sf(j,k,l) = rhs_vf(momxb)%sf(j,k,l) - 0.5_wp*mu_L(-1)*((47._wp/60._wp)*(dvel1+dvel2))*(1._wp/dy(k))
                                !$acc atomic
                                rhs_vf(momxb)%sf(j,k-1,l) = rhs_vf(momxb)%sf(j,k-1,l) - 0.5_wp*mu_L(-2)*((-13._wp/60._wp)*(dvel1+dvel2))*(1._wp/dy(k-1))
                                !$acc atomic
                                rhs_vf(momxb)%sf(j,k-2,l) = rhs_vf(momxb)%sf(j,k-2,l) - 0.5_wp*mu_L(-3)*((2._wp/60._wp)*(dvel1+dvel2))*(1._wp/dy(k-2))

                                !$acc atomic
                                rhs_vf(E_idx)%sf(j,k+2,l) = rhs_vf(E_idx)%sf(j,k+2,l) - 0.5_wp*mu_L(1)*((-3._wp/60._wp)*(dvel1+dvel2))*vel_L(1,1)*(1._wp/dy(k+2))
                                !$acc atomic
                                rhs_vf(E_idx)%sf(j,k+1,l) = rhs_vf(E_idx)%sf(j,k+1,l) - 0.5_wp*mu_L(0)*((27._wp/60._wp)*(dvel1+dvel2))*vel_L(1,0)*(1._wp/dy(k+1))
                                !$acc atomic
                                rhs_vf(E_idx)%sf(j,k,l) = rhs_vf(E_idx)%sf(j,k,l) - 0.5_wp*mu_L(-1)*((47._wp/60._wp)*(dvel1+dvel2))*vel_L(1,-1)*(1._wp/dy(k))
                                !$acc atomic
                                rhs_vf(E_idx)%sf(j,k-1,l) = rhs_vf(E_idx)%sf(j,k-1,l) - 0.5_wp*mu_L(-2)*((-13._wp/60._wp)*(dvel1+dvel2))*vel_L(1,-2)*(1._wp/dy(k-1))
                                !$acc atomic
                                rhs_vf(E_idx)%sf(j,k-2,l) = rhs_vf(E_idx)%sf(j,k-2,l) - 0.5_wp*mu_L(-3)*((2._wp/60._wp)*(dvel1+dvel2))*vel_L(1,-3)*(1._wp/dy(k-2))

                                !$acc atomic
                                rhs_vf(momxb)%sf(j,k+1,l) = rhs_vf(momxb)%sf(j,k+1,l) + 0.5_wp*mu_L(1)*((-3._wp/60._wp)*(dvel1+dvel2))*(1._wp/dy(k+1))
                                !$acc atomic
                                rhs_vf(momxb)%sf(j,k,l) = rhs_vf(momxb)%sf(j,k,l) + 0.5_wp*mu_L(0)*((27._wp/60._wp)*(dvel1+dvel2))*(1._wp/dy(k))
                                !$acc atomic
                                rhs_vf(momxb)%sf(j,k-1,l) = rhs_vf(momxb)%sf(j,k-1,l) + 0.5_wp*mu_L(-1)*((47._wp/60._wp)*(dvel1+dvel2))*(1._wp/dy(k-1))
                                !$acc atomic
                                rhs_vf(momxb)%sf(j,k-2,l) = rhs_vf(momxb)%sf(j,k-2,l) + 0.5_wp*mu_L(-2)*((-13._wp/60._wp)*(dvel1+dvel2))*(1._wp/dy(k-2))
                                !$acc atomic
                                rhs_vf(momxb)%sf(j,k-3,l) = rhs_vf(momxb)%sf(j,k-3,l) + 0.5_wp*mu_L(-3)*((2._wp/60._wp)*(dvel1+dvel2))*(1._wp/dy(k-3))

                                !$acc atomic
                                rhs_vf(E_idx)%sf(j,k+1,l) = rhs_vf(E_idx)%sf(j,k+1,l) + 0.5_wp*mu_L(1)*((-3._wp/60._wp)*(dvel1+dvel2))*vel_L(1,1)*(1._wp/dy(k+1))
                                !$acc atomic
                                rhs_vf(E_idx)%sf(j,k,l) = rhs_vf(E_idx)%sf(j,k,l) + 0.5_wp*mu_L(0)*((27._wp/60._wp)*(dvel1+dvel2))*vel_L(1,0)*(1._wp/dy(k))
                                !$acc atomic
                                rhs_vf(E_idx)%sf(j,k-1,l) = rhs_vf(E_idx)%sf(j,k-1,l) + 0.5_wp*mu_L(-1)*((47._wp/60._wp)*(dvel1+dvel2))*vel_L(1,-1)*(1._wp/dy(k-1))
                                !$acc atomic
                                rhs_vf(E_idx)%sf(j,k-2,l) = rhs_vf(E_idx)%sf(j,k-2,l) + 0.5_wp*mu_L(-2)*((-13._wp/60._wp)*(dvel1+dvel2))*vel_L(1,-2)*(1._wp/dy(k-2))
                                !$acc atomic
                                rhs_vf(E_idx)%sf(j,k-3,l) = rhs_vf(E_idx)%sf(j,k-3,l) + 0.5_wp*mu_L(-3)*((2._wp/60._wp)*(dvel1+dvel2))*vel_L(1,-3)*(1._wp/dy(k-3))
                               
                                !$acc atomic
                                rhs_vf(momxb)%sf(j,k-1,l) = rhs_vf(momxb)%sf(j,k-1,l) - 0.5_wp*mu_R(-2)*((-3._wp/60._wp)*(dvel1+dvel2))*(1._wp/dy(k-1))
                                !$acc atomic
                                rhs_vf(momxb)%sf(j,k,l) = rhs_vf(momxb)%sf(j,k,l) - 0.5_wp*mu_R(-1)*((27._wp/60._wp)*(dvel1+dvel2))*(1._wp/dy(k))
                                !$acc atomic
                                rhs_vf(momxb)%sf(j,k+1,l) = rhs_vf(momxb)%sf(j,k+1,l) - 0.5_wp*mu_R(0)*((47._wp/60._wp)*(dvel1+dvel2))*(1._wp/dy(k+1))
                                !$acc atomic
                                rhs_vf(momxb)%sf(j,k+2,l) = rhs_vf(momxb)%sf(j,k+2,l) - 0.5_wp*mu_R(1)*((-13._wp/60._wp)*(dvel1+dvel2))*(1._wp/dy(k+2))
                                !$acc atomic
                                rhs_vf(momxb)%sf(j,k+3,l) = rhs_vf(momxb)%sf(j,k+3,l) - 0.5_wp*mu_R(2)*((2._wp/60._wp)*(dvel1+dvel2))*(1._wp/dy(k+3))

                                !$acc atomic
                                rhs_vf(E_idx)%sf(j,k-1,l) = rhs_vf(E_idx)%sf(j,k-1,l) - 0.5_wp*mu_R(-2)*((-3._wp/60._wp)*(dvel1+dvel2))*vel_R(1,-2)*(1._wp/dy(k-1))
                                !$acc atomic
                                rhs_vf(E_idx)%sf(j,k,l) = rhs_vf(E_idx)%sf(j,k,l) - 0.5_wp*mu_R(-1)*((27._wp/60._wp)*(dvel1+dvel2))*vel_R(1,-1)*(1._wp/dy(k))
                                !$acc atomic
                                rhs_vf(E_idx)%sf(j,k+1,l) = rhs_vf(E_idx)%sf(j,k+1,l) - 0.5_wp*mu_R(0)*((47._wp/60._wp)*(dvel1+dvel2))*vel_R(1,0)*(1._wp/dy(k+1))
                                !$acc atomic
                                rhs_vf(E_idx)%sf(j,k+2,l) = rhs_vf(E_idx)%sf(j,k+2,l) - 0.5_wp*mu_R(1)*((-13._wp/60._wp)*(dvel1+dvel2))*vel_R(1,1)*(1._wp/dy(k+2))
                                !$acc atomic
                                rhs_vf(E_idx)%sf(j,k+3,l) = rhs_vf(E_idx)%sf(j,k+3,l) - 0.5_wp*mu_R(2)*((2._wp/60._wp)*(dvel1+dvel2))*vel_R(1,2)*(1._wp/dy(k+3)) 

                                !$acc atomic
                                rhs_vf(momxb)%sf(j,k-2,l) = rhs_vf(momxb)%sf(j,k-2,l) + 0.5_wp*mu_R(-2)*((-3._wp/60._wp)*(dvel1+dvel2))*(1._wp/dy(k-2))
                                !$acc atomic
                                rhs_vf(momxb)%sf(j,k-1,l) = rhs_vf(momxb)%sf(j,k-1,l) + 0.5_wp*mu_R(-1)*((27._wp/60._wp)*(dvel1+dvel2))*(1._wp/dy(k-1))
                                !$acc atomic
                                rhs_vf(momxb)%sf(j,k,l) = rhs_vf(momxb)%sf(j,k,l) + 0.5_wp*mu_R(0)*((47._wp/60._wp)*(dvel1+dvel2))*(1._wp/dy(k))
                                !$acc atomic
                                rhs_vf(momxb)%sf(j,k+1,l) = rhs_vf(momxb)%sf(j,k+1,l) + 0.5_wp*mu_R(1)*((-13._wp/60._wp)*(dvel1+dvel2))*(1._wp/dy(k+1))
                                !$acc atomic
                                rhs_vf(momxb)%sf(j,k+2,l) = rhs_vf(momxb)%sf(j,k+2,l) + 0.5_wp*mu_R(2)*((2._wp/60._wp)*(dvel1+dvel2))*(1._wp/dy(k+2))

                                !$acc atomic
                                rhs_vf(E_idx)%sf(j,k-2,l) = rhs_vf(E_idx)%sf(j,k-2,l) + 0.5_wp*mu_R(-2)*((-3._wp/60._wp)*(dvel1+dvel2))*vel_R(1,-2)*(1._wp/dy(k-2))
                                !$acc atomic
                                rhs_vf(E_idx)%sf(j,k-1,l) = rhs_vf(E_idx)%sf(j,k-1,l) + 0.5_wp*mu_R(-1)*((27._wp/60._wp)*(dvel1+dvel2))*vel_R(1,-1)*(1._wp/dy(k-1))
                                !$acc atomic
                                rhs_vf(E_idx)%sf(j,k,l) = rhs_vf(E_idx)%sf(j,k,l) + 0.5_wp*mu_R(0)*((47._wp/60._wp)*(dvel1+dvel2))*vel_R(1,0)*(1._wp/dy(k))
                                !$acc atomic
                                rhs_vf(E_idx)%sf(j,k+1,l) = rhs_vf(E_idx)%sf(j,k+1,l) + 0.5_wp*mu_R(1)*((-13._wp/60._wp)*(dvel1+dvel2))*vel_R(1,1)*(1._wp/dy(k+1))
                                !$acc atomic
                                rhs_vf(E_idx)%sf(j,k+2,l) = rhs_vf(E_idx)%sf(j,k+2,l) + 0.5_wp*mu_R(2)*((2._wp/60._wp)*(dvel1+dvel2))*vel_R(1,2)*(1._wp/dy(k+2))

                            end if

                            !! dwy & dvz
                            dvel1 = (1/(12._wp*dy(k))) * ( &
                            8._wp*q_prim_vf(momxb+2)%sf(j,k+1,l)/rho_sf(1,2) - &
                            8._wp*q_prim_vf(momxb+2)%sf(j,k-1,l)/rho_sf(-1,2) + &
                            q_prim_vf(momxb+2)%sf(j,k-2,l)/rho_sf(-2,2) - &
                            q_prim_vf(momxb+2)%sf(j,k+2,l)/rho_sf(2,2) )

                            dvel2 = (1/(12._wp*dz(l))) * ( &
                            8._wp*q_prim_vf(momxb+1)%sf(j,k,l+1)/rho_sf(1,3)  - &
                            8._wp*q_prim_vf(momxb+1)%sf(j,k,l-1)/rho_sf(-1,3) + &
                            q_prim_vf(momxb+1)%sf(j,k,l-2)/rho_sf(-2,3) - &
                            q_prim_vf(momxb+1)%sf(j,k,l+2)/rho_sf(2,3) )

                            if(viscous) then

                                !$acc atomic
                                rhs_vf(momxb+2)%sf(j,k+2,l) = rhs_vf(momxb+2)%sf(j,k+2,l) - 0.5_wp*mu_L(1)*((-3._wp/60._wp)*(dvel1+dvel2))*(1._wp/dy(k+2))
                                !$acc atomic
                                rhs_vf(momxb+2)%sf(j,k+1,l) = rhs_vf(momxb+2)%sf(j,k+1,l) - 0.5_wp*mu_L(0)*((27._wp/60._wp)*(dvel1+dvel2))*(1._wp/dy(k+1))
                                !$acc atomic
                                rhs_vf(momxb+2)%sf(j,k,l) = rhs_vf(momxb+2)%sf(j,k,l) - 0.5_wp*mu_L(-1)*((47._wp/60._wp)*(dvel1+dvel2))*(1._wp/dy(k))
                                !$acc atomic
                                rhs_vf(momxb+2)%sf(j,k-1,l) = rhs_vf(momxb+2)%sf(j,k-1,l) - 0.5_wp*mu_L(-2)*((-13._wp/60._wp)*(dvel1+dvel2))*(1._wp/dy(k-1))
                                !$acc atomic
                                rhs_vf(momxb+2)%sf(j,k-2,l) = rhs_vf(momxb+2)%sf(j,k-2,l) - 0.5_wp*mu_L(-3)*((2._wp/60._wp)*(dvel1+dvel2))*(1._wp/dy(k-2))

                                !$acc atomic
                                rhs_vf(E_idx)%sf(j,k+2,l) = rhs_vf(E_idx)%sf(j,k+2,l) - 0.5_wp*mu_L(1)*((-3._wp/60._wp)*(dvel1+dvel2))*vel_L(3,1)*(1._wp/dy(k+2))
                                !$acc atomic
                                rhs_vf(E_idx)%sf(j,k+1,l) = rhs_vf(E_idx)%sf(j,k+1,l) - 0.5_wp*mu_L(0)*((27._wp/60._wp)*(dvel1+dvel2))*vel_L(3,0)*(1._wp/dy(k+1))
                                !$acc atomic
                                rhs_vf(E_idx)%sf(j,k,l) = rhs_vf(E_idx)%sf(j,k,l) - 0.5_wp*mu_L(-1)*((47._wp/60._wp)*(dvel1+dvel2))*vel_L(3,-1)*(1._wp/dy(k))
                                !$acc atomic
                                rhs_vf(E_idx)%sf(j,k-1,l) = rhs_vf(E_idx)%sf(j,k-1,l) - 0.5_wp*mu_L(-2)*((-13._wp/60._wp)*(dvel1+dvel2))*vel_L(3,-2)*(1._wp/dy(k-1))
                                !$acc atomic
                                rhs_vf(E_idx)%sf(j,k-2,l) = rhs_vf(E_idx)%sf(j,k-2,l) - 0.5_wp*mu_L(-3)*((2._wp/60._wp)*(dvel1+dvel2))*vel_L(3,-3)*(1._wp/dy(k-2))

                                !$acc atomic
                                rhs_vf(momxb+2)%sf(j,k+1,l) = rhs_vf(momxb+2)%sf(j,k+1,l) + 0.5_wp*mu_L(1)*((-3._wp/60._wp)*(dvel1+dvel2))*(1._wp/dy(k+1))
                                !$acc atomic
                                rhs_vf(momxb+2)%sf(j,k,l) = rhs_vf(momxb+2)%sf(j,k,l) + 0.5_wp*mu_L(0)*((27._wp/60._wp)*(dvel1+dvel2))*(1._wp/dy(k))
                                !$acc atomic
                                rhs_vf(momxb+2)%sf(j,k-1,l) = rhs_vf(momxb+2)%sf(j,k-1,l) + 0.5_wp*mu_L(-1)*((47._wp/60._wp)*(dvel1+dvel2))*(1._wp/dy(k-1))
                                !$acc atomic
                                rhs_vf(momxb+2)%sf(j,k-2,l) = rhs_vf(momxb+2)%sf(j,k-2,l) + 0.5_wp*mu_L(-2)*((-13._wp/60._wp)*(dvel1+dvel2))*(1._wp/dy(k-2))
                                !$acc atomic
                                rhs_vf(momxb+2)%sf(j,k-3,l) = rhs_vf(momxb+2)%sf(j,k-3,l) + 0.5_wp*mu_L(-3)*((2._wp/60._wp)*(dvel1+dvel2))*(1._wp/dy(k-3))

                                !$acc atomic
                                rhs_vf(E_idx)%sf(j,k+1,l) = rhs_vf(E_idx)%sf(j,k+1,l) + 0.5_wp*mu_L(1)*((-3._wp/60._wp)*(dvel1+dvel2))*vel_L(3,1)*(1._wp/dy(k+1))
                                !$acc atomic
                                rhs_vf(E_idx)%sf(j,k,l) = rhs_vf(E_idx)%sf(j,k,l) + 0.5_wp*mu_L(0)*((27._wp/60._wp)*(dvel1+dvel2))*vel_L(3,0)*(1._wp/dy(k))
                                !$acc atomic
                                rhs_vf(E_idx)%sf(j,k-1,l) = rhs_vf(E_idx)%sf(j,k-1,l) + 0.5_wp*mu_L(-1)*((47._wp/60._wp)*(dvel1+dvel2))*vel_L(3,-1)*(1._wp/dy(k-1))
                                !$acc atomic
                                rhs_vf(E_idx)%sf(j,k-2,l) = rhs_vf(E_idx)%sf(j,k-2,l) + 0.5_wp*mu_L(-2)*((-13._wp/60._wp)*(dvel1+dvel2))*vel_L(3,-2)*(1._wp/dy(k-2))
                                !$acc atomic
                                rhs_vf(E_idx)%sf(j,k-3,l) = rhs_vf(E_idx)%sf(j,k-3,l) + 0.5_wp*mu_L(-3)*((2._wp/60._wp)*(dvel1+dvel2))*vel_L(3,-3)*(1._wp/dy(k-3))

                                !$acc atomic
                                rhs_vf(momxb+2)%sf(j,k-1,l) = rhs_vf(momxb+2)%sf(j,k-1,l) - 0.5_wp*mu_R(-2)*((-3._wp/60._wp)*(dvel1+dvel2))*(1._wp/dy(k-1))
                                !$acc atomic
                                rhs_vf(momxb+2)%sf(j,k,l) = rhs_vf(momxb+2)%sf(j,k,l) - 0.5_wp*mu_R(-1)*((27._wp/60._wp)*(dvel1+dvel2))*(1._wp/dy(k))
                                !$acc atomic
                                rhs_vf(momxb+2)%sf(j,k+1,l) = rhs_vf(momxb+2)%sf(j,k+1,l) - 0.5_wp*mu_R(0)*((47._wp/60._wp)*(dvel1+dvel2))*(1._wp/dy(k+1))
                                !$acc atomic
                                rhs_vf(momxb+2)%sf(j,k+2,l) = rhs_vf(momxb+2)%sf(j,k+2,l) - 0.5_wp*mu_R(1)*((-13._wp/60._wp)*(dvel1+dvel2))*(1._wp/dy(k+2))
                                !$acc atomic
                                rhs_vf(momxb+2)%sf(j,k+3,l) = rhs_vf(momxb+2)%sf(j,k+3,l) - 0.5_wp*mu_R(2)*((2._wp/60._wp)*(dvel1+dvel2))*(1._wp/dy(k+3))

                                !$acc atomic
                                rhs_vf(E_idx)%sf(j,k-1,l) = rhs_vf(E_idx)%sf(j,k-1,l) - 0.5_wp*mu_R(-2)*((-3._wp/60._wp)*(dvel1+dvel2))*vel_R(3,-2)*(1._wp/dy(k-1))
                                !$acc atomic
                                rhs_vf(E_idx)%sf(j,k,l) = rhs_vf(E_idx)%sf(j,k,l) - 0.5_wp*mu_R(-1)*((27._wp/60._wp)*(dvel1+dvel2))*vel_R(3,-1)*(1._wp/dy(k))
                                !$acc atomic
                                rhs_vf(E_idx)%sf(j,k+1,l) = rhs_vf(E_idx)%sf(j,k+1,l) - 0.5_wp*mu_R(0)*((47._wp/60._wp)*(dvel1+dvel2))*vel_R(3,0)*(1._wp/dy(k+1))
                                !$acc atomic
                                rhs_vf(E_idx)%sf(j,k+2,l) = rhs_vf(E_idx)%sf(j,k+2,l) - 0.5_wp*mu_R(1)*((-13._wp/60._wp)*(dvel1+dvel2))*vel_R(3,1)*(1._wp/dy(k+2))
                                !$acc atomic
                                rhs_vf(E_idx)%sf(j,k+3,l) = rhs_vf(E_idx)%sf(j,k+3,l) - 0.5_wp*mu_R(2)*((2._wp/60._wp)*(dvel1+dvel2))*vel_R(3,2)*(1._wp/dy(k+3)) 

                                !$acc atomic
                                rhs_vf(momxb+2)%sf(j,k-2,l) = rhs_vf(momxb+2)%sf(j,k-2,l) + 0.5_wp*mu_R(-2)*((-3._wp/60._wp)*(dvel1+dvel2))*(1._wp/dy(k-2))
                                !$acc atomic
                                rhs_vf(momxb+2)%sf(j,k-1,l) = rhs_vf(momxb+2)%sf(j,k-1,l) + 0.5_wp*mu_R(-1)*((27._wp/60._wp)*(dvel1+dvel2))*(1._wp/dy(k-1))
                                !$acc atomic
                                rhs_vf(momxb+2)%sf(j,k,l) = rhs_vf(momxb+2)%sf(j,k,l) + 0.5_wp*mu_R(0)*((47._wp/60._wp)*(dvel1+dvel2))*(1._wp/dy(k))
                                !$acc atomic
                                rhs_vf(momxb+2)%sf(j,k+1,l) = rhs_vf(momxb+2)%sf(j,k+1,l) + 0.5_wp*mu_R(1)*((-13._wp/60._wp)*(dvel1+dvel2))*(1._wp/dy(k+1))
                                !$acc atomic
                                rhs_vf(momxb+2)%sf(j,k+2,l) = rhs_vf(momxb+2)%sf(j,k+2,l) + 0.5_wp*mu_R(2)*((2._wp/60._wp)*(dvel1+dvel2))*(1._wp/dy(k+2))

                                !$acc atomic
                                rhs_vf(E_idx)%sf(j,k-2,l) = rhs_vf(E_idx)%sf(j,k-2,l) + 0.5_wp*mu_R(-2)*((-3._wp/60._wp)*(dvel1+dvel2))*vel_R(3,-2)*(1._wp/dy(k-2))
                                !$acc atomic
                                rhs_vf(E_idx)%sf(j,k-1,l) = rhs_vf(E_idx)%sf(j,k-1,l) + 0.5_wp*mu_R(-1)*((27._wp/60._wp)*(dvel1+dvel2))*vel_R(3,-1)*(1._wp/dy(k-1))
                                !$acc atomic
                                rhs_vf(E_idx)%sf(j,k,l) = rhs_vf(E_idx)%sf(j,k,l) + 0.5_wp*mu_R(0)*((47._wp/60._wp)*(dvel1+dvel2))*vel_R(3,0)*(1._wp/dy(k))
                                !$acc atomic
                                rhs_vf(E_idx)%sf(j,k+1,l) = rhs_vf(E_idx)%sf(j,k+1,l) + 0.5_wp*mu_R(1)*((-13._wp/60._wp)*(dvel1+dvel2))*vel_R(3,1)*(1._wp/dy(k+1))
                                !$acc atomic
                                rhs_vf(E_idx)%sf(j,k+2,l) = rhs_vf(E_idx)%sf(j,k+2,l) + 0.5_wp*mu_R(2)*((2._wp/60._wp)*(dvel1+dvel2))*vel_R(3,2)*(1._wp/dy(k+2))

                            end if

                            !! dvy & dux
                            dvel1 = (1/(12._wp*dy(k))) * ( &
                            8._wp*q_prim_vf(momxb+1)%sf(j,k+1,l)/rho_sf(1,2) - &
                            8._wp*q_prim_vf(momxb+1)%sf(j,k-1,l)/rho_sf(-1,2) + &
                            q_prim_vf(momxb+1)%sf(j,k-2,l)/rho_sf(-2,2) - &
                            q_prim_vf(momxb+1)%sf(j,k+2,l)/rho_sf(2,2) )

                            dvel2 = (1/(12._wp*dx(j))) * ( &
                            8._wp*q_prim_vf(momxb)%sf(j+1,k,l)/rho_sf(1,1) - &
                            8._wp*q_prim_vf(momxb)%sf(j-1,k,l)/rho_sf(-1,1) + &
                            q_prim_vf(momxb)%sf(j-2,k,l)/rho_sf(-2,1) - &
                            q_prim_vf(momxb)%sf(j+2,k,l)/rho_sf(2,1) )

                            if(viscous) then

                                !$acc atomic
                                rhs_vf(momxb+1)%sf(j,k+2,l) = rhs_vf(momxb+1)%sf(j,k+2,l) - 0.5_wp*mu_L(1)*((-3._wp/60._wp)*(4._wp*dvel1-2._wp*dvel2)/3._wp)*(1._wp/dy(k+2))
                                !$acc atomic
                                rhs_vf(momxb+1)%sf(j,k+1,l) = rhs_vf(momxb+1)%sf(j,k+1,l) - 0.5_wp*mu_L(0)*((27._wp/60._wp)*(4._wp*dvel1-2._wp*dvel2)/3._wp)*(1._wp/dy(k+1))
                                !$acc atomic
                                rhs_vf(momxb+1)%sf(j,k,l) = rhs_vf(momxb+1)%sf(j,k,l) - 0.5_wp*mu_L(-1)*((47._wp/60._wp)*(4._wp*dvel1-2._wp*dvel2)/3._wp)*(1._wp/dy(k))
                                !$acc atomic
                                rhs_vf(momxb+1)%sf(j,k-1,l) = rhs_vf(momxb+1)%sf(j,k-1,l) - 0.5_wp*mu_L(-2)*((-13._wp/60._wp)*(4._wp*dvel1-2._wp*dvel2)/3._wp)*(1._wp/dy(k-1))
                                !$acc atomic
                                rhs_vf(momxb+1)%sf(j,k-2,l) = rhs_vf(momxb+1)%sf(j,k-2,l) - 0.5_wp*mu_L(-3)*((2._wp/60._wp)*(4._wp*dvel1-2._wp*dvel2)/3._wp)*(1._wp/dy(k-2))

                                !$acc atomic
                                rhs_vf(E_idx)%sf(j,k+2,l) = rhs_vf(E_idx)%sf(j,k+2,l) - 0.5_wp*mu_L(1)*((-3._wp/60._wp)*(4._wp*dvel1-2._wp*dvel2)/3._wp)*vel_L(2,1)*(1._wp/dy(k+2))
                                !$acc atomic
                                rhs_vf(E_idx)%sf(j,k+1,l) = rhs_vf(E_idx)%sf(j,k+1,l) - 0.5_wp*mu_L(0)*((27._wp/60._wp)*(4._wp*dvel1-2._wp*dvel2)/3._wp)*vel_L(2,0)*(1._wp/dy(k+1))
                                !$acc atomic
                                rhs_vf(E_idx)%sf(j,k,l) = rhs_vf(E_idx)%sf(j,k,l) - 0.5_wp*mu_L(-1)*((47._wp/60._wp)*(4._wp*dvel1-2._wp*dvel2)/3._wp)*vel_L(2,-1)*(1._wp/dy(k))
                                !$acc atomic
                                rhs_vf(E_idx)%sf(j,k-1,l) = rhs_vf(E_idx)%sf(j,k-1,l) - 0.5_wp*mu_L(-2)*((-13._wp/60._wp)*(4._wp*dvel1-2._wp*dvel2)/3._wp)*vel_L(2,-2)*(1._wp/dy(k-1))
                                !$acc atomic
                                rhs_vf(E_idx)%sf(j,k-2,l) = rhs_vf(E_idx)%sf(j,k-2,l) - 0.5_wp*mu_L(-3)*((2._wp/60._wp)*(4._wp*dvel1-2._wp*dvel2)/3._wp)*vel_L(2,-3)*(1._wp/dy(k-2))

                                !$acc atomic
                                rhs_vf(momxb+1)%sf(j,k+1,l) = rhs_vf(momxb+1)%sf(j,k+1,l) + 0.5_wp*mu_L(1)*((-3._wp/60._wp)*(4._wp*dvel1-2._wp*dvel2)/3._wp)*(1._wp/dy(k+1))
                                !$acc atomic
                                rhs_vf(momxb+1)%sf(j,k,l) = rhs_vf(momxb+1)%sf(j,k,l) + 0.5_wp*mu_L(0)*((27._wp/60._wp)*(4._wp*dvel1-2._wp*dvel2)/3._wp)*(1._wp/dy(k))
                                !$acc atomic
                                rhs_vf(momxb+1)%sf(j,k-1,l) = rhs_vf(momxb+1)%sf(j,k-1,l) + 0.5_wp*mu_L(-1)*((47._wp/60._wp)*(4._wp*dvel1-2._wp*dvel2)/3._wp)*(1._wp/dy(k-1))
                                !$acc atomic
                                rhs_vf(momxb+1)%sf(j,k-2,l) = rhs_vf(momxb+1)%sf(j,k-2,l) + 0.5_wp*mu_L(-2)*((-13._wp/60._wp)*(4._wp*dvel1-2._wp*dvel2)/3._wp)*(1._wp/dy(k-2))
                                !$acc atomic
                                rhs_vf(momxb+1)%sf(j,k-3,l) = rhs_vf(momxb+1)%sf(j,k-3,l) + 0.5_wp*mu_L(-3)*((2._wp/60._wp)*(4._wp*dvel1-2._wp*dvel2)/3._wp)*(1._wp/dy(k-3))

                                !$acc atomic
                                rhs_vf(E_idx)%sf(j,k+1,l) = rhs_vf(E_idx)%sf(j,k+1,l) + 0.5_wp*mu_L(1)*((-3._wp/60._wp)*(4._wp*dvel1-2._wp*dvel2)/3._wp)*vel_L(2,1)*(1._wp/dy(k+1))
                                !$acc atomic
                                rhs_vf(E_idx)%sf(j,k,l) = rhs_vf(E_idx)%sf(j,k,l) + 0.5_wp*mu_L(0)*((27._wp/60._wp)*(4._wp*dvel1-2._wp*dvel2)/3._wp)*vel_L(2,0)*(1._wp/dy(k))
                                !$acc atomic
                                rhs_vf(E_idx)%sf(j,k-1,l) = rhs_vf(E_idx)%sf(j,k-1,l) + 0.5_wp*mu_L(-1)*((47._wp/60._wp)*(4._wp*dvel1-2._wp*dvel2)/3._wp)*vel_L(2,-1)*(1._wp/dy(k-1))
                                !$acc atomic
                                rhs_vf(E_idx)%sf(j,k-2,l) = rhs_vf(E_idx)%sf(j,k-2,l) + 0.5_wp*mu_L(-2)*((-13._wp/60._wp)*(4._wp*dvel1-2._wp*dvel2)/3._wp)*vel_L(2,-2)*(1._wp/dy(k-2))
                                !$acc atomic
                                rhs_vf(E_idx)%sf(j,k-3,l) = rhs_vf(E_idx)%sf(j,k-3,l) + 0.5_wp*mu_L(-3)*((2._wp/60._wp)*(4._wp*dvel1-2._wp*dvel2)/3._wp)*vel_L(2,-3)*(1._wp/dy(k-3))
                               
                                !$acc atomic
                                rhs_vf(momxb+1)%sf(j,k-1,l) = rhs_vf(momxb+1)%sf(j,k-1,l) - 0.5_wp*mu_R(-2)*((-3._wp/60._wp)*(4._wp*dvel1-2._wp*dvel2)/3._wp)*(1._wp/dy(k-1))
                                !$acc atomic
                                rhs_vf(momxb+1)%sf(j,k,l) = rhs_vf(momxb+1)%sf(j,k,l) - 0.5_wp*mu_R(-1)*((27._wp/60._wp)*(4._wp*dvel1-2._wp*dvel2)/3._wp)*(1._wp/dy(k))
                                !$acc atomic
                                rhs_vf(momxb+1)%sf(j,k+1,l) = rhs_vf(momxb+1)%sf(j,k+1,l) - 0.5_wp*mu_R(0)*((47._wp/60._wp)*(4._wp*dvel1-2._wp*dvel2)/3._wp)*(1._wp/dy(k+1))
                                !$acc atomic
                                rhs_vf(momxb+1)%sf(j,k+2,l) = rhs_vf(momxb+1)%sf(j,k+2,l) - 0.5_wp*mu_R(1)*((-13._wp/60._wp)*(4._wp*dvel1-2._wp*dvel2)/3._wp)*(1._wp/dy(k+2))
                                !$acc atomic
                                rhs_vf(momxb+1)%sf(j,k+3,l) = rhs_vf(momxb+1)%sf(j,k+3,l) - 0.5_wp*mu_R(2)*((2._wp/60._wp)*(4._wp*dvel1-2._wp*dvel2)/3._wp)*(1._wp/dy(k+3))

                                !$acc atomic
                                rhs_vf(E_idx)%sf(j,k-1,l) = rhs_vf(E_idx)%sf(j,k-1,l) - 0.5_wp*mu_R(-2)*((-3._wp/60._wp)*(4._wp*dvel1-2._wp*dvel2)/3._wp)*vel_R(2,-2)*(1._wp/dy(k-1))
                                !$acc atomic
                                rhs_vf(E_idx)%sf(j,k,l) = rhs_vf(E_idx)%sf(j,k,l) - 0.5_wp*mu_R(-1)*((27._wp/60._wp)*(4._wp*dvel1-2._wp*dvel2)/3._wp)*vel_R(2,-1)*(1._wp/dy(k))
                                !$acc atomic
                                rhs_vf(E_idx)%sf(j,k+1,l) = rhs_vf(E_idx)%sf(j,k+1,l) - 0.5_wp*mu_R(0)*((47._wp/60._wp)*(4._wp*dvel1-2._wp*dvel2)/3._wp)*vel_R(2,0)*(1._wp/dy(k+1))
                                !$acc atomic
                                rhs_vf(E_idx)%sf(j,k+2,l) = rhs_vf(E_idx)%sf(j,k+2,l) - 0.5_wp*mu_R(1)*((-13._wp/60._wp)*(4._wp*dvel1-2._wp*dvel2)/3._wp)*vel_R(2,1)*(1._wp/dy(k+2))
                                !$acc atomic
                                rhs_vf(E_idx)%sf(j,k+3,l) = rhs_vf(E_idx)%sf(j,k+3,l) - 0.5_wp*mu_R(2)*((2._wp/60._wp)*(4._wp*dvel1-2._wp*dvel2)/3._wp)*vel_R(2,2)*(1._wp/dy(k+3)) 

                                !$acc atomic
                                rhs_vf(momxb+1)%sf(j,k-2,l) = rhs_vf(momxb+1)%sf(j,k-2,l) + 0.5_wp*mu_R(-2)*((-3._wp/60._wp)*(4._wp*dvel1-2._wp*dvel2)/3._wp)*(1._wp/dy(k-2))
                                !$acc atomic
                                rhs_vf(momxb+1)%sf(j,k-1,l) = rhs_vf(momxb+1)%sf(j,k-1,l) + 0.5_wp*mu_R(-1)*((27._wp/60._wp)*(4._wp*dvel1-2._wp*dvel2)/3._wp)*(1._wp/dy(k-1))
                                !$acc atomic
                                rhs_vf(momxb+1)%sf(j,k,l) = rhs_vf(momxb+1)%sf(j,k,l) + 0.5_wp*mu_R(0)*((47._wp/60._wp)*(4._wp*dvel1-2._wp*dvel2)/3._wp)*(1._wp/dy(k))
                                !$acc atomic
                                rhs_vf(momxb+1)%sf(j,k+1,l) = rhs_vf(momxb+1)%sf(j,k+1,l) + 0.5_wp*mu_R(1)*((-13._wp/60._wp)*(4._wp*dvel1-2._wp*dvel2)/3._wp)*(1._wp/dy(k+1))
                                !$acc atomic
                                rhs_vf(momxb+1)%sf(j,k+2,l) = rhs_vf(momxb+1)%sf(j,k+2,l) + 0.5_wp*mu_R(2)*((2._wp/60._wp)*(4._wp*dvel1-2._wp*dvel2)/3._wp)*(1._wp/dy(k+2))

                                !$acc atomic
                                rhs_vf(E_idx)%sf(j,k-2,l) = rhs_vf(E_idx)%sf(j,k-2,l) + 0.5_wp*mu_R(-2)*((-3._wp/60._wp)*(4._wp*dvel1-2._wp*dvel2)/3._wp)*vel_R(2,-2)*(1._wp/dy(k-2))
                                !$acc atomic
                                rhs_vf(E_idx)%sf(j,k-1,l) = rhs_vf(E_idx)%sf(j,k-1,l) + 0.5_wp*mu_R(-1)*((27._wp/60._wp)*(4._wp*dvel1-2._wp*dvel2)/3._wp)*vel_R(2,-1)*(1._wp/dy(k-1))
                                !$acc atomic
                                rhs_vf(E_idx)%sf(j,k,l) = rhs_vf(E_idx)%sf(j,k,l) + 0.5_wp*mu_R(0)*((47._wp/60._wp)*(4._wp*dvel1-2._wp*dvel2)/3._wp)*vel_R(2,0)*(1._wp/dy(k))
                                !$acc atomic
                                rhs_vf(E_idx)%sf(j,k+1,l) = rhs_vf(E_idx)%sf(j,k+1,l) + 0.5_wp*mu_R(1)*((-13._wp/60._wp)*(4._wp*dvel1-2._wp*dvel2)/3._wp)*vel_R(2,1)*(1._wp/dy(k+1))
                                !$acc atomic
                                rhs_vf(E_idx)%sf(j,k+2,l) = rhs_vf(E_idx)%sf(j,k+2,l) + 0.5_wp*mu_R(2)*((2._wp/60._wp)*(4._wp*dvel1-2._wp*dvel2)/3._wp)*vel_R(2,2)*(1._wp/dy(k+2))

                            end if

                            !!dwz
                            dvel2 = (1/(12._wp*dz(l))) * ( &
                            8._wp*q_prim_vf(momxb+2)%sf(j,k,l+1)/rho_sf(1,3) - &
                            8._wp*q_prim_vf(momxb+2)%sf(j,k,l-1)/rho_sf(-1,3) + &
                            q_prim_vf(momxb+2)%sf(j,k,l-2)/rho_sf(-2,3) - &
                            q_prim_vf(momxb+2)%sf(j,k,l+2)/rho_sf(2,3) )

                            if(viscous) then

                                !$acc atomic
                                rhs_vf(momxb+1)%sf(j,k+2,l) = rhs_vf(momxb+1)%sf(j,k+2,l) - 0.5_wp*mu_L(1)*((-3._wp/60._wp)*(-2._wp*dvel2)/3._wp)*(1._wp/dy(k+2))
                                !$acc atomic
                                rhs_vf(momxb+1)%sf(j,k+1,l) = rhs_vf(momxb+1)%sf(j,k+1,l) - 0.5_wp*mu_L(0)*((27._wp/60._wp)*(-2._wp*dvel2)/3._wp)*(1._wp/dy(k+1))
                                !$acc atomic
                                rhs_vf(momxb+1)%sf(j,k,l) = rhs_vf(momxb+1)%sf(j,k,l) - 0.5_wp*mu_L(-1)*((47._wp/60._wp)*(-2._wp*dvel2)/3._wp)*(1._wp/dy(k))
                                !$acc atomic
                                rhs_vf(momxb+1)%sf(j,k-1,l) = rhs_vf(momxb+1)%sf(j,k-1,l) - 0.5_wp*mu_L(-2)*((-13._wp/60._wp)*(-2._wp*dvel2)/3._wp)*(1._wp/dy(k-1))
                                !$acc atomic
                                rhs_vf(momxb+1)%sf(j,k-2,l) = rhs_vf(momxb+1)%sf(j,k-2,l) - 0.5_wp*mu_L(-3)*((2._wp/60._wp)*(-2._wp*dvel2)/3._wp)*(1._wp/dy(k-2))

                                !$acc atomic
                                rhs_vf(E_idx)%sf(j,k+2,l) = rhs_vf(E_idx)%sf(j,k+2,l) - 0.5_wp*mu_L(1)*((-3._wp/60._wp)*(-2._wp*dvel2)/3._wp)*vel_L(2,1)*(1._wp/dy(k+2))
                                !$acc atomic
                                rhs_vf(E_idx)%sf(j,k+1,l) = rhs_vf(E_idx)%sf(j,k+1,l) - 0.5_wp*mu_L(0)*((27._wp/60._wp)*(-2._wp*dvel2)/3._wp)*vel_L(2,0)*(1._wp/dy(k+1))
                                !$acc atomic
                                rhs_vf(E_idx)%sf(j,k,l) = rhs_vf(E_idx)%sf(j,k,l) - 0.5_wp*mu_L(-1)*((47._wp/60._wp)*(-2._wp*dvel2)/3._wp)*vel_L(2,-1)*(1._wp/dy(k))
                                !$acc atomic
                                rhs_vf(E_idx)%sf(j,k-1,l) = rhs_vf(E_idx)%sf(j,k-1,l) - 0.5_wp*mu_L(-2)*((-13._wp/60._wp)*(-2._wp*dvel2)/3._wp)*vel_L(2,-2)*(1._wp/dy(k-1))
                                !$acc atomic
                                rhs_vf(E_idx)%sf(j,k-2,l) = rhs_vf(E_idx)%sf(j,k-2,l) - 0.5_wp*mu_L(-3)*((2._wp/60._wp)*(-2._wp*dvel2)/3._wp)*vel_L(2,-3)*(1._wp/dy(k-2))

                                !$acc atomic
                                rhs_vf(momxb+1)%sf(j,k+1,l) = rhs_vf(momxb+1)%sf(j,k+1,l) + 0.5_wp*mu_L(1)*((-3._wp/60._wp)*(-2._wp*dvel2)/3._wp)*(1._wp/dy(k+1))
                                !$acc atomic
                                rhs_vf(momxb+1)%sf(j,k,l) = rhs_vf(momxb+1)%sf(j,k,l) + 0.5_wp*mu_L(0)*((27._wp/60._wp)*(-2._wp*dvel2)/3._wp)*(1._wp/dy(k))
                                !$acc atomic
                                rhs_vf(momxb+1)%sf(j,k-1,l) = rhs_vf(momxb+1)%sf(j,k-1,l) + 0.5_wp*mu_L(-1)*((47._wp/60._wp)*(-2._wp*dvel2)/3._wp)*(1._wp/dy(k-1))
                                !$acc atomic
                                rhs_vf(momxb+1)%sf(j,k-2,l) = rhs_vf(momxb+1)%sf(j,k-2,l) + 0.5_wp*mu_L(-2)*((-13._wp/60._wp)*(-2._wp*dvel2)/3._wp)*(1._wp/dy(k-2))
                                !$acc atomic
                                rhs_vf(momxb+1)%sf(j,k-3,l) = rhs_vf(momxb+1)%sf(j,k-3,l) + 0.5_wp*mu_L(-3)*((2._wp/60._wp)*(-2._wp*dvel2)/3._wp)*(1._wp/dy(k-3))

                                !$acc atomic
                                rhs_vf(E_idx)%sf(j,k+1,l) = rhs_vf(E_idx)%sf(j,k+1,l) + 0.5_wp*mu_L(1)*((-3._wp/60._wp)*(-2._wp*dvel2)/3._wp)*vel_L(2,1)*(1._wp/dy(k+1))
                                !$acc atomic
                                rhs_vf(E_idx)%sf(j,k,l) = rhs_vf(E_idx)%sf(j,k,l) + 0.5_wp*mu_L(0)*((27._wp/60._wp)*(-2._wp*dvel2)/3._wp)*vel_L(2,0)*(1._wp/dy(k))
                                !$acc atomic
                                rhs_vf(E_idx)%sf(j,k-1,l) = rhs_vf(E_idx)%sf(j,k-1,l) + 0.5_wp*mu_L(-1)*((47._wp/60._wp)*(-2._wp*dvel2)/3._wp)*vel_L(2,-1)*(1._wp/dy(k-1))
                                !$acc atomic
                                rhs_vf(E_idx)%sf(j,k-2,l) = rhs_vf(E_idx)%sf(j,k-2,l) + 0.5_wp*mu_L(-2)*((-13._wp/60._wp)*(-2._wp*dvel2)/3._wp)*vel_L(2,-2)*(1._wp/dy(k-2))
                                !$acc atomic
                                rhs_vf(E_idx)%sf(j,k-3,l) = rhs_vf(E_idx)%sf(j,k-3,l) + 0.5_wp*mu_L(-3)*((2._wp/60._wp)*(-2._wp*dvel2)/3._wp)*vel_L(2,-3)*(1._wp/dy(k-3))
                               
                                !$acc atomic
                                rhs_vf(momxb+1)%sf(j,k-1,l) = rhs_vf(momxb+1)%sf(j,k-1,l) - 0.5_wp*mu_R(-2)*((-3._wp/60._wp)*(-2._wp*dvel2)/3._wp)*(1._wp/dy(k-1))
                                !$acc atomic
                                rhs_vf(momxb+1)%sf(j,k,l) = rhs_vf(momxb+1)%sf(j,k,l) - 0.5_wp*mu_R(-1)*((27._wp/60._wp)*(-2._wp*dvel2)/3._wp)*(1._wp/dy(k))
                                !$acc atomic
                                rhs_vf(momxb+1)%sf(j,k+1,l) = rhs_vf(momxb+1)%sf(j,k+1,l) - 0.5_wp*mu_R(0)*((47._wp/60._wp)*(-2._wp*dvel2)/3._wp)*(1._wp/dy(k+1))
                                !$acc atomic
                                rhs_vf(momxb+1)%sf(j,k+2,l) = rhs_vf(momxb+1)%sf(j,k+2,l) - 0.5_wp*mu_R(1)*((-13._wp/60._wp)*(-2._wp*dvel2)/3._wp)*(1._wp/dy(k+2))
                                !$acc atomic
                                rhs_vf(momxb+1)%sf(j,k+3,l) = rhs_vf(momxb+1)%sf(j,k+3,l) - 0.5_wp*mu_R(2)*((2._wp/60._wp)*(-2._wp*dvel2)/3._wp)*(1._wp/dy(k+3))

                                !$acc atomic
                                rhs_vf(E_idx)%sf(j,k-1,l) = rhs_vf(E_idx)%sf(j,k-1,l) - 0.5_wp*mu_R(-2)*((-3._wp/60._wp)*(-2._wp*dvel2)/3._wp)*vel_R(2,-2)*(1._wp/dy(k-1))
                                !$acc atomic
                                rhs_vf(E_idx)%sf(j,k,l) = rhs_vf(E_idx)%sf(j,k,l) - 0.5_wp*mu_R(-1)*((27._wp/60._wp)*(-2._wp*dvel2)/3._wp)*vel_R(2,-1)*(1._wp/dy(k))
                                !$acc atomic
                                rhs_vf(E_idx)%sf(j,k+1,l) = rhs_vf(E_idx)%sf(j,k+1,l) - 0.5_wp*mu_R(0)*((47._wp/60._wp)*(-2._wp*dvel2)/3._wp)*vel_R(2,0)*(1._wp/dy(k+1))
                                !$acc atomic
                                rhs_vf(E_idx)%sf(j,k+2,l) = rhs_vf(E_idx)%sf(j,k+2,l) - 0.5_wp*mu_R(1)*((-13._wp/60._wp)*(-2._wp*dvel2)/3._wp)*vel_R(2,1)*(1._wp/dy(k+2))
                                !$acc atomic
                                rhs_vf(E_idx)%sf(j,k+3,l) = rhs_vf(E_idx)%sf(j,k+3,l) - 0.5_wp*mu_R(2)*((2._wp/60._wp)*(-2._wp*dvel2)/3._wp)*vel_R(2,2)*(1._wp/dy(k+3)) 

                                !$acc atomic
                                rhs_vf(momxb+1)%sf(j,k-2,l) = rhs_vf(momxb+1)%sf(j,k-2,l) + 0.5_wp*mu_R(-2)*((-3._wp/60._wp)*(-2._wp*dvel2)/3._wp)*(1._wp/dy(k-2))
                                !$acc atomic
                                rhs_vf(momxb+1)%sf(j,k-1,l) = rhs_vf(momxb+1)%sf(j,k-1,l) + 0.5_wp*mu_R(-1)*((27._wp/60._wp)*(-2._wp*dvel2)/3._wp)*(1._wp/dy(k-1))
                                !$acc atomic
                                rhs_vf(momxb+1)%sf(j,k,l) = rhs_vf(momxb+1)%sf(j,k,l) + 0.5_wp*mu_R(0)*((47._wp/60._wp)*(-2._wp*dvel2)/3._wp)*(1._wp/dy(k))
                                !$acc atomic
                                rhs_vf(momxb+1)%sf(j,k+1,l) = rhs_vf(momxb+1)%sf(j,k+1,l) + 0.5_wp*mu_R(1)*((-13._wp/60._wp)*(-2._wp*dvel2)/3._wp)*(1._wp/dy(k+1))
                                !$acc atomic
                                rhs_vf(momxb+1)%sf(j,k+2,l) = rhs_vf(momxb+1)%sf(j,k+2,l) + 0.5_wp*mu_R(2)*((2._wp/60._wp)*(-2._wp*dvel2)/3._wp)*(1._wp/dy(k+2))

                                !$acc atomic
                                rhs_vf(E_idx)%sf(j,k-2,l) = rhs_vf(E_idx)%sf(j,k-2,l) + 0.5_wp*mu_R(-2)*((-3._wp/60._wp)*(-2._wp*dvel2)/3._wp)*vel_R(2,-2)*(1._wp/dy(k-2))
                                !$acc atomic
                                rhs_vf(E_idx)%sf(j,k-1,l) = rhs_vf(E_idx)%sf(j,k-1,l) + 0.5_wp*mu_R(-1)*((27._wp/60._wp)*(-2._wp*dvel2)/3._wp)*vel_R(2,-1)*(1._wp/dy(k-1))
                                !$acc atomic
                                rhs_vf(E_idx)%sf(j,k,l) = rhs_vf(E_idx)%sf(j,k,l) + 0.5_wp*mu_R(0)*((47._wp/60._wp)*(-2._wp*dvel2)/3._wp)*vel_R(2,0)*(1._wp/dy(k))
                                !$acc atomic
                                rhs_vf(E_idx)%sf(j,k+1,l) = rhs_vf(E_idx)%sf(j,k+1,l) + 0.5_wp*mu_R(1)*((-13._wp/60._wp)*(-2._wp*dvel2)/3._wp)*vel_R(2,1)*(1._wp/dy(k+1))
                                !$acc atomic
                                rhs_vf(E_idx)%sf(j,k+2,l) = rhs_vf(E_idx)%sf(j,k+2,l) + 0.5_wp*mu_R(2)*((2._wp/60._wp)*(-2._wp*dvel2)/3._wp)*vel_R(2,2)*(1._wp/dy(k+2))

                            end if


                            F_L = (1._wp/60._wp) * (-3._wp * jac(j, k+2, l) + &
                                                27._wp * jac(j, k+1, l) + &
                                                47._wp * jac(j, k, l) -   &
                                                13._wp * jac(j, k-1, l) + &
                                                2._wp * jac(j, k-2, l))

                            !$acc loop seq
                            do i = 1, num_fluids
                                !$acc atomic
                                rhs_vf(i)%sf(j,k+1,l) = rhs_vf(i)%sf(j,k+1,l) + &
                                    (0.5_wp * (alpha_rho_R(i) * &
                                    vel_R(2,0))*(1._wp/dy(k+1)) + &
                                    0.5_wp*cfl * (alpha_rho_R(i))*(1._wp/dy(k+1)))

                                !$acc atomic
                                rhs_vf(i)%sf(j,k,l) = rhs_vf(i)%sf(j,k,l) - &
                                    (0.5_wp * (alpha_rho_R(i) * &
                                    vel_R(2,0))*(1._wp/dy(k)) + &
                                    0.5_wp*cfl * (alpha_rho_R(i))*(1._wp/dy(k))) 
                            end do

                            if(num_fluids > 1) then 
                                !$acc loop seq
                                do i = 1, num_fluids
                                    !$acc atomic
                                    rhs_vf(advxb+i-1)%sf(j,k+1,l) = rhs_vf(advxb+i-1)%sf(j,k+1,l) + &
                                        (0.5_wp * (alpha_R(i) * &
                                        vel_R(2,0))*(1._wp/dy(k+1)) + &
                                        0.5_wp*cfl*(alpha_R(i))*(1._wp/dy(k+1)))

                                    !$acc atomic
                                    rhs_vf(advxb+i-1)%sf(j,k+1,l) = rhs_vf(advxb+i-1)%sf(j,k+1,l) &
                                    - (0.5_wp * q_prim_vf(advxb+i-1)%sf(j,k+1,l)  * vel_R(2,0)*(1._wp/dy(k+1)))

                                    !$acc atomic
                                    rhs_vf(advxb+i-1)%sf(j,k,l) = rhs_vf(advxb+i-1)%sf(j,k,l) - &
                                        (0.5_wp * (alpha_R(i) * &
                                        vel_R(2,0))*(1._wp/dy(k)) + &
                                        0.5_wp*cfl*(alpha_R(i))*(1._wp/dy(k)))

                                    !$acc atomic
                                    rhs_vf(advxb+i-1)%sf(j,k,l) = rhs_vf(advxb+i-1)%sf(j,k,l) &
                                    + (0.5_wp * q_prim_vf(advxb+i-1)%sf(j,k,l)  * vel_R(2,0)*(1._wp/dy(k)))
                                end do
                            end if

                            !$acc atomic
                            rhs_vf(momxb+1)%sf(j,k+1,l) = rhs_vf(momxb+1)%sf(j,k+1,l) + &
                                 (0.5_wp* (rho_R * (vel_R(2,0))**2.0 + &
                                 pres_R+F_L)*(1._wp/dy(k+1)) + &
                                 0.5_wp*cfl * (rho_R*vel_R(2,0))*(1._wp/dy(k+1)))

                            !$acc atomic
                            rhs_vf(momxb)%sf(j, k+1, l) =  rhs_vf(momxb)%sf(j, k+1, l) + &
                                (0.5_wp * rho_R * vel_R(2,0)*vel_R(1,0)*(1._wp/dy(k+1)) + &
                                0.5_wp*cfl * (rho_R*vel_R(1,0))*(1._wp/dy(k+1)))

                            !$acc atomic
                            rhs_vf(momxb+2)%sf(j, k+1, l) =  rhs_vf(momxb+2)%sf(j, k+1, l) + &
                                (0.5_wp * rho_R * vel_R(2,0)*vel_R(3,0)*(1._wp/dy(k+1)) + &
                                0.5_wp*cfl * (rho_R*vel_R(3,0))*(1._wp/dy(k+1)))

                            !$acc atomic
                             rhs_vf(E_idx)%sf(j, k+1, l) = rhs_vf(E_idx)%sf(j, k+1, l) +  &
                                (0.5_wp * (vel_R(2,0) * (E_R+ &
                                pres_R+F_L ) )*(1._wp/dy(k+1)) + &
                                0.5_wp*cfl * (E_R)*(1._wp/dy(k+1)))

                            !$acc atomic
                            rhs_vf(momxb+1)%sf(j,k,l) = rhs_vf(momxb+1)%sf(j,k,l) - &
                                 (0.5_wp* (rho_R * (vel_R(2,0))**2.0 + &
                                 pres_R+F_L)*(1._wp/dy(k)) + &
                                 0.5_wp*cfl * (rho_R*vel_R(2,0))*(1._wp/dy(k)))

                            !$acc atomic
                            rhs_vf(momxb)%sf(j, k, l) =  rhs_vf(momxb)%sf(j, k, l) - &
                                (0.5_wp * rho_R * vel_R(2,0)*vel_R(1,0)*(1._wp/dy(k)) + &
                                0.5_wp*cfl * (rho_R*vel_R(1,0))*(1._wp/dy(k)))

                            !$acc atomic
                            rhs_vf(momxb+2)%sf(j, k, l) =  rhs_vf(momxb+2)%sf(j, k, l) - &
                                (0.5_wp * rho_R * vel_R(2,0)*vel_R(3,0)*(1._wp/dy(k)) + &
                                0.5_wp*cfl * (rho_R*vel_R(3,0))*(1._wp/dy(k)))

                            !$acc atomic
                             rhs_vf(E_idx)%sf(j, k, l) = rhs_vf(E_idx)%sf(j, k, l) -  &
                                (0.5_wp * (vel_R(2,0) * (E_R + &
                                pres_R+F_L) )*(1._wp/dy(k)) + &
                                0.5_wp*cfl * (E_R)*(1._wp/dy(k)))

                        end do
                    end do
                end do
            end if
        elseif (idir == 3) then
                !$acc parallel loop collapse(3) gang vector default(present) private(rho_L,gamma_L,pi_inf_L,mu_L,a_L,vel_L,vel_R, pres_L,alpha_L,alpha_R,alpha_rho_L,cfl,dvel1,dvel2,F_L,E_L,mu_R,rho_sf,alpha_rho_R)
                do l = -buff_size+5, p+buff_size-5
                    do k = 0, n
                        do j = 0, m

                            !$acc loop seq 
                            do i = 1, num_fluids
                                alpha_rho_L(i) = (1._wp/60._wp) * (-3._wp * q_prim_vf(i)%sf(j, k, l-1) + &
                                                27._wp * q_prim_vf(i)%sf(j, k, l) + &
                                                47._wp * q_prim_vf(i)%sf(j, k, l+1) -   &
                                                13._wp * q_prim_vf(i)%sf(j, k, l+2) + &
                                                2._wp * q_prim_vf(i)%sf(j, k, l+3))

                                alpha_rho_R(i) = (1._wp/60._wp) * (-3._wp * q_prim_vf(i)%sf(j, k, l+2) + &
                                                27._wp * q_prim_vf(i)%sf(j, k, l+1) + &
                                                47._wp * q_prim_vf(i)%sf(j, k, l) -   &
                                                13._wp * q_prim_vf(i)%sf(j, k, l-1) + &
                                                2._wp * q_prim_vf(i)%sf(j, k, l-2))

                                alpha_L(i) =  (1._wp/60._wp) * (-3._wp * q_prim_vf(E_idx+i)%sf(j, k, l-1) + &
                                                27._wp * q_prim_vf(E_idx+i)%sf(j, k, l) + &
                                                47._wp * q_prim_vf(E_idx+i)%sf(j, k, l+1) -   &
                                                13._wp * q_prim_vf(E_idx+i)%sf(j, k, l+2) + &
                                                2._wp * q_prim_vf(E_idx+i)%sf(j, k, l+3))
                            end do 

                            if(num_fluids > 1) then  
                                !$acc loop seq 
                                do i = 1,num_fluids                           
                                    alpha_R(i) = (1._wp/60._wp) * (-3._wp * q_prim_vf(E_idx+i)%sf(j, k, l+2) + &
                                                        27._wp * q_prim_vf(E_idx+i)%sf(j, k, l+1) + &
                                                        47._wp * q_prim_vf(E_idx+i)%sf(j, k, l) -   &
                                                        13._wp * q_prim_vf(E_idx+i)%sf(j, k, l-1) + &
                                                        2._wp * q_prim_vf(E_idx+i)%sf(j, k, l-2))
                                end do 
                            else 
                                alpha_R(1) = 1._wp
                            end if

                            F_L = (1._wp/60._wp) * (-3._wp * jac(j, k, l-1) + &
                                                27._wp * jac(j, k, l) + &
                                                47._wp * jac(j, k, l+1) -   &
                                                13._wp * jac(j, k, l+2) + &
                                                2._wp * jac(j, k, l+3))
                            !$acc loop seq 
                            do i = -2, 2
                                rho_L = 0._wp
                                !$acc loop seq 
                                do q = 1, num_fluids
                                    rho_L = rho_L + q_prim_vf(q)%sf(j+i,k,l)
                                end do
                                rho_sf(i,1) = rho_L
                            end do

                            !$acc loop seq 
                            do i = -2, 2
                                rho_L = 0._wp
                                !$acc loop seq 
                                do q = 1, num_fluids
                                    rho_L = rho_L + q_prim_vf(q)%sf(j,k+i,l)
                                end do
                                rho_sf(i,2) = rho_L
                            end do

                            !$acc loop seq 
                            do i = -2, 2
                                rho_L = 0._wp
                                !$acc loop seq 
                                do q = 1, num_fluids
                                    rho_L = rho_L + q_prim_vf(q)%sf(j,k,l+i)
                                end do
                                rho_sf(i,3) = rho_L
                            end do

                            rho_L = sum(alpha_rho_L)
                            gamma_L = sum(alpha_L*gammas)
                            pi_inf_L = sum(alpha_L*pi_infs)

                            rho_R = sum(alpha_rho_R)
                            gamma_R = sum(alpha_R*gammas)
                            pi_inf_R = sum(alpha_R*pi_infs)

                            !$acc loop seq 
                            do q = -3, 2
                                !$acc loop seq 
                                do i = 1, num_dims
                                    vel_L(i,q) = (1._wp/60._wp) * (-3._wp * q_prim_vf(momxb+i-1)%sf(j,k,l-1+q) + &
                                                    27._wp * q_prim_vf(momxb+i-1)%sf(j,k,l+q) + &
                                                    47._wp * q_prim_vf(momxb+i-1)%sf(j,k,l+1+q) -   &
                                                    13._wp * q_prim_vf(momxb+i-1)%sf(j,k,l+2+q) + &
                                                    2._wp * q_prim_vf(momxb+i-1)%sf(j,k,l+3+q)) / rho_L
                                end do 

                                !$acc loop seq 
                                do i = 1, num_dims
                                    vel_R(i,q) = (1._wp/60._wp) * (-3._wp * q_prim_vf(momxb+i-1)%sf(j,k,l+2+q) + &
                                                    27._wp * q_prim_vf(momxb+i-1)%sf(j,k,l+1+q) + &
                                                    47._wp * q_prim_vf(momxb+i-1)%sf(j,k,l+q) -   &
                                                    13._wp * q_prim_vf(momxb+i-1)%sf(j,k,l-1+q) + &
                                                    2._wp * q_prim_vf(momxb+i-1)%sf(j,k,l-2+q)) / rho_R
                                end do
                            end do

                            E_L = (1._wp/60._wp) * (-3._wp * q_prim_vf(E_idx)%sf(j,k,l-1) + &
                                                27._wp * q_prim_vf(E_idx)%sf(j,k,l) + &
                                                47._wp * q_prim_vf(E_idx)%sf(j,k,l+1) -   &
                                                13._wp * q_prim_vf(E_idx)%sf(j,k,l+2) + &
                                                2._wp * q_prim_vf(E_idx)%sf(j,k,l+3))

                            E_R = (1._wp/60._wp) * (-3._wp * q_prim_vf(E_idx)%sf(j,k,l+2) + &
                                                27._wp * q_prim_vf(E_idx)%sf(j,k,l+1) + &
                                                47._wp * q_prim_vf(E_idx)%sf(j,k,l) -   &
                                                13._wp * q_prim_vf(E_idx)%sf(j,k,l-1) + &
                                                2._wp * q_prim_vf(E_idx)%sf(j,k,l-2))       

                            pres_L = (E_L - pi_inf_L - 0.5_wp*rho_L*(vel_L(1,0)**2._wp + vel_L(2,0)**2._wp + vel_L(3,0)**2._wp ))/gamma_L                    

                            pres_R = (E_R - pi_inf_R - 0.5_wp*rho_R*(vel_R(1,0)**2._wp + vel_R(2,0)**2._wp + vel_R(3,0)**2._wp ))/gamma_R 

                            a_L = sqrt((pres_L*(1._wp/gamma_L + 1._wp) + pi_inf_L / gamma_L) / rho_L)
                            a_R = sqrt((pres_R*(1._wp/gamma_R + 1._wp) + pi_inf_R / gamma_R) / rho_R)

                            cfl = (max(sqrt(vel_L(1,0)**2._wp + vel_L(2,0)**2._wp + vel_L(3,0)**2._wp),sqrt(vel_R(1,0)**2._wp + vel_R(2,0)**2._wp + vel_R(3,0)**2._wp)) + max(a_L,a_R))


                            !$acc loop seq
                            do i = 1, num_fluids
                                !$acc atomic
                                rhs_vf(i)%sf(j,k,l+1) = rhs_vf(i)%sf(j,k,l+1) + &
                                    (0.5_wp * (alpha_rho_L(i) * &
                                    vel_L(3,0))*(1._wp/dz(l+1)) - &
                                    0.5_wp*cfl * (alpha_rho_L(i))*(1._wp/dz(l+1)))

                                !$acc atomic
                                rhs_vf(i)%sf(j,k,l) = rhs_vf(i)%sf(j,k,l) - &
                                    (0.5_wp * (alpha_rho_L(i) * &
                                    vel_L(3,0))*(1._wp/dz(l)) - &
                                    0.5_wp*cfl * (alpha_rho_L(i))*(1._wp/dz(l))) 
                            end do

                            if(num_fluids > 1) then 
                                !$acc loop seq
                                do i = 1, num_fluids
                                    !$acc atomic
                                    rhs_vf(advxb+i-1)%sf(j,k,l+1) = rhs_vf(advxb+i-1)%sf(j,k,l+1) + &
                                        (0.5_wp * (alpha_L(i) * &
                                        vel_L(3,0))*(1._wp/dz(l+1)) - &
                                        0.5_wp*cfl*(alpha_L(i))*(1._wp/dz(l+1)))

                                    !$acc atomic
                                    rhs_vf(advxb+i-1)%sf(j,k,l+1) = rhs_vf(advxb+i-1)%sf(j,k,l+1) &
                                    - (0.5_wp * q_prim_vf(advxb+i-1)%sf(j,k,l+1)  * vel_L(3,0)*(1._wp/dz(l+1)))

                                    !$acc atomic
                                    rhs_vf(advxb+i-1)%sf(j,k,l) = rhs_vf(advxb+i-1)%sf(j,k,l) - &
                                        (0.5_wp * (alpha_L(i) * &
                                        vel_L(3,0))*(1._wp/dz(l)) - &
                                        0.5_wp*cfl*(alpha_L(i))*(1._wp/dz(l)))

                                    !$acc atomic
                                    rhs_vf(advxb+i-1)%sf(j,k,l) = rhs_vf(advxb+i-1)%sf(j,k,l) &
                                    + (0.5_wp * q_prim_vf(advxb+i-1)%sf(j,k,l)  * vel_L(3,0)*(1._wp/dz(l)))
                                end do
                            end if

                            !$acc atomic
                            rhs_vf(momxb+2)%sf(j,k,l+1) = rhs_vf(momxb+2)%sf(j,k,l+1) + &
                                 (0.5_wp* (rho_L * (vel_L(3,0))**2.0 + &
                                 pres_L+F_L)*(1._wp/dz(l+1)) - &
                                 0.5_wp*cfl * (rho_L*vel_L(3,0))*(1._wp/dz(l+1)))

                            !$acc atomic
                            rhs_vf(momxb)%sf(j, k,l+1) =  rhs_vf(momxb)%sf(j, k,l+1) + &
                                (0.5_wp * rho_L * vel_L(1,0)*vel_L(3,0)*(1._wp/dz(l+1)) - &
                                0.5_wp*cfl * (rho_L*vel_L(1,0))*(1._wp/dz(l+1)))

                            !$acc atomic
                            rhs_vf(momxb+1)%sf(j,k,l+1) =  rhs_vf(momxb+1)%sf(j,k,l+1) + &
                                (0.5_wp * rho_L * vel_L(2,0)*vel_L(3,0)*(1._wp/dz(l+1)) - &
                                0.5_wp*cfl * (rho_L*vel_L(2,0))*(1._wp/dz(l+1)))

                            !$acc atomic
                            rhs_vf(E_idx)%sf(j,k,l+1) = rhs_vf(E_idx)%sf(j,k,l+1) +  &
                                (0.5_wp * (vel_L(3,0) * (E_L + &
                                pres_L+F_L) )*(1._wp/dz(l+1)) - &
                                0.5_wp*cfl * (E_L)*(1._wp/dz(l+1)))

                            !$acc atomic
                            rhs_vf(momxb+2)%sf(j,k,l) = rhs_vf(momxb+2)%sf(j,k,l) - &
                                 (0.5_wp* (rho_L * (vel_L(3,0))**2.0 + &
                                 pres_L+F_L)*(1._wp/dz(l)) - &
                                 0.5_wp*cfl * (rho_L*vel_L(3,0))*(1._wp/dz(l)))

                            !$acc atomic
                            rhs_vf(momxb)%sf(j, k, l) =  rhs_vf(momxb)%sf(j, k, l) - &
                                (0.5_wp * rho_L * vel_L(1,0)*vel_L(3,0)*(1._wp/dz(l)) - &
                                0.5_wp*cfl * (rho_L*vel_L(1,0))*(1._wp/dz(l)))

                            !$acc atomic
                            rhs_vf(momxb+1)%sf(j, k, l) =  rhs_vf(momxb+1)%sf(j, k, l) - &
                                (0.5_wp * rho_L * vel_L(2,0)*vel_L(3,0)*(1._wp/dz(l)) - &
                                0.5_wp*cfl * (rho_L*vel_L(2,0))*(1._wp/dz(l)))

                            !$acc atomic
                            rhs_vf(E_idx)%sf(j, k, l) = rhs_vf(E_idx)%sf(j, k, l) -  &
                                (0.5_wp * (vel_L(3,0) * (E_L + &
                                pres_L+F_L) )*(1._wp/dz(l)) - &
                                0.5_wp*cfl * (E_L)*(1._wp/dz(l)))

                            !! dwx & duz
                            dvel1 = (1/(12._wp*dz(l))) * ( &
                            8._wp*q_prim_vf(momxb)%sf(j,k,l+1)/rho_sf(1,3) - &
                            8._wp*q_prim_vf(momxb)%sf(j,k,l-1)/rho_sf(-1,3)  + &
                            q_prim_vf(momxb)%sf(j,k,l-2)/rho_sf(-2,3)  - &
                            q_prim_vf(momxb)%sf(j,k,l+2)/rho_sf(2,3)  )


                            dvel2 = (1/(12._wp*dx(j))) * ( &
                            8._wp*q_prim_vf(momxb+2)%sf(j+1,k,l)/rho_sf(1,1) - &
                            8._wp*q_prim_vf(momxb+2)%sf(j-1,k,l)/rho_sf(-1,1) + &
                            q_prim_vf(momxb+2)%sf(j-2,k,l)/rho_sf(-2,1) - &
                            q_prim_vf(momxb+2)%sf(j+2,k,l)/rho_sf(2,1) )

                            if(num_fluids > 1) then 
                                !$acc loop seq
                                do q = -3, 2 
                                    mu_L(q) = 0._wp
                                    mu_R(q) = 0._wp
                                    !$acc loop seq 
                                    do i = 1, num_fluids
                                        mu_L(q) = (1._wp/60._wp) * (-3._wp * q_prim_vf(E_idx+i)%sf(j, k, l-1+q) + &
                                                27._wp * q_prim_vf(E_idx+i)%sf(j, k, l+q) + &
                                                47._wp * q_prim_vf(E_idx+i)%sf(j, k, l+1+q) -   &
                                                13._wp * q_prim_vf(E_idx+i)%sf(j, k, l+2+q) + &
                                                2._wp * q_prim_vf(E_idx+i)%sf(j, k, l+3+q)) / Res(1, i) + mu_L(q)
                                        
                                        mu_R(q) = (1._wp/60._wp) * (-3._wp * q_prim_vf(E_idx+i)%sf(j, k, l+2+q) + &
                                                27._wp * q_prim_vf(E_idx+i)%sf(j, k, l+1+q) + &
                                                47._wp * q_prim_vf(E_idx+i)%sf(j, k, l+q) -   &
                                                13._wp * q_prim_vf(E_idx+i)%sf(j, k, l-1+q) + &
                                                2._wp * q_prim_vf(E_idx+i)%sf(j, k, l-2+q)) / Res(1, i) + mu_R(q)
                                    end do
                                end do
                            else
                                !$acc loop seq
                                do q = -3, 2
                                    mu_L(q) = 1._wp / Res(1, 1) 
                                    mu_R(q) = 1._wp / Res(1, 1) 
                                end do
                            end if
                            if(viscous) then

                                !$acc atomic
                                rhs_vf(momxb)%sf(j,k,l+2) = rhs_vf(momxb)%sf(j,k,l+2) - 0.5_wp*mu_L(1)*((-3._wp/60._wp)*(dvel1+dvel2))*(1._wp/dz(l+2))
                                !$acc atomic
                                rhs_vf(momxb)%sf(j,k,l+1) = rhs_vf(momxb)%sf(j,k,l+1) - 0.5_wp*mu_L(0)*((27._wp/60._wp)*(dvel1+dvel2))*(1._wp/dz(l+1))
                                !$acc atomic
                                rhs_vf(momxb)%sf(j,k,l) = rhs_vf(momxb)%sf(j,k,l) - 0.5_wp*mu_L(-1)*((47._wp/60._wp)*(dvel1+dvel2))*(1._wp/dz(l))
                                !$acc atomic
                                rhs_vf(momxb)%sf(j,k,l-1) = rhs_vf(momxb)%sf(j,k,l-1) - 0.5_wp*mu_L(-2)*((-13._wp/60._wp)*(dvel1+dvel2))*(1._wp/dz(l-1))
                                !$acc atomic
                                rhs_vf(momxb)%sf(j,k,l-2) = rhs_vf(momxb)%sf(j,k,l-2) - 0.5_wp*mu_L(-3)*((2._wp/60._wp)*(dvel1+dvel2))*(1._wp/dz(l-2))

                                !$acc atomic
                                rhs_vf(E_idx)%sf(j,k,l+2) = rhs_vf(E_idx)%sf(j,k,l+2) - 0.5_wp*mu_L(1)*((-3._wp/60._wp)*(dvel1+dvel2))*vel_L(1,1)*(1._wp/dz(l+2))
                                !$acc atomic
                                rhs_vf(E_idx)%sf(j,k,l+1) = rhs_vf(E_idx)%sf(j,k,l+1) - 0.5_wp*mu_L(0)*((27._wp/60._wp)*(dvel1+dvel2))*vel_L(1,0)*(1._wp/dz(l+1))
                                !$acc atomic
                                rhs_vf(E_idx)%sf(j,k,l) = rhs_vf(E_idx)%sf(j,k,l) - 0.5_wp*mu_L(-1)*((47._wp/60._wp)*(dvel1+dvel2))*vel_L(1,-1)*(1._wp/dz(l))
                                !$acc atomic
                                rhs_vf(E_idx)%sf(j,k,l-1) = rhs_vf(E_idx)%sf(j,k,l-1) - 0.5_wp*mu_L(-2)*((-13._wp/60._wp)*(dvel1+dvel2))*vel_L(1,-2)*(1._wp/dz(l-1))
                                !$acc atomic
                                rhs_vf(E_idx)%sf(j,k,l-2) = rhs_vf(E_idx)%sf(j,k,l-2) - 0.5_wp*mu_L(-3)*((2._wp/60._wp)*(dvel1+dvel2))*vel_L(1,-3)*(1._wp/dz(l-2))

                                !$acc atomic
                                rhs_vf(momxb)%sf(j,k,l+1) = rhs_vf(momxb)%sf(j,k,l+1) + 0.5_wp*mu_L(1)*((-3._wp/60._wp)*(dvel1+dvel2))*(1._wp/dz(l+1))
                                !$acc atomic
                                rhs_vf(momxb)%sf(j,k,l) = rhs_vf(momxb)%sf(j,k,l) + 0.5_wp*mu_L(0)*((27._wp/60._wp)*(dvel1+dvel2))*(1._wp/dz(l))
                                !$acc atomic
                                rhs_vf(momxb)%sf(j,k,l-1) = rhs_vf(momxb)%sf(j,k,l-1) + 0.5_wp*mu_L(-1)*((47._wp/60._wp)*(dvel1+dvel2))*(1._wp/dz(l-1))
                                !$acc atomic
                                rhs_vf(momxb)%sf(j,k,l-2) = rhs_vf(momxb)%sf(j,k,l-2) + 0.5_wp*mu_L(-2)*((-13._wp/60._wp)*(dvel1+dvel2))*(1._wp/dz(l-2))
                                !$acc atomic
                                rhs_vf(momxb)%sf(j,k,l-3) = rhs_vf(momxb)%sf(j,k,l-3) + 0.5_wp*mu_L(-3)*((2._wp/60._wp)*(dvel1+dvel2))*(1._wp/dz(l-3))

                                !$acc atomic
                                rhs_vf(E_idx)%sf(j,k,l+1) = rhs_vf(E_idx)%sf(j,k,l+1) + 0.5_wp*mu_L(1)*((-3._wp/60._wp)*(dvel1+dvel2))*vel_L(1,1)*(1._wp/dz(l+1))
                                !$acc atomic
                                rhs_vf(E_idx)%sf(j,k,l) = rhs_vf(E_idx)%sf(j,k,l) + 0.5_wp*mu_L(0)*((27._wp/60._wp)*(dvel1+dvel2))*vel_L(1,0)*(1._wp/dz(l))
                                !$acc atomic
                                rhs_vf(E_idx)%sf(j,k,l-1) = rhs_vf(E_idx)%sf(j,k,l-1) + 0.5_wp*mu_L(-1)*((47._wp/60._wp)*(dvel1+dvel2))*vel_L(1,-1)*(1._wp/dz(l-1))
                                !$acc atomic
                                rhs_vf(E_idx)%sf(j,k,l-2) = rhs_vf(E_idx)%sf(j,k,l-2) + 0.5_wp*mu_L(-2)*((-13._wp/60._wp)*(dvel1+dvel2))*vel_L(1,-2)*(1._wp/dz(l-2))
                                !$acc atomic
                                rhs_vf(E_idx)%sf(j,k,l-3) = rhs_vf(E_idx)%sf(j,k,l-3) + 0.5_wp*mu_L(-3)*((2._wp/60._wp)*(dvel1+dvel2))*vel_L(1,-3)*(1._wp/dz(l-3))

                                !$acc atomic
                                rhs_vf(momxb)%sf(j,k,l-1) = rhs_vf(momxb)%sf(j,k,l-1) - 0.5_wp*mu_R(-2)*((-3._wp/60._wp)*(dvel1+dvel2))*(1._wp/dz(l-1))
                                !$acc atomic
                                rhs_vf(momxb)%sf(j,k,l) = rhs_vf(momxb)%sf(j,k,l) - 0.5_wp*mu_R(-1)*((27._wp/60._wp)*(dvel1+dvel2))*(1._wp/dz(l))
                                !$acc atomic
                                rhs_vf(momxb)%sf(j,k,l+1) = rhs_vf(momxb)%sf(j,k,l+1) - 0.5_wp*mu_R(0)*((47._wp/60._wp)*(dvel1+dvel2))*(1._wp/dz(l+1))
                                !$acc atomic
                                rhs_vf(momxb)%sf(j,k,l+2) = rhs_vf(momxb)%sf(j,k,l+2) - 0.5_wp*mu_R(1)*((-13._wp/60._wp)*(dvel1+dvel2))*(1._wp/dz(l+2))
                                !$acc atomic
                                rhs_vf(momxb)%sf(j,k,l+3) = rhs_vf(momxb)%sf(j,k,l+3) - 0.5_wp*mu_R(2)*((2._wp/60._wp)*(dvel1+dvel2))*(1._wp/dz(l+3))

                                !$acc atomic
                                rhs_vf(E_idx)%sf(j,k,l-1) = rhs_vf(E_idx)%sf(j,k,l-1) - 0.5_wp*mu_R(-2)*((-3._wp/60._wp)*(dvel1+dvel2))*vel_R(1,-2)*(1._wp/dz(l-1))
                                !$acc atomic
                                rhs_vf(E_idx)%sf(j,k,l) = rhs_vf(E_idx)%sf(j,k,l) - 0.5_wp*mu_R(-1)*((27._wp/60._wp)*(dvel1+dvel2))*vel_R(1,-1)*(1._wp/dz(l))
                                !$acc atomic
                                rhs_vf(E_idx)%sf(j,k,l+1) = rhs_vf(E_idx)%sf(j,k,l+1) - 0.5_wp*mu_R(0)*((47._wp/60._wp)*(dvel1+dvel2))*vel_R(1,0)*(1._wp/dz(l+1))
                                !$acc atomic
                                rhs_vf(E_idx)%sf(j,k,l+2) = rhs_vf(E_idx)%sf(j,k,l+2) - 0.5_wp*mu_R(1)*((-13._wp/60._wp)*(dvel1+dvel2))*vel_R(1,1)*(1._wp/dz(l+2))
                                !$acc atomic
                                rhs_vf(E_idx)%sf(j,k,l+3) = rhs_vf(E_idx)%sf(j,k,l+3) - 0.5_wp*mu_R(2)*((2._wp/60._wp)*(dvel1+dvel2))*vel_R(1,2)*(1._wp/dz(l+3)) 

                                !$acc atomic
                                rhs_vf(momxb)%sf(j,k,l-2) = rhs_vf(momxb)%sf(j,k,l-2) + 0.5_wp*mu_R(-2)*((-3._wp/60._wp)*(dvel1+dvel2))*(1._wp/dz(l-2))
                                !$acc atomic
                                rhs_vf(momxb)%sf(j,k,l-1) = rhs_vf(momxb)%sf(j,k,l-1) + 0.5_wp*mu_R(-1)*((27._wp/60._wp)*(dvel1+dvel2))*(1._wp/dz(l-1))
                                !$acc atomic
                                rhs_vf(momxb)%sf(j,k,l) = rhs_vf(momxb)%sf(j,k,l) + 0.5_wp*mu_R(0)*((47._wp/60._wp)*(dvel1+dvel2))*(1._wp/dz(l))
                                !$acc atomic
                                rhs_vf(momxb)%sf(j,k,l+1) = rhs_vf(momxb)%sf(j,k,l+1) + 0.5_wp*mu_R(1)*((-13._wp/60._wp)*(dvel1+dvel2))*(1._wp/dz(l+1))
                                !$acc atomic
                                rhs_vf(momxb)%sf(j,k,l+2) = rhs_vf(momxb)%sf(j,k,l+2) + 0.5_wp*mu_R(2)*((2._wp/60._wp)*(dvel1+dvel2))*(1._wp/dz(l+2))

                                !$acc atomic
                                rhs_vf(E_idx)%sf(j,k,l-2) = rhs_vf(E_idx)%sf(j,k,l-2) + 0.5_wp*mu_R(-2)*((-3._wp/60._wp)*(dvel1+dvel2))*vel_R(1,-2)*(1._wp/dz(l-2))
                                !$acc atomic
                                rhs_vf(E_idx)%sf(j,k,l-1) = rhs_vf(E_idx)%sf(j,k,l-1) + 0.5_wp*mu_R(-1)*((27._wp/60._wp)*(dvel1+dvel2))*vel_R(1,-1)*(1._wp/dz(l-1))
                                !$acc atomic
                                rhs_vf(E_idx)%sf(j,k,l) = rhs_vf(E_idx)%sf(j,k,l) + 0.5_wp*mu_R(0)*((47._wp/60._wp)*(dvel1+dvel2))*vel_R(1,0)*(1._wp/dz(l))
                                !$acc atomic
                                rhs_vf(E_idx)%sf(j,k,l+1) = rhs_vf(E_idx)%sf(j,k,l+1) + 0.5_wp*mu_R(1)*((-13._wp/60._wp)*(dvel1+dvel2))*vel_R(1,1)*(1._wp/dz(l+1))
                                !$acc atomic
                                rhs_vf(E_idx)%sf(j,k,l+2) = rhs_vf(E_idx)%sf(j,k,l+2) + 0.5_wp*mu_R(2)*((2._wp/60._wp)*(dvel1+dvel2))*vel_R(1,2)*(1._wp/dz(l+2))

                            end if

                            !! dwy & dvz
                            dvel1 = (1/(12._wp*dy(k))) * ( &
                            8._wp*q_prim_vf(momxb+2)%sf(j,k+1,l)/rho_sf(1,2) - &
                            8._wp*q_prim_vf(momxb+2)%sf(j,k-1,l)/rho_sf(-1,2)  + &
                            q_prim_vf(momxb+2)%sf(j,k-2,l)/rho_sf(-2,2)  - &
                            q_prim_vf(momxb+2)%sf(j,k+2,l)/rho_sf(2,2)  )

                            dvel2 = (1/(12._wp*dz(l))) * ( &
                            8._wp*q_prim_vf(momxb+1)%sf(j,k,l+1)/rho_sf(1,3) - &
                            8._wp*q_prim_vf(momxb+1)%sf(j,k,l-1)/rho_sf(-1,3)  + &
                            q_prim_vf(momxb+1)%sf(j,k,l-2)/rho_sf(-2,3)  - &
                            q_prim_vf(momxb+1)%sf(j,k,l+2)/rho_sf(2,3)  )

                            if(viscous) then

                                !$acc atomic
                                rhs_vf(momxb+1)%sf(j,k,l+2) = rhs_vf(momxb+1)%sf(j,k,l+2) - 0.5_wp*mu_L(1)*((-3._wp/60._wp)*(dvel1+dvel2))*(1._wp/dz(l+2))
                                !$acc atomic
                                rhs_vf(momxb+1)%sf(j,k,l+1) = rhs_vf(momxb+1)%sf(j,k,l+1) - 0.5_wp*mu_L(0)*((27._wp/60._wp)*(dvel1+dvel2))*(1._wp/dz(l+1))
                                !$acc atomic
                                rhs_vf(momxb+1)%sf(j,k,l) = rhs_vf(momxb+1)%sf(j,k,l) - 0.5_wp*mu_L(-1)*((47._wp/60._wp)*(dvel1+dvel2))*(1._wp/dz(l))
                                !$acc atomic
                                rhs_vf(momxb+1)%sf(j,k,l-1) = rhs_vf(momxb+1)%sf(j,k,l-1) - 0.5_wp*mu_L(-2)*((-13._wp/60._wp)*(dvel1+dvel2))*(1._wp/dz(l-1))
                                !$acc atomic
                                rhs_vf(momxb+1)%sf(j,k,l-2) = rhs_vf(momxb+1)%sf(j,k,l-2) - 0.5_wp*mu_L(-3)*((2._wp/60._wp)*(dvel1+dvel2))*(1._wp/dz(l-2))

                                !$acc atomic
                                rhs_vf(E_idx)%sf(j,k,l+2) = rhs_vf(E_idx)%sf(j,k,l+2) - 0.5_wp*mu_L(1)*((-3._wp/60._wp)*(dvel1+dvel2))*vel_L(2,1)*(1._wp/dz(l+2))
                                !$acc atomic
                                rhs_vf(E_idx)%sf(j,k,l+1) = rhs_vf(E_idx)%sf(j,k,l+1) - 0.5_wp*mu_L(0)*((27._wp/60._wp)*(dvel1+dvel2))*vel_L(2,0)*(1._wp/dz(l+1))
                                !$acc atomic
                                rhs_vf(E_idx)%sf(j,k,l) = rhs_vf(E_idx)%sf(j,k,l) - 0.5_wp*mu_L(-1)*((47._wp/60._wp)*(dvel1+dvel2))*vel_L(2,-1)*(1._wp/dz(l))
                                !$acc atomic
                                rhs_vf(E_idx)%sf(j,k,l-1) = rhs_vf(E_idx)%sf(j,k,l-1) - 0.5_wp*mu_L(-2)*((-13._wp/60._wp)*(dvel1+dvel2))*vel_L(2,-2)*(1._wp/dz(l-1))
                                !$acc atomic
                                rhs_vf(E_idx)%sf(j,k,l-2) = rhs_vf(E_idx)%sf(j,k,l-2) - 0.5_wp*mu_L(-3)*((2._wp/60._wp)*(dvel1+dvel2))*vel_L(2,-3)*(1._wp/dz(l-2))

                                !$acc atomic
                                rhs_vf(momxb+1)%sf(j,k,l+1) = rhs_vf(momxb+1)%sf(j,k,l+1) + 0.5_wp*mu_L(1)*((-3._wp/60._wp)*(dvel1+dvel2))*(1._wp/dz(l+1))
                                !$acc atomic
                                rhs_vf(momxb+1)%sf(j,k,l) = rhs_vf(momxb+1)%sf(j,k,l) + 0.5_wp*mu_L(0)*((27._wp/60._wp)*(dvel1+dvel2))*(1._wp/dz(l))
                                !$acc atomic
                                rhs_vf(momxb+1)%sf(j,k,l-1) = rhs_vf(momxb+1)%sf(j,k,l-1) + 0.5_wp*mu_L(-1)*((47._wp/60._wp)*(dvel1+dvel2))*(1._wp/dz(l-1))
                                !$acc atomic
                                rhs_vf(momxb+1)%sf(j,k,l-2) = rhs_vf(momxb+1)%sf(j,k,l-2) + 0.5_wp*mu_L(-2)*((-13._wp/60._wp)*(dvel1+dvel2))*(1._wp/dz(l-2))
                                !$acc atomic
                                rhs_vf(momxb+1)%sf(j,k,l-3) = rhs_vf(momxb+1)%sf(j,k,l-3) + 0.5_wp*mu_L(-3)*((2._wp/60._wp)*(dvel1+dvel2))*(1._wp/dz(l-3))

                                !$acc atomic
                                rhs_vf(E_idx)%sf(j,k,l+1) = rhs_vf(E_idx)%sf(j,k,l+1) + 0.5_wp*mu_L(1)*((-3._wp/60._wp)*(dvel1+dvel2))*vel_L(2,1)*(1._wp/dz(l+1))
                                !$acc atomic
                                rhs_vf(E_idx)%sf(j,k,l) = rhs_vf(E_idx)%sf(j,k,l) + 0.5_wp*mu_L(0)*((27._wp/60._wp)*(dvel1+dvel2))*vel_L(2,0)*(1._wp/dz(l))
                                !$acc atomic
                                rhs_vf(E_idx)%sf(j,k,l-1) = rhs_vf(E_idx)%sf(j,k,l-1) + 0.5_wp*mu_L(-1)*((47._wp/60._wp)*(dvel1+dvel2))*vel_L(2,-1)*(1._wp/dz(l-1))
                                !$acc atomic
                                rhs_vf(E_idx)%sf(j,k,l-2) = rhs_vf(E_idx)%sf(j,k,l-2) + 0.5_wp*mu_L(-2)*((-13._wp/60._wp)*(dvel1+dvel2))*vel_L(2,-2)*(1._wp/dz(l-2))
                                !$acc atomic
                                rhs_vf(E_idx)%sf(j,k,l-3) = rhs_vf(E_idx)%sf(j,k,l-3) + 0.5_wp*mu_L(-3)*((2._wp/60._wp)*(dvel1+dvel2))*vel_L(2,-3)*(1._wp/dz(l-3))
                               
                                !$acc atomic
                                rhs_vf(momxb+1)%sf(j,k,l-1) = rhs_vf(momxb+1)%sf(j,k,l-1) - 0.5_wp*mu_R(-2)*((-3._wp/60._wp)*(dvel1+dvel2))*(1._wp/dz(l-1))
                                !$acc atomic
                                rhs_vf(momxb+1)%sf(j,k,l) = rhs_vf(momxb+1)%sf(j,k,l) - 0.5_wp*mu_R(-1)*((27._wp/60._wp)*(dvel1+dvel2))*(1._wp/dz(l))
                                !$acc atomic
                                rhs_vf(momxb+1)%sf(j,k,l+1) = rhs_vf(momxb+1)%sf(j,k,l+1) - 0.5_wp*mu_R(0)*((47._wp/60._wp)*(dvel1+dvel2))*(1._wp/dz(l+1))
                                !$acc atomic
                                rhs_vf(momxb+1)%sf(j,k,l+2) = rhs_vf(momxb+1)%sf(j,k,l+2) - 0.5_wp*mu_R(1)*((-13._wp/60._wp)*(dvel1+dvel2))*(1._wp/dz(l+2))
                                !$acc atomic
                                rhs_vf(momxb+1)%sf(j,k,l+3) = rhs_vf(momxb+1)%sf(j,k,l+3) - 0.5_wp*mu_R(2)*((2._wp/60._wp)*(dvel1+dvel2))*(1._wp/dz(l+3))

                                !$acc atomic
                                rhs_vf(E_idx)%sf(j,k,l-1) = rhs_vf(E_idx)%sf(j,k,l-1) - 0.5_wp*mu_R(-2)*((-3._wp/60._wp)*(dvel1+dvel2))*vel_R(2,-2)*(1._wp/dz(l-1))
                                !$acc atomic
                                rhs_vf(E_idx)%sf(j,k,l) = rhs_vf(E_idx)%sf(j,k,l) - 0.5_wp*mu_R(-1)*((27._wp/60._wp)*(dvel1+dvel2))*vel_R(2,-1)*(1._wp/dz(l))
                                !$acc atomic
                                rhs_vf(E_idx)%sf(j,k,l+1) = rhs_vf(E_idx)%sf(j,k,l+1) - 0.5_wp*mu_R(0)*((47._wp/60._wp)*(dvel1+dvel2))*vel_R(2,0)*(1._wp/dz(l+1))
                                !$acc atomic
                                rhs_vf(E_idx)%sf(j,k,l+2) = rhs_vf(E_idx)%sf(j,k,l+2) - 0.5_wp*mu_R(1)*((-13._wp/60._wp)*(dvel1+dvel2))*vel_R(2,1)*(1._wp/dz(l+2))
                                !$acc atomic
                                rhs_vf(E_idx)%sf(j,k,l+3) = rhs_vf(E_idx)%sf(j,k,l+3) - 0.5_wp*mu_R(2)*((2._wp/60._wp)*(dvel1+dvel2))*vel_R(2,2)*(1._wp/dz(l+3)) 

                                !$acc atomic
                                rhs_vf(momxb+1)%sf(j,k,l-2) = rhs_vf(momxb+1)%sf(j,k,l-2) + 0.5_wp*mu_R(-2)*((-3._wp/60._wp)*(dvel1+dvel2))*(1._wp/dz(l-2))
                                !$acc atomic
                                rhs_vf(momxb+1)%sf(j,k,l-1) = rhs_vf(momxb+1)%sf(j,k,l-1) + 0.5_wp*mu_R(-1)*((27._wp/60._wp)*(dvel1+dvel2))*(1._wp/dz(l-1))
                                !$acc atomic
                                rhs_vf(momxb+1)%sf(j,k,l) = rhs_vf(momxb+1)%sf(j,k,l) + 0.5_wp*mu_R(0)*((47._wp/60._wp)*(dvel1+dvel2))*(1._wp/dz(l))
                                !$acc atomic
                                rhs_vf(momxb+1)%sf(j,k,l+1) = rhs_vf(momxb+1)%sf(j,k,l+1) + 0.5_wp*mu_R(1)*((-13._wp/60._wp)*(dvel1+dvel2))*(1._wp/dz(l+1))
                                !$acc atomic
                                rhs_vf(momxb+1)%sf(j,k,l+2) = rhs_vf(momxb+1)%sf(j,k,l+2) + 0.5_wp*mu_R(2)*((2._wp/60._wp)*(dvel1+dvel2))*(1._wp/dz(l+2))

                                !$acc atomic
                                rhs_vf(E_idx)%sf(j,k,l-2) = rhs_vf(E_idx)%sf(j,k,l-2) + 0.5_wp*mu_R(-2)*((-3._wp/60._wp)*(dvel1+dvel2))*vel_R(2,-2)*(1._wp/dz(l-2))
                                !$acc atomic
                                rhs_vf(E_idx)%sf(j,k,l-1) = rhs_vf(E_idx)%sf(j,k,l-1) + 0.5_wp*mu_R(-1)*((27._wp/60._wp)*(dvel1+dvel2))*vel_R(2,-1)*(1._wp/dz(l-1))
                                !$acc atomic
                                rhs_vf(E_idx)%sf(j,k,l) = rhs_vf(E_idx)%sf(j,k,l) + 0.5_wp*mu_R(0)*((47._wp/60._wp)*(dvel1+dvel2))*vel_R(2,0)*(1._wp/dz(l))
                                !$acc atomic
                                rhs_vf(E_idx)%sf(j,k,l+1) = rhs_vf(E_idx)%sf(j,k,l+1) + 0.5_wp*mu_R(1)*((-13._wp/60._wp)*(dvel1+dvel2))*vel_R(2,1)*(1._wp/dz(l+1))
                                !$acc atomic
                                rhs_vf(E_idx)%sf(j,k,l+2) = rhs_vf(E_idx)%sf(j,k,l+2) + 0.5_wp*mu_R(2)*((2._wp/60._wp)*(dvel1+dvel2))*vel_R(2,2)*(1._wp/dz(l+2))

                            end if

                            !! dwz & dux
                            dvel1 = (1/(12._wp*dz(l))) * ( &
                            8._wp*q_prim_vf(momxb+2)%sf(j,k,l+1)/rho_sf(1,3) - &
                            8._wp*q_prim_vf(momxb+2)%sf(j,k,l-1)/rho_sf(-1,3)  + &
                            q_prim_vf(momxb+2)%sf(j,k,l-2)/rho_sf(-2,3)  - &
                            q_prim_vf(momxb+2)%sf(j,k,l+2)/rho_sf(2,3)  )


                            dvel2 = (1/(12._wp*dx(j))) * ( &
                            8._wp*q_prim_vf(momxb)%sf(j+1,k,l)/rho_sf(1,1) - &
                            8._wp*q_prim_vf(momxb)%sf(j-1,k,l)/rho_sf(-1,1) + &
                            q_prim_vf(momxb)%sf(j-2,k,l)/rho_sf(-2,1) - &
                            q_prim_vf(momxb)%sf(j+2,k,l)/rho_sf(2,1) )

                            if(viscous) then

                                !$acc atomic
                                rhs_vf(momxb+2)%sf(j,k,l+2) = rhs_vf(momxb+2)%sf(j,k,l+2) - 0.5_wp*mu_L(1)*((-3._wp/60._wp)*(4_wp*dvel1-2_wp*dvel2)/3_wp)*(1._wp/dz(l+2))
                                !$acc atomic
                                rhs_vf(momxb+2)%sf(j,k,l+1) = rhs_vf(momxb+2)%sf(j,k,l+1) - 0.5_wp*mu_L(0)*((27._wp/60._wp)*(4_wp*dvel1-2_wp*dvel2)/3_wp)*(1._wp/dz(l+1))
                                !$acc atomic
                                rhs_vf(momxb+2)%sf(j,k,l) = rhs_vf(momxb+2)%sf(j,k,l) - 0.5_wp*mu_L(-1)*((47._wp/60._wp)*(4_wp*dvel1-2_wp*dvel2)/3_wp)*(1._wp/dz(l))
                                !$acc atomic
                                rhs_vf(momxb+2)%sf(j,k,l-1) = rhs_vf(momxb+2)%sf(j,k,l-1) - 0.5_wp*mu_L(-2)*((-13._wp/60._wp)*(4_wp*dvel1-2_wp*dvel2)/3_wp)*(1._wp/dz(l-1))
                                !$acc atomic
                                rhs_vf(momxb+2)%sf(j,k,l-2) = rhs_vf(momxb+2)%sf(j,k,l-2) - 0.5_wp*mu_L(-3)*((2._wp/60._wp)*(4_wp*dvel1-2_wp*dvel2)/3_wp)*(1._wp/dz(l-2))

                                !$acc atomic
                                rhs_vf(E_idx)%sf(j,k,l+2) = rhs_vf(E_idx)%sf(j,k,l+2) - 0.5_wp*mu_L(1)*((-3._wp/60._wp)*(4_wp*dvel1-2_wp*dvel2)/3_wp)*vel_L(3,1)*(1._wp/dz(l+2))
                                !$acc atomic
                                rhs_vf(E_idx)%sf(j,k,l+1) = rhs_vf(E_idx)%sf(j,k,l+1) - 0.5_wp*mu_L(0)*((27._wp/60._wp)*(4_wp*dvel1-2_wp*dvel2)/3_wp)*vel_L(3,0)*(1._wp/dz(l+1))
                                !$acc atomic
                                rhs_vf(E_idx)%sf(j,k,l) = rhs_vf(E_idx)%sf(j,k,l) - 0.5_wp*mu_L(-1)*((47._wp/60._wp)*(4_wp*dvel1-2_wp*dvel2)/3_wp)*vel_L(3,-1)*(1._wp/dz(l))
                                !$acc atomic
                                rhs_vf(E_idx)%sf(j,k,l-1) = rhs_vf(E_idx)%sf(j,k,l-1) - 0.5_wp*mu_L(-2)*((-13._wp/60._wp)*(4_wp*dvel1-2_wp*dvel2)/3_wp)*vel_L(3,-2)*(1._wp/dz(l-1))
                                !$acc atomic
                                rhs_vf(E_idx)%sf(j,k,l-2) = rhs_vf(E_idx)%sf(j,k,l-2) - 0.5_wp*mu_L(-3)*((2._wp/60._wp)*(4_wp*dvel1-2_wp*dvel2)/3_wp)*vel_L(3,-3)*(1._wp/dz(l-2))

                                !$acc atomic
                                rhs_vf(momxb+2)%sf(j,k,l+1) = rhs_vf(momxb+2)%sf(j,k,l+1) + 0.5_wp*mu_L(1)*((-3._wp/60._wp)*(4_wp*dvel1-2_wp*dvel2)/3_wp)*(1._wp/dz(l+1))
                                !$acc atomic
                                rhs_vf(momxb+2)%sf(j,k,l) = rhs_vf(momxb+2)%sf(j,k,l) + 0.5_wp*mu_L(0)*((27._wp/60._wp)*(4_wp*dvel1-2_wp*dvel2)/3_wp)*(1._wp/dz(l))
                                !$acc atomic
                                rhs_vf(momxb+2)%sf(j,k,l-1) = rhs_vf(momxb+2)%sf(j,k,l-1) + 0.5_wp*mu_L(-1)*((47._wp/60._wp)*(4_wp*dvel1-2_wp*dvel2)/3_wp)*(1._wp/dz(l-1))
                                !$acc atomic
                                rhs_vf(momxb+2)%sf(j,k,l-2) = rhs_vf(momxb+2)%sf(j,k,l-2) + 0.5_wp*mu_L(-2)*((-13._wp/60._wp)*(4_wp*dvel1-2_wp*dvel2)/3_wp)*(1._wp/dz(l-2))
                                !$acc atomic
                                rhs_vf(momxb+2)%sf(j,k,l-3) = rhs_vf(momxb+2)%sf(j,k,l-3) + 0.5_wp*mu_L(-3)*((2._wp/60._wp)*(4_wp*dvel1-2_wp*dvel2)/3_wp)*(1._wp/dz(l-3))

                                !$acc atomic
                                rhs_vf(E_idx)%sf(j,k,l+1) = rhs_vf(E_idx)%sf(j,k,l+1) + 0.5_wp*mu_L(1)*((-3._wp/60._wp)*(4_wp*dvel1-2_wp*dvel2)/3_wp)*vel_L(3,1)*(1._wp/dz(l+1))
                                !$acc atomic
                                rhs_vf(E_idx)%sf(j,k,l) = rhs_vf(E_idx)%sf(j,k,l) + 0.5_wp*mu_L(0)*((27._wp/60._wp)*(4_wp*dvel1-2_wp*dvel2)/3_wp)*vel_L(3,0)*(1._wp/dz(l))
                                !$acc atomic
                                rhs_vf(E_idx)%sf(j,k,l-1) = rhs_vf(E_idx)%sf(j,k,l-1) + 0.5_wp*mu_L(-1)*((47._wp/60._wp)*(4_wp*dvel1-2_wp*dvel2)/3_wp)*vel_L(3,-1)*(1._wp/dz(l-1))
                                !$acc atomic
                                rhs_vf(E_idx)%sf(j,k,l-2) = rhs_vf(E_idx)%sf(j,k,l-2) + 0.5_wp*mu_L(-2)*((-13._wp/60._wp)*(4_wp*dvel1-2_wp*dvel2)/3_wp)*vel_L(3,-2)*(1._wp/dz(l-2))
                                !$acc atomic
                                rhs_vf(E_idx)%sf(j,k,l-3) = rhs_vf(E_idx)%sf(j,k,l-3) + 0.5_wp*mu_L(-3)*((2._wp/60._wp)*(4_wp*dvel1-2_wp*dvel2)/3_wp)*vel_L(3,-3)*(1._wp/dz(l-3))
                               
                                !$acc atomic
                                rhs_vf(momxb+2)%sf(j,k,l-1) = rhs_vf(momxb+2)%sf(j,k,l-1) - 0.5_wp*mu_R(-2)*((-3._wp/60._wp)*(4_wp*dvel1-2_wp*dvel2)/3_wp)*(1._wp/dz(l-1))
                                !$acc atomic
                                rhs_vf(momxb+2)%sf(j,k,l) = rhs_vf(momxb+2)%sf(j,k,l) - 0.5_wp*mu_R(-1)*((27._wp/60._wp)*(4_wp*dvel1-2_wp*dvel2)/3_wp)*(1._wp/dz(l))
                                !$acc atomic
                                rhs_vf(momxb+2)%sf(j,k,l+1) = rhs_vf(momxb+2)%sf(j,k,l+1) - 0.5_wp*mu_R(0)*((47._wp/60._wp)*(4_wp*dvel1-2_wp*dvel2)/3_wp)*(1._wp/dz(l+1))
                                !$acc atomic
                                rhs_vf(momxb+2)%sf(j,k,l+2) = rhs_vf(momxb+2)%sf(j,k,l+2) - 0.5_wp*mu_R(1)*((-13._wp/60._wp)*(4_wp*dvel1-2_wp*dvel2)/3_wp)*(1._wp/dz(l+2))
                                !$acc atomic
                                rhs_vf(momxb+2)%sf(j,k,l+3) = rhs_vf(momxb+2)%sf(j,k,l+3) - 0.5_wp*mu_R(2)*((2._wp/60._wp)*(4_wp*dvel1-2_wp*dvel2)/3_wp)*(1._wp/dz(l+3))

                                !$acc atomic
                                rhs_vf(E_idx)%sf(j,k,l-1) = rhs_vf(E_idx)%sf(j,k,l-1) - 0.5_wp*mu_R(-2)*((-3._wp/60._wp)*(4_wp*dvel1-2_wp*dvel2)/3_wp)*vel_R(3,-2)*(1._wp/dz(l-1))
                                !$acc atomic
                                rhs_vf(E_idx)%sf(j,k,l) = rhs_vf(E_idx)%sf(j,k,l) - 0.5_wp*mu_R(-1)*((27._wp/60._wp)*(4_wp*dvel1-2_wp*dvel2)/3_wp)*vel_R(3,-1)*(1._wp/dz(l))
                                !$acc atomic
                                rhs_vf(E_idx)%sf(j,k,l+1) = rhs_vf(E_idx)%sf(j,k,l+1) - 0.5_wp*mu_R(0)*((47._wp/60._wp)*(4_wp*dvel1-2_wp*dvel2)/3_wp)*vel_R(3,0)*(1._wp/dz(l+1))
                                !$acc atomic
                                rhs_vf(E_idx)%sf(j,k,l+2) = rhs_vf(E_idx)%sf(j,k,l+2) - 0.5_wp*mu_R(1)*((-13._wp/60._wp)*(4_wp*dvel1-2_wp*dvel2)/3_wp)*vel_R(3,1)*(1._wp/dz(l+2))
                                !$acc atomic
                                rhs_vf(E_idx)%sf(j,k,l+3) = rhs_vf(E_idx)%sf(j,k,l+3) - 0.5_wp*mu_R(2)*((2._wp/60._wp)*(4_wp*dvel1-2_wp*dvel2)/3_wp)*vel_R(3,2)*(1._wp/dz(l+3)) 

                                !$acc atomic
                                rhs_vf(momxb+2)%sf(j,k,l-2) = rhs_vf(momxb+2)%sf(j,k,l-2) + 0.5_wp*mu_R(-2)*((-3._wp/60._wp)*(4_wp*dvel1-2_wp*dvel2)/3_wp)*(1._wp/dz(l-2))
                                !$acc atomic
                                rhs_vf(momxb+2)%sf(j,k,l-1) = rhs_vf(momxb+2)%sf(j,k,l-1) + 0.5_wp*mu_R(-1)*((27._wp/60._wp)*(4_wp*dvel1-2_wp*dvel2)/3_wp)*(1._wp/dz(l-1))
                                !$acc atomic
                                rhs_vf(momxb+2)%sf(j,k,l) = rhs_vf(momxb+2)%sf(j,k,l) + 0.5_wp*mu_R(0)*((47._wp/60._wp)*(4_wp*dvel1-2_wp*dvel2)/3_wp)*(1._wp/dz(l))
                                !$acc atomic
                                rhs_vf(momxb+2)%sf(j,k,l+1) = rhs_vf(momxb+2)%sf(j,k,l+1) + 0.5_wp*mu_R(1)*((-13._wp/60._wp)*(4_wp*dvel1-2_wp*dvel2)/3_wp)*(1._wp/dz(l+1))
                                !$acc atomic
                                rhs_vf(momxb+2)%sf(j,k,l+2) = rhs_vf(momxb+2)%sf(j,k,l+2) + 0.5_wp*mu_R(2)*((2._wp/60._wp)*(4_wp*dvel1-2_wp*dvel2)/3_wp)*(1._wp/dz(l+2))

                                !$acc atomic
                                rhs_vf(E_idx)%sf(j,k,l-2) = rhs_vf(E_idx)%sf(j,k,l-2) + 0.5_wp*mu_R(-2)*((-3._wp/60._wp)*(4_wp*dvel1-2_wp*dvel2)/3_wp)*vel_R(3,-2)*(1._wp/dz(l-2))
                                !$acc atomic
                                rhs_vf(E_idx)%sf(j,k,l-1) = rhs_vf(E_idx)%sf(j,k,l-1) + 0.5_wp*mu_R(-1)*((27._wp/60._wp)*(4_wp*dvel1-2_wp*dvel2)/3_wp)*vel_R(3,-1)*(1._wp/dz(l-1))
                                !$acc atomic
                                rhs_vf(E_idx)%sf(j,k,l) = rhs_vf(E_idx)%sf(j,k,l) + 0.5_wp*mu_R(0)*((47._wp/60._wp)*(4_wp*dvel1-2_wp*dvel2)/3_wp)*vel_R(3,0)*(1._wp/dz(l))
                                !$acc atomic
                                rhs_vf(E_idx)%sf(j,k,l+1) = rhs_vf(E_idx)%sf(j,k,l+1) + 0.5_wp*mu_R(1)*((-13._wp/60._wp)*(4_wp*dvel1-2_wp*dvel2)/3_wp)*vel_R(3,1)*(1._wp/dz(l+1))
                                !$acc atomic
                                rhs_vf(E_idx)%sf(j,k,l+2) = rhs_vf(E_idx)%sf(j,k,l+2) + 0.5_wp*mu_R(2)*((2._wp/60._wp)*(4_wp*dvel1-2_wp*dvel2)/3_wp)*vel_R(3,2)*(1._wp/dz(l+2))

                            end if


                            !!dvy
                            dvel2 = (1/(12._wp*dy(k))) * ( &
                            8._wp*q_prim_vf(momxb+1)%sf(j,k+1,l)/rho_sf(1,2) - &
                            8._wp*q_prim_vf(momxb+1)%sf(j,k-1,l)/rho_sf(-1,2) + &
                            q_prim_vf(momxb+1)%sf(j,k-2,l)/rho_sf(-2,2) - &
                            q_prim_vf(momxb+1)%sf(j,k+2,l)/rho_sf(2,2) )

                            if(viscous) then

                                !$acc atomic
                                rhs_vf(momxb+2)%sf(j,k,l+2) = rhs_vf(momxb+2)%sf(j,k,l+2) - 0.5_wp*mu_L(1)*((-3._wp/60._wp)*(-2_wp*dvel2)/3_wp)*(1._wp/dz(l+2))
                                !$acc atomic
                                rhs_vf(momxb+2)%sf(j,k,l+1) = rhs_vf(momxb+2)%sf(j,k,l+1) - 0.5_wp*mu_L(0)*((27._wp/60._wp)*(-2_wp*dvel2)/3_wp)*(1._wp/dz(l+1))
                                !$acc atomic
                                rhs_vf(momxb+2)%sf(j,k,l) = rhs_vf(momxb+2)%sf(j,k,l) - 0.5_wp*mu_L(-1)*((47._wp/60._wp)*(-2_wp*dvel2)/3_wp)*(1._wp/dz(l))
                                !$acc atomic
                                rhs_vf(momxb+2)%sf(j,k,l-1) = rhs_vf(momxb+2)%sf(j,k,l-1) - 0.5_wp*mu_L(-2)*((-13._wp/60._wp)*(-2_wp*dvel2)/3_wp)*(1._wp/dz(l-1))
                                !$acc atomic
                                rhs_vf(momxb+2)%sf(j,k,l-2) = rhs_vf(momxb+2)%sf(j,k,l-2) - 0.5_wp*mu_L(-3)*((2._wp/60._wp)*(-2_wp*dvel2)/3_wp)*(1._wp/dz(l-2))

                                !$acc atomic
                                rhs_vf(E_idx)%sf(j,k,l+2) = rhs_vf(E_idx)%sf(j,k,l+2) - 0.5_wp*mu_L(1)*((-3._wp/60._wp)*(-2_wp*dvel2)/3_wp)*vel_L(3,1)*(1._wp/dz(l+2))
                                !$acc atomic
                                rhs_vf(E_idx)%sf(j,k,l+1) = rhs_vf(E_idx)%sf(j,k,l+1) - 0.5_wp*mu_L(0)*((27._wp/60._wp)*(-2_wp*dvel2)/3_wp)*vel_L(3,0)*(1._wp/dz(l+1))
                                !$acc atomic
                                rhs_vf(E_idx)%sf(j,k,l) = rhs_vf(E_idx)%sf(j,k,l) - 0.5_wp*mu_L(-1)*((47._wp/60._wp)*(-2_wp*dvel2)/3_wp)*vel_L(3,-1)*(1._wp/dz(l))
                                !$acc atomic
                                rhs_vf(E_idx)%sf(j,k,l-1) = rhs_vf(E_idx)%sf(j,k,l-1) - 0.5_wp*mu_L(-2)*((-13._wp/60._wp)*(-2_wp*dvel2)/3_wp)*vel_L(3,-2)*(1._wp/dz(l-1))
                                !$acc atomic
                                rhs_vf(E_idx)%sf(j,k,l-2) = rhs_vf(E_idx)%sf(j,k,l-2) - 0.5_wp*mu_L(-3)*((2._wp/60._wp)*(-2_wp*dvel2)/3_wp)*vel_L(3,-3)*(1._wp/dz(l-2))

                                !$acc atomic
                                rhs_vf(momxb+2)%sf(j,k,l+1) = rhs_vf(momxb+2)%sf(j,k,l+1) + 0.5_wp*mu_L(1)*((-3._wp/60._wp)*(-2_wp*dvel2)/3_wp)*(1._wp/dz(l+1))
                                !$acc atomic
                                rhs_vf(momxb+2)%sf(j,k,l) = rhs_vf(momxb+2)%sf(j,k,l) + 0.5_wp*mu_L(0)*((27._wp/60._wp)*(-2_wp*dvel2)/3_wp)*(1._wp/dz(l))
                                !$acc atomic
                                rhs_vf(momxb+2)%sf(j,k,l-1) = rhs_vf(momxb+2)%sf(j,k,l-1) + 0.5_wp*mu_L(-1)*((47._wp/60._wp)*(-2_wp*dvel2)/3_wp)*(1._wp/dz(l-1))
                                !$acc atomic
                                rhs_vf(momxb+2)%sf(j,k,l-2) = rhs_vf(momxb+2)%sf(j,k,l-2) + 0.5_wp*mu_L(-2)*((-13._wp/60._wp)*(-2_wp*dvel2)/3_wp)*(1._wp/dz(l-2))
                                !$acc atomic
                                rhs_vf(momxb+2)%sf(j,k,l-3) = rhs_vf(momxb+2)%sf(j,k,l-3) + 0.5_wp*mu_L(-3)*((2._wp/60._wp)*(-2_wp*dvel2)/3_wp)*(1._wp/dz(l-3))

                                !$acc atomic
                                rhs_vf(E_idx)%sf(j,k,l+1) = rhs_vf(E_idx)%sf(j,k,l+1) + 0.5_wp*mu_L(1)*((-3._wp/60._wp)*(-2_wp*dvel2)/3_wp)*vel_L(3,1)*(1._wp/dz(l+1))
                                !$acc atomic
                                rhs_vf(E_idx)%sf(j,k,l) = rhs_vf(E_idx)%sf(j,k,l) + 0.5_wp*mu_L(0)*((27._wp/60._wp)*(-2_wp*dvel2)/3_wp)*vel_L(3,0)*(1._wp/dz(l))
                                !$acc atomic
                                rhs_vf(E_idx)%sf(j,k,l-1) = rhs_vf(E_idx)%sf(j,k,l-1) + 0.5_wp*mu_L(-1)*((47._wp/60._wp)*(-2_wp*dvel2)/3_wp)*vel_L(3,-1)*(1._wp/dz(l-1))
                                !$acc atomic
                                rhs_vf(E_idx)%sf(j,k,l-2) = rhs_vf(E_idx)%sf(j,k,l-2) + 0.5_wp*mu_L(-2)*((-13._wp/60._wp)*(-2_wp*dvel2)/3_wp)*vel_L(3,-2)*(1._wp/dz(l-2))
                                !$acc atomic
                                rhs_vf(E_idx)%sf(j,k,l-3) = rhs_vf(E_idx)%sf(j,k,l-3) + 0.5_wp*mu_L(-3)*((2._wp/60._wp)*(-2_wp*dvel2)/3_wp)*vel_L(3,-3)*(1._wp/dz(l-3))
                               
                                !$acc atomic
                                rhs_vf(momxb+2)%sf(j,k,l-1) = rhs_vf(momxb+2)%sf(j,k,l-1) - 0.5_wp*mu_R(-2)*((-3._wp/60._wp)*(-2_wp*dvel2)/3_wp)*(1._wp/dz(l-1))
                                !$acc atomic
                                rhs_vf(momxb+2)%sf(j,k,l) = rhs_vf(momxb+2)%sf(j,k,l) - 0.5_wp*mu_R(-1)*((27._wp/60._wp)*(-2_wp*dvel2)/3_wp)*(1._wp/dz(l))
                                !$acc atomic
                                rhs_vf(momxb+2)%sf(j,k,l+1) = rhs_vf(momxb+2)%sf(j,k,l+1) - 0.5_wp*mu_R(0)*((47._wp/60._wp)*(-2_wp*dvel2)/3_wp)*(1._wp/dz(l+1))
                                !$acc atomic
                                rhs_vf(momxb+2)%sf(j,k,l+2) = rhs_vf(momxb+2)%sf(j,k,l+2) - 0.5_wp*mu_R(1)*((-13._wp/60._wp)*(-2_wp*dvel2)/3_wp)*(1._wp/dz(l+2))
                                !$acc atomic
                                rhs_vf(momxb+2)%sf(j,k,l+3) = rhs_vf(momxb+2)%sf(j,k,l+3) - 0.5_wp*mu_R(2)*((2._wp/60._wp)*(-2_wp*dvel2)/3_wp)*(1._wp/dz(l+3))

                                !$acc atomic
                                rhs_vf(E_idx)%sf(j,k,l-1) = rhs_vf(E_idx)%sf(j,k,l-1) - 0.5_wp*mu_R(-2)*((-3._wp/60._wp)*(-2_wp*dvel2)/3_wp)*vel_R(3,-2)*(1._wp/dz(l-1))
                                !$acc atomic
                                rhs_vf(E_idx)%sf(j,k,l) = rhs_vf(E_idx)%sf(j,k,l) - 0.5_wp*mu_R(-1)*((27._wp/60._wp)*(-2_wp*dvel2)/3_wp)*vel_R(3,-1)*(1._wp/dz(l))
                                !$acc atomic
                                rhs_vf(E_idx)%sf(j,k,l+1) = rhs_vf(E_idx)%sf(j,k,l+1) - 0.5_wp*mu_R(0)*((47._wp/60._wp)*(-2_wp*dvel2)/3_wp)*vel_R(3,0)*(1._wp/dz(l+1))
                                !$acc atomic
                                rhs_vf(E_idx)%sf(j,k,l+2) = rhs_vf(E_idx)%sf(j,k,l+2) - 0.5_wp*mu_R(1)*((-13._wp/60._wp)*(-2_wp*dvel2)/3_wp)*vel_R(3,1)*(1._wp/dz(l+2))
                                !$acc atomic
                                rhs_vf(E_idx)%sf(j,k,l+3) = rhs_vf(E_idx)%sf(j,k,l+3) - 0.5_wp*mu_R(2)*((2._wp/60._wp)*(-2_wp*dvel2)/3_wp)*vel_R(3,2)*(1._wp/dz(l+3)) 

                                !$acc atomic
                                rhs_vf(momxb+2)%sf(j,k,l-2) = rhs_vf(momxb+2)%sf(j,k,l-2) + 0.5_wp*mu_R(-2)*((-3._wp/60._wp)*(-2_wp*dvel2)/3_wp)*(1._wp/dz(l-2))
                                !$acc atomic
                                rhs_vf(momxb+2)%sf(j,k,l-1) = rhs_vf(momxb+2)%sf(j,k,l-1) + 0.5_wp*mu_R(-1)*((27._wp/60._wp)*(-2_wp*dvel2)/3_wp)*(1._wp/dz(l-1))
                                !$acc atomic
                                rhs_vf(momxb+2)%sf(j,k,l) = rhs_vf(momxb+2)%sf(j,k,l) + 0.5_wp*mu_R(0)*((47._wp/60._wp)*(-2_wp*dvel2)/3_wp)*(1._wp/dz(l))
                                !$acc atomic
                                rhs_vf(momxb+2)%sf(j,k,l+1) = rhs_vf(momxb+2)%sf(j,k,l+1) + 0.5_wp*mu_R(1)*((-13._wp/60._wp)*(-2_wp*dvel2)/3_wp)*(1._wp/dz(l+1))
                                !$acc atomic
                                rhs_vf(momxb+2)%sf(j,k,l+2) = rhs_vf(momxb+2)%sf(j,k,l+2) + 0.5_wp*mu_R(2)*((2._wp/60._wp)*(-2_wp*dvel2)/3_wp)*(1._wp/dz(l+2))

                                !$acc atomic
                                rhs_vf(E_idx)%sf(j,k,l-2) = rhs_vf(E_idx)%sf(j,k,l-2) + 0.5_wp*mu_R(-2)*((-3._wp/60._wp)*(-2_wp*dvel2)/3_wp)*vel_R(3,-2)*(1._wp/dz(l-2))
                                !$acc atomic
                                rhs_vf(E_idx)%sf(j,k,l-1) = rhs_vf(E_idx)%sf(j,k,l-1) + 0.5_wp*mu_R(-1)*((27._wp/60._wp)*(-2_wp*dvel2)/3_wp)*vel_R(3,-1)*(1._wp/dz(l-1))
                                !$acc atomic
                                rhs_vf(E_idx)%sf(j,k,l) = rhs_vf(E_idx)%sf(j,k,l) + 0.5_wp*mu_R(0)*((47._wp/60._wp)*(-2_wp*dvel2)/3_wp)*vel_R(3,0)*(1._wp/dz(l))
                                !$acc atomic
                                rhs_vf(E_idx)%sf(j,k,l+1) = rhs_vf(E_idx)%sf(j,k,l+1) + 0.5_wp*mu_R(1)*((-13._wp/60._wp)*(-2_wp*dvel2)/3_wp)*vel_R(3,1)*(1._wp/dz(l+1))
                                !$acc atomic
                                rhs_vf(E_idx)%sf(j,k,l+2) = rhs_vf(E_idx)%sf(j,k,l+2) + 0.5_wp*mu_R(2)*((2._wp/60._wp)*(-2_wp*dvel2)/3_wp)*vel_R(3,2)*(1._wp/dz(l+2))

                            end if

                            F_L = (1._wp/60._wp) * (-3._wp * jac(j, k, l+2) + &
                                                27._wp * jac(j, k, l+1) + &
                                                47._wp * jac(j, k, l) -   &
                                                13._wp * jac(j, k, l-1) + &
                                                2._wp * jac(j, k, l-2))
                            !$acc loop seq
                            do i = 1, num_fluids
                                !$acc atomic
                                rhs_vf(i)%sf(j,k,l+1) = rhs_vf(i)%sf(j,k,l+1) + &
                                    (0.5_wp * (alpha_rho_R(i) * &
                                    vel_R(3,0))*(1._wp/dz(l+1)) + &
                                    0.5_wp*cfl * (alpha_rho_R(i))*(1._wp/dz(l+1)))

                                !$acc atomic
                                rhs_vf(i)%sf(j,k,l) = rhs_vf(i)%sf(j,k,l) - &
                                    (0.5_wp * (alpha_rho_R(i) * &
                                    vel_R(3,0))*(1._wp/dz(l)) + &
                                    0.5_wp*cfl * (alpha_rho_R(i))*(1._wp/dz(l)))
                            end do

                            if(num_fluids > 1) then 
                                !$acc loop seq
                                do i = 1, num_fluids
                                    !$acc atomic
                                    rhs_vf(advxb+i-1)%sf(j,k,l+1) = rhs_vf(advxb+i-1)%sf(j,k,l+1) + &
                                        (0.5_wp * (alpha_R(i) * &
                                        vel_R(3,0))*(1._wp/dz(l+1)) + &
                                        0.5_wp*cfl*(alpha_R(i))*(1._wp/dz(l+1)))

                                    !$acc atomic
                                    rhs_vf(advxb+i-1)%sf(j,k,l+1) = rhs_vf(advxb+i-1)%sf(j,k,l+1) &
                                    - (0.5_wp * q_prim_vf(advxb+i-1)%sf(j,k,l+1)  * vel_R(3,0)*(1._wp/dz(l+1)))

                                    !$acc atomic
                                    rhs_vf(advxb+i-1)%sf(j,k,l) = rhs_vf(advxb+i-1)%sf(j,k,l) - &
                                        (0.5_wp * (alpha_R(i) * &
                                        vel_R(3,0))*(1._wp/dz(l)) + &
                                        0.5_wp*cfl*(alpha_R(i))*(1._wp/dz(l)))

                                    !$acc atomic
                                    rhs_vf(advxb+i-1)%sf(j,k,l) = rhs_vf(advxb+i-1)%sf(j,k,l) &
                                    + (0.5_wp * q_prim_vf(advxb+i-1)%sf(j,k,l)  * vel_R(3,0)*(1._wp/dz(l)))
                                end do
                            end if

                            !$acc atomic
                            rhs_vf(momxb+2)%sf(j,k,l+1) = rhs_vf(momxb+2)%sf(j,k,l+1) + &
                                 (0.5_wp* (rho_R * (vel_R(3,0))**2.0 + &
                                 pres_R+F_L)*(1._wp/dz(l+1)) + &
                                 0.5_wp*cfl * (rho_R*vel_R(3,0))*(1._wp/dz(l+1)))

                            !$acc atomic
                            rhs_vf(momxb)%sf(j, k,l+1) =  rhs_vf(momxb)%sf(j, k,l+1) + &
                                (0.5_wp * rho_R * vel_R(1,0)*vel_R(3,0)*(1._wp/dz(l+1)) + &
                                0.5_wp*cfl * (rho_R*vel_R(1,0))*(1._wp/dz(l+1)))

                            !$acc atomic
                            rhs_vf(momxb+1)%sf(j,k,l+1) =  rhs_vf(momxb+1)%sf(j,k,l+1) + &
                                (0.5_wp * rho_R * vel_R(2,0)*vel_R(3,0)*(1._wp/dz(l+1)) + &
                                0.5_wp*cfl * (rho_R*vel_R(2,0))*(1._wp/dz(l+1)))

                            !$acc atomic
                            rhs_vf(E_idx)%sf(j,k,l+1) = rhs_vf(E_idx)%sf(j,k,l+1) +  &
                                (0.5_wp * (vel_R(3,0) * (E_R + &
                                pres_R+F_L) )*(1._wp/dz(l+1)) + &
                                0.5_wp*cfl * (E_R)*(1._wp/dz(l+1)))

                            !$acc atomic
                            rhs_vf(momxb+2)%sf(j,k,l) = rhs_vf(momxb+2)%sf(j,k,l) - &
                                 (0.5_wp* (rho_R * (vel_R(3,0))**2.0 + &
                                 pres_R+F_L)*(1._wp/dz(l)) + &
                                 0.5_wp*cfl * (rho_R*vel_R(3,0))*(1._wp/dz(l)))

                            !$acc atomic
                            rhs_vf(momxb)%sf(j, k, l) =  rhs_vf(momxb)%sf(j, k, l) - &
                                (0.5_wp * rho_R * vel_R(1,0)*vel_R(3,0)*(1._wp/dz(l)) + &
                                0.5_wp*cfl * (rho_R*vel_R(1,0))*(1._wp/dz(l)))

                            !$acc atomic
                            rhs_vf(momxb+1)%sf(j, k, l) =  rhs_vf(momxb+1)%sf(j, k, l) - &
                                (0.5_wp * rho_R * vel_R(2,0)*vel_R(3,0)*(1._wp/dz(l)) + &
                                0.5_wp*cfl * (rho_R*vel_R(2,0))*(1._wp/dz(l)))

                            !$acc atomic
                            rhs_vf(E_idx)%sf(j, k, l) = rhs_vf(E_idx)%sf(j, k, l) -  &
                                (0.5_wp * (vel_R(3,0) * (E_R + &
                                pres_R+F_L) )*(1._wp/dz(l)) + &
                                0.5_wp*cfl * (E_R)*(1._wp/dz(l)))

                        end do
                    end do
                end do
        end if

    end subroutine s_igr_riemann_solver

    subroutine s_initialize_igr(q_prim_vf, rhs_vf)

        type(scalar_field), dimension(sys_size), intent(in) :: q_prim_vf, rhs_vf

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
                    do i = 1, sys_size
                        rhs_vf(i)%sf(j,k,l) = 0._wp 
                    end do
                end do
            end do
        end do

        if(num_fluids > 1) then 
            !$acc parallel loop collapse(3) gang vector default(present)
            do l = idwbuff(3)%beg, idwbuff(3)%end
                do k = idwbuff(2)%beg, idwbuff(2)%end
                    do j = idwbuff(1)%beg, idwbuff(1)%end
                        rho_igr(j,k,l) = 0._wp
                        do i = 1, num_fluids
                            rho_igr(j,k,l) = rho_igr(j,k,l) + q_prim_vf(contxb + i - 1)%sf(j,k,l)
                        end do
                    end do
                end do
            end do
        end if

    end subroutine s_initialize_igr

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
        elseif (idir == 3) then
            !$acc parallel loop collapse(4) gang vector default(present)
            do i = 1, sys_size
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
        end if

    end subroutine s_igr_flux_add

end module m_igr
