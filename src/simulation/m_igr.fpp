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
    real(wp), allocatable, dimension(:, :, :) :: jac,jac_rhs,rho_igr,cfl
    real(wp), allocatable, dimension(:, :, :) :: dvel1, dvel2, dvel_L, FL
    !$acc declare create(fd_coeff,jac, jac_rhs,rho_igr,cfl)
    !$acc declare create(dvel1, dvel2, dvel_L, FL)

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

            @:ALLOCATE(qL_rs_vf(idwbuff(1)%beg:idwbuff(1)%end, &
                            idwbuff(2)%beg:idwbuff(2)%end, &
                            idwbuff(3)%beg:idwbuff(3)%end, &
                                1:sys_size))
            @:ALLOCATE(qR_rs_vf(idwbuff(1)%beg:idwbuff(1)%end, &
                                idwbuff(2)%beg:idwbuff(2)%end, &
                                idwbuff(3)%beg:idwbuff(3)%end, &
                                1:1))
            if(num_fluids > 1) then 
                @:ALLOCATE(qR_rs_vf(idwbuff(1)%beg:idwbuff(1)%end, &
                    idwbuff(2)%beg:idwbuff(2)%end, &
                    idwbuff(3)%beg:idwbuff(3)%end, &
                    advxb:advxe))
                 @:ALLOCATE(rho_igr(idwbuff(1)%beg:idwbuff(1)%end, &
                    idwbuff(2)%beg:idwbuff(2)%end, &
                    idwbuff(3)%beg:idwbuff(3)%end))
            end if
            #:for VAR in [ 'jac','jac_rhs','fd_coeff', 'dvel1', 'dvel2','cfl']
                @:ALLOCATE(${VAR}$(idwbuff(1)%beg:idwbuff(1)%end, &
                             idwbuff(2)%beg:idwbuff(2)%end, &
                             idwbuff(3)%beg:idwbuff(3)%end))
            #:endfor

            if(viscous) then
                #:for VAR in ['dvel_L']
                    @:ALLOCATE(${VAR}$(idwbuff(1)%beg:idwbuff(1)%end, &
                                      idwbuff(2)%beg:idwbuff(2)%end, &
                                      idwbuff(3)%beg:idwbuff(3)%end))
                #:endfor

            end if

            @:ALLOCATE(FL(idwbuff(1)%beg:idwbuff(1)%end, &
                          idwbuff(2)%beg:idwbuff(2)%end, &
                          idwbuff(3)%beg:idwbuff(3)%end))

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

    subroutine s_igr_sigma(flux_vf, idir)

        type(scalar_field), &
            dimension(sys_size), &
            intent(inout) :: flux_vf
        integer, intent(in) :: idir

        if(idir == 1) then 
            if(p == 0) then 
                !$acc parallel loop collapse(3) gang vector default(present) 
                do l = 0, p
                    do k = -buff_size+2, n+buff_size-3
                        do j = -buff_size+2, m+buff_size-3

                            FL(j+1, k, l) = (1._wp/60._wp) * (-3._wp * jac(j-1, k, l) + &
                                                27._wp * jac(j, k, l) + &
                                                47._wp * jac(j+1, k, l) -   &
                                                13._wp * jac(j+2, k, l) + &
                                                2._wp * jac(j+3, k, l))

                            flux_vf(momxb)%sf(j,k,l) = flux_vf(momxb)%sf(j,k,l) + &
                                                      0.5_wp * FL(j+1,k,l)

                            flux_vf(E_idx)%sf(j,k,l) = flux_vf(E_idx)%sf(j,k,l) + &
                                                      0.5_wp * qL_rs_vf(j+1, k, l, momxb) * FL(j+1,k,l)


                            FL(j, k, l) = (1._wp/60._wp) * (-3._wp * jac(j+2, k, l) + &
                                                27._wp * jac(j+1, k, l) + &
                                                47._wp * jac(j, k, l) -   &
                                                13._wp * jac(j-1, k, l) + &
                                                2._wp * jac(j-2, k, l))

                            flux_vf(momxb)%sf(j,k,l) = flux_vf(momxb)%sf(j,k,l) + &
                                                      0.5_wp * FL(j,k,l)

                            flux_vf(E_idx)%sf(j,k,l) = flux_vf(E_idx)%sf(j,k,l) + &
                                                      0.5_wp * qR_rs_vf(j, k, l, 1) * FL(j,k,l)                         


                        end do 
                    end do 
                end do
            else 
                !$acc parallel loop collapse(3) gang vector default(present) 
                do l = -buff_size+2, p+buff_size-3
                    do k = -buff_size+2, n+buff_size-3
                        do j = -buff_size+2, m+buff_size-3

                            FL(j+1, k, l) = (1._wp/60._wp) * (-3._wp * jac(j-1, k, l) + &
                                                27._wp * jac(j, k, l) + &
                                                47._wp * jac(j+1, k, l) -   &
                                                13._wp * jac(j+2, k, l) + &
                                                2._wp * jac(j+3, k, l))

                            flux_vf(momxb)%sf(j,k,l) = flux_vf(momxb)%sf(j,k,l) + &
                                                      0.5_wp * FL(j+1,k,l)

                            flux_vf(E_idx)%sf(j,k,l) = flux_vf(E_idx)%sf(j,k,l) + &
                                                      0.5_wp * qL_rs_vf(j+1, k, l, momxb) * FL(j+1,k,l)

                            FL(j, k, l) = (1._wp/60._wp) * (-3._wp * jac(j+2, k, l) + &
                                                27._wp * jac(j+1, k, l) + &
                                                47._wp * jac(j, k, l) -   &
                                                13._wp * jac(j-1, k, l) + &
                                                2._wp * jac(j-2, k, l))

                            flux_vf(momxb)%sf(j,k,l) = flux_vf(momxb)%sf(j,k,l) + &
                                                      0.5_wp * FL(j,k,l)

                            flux_vf(E_idx)%sf(j,k,l) = flux_vf(E_idx)%sf(j,k,l) + &
                                                      0.5_wp * qR_rs_vf(j, k, l, 1) * FL(j,k,l)  

                        end do 
                    end do 
                end do                
            end if
        else if(idir == 2) then 
            if(p == 0) then 
                !$acc parallel loop collapse(3) gang vector default(present) 
                do l = 0, p
                    do k = -buff_size+2, n+buff_size-3
                        do j = -buff_size+2, m+buff_size-3

                            FL(j, k+1, l) = (1._wp/60._wp) * (-3._wp * jac(j, k-1, l) + &
                                                27._wp * jac(j, k, l) + &
                                                47._wp * jac(j, k+1, l) -   &
                                                13._wp * jac(j, k+2, l) + &
                                                2._wp * jac(j, k+3, l))

                            flux_vf(momxb+1)%sf(j,k,l) = flux_vf(momxb+1)%sf(j,k,l) + &
                                                      0.5_wp * FL(j,k+1,l)

                            flux_vf(E_idx)%sf(j,k,l) = flux_vf(E_idx)%sf(j,k,l) + &
                                                      0.5_wp * qL_rs_vf(j, k+1, l, momxb+1) * FL(j,k+1,l)


                            FL(j, k, l) = (1._wp/60._wp) * (-3._wp * jac(j, k+2, l) + &
                                                27._wp * jac(j, k+1, l) + &
                                                47._wp * jac(j, k, l) -   &
                                                13._wp * jac(j, k-1, l) + &
                                                2._wp * jac(j, k-2, l))

                            flux_vf(momxb+1)%sf(j,k,l) = flux_vf(momxb+1)%sf(j,k,l) + &
                                                      0.5_wp * FL(j,k,l)

                            flux_vf(E_idx)%sf(j,k,l) = flux_vf(E_idx)%sf(j,k,l) + &
                                                      0.5_wp * qR_rs_vf(j, k, l, 1) * FL(j,k,l)                            

                        end do 
                    end do 
                end do
            else 
                !$acc parallel loop collapse(3) gang vector default(present) 
                do l = -buff_size+2, p+buff_size-3
                    do k = -buff_size+2, n+buff_size-3
                        do j = -buff_size+2, m+buff_size-3

                            FL(j, k+1, l) = (1._wp/60._wp) * (-3._wp * jac(j, k-1, l) + &
                                                27._wp * jac(j, k, l) + &
                                                47._wp * jac(j, k+1, l) -   &
                                                13._wp * jac(j, k+2, l) + &
                                                2._wp * jac(j, k+3, l))

                            flux_vf(momxb+1)%sf(j,k,l) = flux_vf(momxb+1)%sf(j,k,l) + &
                                                      0.5_wp * FL(j,k+1,l)

                            flux_vf(E_idx)%sf(j,k,l) = flux_vf(E_idx)%sf(j,k,l) + &
                                                      0.5_wp * qL_rs_vf(j, k+1, l, momxb+1) * FL(j,k+1,l)


                            FL(j, k, l) = (1._wp/60._wp) * (-3._wp * jac(j, k+2, l) + &
                                                27._wp * jac(j, k+1, l) + &
                                                47._wp * jac(j, k, l) -   &
                                                13._wp * jac(j, k-1, l) + &
                                                2._wp * jac(j, k-2, l))

                            flux_vf(momxb+1)%sf(j,k,l) = flux_vf(momxb+1)%sf(j,k,l) + &
                                                      0.5_wp * FL(j,k,l)

                            flux_vf(E_idx)%sf(j,k,l) = flux_vf(E_idx)%sf(j,k,l) + &
                                                      0.5_wp * qR_rs_vf(j, k, l, 1) * FL(j,k,l)    

                        end do 
                    end do 
                end do                
            end if  
        else if(idir == 3) then 
            !$acc parallel loop collapse(3) gang vector default(present) 
            do l = -buff_size+2, p+buff_size-3
                do k = -buff_size+2, n+buff_size-3
                    do j = -buff_size+2, m+buff_size-3

                        FL(j, k, l+1) = (1._wp/60._wp) * (-3._wp * jac(j, k, l-1) + &
                                            27._wp * jac(j, k, l) + &
                                            47._wp * jac(j, k, l+1) -   &
                                            13._wp * jac(j, k, l+2) + &
                                            2._wp * jac(j, k, l+3))

                        flux_vf(momxb+2)%sf(j,k,l) = flux_vf(momxb+2)%sf(j,k,l) + &
                                                  0.5_wp * FL(j,k,l+1)

                        flux_vf(E_idx)%sf(j,k,l) = flux_vf(E_idx)%sf(j,k,l) + &
                                                  0.5_wp * qL_rs_vf(j, k, l+1, momxb+2) * FL(j,k,l+1)


                        FL(j, k, l) = (1._wp/60._wp) * (-3._wp * jac(j, k, l+2) + &
                                            27._wp * jac(j, k, l+1) + &
                                            47._wp * jac(j, k, l) -   &
                                            13._wp * jac(j, k, l-1) + &
                                            2._wp * jac(j, k, l-2))

                        flux_vf(momxb+2)%sf(j,k,l) = flux_vf(momxb+2)%sf(j,k,l) + &
                                                  0.5_wp * FL(j,k,l)

                        flux_vf(E_idx)%sf(j,k,l) = flux_vf(E_idx)%sf(j,k,l) + &
                                                  0.5_wp * qR_rs_vf(j, k, l, 1) * FL(j,k,l)    

                    end do 
                end do 
            end do
        end if

    end subroutine s_igr_sigma

    subroutine s_igr_riemann_solver(q_prim_vf, flux_vf, idir)

        type(scalar_field), &
            dimension(sys_size), &
            intent(inout) :: flux_vf
        type(scalar_field), &
            dimension(sys_size), &
            intent(inout) :: q_prim_vf
        integer, intent(in) :: idir

        real(wp) :: rho_L, gamma_L, pi_inf_L, mu_L, a_L

        if (idir == 1) then
            if(p == 0) then
                !$acc parallel loop collapse(3) gang vector default(present) private(rho_L,gamma_L,pi_inf_L,mu_L,a_L)
                do l = 0, p
                    do k = -buff_size+2, n+buff_size-3
                        do j = -buff_size+2, m+buff_size-3

                            !$acc loop seq
                            do i = 1, sys_size
                                qL_rs_vf(j+1, k, l, i) = (1._wp/60._wp) * (-3._wp * q_prim_vf(i)%sf(j-1, k, l) + &
                                                27._wp * q_prim_vf(i)%sf(j, k, l) + &
                                                47._wp * q_prim_vf(i)%sf(j+1, k, l) -   &
                                                13._wp * q_prim_vf(i)%sf(j+2, k, l) + &
                                                2._wp * q_prim_vf(i)%sf(j+3, k, l))
                            end do

                            qR_rs_vf(j, k, l, 1) = (1._wp/60._wp) * (-3._wp * q_prim_vf(momxb)%sf(j+2, k, l) + &
                                                27._wp * q_prim_vf(momxb)%sf(j+1, k, l) + &
                                                47._wp * q_prim_vf(momxb)%sf(j, k, l) -   &
                                                13._wp * q_prim_vf(momxb)%sf(j-1, k, l) + &
                                                2._wp * q_prim_vf(momxb)%sf(j-2, k, l))
                            
                            if(num_fluids > 1) then  
                                !$acc loop seq 
                                do i = advxb, advxe                           
                                    qR_rs_vf(j, k, l,i) = (1._wp/60._wp) * (-3._wp * q_prim_vf(i)%sf(j+2, k, l) + &
                                                        27._wp * q_prim_vf(i)%sf(j+1, k, l) + &
                                                        47._wp * q_prim_vf(i)%sf(j, k, l) -   &
                                                        13._wp * q_prim_vf(i)%sf(j-1, k, l) + &
                                                        2._wp * q_prim_vf(i)%sf(j-2, k, l))
                                end do 
                            end if

                            rho_L = sum(qL_rs_vf(j+1,k,l,1:num_fluids))
                            gamma_L = sum(qL_rs_vf(j+1,k,l,advxb:advxb+num_fluids-1)*gammas)
                            pi_inf_L = sum(qL_rs_vf(j+1,k,l,advxb:advxb+num_fluids-1)*pi_infs)

                            a_L = sqrt((qL_rs_vf(j+1,k,l,E_idx)*(1._wp/gamma_L + 1._wp) + pi_inf_L / gamma_L) / rho_L)

                            cfl(j,k,l) = 2_wp*(sqrt(qL_rs_vf(j+1,k,l,momxb)**2._wp + qL_rs_vf(j+1,k,l,momxb+1)**2._wp) + a_L)

                            !$acc loop seq
                            do i = 1, num_fluids
                                flux_vf(i)%sf(j,k,l) = &
                                    0.5_wp * (qL_rs_vf(j+1,k,l,i) * &
                                    qL_rs_vf(j+1,k,l,momxb)) - &
                                    0.5_wp*cfl(j,k,l) * (qL_rs_vf(j+1,k,l,i))
                            end do

                            if(num_fluids > 1) then 
                                !$acc loop seq
                                do i = 1, num_fluids
                                    flux_vf(advxb+i-1)%sf(j,k,l) = &
                                        0.5_wp * (qL_rs_vf(j+1,k,l,advxb+i-1) * &
                                        qL_rs_vf(j+1,k,l,momxb )) - &
                                        0.5_wp*cfl(j,k,l) * (qL_rs_vf(j+1,k,l,advxb+i-1))

                                    flux_vf(advxb+i-1)%sf(j,k,l) = flux_vf(advxb+i-1)%sf(j,k,l) &
                                    - 0.5_wp * q_prim_vf(advxb+i-1)%sf(j,k,l) * qL_rs_vf(j+1,k,l,momxb)
                                end do
                            end if

                            ! Momentum -> rho*u^2 + p + [[F_igr]]
                            flux_vf(momxb)%sf(j,k,l) = &
                                 0.5_wp* (rho_L * (qL_rs_vf(j+1,k,l,momxb))**2.0 + &
                                 qL_rs_vf(j+1,k,l,E_idx)) - &
                                 0.5_wp*cfl(j,k,l) * (qL_rs_vf(j+1,k,l,momxb))

                            flux_vf(momxb+1)%sf(j, k, l) =  &
                                0.5_wp * rho_L * qL_rs_vf(j+1,k,l,momxb)*qL_rs_vf(j+1,k,l,momxb+1) - &
                                0.5_wp*cfl(j,k,l) * (qL_rs_vf(j+1,k,l,momxb+1))

                             flux_vf(E_idx)%sf(j, k, l) = &
                                0.5_wp * (qL_rs_vf(j+1,k,l,momxb) * (qL_rs_vf(j+1,k,l,E_idx)*gamma_L + pi_inf_L + &
                                0.5_wp * rho_L * (qL_rs_vf(j+1,k,l,momxb)**2._wp + qL_rs_vf(j+1,k,l,momxb+1)**2._wp ) + &
                                qL_rs_vf(j+1,k,l,E_idx)) ) - &
                                0.5_wp*cfl(j,k,l) * (qL_rs_vf(j+1,k,l,E_idx))

                            !! duy & dvx
                            dvel1(j,k,l) = (1/(12._wp*dy(k))) * ( &
                            8._wp*q_prim_vf(momxb)%sf(j,k+1,l) - &
                            8._wp*q_prim_vf(momxb)%sf(j,k-1,l) + &
                            q_prim_vf(momxb)%sf(j,k-2,l) - &
                            q_prim_vf(momxb)%sf(j,k+2,l) )

                            dvel2(j,k,l) = (1/(12._wp*dx(j))) * ( &
                            8._wp*q_prim_vf(momxb+1)%sf(j+1,k,l) - &
                            8._wp*q_prim_vf(momxb+1)%sf(j-1,k,l) + &
                            q_prim_vf(momxb+1)%sf(j-2,k,l) - &
                            q_prim_vf(momxb+1)%sf(j+2,k,l) )

                            jac_rhs(j, k, l) = alf_igr* (2_wp * dvel1(j,k,l)* dvel2(j,k,l))

                        end do 
                    end do 
                end do

                !$acc parallel loop collapse(3) gang vector default(present) private(rho_L,gamma_L,pi_inf_L,mu_L,a_L)
                do l = 0, p
                    do k = -buff_size+2, n+buff_size-3
                        do j = -buff_size+2, m+buff_size-3

                            if(viscous) then
                                dvel_L(j+1, k, l) = (1._wp/60._wp) * (-3._wp * dvel1(j-1, k, l) + &
                                                        27._wp * dvel1(j, k, l) + &
                                                        47._wp * dvel1(j+1, k, l) -   &
                                                        13._wp * dvel1(j+2, k, l) + &
                                                        2._wp * dvel1(j+3, k, l))

                                mu_L = 0._wp
                                !$acc loop seq
                                do q = 1, Re_size(1)
                                    mu_L =  qL_rs_vf(j+1,k,l,advxb+Re_idx(1, q)-1) / Res(1, q) &
                                              + mu_L
                                end do

                                flux_vf(momxb+1)%sf(j,k,l) = flux_vf(momxb+1)%sf(j,k,l) - 0.5_wp*mu_L*(dvel_L(j+1,k,l))
                                flux_vf(E_idx)%sf(j,k,l) = flux_vf(E_idx)%sf(j,k,l) - 0.5_wp*mu_L*qL_rs_vf(j+1,k,l,momxb)*dvel_L(j+1,k,l) 

                                dvel_L(j+1, k, l) = (1._wp/60._wp) * (-3._wp * dvel2(j-1, k, l) + &
                                                        27._wp * dvel2(j, k, l) + &
                                                        47._wp * dvel2(j+1, k, l) -   &
                                                        13._wp * dvel2(j+2, k, l) + &
                                                        2._wp * dvel2(j+3, k, l))

                                flux_vf(momxb+1)%sf(j,k,l) = flux_vf(momxb+1)%sf(j,k,l) - 0.5_wp*mu_L*(dvel_L(j+1,k,l))
                                flux_vf(E_idx)%sf(j,k,l) = flux_vf(E_idx)%sf(j,k,l) - 0.5_wp*mu_L*qL_rs_vf(j+1,k,l,momxb)*dvel_L(j+1,k,l)

                                if(num_fluids > 1) then 
                                    mu_L = 0._wp
                                    !$acc loop seq
                                    do q = 1, Re_size(1)
                                        mu_L =  qR_rs_vf(j,k,l,advxb+Re_idx(1, q)-1) / Res(1, q) &
                                                  + mu_L
                                    end do
                                end if

                                dvel_L(j, k, l) = (1._wp/60._wp) * (-3._wp * dvel1(j+2, k, l) + &
                                                        27._wp * dvel1(j+1, k, l) + &
                                                        47._wp * dvel1(j, k, l) -   &
                                                        13._wp * dvel1(j-1, k, l) + &
                                                        2._wp * dvel1(j-2, k, l))
                                

                                flux_vf(momxb+1)%sf(j,k,l) = flux_vf(momxb+1)%sf(j,k,l) - 0.5_wp*mu_L*(dvel_L(j,k,l))
                                flux_vf(E_idx)%sf(j,k,l) = flux_vf(E_idx)%sf(j,k,l) - 0.5_wp*mu_L*qR_rs_vf(j,k,l,1)*dvel_L(j,k,l) 

                                dvel_L(j, k, l) = (1._wp/60._wp) * (-3._wp * dvel2(j+2, k, l) + &
                                                        27._wp * dvel2(j+1, k, l) + &
                                                        47._wp * dvel2(j, k, l) -   &
                                                        13._wp * dvel2(j-1, k, l) + &
                                                        2._wp * dvel2(j-2, k, l))

                                flux_vf(momxb+1)%sf(j,k,l) = flux_vf(momxb+1)%sf(j,k,l) - 0.5_wp*mu_L*(dvel_L(j,k,l))
                                flux_vf(E_idx)%sf(j,k,l) = flux_vf(E_idx)%sf(j,k,l) - 0.5_wp*mu_L*qR_rs_vf(j,k,l,1)*dvel_L(j,k,l)

                            end if

                            !! dux & dvy
                            dvel1(j,k,l) = (1/(12._wp*dx(j))) * ( &
                            8._wp*q_prim_vf(momxb)%sf(j+1,k,l) - &
                            8._wp*q_prim_vf(momxb)%sf(j-1,k,l) + &
                            q_prim_vf(momxb)%sf(j-2,k,l) - &
                            q_prim_vf(momxb)%sf(j+2,k,l) )

                            dvel2(j,k,l) = (1/(12._wp*dy(k))) * ( &
                            8._wp*q_prim_vf(momxb+1)%sf(j,k+1,l) - &
                            8._wp*q_prim_vf(momxb+1)%sf(j,k-1,l) + &
                            q_prim_vf(momxb+1)%sf(j,k-2,l) - &
                            q_prim_vf(momxb+1)%sf(j,k+2,l) )

                            jac_rhs(j, k, l) = jac_rhs(j,k,l) + alf_igr* (dvel1(j,k,l)**2_wp + dvel2(j,k,l)**2_wp + (dvel1(j,k,l) + dvel2(j,k,l))**2_wp)

                        end do 
                    end do 
                end do

                !$acc parallel loop collapse(3) gang vector default(present) private(rho_L,gamma_L,pi_inf_L,mu_L,a_L)
                do l = 0, p
                    do k = -buff_size+2, n+buff_size-3
                        do j = -buff_size+2, m+buff_size-3
                            if(viscous) then

                                mu_L = 0._wp
                                !$acc loop seq
                                do q = 1, Re_size(1)
                                    mu_L =  qL_rs_vf(j+1,k,l,advxb+Re_idx(1, q)-1) / Res(1, q) &
                                              + mu_L
                                end do

                                dvel_L(j+1, k, l) = (1._wp/60._wp) * (-3._wp * dvel1(j-1, k, l) + &
                                                        27._wp * dvel1(j, k, l) + &
                                                        47._wp * dvel1(j+1, k, l) -   &
                                                        13._wp * dvel1(j+2, k, l) + &
                                                        2._wp * dvel1(j+3, k, l))

                                flux_vf(momxb)%sf(j,k,l) = flux_vf(momxb)%sf(j,k,l) - 0.5_wp*mu_L*((4._wp/3._wp)*dvel_L(j+1,k,l))
                                flux_vf(E_idx)%sf(j,k,l) = flux_vf(E_idx)%sf(j,k,l) - 0.5_wp*mu_L*qL_rs_vf(j+1,k,l,momxb)*(4._wp/3._wp)*dvel_L(j+1,k,l)

                                dvel_L(j+1, k, l) = (1._wp/60._wp) * (-3._wp * dvel2(j-1, k, l) + &
                                                        27._wp * dvel2(j, k, l) + &
                                                        47._wp * dvel2(j+1, k, l) -   &
                                                        13._wp * dvel2(j+2, k, l) + &
                                                        2._wp * dvel2(j+3, k, l))

                                flux_vf(momxb)%sf(j,k,l) = flux_vf(momxb)%sf(j,k,l) + 0.5_wp*mu_L*((2._wp/3._wp)*dvel_L(j+1,k,l))
                                flux_vf(E_idx)%sf(j,k,l) = flux_vf(E_idx)%sf(j,k,l) + 0.5_wp*mu_L*qL_rs_vf(j+1,k,l,momxb)*((2._wp/3._wp)*dvel_L(j+1,k,l))

                                if(num_fluids > 1) then 
                                    mu_L = 0._wp
                                    !$acc loop seq
                                    do q = 1, Re_size(1)
                                        mu_L =  qR_rs_vf(j,k,l,advxb+Re_idx(1, q)-1) / Res(1, q) &
                                                  + mu_L
                                    end do
                                end if

                                dvel_L(j, k, l) = (1._wp/60._wp) * (-3._wp * dvel1(j+2, k, l) + &
                                        27._wp * dvel1(j+1, k, l) + &
                                        47._wp * dvel1(j, k, l) -   &
                                        13._wp * dvel1(j-1, k, l) + &
                                        2._wp * dvel1(j-2, k, l))
                                

                                flux_vf(momxb)%sf(j,k,l) = flux_vf(momxb)%sf(j,k,l) - 0.5_wp*mu_L*((4._wp/3._wp)*dvel_L(j,k,l))
                                flux_vf(E_idx)%sf(j,k,l) = flux_vf(E_idx)%sf(j,k,l) - 0.5_wp*mu_L*qR_rs_vf(j,k,l,1)*(4._wp/3._wp)*dvel_L(j,k,l) 

                                dvel_L(j, k, l) = (1._wp/60._wp) * (-3._wp * dvel2(j+2, k, l) + &
                                                        27._wp * dvel2(j+1, k, l) + &
                                                        47._wp * dvel2(j, k, l) -   &
                                                        13._wp * dvel2(j-1, k, l) + &
                                                        2._wp * dvel2(j-2, k, l))

                                flux_vf(momxb)%sf(j,k,l) = flux_vf(momxb)%sf(j,k,l) + 0.5_wp*mu_L*((2._wp/3._wp)*dvel_L(j,k,l))
                                flux_vf(E_idx)%sf(j,k,l) = flux_vf(E_idx)%sf(j,k,l) + 0.5_wp*mu_L*qR_rs_vf(j,k,l,1)*(2._wp/3._wp)*dvel_L(j,k,l)

                            end if

                            !$acc loop seq 
                            do i = 1, momxb - 1
                                qL_rs_vf(j, k, l, i) = (1._wp/60._wp) * (-3._wp * q_prim_vf(i)%sf(j+2, k, l) + &
                                                27._wp * q_prim_vf(i)%sf(j+1, k, l) + &
                                                47._wp * q_prim_vf(i)%sf(j, k, l) -   &
                                                13._wp * q_prim_vf(i)%sf(j-1, k, l) + &
                                                2._wp * q_prim_vf(i)%sf(j-2, k, l))
                            end do

                            !$acc loop seq 
                            do i = momxb+1,advxe
                                qL_rs_vf(j, k, l, i) = (1._wp/60._wp) * (-3._wp * q_prim_vf(i)%sf(j+2, k, l) + &
                                                27._wp * q_prim_vf(i)%sf(j+1, k, l) + &
                                                47._wp * q_prim_vf(i)%sf(j, k, l) -   &
                                                13._wp * q_prim_vf(i)%sf(j-1, k, l) + &
                                                2._wp * q_prim_vf(i)%sf(j-2, k, l))
                            end do 

                            if(sys_size > advxe) then 
                                do i = advxe + 1, sys_size
                                  qL_rs_vf(j, k, l, i) = (1._wp/60._wp) * (-3._wp * q_prim_vf(i)%sf(j+2, k, l) + &
                                                27._wp * q_prim_vf(i)%sf(j+1, k, l) + &
                                                47._wp * q_prim_vf(i)%sf(j, k, l) -   &
                                                13._wp * q_prim_vf(i)%sf(j-1, k, l) + &
                                                2._wp * q_prim_vf(i)%sf(j-2, k, l))
                                end do
                            end if  

                            rho_L = sum(qL_rs_vf(j,k,l,1:num_fluids))
                            gamma_L = sum(qL_rs_vf(j,k,l,advxb:advxb+num_fluids-1)*gammas)
                            pi_inf_L = sum(qL_rs_vf(j,k,l,advxb:advxb+num_fluids-1)*pi_infs)

                            !$acc loop seq
                            do i = 1, num_fluids
                                flux_vf(i)%sf(j,k,l) = flux_vf(i)%sf(j,k,l) + &
                                    0.5_wp * (qL_rs_vf(j,k,l,i) * &
                                    qR_rs_vf(j,k,l,1)) + &
                                    0.5_wp*cfl(j,k,l) * (qL_rs_vf(j,k,l,i))
                            end do

                            if(num_fluids > 1) then 
                                !$acc loop seq
                                do i = 1, num_fluids
                                    flux_vf(advxb+i-1)%sf(j,k,l) = &
                                        0.5_wp * (qR_rs_vf(j,k,l,advxb+i-1) * &
                                        qR_rs_vf(j,k,l,1)) + &
                                        0.5_wp*cfl(j,k,l) * (qR_rs_vf(j,k,l,advxb+i-1))

                                    flux_vf(advxb+i-1)%sf(j,k,l) = flux_vf(advxb+i-1)%sf(j,k,l) &
                                    - 0.5_wp * q_prim_vf(advxb+i-1)%sf(j,k,l) * qR_rs_vf(j,k,l,1)
                                end do
                            end if

                            ! Momentum -> rho*u^2 + p + [[F_igr]]
                            flux_vf(momxb)%sf(j,k,l) = flux_vf(momxb)%sf(j,k,l) + &
                                 0.5_wp* (rho_L * (qR_rs_vf(j,k,l,1))**2.0 + &
                                 qL_rs_vf(j,k,l,E_idx)) + &
                                 0.5_wp*cfl(j,k,l) * (qR_rs_vf(j,k,l,1))

                            flux_vf(momxb+1)%sf(j, k, l) = flux_vf(momxb+1)%sf(j, k, l) + &
                                0.5_wp * rho_L * qR_rs_vf(j,k,l,1)*qL_rs_vf(j,k,l,momxb+1) + &
                                0.5_wp*cfl(j,k,l) * (qL_rs_vf(j,k,l,momxb+1))

                             flux_vf(E_idx)%sf(j, k, l) = flux_vf(E_idx)%sf(j, k, l) +  &
                                0.5_wp * (qR_rs_vf(j,k,l,1) * (qL_rs_vf(j,k,l,E_idx)*gamma_L + pi_inf_L + &
                                0.5_wp * rho_L * (qR_rs_vf(j,k,l,1)**2._wp + qL_rs_vf(j,k,l,momxb+1)**2._wp  ) + &
                                qL_rs_vf(j,k,l,E_idx)) ) + &
                                0.5_wp*cfl(j,k,l) * (qL_rs_vf(j,k,l,E_idx))

                        end do
                    end do
                end do
            else
                !$acc parallel loop collapse(3) gang vector default(present) private(rho_L,gamma_L,pi_inf_L,mu_L,a_L)
                do l = -buff_size+2,p+buff_size-3
                    do k = -buff_size+2, n+buff_size-3
                        do j = -buff_size+2, m+buff_size-3

                            !$acc loop seq
                            do i = 1, sys_size
                                qL_rs_vf(j+1, k, l, i) = (1._wp/60._wp) * (-3._wp * q_prim_vf(i)%sf(j-1, k, l) + &
                                                27._wp * q_prim_vf(i)%sf(j, k, l) + &
                                                47._wp * q_prim_vf(i)%sf(j+1, k, l) -   &
                                                13._wp * q_prim_vf(i)%sf(j+2, k, l) + &
                                                2._wp * q_prim_vf(i)%sf(j+3, k, l))
                            end do

                            qR_rs_vf(j, k, l, 1) = (1._wp/60._wp) * (-3._wp * q_prim_vf(momxb)%sf(j+2, k, l) + &
                                                27._wp * q_prim_vf(momxb)%sf(j+1, k, l) + &
                                                47._wp * q_prim_vf(momxb)%sf(j, k, l) -   &
                                                13._wp * q_prim_vf(momxb)%sf(j-1, k, l) + &
                                                2._wp * q_prim_vf(momxb)%sf(j-2, k, l))
                            
                            if(num_fluids > 1) then  
                                !$acc loop seq 
                                do i = advxb, advxe                           
                                    qR_rs_vf(j, k, l,i) = (1._wp/60._wp) * (-3._wp * q_prim_vf(i)%sf(j+2, k, l) + &
                                                        27._wp * q_prim_vf(i)%sf(j+1, k, l) + &
                                                        47._wp * q_prim_vf(i)%sf(j, k, l) -   &
                                                        13._wp * q_prim_vf(i)%sf(j-1, k, l) + &
                                                        2._wp * q_prim_vf(i)%sf(j-2, k, l))
                                end do 
                            end if

                            rho_L = sum(qL_rs_vf(j+1,k,l,1:num_fluids))
                            gamma_L = sum(qL_rs_vf(j+1,k,l,advxb:advxb+num_fluids-1)*gammas)
                            pi_inf_L = sum(qL_rs_vf(j+1,k,l,advxb:advxb+num_fluids-1)*pi_infs)

                            a_L = sqrt((qL_rs_vf(j+1,k,l,E_idx)*(1._wp/gamma_L + 1._wp) + pi_inf_L / gamma_L) / rho_L)

                            cfl(j,k,l) = (sqrt(qL_rs_vf(j+1,k,l,momxb)**2._wp + qL_rs_vf(j+1,k,l,momxb+1)**2._wp + qL_rs_vf(j+1,k,l,momxb+2)**2._wp) + a_L)

                            !$acc loop seq
                            do i = 1, num_fluids
                                flux_vf(i)%sf(j,k,l) = &
                                    0.5_wp * (qL_rs_vf(j+1,k,l,i) * &
                                    qL_rs_vf(j+1,k,l,momxb)) - &
                                    0.5_wp*cfl(j,k,l) * (qL_rs_vf(j+1,k,l,i))
                            end do

                            if(num_fluids > 1) then 
                                !$acc loop seq
                                do i = 1, num_fluids
                                    flux_vf(advxb+i-1)%sf(j,k,l) = &
                                        0.5_wp * (qL_rs_vf(j+1,k,l,advxb+i-1) * &
                                        qL_rs_vf(j+1,k,l,momxb )) - &
                                        0.5_wp*cfl(j,k,l) * (qL_rs_vf(j+1,k,l,advxb+i-1))

                                    flux_vf(advxb+i-1)%sf(j,k,l) = flux_vf(advxb+i-1)%sf(j,k,l) &
                                    - 0.5_wp * q_prim_vf(advxb+i-1)%sf(j,k,l) * qL_rs_vf(j+1,k,l,momxb)
                                end do
                            end if

                            ! Momentum -> rho*u^2 + p + [[F_igr]]
                            flux_vf(momxb)%sf(j,k,l) = &
                                 0.5_wp* (rho_L * (qL_rs_vf(j+1,k,l,momxb))**2.0 + &
                                 qL_rs_vf(j+1,k,l,E_idx)) - &
                                 0.5_wp*cfl(j,k,l) * (qL_rs_vf(j+1,k,l,momxb))

                            flux_vf(momxb+1)%sf(j, k, l) =  &
                                0.5_wp * rho_L * qL_rs_vf(j+1,k,l,momxb)*qL_rs_vf(j+1,k,l,momxb+1) - &
                                0.5_wp*cfl(j,k,l) * (qL_rs_vf(j+1,k,l,momxb+1))

                            flux_vf(momxb+2)%sf(j, k, l) = &
                                0.5_wp * rho_L * qL_rs_vf(j+1,k,l,momxb)*qL_rs_vf(j+1,k,l,momxb+2) - &
                                0.5_wp*cfl(j,k,l) * (qL_rs_vf(j+1,k,l,momxb+2))

                             flux_vf(E_idx)%sf(j, k, l) = &
                                0.5_wp * (qL_rs_vf(j+1,k,l,momxb) * (qL_rs_vf(j+1,k,l,E_idx)*gamma_L + pi_inf_L + &
                                0.5_wp * rho_L * (qL_rs_vf(j+1,k,l,momxb)**2._wp + qL_rs_vf(j+1,k,l,momxb+1)**2._wp  + qL_rs_vf(j+1,k,l,momxb+2)**2._wp ) + &
                                qL_rs_vf(j+1,k,l,E_idx)) ) - &
                                0.5_wp*cfl(j,k,l) * (qL_rs_vf(j+1,k,l,E_idx))

                            !! duz & dwx
                            dvel1(j,k,l) = (1/(12._wp*dz(l))) * ( &
                            8._wp*q_prim_vf(momxb)%sf(j,k,l+1) - &
                            8._wp*q_prim_vf(momxb)%sf(j,k,l-1) + &
                            q_prim_vf(momxb)%sf(j,k,l-2) - &
                            q_prim_vf(momxb)%sf(j,k,l+2) )

                            dvel2(j,k,l) = (1/(12._wp*dx(j))) * ( &
                            8._wp*q_prim_vf(momxb+2)%sf(j+1,k,l) - &
                            8._wp*q_prim_vf(momxb+2)%sf(j-1,k,l) + &
                            q_prim_vf(momxb+2)%sf(j-2,k,l) - &
                            q_prim_vf(momxb+2)%sf(j+2,k,l) )

                            jac_rhs(j, k, l) =  alf_igr* (2_wp * dvel1(j,k,l)* dvel2(j,k,l))
                        end do 
                    end do 
                end do 

                !$acc parallel loop collapse(3) gang vector default(present) private(rho_L,gamma_L,pi_inf_L,mu_L,a_L)
                do l = -buff_size+2,p+buff_size-3
                    do k = -buff_size+2, n+buff_size-3
                        do j = -buff_size+2, m+buff_size-3

                            if(viscous) then
                                dvel_L(j+1, k, l) = (1._wp/60._wp) * (-3._wp * dvel1(j-1, k, l) + &
                                                        27._wp * dvel1(j, k, l) + &
                                                        47._wp * dvel1(j+1, k, l) -   &
                                                        13._wp * dvel1(j+2, k, l) + &
                                                        2._wp * dvel1(j+3, k, l))

                                mu_L = 0._wp
                                !$acc loop seq
                                do q = 1, Re_size(1)
                                    mu_L =  qL_rs_vf(j+1,k,l,advxb+Re_idx(1, q)-1) / Res(1, q) &
                                              + mu_L
                                end do

                                flux_vf(momxb+2)%sf(j,k,l) = flux_vf(momxb+2)%sf(j,k,l) - 0.5_wp*mu_L*(dvel_L(j+1,k,l))
                                flux_vf(E_idx)%sf(j,k,l) = flux_vf(E_idx)%sf(j,k,l) - 0.5_wp*mu_L*qL_rs_vf(j+1,k,l,momxb)*dvel_L(j+1,k,l) 

                                dvel_L(j+1, k, l) = (1._wp/60._wp) * (-3._wp * dvel2(j-1, k, l) + &
                                                        27._wp * dvel2(j, k, l) + &
                                                        47._wp * dvel2(j+1, k, l) -   &
                                                        13._wp * dvel2(j+2, k, l) + &
                                                        2._wp * dvel2(j+3, k, l))

                                flux_vf(momxb+2)%sf(j,k,l) = flux_vf(momxb+2)%sf(j,k,l) - 0.5_wp*mu_L*(dvel_L(j+1,k,l))
                                flux_vf(E_idx)%sf(j,k,l) = flux_vf(E_idx)%sf(j,k,l) - 0.5_wp*mu_L*qL_rs_vf(j+1,k,l,momxb)*dvel_L(j+1,k,l)

                                if(num_fluids > 1) then 
                                    mu_L = 0._wp
                                    !$acc loop seq
                                    do q = 1, Re_size(1)
                                        mu_L =  qR_rs_vf(j,k,l,advxb+Re_idx(1, q)-1) / Res(1, q) &
                                                  + mu_L
                                    end do
                                end if

                                dvel_L(j, k, l) = (1._wp/60._wp) * (-3._wp * dvel1(j+2, k, l) + &
                                                        27._wp * dvel1(j+1, k, l) + &
                                                        47._wp * dvel1(j, k, l) -   &
                                                        13._wp * dvel1(j-1, k, l) + &
                                                        2._wp * dvel1(j-2, k, l))
                                

                                flux_vf(momxb+2)%sf(j,k,l) = flux_vf(momxb+2)%sf(j,k,l) - 0.5_wp*mu_L*(dvel_L(j,k,l))
                                flux_vf(E_idx)%sf(j,k,l) = flux_vf(E_idx)%sf(j,k,l) - 0.5_wp*mu_L*qR_rs_vf(j,k,l,1)*dvel_L(j,k,l) 

                                dvel_L(j, k, l) = (1._wp/60._wp) * (-3._wp * dvel2(j+2, k, l) + &
                                                        27._wp * dvel2(j+1, k, l) + &
                                                        47._wp * dvel2(j, k, l) -   &
                                                        13._wp * dvel2(j-1, k, l) + &
                                                        2._wp * dvel2(j-2, k, l))

                                flux_vf(momxb+2)%sf(j,k,l) = flux_vf(momxb+2)%sf(j,k,l) - 0.5_wp*mu_L*(dvel_L(j,k,l))
                                flux_vf(E_idx)%sf(j,k,l) = flux_vf(E_idx)%sf(j,k,l) - 0.5_wp*mu_L*qR_rs_vf(j,k,l,1)*dvel_L(j,k,l)

                            end if

                            !! duy & dvx
                            dvel1(j,k,l) = (1/(12._wp*dy(k))) * ( &
                            8._wp*q_prim_vf(momxb)%sf(j,k+1,l) - &
                            8._wp*q_prim_vf(momxb)%sf(j,k-1,l) + &
                            q_prim_vf(momxb)%sf(j,k-2,l) - &
                            q_prim_vf(momxb)%sf(j,k+2,l) )

                            dvel2(j,k,l) = (1/(12._wp*dx(j))) * ( &
                            8._wp*q_prim_vf(momxb+1)%sf(j+1,k,l) - &
                            8._wp*q_prim_vf(momxb+1)%sf(j-1,k,l) + &
                            q_prim_vf(momxb+1)%sf(j-2,k,l) - &
                            q_prim_vf(momxb+1)%sf(j+2,k,l) )

                            jac_rhs(j, k, l) = jac_rhs(j,k,l) + alf_igr* (2_wp * dvel1(j,k,l)* dvel2(j,k,l))
                        end do 
                    end do 
                end do 

                !$acc parallel loop collapse(3) gang vector default(present) private(rho_L,gamma_L,pi_inf_L,mu_L,a_L)
                do l = -buff_size+2,p+buff_size-3
                    do k = -buff_size+2, n+buff_size-3
                        do j = -buff_size+2, m+buff_size-3

                            if(viscous) then
                                dvel_L(j+1, k, l) = (1._wp/60._wp) * (-3._wp * dvel1(j-1, k, l) + &
                                                        27._wp * dvel1(j, k, l) + &
                                                        47._wp * dvel1(j+1, k, l) -   &
                                                        13._wp * dvel1(j+2, k, l) + &
                                                        2._wp * dvel1(j+3, k, l))

                                mu_L = 0._wp
                                !$acc loop seq
                                do q = 1, Re_size(1)
                                    mu_L =  qL_rs_vf(j+1,k,l,advxb+Re_idx(1, q)-1) / Res(1, q) &
                                              + mu_L
                                end do

                                flux_vf(momxb+1)%sf(j,k,l) = flux_vf(momxb+1)%sf(j,k,l) - 0.5_wp*mu_L*(dvel_L(j+1,k,l))
                                flux_vf(E_idx)%sf(j,k,l) = flux_vf(E_idx)%sf(j,k,l) - 0.5_wp*mu_L*qL_rs_vf(j+1,k,l,momxb)*dvel_L(j+1,k,l) 

                                dvel_L(j+1, k, l) = (1._wp/60._wp) * (-3._wp * dvel2(j-1, k, l) + &
                                                        27._wp * dvel2(j, k, l) + &
                                                        47._wp * dvel2(j+1, k, l) -   &
                                                        13._wp * dvel2(j+2, k, l) + &
                                                        2._wp * dvel2(j+3, k, l))

                                flux_vf(momxb+1)%sf(j,k,l) = flux_vf(momxb+1)%sf(j,k,l) - 0.5_wp*mu_L*(dvel_L(j+1,k,l))
                                flux_vf(E_idx)%sf(j,k,l) = flux_vf(E_idx)%sf(j,k,l) - 0.5_wp*mu_L*qL_rs_vf(j+1,k,l,momxb)*dvel_L(j+1,k,l)

                                if(num_fluids > 1) then 
                                    mu_L = 0._wp
                                    !$acc loop seq
                                    do q = 1, Re_size(1)
                                        mu_L =  qR_rs_vf(j,k,l,advxb+Re_idx(1, q)-1) / Res(1, q) &
                                                  + mu_L
                                    end do
                                end if

                                dvel_L(j, k, l) = (1._wp/60._wp) * (-3._wp * dvel1(j+2, k, l) + &
                                                        27._wp * dvel1(j+1, k, l) + &
                                                        47._wp * dvel1(j, k, l) -   &
                                                        13._wp * dvel1(j-1, k, l) + &
                                                        2._wp * dvel1(j-2, k, l))
                                

                                flux_vf(momxb+1)%sf(j,k,l) = flux_vf(momxb+1)%sf(j,k,l) - 0.5_wp*mu_L*(dvel_L(j,k,l))
                                flux_vf(E_idx)%sf(j,k,l) = flux_vf(E_idx)%sf(j,k,l) - 0.5_wp*mu_L*qR_rs_vf(j,k,l,1)*dvel_L(j,k,l) 

                                dvel_L(j, k, l) = (1._wp/60._wp) * (-3._wp * dvel2(j+2, k, l) + &
                                                        27._wp * dvel2(j+1, k, l) + &
                                                        47._wp * dvel2(j, k, l) -   &
                                                        13._wp * dvel2(j-1, k, l) + &
                                                        2._wp * dvel2(j-2, k, l))

                                flux_vf(momxb+1)%sf(j,k,l) = flux_vf(momxb+1)%sf(j,k,l) - 0.5_wp*mu_L*(dvel_L(j,k,l))
                                flux_vf(E_idx)%sf(j,k,l) = flux_vf(E_idx)%sf(j,k,l) - 0.5_wp*mu_L*qR_rs_vf(j,k,l,1)*dvel_L(j,k,l)

                            end if

                            !! dvz & dwy
                            dvel1(j,k,l) = (1/(12._wp*dz(l))) * ( &
                            8._wp*q_prim_vf(momxb+1)%sf(j,k,l+1) - &
                            8._wp*q_prim_vf(momxb+1)%sf(j,k,l-1) + &
                            q_prim_vf(momxb+1)%sf(j,k,l-2) - &
                            q_prim_vf(momxb+1)%sf(j,k,l+2) )

                            dvel2(j,k,l) = (1/(12._wp*dy(k))) * ( &
                            8._wp*q_prim_vf(momxb+2)%sf(j,k+1,l) - &
                            8._wp*q_prim_vf(momxb+2)%sf(j,k-1,l) + &
                            q_prim_vf(momxb+2)%sf(j,k-2,l) - &
                            q_prim_vf(momxb+2)%sf(j,k+2,l))

                            jac_rhs(j,k,l) = jac_rhs(j,k,l) + alf_igr*(2_wp*dvel1(j,k,l)*dvel2(j,k,l))

                            !! dux & dvy
                            dvel1(j,k,l) = (1/(12._wp*dx(j))) * ( &
                            8._wp*q_prim_vf(momxb)%sf(j+1,k,l) - &
                            8._wp*q_prim_vf(momxb)%sf(j-1,k,l) + &
                            q_prim_vf(momxb)%sf(j-2,k,l) - &
                            q_prim_vf(momxb)%sf(j+2,k,l) )

                            dvel2(j,k,l) = (1/(12._wp*dy(k))) * ( &
                            8._wp*q_prim_vf(momxb+1)%sf(j,k+1,l) - &
                            8._wp*q_prim_vf(momxb+1)%sf(j,k-1,l) + &
                            q_prim_vf(momxb+1)%sf(j,k-2,l) - &
                            q_prim_vf(momxb+1)%sf(j,k+2,l) )

                            jac_rhs(j, k, l) = jac_rhs(j,k,l) + alf_igr* (dvel1(j,k,l)**2_wp + dvel2(j,k,l)**2_wp)
                        end do 
                    end do 
                end do 

                !$acc parallel loop collapse(3) gang vector default(present) private(rho_L,gamma_L,pi_inf_L,mu_L,a_L)
                do l = -buff_size+2,p+buff_size-3
                    do k = -buff_size+2, n+buff_size-3
                        do j = -buff_size+2, m+buff_size-3

                            if(viscous) then

                                mu_L = 0._wp
                                !$acc loop seq
                                do q = 1, Re_size(1)
                                    mu_L =  qL_rs_vf(j+1,k,l,advxb+Re_idx(1, q)-1) / Res(1, q) &
                                              + mu_L
                                end do

                                dvel_L(j+1, k, l) = (1._wp/60._wp) * (-3._wp * dvel1(j-1, k, l) + &
                                                        27._wp * dvel1(j, k, l) + &
                                                        47._wp * dvel1(j+1, k, l) -   &
                                                        13._wp * dvel1(j+2, k, l) + &
                                                        2._wp * dvel1(j+3, k, l))

                                flux_vf(momxb)%sf(j,k,l) = flux_vf(momxb)%sf(j,k,l) - 0.5_wp*mu_L*((4._wp/3._wp)*dvel_L(j+1,k,l))
                                flux_vf(E_idx)%sf(j,k,l) = flux_vf(E_idx)%sf(j,k,l) - 0.5_wp*mu_L*qL_rs_vf(j+1,k,l,momxb)*(4._wp/3._wp)*dvel_L(j+1,k,l)

                                dvel_L(j+1, k, l) = (1._wp/60._wp) * (-3._wp * dvel2(j-1, k, l) + &
                                                        27._wp * dvel2(j, k, l) + &
                                                        47._wp * dvel2(j+1, k, l) -   &
                                                        13._wp * dvel2(j+2, k, l) + &
                                                        2._wp * dvel2(j+3, k, l))

                                flux_vf(momxb)%sf(j,k,l) = flux_vf(momxb)%sf(j,k,l) + 0.5_wp*mu_L*((2._wp/3._wp)*dvel_L(j+1,k,l))
                                flux_vf(E_idx)%sf(j,k,l) = flux_vf(E_idx)%sf(j,k,l) + 0.5_wp*mu_L*qL_rs_vf(j+1,k,l,momxb)*((2._wp/3._wp)*dvel_L(j+1,k,l))

                                if(num_fluids > 1) then 
                                    mu_L = 0._wp
                                    !$acc loop seq
                                    do q = 1, Re_size(1)
                                        mu_L =  qR_rs_vf(j,k,l,advxb+Re_idx(1, q)-1) / Res(1, q) &
                                                  + mu_L
                                    end do
                                end if

                                dvel_L(j, k, l) = (1._wp/60._wp) * (-3._wp * dvel1(j+2, k, l) + &
                                        27._wp * dvel1(j+1, k, l) + &
                                        47._wp * dvel1(j, k, l) -   &
                                        13._wp * dvel1(j-1, k, l) + &
                                        2._wp * dvel1(j-2, k, l))
                                

                                flux_vf(momxb)%sf(j,k,l) = flux_vf(momxb)%sf(j,k,l) - 0.5_wp*mu_L*((4._wp/3._wp)*dvel_L(j,k,l))
                                flux_vf(E_idx)%sf(j,k,l) = flux_vf(E_idx)%sf(j,k,l) - 0.5_wp*mu_L*qR_rs_vf(j,k,l,1)*(4._wp/3._wp)*dvel_L(j,k,l) 

                                dvel_L(j, k, l) = (1._wp/60._wp) * (-3._wp * dvel2(j+2, k, l) + &
                                                        27._wp * dvel2(j+1, k, l) + &
                                                        47._wp * dvel2(j, k, l) -   &
                                                        13._wp * dvel2(j-1, k, l) + &
                                                        2._wp * dvel2(j-2, k, l))

                                flux_vf(momxb)%sf(j,k,l) = flux_vf(momxb)%sf(j,k,l) + 0.5_wp*mu_L*((2._wp/3._wp)*dvel_L(j,k,l))
                                flux_vf(E_idx)%sf(j,k,l) = flux_vf(E_idx)%sf(j,k,l) + 0.5_wp*mu_L*qR_rs_vf(j,k,l,1)*(2._wp/3._wp)*dvel_L(j,k,l)

                            end if

                            !! dux + dvy, dwz

                            dvel1(j,k,l) = dvel1(j,k,l) + dvel2(j,k,l)

                            dvel2(j,k,l) = (1/(12._wp*dz(l))) * ( &
                            8._wp*q_prim_vf(momxb+2)%sf(j,k,l+1) - &
                            8._wp*q_prim_vf(momxb+2)%sf(j,k,l-1) + &
                            q_prim_vf(momxb+2)%sf(j,k,l-2) - &
                            q_prim_vf(momxb+2)%sf(j,k,l+2) )

                            jac_rhs(j,k,l) = jac_rhs(j,k,l) + alf_igr * (dvel2(j,k,l) **2_wp + (dvel1(j,k,l) + dvel2(j,k,l))**2_wp)
                        end do 
                    end do 
                end do 

                !$acc parallel loop collapse(3) gang vector default(present) private(rho_L,gamma_L,pi_inf_L,mu_L,a_L)
                do l = -buff_size+2,p+buff_size-3
                    do k = -buff_size+2, n+buff_size-3
                        do j = -buff_size+2, m+buff_size-3

                            if(viscous) then

                                mu_L = 0._wp
                                !$acc loop seq
                                do q = 1, Re_size(1)
                                    mu_L =  qL_rs_vf(j+1,k,l,advxb+Re_idx(1, q)-1) / Res(1, q) &
                                              + mu_L
                                end do

                                dvel_L(j+1, k, l) = (1._wp/60._wp) * (-3._wp * dvel2(j-1, k, l) + &
                                                        27._wp * dvel2(j, k, l) + &
                                                        47._wp * dvel2(j+1, k, l) -   &
                                                        13._wp * dvel2(j+2, k, l) + &
                                                        2._wp * dvel2(j+3, k, l))

                                flux_vf(momxb)%sf(j,k,l) = flux_vf(momxb)%sf(j,k,l) + 0.5_wp*mu_L*((2._wp/3._wp)*dvel_L(j+1,k,l))
                                flux_vf(E_idx)%sf(j,k,l) = flux_vf(E_idx)%sf(j,k,l) + 0.5_wp*mu_L*qL_rs_vf(j+1,k,l,momxb)*((2._wp/3._wp)*dvel_L(j+1,k,l))

                                if(num_fluids > 1) then 
                                    mu_L = 0._wp
                                    !$acc loop seq
                                    do q = 1, Re_size(1)
                                        mu_L =  qR_rs_vf(j,k,l,advxb+Re_idx(1, q)-1) / Res(1, q) &
                                                  + mu_L
                                    end do
                                end if

                                dvel_L(j, k, l) = (1._wp/60._wp) * (-3._wp * dvel2(j+2, k, l) + &
                                                        27._wp * dvel2(j+1, k, l) + &
                                                        47._wp * dvel2(j, k, l) -   &
                                                        13._wp * dvel2(j-1, k, l) + &
                                                        2._wp * dvel2(j-2, k, l))

                                flux_vf(momxb)%sf(j,k,l) = flux_vf(momxb)%sf(j,k,l) + 0.5_wp*mu_L*((2._wp/3._wp)*dvel_L(j,k,l))
                                flux_vf(E_idx)%sf(j,k,l) = flux_vf(E_idx)%sf(j,k,l) + 0.5_wp*mu_L*qR_rs_vf(j,k,l,1)*(2._wp/3._wp)*dvel_L(j,k,l)

                            end if

                            !$acc loop seq 
                            do i = 1, momxb - 1
                                qL_rs_vf(j, k, l, i) = (1._wp/60._wp) * (-3._wp * q_prim_vf(i)%sf(j+2, k, l) + &
                                                27._wp * q_prim_vf(i)%sf(j+1, k, l) + &
                                                47._wp * q_prim_vf(i)%sf(j, k, l) -   &
                                                13._wp * q_prim_vf(i)%sf(j-1, k, l) + &
                                                2._wp * q_prim_vf(i)%sf(j-2, k, l))
                            end do

                            !$acc loop seq 
                            do i = momxb+1,advxe
                                qL_rs_vf(j, k, l, i) = (1._wp/60._wp) * (-3._wp * q_prim_vf(i)%sf(j+2, k, l) + &
                                                27._wp * q_prim_vf(i)%sf(j+1, k, l) + &
                                                47._wp * q_prim_vf(i)%sf(j, k, l) -   &
                                                13._wp * q_prim_vf(i)%sf(j-1, k, l) + &
                                                2._wp * q_prim_vf(i)%sf(j-2, k, l))
                            end do 

                            if(sys_size > advxe) then 
                                do i = advxe + 1, sys_size
                                  qL_rs_vf(j, k, l, i) = (1._wp/60._wp) * (-3._wp * q_prim_vf(i)%sf(j+2, k, l) + &
                                                27._wp * q_prim_vf(i)%sf(j+1, k, l) + &
                                                47._wp * q_prim_vf(i)%sf(j, k, l) -   &
                                                13._wp * q_prim_vf(i)%sf(j-1, k, l) + &
                                                2._wp * q_prim_vf(i)%sf(j-2, k, l))
                                end do
                            end if  

                            rho_L = sum(qL_rs_vf(j,k,l,1:num_fluids))
                            gamma_L = sum(qL_rs_vf(j,k,l,advxb:advxb+num_fluids-1)*gammas)
                            pi_inf_L = sum(qL_rs_vf(j,k,l,advxb:advxb+num_fluids-1)*pi_infs)

                            !$acc loop seq
                            do i = 1, num_fluids
                                flux_vf(i)%sf(j,k,l) = flux_vf(i)%sf(j,k,l) + &
                                    0.5_wp * (qL_rs_vf(j,k,l,i) * &
                                    qR_rs_vf(j,k,l,1)) + &
                                    0.5_wp*cfl(j,k,l) * (qL_rs_vf(j,k,l,i))
                            end do

                            if(num_fluids > 1) then 
                                !$acc loop seq
                                do i = 1, num_fluids
                                    flux_vf(advxb+i-1)%sf(j,k,l) = &
                                        0.5_wp * (qR_rs_vf(j,k,l,advxb+i-1) * &
                                        qR_rs_vf(j,k,l,1)) + &
                                        0.5_wp*cfl(j,k,l) * (qR_rs_vf(j,k,l,advxb+i-1))

                                    flux_vf(advxb+i-1)%sf(j,k,l) = flux_vf(advxb+i-1)%sf(j,k,l) &
                                    - 0.5_wp * q_prim_vf(advxb+i-1)%sf(j,k,l) * qR_rs_vf(j,k,l,1)
                                end do
                            end if


                            ! Momentum -> rho*u^2 + p + [[F_igr]]
                            flux_vf(momxb)%sf(j,k,l) = flux_vf(momxb)%sf(j,k,l) + &
                                 0.5_wp* (rho_L * (qR_rs_vf(j,k,l,1))**2.0 + &
                                 qL_rs_vf(j,k,l,E_idx)) + &
                                 0.5_wp*cfl(j,k,l) * (qR_rs_vf(j,k,l,1))

                            flux_vf(momxb+1)%sf(j, k, l) = flux_vf(momxb+1)%sf(j, k, l) + &
                                0.5_wp * rho_L * qR_rs_vf(j,k,l,1)*qL_rs_vf(j,k,l,momxb+1) + &
                                0.5_wp*cfl(j,k,l) * (qL_rs_vf(j,k,l,momxb+1))

                            flux_vf(momxb+2)%sf(j, k, l) = flux_vf(momxb+2)%sf(j, k, l) + &
                                0.5_wp * rho_L * qR_rs_vf(j,k,l,1)*qL_rs_vf(j,k,l,momxb+2) + &
                                0.5_wp*cfl(j,k,l) * (qL_rs_vf(j,k,l,momxb+2))

                             flux_vf(E_idx)%sf(j, k, l) = flux_vf(E_idx)%sf(j, k, l) +  &
                                0.5_wp * (qR_rs_vf(j,k,l,1) * (qL_rs_vf(j,k,l,E_idx)*gamma_L + pi_inf_L + &
                                0.5_wp * rho_L * (qR_rs_vf(j,k,l,1)**2._wp + qL_rs_vf(j,k,l,momxb+1)**2._wp + qL_rs_vf(j,k,l,momxb+2)**2._wp ) + &
                                qL_rs_vf(j,k,l,E_idx)) ) + &
                                0.5_wp*cfl(j,k,l) * (qL_rs_vf(j,k,l,E_idx))

                        end do
                    end do
                end do
            end if
        else if (idir == 2) then
            if(p == 0) then
                !$acc parallel loop collapse(3) gang vector default(present) private(rho_L,gamma_L,pi_inf_L,mu_L,a_L)
                do l = 0, p
                    do k = -buff_size+2, n+buff_size-3
                        do j = -buff_size+2, m+buff_size-3

                            !$acc loop seq
                            do i = 1, sys_size
                                qL_rs_vf(j, k+1, l, i) = (1._wp/60._wp) * (-3._wp * q_prim_vf(i)%sf(j, k-1, l) + &
                                                27._wp * q_prim_vf(i)%sf(j, k, l) + &
                                                47._wp * q_prim_vf(i)%sf(j, k+1, l) -   &
                                                13._wp * q_prim_vf(i)%sf(j, k+2, l) + &
                                                2._wp * q_prim_vf(i)%sf(j, k+3, l))
                            end do

                            qR_rs_vf(j, k, l, 1) = (1._wp/60._wp) * (-3._wp * q_prim_vf(momxb+1)%sf(j, k+2, l) + &
                                                27._wp * q_prim_vf(momxb+1)%sf(j, k+1, l) + &
                                                47._wp * q_prim_vf(momxb+1)%sf(j, k, l) -   &
                                                13._wp * q_prim_vf(momxb+1)%sf(j, k-1, l) + &
                                                2._wp * q_prim_vf(momxb+1)%sf(j, k-2, l))
                            
                            if(num_fluids > 1) then  
                                !$acc loop seq 
                                do i = advxb, advxe                           
                                    qR_rs_vf(j, k, l,i) = (1._wp/60._wp) * (-3._wp * q_prim_vf(i)%sf(j, k+2, l) + &
                                                        27._wp * q_prim_vf(i)%sf(j, k+1, l) + &
                                                        47._wp * q_prim_vf(i)%sf(j, k, l) -   &
                                                        13._wp * q_prim_vf(i)%sf(j, k-1, l) + &
                                                        2._wp * q_prim_vf(i)%sf(j, k-2, l))
                                end do 
                            end if

                            rho_L = sum(qL_rs_vf(j,k+1,l,1:num_fluids))
                            gamma_L = sum(qL_rs_vf(j,k+1,l,advxb:advxb+num_fluids-1)*gammas)
                            pi_inf_L = sum(qL_rs_vf(j,k+1,l,advxb:advxb+num_fluids-1)*pi_infs)

                            a_L = sqrt((qL_rs_vf(j,k+1,l,E_idx)*(1._wp/gamma_L + 1._wp) + pi_inf_L / gamma_L) / rho_L)

                            cfl(j,k,l) = 2_wp*(sqrt(qL_rs_vf(j,k+1,l,momxb)**2._wp + qL_rs_vf(j,k+1,l,momxb+1)**2._wp) + a_L)

                            !$acc loop seq
                            do i = 1, num_fluids
                                flux_vf(i)%sf(j,k,l) = &
                                    0.5_wp * (qL_rs_vf(j,k+1,l,i) * &
                                    qL_rs_vf(j,k+1,l,momxb+1)) - &
                                    0.5_wp*cfl(j,k,l) * (qL_rs_vf(j,k+1,l,i))
                            end do

                            if(num_fluids > 1) then 
                                !$acc loop seq
                                do i = 1, num_fluids
                                    flux_vf(advxb+i-1)%sf(j,k,l) = &
                                        0.5_wp * (qL_rs_vf(j,k+1,l,advxb+i-1) * &
                                        qL_rs_vf(j,k+1,l,momxb+1 )) - &
                                        0.5_wp*cfl(j,k,l) * (qL_rs_vf(j,k+1,l,advxb+i-1))

                                    flux_vf(advxb+i-1)%sf(j,k,l) = flux_vf(advxb+i-1)%sf(j,k,l) &
                                    - 0.5_wp * q_prim_vf(advxb+i-1)%sf(j,k,l) * qL_rs_vf(j,k+1,l,momxb+1)
                                end do
                            end if

                            ! Momentum -> rho*u^2 + p + [[F_igr]]
                            flux_vf(momxb+1)%sf(j,k,l) = &
                                 0.5_wp* (rho_L * (qL_rs_vf(j,k+1,l,momxb+1))**2.0 + &
                                 qL_rs_vf(j,k+1,l,E_idx)) - &
                                 0.5_wp*cfl(j,k,l) * (qL_rs_vf(j,k+1,l,momxb+1))

                            flux_vf(momxb)%sf(j, k, l) =  &
                                0.5_wp * rho_L * qL_rs_vf(j,k+1,l,momxb)*qL_rs_vf(j,k+1,l,momxb+1) - &
                                0.5_wp*cfl(j,k,l) * (qL_rs_vf(j,k+1,l,momxb))

                             flux_vf(E_idx)%sf(j, k, l) = &
                                0.5_wp * (qL_rs_vf(j,k+1,l,momxb+1) * (qL_rs_vf(j,k+1,l,E_idx)*gamma_L + pi_inf_L + &
                                0.5_wp * rho_L * (qL_rs_vf(j,k+1,l,momxb)**2._wp + qL_rs_vf(j,k+1,l,momxb+1)**2._wp ) + &
                                qL_rs_vf(j,k+1,l,E_idx)) ) - &
                                0.5_wp*cfl(j,k,l) * (qL_rs_vf(j,k+1,l,E_idx))

                            !! duy & dvx
                            dvel1(j,k,l) = (1/(12._wp*dy(k))) * ( &
                            8._wp*q_prim_vf(momxb)%sf(j,k+1,l) - &
                            8._wp*q_prim_vf(momxb)%sf(j,k-1,l) + &
                            q_prim_vf(momxb)%sf(j,k-2,l) - &
                            q_prim_vf(momxb)%sf(j,k+2,l) )

                            dvel2(j,k,l) = (1/(12._wp*dx(j))) * ( &
                            8._wp*q_prim_vf(momxb+1)%sf(j+1,k,l) - &
                            8._wp*q_prim_vf(momxb+1)%sf(j-1,k,l) + &
                            q_prim_vf(momxb+1)%sf(j-2,k,l) - &
                            q_prim_vf(momxb+1)%sf(j+2,k,l) )
                        end do 
                    end do 
                end do

                !$acc parallel loop collapse(3) gang vector default(present) private(rho_L,gamma_L,pi_inf_L,mu_L,a_L)
                do l = 0, p
                    do k = -buff_size+2, n+buff_size-3
                        do j = -buff_size+2, m+buff_size-3

                            if(viscous) then
                                dvel_L(j, k+1, l) = (1._wp/60._wp) * (-3._wp * dvel1(j, k-1, l) + &
                                                        27._wp * dvel1(j, k, l) + &
                                                        47._wp * dvel1(j, k+1, l) -   &
                                                        13._wp * dvel1(j, k+2, l) + &
                                                        2._wp * dvel1(j, k+3, l))

                                mu_L = 0._wp
                                !$acc loop seq
                                do q = 1, Re_size(1)
                                    mu_L =  qL_rs_vf(j,k+1,l,advxb+Re_idx(1, q)-1) / Res(1, q) &
                                              + mu_L
                                end do

                                flux_vf(momxb)%sf(j,k,l) = flux_vf(momxb)%sf(j,k,l) - 0.5_wp*mu_L*(dvel_L(j,k+1,l))
                                flux_vf(E_idx)%sf(j,k,l) = flux_vf(E_idx)%sf(j,k,l) - 0.5_wp*mu_L*qL_rs_vf(j,k+1,l,momxb+1)*dvel_L(j,k+1,l) 

                                dvel_L(j, k+1, l) = (1._wp/60._wp) * (-3._wp * dvel2(j, k-1, l) + &
                                                        27._wp * dvel2(j, k, l) + &
                                                        47._wp * dvel2(j, k+1, l) -   &
                                                        13._wp * dvel2(j, k+2, l) + &
                                                        2._wp * dvel2(j, k+3, l))

                                flux_vf(momxb)%sf(j,k,l) = flux_vf(momxb)%sf(j,k,l) - 0.5_wp*mu_L*(dvel_L(j,k+1,l))
                                flux_vf(E_idx)%sf(j,k,l) = flux_vf(E_idx)%sf(j,k,l) - 0.5_wp*mu_L*qL_rs_vf(j,k+1,l,momxb+1)*dvel_L(j,k+1,l)

                                if(num_fluids > 1) then 
                                    mu_L = 0._wp
                                    !$acc loop seq
                                    do q = 1, Re_size(1)
                                        mu_L =  qR_rs_vf(j,k,l,advxb+Re_idx(1, q)-1) / Res(1, q) &
                                                  + mu_L
                                    end do
                                end if

                                dvel_L(j, k, l) = (1._wp/60._wp) * (-3._wp * dvel1(j, k+2, l) + &
                                                        27._wp * dvel1(j, k+1, l) + &
                                                        47._wp * dvel1(j, k, l) -   &
                                                        13._wp * dvel1(j, k-1, l) + &
                                                        2._wp * dvel1(j, k-2, l))
                                

                                flux_vf(momxb)%sf(j,k,l) = flux_vf(momxb)%sf(j,k,l) - 0.5_wp*mu_L*(dvel_L(j,k,l))
                                flux_vf(E_idx)%sf(j,k,l) = flux_vf(E_idx)%sf(j,k,l) - 0.5_wp*mu_L*qR_rs_vf(j,k,l,1)*dvel_L(j,k,l) 

                                dvel_L(j, k, l) = (1._wp/60._wp) * (-3._wp * dvel2(j, k+2, l) + &
                                                        27._wp * dvel2(j, k+1, l) + &
                                                        47._wp * dvel2(j, k, l) -   &
                                                        13._wp * dvel2(j, k-1, l) + &
                                                        2._wp * dvel2(j, k-2, l))

                                flux_vf(momxb)%sf(j,k,l) = flux_vf(momxb)%sf(j,k,l) - 0.5_wp*mu_L*(dvel_L(j,k,l))
                                flux_vf(E_idx)%sf(j,k,l) = flux_vf(E_idx)%sf(j,k,l) - 0.5_wp*mu_L*qR_rs_vf(j,k,l,1)*dvel_L(j,k,l)

                            end if


                            !! dvy & dux
                            dvel1(j,k,l) = (1/(12._wp*dy(k))) * ( &
                            8._wp*q_prim_vf(momxb+1)%sf(j,k+1,l) - &
                            8._wp*q_prim_vf(momxb+1)%sf(j,k-1,l) + &
                            q_prim_vf(momxb+1)%sf(j,k-2,l) - &
                            q_prim_vf(momxb+1)%sf(j,k+2,l) )

                            dvel2(j,k,l) = (1/(12._wp*dx(j))) * ( &
                            8._wp*q_prim_vf(momxb)%sf(j+1,k,l) - &
                            8._wp*q_prim_vf(momxb)%sf(j-1,k,l) + &
                            q_prim_vf(momxb)%sf(j-2,k,l) - &
                            q_prim_vf(momxb)%sf(j+2,k,l) )

                        end do 
                    end do 
                end do

                !$acc parallel loop collapse(3) gang vector default(present) private(rho_L,gamma_L,pi_inf_L,mu_L,a_L)
                do l = 0, p
                    do k = -buff_size+2, n+buff_size-3
                        do j = -buff_size+2, m+buff_size-3

                            if(viscous) then
                                dvel_L(j, k+1, l) = (1._wp/60._wp) * (-3._wp * dvel1(j, k-1, l) + &
                                                        27._wp * dvel1(j, k, l) + &
                                                        47._wp * dvel1(j, k+1, l) -   &
                                                        13._wp * dvel1(j, k+2, l) + &
                                                        2._wp * dvel1(j, k+3, l))


                                mu_L = 0._wp
                                !$acc loop seq
                                do q = 1, Re_size(1)
                                    mu_L =  qL_rs_vf(j,k+1,l,advxb+Re_idx(1, q)-1) / Res(1, q) &
                                              + mu_L
                                end do
                                

                                flux_vf(momxb+1)%sf(j,k,l) = flux_vf(momxb+1)%sf(j,k,l) - 0.5_wp*mu_L*((4._wp/3._wp)*dvel_L(j,k+1,l))
                                flux_vf(E_idx)%sf(j,k,l) = flux_vf(E_idx)%sf(j,k,l) - 0.5_wp*mu_L*(4._wp/3._wp)*qL_rs_vf(j,k+1,l,momxb+1)*dvel_L(j,k+1,l) 

                                dvel_L(j, k+1, l) = (1._wp/60._wp) * (-3._wp * dvel2(j, k-1, l) + &
                                                        27._wp * dvel2(j, k, l) + &
                                                        47._wp * dvel2(j, k+1, l) -   &
                                                        13._wp * dvel2(j, k+2, l) + &
                                                        2._wp * dvel2(j, k+3, l))

                                flux_vf(momxb+1)%sf(j,k,l) = flux_vf(momxb+1)%sf(j,k,l) + 0.5_wp*mu_L*((2._wp/3._wp)*dvel_L(j,k+1,l))
                                flux_vf(E_idx)%sf(j,k,l) = flux_vf(E_idx)%sf(j,k,l) + 0.5_wp*mu_L*(2._wp/3._wp)*qL_rs_vf(j,k+1,l,momxb+1)*dvel_L(j,k+1,l) 

                                if(num_fluids > 1) then 
                                    mu_L = 0._wp
                                    !$acc loop seq
                                    do q = 1, Re_size(1)
                                        mu_L =  qR_rs_vf(j,k,l,advxb+Re_idx(1, q)-1) / Res(1, q) &
                                                  + mu_L
                                    end do
                                end if

                                dvel_L(j, k, l) = (1._wp/60._wp) * (-3._wp * dvel1(j, k+2, l) + &
                                                        27._wp * dvel1(j, k+1, l) + &
                                                        47._wp * dvel1(j, k, l) -   &
                                                        13._wp * dvel1(j, k-1, l) + &
                                                        2._wp * dvel1(j, k-2, l))
                                

                                flux_vf(momxb+1)%sf(j,k,l) = flux_vf(momxb+1)%sf(j,k,l) - 0.5_wp*mu_L*(4._wp/3._wp)*(dvel_L(j,k,l))
                                flux_vf(E_idx)%sf(j,k,l) = flux_vf(E_idx)%sf(j,k,l) - 0.5_wp*mu_L*(4._wp/3._wp)*qR_rs_vf(j,k,l,1)*dvel_L(j,k,l) 

                                dvel_L(j, k, l) = (1._wp/60._wp) * (-3._wp * dvel2(j, k+2, l) + &
                                                        27._wp * dvel2(j, k+1, l) + &
                                                        47._wp * dvel2(j, k, l) -   &
                                                        13._wp * dvel2(j, k-1, l) + &
                                                        2._wp * dvel2(j, k-2, l))

                                flux_vf(momxb+1)%sf(j,k,l) = flux_vf(momxb+1)%sf(j,k,l) + 0.5_wp*mu_L*((2._wp/3._wp)*dvel_L(j,k,l))
                                flux_vf(E_idx)%sf(j,k,l) = flux_vf(E_idx)%sf(j,k,l) + 0.5_wp*mu_L*(2._wp/3._wp)*qR_rs_vf(j,k,l,1)*dvel_L(j,k,l)

                            end if

                            !$acc loop seq 
                            do i = 1, momxb 
                                qL_rs_vf(j, k, l, i) = (1._wp/60._wp) * (-3._wp * q_prim_vf(i)%sf(j, k+2, l) + &
                                                27._wp * q_prim_vf(i)%sf(j, k+1, l) + &
                                                47._wp * q_prim_vf(i)%sf(j, k, l) -   &
                                                13._wp * q_prim_vf(i)%sf(j, k-1, l) + &
                                                2._wp * q_prim_vf(i)%sf(j, k-2, l))
                            end do

                            !$acc loop seq 
                            do i = momxb+2,advxe
                                qL_rs_vf(j, k, l, i) = (1._wp/60._wp) * (-3._wp * q_prim_vf(i)%sf(j, k+2, l) + &
                                                27._wp * q_prim_vf(i)%sf(j, k+1, l) + &
                                                47._wp * q_prim_vf(i)%sf(j, k, l) -   &
                                                13._wp * q_prim_vf(i)%sf(j, k-1, l) + &
                                                2._wp * q_prim_vf(i)%sf(j, k-2, l))
                            end do 

                            if(sys_size > advxe) then 
                                do i = advxe + 1, sys_size
                                    qL_rs_vf(j, k, l, i) = (1._wp/60._wp) * (-3._wp * q_prim_vf(i)%sf(j, k+2, l) + &
                                                    27._wp * q_prim_vf(i)%sf(j, k+1, l) + &
                                                    47._wp * q_prim_vf(i)%sf(j, k, l) -   &
                                                    13._wp * q_prim_vf(i)%sf(j, k-1, l) + &
                                                    2._wp * q_prim_vf(i)%sf(j, k-2, l))
                                end do
                            end if  

                            !$acc loop seq
                            do i = 1, num_fluids
                                flux_vf(i)%sf(j,k,l) = flux_vf(i)%sf(j,k,l) + &
                                    0.5_wp * (qL_rs_vf(j,k,l,i) * &
                                    qR_rs_vf(j,k,l,1)) + &
                                    0.5_wp*cfl(j,k,l) * (qL_rs_vf(j,k,l,i))
                            end do

                            rho_L = sum(qL_rs_vf(j,k,l,1:num_fluids))
                            gamma_L = sum(qL_rs_vf(j,k,l,advxb:advxb+num_fluids-1)*gammas)
                            pi_inf_L = sum(qL_rs_vf(j,k,l,advxb:advxb+num_fluids-1)*pi_infs)

                            if(num_fluids > 1) then 
                                !$acc loop seq
                                do i = 1, num_fluids
                                    flux_vf(advxb+i-1)%sf(j,k,l) = &
                                        0.5_wp * (qR_rs_vf(j,k,l,advxb+i-1) * &
                                        qR_rs_vf(j,k,l,1)) + &
                                        0.5_wp*cfl(j,k,l) * (qR_rs_vf(j,k,l,advxb+i-1))

                                    flux_vf(advxb+i-1)%sf(j,k,l) = flux_vf(advxb+i-1)%sf(j,k,l) &
                                    - 0.5_wp * q_prim_vf(advxb+i-1)%sf(j,k,l) * qR_rs_vf(j,k,l,1)
                                end do
                            end if

                            ! Momentum -> rho*u^2 + p + [[F_igr]]
                            flux_vf(momxb+1)%sf(j,k,l) = flux_vf(momxb+1)%sf(j,k,l) + &
                                 0.5_wp* (rho_L * (qR_rs_vf(j,k,l,1))**2.0 + &
                                 qL_rs_vf(j,k,l,E_idx)) + &
                                 0.5_wp*cfl(j,k,l) * (qR_rs_vf(j,k,l,1))

                            flux_vf(momxb)%sf(j, k, l) = flux_vf(momxb)%sf(j, k, l) +  &
                                0.5_wp * rho_L * qL_rs_vf(j,k,l,momxb)*qR_rs_vf(j,k,l,1) + &
                                0.5_wp*cfl(j,k,l) * (qL_rs_vf(j,k,l,momxb))

                             flux_vf(E_idx)%sf(j, k, l) = flux_vf(E_idx)%sf(j, k, l) + &
                                0.5_wp * (qR_rs_vf(j,k,l,1) * (qL_rs_vf(j,k,l,E_idx)*gamma_L + pi_inf_L + &
                                0.5_wp * rho_L * (qL_rs_vf(j,k,l,momxb)**2._wp + qR_rs_vf(j,k,l,1)**2._wp ) + &
                                qL_rs_vf(j,k,l,E_idx)) ) + &
                                0.5_wp*cfl(j,k,l) * (qL_rs_vf(j,k,l,E_idx))

                        end do
                    end do
                end do
            else
                !$acc parallel loop collapse(3) gang vector default(present) private(rho_L,gamma_L,pi_inf_L,mu_L,a_L)
                do l = -buff_size+2, p+buff_size-3
                    do k = -buff_size+2, n+buff_size-3
                        do j = -buff_size+2, m+buff_size-3

                            !$acc loop seq
                            do i = 1, sys_size
                                qL_rs_vf(j, k+1, l, i) = (1._wp/60._wp) * (-3._wp * q_prim_vf(i)%sf(j, k-1, l) + &
                                                27._wp * q_prim_vf(i)%sf(j, k, l) + &
                                                47._wp * q_prim_vf(i)%sf(j, k+1, l) -   &
                                                13._wp * q_prim_vf(i)%sf(j, k+2, l) + &
                                                2._wp * q_prim_vf(i)%sf(j, k+3, l))
                            end do

                            qR_rs_vf(j, k, l, 1) = (1._wp/60._wp) * (-3._wp * q_prim_vf(momxb+1)%sf(j, k+2, l) + &
                                                27._wp * q_prim_vf(momxb+1)%sf(j, k+1, l) + &
                                                47._wp * q_prim_vf(momxb+1)%sf(j, k, l) -   &
                                                13._wp * q_prim_vf(momxb+1)%sf(j, k-1, l) + &
                                                2._wp * q_prim_vf(momxb+1)%sf(j, k-2, l))
                            
                            if(num_fluids > 1) then  
                                !$acc loop seq 
                                do i = advxb, advxe                           
                                    qR_rs_vf(j, k, l,i) = (1._wp/60._wp) * (-3._wp * q_prim_vf(i)%sf(j, k+2, l) + &
                                                        27._wp * q_prim_vf(i)%sf(j, k+1, l) + &
                                                        47._wp * q_prim_vf(i)%sf(j, k, l) -   &
                                                        13._wp * q_prim_vf(i)%sf(j, k-1, l) + &
                                                        2._wp * q_prim_vf(i)%sf(j, k-2, l))
                                end do 
                            end if

                            rho_L = sum(qL_rs_vf(j,k+1,l,1:num_fluids))
                            gamma_L = sum(qL_rs_vf(j,k+1,l,advxb:advxb+num_fluids-1)*gammas)
                            pi_inf_L = sum(qL_rs_vf(j,k+1,l,advxb:advxb+num_fluids-1)*pi_infs)

                            a_L = sqrt((qL_rs_vf(j,k+1,l,E_idx)*(1._wp/gamma_L + 1._wp) + pi_inf_L / gamma_L) / rho_L)

                            cfl(j,k,l) = (sqrt(qL_rs_vf(j,k+1,l,momxb)**2._wp + qL_rs_vf(j,k+1,l,momxb+1)**2._wp + qL_rs_vf(j,k+1,l,momxb+2)**2._wp) + a_L)

                            !$acc loop seq
                            do i = 1, num_fluids
                                flux_vf(i)%sf(j,k,l) = &
                                    0.5_wp * (qL_rs_vf(j,k+1,l,i) * &
                                    qL_rs_vf(j,k+1,l,momxb+1)) - &
                                    0.5_wp*cfl(j,k,l) * (qL_rs_vf(j,k+1,l,i))
                            end do

                            if(num_fluids > 1) then 
                                !$acc loop seq
                                do i = 1, num_fluids
                                    flux_vf(advxb+i-1)%sf(j,k,l) = &
                                        0.5_wp * (qL_rs_vf(j,k+1,l,advxb+i-1) * &
                                        qL_rs_vf(j,k+1,l,momxb+1 )) - &
                                        0.5_wp*cfl(j,k,l) * (qL_rs_vf(j,k+1,l,advxb+i-1))

                                    flux_vf(advxb+i-1)%sf(j,k,l) = flux_vf(advxb+i-1)%sf(j,k,l) &
                                    - 0.5_wp * q_prim_vf(advxb+i-1)%sf(j,k,l) * qL_rs_vf(j,k+1,l,momxb+1)
                                end do
                            end if

                            ! Momentum -> rho*u^2 + p + [[F_igr]]
                            flux_vf(momxb+1)%sf(j,k,l) = &
                                 0.5_wp* (rho_L * (qL_rs_vf(j,k+1,l,momxb+1))**2.0 + &
                                 qL_rs_vf(j,k+1,l,E_idx)) - &
                                 0.5_wp*cfl(j,k,l) * (qL_rs_vf(j,k+1,l,momxb+1))

                            flux_vf(momxb)%sf(j, k, l) =  &
                                0.5_wp * rho_L * qL_rs_vf(j,k+1,l,momxb)*qL_rs_vf(j,k+1,l,momxb+1) - &
                                0.5_wp*cfl(j,k,l) * (qL_rs_vf(j,k+1,l,momxb))

                            flux_vf(momxb+2)%sf(j, k, l) =  &
                                0.5_wp * rho_L * qL_rs_vf(j,k+1,l,momxb+2)*qL_rs_vf(j,k+1,l,momxb+1) - &
                                0.5_wp*cfl(j,k,l) * (qL_rs_vf(j,k+1,l,momxb+2))

                             flux_vf(E_idx)%sf(j, k, l) = &
                                0.5_wp * (qL_rs_vf(j,k+1,l,momxb+1) * (qL_rs_vf(j,k+1,l,E_idx)*gamma_L + pi_inf_L + &
                                0.5_wp * rho_L * (qL_rs_vf(j,k+1,l,momxb)**2._wp + qL_rs_vf(j,k+1,l,momxb+1)**2._wp +qL_rs_vf(j,k+1,l,momxb+2)**2._wp) + &
                                qL_rs_vf(j,k+1,l,E_idx)) ) - &
                                0.5_wp*cfl(j,k,l) * (qL_rs_vf(j,k+1,l,E_idx))

                            !! duy & dvx
                            dvel1(j,k,l) = (1/(12._wp*dy(k))) * ( &
                            8._wp*q_prim_vf(momxb)%sf(j,k+1,l) - &
                            8._wp*q_prim_vf(momxb)%sf(j,k-1,l) + &
                            q_prim_vf(momxb)%sf(j,k-2,l) - &
                            q_prim_vf(momxb)%sf(j,k+2,l) )

                            dvel2(j,k,l) = (1/(12._wp*dx(j))) * ( &
                            8._wp*q_prim_vf(momxb+1)%sf(j+1,k,l) - &
                            8._wp*q_prim_vf(momxb+1)%sf(j-1,k,l) + &
                            q_prim_vf(momxb+1)%sf(j-2,k,l) - &
                            q_prim_vf(momxb+1)%sf(j+2,k,l) )
                        end do 
                    end do 
                end do 

                !$acc parallel loop collapse(3) gang vector default(present) private(rho_L,gamma_L,pi_inf_L,mu_L,a_L)
                do l = -buff_size+2,p+buff_size-3
                    do k = -buff_size+2, n+buff_size-3
                        do j = -buff_size+2, m+buff_size-3

                            if(viscous) then
                                dvel_L(j, k+1, l) = (1._wp/60._wp) * (-3._wp * dvel1(j, k-1, l) + &
                                                        27._wp * dvel1(j, k, l) + &
                                                        47._wp * dvel1(j, k+1, l) -   &
                                                        13._wp * dvel1(j, k+2, l) + &
                                                        2._wp * dvel1(j, k+3, l))

                                mu_L = 0._wp
                                !$acc loop seq
                                do q = 1, Re_size(1)
                                    mu_L =  qL_rs_vf(j,k+1,l,advxb+Re_idx(1, q)-1) / Res(1, q) &
                                              + mu_L
                                end do
                                

                                flux_vf(momxb)%sf(j,k,l) = flux_vf(momxb)%sf(j,k,l) - 0.5_wp*mu_L*(dvel_L(j,k+1,l))
                                flux_vf(E_idx)%sf(j,k,l) = flux_vf(E_idx)%sf(j,k,l) - 0.5_wp*mu_L*qL_rs_vf(j,k+1,l,momxb+1)*dvel_L(j,k+1,l) 

                                dvel_L(j, k+1, l) = (1._wp/60._wp) * (-3._wp * dvel2(j, k-1, l) + &
                                                        27._wp * dvel2(j, k, l) + &
                                                        47._wp * dvel2(j, k+1, l) -   &
                                                        13._wp * dvel2(j, k+2, l) + &
                                                        2._wp * dvel2(j, k+3, l))

                                flux_vf(momxb)%sf(j,k,l) = flux_vf(momxb)%sf(j,k,l) - 0.5_wp*mu_L*(dvel_L(j,k+1,l))
                                flux_vf(E_idx)%sf(j,k,l) = flux_vf(E_idx)%sf(j,k,l) - 0.5_wp*mu_L*qL_rs_vf(j,k+1,l,momxb+1)*dvel_L(j,k+1,l)

                                if(num_fluids > 1) then 
                                    mu_L = 0._wp
                                    !$acc loop seq
                                    do q = 1, Re_size(1)
                                        mu_L =  qR_rs_vf(j,k,l,advxb+Re_idx(1, q)-1) / Res(1, q) &
                                                  + mu_L
                                    end do
                                end if

                                dvel_L(j, k, l) = (1._wp/60._wp) * (-3._wp * dvel1(j, k+2, l) + &
                                                        27._wp * dvel1(j, k+1, l) + &
                                                        47._wp * dvel1(j, k, l) -   &
                                                        13._wp * dvel1(j, k-1, l) + &
                                                        2._wp * dvel1(j, k-2, l))
                                

                                flux_vf(momxb)%sf(j,k,l) = flux_vf(momxb)%sf(j,k,l) - 0.5_wp*mu_L*(dvel_L(j,k,l))
                                flux_vf(E_idx)%sf(j,k,l) = flux_vf(E_idx)%sf(j,k,l) - 0.5_wp*mu_L*qR_rs_vf(j,k,l,1)*dvel_L(j,k,l) 

                                dvel_L(j, k, l) = (1._wp/60._wp) * (-3._wp * dvel2(j, k+2, l) + &
                                                        27._wp * dvel2(j, k+1, l) + &
                                                        47._wp * dvel2(j, k, l) -   &
                                                        13._wp * dvel2(j, k-1, l) + &
                                                        2._wp * dvel2(j, k-2, l))

                                flux_vf(momxb)%sf(j,k,l) = flux_vf(momxb)%sf(j,k,l) - 0.5_wp*mu_L*(dvel_L(j,k,l))
                                flux_vf(E_idx)%sf(j,k,l) = flux_vf(E_idx)%sf(j,k,l) - 0.5_wp*mu_L*qR_rs_vf(j,k,l,1)*dvel_L(j,k,l)

                            end if

                            !! dwy & dvz
                            dvel1(j,k,l) = (1/(12._wp*dy(k))) * ( &
                            8._wp*q_prim_vf(momxb+2)%sf(j,k+1,l) - &
                            8._wp*q_prim_vf(momxb+2)%sf(j,k-1,l) + &
                            q_prim_vf(momxb+2)%sf(j,k-2,l) - &
                            q_prim_vf(momxb+2)%sf(j,k+2,l) )

                            dvel2(j,k,l) = (1/(12._wp*dz(l))) * ( &
                            8._wp*q_prim_vf(momxb+1)%sf(j,k,l+1) - &
                            8._wp*q_prim_vf(momxb+1)%sf(j,k,l-1) + &
                            q_prim_vf(momxb+1)%sf(j,k,l-2) - &
                            q_prim_vf(momxb+1)%sf(j,k,l+2) )
                        end do 
                    end do 
                end do 

                !$acc parallel loop collapse(3) gang vector default(present) private(rho_L,gamma_L,pi_inf_L,mu_L,a_L)
                do l = -buff_size+2,p+buff_size-3
                    do k = -buff_size+2, n+buff_size-3
                        do j = -buff_size+2, m+buff_size-3

                            if(viscous) then

                                dvel_L(j, k+1, l) = (1._wp/60._wp) * (-3._wp * dvel1(j, k-1, l) + &
                                                        27._wp * dvel1(j, k, l) + &
                                                        47._wp * dvel1(j, k+1, l) -   &
                                                        13._wp * dvel1(j, k+2, l) + &
                                                        2._wp * dvel1(j, k+3, l))

                                mu_L = 0._wp
                                !$acc loop seq
                                do q = 1, Re_size(1)
                                    mu_L =  qL_rs_vf(j,k+1,l,advxb+Re_idx(1, q)-1) / Res(1, q) &
                                              + mu_L
                                end do
                                

                                flux_vf(momxb+2)%sf(j,k,l) = flux_vf(momxb+2)%sf(j,k,l) - 0.5_wp*mu_L*(dvel_L(j,k+1,l))
                                flux_vf(E_idx)%sf(j,k,l) = flux_vf(E_idx)%sf(j,k,l) - 0.5_wp*mu_L*qL_rs_vf(j,k+1,l,momxb+1)*dvel_L(j,k+1,l) 

                                dvel_L(j, k+1, l) = (1._wp/60._wp) * (-3._wp * dvel2(j, k-1, l) + &
                                                        27._wp * dvel2(j, k, l) + &
                                                        47._wp * dvel2(j, k+1, l) -   &
                                                        13._wp * dvel2(j, k+2, l) + &
                                                        2._wp * dvel2(j, k+3, l))

                                flux_vf(momxb+2)%sf(j,k,l) = flux_vf(momxb+2)%sf(j,k,l) - 0.5_wp*mu_L*(dvel_L(j,k+1,l))
                                flux_vf(E_idx)%sf(j,k,l) = flux_vf(E_idx)%sf(j,k,l) - 0.5_wp*mu_L*qL_rs_vf(j,k+1,l,momxb+1)*dvel_L(j,k+1,l)
                                

                                if(num_fluids > 1) then 
                                    mu_L = 0._wp
                                    !$acc loop seq
                                    do q = 1, Re_size(1)
                                        mu_L =  qR_rs_vf(j,k,l,advxb+Re_idx(1, q)-1) / Res(1, q) &
                                                  + mu_L
                                    end do
                                end if

                                dvel_L(j, k, l) = (1._wp/60._wp) * (-3._wp * dvel1(j, k+2, l) + &
                                                        27._wp * dvel1(j, k+1, l) + &
                                                        47._wp * dvel1(j, k, l) -   &
                                                        13._wp * dvel1(j, k-1, l) + &
                                                        2._wp * dvel1(j, k-2, l))
                                

                                flux_vf(momxb+2)%sf(j,k,l) = flux_vf(momxb+2)%sf(j,k,l) - 0.5_wp*mu_L*(dvel_L(j,k,l))
                                flux_vf(E_idx)%sf(j,k,l) = flux_vf(E_idx)%sf(j,k,l) - 0.5_wp*mu_L*qR_rs_vf(j,k,l,1)*dvel_L(j,k,l) 

                                dvel_L(j, k, l) = (1._wp/60._wp) * (-3._wp * dvel2(j, k+2, l) + &
                                                        27._wp * dvel2(j, k+1, l) + &
                                                        47._wp * dvel2(j, k, l) -   &
                                                        13._wp * dvel2(j, k-1, l) + &
                                                        2._wp * dvel2(j, k-2, l))

                                flux_vf(momxb+2)%sf(j,k,l) = flux_vf(momxb+2)%sf(j,k,l) - 0.5_wp*mu_L*(dvel_L(j,k,l))
                                flux_vf(E_idx)%sf(j,k,l) = flux_vf(E_idx)%sf(j,k,l) - 0.5_wp*mu_L*qR_rs_vf(j,k,l,1)*dvel_L(j,k,l)

                            end if


                            !! dvy & dux
                            dvel1(j,k,l) = (1/(12._wp*dy(k))) * ( &
                            8._wp*q_prim_vf(momxb+1)%sf(j,k+1,l) - &
                            8._wp*q_prim_vf(momxb+1)%sf(j,k-1,l) + &
                            q_prim_vf(momxb+1)%sf(j,k-2,l) - &
                            q_prim_vf(momxb+1)%sf(j,k+2,l) )

                            dvel2(j,k,l) = (1/(12._wp*dx(j))) * ( &
                            8._wp*q_prim_vf(momxb)%sf(j+1,k,l) - &
                            8._wp*q_prim_vf(momxb)%sf(j-1,k,l) + &
                            q_prim_vf(momxb)%sf(j-2,k,l) - &
                            q_prim_vf(momxb)%sf(j+2,k,l) )
                        end do 
                    end do 
                end do 

                !$acc parallel loop collapse(3) gang vector default(present) private(rho_L,gamma_L,pi_inf_L,mu_L,a_L)
                do l = -buff_size+2,p+buff_size-3
                    do k = -buff_size+2, n+buff_size-3
                        do j = -buff_size+2, m+buff_size-3

                            if(viscous) then
                                dvel_L(j, k+1, l) = (1._wp/60._wp) * (-3._wp * dvel1(j, k-1, l) + &
                                                        27._wp * dvel1(j, k, l) + &
                                                        47._wp * dvel1(j, k+1, l) -   &
                                                        13._wp * dvel1(j, k+2, l) + &
                                                        2._wp * dvel1(j, k+3, l))

                                mu_L = 0._wp
                                !$acc loop seq
                                do q = 1, Re_size(1)
                                    mu_L =  qL_rs_vf(j,k+1,l,advxb+Re_idx(1, q)-1) / Res(1, q) &
                                              + mu_L
                                end do
                                

                                flux_vf(momxb+1)%sf(j,k,l) = flux_vf(momxb+1)%sf(j,k,l) - 0.5_wp*mu_L*((4._wp/3._wp)*dvel_L(j,k+1,l))
                                flux_vf(E_idx)%sf(j,k,l) = flux_vf(E_idx)%sf(j,k,l) - 0.5_wp*mu_L*(4._wp/3._wp)*qL_rs_vf(j,k+1,l,momxb+1)*dvel_L(j,k+1,l) 

                                dvel_L(j, k+1, l) = (1._wp/60._wp) * (-3._wp * dvel2(j, k-1, l) + &
                                                        27._wp * dvel2(j, k, l) + &
                                                        47._wp * dvel2(j, k+1, l) -   &
                                                        13._wp * dvel2(j, k+2, l) + &
                                                        2._wp * dvel2(j, k+3, l))

                                flux_vf(momxb+1)%sf(j,k,l) = flux_vf(momxb+1)%sf(j,k,l) + 0.5_wp*mu_L*((2._wp/3._wp)*dvel_L(j,k+1,l))
                                flux_vf(E_idx)%sf(j,k,l) = flux_vf(E_idx)%sf(j,k,l) + 0.5_wp*mu_L*(2._wp/3._wp)*qL_rs_vf(j,k+1,l,momxb+1)*dvel_L(j,k+1,l) 

                                if(num_fluids > 1) then 
                                    mu_L = 0._wp
                                    !$acc loop seq
                                    do q = 1, Re_size(1)
                                        mu_L =  qR_rs_vf(j,k,l,advxb+Re_idx(1, q)-1) / Res(1, q) &
                                                  + mu_L
                                    end do
                                end if

                                dvel_L(j, k, l) = (1._wp/60._wp) * (-3._wp * dvel1(j, k+2, l) + &
                                                        27._wp * dvel1(j, k+1, l) + &
                                                        47._wp * dvel1(j, k, l) -   &
                                                        13._wp * dvel1(j, k-1, l) + &
                                                        2._wp * dvel1(j, k-2, l))
                                

                                flux_vf(momxb+1)%sf(j,k,l) = flux_vf(momxb+1)%sf(j,k,l) - 0.5_wp*mu_L*(4._wp/3._wp)*(dvel_L(j,k,l))
                                flux_vf(E_idx)%sf(j,k,l) = flux_vf(E_idx)%sf(j,k,l) - 0.5_wp*mu_L*(4._wp/3._wp)*qR_rs_vf(j,k,l,1)*dvel_L(j,k,l) 

                                dvel_L(j, k, l) = (1._wp/60._wp) * (-3._wp * dvel2(j, k+2, l) + &
                                                        27._wp * dvel2(j, k+1, l) + &
                                                        47._wp * dvel2(j, k, l) -   &
                                                        13._wp * dvel2(j, k-1, l) + &
                                                        2._wp * dvel2(j, k-2, l))

                                flux_vf(momxb+1)%sf(j,k,l) = flux_vf(momxb+1)%sf(j,k,l) + 0.5_wp*mu_L*((2._wp/3._wp)*dvel_L(j,k,l))
                                flux_vf(E_idx)%sf(j,k,l) = flux_vf(E_idx)%sf(j,k,l) + 0.5_wp*mu_L*(2._wp/3._wp)*qR_rs_vf(j,k,l,1)*dvel_L(j,k,l)

                            end if

                            !! dwz
                            dvel2(j,k,l) = (1/(12._wp*dz(l))) * ( &
                            8._wp*q_prim_vf(momxb+2)%sf(j,k,l+1) - &
                            8._wp*q_prim_vf(momxb+2)%sf(j,k,l-1) + &
                            q_prim_vf(momxb+2)%sf(j,k,l-2) - &
                            q_prim_vf(momxb+2)%sf(j,k,l+2) )
                        end do 
                    end do 
                end do 

                !$acc parallel loop collapse(3) gang vector default(present) private(rho_L,gamma_L,pi_inf_L,mu_L,a_L)
                do l = -buff_size+2,p+buff_size-3
                    do k = -buff_size+2, n+buff_size-3
                        do j = -buff_size+2, m+buff_size-3

                            if(viscous) then

                                mu_L = 0._wp
                                !$acc loop seq
                                do q = 1, Re_size(1)
                                    mu_L =  qL_rs_vf(j,k+1,l,advxb+Re_idx(1, q)-1) / Res(1, q) &
                                              + mu_L
                                end do

                                dvel_L(j, k+1, l) = (1._wp/60._wp) * (-3._wp * dvel2(j, k-1, l) + &
                                                        27._wp * dvel2(j, k, l) + &
                                                        47._wp * dvel2(j, k+1, l) -   &
                                                        13._wp * dvel2(j, k+2, l) + &
                                                        2._wp * dvel2(j, k+3, l))

                                flux_vf(momxb+1)%sf(j,k,l) = flux_vf(momxb+1)%sf(j,k,l) + 0.5_wp*mu_L*((2._wp/3._wp)*dvel_L(j,k+1,l))
                                flux_vf(E_idx)%sf(j,k,l) = flux_vf(E_idx)%sf(j,k,l) + 0.5_wp*mu_L*(2._wp/3._wp)*qL_rs_vf(j,k+1,l,momxb+1)*dvel_L(j,k+1,l) 

                                if(num_fluids > 1) then 
                                    mu_L = 0._wp
                                    !$acc loop seq
                                    do q = 1, Re_size(1)
                                        mu_L =  qR_rs_vf(j,k,l,advxb+Re_idx(1, q)-1) / Res(1, q) &
                                                  + mu_L
                                    end do
                                end if

                                dvel_L(j, k, l) = (1._wp/60._wp) * (-3._wp * dvel2(j, k+2, l) + &
                                                        27._wp * dvel2(j, k+1, l) + &
                                                        47._wp * dvel2(j, k, l) -   &
                                                        13._wp * dvel2(j, k-1, l) + &
                                                        2._wp * dvel2(j, k-2, l))

                                flux_vf(momxb+1)%sf(j,k,l) = flux_vf(momxb+1)%sf(j,k,l) + 0.5_wp*mu_L*((2._wp/3._wp)*dvel_L(j,k,l))
                                flux_vf(E_idx)%sf(j,k,l) = flux_vf(E_idx)%sf(j,k,l) + 0.5_wp*mu_L*(2._wp/3._wp)*qR_rs_vf(j,k,l,1)*dvel_L(j,k,l)

                            end if

                            !$acc loop seq 
                            do i = 1, momxb 
                                qL_rs_vf(j, k, l, i) = (1._wp/60._wp) * (-3._wp * q_prim_vf(i)%sf(j, k+2, l) + &
                                                27._wp * q_prim_vf(i)%sf(j, k+1, l) + &
                                                47._wp * q_prim_vf(i)%sf(j, k, l) -   &
                                                13._wp * q_prim_vf(i)%sf(j, k-1, l) + &
                                                2._wp * q_prim_vf(i)%sf(j, k-2, l))
                            end do

                            !$acc loop seq 
                            do i = momxb+2,advxe
                                qL_rs_vf(j, k, l, i) = (1._wp/60._wp) * (-3._wp * q_prim_vf(i)%sf(j, k+2, l) + &
                                                27._wp * q_prim_vf(i)%sf(j, k+1, l) + &
                                                47._wp * q_prim_vf(i)%sf(j, k, l) -   &
                                                13._wp * q_prim_vf(i)%sf(j, k-1, l) + &
                                                2._wp * q_prim_vf(i)%sf(j, k-2, l))
                            end do 

                            if(sys_size > advxe) then 
                                do i = advxe + 1, sys_size
                                    qL_rs_vf(j, k, l, i) = (1._wp/60._wp) * (-3._wp * q_prim_vf(i)%sf(j, k+2, l) + &
                                                    27._wp * q_prim_vf(i)%sf(j, k+1, l) + &
                                                    47._wp * q_prim_vf(i)%sf(j, k, l) -   &
                                                    13._wp * q_prim_vf(i)%sf(j, k-1, l) + &
                                                    2._wp * q_prim_vf(i)%sf(j, k-2, l))
                                end do
                            end if  

                            rho_L = sum(qL_rs_vf(j,k,l,1:num_fluids))
                            gamma_L = sum(qL_rs_vf(j,k,l,advxb:advxb+num_fluids-1)*gammas)
                            pi_inf_L = sum(qL_rs_vf(j,k,l,advxb:advxb+num_fluids-1)*pi_infs)

                            !$acc loop seq
                            do i = 1, num_fluids
                                flux_vf(i)%sf(j,k,l) = flux_vf(i)%sf(j,k,l) + &
                                    0.5_wp * (qL_rs_vf(j,k,l,i) * &
                                    qR_rs_vf(j,k,l,1)) + &
                                    0.5_wp*cfl(j,k,l) * (qL_rs_vf(j,k,l,i))
                            end do

                            if(num_fluids > 1) then 
                                !$acc loop seq
                                do i = 1, num_fluids
                                    flux_vf(advxb+i-1)%sf(j,k,l) = &
                                        0.5_wp * (qR_rs_vf(j,k,l,advxb+i-1) * &
                                        qR_rs_vf(j,k,l,1)) + &
                                        0.5_wp*cfl(j,k,l) * (qR_rs_vf(j,k,l,advxb+i-1))

                                    flux_vf(advxb+i-1)%sf(j,k,l) = flux_vf(advxb+i-1)%sf(j,k,l) &
                                    - 0.5_wp * q_prim_vf(advxb+i-1)%sf(j,k,l) * qR_rs_vf(j,k,l,1)
                                end do
                            end if

                            ! Momentum -> rho*u^2 + p + [[F_igr]]
                            flux_vf(momxb+1)%sf(j,k,l) = flux_vf(momxb+1)%sf(j,k,l) + &
                                 0.5_wp* (rho_L * (qR_rs_vf(j,k,l,1))**2.0 + &
                                 qL_rs_vf(j,k,l,E_idx)) + &
                                 0.5_wp*cfl(j,k,l) * (qR_rs_vf(j,k,l,1))

                            flux_vf(momxb)%sf(j, k, l) = flux_vf(momxb)%sf(j, k, l) +  &
                                0.5_wp * rho_L * qL_rs_vf(j,k,l,momxb)*qR_rs_vf(j,k,l,1) + &
                                0.5_wp*cfl(j,k,l) * (qL_rs_vf(j,k,l,momxb))

                            flux_vf(momxb+2)%sf(j, k, l) = flux_vf(momxb+2)%sf(j, k, l) +  &
                                0.5_wp * rho_L * qL_rs_vf(j,k,l,momxb+2)*qR_rs_vf(j,k,l,1) + &
                                0.5_wp*cfl(j,k,l) * (qL_rs_vf(j,k,l,momxb+2))

                             flux_vf(E_idx)%sf(j, k, l) = flux_vf(E_idx)%sf(j, k, l) + &
                                0.5_wp * (qR_rs_vf(j,k,l,1) * (qL_rs_vf(j,k,l,E_idx)*gamma_L + pi_inf_L + &
                                0.5_wp * rho_L * (qL_rs_vf(j,k,l,momxb)**2._wp + qR_rs_vf(j,k,l,1)**2._wp + qL_rs_vf(j,k,l,momxb+2)**2._wp ) + &
                                qL_rs_vf(j,k,l,E_idx)) ) + &
                                0.5_wp*cfl(j,k,l) * (qL_rs_vf(j,k,l,E_idx))

                        end do
                    end do
                end do
            end if
        elseif (idir == 3) then
           !$acc parallel loop collapse(3) gang vector default(present) private(rho_L,gamma_L,pi_inf_L,mu_L,a_L)
            do l = -buff_size+2, p+buff_size-3
                do k = -buff_size+2, n+buff_size-3
                    do j = -buff_size+2, m+buff_size-3

                        !$acc loop seq
                        do i = 1, sys_size
                            qL_rs_vf(j, k, l+1, i) = (1._wp/60._wp) * (-3._wp * q_prim_vf(i)%sf(j, k, l-1) + &
                                            27._wp * q_prim_vf(i)%sf(j, k, l) + &
                                            47._wp * q_prim_vf(i)%sf(j, k, l+1) -   &
                                            13._wp * q_prim_vf(i)%sf(j, k, l+2) + &
                                            2._wp * q_prim_vf(i)%sf(j, k, l+3))
                        end do

                        qR_rs_vf(j, k, l, 1) = (1._wp/60._wp) * (-3._wp * q_prim_vf(momxb+2)%sf(j, k, l+2) + &
                                            27._wp * q_prim_vf(momxb+2)%sf(j, k, l+1) + &
                                            47._wp * q_prim_vf(momxb+2)%sf(j, k, l) -   &
                                            13._wp * q_prim_vf(momxb+2)%sf(j, k, l-1) + &
                                            2._wp * q_prim_vf(momxb+2)%sf(j, k, l-2))
                        
                        if(num_fluids > 1) then  
                            !$acc loop seq 
                            do i = advxb, advxe                           
                                qR_rs_vf(j, k, l, 1) = (1._wp/60._wp) * (-3._wp * q_prim_vf(i)%sf(j, k, l+2) + &
                                            27._wp * q_prim_vf(i)%sf(j, k, l+1) + &
                                            47._wp * q_prim_vf(i)%sf(j, k, l) -   &
                                            13._wp * q_prim_vf(i)%sf(j, k, l-1) + &
                                            2._wp * q_prim_vf(i)%sf(j, k, l-2))
                            end do 
                        end if

                        rho_L = sum(qL_rs_vf(j,k,l+1,1:num_fluids))
                        gamma_L = sum(qL_rs_vf(j,k,l+1,advxb:advxb+num_fluids-1)*gammas)
                        pi_inf_L = sum(qL_rs_vf(j,k,l+1,advxb:advxb+num_fluids-1)*pi_infs)

                        a_L = sqrt((qL_rs_vf(j,k,l+1,E_idx)*(1._wp/gamma_L + 1._wp) + pi_inf_L / gamma_L) / rho_L)

                        cfl(j,k,l) = (sqrt(qL_rs_vf(j,k,l+1,momxb)**2._wp + qL_rs_vf(j,k,l+1,momxb+1)**2._wp + qL_rs_vf(j,k,l+1,momxb+2)**2._wp) + a_L)

                        !$acc loop seq
                        do i = 1, num_fluids
                            flux_vf(i)%sf(j,k,l) = &
                                0.5_wp * (qL_rs_vf(j,k,l+1,i) * &
                                qL_rs_vf(j,k,l+1,momxb+2)) - &
                                0.5_wp*cfl(j,k,l) * (qL_rs_vf(j,k,l+1,i))
                        end do

                        if(num_fluids > 1) then 
                            !$acc loop seq
                            do i = 1, num_fluids
                                flux_vf(advxb+i-1)%sf(j,k,l) = &
                                    0.5_wp * (qL_rs_vf(j,k,l+1,advxb+i-1) * &
                                    qL_rs_vf(j,k,l+1,momxb+2 )) - &
                                    0.5_wp*cfl(j,k,l) * (qL_rs_vf(j,k,l+1,advxb+i-1))

                                flux_vf(advxb+i-1)%sf(j,k,l) = flux_vf(advxb+i-1)%sf(j,k,l) &
                                - 0.5_wp * q_prim_vf(advxb+i-1)%sf(j,k,l) * qL_rs_vf(j,k,l+1,momxb+2)
                            end do
                        end if

                        ! Momentum -> rho*u^2 + p + [[F_igr]]
                        flux_vf(momxb+2)%sf(j,k,l) = &
                             0.5_wp* (rho_L * (qL_rs_vf(j,k,l+1,momxb+2))**2.0 + &
                             qL_rs_vf(j,k,l+1,E_idx)) - &
                             0.5_wp*cfl(j,k,l) * (qL_rs_vf(j,k,l+1,momxb+2))

                        flux_vf(momxb)%sf(j, k, l) =  &
                            0.5_wp * rho_L * qL_rs_vf(j,k,l+1,momxb)*qL_rs_vf(j,k,l+1,momxb+2) - &
                            0.5_wp*cfl(j,k,l) * (qL_rs_vf(j,k,l+1,momxb))

                        flux_vf(momxb+1)%sf(j, k, l) =  &
                            0.5_wp * rho_L * qL_rs_vf(j,k,l+1,momxb+2)*qL_rs_vf(j,k,l+1,momxb+1) - &
                            0.5_wp*cfl(j,k,l) * (qL_rs_vf(j,k,l+1,momxb+1))

                         flux_vf(E_idx)%sf(j, k, l) = &
                            0.5_wp * (qL_rs_vf(j,k,l+1,momxb+2) * (qL_rs_vf(j,k,l+1,E_idx)*gamma_L + pi_inf_L + &
                            0.5_wp * rho_L * (qL_rs_vf(j,k,l+1,momxb)**2._wp + qL_rs_vf(j,k,l+1,momxb+1)**2._wp +qL_rs_vf(j,k,l+1,momxb+2)**2._wp) + &
                            qL_rs_vf(j,k,l+1,E_idx)) ) - &
                            0.5_wp*cfl(j,k,l) * (qL_rs_vf(j,k,l+1,E_idx))

                        !! duz & dwx
                        dvel1(j,k,l) = (1/(12._wp*dz(l))) * ( &
                        8._wp*q_prim_vf(momxb)%sf(j,k,l+1) - &
                        8._wp*q_prim_vf(momxb)%sf(j,k,l-1) + &
                        q_prim_vf(momxb)%sf(j,k,l-2) - &
                        q_prim_vf(momxb)%sf(j,k,l+2) )

                        dvel2(j,k,l) = (1/(12._wp*dx(j))) * ( &
                        8._wp*q_prim_vf(momxb+2)%sf(j+1,k,l) - &
                        8._wp*q_prim_vf(momxb+2)%sf(j-1,k,l) + &
                        q_prim_vf(momxb+2)%sf(j-2,k,l) - &
                        q_prim_vf(momxb+2)%sf(j+2,k,l) )
                    end do 
                end do 
            end do 

            !$acc parallel loop collapse(3) gang vector default(present) private(rho_L,gamma_L,pi_inf_L,mu_L,a_L)
            do l = -buff_size+2,p+buff_size-3
                do k = -buff_size+2, n+buff_size-3
                    do j = -buff_size+2, m+buff_size-3

                        if(viscous) then
                            dvel_L(j, k, l+1) = (1._wp/60._wp) * (-3._wp * dvel1(j, k, l-1) + &
                                                    27._wp * dvel1(j, k, l) + &
                                                    47._wp * dvel1(j, k, l+1) -   &
                                                    13._wp * dvel1(j, k, l+2) + &
                                                    2._wp * dvel1(j, k, l+3))

                            mu_L = 0._wp
                            !$acc loop seq
                            do q = 1, Re_size(1)
                                mu_L =  qL_rs_vf(j,k,l+1,advxb+Re_idx(1, q)-1) / Res(1, q) &
                                          + mu_L
                            end do
                            

                            flux_vf(momxb)%sf(j,k,l) = flux_vf(momxb)%sf(j,k,l) - 0.5_wp*mu_L*(dvel_L(j,k,l+1))
                            flux_vf(E_idx)%sf(j,k,l) = flux_vf(E_idx)%sf(j,k,l) - 0.5_wp*mu_L*qL_rs_vf(j,k,l+1,momxb+2)*dvel_L(j,k,l+1) 

                            dvel_L(j, k, l+1) = (1._wp/60._wp) * (-3._wp * dvel2(j, k, l-1) + &
                                                    27._wp * dvel2(j, k, l) + &
                                                    47._wp * dvel2(j, k, l+1) -   &
                                                    13._wp * dvel2(j, k, l+2) + &
                                                    2._wp * dvel2(j, k, l+3))

                            flux_vf(momxb)%sf(j,k,l) = flux_vf(momxb)%sf(j,k,l) - 0.5_wp*mu_L*(dvel_L(j,k,l+1))
                            flux_vf(E_idx)%sf(j,k,l) = flux_vf(E_idx)%sf(j,k,l) - 0.5_wp*mu_L*qL_rs_vf(j,k,l+1,momxb+2)*dvel_L(j,k,l+1)

                            if(num_fluids > 1) then 
                                mu_L = 0._wp
                                !$acc loop seq
                                do q = 1, Re_size(1)
                                    mu_L =  qR_rs_vf(j,k,l,advxb+Re_idx(1, q)-1) / Res(1, q) &
                                              + mu_L
                                end do
                            end if

                            dvel_L(j, k, l) = (1._wp/60._wp) * (-3._wp * dvel1(j, k, l+2) + &
                                                    27._wp * dvel1(j, k, l+1) + &
                                                    47._wp * dvel1(j, k, l) -   &
                                                    13._wp * dvel1(j, k, l-1) + &
                                                    2._wp * dvel1(j, k, l-2))
                            

                            flux_vf(momxb)%sf(j,k,l) = flux_vf(momxb)%sf(j,k,l) - 0.5_wp*mu_L*(dvel_L(j,k,l))
                            flux_vf(E_idx)%sf(j,k,l) = flux_vf(E_idx)%sf(j,k,l) - 0.5_wp*mu_L*qR_rs_vf(j,k,l,1)*dvel_L(j,k,l) 

                            dvel_L(j, k, l) = (1._wp/60._wp) * (-3._wp * dvel2(j, k, l+2) + &
                                                    27._wp * dvel2(j, k, l+1) + &
                                                    47._wp * dvel2(j, k, l) -   &
                                                    13._wp * dvel2(j, k, l-1) + &
                                                    2._wp * dvel2(j, k, l-2))

                            flux_vf(momxb)%sf(j,k,l) = flux_vf(momxb)%sf(j,k,l) - 0.5_wp*mu_L*(dvel_L(j,k,l))
                            flux_vf(E_idx)%sf(j,k,l) = flux_vf(E_idx)%sf(j,k,l) - 0.5_wp*mu_L*qR_rs_vf(j,k,l,1)*dvel_L(j,k,l)

                        end if

                        !! dvz & dwy
                        dvel1(j,k,l) = (1/(12._wp*dz(l))) * ( &
                        8._wp*q_prim_vf(momxb+1)%sf(j,k,l+1) - &
                        8._wp*q_prim_vf(momxb+1)%sf(j,k,l-1) + &
                        q_prim_vf(momxb+1)%sf(j,k,l-2) - &
                        q_prim_vf(momxb+1)%sf(j,k,l+2) )

                        dvel2(j,k,l) = (1/(12._wp*dy(k))) * ( &
                        8._wp*q_prim_vf(momxb+2)%sf(j,k+1,l) - &
                        8._wp*q_prim_vf(momxb+2)%sf(j,k-1,l) + &
                        q_prim_vf(momxb+2)%sf(j,k-2,l) - &
                        q_prim_vf(momxb+2)%sf(j,k+2,l) )
                    end do 
                end do 
            end do 

            !$acc parallel loop collapse(3) gang vector default(present) private(rho_L,gamma_L,pi_inf_L,mu_L,a_L)
            do l = -buff_size+2,p+buff_size-3
                do k = -buff_size+2, n+buff_size-3
                    do j = -buff_size+2, m+buff_size-3

                        if(viscous) then
                            dvel_L(j, k, l+1) = (1._wp/60._wp) * (-3._wp * dvel1(j, k, l-1) + &
                                                    27._wp * dvel1(j, k, l) + &
                                                    47._wp * dvel1(j, k, l+1) -   &
                                                    13._wp * dvel1(j, k, l+2) + &
                                                    2._wp * dvel1(j, k, l+3))

                            mu_L = 0._wp
                            !$acc loop seq
                            do q = 1, Re_size(1)
                                mu_L =  qL_rs_vf(j,k,l+1,advxb+Re_idx(1, q)-1) / Res(1, q) &
                                          + mu_L
                            end do
                            

                            flux_vf(momxb+1)%sf(j,k,l) = flux_vf(momxb+1)%sf(j,k,l) - 0.5_wp*mu_L*(dvel_L(j,k,l+1))
                            flux_vf(E_idx)%sf(j,k,l) = flux_vf(E_idx)%sf(j,k,l) - 0.5_wp*mu_L*qL_rs_vf(j,k,l+1,momxb+2)*dvel_L(j,k,l+1) 

                            dvel_L(j, k, l+1) = (1._wp/60._wp) * (-3._wp * dvel2(j, k, l-1) + &
                                                    27._wp * dvel2(j, k, l) + &
                                                    47._wp * dvel2(j, k, l+1) -   &
                                                    13._wp * dvel2(j, k, l+2) + &
                                                    2._wp * dvel2(j, k, l+3))

                            flux_vf(momxb+1)%sf(j,k,l) = flux_vf(momxb+1)%sf(j,k,l) - 0.5_wp*mu_L*(dvel_L(j,k,l+1))
                            flux_vf(E_idx)%sf(j,k,l) = flux_vf(E_idx)%sf(j,k,l) - 0.5_wp*mu_L*qL_rs_vf(j,k,l+1,momxb+2)*dvel_L(j,k,l+1)

                            if(num_fluids > 1) then 
                                mu_L = 0._wp
                                !$acc loop seq
                                do q = 1, Re_size(1)
                                    mu_L =  qR_rs_vf(j,k,l,advxb+Re_idx(1, q)-1) / Res(1, q) &
                                              + mu_L
                                end do
                            end if

                            dvel_L(j, k, l) = (1._wp/60._wp) * (-3._wp * dvel1(j, k, l+2) + &
                                                    27._wp * dvel1(j, k, l+1) + &
                                                    47._wp * dvel1(j, k, l) -   &
                                                    13._wp * dvel1(j, k, l-1) + &
                                                    2._wp * dvel1(j, k, l-2))
                            

                            flux_vf(momxb+1)%sf(j,k,l) = flux_vf(momxb+1)%sf(j,k,l) - 0.5_wp*mu_L*(dvel_L(j,k,l))
                            flux_vf(E_idx)%sf(j,k,l) = flux_vf(E_idx)%sf(j,k,l) - 0.5_wp*mu_L*qR_rs_vf(j,k,l,1)*dvel_L(j,k,l) 

                            dvel_L(j, k, l) = (1._wp/60._wp) * (-3._wp * dvel2(j, k, l+2) + &
                                                    27._wp * dvel2(j, k, l+1) + &
                                                    47._wp * dvel2(j, k, l) -   &
                                                    13._wp * dvel2(j, k, l-1) + &
                                                    2._wp * dvel2(j, k, l-2))

                            flux_vf(momxb+1)%sf(j,k,l) = flux_vf(momxb+1)%sf(j,k,l) - 0.5_wp*mu_L*(dvel_L(j,k,l))
                            flux_vf(E_idx)%sf(j,k,l) = flux_vf(E_idx)%sf(j,k,l) - 0.5_wp*mu_L*qR_rs_vf(j,k,l,1)*dvel_L(j,k,l)

                        end if


                        !! dwz & dux
                        dvel1(j,k,l) = (1/(12._wp*dz(l))) * ( &
                        8._wp*q_prim_vf(momxb+2)%sf(j,k,l+1) - &
                        8._wp*q_prim_vf(momxb+2)%sf(j,k,l-1) + &
                        q_prim_vf(momxb+2)%sf(j,k,l-2) - &
                        q_prim_vf(momxb+2)%sf(j,k,l+2) )

                        dvel2(j,k,l) = (1/(12._wp*dx(j))) * ( &
                        8._wp*q_prim_vf(momxb)%sf(j+1,k,l) - &
                        8._wp*q_prim_vf(momxb)%sf(j-1,k,l) + &
                        q_prim_vf(momxb)%sf(j-2,k,l) - &
                        q_prim_vf(momxb)%sf(j+2,k,l) )
                    end do 
                end do 
            end do 

            !$acc parallel loop collapse(3) gang vector default(present) private(rho_L,gamma_L,pi_inf_L,mu_L,a_L)
            do l = -buff_size+2,p+buff_size-3
                do k = -buff_size+2, n+buff_size-3
                    do j = -buff_size+2, m+buff_size-3

                        if(viscous) then
                            dvel_L(j, k, l+1) = (1._wp/60._wp) * (-3._wp * dvel1(j, k, l-1) + &
                                                    27._wp * dvel1(j, k, l) + &
                                                    47._wp * dvel1(j, k, l+1) -   &
                                                    13._wp * dvel1(j, k, l+2) + &
                                                    2._wp * dvel1(j, k, l+3))
 
                            mu_L = 0._wp
                            !$acc loop seq
                            do q = 1, Re_size(1)
                                mu_L =  qL_rs_vf(j,k,l+1,advxb+Re_idx(1, q)-1) / Res(1, q) &
                                          + mu_L
                            end do
                            

                            flux_vf(momxb+2)%sf(j,k,l) = flux_vf(momxb+2)%sf(j,k,l) - 0.5_wp*mu_L*((4._wp/3._wp)*dvel_L(j,k,l+1))
                            flux_vf(E_idx)%sf(j,k,l) = flux_vf(E_idx)%sf(j,k,l) - 0.5_wp*mu_L*(4._wp/3._wp)*qL_rs_vf(j,k,l+1,momxb+2)*dvel_L(j,k,l+1) 

                            dvel_L(j, k, l+1) = (1._wp/60._wp) * (-3._wp * dvel2(j, k, l-1) + &
                                                    27._wp * dvel2(j, k, l) + &
                                                    47._wp * dvel2(j, k, l+1) -   &
                                                    13._wp * dvel2(j, k, l+2) + &
                                                    2._wp * dvel2(j, k, l+3))

                            flux_vf(momxb+2)%sf(j,k,l) = flux_vf(momxb+2)%sf(j,k,l) + 0.5_wp*mu_L*((2._wp/3._wp)*dvel_L(j,k,l+1))
                            flux_vf(E_idx)%sf(j,k,l) = flux_vf(E_idx)%sf(j,k,l) + 0.5_wp*mu_L*(2._wp/3._wp)*qL_rs_vf(j,k,l+1,momxb+2)*dvel_L(j,k,l+1) 

                            if(num_fluids > 1) then 
                                mu_L = 0._wp
                                !$acc loop seq
                                do q = 1, Re_size(1)
                                    mu_L =  qR_rs_vf(j,k,l,advxb+Re_idx(1, q)-1) / Res(1, q) &
                                              + mu_L
                                end do
                            end if

                            dvel_L(j, k, l) = (1._wp/60._wp) * (-3._wp * dvel1(j, k, l+2) + &
                                                    27._wp * dvel1(j, k, l+1) + &
                                                    47._wp * dvel1(j, k, l) -   &
                                                    13._wp * dvel1(j, k, l-1) + &
                                                    2._wp * dvel1(j, k, l-2))
                            

                            flux_vf(momxb+2)%sf(j,k,l) = flux_vf(momxb+2)%sf(j,k,l) - 0.5_wp*mu_L*(4._wp/3._wp)*(dvel_L(j,k,l))
                            flux_vf(E_idx)%sf(j,k,l) = flux_vf(E_idx)%sf(j,k,l) - 0.5_wp*mu_L*(4._wp/3._wp)*qR_rs_vf(j,k,l,1)*dvel_L(j,k,l) 

                            dvel_L(j, k, l) = (1._wp/60._wp) * (-3._wp * dvel2(j, k, l+2) + &
                                                    27._wp * dvel2(j, k, l+1) + &
                                                    47._wp * dvel2(j, k, l) -   &
                                                    13._wp * dvel2(j, k, l-1) + &
                                                    2._wp * dvel2(j, k, l-2))

                            flux_vf(momxb+2)%sf(j,k,l) = flux_vf(momxb+2)%sf(j,k,l) + 0.5_wp*mu_L*((2._wp/3._wp)*dvel_L(j,k,l))
                            flux_vf(E_idx)%sf(j,k,l) = flux_vf(E_idx)%sf(j,k,l) + 0.5_wp*mu_L*(2._wp/3._wp)*qR_rs_vf(j,k,l,1)*dvel_L(j,k,l)

                        end if

                        !! dvy
                        dvel2(j,k,l) = (1/(12._wp*dy(k))) * ( &
                        8._wp*q_prim_vf(momxb+1)%sf(j,k+1,l) - &
                        8._wp*q_prim_vf(momxb+1)%sf(j,k-1,l) + &
                        q_prim_vf(momxb+1)%sf(j,k-2,l) - &
                        q_prim_vf(momxb+1)%sf(j,k+2,l) )
                    end do 
                end do 
            end do 

            !$acc parallel loop collapse(3) gang vector default(present) private(rho_L,gamma_L,pi_inf_L,mu_L,a_L)
            do l = -buff_size+2,p+buff_size-3
                do k = -buff_size+2, n+buff_size-3
                    do j = -buff_size+2, m+buff_size-3

                        if(viscous) then

                            mu_L = 0._wp
                            !$acc loop seq
                            do q = 1, Re_size(1)
                                mu_L =  qL_rs_vf(j,k,l+1,advxb+Re_idx(1, q)-1) / Res(1, q) &
                                          + mu_L
                            end do

                            dvel_L(j, k, l+1) = (1._wp/60._wp) * (-3._wp * dvel2(j, k, l-1) + &
                                                    27._wp * dvel2(j, k, l) + &
                                                    47._wp * dvel2(j, k, l+1) -   &
                                                    13._wp * dvel2(j, k, l+2) + &
                                                    2._wp * dvel2(j, k, l+3))

                            flux_vf(momxb+2)%sf(j,k,l) = flux_vf(momxb+2)%sf(j,k,l) + 0.5_wp*mu_L*((2._wp/3._wp)*dvel_L(j,k,l+1))
                            flux_vf(E_idx)%sf(j,k,l) = flux_vf(E_idx)%sf(j,k,l) + 0.5_wp*mu_L*(2._wp/3._wp)*qL_rs_vf(j,k,l+1,momxb+2)*dvel_L(j,k,l+1) 

                            if(num_fluids > 1) then 
                                mu_L = 0._wp
                                !$acc loop seq
                                do q = 1, Re_size(1)
                                    mu_L =  qR_rs_vf(j,k,l,advxb+Re_idx(1, q)-1) / Res(1, q) &
                                              + mu_L
                                end do
                            end if

                            dvel_L(j, k, l) = (1._wp/60._wp) * (-3._wp * dvel2(j, k, l+2) + &
                                                    27._wp * dvel2(j, k, l+1) + &
                                                    47._wp * dvel2(j, k, l) -   &
                                                    13._wp * dvel2(j, k, l-1) + &
                                                    2._wp * dvel2(j, k, l-2))

                            flux_vf(momxb+2)%sf(j,k,l) = flux_vf(momxb+2)%sf(j,k,l) + 0.5_wp*mu_L*((2._wp/3._wp)*dvel_L(j,k,l))
                            flux_vf(E_idx)%sf(j,k,l) = flux_vf(E_idx)%sf(j,k,l) + 0.5_wp*mu_L*(2._wp/3._wp)*qR_rs_vf(j,k,l,1)*dvel_L(j,k,l)

                        end if

                        !$acc loop seq 
                        do i = 1, momxb+1 
                            qL_rs_vf(j, k, l, i) = (1._wp/60._wp) * (-3._wp * q_prim_vf(i)%sf(j, k, l+2) + &
                                            27._wp * q_prim_vf(i)%sf(j, k, l+1) + &
                                            47._wp * q_prim_vf(i)%sf(j, k, l) -   &
                                            13._wp * q_prim_vf(i)%sf(j, k, l-1) + &
                                            2._wp * q_prim_vf(i)%sf(j, k, l-2))
                        end do

                        !$acc loop seq 
                        do i = momxb+3,advxe
                            qL_rs_vf(j, k, l, i) = (1._wp/60._wp) * (-3._wp * q_prim_vf(i)%sf(j, k, l+2) + &
                                            27._wp * q_prim_vf(i)%sf(j, k, l+1) + &
                                            47._wp * q_prim_vf(i)%sf(j, k, l) -   &
                                            13._wp * q_prim_vf(i)%sf(j, k, l-1) + &
                                            2._wp * q_prim_vf(i)%sf(j, k, l-2))
                        end do 

                        if(sys_size > advxe) then 
                            do i = advxe + 1, sys_size
                                qL_rs_vf(j, k, l, i) = (1._wp/60._wp) * (-3._wp * q_prim_vf(i)%sf(j, k, l+2) + &
                                            27._wp * q_prim_vf(i)%sf(j, k, l+1) + &
                                            47._wp * q_prim_vf(i)%sf(j, k, l) -   &
                                            13._wp * q_prim_vf(i)%sf(j, k, l-1) + &
                                            2._wp * q_prim_vf(i)%sf(j, k, l-2))
                            end do
                        end if  

                        rho_L = sum(qL_rs_vf(j,k,l,1:num_fluids))
                        gamma_L = sum(qL_rs_vf(j,k,l,advxb:advxb+num_fluids-1)*gammas)
                        pi_inf_L = sum(qL_rs_vf(j,k,l,advxb:advxb+num_fluids-1)*pi_infs)

                        !$acc loop seq
                        do i = 1, num_fluids
                            flux_vf(i)%sf(j,k,l) = flux_vf(i)%sf(j,k,l) + &
                                0.5_wp * (qL_rs_vf(j,k,l,i) * &
                                qR_rs_vf(j,k,l,1)) + &
                                0.5_wp*cfl(j,k,l) * (qL_rs_vf(j,k,l,i))
                        end do

                        if(num_fluids > 1) then 
                            !$acc loop seq
                            do i = 1, num_fluids
                                flux_vf(advxb+i-1)%sf(j,k,l) = &
                                    0.5_wp * (qR_rs_vf(j,k,l,advxb+i-1) * &
                                    qR_rs_vf(j,k,l,1)) + &
                                    0.5_wp*cfl(j,k,l) * (qR_rs_vf(j,k,l,advxb+i-1))

                                flux_vf(advxb+i-1)%sf(j,k,l) = flux_vf(advxb+i-1)%sf(j,k,l) &
                                - 0.5_wp * q_prim_vf(advxb+i-1)%sf(j,k,l) * qR_rs_vf(j,k,l,1)
                            end do
                        end if

                        ! Momentum -> rho*u^2 + p + [[F_igr]]
                        flux_vf(momxb+2)%sf(j,k,l) = flux_vf(momxb+2)%sf(j,k,l) + &
                             0.5_wp* (rho_L * (qR_rs_vf(j,k,l,1))**2.0 + &
                             qL_rs_vf(j,k,l,E_idx)) + &
                             0.5_wp*cfl(j,k,l) * (qR_rs_vf(j,k,l,1))

                        flux_vf(momxb)%sf(j, k, l) = flux_vf(momxb)%sf(j, k, l) +  &
                            0.5_wp * rho_L * qL_rs_vf(j,k,l,momxb)*qR_rs_vf(j,k,l,1) + &
                            0.5_wp*cfl(j,k,l) * (qL_rs_vf(j,k,l,momxb))

                        flux_vf(momxb+1)%sf(j, k, l) = flux_vf(momxb+1)%sf(j, k, l) +  &
                            0.5_wp * rho_L * qL_rs_vf(j,k,l,momxb+1)*qR_rs_vf(j,k,l,1) + &
                            0.5_wp*cfl(j,k,l) * (qL_rs_vf(j,k,l,momxb+1))

                         flux_vf(E_idx)%sf(j, k, l) = flux_vf(E_idx)%sf(j, k, l) + &
                            0.5_wp * (qR_rs_vf(j,k,l,1) * (qL_rs_vf(j,k,l,E_idx)*gamma_L + pi_inf_L + &
                            0.5_wp * rho_L * (qL_rs_vf(j,k,l,momxb)**2._wp + qL_rs_vf(j,k,l,momxb+1)**2._wp + qR_rs_vf(j,k,l,1)**2._wp) + &
                            qL_rs_vf(j,k,l,E_idx)) ) + &
                            0.5_wp*cfl(j,k,l) * (qL_rs_vf(j,k,l,E_idx))

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
