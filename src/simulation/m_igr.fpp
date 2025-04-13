#:include 'macros.fpp'

module m_igr

    use m_derived_types        !< Definitions of the derived types

    use m_global_parameters

    use m_variables_conversion

    use m_mpi_proxy

    use m_helper

    use atomic_mod
    
    implicit none

    private; public :: s_initialize_igr_module, &
        s_igr_jacobi_iteration, &
        s_igr_riemann_solver, &
        s_igr_sigma, &
        s_initialize_igr, &
        s_igr_flux_add

#ifdef __NVCOMPILER_GPU_UNIFIED_MEM
    real(2), pointer, contiguous, dimension(:, :, :) :: jac,jac_rhs
#else
    real(2), allocatable, dimension(:, :, :) :: jac,jac_rhs
    !$acc declare create(jac, jac_rhs)
#endif

    real(wp) :: alf_igr, omega, mu, bcxb, bcxe, bcyb, bcye, bczb, bcze
    !$acc declare create(alf_igr, omega, mu, bcxb, bcxe, bcyb, bcye, bczb, bcze)

    type(int_bounds_info) :: ix, iy, iz
    !$acc declare create(ix, iy, iz)

    real(wp), allocatable, dimension(:, :) :: Res
    !$acc declare create(Res)

#ifdef _CRAYFTN
    real(wp), allocatable, dimension(:) :: coeff_L, coeff_R
    !$acc declare create(coeff_L, coeff_R)
#else
    real(wp), parameter :: coeff_L(-1:3) = [ &
        -3._wp/60._wp,   &  ! Index -1
        27._wp/60._wp,   &  ! Index 0
        47._wp/60._wp,   &  ! Index 1
        -13._wp/60._wp,  &  ! Index 2
        2._wp/60._wp     &  ! Index 3
        ]

    real(wp), parameter :: coeff_R(-2:2) = [ &
        2._wp/60._wp,    &  ! Index -2
        -13._wp/60._wp,  &  ! Index -1
        47._wp/60._wp,   &  ! Index 0
        27._wp/60._wp,   &  ! Index 1
        -3._wp/60._wp    &  ! Index 2
        ]
#endif

    integer :: i, j, k, l, q, r

#ifdef __NVCOMPILER_GPU_UNIFIED_MEM
    real(wp), allocatable, dimension(:, :, :, :), pinned, target :: m_igr_pool_host
#endif

contains

    subroutine s_initialize_igr_module()

        integer :: igr_temps_on_gpu = 3
        integer :: igr_temps_on_cpu = 0
        integer :: pool_idx = 1
        character(len=10) :: igr_temps_on_gpu_str

#ifdef __NVCOMPILER_GPU_UNIFIED_MEM
        call get_environment_variable("NVIDIA_IGR_TEMPS_ON_GPU", igr_temps_on_gpu_str)

        if (trim(igr_temps_on_gpu_str) == "0") then
            igr_temps_on_gpu = 0 ! jac, jac_rhs and jac_old on CPU
        elseif (trim(igr_temps_on_gpu_str) == "1") then
            igr_temps_on_gpu = 1 ! jac on GPU, jac_rhs on CPU, jac_old on CPU
        elseif (trim(igr_temps_on_gpu_str) == "2") then
            igr_temps_on_gpu = 2 ! jac and jac_rhs on GPU, jac_old on CPU
        elseif (trim(igr_temps_on_gpu_str) == "3") then
            igr_temps_on_gpu = 3 ! jac, jac_rhs and jac_old on GPU
        else ! default on GPU
            igr_temps_on_gpu = 3
        endif
#endif

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
            @:PREFER_GPU(Res)
            @:PREFER_GPU(Re_idx)
        end if

        if(igr) then 
#ifdef __NVCOMPILER_GPU_UNIFIED_MEM
           igr_temps_on_cpu = 3 - igr_temps_on_gpu
           if ( igr_temps_on_cpu >= 1 ) then
               allocate(m_igr_pool_host(idwbuff(1)%beg:idwbuff(1)%end, &
                                        idwbuff(2)%beg:idwbuff(2)%end, &
                                        idwbuff(3)%beg:idwbuff(3)%end, &
                                        1:igr_temps_on_cpu))
               pool_idx = 1
               if ( igr_temps_on_cpu >= 1 ) then
                   !print*, 'jac_old on CPU'
                   !jac_old(idwbuff(1)%beg:idwbuff(1)%end, &
                   !    idwbuff(2)%beg:idwbuff(2)%end, &
                   !    idwbuff(3)%beg:idwbuff(3)%end) => m_igr_pool_host(:,:,:,pool_idx)
                   !pool_idx = pool_idx + 1
               end if
               if ( igr_temps_on_cpu >= 2 ) then
                   !print*, 'jac_rhs on CPU'
                   jac_rhs(idwbuff(1)%beg:idwbuff(1)%end, &
                       idwbuff(2)%beg:idwbuff(2)%end, &
                       idwbuff(3)%beg:idwbuff(3)%end) => m_igr_pool_host(:,:,:,pool_idx)
                   pool_idx = pool_idx + 1
               end if
               if ( igr_temps_on_cpu >= 3 ) then
                   !print*, 'jac on CPU'
                   jac(idwbuff(1)%beg:idwbuff(1)%end, &
                       idwbuff(2)%beg:idwbuff(2)%end, &
                       idwbuff(3)%beg:idwbuff(3)%end) => m_igr_pool_host(:,:,:,pool_idx)
                   pool_idx = pool_idx + 1
               end if
           end if
           if ( igr_temps_on_gpu >= 1 ) then
                !print*, 'jac on GPU'
                @:ALLOCATE(jac(idwbuff(1)%beg:idwbuff(1)%end, &
                             idwbuff(2)%beg:idwbuff(2)%end, &
                             idwbuff(3)%beg:idwbuff(3)%end))
                @:PREFER_GPU(jac)
           endif
           if ( igr_temps_on_gpu >= 2 ) then
                !print*, 'jac_rhs on GPU'
                @:ALLOCATE(jac_rhs(idwbuff(1)%beg:idwbuff(1)%end, &
                             idwbuff(2)%beg:idwbuff(2)%end, &
                             idwbuff(3)%beg:idwbuff(3)%end))
                @:PREFER_GPU(jac_rhs)
           endif
           if ( igr_temps_on_gpu >= 3 ) then
                !print*, 'jac_old on GPU'
                !@:ALLOCATE(jac_old(idwbuff(1)%beg:idwbuff(1)%end, &
                !             idwbuff(2)%beg:idwbuff(2)%end, &
                !             idwbuff(3)%beg:idwbuff(3)%end))
                !@:PREFER_GPU(jac_old)
           endif
#else
            #:for VAR in [ 'jac','jac_rhs']
                @:ALLOCATE(${VAR}$(idwbuff(1)%beg:idwbuff(1)%end, &
                             idwbuff(2)%beg:idwbuff(2)%end, &
                             idwbuff(3)%beg:idwbuff(3)%end))
            #:endfor
            if ( igr_temps_on_gpu >= 1 ) then
                @:PREFER_GPU(jac)
            end if
            if ( igr_temps_on_gpu >= 2 ) then
                @:PREFER_GPU(jac_rhs)
            end if
            if ( igr_temps_on_gpu >= 3) then 
                !@:PREFER_GPU(jac_old)
            end if 
#endif

            !$acc parallel loop collapse(3) gang vector default(present)
            do l = idwbuff(3)%beg, idwbuff(3)%end
                do k = idwbuff(2)%beg, idwbuff(2)%end
                    do j = idwbuff(1)%beg, idwbuff(1)%end
                        jac(j, k, l) = 0._wp
                   end do
                end do
            end do
#ifdef _CRAYFTN
            @:ALLOCATE(coeff_L(-1:3))
            coeff_L(-1) = (-3._wp/60._wp)
            coeff_L(0) = (27._wp/60._wp)
            coeff_L(1) = (47._wp/60._wp)
            coeff_L(2) = (-13._wp/60._wp)
            coeff_L(3) = (2._wp/60._wp)
            !$acc update device(coeff_L)

            @:ALLOCATE(coeff_R(-2:2))
            coeff_R(2) = (-3._wp/60._wp)
            coeff_R(1) = (27._wp/60._wp)
            coeff_R(0) = (47._wp/60._wp)
            coeff_R(-1) = (-13._wp/60._wp)
            coeff_R(-2) = (2._wp/60._wp)
            !$acc update device(coeff_R)
#endif
        end if

    end subroutine s_initialize_igr_module

    subroutine s_igr_jacobi_iteration(q_prim_vf, t_step)

        type(scalar_field_half), &
            dimension(sys_size), &
            intent(inout) :: q_prim_vf
        real(wp) :: rho_rx, rho_ry, rho_rz, rho_lx, rho_ly, rho_lz
        real(wp) :: resid, fd_coeff
        integer :: num_iters, t_step


        if (t_step  == t_step_start) then
            num_iters = 100
        else
            num_iters = num_igr_iters
        end if

        do q = 1, num_iters
            !$acc parallel loop collapse(3) gang vector default(present) private(rho_lx, rho_rx, rho_ly, rho_ry, rho_lz, rho_rz, fd_coeff)
            do l = 0, p
                do k = 0, n
                    do j = 0, m
                        rho_lx = 0._wp 
                        rho_rx = 0._wp 
                        rho_ly = 0._wp 
                        rho_ry = 0._wp
                        rho_lz = 0._wp 
                        rho_rz = 0._wp 
                        fd_coeff = 0._wp

                        !$acc loop seq
                        do i = 1, num_fluids
                            rho_lx = rho_lx + (q_prim_vf(i)%sf(j,k,l) + q_prim_vf(i)%sf(j-1,k,l)) / 2._wp
                            rho_rx = rho_rx + (q_prim_vf(i)%sf(j,k,l) + q_prim_vf(i)%sf(j+1,k,l)) / 2._wp
                            rho_ly = rho_ly + (q_prim_vf(i)%sf(j,k,l) + q_prim_vf(i)%sf(j,k-1,l)) / 2._wp
                            rho_ry = rho_ry + (q_prim_vf(i)%sf(j,k,l) + q_prim_vf(i)%sf(j,k+1,l)) / 2._wp
                            if(p > 0) then 
                                rho_lz = rho_lz + (q_prim_vf(i)%sf(j,k,l) + q_prim_vf(i)%sf(j,k,l-1)) / 2._wp
                                rho_rz = rho_rz + (q_prim_vf(i)%sf(j,k,l) + q_prim_vf(i)%sf(j,k,l+1)) / 2._wp
                            end if
                            fd_coeff = fd_coeff + q_prim_vf(i)%sf(j,k,l)
                        end do

                        fd_coeff = 1._wp / fd_coeff + alf_igr * ( (1._wp / dx(j)**2._wp) * (1._wp/rho_lx + 1._wp/rho_rx) + (1._wp / dy(k)**2._wp) *(1._wp/rho_ly + 1._wp/rho_ry) )

                        if(p > 0) then 
                            fd_coeff = fd_coeff + alf_igr* (1._wp / dz(l)**2._wp) * (1._wp/rho_lz + 1._wp /rho_rz)
                        end if

                        if(p > 0) then
                            jac(j, k, l) = (alf_igr / fd_coeff) * ( (1._wp / dx(j)**2._wp) * ( jac(j-1,k,l)/rho_lx + jac(j+1,k,l)/rho_rx) + &
                                                    (1._wp / dy(k)**2._wp) * (jac(j,k-1,l)/rho_ly + jac(j,k+1,l)/rho_ry) + &
                                                    (1._wp / dz(l)**2._wp) * (jac(j,k,l-1)/rho_lz + jac(j,k,l+1)/rho_rz) ) + &
                                                    jac_rhs(j,k,l) / fd_coeff
                       else 
                            jac(j, k, l) = (alf_igr / fd_coeff) * ( (1._wp / dx(j)**2._wp) * ( jac(j-1,k,l)/rho_lx + jac(j+1,k,l)/rho_rx) + &
                                                    (1._wp / dy(k)**2._wp) * (jac(j,k-1,l)/rho_ly + jac(j,k+1,l)/rho_ry)) + &
                                                    jac_rhs(j,k,l) / fd_coeff                             
                       end if
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
        end do
    end subroutine s_igr_jacobi_iteration

    subroutine s_igr_sigma(q_prim_vf, rhs_vf, idir)

        type(scalar_field_half), &
            dimension(sys_size), &
            intent(inout) :: rhs_vf
        type(scalar_field_half), &
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

                            @:ATOMIC_ADD(rhs_vf(momxb)%sf(j+1,k,l), (0.5_wp/dx(j+1) * F_L * dt))
                            ! !$omp atomic update
                            ! !$acc atomic
                            ! rhs_vf(momxb)%sf(j+1,k,l) = rhs_vf(momxb)%sf(j+1,k,l) + &
                            !                           0.5_wp * F_L * (1._wp/dx(j+1))

                            @:ATOMIC_ADD(rhs_vf(E_idx)%sf(j+1,k,l), (0.5_wp/dx(j+1) * vel_L * F_L * dt))
                            ! !$omp atomic update
                            ! !$acc atomic
                            ! rhs_vf(E_idx)%sf(j+1,k,l) = rhs_vf(E_idx)%sf(j+1,k,l) + &
                            !                           0.5_wp * vel_L * F_L * (1._wp/dx(j+1))

                            @:ATOMIC_ADD(rhs_vf(momxb)%sf(j,k,l), (-0.5_wp/dx(j) * F_L * dt))
                            ! !$omp atomic update
                            ! !$acc atomic
                            ! rhs_vf(momxb)%sf(j,k,l) = rhs_vf(momxb)%sf(j,k,l) - &
                            !                           0.5_wp * F_L * (1._wp/dx(j))

                            @:ATOMIC_ADD(rhs_vf(E_idx)%sf(j,k,l),(-0.5_wp/dx(j) * vel_L * F_L * dt))
                            ! !$omp atomic update
                            ! !$acc atomic
                            ! rhs_vf(E_idx)%sf(j,k,l) = rhs_vf(E_idx)%sf(j,k,l) - &
                            !                           0.5_wp * vel_L * F_L * (1._wp/dx(j))

                            F_L = (1._wp/60._wp) * (-3._wp * jac(j+2, k, l) + &
                                                27._wp * jac(j+1, k, l) + &
                                                47._wp * jac(j, k, l) -   &
                                                13._wp * jac(j-1, k, l) + &
                                                2._wp * jac(j-2, k, l))

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

                            @:ATOMIC_ADD(rhs_vf(momxb)%sf(j+1,k,l),(0.5_wp/dx(j+1) * F_L * dt))
                            ! !$omp atomic update
                            ! !$acc atomic
                            ! rhs_vf(momxb)%sf(j+1,k,l) = rhs_vf(momxb)%sf(j+1,k,l) + &
                            !                           0.5_wp * F_L * (1._wp/dx(j+1))

                            @:ATOMIC_ADD(rhs_vf(E_idx)%sf(j+1,k,l),(0.5_wp/dx(j+1) * vel_L * F_L * dt))
                            ! !$omp atomic update
                            ! !$acc atomic
                            ! rhs_vf(E_idx)%sf(j+1,k,l) = rhs_vf(E_idx)%sf(j+1,k,l) + &
                            !                           0.5_wp * vel_L * F_L * (1._wp/dx(j+1))

                            @:ATOMIC_ADD(rhs_vf(momxb)%sf(j,k,l),(-0.5_wp/dx(j) * F_L * dt))
                            ! !$omp atomic update
                            ! !$acc atomic
                            ! rhs_vf(momxb)%sf(j,k,l) = rhs_vf(momxb)%sf(j,k,l) - &
                            !                           0.5_wp * F_L * (1._wp/dx(j))

                            @:ATOMIC_ADD(rhs_vf(E_idx)%sf(j,k,l),(-0.5_wp/dx(j) * vel_L * F_L * dt))
                            ! !$omp atomic update
                            ! !$acc atomic
                            ! rhs_vf(E_idx)%sf(j,k,l) = rhs_vf(E_idx)%sf(j,k,l) - &
                            !                           0.5_wp * vel_L * F_L * (1._wp/dx(j))                   

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

                            @:ATOMIC_ADD(rhs_vf(momxb)%sf(j+1,k,l), (0.5_wp/dx(j+1) * F_L * dt))
                            ! !$omp atomic update
                            ! !$acc atomic
                            ! rhs_vf(momxb)%sf(j+1,k,l) = rhs_vf(momxb)%sf(j+1,k,l) + &
                            !                           0.5_wp * F_L * (1._wp/dx(j+1))

                            @:ATOMIC_ADD(rhs_vf(E_idx)%sf(j+1,k,l), (0.5_wp/dx(j+1) * vel_L * F_L * dt))
                            ! !$omp atomic update
                            ! !$acc atomic
                            ! rhs_vf(E_idx)%sf(j+1,k,l) = rhs_vf(E_idx)%sf(j+1,k,l) + &
                            !                           0.5_wp * vel_L * F_L * (1._wp/dx(j+1))

                            @:ATOMIC_ADD(rhs_vf(momxb)%sf(j,k,l), (-0.5_wp/dx(j) * F_L * dt))
                            ! !$omp atomic update
                            ! !$acc atomic
                            ! rhs_vf(momxb)%sf(j,k,l) = rhs_vf(momxb)%sf(j,k,l) - &
                            !                           0.5_wp * F_L * (1._wp/dx(j))

                            @:ATOMIC_ADD(rhs_vf(E_idx)%sf(j,k,l),(-0.5_wp/dx(j) * vel_L * F_L * dt))
                            ! !$omp atomic update
                            ! !$acc atomic
                            ! rhs_vf(E_idx)%sf(j,k,l) = rhs_vf(E_idx)%sf(j,k,l) - &
                            !                           0.5_wp * vel_L * F_L * (1._wp/dx(j))

                            F_L = (1._wp/60._wp) * (-3._wp * jac(j+2, k, l) + &
                                                27._wp * jac(j+1, k, l) + &
                                                47._wp * jac(j, k, l) -   &
                                                13._wp * jac(j-1, k, l) + &
                                                2._wp * jac(j-2, k, l))

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

                            @:ATOMIC_ADD(rhs_vf(momxb)%sf(j+1,k,l),(0.5_wp/dx(j+1) * F_L * dt))
                            ! !$omp atomic update
                            ! !$acc atomic
                            ! rhs_vf(momxb)%sf(j+1,k,l) = rhs_vf(momxb)%sf(j+1,k,l) + &
                            !                           0.5_wp * F_L * (1._wp/dx(j+1))

                            @:ATOMIC_ADD(rhs_vf(E_idx)%sf(j+1,k,l),(0.5_wp/dx(j+1) * vel_L * F_L * dt))
                            ! !$omp atomic update
                            ! !$acc atomic
                            ! rhs_vf(E_idx)%sf(j+1,k,l) = rhs_vf(E_idx)%sf(j+1,k,l) + &
                            !                           0.5_wp * vel_L * F_L * (1._wp/dx(j+1))

                            @:ATOMIC_ADD(rhs_vf(momxb)%sf(j,k,l),(-0.5_wp/dx(j) * F_L * dt))
                            ! !$omp atomic update
                            ! !$acc atomic
                            ! rhs_vf(momxb)%sf(j,k,l) = rhs_vf(momxb)%sf(j,k,l) - &
                            !                           0.5_wp * F_L * (1._wp/dx(j))

                            @:ATOMIC_ADD(rhs_vf(E_idx)%sf(j,k,l),(-0.5_wp/dx(j) * vel_L * F_L * dt))
                            ! !$omp atomic update
                            ! !$acc atomic
                            ! rhs_vf(E_idx)%sf(j,k,l) = rhs_vf(E_idx)%sf(j,k,l) - &
                            !                           0.5_wp * vel_L * F_L * (1._wp/dx(j)) 

                        end do 
                    end do 
                end do                
            end if
        end if

    end subroutine s_igr_sigma

    subroutine s_igr_riemann_solver(q_prim_vf, rhs_vf, idir)

        type(scalar_field_half), &
            dimension(sys_size), &
            intent(inout) :: rhs_vf
        type(scalar_field_half), &
            dimension(sys_size), &
            intent(inout) :: q_prim_vf
        integer, intent(in) :: idir

        real(wp) :: rho_L, gamma_L, pi_inf_L, a_L, cfl, pres_L, F_L, E_L, pres_R, a_R, rho_R, gamma_R, pi_inf_R, E_R, mu_L, mu_R
        real(wp), dimension(num_fluids) :: alpha_rho_L, alpha_L, alpha_R, alpha_rho_R
        real(wp), dimension(num_dims) :: vel_L, vel_R
        real(wp), dimension(-1:1) :: rho_sf_small
        real(wp), dimension(num_dims,num_dims) :: dvel
        real(wp), dimension(3) :: vflux_L_arr, vflux_R_arr
        real(wp), dimension(num_dims) :: dvel_small

        if (idir == 1) then
            if(p == 0) then
                !$omp target teams loop collapse(3) private(rho_L,gamma_L,pi_inf_L,mu_L,a_L,vel_L,vel_R,pres_L,alpha_L,alpha_R,alpha_rho_L,cfl,dvel,F_L,E_L,mu_R,rho_sf_small,alpha_rho_R,vflux_L_arr,vflux_R_arr,dvel_small)
                !$acc parallel loop collapse(3) gang vector default(present) private(rho_L,gamma_L,pi_inf_L,mu_L,a_L,vel_L,vel_R,pres_L,alpha_L,alpha_R,alpha_rho_L,cfl,dvel,F_L,E_L,mu_R,rho_sf_small,alpha_rho_R,vflux_L_arr,vflux_R_arr,dvel_small)                
                do l = 0, p
                    do k = 0, n
                        do j = -buff_size+5, m+buff_size-5

                            dvel = 0._wp
                            vflux_L_arr = 0._wp
                            vflux_R_arr = 0._wp

                            !DIR$ unroll 6
                            !$acc loop seq 
                            do q = -2, 3
                                dvel_small = 0._wp
                                !x-direction contributions
                                !$acc loop seq 
                                do i = -1, 1
                                    rho_L = 0._wp
                                    !$acc loop seq 
                                    do r = 1, num_fluids
                                        rho_L = rho_L + q_prim_vf(r)%sf(j+i+q,k,l)
                                    end do
                                    rho_sf_small(i) = rho_L
                                end do

                                dvel_small(1) = (1/(2._wp*dx(j))) * ( &
                                    1._wp*q_prim_vf(momxb)%sf(j+1+q,k,l)/rho_sf_small(1) - &
                                    1._wp*q_prim_vf(momxb)%sf(j-1+q,k,l)/rho_sf_small(-1))
                                dvel_small(2) = (1/(2._wp*dx(j))) * ( &
                                    q_prim_vf(momxb+1)%sf(j+1+q,k,l)/rho_sf_small(1) - &
                                    q_prim_vf(momxb+1)%sf(j-1+q,k,l)/rho_sf_small(-1))

                                if (q == 0) dvel(:,1) = dvel_small
                                if (q > -2 .and. viscous) then
                                    vflux_L_arr(1) = vflux_L_arr(1) + coeff_L(q)*(dvel_small(2))
                                    vflux_L_arr(3) = vflux_L_arr(3) + coeff_L(q)*(4._wp*dvel_small(1))/3._wp
                                end if
                                if (q < 3 .and. viscous) then
                                    vflux_R_arr(1) = vflux_R_arr(1) + coeff_R(q)*(dvel_small(2))
                                    vflux_R_arr(3) = vflux_R_arr(3) + coeff_R(q)*(4._wp*dvel_small(1))/3._wp
                                end if

                                !y-direction contributions
                                !$acc loop seq 
                                do i = -1, 1
                                    rho_L = 0._wp
                                    !$acc loop seq 
                                    do r = 1, num_fluids
                                        rho_L = rho_L + q_prim_vf(r)%sf(j+q,k+i,l)
                                    end do
                                    rho_sf_small(i) = rho_L
                                end do

                                dvel_small(1) = (1/(2._wp*dy(k))) * ( &
                                    q_prim_vf(momxb)%sf(j+q,k+1,l)/rho_sf_small(1) - &
                                    q_prim_vf(momxb)%sf(j+q,k-1,l)/rho_sf_small(-1))
                                dvel_small(2) = (1/(2._wp*dy(k))) * ( &
                                    q_prim_vf(momxb+1)%sf(j+q,k+1,l)/rho_sf_small(1) - &
                                    q_prim_vf(momxb+1)%sf(j+q,k-1,l)/rho_sf_small(-1))


                                if (q == 0) dvel(:,2) = dvel_small

                                if (q > -2 .and. viscous) then
                                    vflux_L_arr(1) = vflux_L_arr(1) + coeff_L(q)*(dvel_small(1))
                                    vflux_L_arr(3) = vflux_L_arr(3) + coeff_L(q)*(-2._wp*dvel_small(2))/3._wp
                                end if
                                if (q < 3 .and. viscous) then
                                    vflux_R_arr(1) = vflux_R_arr(1) + coeff_R(q)*(dvel_small(1))
                                    vflux_R_arr(3) = vflux_R_arr(3) + coeff_R(q)*(-2._wp*dvel_small(2))/3._wp
                                end if

                                if (q == 0) then
                                    jac_rhs(j,k,l) = alf_igr * (2._wp*(dvel(1,2)*dvel(2,1)) &
                                        + dvel(1,1)**2_wp + dvel(2,2)**2_wp &
                                        + (dvel(1,1) + dvel(2,2))**2_wp)
                                end if
                            end do

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
                            end do 

                            if(num_fluids > 1) then 
                                !$acc loop seq 
                                do i = 1, num_fluids - 1
                                    alpha_L(i) =  (1._wp/60._wp) * (-3._wp * q_prim_vf(E_idx+i)%sf(j-1, k, l) + &
                                                    27._wp * q_prim_vf(E_idx+i)%sf(j, k, l) + &
                                                    47._wp * q_prim_vf(E_idx+i)%sf(j+1, k, l) -   &
                                                    13._wp * q_prim_vf(E_idx+i)%sf(j+2, k, l) + &
                                                    2._wp * q_prim_vf(E_idx+i)%sf(j+3, k, l))
                                    alpha_R(i) = (1._wp/60._wp) * (-3._wp * q_prim_vf(E_idx+i)%sf(j+2, k, l) + &
                                                            27._wp * q_prim_vf(E_idx+i)%sf(j+1, k, l) + &
                                                            47._wp * q_prim_vf(E_idx+i)%sf(j, k, l) -   &
                                                            13._wp * q_prim_vf(E_idx+i)%sf(j-1, k, l) + &
                                                            2._wp * q_prim_vf(E_idx+i)%sf(j-2, k, l))
                                end do

                                alpha_L(num_fluids) = 1._wp - sum(alpha_L(1:num_fluids - 1))
                                alpha_R(num_fluids) = 1._wp - sum(alpha_R(1:num_fluids - 1))    
                            else 
                                alpha_L(1) = 1._wp 
                                alpha_R(1) = 1._wp
                            end if

                            rho_L = sum(alpha_rho_L)
                            gamma_L = sum(alpha_L*gammas)
                            pi_inf_L = sum(alpha_L*pi_infs)

                            rho_R = sum(alpha_rho_R)
                            gamma_R = sum(alpha_R*gammas)
                            pi_inf_R = sum(alpha_R*pi_infs)

                            !$acc loop seq 
                            do i = 1, num_dims
                                vel_L(i) = (1._wp/60._wp) * (-3._wp * q_prim_vf(momxb+i-1)%sf(j-1,k,l) + &
                                                27._wp * q_prim_vf(momxb+i-1)%sf(j,k,l) + &
                                                47._wp * q_prim_vf(momxb+i-1)%sf(j+1,k,l) -   &
                                                13._wp * q_prim_vf(momxb+i-1)%sf(j+2,k,l) + &
                                                2._wp * q_prim_vf(momxb+i-1)%sf(j+3,k,l)) / rho_L
                            end do 

                            !$acc loop seq 
                            do i = 1, num_dims
                                vel_R(i) = (1._wp/60._wp) * (-3._wp * q_prim_vf(momxb+i-1)%sf(j+2,k,l) + &
                                                27._wp * q_prim_vf(momxb+i-1)%sf(j+1,k,l) + &
                                                47._wp * q_prim_vf(momxb+i-1)%sf(j,k,l) -   &
                                                13._wp * q_prim_vf(momxb+i-1)%sf(j-1,k,l) + &
                                                2._wp * q_prim_vf(momxb+i-1)%sf(j-2,k,l)) / rho_R
                            end do

                            if(viscous) then

                                if(num_fluids > 1) then 
                                    mu_L = 0._wp
                                    mu_R = 0._wp
                                    !$acc loop seq 
                                    do i = 1, num_fluids - 1
                                        mu_L = (1._wp/60._wp) * (-3._wp * q_prim_vf(E_idx+i)%sf(j-1, k, l) + &
                                                27._wp * q_prim_vf(E_idx+i)%sf(j, k, l) + &
                                                47._wp * q_prim_vf(E_idx+i)%sf(j+1, k, l) -   &
                                                13._wp * q_prim_vf(E_idx+i)%sf(j+2, k, l) + &
                                                2._wp * q_prim_vf(E_idx+i)%sf(j+3, k, l)) / Res(1, i) + mu_L

                                        mu_R = (1._wp/60._wp) * (-3._wp * q_prim_vf(E_idx+i)%sf(j+2, k, l) + &
                                                27._wp * q_prim_vf(E_idx+i)%sf(j+1, k, l) + &
                                                47._wp * q_prim_vf(E_idx+i)%sf(j, k, l) -   &
                                                13._wp * q_prim_vf(E_idx+i)%sf(j-1, k, l) + &
                                                2._wp * q_prim_vf(E_idx+i)%sf(j-2, k, l)) / Res(1, i) + mu_R
                                    end do

                                    mu_L = 1._wp / Res(1, num_fluids) + mu_L
                                    mu_R = 1._wp / Res(1, num_fluids) + mu_R

                                    !$acc loop seq 
                                    do i = 1, num_fluids - 1
                                        mu_L = -(1._wp/60._wp) * (-3._wp * q_prim_vf(E_idx+i)%sf(j-1, k, l) + &
                                                27._wp * q_prim_vf(E_idx+i)%sf(j, k, l) + &
                                                47._wp * q_prim_vf(E_idx+i)%sf(j+1, k, l) -   &
                                                13._wp * q_prim_vf(E_idx+i)%sf(j+2, k, l) + &
                                                2._wp * q_prim_vf(E_idx+i)%sf(j+3, k, l)) / Res(1, num_fluids) + mu_L

                                        mu_R = -(1._wp/60._wp) * (-3._wp * q_prim_vf(E_idx+i)%sf(j+2, k, l) + &
                                                27._wp * q_prim_vf(E_idx+i)%sf(j+1, k, l) + &
                                                47._wp * q_prim_vf(E_idx+i)%sf(j, k, l) -   &
                                                13._wp * q_prim_vf(E_idx+i)%sf(j-1, k, l) + &
                                                2._wp * q_prim_vf(E_idx+i)%sf(j-2, k, l)) / Res(1, num_fluids) + mu_R
                                    end do
                                else
                                    mu_L = 1._wp / Res(1, 1) 
                                    mu_R = 1._wp / Res(1, 1) 
                                end if

                                @:ATOMIC_ADD(rhs_vf(momxb+1)%sf(j+1,k,l),(- 0.5_wp*mu_L*vflux_L_arr(1)*(1._wp/dx(j+1))*dt))
                                !!$omp atomic update
                                !!$acc atomic
                                !rhs_vf(momxb+1)%sf(j+1,k,l) = rhs_vf(momxb+1)%sf(j+1,k,l) - 0.5_wp*mu_L*vflux_L_arr(1)*(1._wp/dx(j+1))
                                @:ATOMIC_ADD(rhs_vf(E_idx)%sf(j+1,k,l),(- 0.5_wp*mu_L*vflux_L_arr(1)*vel_L(2)*(1._wp/dx(j+1))*dt))
                                !!$omp atomic update
                                !!$acc atomic
                                !rhs_vf(E_idx)%sf(j+1,k,l) = rhs_vf(E_idx)%sf(j+1,k,l) - 0.5_wp*mu_L*vflux_L_arr(1)*vel_L(2)*(1._wp/dx(j+1))
                        
                                @:ATOMIC_ADD(rhs_vf(momxb+1)%sf(j,k,l),(0.5_wp*mu_L*vflux_L_arr(1)*(1._wp/dx(j))*dt))
                                !!$omp atomic update
                                !!$acc atomic
                                !rhs_vf(momxb+1)%sf(j,k,l) = rhs_vf(momxb+1)%sf(j,k,l) + 0.5_wp*mu_L*vflux_L_arr(1)*(1._wp/dx(j))
                                @:ATOMIC_ADD(rhs_vf(E_idx)%sf(j,k,l),(0.5_wp*mu_L*vflux_L_arr(1)*vel_L(2)*(1._wp/dx(j))*dt))
                                !!$omp atomic update
                                !!$acc atomic
                                !rhs_vf(E_idx)%sf(j,k,l) = rhs_vf(E_idx)%sf(j,k,l) + 0.5_wp*mu_L*vflux_L_arr(1)*vel_L(2)*(1._wp/dx(j))                    
                             
                                @:ATOMIC_ADD(rhs_vf(momxb+1)%sf(j+1,k,l),(- 0.5_wp*mu_R*vflux_R_arr(1)*(1._wp/dx(j+1))*dt))
                                !!$omp atomic update
                                !!$acc atomic
                                !rhs_vf(momxb+1)%sf(j+1,k,l) = rhs_vf(momxb+1)%sf(j+1,k,l) - 0.5_wp*mu_R*vflux_R_arr(1)*(1._wp/dx(j+1))
                                @:ATOMIC_ADD(rhs_vf(E_idx)%sf(j+1,k,l),(- 0.5_wp*mu_R*vflux_R_arr(1)*vel_R(2)*(1._wp/dx(j+1))*dt))
                                !!$omp atomic update
                                !!$acc atomic
                                !rhs_vf(E_idx)%sf(j+1,k,l) = rhs_vf(E_idx)%sf(j+1,k,l) - 0.5_wp*mu_R*vflux_R_arr(1)*vel_R(2)*(1._wp/dx(j+1))
                             
                                @:ATOMIC_ADD(rhs_vf(momxb+1)%sf(j,k,l),(0.5_wp*mu_R*vflux_R_arr(1)*(1._wp/dx(j))*dt))
                                !!$omp atomic update
                                !!$acc atomic
                                !rhs_vf(momxb+1)%sf(j,k,l) = rhs_vf(momxb+1)%sf(j,k,l) + 0.5_wp*mu_R*vflux_R_arr(1)*(1._wp/dx(j))
                                @:ATOMIC_ADD(rhs_vf(E_idx)%sf(j,k,l),(0.5_wp*mu_R*vflux_R_arr(1)*vel_R(2)*(1._wp/dx(j))*dt))
                                !!$omp atomic update
                                !!$acc atomic
                                !rhs_vf(E_idx)%sf(j,k,l) = rhs_vf(E_idx)%sf(j,k,l) + 0.5_wp*mu_R*vflux_R_arr(1)*vel_R(2)*(1._wp/dx(j))

                                @:ATOMIC_ADD(rhs_vf(momxb)%sf(j+1,k,l),(- 0.5_wp*mu_L*vflux_L_arr(3)*(1._wp/dx(j+1))*dt))
                                !!$omp atomic update
                                !!$acc atomic
                                !rhs_vf(momxb)%sf(j+1,k,l) = rhs_vf(momxb)%sf(j+1,k,l) - 0.5_wp*mu_L*vflux_L_arr(3)*(1._wp/dx(j+1))
                                @:ATOMIC_ADD(rhs_vf(E_idx)%sf(j+1,k,l),(- 0.5_wp*mu_L*vflux_L_arr(3)*vel_L(1)*(1._wp/dx(j+1))*dt))
                                !!$omp atomic update
                                !!$acc atomic
                                !rhs_vf(E_idx)%sf(j+1,k,l) = rhs_vf(E_idx)%sf(j+1,k,l) - 0.5_wp*mu_L*vflux_L_arr(3)*vel_L(1)*(1._wp/dx(j+1)) 

                                @:ATOMIC_ADD(rhs_vf(momxb)%sf(j,k,l),(0.5_wp*mu_L*vflux_L_arr(3)*(1._wp/dx(j))*dt))
                                !!$omp atomic update
                                !!$acc atomic
                                !rhs_vf(momxb)%sf(j,k,l) = rhs_vf(momxb)%sf(j,k,l) + 0.5_wp*mu_L*vflux_L_arr(3)*(1._wp/dx(j))
                                @:ATOMIC_ADD(rhs_vf(E_idx)%sf(j,k,l),(0.5_wp*mu_L*vflux_L_arr(3)*vel_L(1)*(1._wp/dx(j))*dt))
                                !!$omp atomic update
                                !!$acc atomic
                                !rhs_vf(E_idx)%sf(j,k,l) = rhs_vf(E_idx)%sf(j,k,l) + 0.5_wp*mu_L*vflux_L_arr(3)*vel_L(1)*(1._wp/dx(j))
                            
                                @:ATOMIC_ADD(rhs_vf(momxb)%sf(j+1,k,l),(- 0.5_wp*mu_R*vflux_R_arr(3)*(1._wp/dx(j+1))*dt))
                                !!$omp atomic update
                                !!$acc atomic
                                !rhs_vf(momxb)%sf(j+1,k,l) = rhs_vf(momxb)%sf(j+1,k,l) - 0.5_wp*mu_R*vflux_R_arr(3)*(1._wp/dx(j+1))
                                @:ATOMIC_ADD(rhs_vf(E_idx)%sf(j+1,k,l),(- 0.5_wp*mu_R*vflux_R_arr(3)*vel_R(1)*(1._wp/dx(j+1))*dt))
                                !!$omp atomic update
                                !!$acc atomic
                                !rhs_vf(E_idx)%sf(j+1,k,l) = rhs_vf(E_idx)%sf(j+1,k,l) - 0.5_wp*mu_R*vflux_R_arr(3)*vel_R(1)*(1._wp/dx(j+1))
                            
                                @:ATOMIC_ADD(rhs_vf(momxb)%sf(j,k,l),(0.5_wp*mu_R*vflux_R_arr(3)*(1._wp/dx(j))*dt))
                                !!$omp atomic update
                                !!$acc atomic
                                !rhs_vf(momxb)%sf(j,k,l) = rhs_vf(momxb)%sf(j,k,l) + 0.5_wp*mu_R*vflux_R_arr(3)*(1._wp/dx(j))
                                @:ATOMIC_ADD(rhs_vf(E_idx)%sf(j,k,l),(0.5_wp*mu_R*vflux_R_arr(3)*vel_R(1)*(1._wp/dx(j))*dt))
                                !!$omp atomic update
                                !!$acc atomic
                                !rhs_vf(E_idx)%sf(j,k,l) = rhs_vf(E_idx)%sf(j,k,l) + 0.5_wp*mu_R*vflux_R_arr(3)*vel_R(1)*(1._wp/dx(j))
                            
                            end if

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

                            pres_L = (E_L - pi_inf_L - 0.5_wp*rho_L*(vel_L(1)**2._wp + vel_L(2)**2._wp))/gamma_L                    

                            pres_R = (E_R - pi_inf_R - 0.5_wp*rho_R*(vel_R(1)**2._wp + vel_R(2)**2._wp))/gamma_R 

                            a_L = sqrt((pres_L*(1._wp/gamma_L + 1._wp) + pi_inf_L / gamma_L) / rho_L)
                            a_R = sqrt((pres_R*(1._wp/gamma_R + 1._wp) + pi_inf_R / gamma_R) / rho_R)

                            cfl = (max(sqrt(vel_L(1)**2._wp + vel_L(2)**2._wp),sqrt(vel_R(1)**2._wp + vel_R(2)**2._wp)) + max(a_L,a_R))

                            !$acc loop seq
                            do i = 1, num_fluids
                                @:ATOMIC_ADD(rhs_vf(i)%sf(j+1,k,l),(0.5_wp * (alpha_rho_L(i) * &
                                    vel_L(1))*(1._wp/dx(j+1)) - &
                                    0.5_wp*cfl *(alpha_rho_L(i))*(1._wp/dx(j+1)))*dt)
                                ! !$omp atomic update
                                ! !$acc atomic
                                ! rhs_vf(i)%sf(j+1,k,l) = rhs_vf(i)%sf(j+1,k,l) + &
                                !     (0.5_wp * (alpha_rho_L(i) * &
                                !     vel_L(1))*(1._wp/dx(j+1)) - &
                                !     0.5_wp*cfl *(alpha_rho_L(i))*(1._wp/dx(j+1)))

                                @:ATOMIC_ADD(rhs_vf(i)%sf(j,k,l), &
                                    (-0.5_wp * (alpha_rho_L(i) * &
                                    vel_L(1))*(1._wp/dx(j)) + &
                                    0.5_wp*cfl * (alpha_rho_L(i))*(1._wp/dx(j)))*dt)
                                ! !$omp atomic update
                                ! !$acc atomic
                                ! rhs_vf(i)%sf(j,k,l) = rhs_vf(i)%sf(j,k,l) - &
                                !     (0.5_wp * (alpha_rho_L(i) * &
                                !     vel_L(1))*(1._wp/dx(j)) - &
                                !     0.5_wp*cfl * (alpha_rho_L(i))*(1._wp/dx(j))) 
                            end do

                            if(num_fluids > 1) then 
                                !$acc loop seq
                                do i = 1, num_fluids - 1
                                    @:ATOMIC_ADD(rhs_vf(advxb+i-1)%sf(j+1,k,l),(0.5_wp * (alpha_L(i) * &
                                        vel_L(1))*(1._wp/dx(j+1)) - &
                                        0.5_wp*cfl*(alpha_L(i))*(1._wp/dx(j+1)))*dt)
                                    ! !$omp atomic update
                                    ! !$acc atomic
                                    ! rhs_vf(advxb+i-1)%sf(j+1,k,l) = rhs_vf(advxb+i-1)%sf(j+1,k,l) + &
                                    !     (0.5_wp * (alpha_L(i) * &
                                    !     vel_L(1))*(1._wp/dx(j+1)) - &
                                    !     0.5_wp*cfl*(alpha_L(i))*(1._wp/dx(j+1)))

                                    @:ATOMIC_ADD(rhs_vf(advxb+i-1)%sf(j+1,k,l),(-0.5_wp * q_prim_vf(advxb+i-1)%sf(j+1,k,l) * vel_L(1)*(1._wp/dx(j+1)))*dt)
                                    ! !$omp atomic update
                                    ! !$acc atomic
                                    ! rhs_vf(advxb+i-1)%sf(j+1,k,l) = rhs_vf(advxb+i-1)%sf(j+1,k,l) &
                                    ! - (0.5_wp * q_prim_vf(advxb+i-1)%sf(j+1,k,l) * vel_L(1)*(1._wp/dx(j+1)))

                                    @:ATOMIC_ADD(rhs_vf(advxb+i-1)%sf(j,k,l),(-0.5_wp * (alpha_L(i) * &
                                        vel_L(1))*(1._wp/dx(j)) + &
                                        0.5_wp*cfl*(alpha_L(i))*(1._wp/dx(j)))*dt)
                                    ! !$omp atomic update
                                    ! !$acc atomic
                                    ! rhs_vf(advxb+i-1)%sf(j,k,l) = rhs_vf(advxb+i-1)%sf(j,k,l) - &
                                    !     (0.5_wp * (alpha_L(i) * &
                                    !     vel_L(1))*(1._wp/dx(j)) - &
                                    !     0.5_wp*cfl*(alpha_L(i))*(1._wp/dx(j)))

                                    @:ATOMIC_ADD(rhs_vf(advxb+i-1)%sf(j,k,l),(0.5_wp * q_prim_vf(advxb+i-1)%sf(j,k,l)  * vel_L(1)*(1._wp/dx(j)))*dt)
                                    ! !$omp atomic update
                                    ! !$acc atomic
                                    ! rhs_vf(advxb+i-1)%sf(j,k,l) = rhs_vf(advxb+i-1)%sf(j,k,l) &
                                    ! + (0.5_wp * q_prim_vf(advxb+i-1)%sf(j,k,l)  * vel_L(1)*(1._wp/dx(j)))
                                end do
                            end if

                            @:ATOMIC_ADD(rhs_vf(momxb)%sf(j+1,k,l),(0.5_wp* (rho_L * (vel_L(1))**2.0 + &
                                 pres_L)*(1._wp/dx(j+1)) - &
                                 0.5_wp*cfl * (rho_L*vel_L(1))*(1._wp/dx(j+1)))*dt)
                            ! !$omp atomic update
                            ! !$acc atomic
                            ! rhs_vf(momxb)%sf(j+1,k,l) = rhs_vf(momxb)%sf(j+1,k,l) + &
                            !      (0.5_wp* (rho_L * (vel_L(1))**2.0 + &
                            !      pres_L)*(1._wp/dx(j+1)) - &
                            !      0.5_wp*cfl * (rho_L*vel_L(1))*(1._wp/dx(j+1)))

                            @:ATOMIC_ADD(rhs_vf(momxb+1)%sf(j+1, k, l),(0.5_wp * rho_L * vel_L(1)*vel_L(2)*(1._wp/dx(j+1)) - &
                                0.5_wp*cfl * (rho_L*vel_L(2))*(1._wp/dx(j+1)))*dt)
                            ! !$omp atomic update
                            ! !$acc atomic
                            ! rhs_vf(momxb+1)%sf(j+1, k, l) =  rhs_vf(momxb+1)%sf(j+1, k, l) + &
                            !     (0.5_wp * rho_L * vel_L(1)*vel_L(2)*(1._wp/dx(j+1)) - &
                            !     0.5_wp*cfl * (rho_L*vel_L(2))*(1._wp/dx(j+1)))

                            @:ATOMIC_ADD(rhs_vf(E_idx)%sf(j+1, k, l),(0.5_wp * (vel_L(1) * (E_L + &
                                pres_L) )*(1._wp/dx(j+1)) - &
                                0.5_wp*cfl * (E_L)*(1._wp/dx(j+1)))*dt)
                            ! !$omp atomic update
                            ! !$acc atomic
                            ! rhs_vf(E_idx)%sf(j+1, k, l) = rhs_vf(E_idx)%sf(j+1, k, l) +  &
                            !     (0.5_wp * (vel_L(1) * (E_L + &
                            !     pres_L) )*(1._wp/dx(j+1)) - &
                            !     0.5_wp*cfl * (E_L)*(1._wp/dx(j+1)))

                            @:ATOMIC_ADD(rhs_vf(momxb)%sf(j,k,l),(-0.5_wp* (rho_L * (vel_L(1))**2.0 + &
                                 pres_L)*(1._wp/dx(j)) + &
                                 0.5_wp*cfl * (rho_L*vel_L(1))*(1._wp/dx(j)))*dt)
                            ! !$omp atomic update
                            ! !$acc atomic
                            ! rhs_vf(momxb)%sf(j,k,l) = rhs_vf(momxb)%sf(j,k,l) - &
                            !      (0.5_wp* (rho_L * (vel_L(1))**2.0 + &
                            !      pres_L)*(1._wp/dx(j)) - &
                            !      0.5_wp*cfl * (rho_L*vel_L(1))*(1._wp/dx(j)))

                            @:ATOMIC_ADD(rhs_vf(momxb+1)%sf(j, k, l),(-0.5_wp * rho_L * vel_L(1)*vel_L(2)*(1._wp/dx(j)) + &
                                0.5_wp*cfl * (rho_L*vel_L(2))*(1._wp/dx(j)))*dt)
                            ! !$omp atomic update
                            ! !$acc atomic
                            ! rhs_vf(momxb+1)%sf(j, k, l) =  rhs_vf(momxb+1)%sf(j, k, l) - &
                            !     (0.5_wp * rho_L * vel_L(1)*vel_L(2)*(1._wp/dx(j)) - &
                            !     0.5_wp*cfl * (rho_L*vel_L(2))*(1._wp/dx(j)))

                            @:ATOMIC_ADD(rhs_vf(E_idx)%sf(j, k, l),(-0.5_wp * (vel_L(1) * (E_L + &
                                pres_L) )*(1._wp/dx(j)) + &
                                0.5_wp*cfl * (E_L)*(1._wp/dx(j)))*dt)
                            ! !$omp atomic update
                            ! !$acc atomic
                            !  rhs_vf(E_idx)%sf(j, k, l) = rhs_vf(E_idx)%sf(j, k, l) -  &
                            !     (0.5_wp * (vel_L(1) * (E_L + &
                            !     pres_L) )*(1._wp/dx(j)) - &
                            !     0.5_wp*cfl * (E_L)*(1._wp/dx(j)))

                            !$acc loop seq
                            do i = 1, num_fluids
                                @:ATOMIC_ADD(rhs_vf(i)%sf(j+1,k,l),(0.5_wp * (alpha_rho_R(i) * &
                                    vel_R(1))*(1._wp/dx(j+1)) + &
                                    0.5_wp*cfl * (alpha_rho_R(i))*(1._wp/dx(j+1)))*dt)
                                ! !$omp atomic update
                                ! !$acc atomic
                                ! rhs_vf(i)%sf(j+1,k,l) = rhs_vf(i)%sf(j+1,k,l) + &
                                !     (0.5_wp * (alpha_rho_R(i) * &
                                !     vel_R(1))*(1._wp/dx(j+1)) + &
                                !     0.5_wp*cfl * (alpha_rho_R(i))*(1._wp/dx(j+1)))

                                @:ATOMIC_ADD(rhs_vf(i)%sf(j,k,l),(-0.5_wp * (alpha_rho_R(i) * &
                                    vel_R(1))*(1._wp/dx(j)) - &
                                    0.5_wp*cfl * (alpha_rho_R(i))*(1._wp/dx(j)))*dt)
                                ! !$omp atomic update
                                ! !$acc atomic
                                ! rhs_vf(i)%sf(j,k,l) = rhs_vf(i)%sf(j,k,l) - &
                                !     (0.5_wp * (alpha_rho_R(i) * &
                                !     vel_R(1))*(1._wp/dx(j)) + &
                                !     0.5_wp*cfl * (alpha_rho_R(i))*(1._wp/dx(j)))
                            end do

                            if(num_fluids > 1) then 
                                !$acc loop seq
                                do i = 1, num_fluids - 1
                                    @:ATOMIC_ADD(rhs_vf(advxb+i-1)%sf(j+1,k,l),(0.5_wp * (alpha_R(i) * &
                                        vel_R(1))*(1._wp/dx(j+1)) + &
                                        0.5_wp*cfl*(alpha_R(i))*(1._wp/dx(j+1)))*dt)
                                    ! !$omp atomic update
                                    ! !$acc atomic
                                    ! rhs_vf(advxb+i-1)%sf(j+1,k,l) = rhs_vf(advxb+i-1)%sf(j+1,k,l) + &
                                    !     (0.5_wp * (alpha_R(i) * &
                                    !     vel_R(1))*(1._wp/dx(j+1)) + &
                                    !     0.5_wp*cfl*(alpha_R(i))*(1._wp/dx(j+1)))

                                    @:ATOMIC_ADD(rhs_vf(advxb+i-1)%sf(j+1,k,l),(-0.5_wp * q_prim_vf(advxb+i-1)%sf(j+1,k,l)  * vel_R(1)*(1._wp/dx(j+1)))*dt)
                                    ! !$omp atomic update
                                    ! !$acc atomic
                                    ! rhs_vf(advxb+i-1)%sf(j+1,k,l) = rhs_vf(advxb+i-1)%sf(j+1,k,l) &
                                    ! - (0.5_wp * q_prim_vf(advxb+i-1)%sf(j+1,k,l)  * vel_R(1)*(1._wp/dx(j+1)))

                                    @:ATOMIC_ADD(rhs_vf(advxb+i-1)%sf(j,k,l),(-0.5_wp * (alpha_R(i) * &
                                        vel_R(1))*(1._wp/dx(j)) - &
                                        0.5_wp*cfl*(alpha_R(i))*(1._wp/dx(j)))*dt)
                                    ! !$omp atomic update
                                    ! !$acc atomic
                                    ! rhs_vf(advxb+i-1)%sf(j,k,l) = rhs_vf(advxb+i-1)%sf(j,k,l) - &
                                    !     (0.5_wp * (alpha_R(i) * &
                                    !     vel_R(1))*(1._wp/dx(j)) + &
                                    !     0.5_wp*cfl*(alpha_R(i))*(1._wp/dx(j)))

                                    @:ATOMIC_ADD(rhs_vf(advxb+i-1)%sf(j,k,l),(0.5_wp * q_prim_vf(advxb+i-1)%sf(j,k,l)  * vel_R(1)*(1._wp/dx(j)))*dt)
                                    ! !$omp atomic update
                                    ! !$acc atomic
                                    ! rhs_vf(advxb+i-1)%sf(j,k,l) = rhs_vf(advxb+i-1)%sf(j,k,l) &
                                    ! + (0.5_wp * q_prim_vf(advxb+i-1)%sf(j,k,l)  * vel_R(1)*(1._wp/dx(j)))
                                end do
                            end if

                            @:ATOMIC_ADD(rhs_vf(momxb)%sf(j+1,k,l),(0.5_wp* (rho_R * (vel_R(1))**2.0 + &
                                 pres_R)*(1._wp/dx(j+1)) + &
                                 0.5_wp*cfl * (rho_R*vel_R(1))*(1._wp/dx(j+1)))*dt)
                            ! !$omp atomic update
                            ! !$acc atomic
                            ! rhs_vf(momxb)%sf(j+1,k,l) = rhs_vf(momxb)%sf(j+1,k,l) + &
                            !      (0.5_wp* (rho_R * (vel_R(1))**2.0 + &
                            !      pres_R)*(1._wp/dx(j+1)) + &
                            !      0.5_wp*cfl * (rho_R*vel_R(1))*(1._wp/dx(j+1)))

                            @:ATOMIC_ADD(rhs_vf(momxb+1)%sf(j+1, k, l),(0.5_wp * rho_R * vel_R(1)*vel_R(2)*(1._wp/dx(j+1)) + &
                                0.5_wp*cfl * (rho_R*vel_R(2))*(1._wp/dx(j+1)))*dt)
                            ! !$omp atomic update
                            ! !$acc atomic
                            ! rhs_vf(momxb+1)%sf(j+1, k, l) =  rhs_vf(momxb+1)%sf(j+1, k, l) + &
                            !     (0.5_wp * rho_R * vel_R(1)*vel_R(2)*(1._wp/dx(j+1)) + &
                            !     0.5_wp*cfl * (rho_R*vel_R(2))*(1._wp/dx(j+1)))

                            @:ATOMIC_ADD(rhs_vf(E_idx)%sf(j+1, k, l),(0.5_wp * (vel_R(1) * (E_R + &
                                pres_R) )*(1._wp/dx(j+1)) + &
                                0.5_wp*cfl * (E_R)*(1._wp/dx(j+1)))*dt)
                            ! !$omp atomic update
                            ! !$acc atomic
                            !  rhs_vf(E_idx)%sf(j+1, k, l) = rhs_vf(E_idx)%sf(j+1, k, l) +  &
                            !     (0.5_wp * (vel_R(1) * (E_R + &
                            !     pres_R) )*(1._wp/dx(j+1)) + &
                            !     0.5_wp*cfl * (E_R)*(1._wp/dx(j+1)))

                            @:ATOMIC_ADD(rhs_vf(momxb)%sf(j,k,l),(-0.5_wp* (rho_R * (vel_R(1))**2.0 + &
                                 pres_R)*(1._wp/dx(j)) - &
                                 0.5_wp*cfl * (rho_R*vel_R(1))*(1._wp/dx(j)))*dt)
                            ! !$omp atomic update
                            ! !$acc atomic
                            ! rhs_vf(momxb)%sf(j,k,l) = rhs_vf(momxb)%sf(j,k,l) - &
                            !      (0.5_wp* (rho_R * (vel_R(1))**2.0 + &
                            !      pres_R)*(1._wp/dx(j)) + &
                            !      0.5_wp*cfl * (rho_R*vel_R(1))*(1._wp/dx(j)))

                            @:ATOMIC_ADD(rhs_vf(momxb+1)%sf(j, k, l),(-0.5_wp * rho_R * vel_R(1)*vel_R(2)*(1._wp/dx(j)) - &
                                0.5_wp*cfl * (rho_R*vel_R(2))*(1._wp/dx(j)))*dt)
                            ! !$omp atomic update
                            ! !$acc atomic
                            ! rhs_vf(momxb+1)%sf(j, k, l) =  rhs_vf(momxb+1)%sf(j, k, l) - &
                            !     (0.5_wp * rho_R * vel_R(1)*vel_R(2)*(1._wp/dx(j)) + &
                            !     0.5_wp*cfl * (rho_R*vel_R(2))*(1._wp/dx(j)))

                            @:ATOMIC_ADD(rhs_vf(E_idx)%sf(j, k, l),(-0.5_wp * (vel_R(1) * (E_R + &
                                pres_R) )*(1._wp/dx(j)) - &
                                0.5_wp*cfl * (E_R)*(1._wp/dx(j)))*dt)
                            ! !$omp atomic update
                            ! !$acc atomic
                            !  rhs_vf(E_idx)%sf(j, k, l) = rhs_vf(E_idx)%sf(j, k, l) -  &
                            !     (0.5_wp * (vel_R(1) * (E_R + &
                            !     pres_R) )*(1._wp/dx(j)) + &
                            !     0.5_wp*cfl * (E_R)*(1._wp/dx(j)))
                        end do
                    end do
                end do
            else
                !$omp target teams loop collapse(3) private(rho_L,gamma_L,pi_inf_L,mu_L,a_L,vel_L,vel_R,pres_L,alpha_L,alpha_R,alpha_rho_L,cfl,dvel,F_L,E_L,mu_R,rho_sf_small,alpha_rho_R,vflux_L_arr,vflux_R_arr,dvel_small)
                !$acc parallel loop collapse(3) gang vector default(present) private(rho_L,gamma_L,pi_inf_L,mu_L,a_L,vel_L,vel_R,pres_L,alpha_L,alpha_R,alpha_rho_L,cfl,dvel,F_L,E_L,mu_R,rho_sf_small,alpha_rho_R,vflux_L_arr,vflux_R_arr,dvel_small)
                do l = 0, p
                    do k = 0, n
                        do j = -buff_size+5, m+buff_size-5
                            
                            dvel = 0._wp
                            vflux_L_arr = 0._wp
                            vflux_R_arr = 0._wp

                            !DIR$ unroll 6
                            !$acc loop seq 
                            do q = -2, 3
                                dvel_small = 0._wp
                                !x-direction contributions
                                !$acc loop seq 
                                do i = -1, 1
                                    rho_L = 0._wp
                                    !$acc loop seq 
                                    do r = 1, num_fluids
                                        rho_L = rho_L + q_prim_vf(r)%sf(j+i+q,k,l)
                                    end do
                                    rho_sf_small(i) = rho_L
                                end do

                                dvel_small(1) = (1/(2._wp*dx(j))) * ( &
                                    q_prim_vf(momxb)%sf(j+1+q,k,l)/rho_sf_small(1) - &
                                    q_prim_vf(momxb)%sf(j-1+q,k,l)/rho_sf_small(-1))
                                dvel_small(2) = (1/(2._wp*dx(j))) * ( &
                                    q_prim_vf(momxb+1)%sf(j+1+q,k,l)/rho_sf_small(1) - &
                                    q_prim_vf(momxb+1)%sf(j-1+q,k,l)/rho_sf_small(-1))
                                dvel_small(3) = (1/(2._wp*dx(j))) * ( &
                                    q_prim_vf(momxb+2)%sf(j+1+q,k,l)/rho_sf_small(1) - &
                                    q_prim_vf(momxb+2)%sf(j-1+q,k,l)/rho_sf_small(-1))

                                if (q == 0) dvel(:,1) = dvel_small
                                if (q > -2 .and. viscous) then
                                    vflux_L_arr(1) = vflux_L_arr(1) + coeff_L(q)*(dvel_small(2))
                                    vflux_L_arr(2) = vflux_L_arr(2) + coeff_L(q)*(dvel_small(3))
                                    vflux_L_arr(3) = vflux_L_arr(3) + coeff_L(q)*(4._wp*dvel_small(1))/3._wp
                                end if
                                if (q < 3 .and. viscous) then
                                    vflux_R_arr(1) = vflux_R_arr(1) + coeff_R(q)*(dvel_small(2))
                                    vflux_R_arr(2) = vflux_R_arr(2) + coeff_R(q)*(dvel_small(3))
                                    vflux_R_arr(3) = vflux_R_arr(3) + coeff_R(q)*(4._wp*dvel_small(1))/3._wp
                                end if

                                !y-direction contributions
                                !$acc loop seq 
                                do i = -1, 1
                                    rho_L = 0._wp
                                    !$acc loop seq 
                                    do r = 1, num_fluids
                                        rho_L = rho_L + q_prim_vf(r)%sf(j+q,k+i,l)
                                    end do
                                    rho_sf_small(i) = rho_L
                                end do

                                dvel_small(1) = (1/(2._wp*dy(k))) * ( &
                                    q_prim_vf(momxb)%sf(j+q,k+1,l)/rho_sf_small(1) - &
                                    q_prim_vf(momxb)%sf(j+q,k-1,l)/rho_sf_small(-1))
                                dvel_small(2) = (1/(2._wp*dy(k))) * ( &
                                    q_prim_vf(momxb+1)%sf(j+q,k+1,l)/rho_sf_small(1) - &
                                    q_prim_vf(momxb+1)%sf(j+q,k-1,l)/rho_sf_small(-1))
                                if (q == 0) dvel_small(3) = (1/(2._wp*dy(k))) * ( &
                                    q_prim_vf(momxb+2)%sf(j+q,k+1,l)/rho_sf_small(1) - &
                                    q_prim_vf(momxb+2)%sf(j+q,k-1,l)/rho_sf_small(-1))
                                if (q == 0) dvel(:,2) = dvel_small

                                if (q > -2 .and. viscous) then
                                    vflux_L_arr(1) = vflux_L_arr(1) + coeff_L(q)*(dvel_small(1))
                                    vflux_L_arr(3) = vflux_L_arr(3) + coeff_L(q)*(-2._wp*dvel_small(2))/3._wp
                                end if
                                if (q < 3 .and. viscous) then
                                    vflux_R_arr(1) = vflux_R_arr(1) + coeff_R(q)*(dvel_small(1))
                                    vflux_R_arr(3) = vflux_R_arr(3) + coeff_R(q)*(-2._wp*dvel_small(2))/3._wp
                                end if

                                !z-direction contributions
                                !$acc loop seq 
                                do i = -1, 1
                                    rho_L = 0._wp
                                    !$acc loop seq 
                                    do r = 1, num_fluids
                                        rho_L = rho_L + q_prim_vf(r)%sf(j+q,k,l+i)
                                    end do
                                    rho_sf_small(i) = rho_L
                                end do

                                dvel_small(1) = (1/(2._wp*dz(l))) * ( &
                                    q_prim_vf(momxb)%sf(j+q,k,l+1)/rho_sf_small(1) - &
                                    q_prim_vf(momxb)%sf(j+q,k,l-1)/rho_sf_small(-1))
                                if (q == 0) dvel_small(2) = (1/(2._wp*dz(l))) * ( &
                                    q_prim_vf(momxb+1)%sf(j+q,k,l+1)/rho_sf_small(1) - &
                                    q_prim_vf(momxb+1)%sf(j+q,k,l-1)/rho_sf_small(-1))
                                dvel_small(3) = (1/(2._wp*dz(l))) * ( &
                                    q_prim_vf(momxb+2)%sf(j+q,k,l+1)/rho_sf_small(1) - &
                                    q_prim_vf(momxb+2)%sf(j+q,k,l-1)/rho_sf_small(-1))
                                if (q == 0) dvel(:,3) = dvel_small
                                if (q > -2 .and. viscous) then
                                    vflux_L_arr(2) = vflux_L_arr(2) + coeff_L(q)*(dvel_small(1))
                                    vflux_L_arr(3) = vflux_L_arr(3) + coeff_L(q)*(-2._wp*dvel_small(3))/3._wp
                                end if
                                if (q < 3 .and. viscous) then
                                    vflux_R_arr(2) = vflux_R_arr(2) + coeff_R(q)*(dvel_small(1))
                                    vflux_R_arr(3) = vflux_R_arr(3) + coeff_R(q)*(-2._wp*dvel_small(3))/3._wp
                                end if

                                if (q == 0) then
                                    jac_rhs(j,k,l) = alf_igr * (2._wp*(dvel(1,2)*dvel(2,1) &
                                        + dvel(1,3)*dvel(3,1) &
                                        + dvel(2,3)*dvel(3,2)) &
                                        + dvel(1,1)**2_wp + dvel(2,2)**2_wp &
                                        + dvel(3,3)**2_wp &
                                        + (dvel(1,1) + dvel(2,2)+ dvel(3,3))**2_wp)
                                end if
                            end do

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
                            end do 

                             if(num_fluids > 1) then 
                                !$acc loop seq 
                                do i = 1, num_fluids - 1
                                    alpha_L(i) =  (1._wp/60._wp) * (-3._wp * q_prim_vf(E_idx+i)%sf(j-1, k, l) + &
                                                    27._wp * q_prim_vf(E_idx+i)%sf(j, k, l) + &
                                                    47._wp * q_prim_vf(E_idx+i)%sf(j+1, k, l) -   &
                                                    13._wp * q_prim_vf(E_idx+i)%sf(j+2, k, l) + &
                                                    2._wp * q_prim_vf(E_idx+i)%sf(j+3, k, l))
                                    alpha_R(i) = (1._wp/60._wp) * (-3._wp * q_prim_vf(E_idx+i)%sf(j+2, k, l) + &
                                                            27._wp * q_prim_vf(E_idx+i)%sf(j+1, k, l) + &
                                                            47._wp * q_prim_vf(E_idx+i)%sf(j, k, l) -   &
                                                            13._wp * q_prim_vf(E_idx+i)%sf(j-1, k, l) + &
                                                            2._wp * q_prim_vf(E_idx+i)%sf(j-2, k, l))
                                end do

                                alpha_L(num_fluids) = 1._wp - sum(alpha_L(1:num_fluids - 1))
                                alpha_R(num_fluids) = 1._wp - sum(alpha_R(1:num_fluids - 1))    
                            else 
                                alpha_L(1) = 1._wp 
                                alpha_R(1) = 1._wp
                            end if

                            rho_L = sum(alpha_rho_L)
                            gamma_L = sum(alpha_L*gammas)
                            pi_inf_L = sum(alpha_L*pi_infs)

                            rho_R = sum(alpha_rho_R)
                            gamma_R = sum(alpha_R*gammas)
                            pi_inf_R = sum(alpha_R*pi_infs)

                            !$acc loop seq 
                            do i = 1, num_dims
                                vel_L(i) = (1._wp/60._wp) * (-3._wp * q_prim_vf(momxb+i-1)%sf(j-1,k,l) + &
                                                27._wp * q_prim_vf(momxb+i-1)%sf(j,k,l) + &
                                                47._wp * q_prim_vf(momxb+i-1)%sf(j+1,k,l) -   &
                                                13._wp * q_prim_vf(momxb+i-1)%sf(j+2,k,l) + &
                                                2._wp * q_prim_vf(momxb+i-1)%sf(j+3,k,l)) / rho_L
                            end do 

                            !$acc loop seq 
                            do i = 1, num_dims
                                vel_R(i) = (1._wp/60._wp) * (-3._wp * q_prim_vf(momxb+i-1)%sf(j+2,k,l) + &
                                                27._wp * q_prim_vf(momxb+i-1)%sf(j+1,k,l) + &
                                                47._wp * q_prim_vf(momxb+i-1)%sf(j,k,l) -   &
                                                13._wp * q_prim_vf(momxb+i-1)%sf(j-1,k,l) + &
                                                2._wp * q_prim_vf(momxb+i-1)%sf(j-2,k,l)) / rho_R
                            end do
                            if(viscous) then

                                if(num_fluids > 1) then 
                                    mu_L = 0._wp
                                    mu_R = 0._wp
                                    !$acc loop seq 
                                    do i = 1, num_fluids - 1
                                        mu_L = (1._wp/60._wp) * (-3._wp * q_prim_vf(E_idx+i)%sf(j-1, k, l) + &
                                                27._wp * q_prim_vf(E_idx+i)%sf(j, k, l) + &
                                                47._wp * q_prim_vf(E_idx+i)%sf(j+1, k, l) -   &
                                                13._wp * q_prim_vf(E_idx+i)%sf(j+2, k, l) + &
                                                2._wp * q_prim_vf(E_idx+i)%sf(j+3, k, l)) / Res(1, i) + mu_L

                                        mu_R = (1._wp/60._wp) * (-3._wp * q_prim_vf(E_idx+i)%sf(j+2, k, l) + &
                                                27._wp * q_prim_vf(E_idx+i)%sf(j+1, k, l) + &
                                                47._wp * q_prim_vf(E_idx+i)%sf(j, k, l) -   &
                                                13._wp * q_prim_vf(E_idx+i)%sf(j-1, k, l) + &
                                                2._wp * q_prim_vf(E_idx+i)%sf(j-2, k, l)) / Res(1, i) + mu_R
                                    end do

                                    mu_L = 1._wp / Res(1, num_fluids) + mu_L
                                    mu_R = 1._wp / Res(1, num_fluids) + mu_R

                                    !$acc loop seq 
                                    do i = 1, num_fluids - 1
                                        mu_L = -(1._wp/60._wp) * (-3._wp * q_prim_vf(E_idx+i)%sf(j-1, k, l) + &
                                                27._wp * q_prim_vf(E_idx+i)%sf(j, k, l) + &
                                                47._wp * q_prim_vf(E_idx+i)%sf(j+1, k, l) -   &
                                                13._wp * q_prim_vf(E_idx+i)%sf(j+2, k, l) + &
                                                2._wp * q_prim_vf(E_idx+i)%sf(j+3, k, l)) / Res(1, num_fluids) + mu_L

                                        mu_R = -(1._wp/60._wp) * (-3._wp * q_prim_vf(E_idx+i)%sf(j+2, k, l) + &
                                                27._wp * q_prim_vf(E_idx+i)%sf(j+1, k, l) + &
                                                47._wp * q_prim_vf(E_idx+i)%sf(j, k, l) -   &
                                                13._wp * q_prim_vf(E_idx+i)%sf(j-1, k, l) + &
                                                2._wp * q_prim_vf(E_idx+i)%sf(j-2, k, l)) / Res(1, num_fluids) + mu_R
                                    end do
                                else
                                    mu_L = 1._wp / Res(1, 1) 
                                    mu_R = 1._wp / Res(1, 1) 
                                end if


                                @:ATOMIC_ADD(rhs_vf(momxb+1)%sf(j+1,k,l),(- 0.5_wp*mu_L*vflux_L_arr(1)*(1._wp/dx(j+1))*dt))
                                !!$omp atomic update
                                !!$acc atomic
                                !rhs_vf(momxb+1)%sf(j+1,k,l) = rhs_vf(momxb+1)%sf(j+1,k,l) - 0.5_wp*mu_L*vflux_L_arr(1)*(1._wp/dx(j+1))
                                @:ATOMIC_ADD(rhs_vf(E_idx)%sf(j+1,k,l),(- 0.5_wp*mu_L*vflux_L_arr(1)*vel_L(2)*(1._wp/dx(j+1))*dt))
                                !!$omp atomic update
                                !!$acc atomic
                                !rhs_vf(E_idx)%sf(j+1,k,l) = rhs_vf(E_idx)%sf(j+1,k,l) - 0.5_wp*mu_L*vflux_L_arr(1)*vel_L(2)*(1._wp/dx(j+1))
                        
                                @:ATOMIC_ADD(rhs_vf(momxb+1)%sf(j,k,l),(0.5_wp*mu_L*vflux_L_arr(1)*(1._wp/dx(j))*dt))
                                !!$omp atomic update
                                !!$acc atomic
                                !rhs_vf(momxb+1)%sf(j,k,l) = rhs_vf(momxb+1)%sf(j,k,l) + 0.5_wp*mu_L*vflux_L_arr(1)*(1._wp/dx(j))
                                @:ATOMIC_ADD(rhs_vf(E_idx)%sf(j,k,l),(0.5_wp*mu_L*vflux_L_arr(1)*vel_L(2)*(1._wp/dx(j))*dt))
                                !!$omp atomic update
                                !!$acc atomic
                                !rhs_vf(E_idx)%sf(j,k,l) = rhs_vf(E_idx)%sf(j,k,l) + 0.5_wp*mu_L*vflux_L_arr(1)*vel_L(2)*(1._wp/dx(j))                    
                             
                                @:ATOMIC_ADD(rhs_vf(momxb+1)%sf(j+1,k,l),(- 0.5_wp*mu_R*vflux_R_arr(1)*(1._wp/dx(j+1))*dt))
                                !!$omp atomic update
                                !!$acc atomic
                                !rhs_vf(momxb+1)%sf(j+1,k,l) = rhs_vf(momxb+1)%sf(j+1,k,l) - 0.5_wp*mu_R*vflux_R_arr(1)*(1._wp/dx(j+1))
                                @:ATOMIC_ADD(rhs_vf(E_idx)%sf(j+1,k,l),(- 0.5_wp*mu_R*vflux_R_arr(1)*vel_R(2)*(1._wp/dx(j+1))*dt))
                                !!$omp atomic update
                                !!$acc atomic
                                !rhs_vf(E_idx)%sf(j+1,k,l) = rhs_vf(E_idx)%sf(j+1,k,l) - 0.5_wp*mu_R*vflux_R_arr(1)*vel_R(2)*(1._wp/dx(j+1))
                             
                                @:ATOMIC_ADD(rhs_vf(momxb+1)%sf(j,k,l),(0.5_wp*mu_R*vflux_R_arr(1)*(1._wp/dx(j))*dt))
                                !!$omp atomic update
                                !!$acc atomic
                                !rhs_vf(momxb+1)%sf(j,k,l) = rhs_vf(momxb+1)%sf(j,k,l) + 0.5_wp*mu_R*vflux_R_arr(1)*(1._wp/dx(j))
                                @:ATOMIC_ADD(rhs_vf(E_idx)%sf(j,k,l),(0.5_wp*mu_R*vflux_R_arr(1)*vel_R(2)*(1._wp/dx(j))*dt))
                                !!$omp atomic update
                                !!$acc atomic
                                !rhs_vf(E_idx)%sf(j,k,l) = rhs_vf(E_idx)%sf(j,k,l) + 0.5_wp*mu_R*vflux_R_arr(1)*vel_R(2)*(1._wp/dx(j))


                                @:ATOMIC_ADD(rhs_vf(momxb+2)%sf(j+1,k,l),(- 0.5_wp*mu_L*vflux_L_arr(2)*(1._wp/dx(j+1))*dt))
                                !!$omp atomic update
                                !!$acc atomic
                                !rhs_vf(momxb+2)%sf(j+1,k,l) = rhs_vf(momxb+2)%sf(j+1,k,l) - 0.5_wp*mu_L*vflux_L_arr(2)*(1._wp/dx(j+1))
                                @:ATOMIC_ADD(rhs_vf(E_idx)%sf(j+1,k,l),(- 0.5_wp*mu_L*vflux_L_arr(2)*vel_L(3)*(1._wp/dx(j+1))*dt))
                                !!$omp atomic update
                                !!$acc atomic
                                !rhs_vf(E_idx)%sf(j+1,k,l) = rhs_vf(E_idx)%sf(j+1,k,l) - 0.5_wp*mu_L*vflux_L_arr(2)*vel_L(3)*(1._wp/dx(j+1))
                            
                                @:ATOMIC_ADD(rhs_vf(momxb+2)%sf(j,k,l),(0.5_wp*mu_L*vflux_L_arr(2)*(1._wp/dx(j))*dt))
                                !!$omp atomic update
                                !!$acc atomic
                                !rhs_vf(momxb+2)%sf(j,k,l) = rhs_vf(momxb+2)%sf(j,k,l) + 0.5_wp*mu_L*vflux_L_arr(2)*(1._wp/dx(j))
                                @:ATOMIC_ADD(rhs_vf(E_idx)%sf(j,k,l),(0.5_wp*mu_L*vflux_L_arr(2)*vel_L(3)*(1._wp/dx(j))*dt))
                                !!$omp atomic update
                                !!$acc atomic
                                !rhs_vf(E_idx)%sf(j,k,l) = rhs_vf(E_idx)%sf(j,k,l) + 0.5_wp*mu_L*vflux_L_arr(2)*vel_L(3)*(1._wp/dx(j))
                            
                                @:ATOMIC_ADD(rhs_vf(momxb+2)%sf(j+1,k,l),(- 0.5_wp*mu_R*vflux_R_arr(2)*(1._wp/dx(j+1))*dt))
                                !!$omp atomic update
                                !!$acc atomic
                                !rhs_vf(momxb+2)%sf(j+1,k,l) = rhs_vf(momxb+2)%sf(j+1,k,l) - 0.5_wp*mu_R*vflux_R_arr(2)*(1._wp/dx(j+1))
                                @:ATOMIC_ADD(rhs_vf(E_idx)%sf(j+1,k,l),(- 0.5_wp*mu_R*vflux_R_arr(2)*vel_R(3)*(1._wp/dx(j+1))*dt))
                                !!$omp atomic update
                                !!$acc atomic
                                !rhs_vf(E_idx)%sf(j+1,k,l) = rhs_vf(E_idx)%sf(j+1,k,l) - 0.5_wp*mu_R*vflux_R_arr(2)*vel_R(3)*(1._wp/dx(j+1))
                            
                                @:ATOMIC_ADD(rhs_vf(momxb+2)%sf(j,k,l),(0.5_wp*mu_R*vflux_R_arr(2)*(1._wp/dx(j))*dt))
                                !!$omp atomic update
                                !!$acc atomic
                                !rhs_vf(momxb+2)%sf(j,k,l) = rhs_vf(momxb+2)%sf(j,k,l) + 0.5_wp*mu_R*vflux_R_arr(2)*(1._wp/dx(j))
                                @:ATOMIC_ADD(rhs_vf(E_idx)%sf(j,k,l),(0.5_wp*mu_R*vflux_R_arr(2)*vel_R(3)*(1._wp/dx(j))*dt))
                                !!$omp atomic update
                                !!$acc atomic
                                !rhs_vf(E_idx)%sf(j,k,l) = rhs_vf(E_idx)%sf(j,k,l) + 0.5_wp*mu_R*vflux_R_arr(2)*vel_R(3)*(1._wp/dx(j))


                                @:ATOMIC_ADD(rhs_vf(momxb)%sf(j+1,k,l),(- 0.5_wp*mu_L*vflux_L_arr(3)*(1._wp/dx(j+1))*dt))
                                !!$omp atomic update
                                !!$acc atomic
                                !rhs_vf(momxb)%sf(j+1,k,l) = rhs_vf(momxb)%sf(j+1,k,l) - 0.5_wp*mu_L*vflux_L_arr(3)*(1._wp/dx(j+1))
                                @:ATOMIC_ADD(rhs_vf(E_idx)%sf(j+1,k,l),(- 0.5_wp*mu_L*vflux_L_arr(3)*vel_L(1)*(1._wp/dx(j+1))*dt))
                                !!$omp atomic update
                                !!$acc atomic
                                !rhs_vf(E_idx)%sf(j+1,k,l) = rhs_vf(E_idx)%sf(j+1,k,l) - 0.5_wp*mu_L*vflux_L_arr(3)*vel_L(1)*(1._wp/dx(j+1)) 

                                @:ATOMIC_ADD(rhs_vf(momxb)%sf(j,k,l),(0.5_wp*mu_L*vflux_L_arr(3)*(1._wp/dx(j))*dt))
                                !!$omp atomic update
                                !!$acc atomic
                                !rhs_vf(momxb)%sf(j,k,l) = rhs_vf(momxb)%sf(j,k,l) + 0.5_wp*mu_L*vflux_L_arr(3)*(1._wp/dx(j))
                                @:ATOMIC_ADD(rhs_vf(E_idx)%sf(j,k,l),(0.5_wp*mu_L*vflux_L_arr(3)*vel_L(1)*(1._wp/dx(j))*dt))
                                !!$omp atomic update
                                !!$acc atomic
                                !rhs_vf(E_idx)%sf(j,k,l) = rhs_vf(E_idx)%sf(j,k,l) + 0.5_wp*mu_L*vflux_L_arr(3)*vel_L(1)*(1._wp/dx(j))
                            
                                @:ATOMIC_ADD(rhs_vf(momxb)%sf(j+1,k,l),(- 0.5_wp*mu_R*vflux_R_arr(3)*(1._wp/dx(j+1))*dt))
                                !!$omp atomic update
                                !!$acc atomic
                                !rhs_vf(momxb)%sf(j+1,k,l) = rhs_vf(momxb)%sf(j+1,k,l) - 0.5_wp*mu_R*vflux_R_arr(3)*(1._wp/dx(j+1))
                                @:ATOMIC_ADD(rhs_vf(E_idx)%sf(j+1,k,l),(- 0.5_wp*mu_R*vflux_R_arr(3)*vel_R(1)*(1._wp/dx(j+1))*dt))
                                !!$omp atomic update
                                !!$acc atomic
                                !rhs_vf(E_idx)%sf(j+1,k,l) = rhs_vf(E_idx)%sf(j+1,k,l) - 0.5_wp*mu_R*vflux_R_arr(3)*vel_R(1)*(1._wp/dx(j+1))
                            
                                @:ATOMIC_ADD(rhs_vf(momxb)%sf(j,k,l),(0.5_wp*mu_R*vflux_R_arr(3)*(1._wp/dx(j))*dt))
                                !!$omp atomic update
                                !!$acc atomic
                                !rhs_vf(momxb)%sf(j,k,l) = rhs_vf(momxb)%sf(j,k,l) + 0.5_wp*mu_R*vflux_R_arr(3)*(1._wp/dx(j))
                                @:ATOMIC_ADD(rhs_vf(E_idx)%sf(j,k,l),(0.5_wp*mu_R*vflux_R_arr(3)*vel_R(1)*(1._wp/dx(j))*dt))
                                !!$omp atomic update
                                !!$acc atomic
                                !rhs_vf(E_idx)%sf(j,k,l) = rhs_vf(E_idx)%sf(j,k,l) + 0.5_wp*mu_R*vflux_R_arr(3)*vel_R(1)*(1._wp/dx(j))  
                            
                                                            
                            endif
                            
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

                            pres_L = (E_L - pi_inf_L - 0.5_wp*rho_L*(vel_L(1)**2._wp + vel_L(2)**2._wp + vel_L(3)**2._wp))/gamma_L                    

                            pres_R = (E_R - pi_inf_R - 0.5_wp*rho_R*(vel_R(1)**2._wp + vel_R(2)**2._wp + vel_R(3)**2._wp))/gamma_R 

                            a_L = sqrt((pres_L*(1._wp/gamma_L + 1._wp) + pi_inf_L / gamma_L) / rho_L)
                            a_R = sqrt((pres_R*(1._wp/gamma_R + 1._wp) + pi_inf_R / gamma_R) / rho_R)

                            cfl = (max(sqrt(vel_L(1)**2._wp + vel_L(2)**2._wp + vel_L(3)**2._wp),sqrt(vel_R(1)**2._wp + vel_R(2)**2._wp + vel_R(3)**2._wp)) + max(a_L,a_R))

                            !$acc loop seq
                            do i = 1, num_fluids
                                @:ATOMIC_ADD(rhs_vf(i)%sf(j+1,k,l),(0.5_wp * (alpha_rho_L(i) * &
                                    vel_L(1))*(1._wp/dx(j+1)) - &
                                    0.5_wp*cfl * (alpha_rho_L(i))*(1._wp/dx(j+1)))*dt)
                                ! !$omp atomic update
                                ! !$acc atomic
                                ! rhs_vf(i)%sf(j+1,k,l) = rhs_vf(i)%sf(j+1,k,l) + &
                                !     (0.5_wp * (alpha_rho_L(i) * &
                                !     vel_L(1))*(1._wp/dx(j+1)) - &
                                !     0.5_wp*cfl * (alpha_rho_L(i))*(1._wp/dx(j+1)))

                                @:ATOMIC_ADD(rhs_vf(i)%sf(j,k,l),(-0.5_wp * (alpha_rho_L(i) * &
                                    vel_L(1))*(1._wp/dx(j)) + &
                                    0.5_wp*cfl * (alpha_rho_L(i))*(1._wp/dx(j)))*dt)
                                ! !$omp atomic update
                                ! !$acc atomic
                                ! rhs_vf(i)%sf(j,k,l) = rhs_vf(i)%sf(j,k,l) - &
                                !     (0.5_wp * (alpha_rho_L(i) * &
                                !     vel_L(1))*(1._wp/dx(j)) - &
                                !     0.5_wp*cfl * (alpha_rho_L(i))*(1._wp/dx(j)))
                            end do

                            if(num_fluids > 1) then 
                                !$acc loop seq
                                do i = 1, num_fluids - 1
                                    @:ATOMIC_ADD(rhs_vf(advxb+i-1)%sf(j+1,k,l),(0.5_wp * (alpha_L(i) * &
                                        vel_L(1))*(1._wp/dx(j+1)) - &
                                        0.5_wp*cfl*(alpha_L(i))*(1._wp/dx(j+1)))*dt)
                                    ! !$omp atomic update
                                    ! !$acc atomic
                                    ! rhs_vf(advxb+i-1)%sf(j+1,k,l) = rhs_vf(advxb+i-1)%sf(j+1,k,l) + &
                                    !     (0.5_wp * (alpha_L(i) * &
                                    !     vel_L(1))*(1._wp/dx(j+1)) - &
                                    !     0.5_wp*cfl*(alpha_L(i))*(1._wp/dx(j+1)))

                                    @:ATOMIC_ADD(rhs_vf(advxb+i-1)%sf(j+1,k,l),(-0.5_wp * q_prim_vf(advxb+i-1)%sf(j+1,k,l)  * vel_L(1)*(1._wp/dx(j+1)))*dt)
                                    ! !$omp atomic update
                                    ! !$acc atomic
                                    ! rhs_vf(advxb+i-1)%sf(j+1,k,l) = rhs_vf(advxb+i-1)%sf(j+1,k,l) &
                                    ! - (0.5_wp * q_prim_vf(advxb+i-1)%sf(j+1,k,l)  * vel_L(1)*(1._wp/dx(j+1)))

                                    @:ATOMIC_ADD(rhs_vf(advxb+i-1)%sf(j,k,l),(-0.5_wp * (alpha_L(i) * &
                                        vel_L(1))*(1._wp/dx(j)) + &
                                        0.5_wp*cfl*(alpha_L(i))*(1._wp/dx(j)))*dt)
                                    ! !$omp atomic update
                                    ! !$acc atomic
                                    ! rhs_vf(advxb+i-1)%sf(j,k,l) = rhs_vf(advxb+i-1)%sf(j,k,l) - &
                                    !     (0.5_wp * (alpha_L(i) * &
                                    !     vel_L(1))*(1._wp/dx(j)) - &
                                    !     0.5_wp*cfl*(alpha_L(i))*(1._wp/dx(j)))

                                    @:ATOMIC_ADD(rhs_vf(advxb+i-1)%sf(j,k,l),(0.5_wp * q_prim_vf(advxb+i-1)%sf(j,k,l)  * vel_L(1)*(1._wp/dx(j)))*dt)
                                    ! !$omp atomic update
                                    ! !$acc atomic
                                    ! rhs_vf(advxb+i-1)%sf(j,k,l) = rhs_vf(advxb+i-1)%sf(j,k,l) &
                                    ! + (0.5_wp * q_prim_vf(advxb+i-1)%sf(j,k,l)  * vel_L(1)*(1._wp/dx(j)))
                                end do
                            end if

                            @:ATOMIC_ADD(rhs_vf(momxb)%sf(j+1,k,l),(0.5_wp* (rho_L * (vel_L(1))**2.0 + &
                                 pres_L)*(1._wp/dx(j+1)) - &
                                 0.5_wp*cfl * (rho_L*vel_L(1))*(1._wp/dx(j+1)))*dt)
                            ! !$omp atomic update
                            ! !$acc atomic
                            ! rhs_vf(momxb)%sf(j+1,k,l) = rhs_vf(momxb)%sf(j+1,k,l) + &
                            !      (0.5_wp* (rho_L * (vel_L(1))**2.0 + &
                            !      pres_L)*(1._wp/dx(j+1)) - &
                            !      0.5_wp*cfl * (rho_L*vel_L(1))*(1._wp/dx(j+1)))

                            @:ATOMIC_ADD(rhs_vf(momxb+1)%sf(j+1, k, l),(0.5_wp * rho_L * vel_L(1)*vel_L(2)*(1._wp/dx(j+1)) - &
                                0.5_wp*cfl * (rho_L*vel_L(2))*(1._wp/dx(j+1)))*dt)
                            ! !$omp atomic update
                            ! !$acc atomic
                            ! rhs_vf(momxb+1)%sf(j+1, k, l) =  rhs_vf(momxb+1)%sf(j+1, k, l) + &
                            !     (0.5_wp * rho_L * vel_L(1)*vel_L(2)*(1._wp/dx(j+1)) - &
                            !     0.5_wp*cfl * (rho_L*vel_L(2))*(1._wp/dx(j+1)))

                            @:ATOMIC_ADD(rhs_vf(momxb+2)%sf(j+1, k, l),(0.5_wp * rho_L * vel_L(1)*vel_L(3)*(1._wp/dx(j+1)) - &
                                0.5_wp*cfl * (rho_L*vel_L(3))*(1._wp/dx(j+1)))*dt)
                            ! !$omp atomic update
                            ! !$acc atomic
                            ! rhs_vf(momxb+2)%sf(j+1, k, l) =  rhs_vf(momxb+2)%sf(j+1, k, l) + &
                            !     (0.5_wp * rho_L * vel_L(1)*vel_L(3)*(1._wp/dx(j+1)) - &
                            !     0.5_wp*cfl * (rho_L*vel_L(3))*(1._wp/dx(j+1)))

                            @:ATOMIC_ADD(rhs_vf(E_idx)%sf(j+1, k, l),(0.5_wp * (vel_L(1) * (E_L + &
                                pres_L) )*(1._wp/dx(j+1)) - &
                                0.5_wp*cfl * (E_L)*(1._wp/dx(j+1)))*dt)
                            ! !$omp atomic update
                            ! !$acc atomic
                            !  rhs_vf(E_idx)%sf(j+1, k, l) = rhs_vf(E_idx)%sf(j+1, k, l) +  &
                            !     (0.5_wp * (vel_L(1) * (E_L + &
                            !     pres_L) )*(1._wp/dx(j+1)) - &
                            !     0.5_wp*cfl * (E_L)*(1._wp/dx(j+1)))

                            @:ATOMIC_ADD(rhs_vf(momxb)%sf(j,k,l),(-0.5_wp* (rho_L * (vel_L(1))**2.0 + &
                                 pres_L)*(1._wp/dx(j)) + &
                                 0.5_wp*cfl * (rho_L*vel_L(1))*(1._wp/dx(j)))*dt)
                            ! !$omp atomic update
                            ! !$acc atomic
                            ! rhs_vf(momxb)%sf(j,k,l) = rhs_vf(momxb)%sf(j,k,l) - &
                            !      (0.5_wp* (rho_L * (vel_L(1))**2.0 + &
                            !      pres_L)*(1._wp/dx(j)) - &
                            !      0.5_wp*cfl * (rho_L*vel_L(1))*(1._wp/dx(j)))

                            @:ATOMIC_ADD(rhs_vf(momxb+1)%sf(j, k, l),(-0.5_wp * rho_L * vel_L(1)*vel_L(2)*(1._wp/dx(j)) + &
                                0.5_wp*cfl * (rho_L*vel_L(2))*(1._wp/dx(j)))*dt)
                            ! !$omp atomic update
                            ! !$acc atomic
                            ! rhs_vf(momxb+1)%sf(j, k, l) =  rhs_vf(momxb+1)%sf(j, k, l) - &
                            !     (0.5_wp * rho_L * vel_L(1)*vel_L(2)*(1._wp/dx(j)) - &
                            !     0.5_wp*cfl * (rho_L*vel_L(2))*(1._wp/dx(j)))

                            @:ATOMIC_ADD(rhs_vf(momxb+2)%sf(j, k, l),(-0.5_wp * rho_L * vel_L(1)*vel_L(3)*(1._wp/dx(j)) + &
                                0.5_wp*cfl * (rho_L*vel_L(3))*(1._wp/dx(j)))*dt)
                            ! !$omp atomic update
                            ! !$acc atomic
                            ! rhs_vf(momxb+2)%sf(j, k, l) =  rhs_vf(momxb+2)%sf(j, k, l) - &
                            !     (0.5_wp * rho_L * vel_L(1)*vel_L(3)*(1._wp/dx(j)) - &
                            !     0.5_wp*cfl * (rho_L*vel_L(3))*(1._wp/dx(j)))

                            @:ATOMIC_ADD(rhs_vf(E_idx)%sf(j, k, l),(-0.5_wp * (vel_L(1) * (E_L + &
                                pres_L) )*(1._wp/dx(j)) + &
                                0.5_wp*cfl * (E_L)*(1._wp/dx(j)))*dt)
                            ! !$omp atomic update
                            ! !$acc atomic
                            !  rhs_vf(E_idx)%sf(j, k, l) = rhs_vf(E_idx)%sf(j, k, l) -  &
                            !     (0.5_wp * (vel_L(1) * (E_L + &
                            !     pres_L) )*(1._wp/dx(j)) - &
                            !     0.5_wp*cfl * (E_L)*(1._wp/dx(j)))

                            !$acc loop seq
                            do i = 1, num_fluids
                                @:ATOMIC_ADD(rhs_vf(i)%sf(j+1,k,l),(0.5_wp * (alpha_rho_R(i) * &
                                    vel_R(1))*(1._wp/dx(j+1)) + &
                                    0.5_wp*cfl * (alpha_rho_R(i))*(1._wp/dx(j+1)))*dt)
                                ! !$omp atomic update
                                ! !$acc atomic
                                ! rhs_vf(i)%sf(j+1,k,l) = rhs_vf(i)%sf(j+1,k,l) + &
                                !     (0.5_wp * (alpha_rho_R(i) * &
                                !     vel_R(1))*(1._wp/dx(j+1)) + &
                                !     0.5_wp*cfl * (alpha_rho_R(i))*(1._wp/dx(j+1)))

                                @:ATOMIC_ADD(rhs_vf(i)%sf(j,k,l),(-0.5_wp * (alpha_rho_R(i) * &
                                    vel_R(1))*(1._wp/dx(j)) - &
                                    0.5_wp*cfl * (alpha_rho_R(i))*(1._wp/dx(j)))*dt)
                                ! !$omp atomic update
                                ! !$acc atomic
                                ! rhs_vf(i)%sf(j,k,l) = rhs_vf(i)%sf(j,k,l) - &
                                !     (0.5_wp * (alpha_rho_R(i) * &
                                !     vel_R(1))*(1._wp/dx(j)) + &
                                !     0.5_wp*cfl * (alpha_rho_R(i))*(1._wp/dx(j))) 
                            end do

                            if(num_fluids > 1) then 
                                !$acc loop seq
                                do i = 1, num_fluids - 1
                                    @:ATOMIC_ADD(rhs_vf(advxb+i-1)%sf(j+1,k,l),(0.5_wp * (alpha_R(i) * &
                                        vel_R(1))*(1._wp/dx(j+1)) + &
                                        0.5_wp*cfl*(alpha_R(i))*(1._wp/dx(j+1)))*dt)
                                    ! !$omp atomic update
                                    ! !$acc atomic
                                    ! rhs_vf(advxb+i-1)%sf(j+1,k,l) = rhs_vf(advxb+i-1)%sf(j+1,k,l) + &
                                    !     (0.5_wp * (alpha_R(i) * &
                                    !     vel_R(1))*(1._wp/dx(j+1)) + &
                                    !     0.5_wp*cfl*(alpha_R(i))*(1._wp/dx(j+1)))

                                    @:ATOMIC_ADD(rhs_vf(advxb+i-1)%sf(j+1,k,l),(-0.5_wp * q_prim_vf(advxb+i-1)%sf(j+1,k,l)  * vel_R(1)*(1._wp/dx(j+1)))*dt)
                                    ! !$omp atomic update
                                    ! !$acc atomic
                                    ! rhs_vf(advxb+i-1)%sf(j+1,k,l) = rhs_vf(advxb+i-1)%sf(j+1,k,l) &
                                    ! - (0.5_wp * q_prim_vf(advxb+i-1)%sf(j+1,k,l)  * vel_R(1)*(1._wp/dx(j+1)))

                                    @:ATOMIC_ADD(rhs_vf(advxb+i-1)%sf(j,k,l),(-0.5_wp * (alpha_R(i) * &
                                        vel_R(1))*(1._wp/dx(j)) - &
                                        0.5_wp*cfl*(alpha_R(i))*(1._wp/dx(j)))*dt)
                                    ! !$omp atomic update
                                    ! !$acc atomic
                                    ! rhs_vf(advxb+i-1)%sf(j,k,l) = rhs_vf(advxb+i-1)%sf(j,k,l) - &
                                    !     (0.5_wp * (alpha_R(i) * &
                                    !     vel_R(1))*(1._wp/dx(j)) + &
                                    !     0.5_wp*cfl*(alpha_R(i))*(1._wp/dx(j)))

                                    @:ATOMIC_ADD(rhs_vf(advxb+i-1)%sf(j,k,l),(0.5_wp * q_prim_vf(advxb+i-1)%sf(j,k,l) * vel_R(1)*(1._wp/dx(j)))*dt)
                                    ! !$omp atomic update
                                    ! !$acc atomic
                                    ! rhs_vf(advxb+i-1)%sf(j,k,l) = rhs_vf(advxb+i-1)%sf(j,k,l) &
                                    ! + (0.5_wp * q_prim_vf(advxb+i-1)%sf(j,k,l) * vel_R(1)*(1._wp/dx(j)))
                                end do
                            end if

                            @:ATOMIC_ADD(rhs_vf(momxb)%sf(j+1,k,l),(0.5_wp* (rho_R * (vel_R(1))**2.0 + &
                                 pres_R)*(1._wp/dx(j+1)) + &
                                 0.5_wp*cfl * (rho_R*vel_R(1))*(1._wp/dx(j+1)))*dt)
                            ! !$omp atomic update
                            ! !$acc atomic
                            ! rhs_vf(momxb)%sf(j+1,k,l) = rhs_vf(momxb)%sf(j+1,k,l) + &
                            !      (0.5_wp* (rho_R * (vel_R(1))**2.0 + &
                            !      pres_R)*(1._wp/dx(j+1)) + &
                            !      0.5_wp*cfl * (rho_R*vel_R(1))*(1._wp/dx(j+1)))

                            @:ATOMIC_ADD(rhs_vf(momxb+1)%sf(j+1, k, l),(0.5_wp * rho_R * vel_R(1)*vel_R(2)*(1._wp/dx(j+1)) + &
                                0.5_wp*cfl * (rho_R*vel_R(2))*(1._wp/dx(j+1)))*dt)
                            ! !$omp atomic update
                            ! !$acc atomic
                            ! rhs_vf(momxb+1)%sf(j+1, k, l) =  rhs_vf(momxb+1)%sf(j+1, k, l) + &
                            !     (0.5_wp * rho_R * vel_R(1)*vel_R(2)*(1._wp/dx(j+1)) + &
                            !     0.5_wp*cfl * (rho_R*vel_R(2))*(1._wp/dx(j+1)))

                            @:ATOMIC_ADD(rhs_vf(momxb+2)%sf(j+1, k, l),(0.5_wp * rho_R * vel_R(1)*vel_R(3)*(1._wp/dx(j+1)) + &
                                0.5_wp*cfl * (rho_R*vel_R(3))*(1._wp/dx(j+1)))*dt)
                            ! !$omp atomic update
                            ! !$acc atomic
                            ! rhs_vf(momxb+2)%sf(j+1, k, l) =  rhs_vf(momxb+2)%sf(j+1, k, l) + &
                            !     (0.5_wp * rho_R * vel_R(1)*vel_R(3)*(1._wp/dx(j+1)) + &
                            !     0.5_wp*cfl * (rho_R*vel_R(3))*(1._wp/dx(j+1)))

                            @:ATOMIC_ADD(rhs_vf(E_idx)%sf(j+1, k, l),(0.5_wp * (vel_R(1) * (E_R + &
                                pres_R) )*(1._wp/dx(j+1)) + &
                                0.5_wp*cfl * (E_R)*(1._wp/dx(j+1)))*dt)
                            ! !$omp atomic update
                            ! !$acc atomic
                            !  rhs_vf(E_idx)%sf(j+1, k, l) = rhs_vf(E_idx)%sf(j+1, k, l) +  &
                            !     (0.5_wp * (vel_R(1) * (E_R + &
                            !     pres_R) )*(1._wp/dx(j+1)) + &
                            !     0.5_wp*cfl * (E_R)*(1._wp/dx(j+1)))

                            @:ATOMIC_ADD(rhs_vf(momxb)%sf(j,k,l),(-0.5_wp* (rho_R * (vel_R(1))**2.0 + &
                                 pres_R)*(1._wp/dx(j)) - &
                                 0.5_wp*cfl * (rho_R*vel_R(1))*(1._wp/dx(j)))*dt)
                            ! !$omp atomic update
                            ! !$acc atomic
                            ! rhs_vf(momxb)%sf(j,k,l) = rhs_vf(momxb)%sf(j,k,l) - &
                            !      (0.5_wp* (rho_R * (vel_R(1))**2.0 + &
                            !      pres_R)*(1._wp/dx(j)) + &
                            !      0.5_wp*cfl * (rho_R*vel_R(1))*(1._wp/dx(j)))

                            @:ATOMIC_ADD(rhs_vf(momxb+1)%sf(j, k, l),(-0.5_wp * rho_R * vel_R(1)*vel_R(2)*(1._wp/dx(j)) - &
                                0.5_wp*cfl * (rho_R*vel_R(2))*(1._wp/dx(j)))*dt)
                            ! !$omp atomic update
                            ! !$acc atomic
                            ! rhs_vf(momxb+1)%sf(j, k, l) =  rhs_vf(momxb+1)%sf(j, k, l) - &
                            !     (0.5_wp * rho_R * vel_R(1)*vel_R(2)*(1._wp/dx(j)) + &
                            !     0.5_wp*cfl * (rho_R*vel_R(2))*(1._wp/dx(j)))

                            @:ATOMIC_ADD(rhs_vf(momxb+2)%sf(j, k, l),(-0.5_wp * rho_R * vel_R(1)*vel_R(3)*(1._wp/dx(j)) - &
                                0.5_wp*cfl * (rho_R*vel_R(3))*(1._wp/dx(j)))*dt)
                            ! !$omp atomic update
                            ! !$acc atomic
                            ! rhs_vf(momxb+2)%sf(j, k, l) =  rhs_vf(momxb+2)%sf(j, k, l) - &
                            !     (0.5_wp * rho_R * vel_R(1)*vel_R(3)*(1._wp/dx(j)) + &
                            !     0.5_wp*cfl * (rho_R*vel_R(3))*(1._wp/dx(j)))

                            @:ATOMIC_ADD(rhs_vf(E_idx)%sf(j, k, l),(-0.5_wp * (vel_R(1) * (E_R + &
                                pres_R) )*(1._wp/dx(j)) - &
                                0.5_wp*cfl * (E_R)*(1._wp/dx(j)))*dt)
                            ! !$omp atomic update
                            ! !$acc atomic
                            !  rhs_vf(E_idx)%sf(j, k, l) = rhs_vf(E_idx)%sf(j, k, l) -  &
                            !     (0.5_wp * (vel_R(1) * (E_R + &
                            !     pres_R) )*(1._wp/dx(j)) + &
                            !     0.5_wp*cfl * (E_R)*(1._wp/dx(j)))
                        end do
                    end do
                end do
            end if
        else if (idir == 2) then
            if(p == 0) then
                !$omp target teams loop collapse(3) private(rho_L,gamma_L,pi_inf_L,mu_L,a_L,vel_L,vel_R,pres_L,alpha_L,alpha_R,alpha_rho_L,cfl,F_L,E_L,mu_R,rho_sf_small,alpha_rho_R,vflux_L_arr,vflux_R_arr,dvel_small)
                !$acc parallel loop collapse(3) gang vector default(present) private(rho_L,gamma_L,pi_inf_L,mu_L,a_L,vel_L,vel_R,pres_L,alpha_L,alpha_R,alpha_rho_L,cfl,F_L,E_L,mu_R,rho_sf_small,alpha_rho_R,vflux_L_arr,vflux_R_arr,dvel_small)
                do l = 0, p
                    do k = -buff_size+5, n+buff_size-5
                        do j = 0, m

                            if(viscous) then 
                                vflux_L_arr = 0._wp
                                vflux_R_arr = 0._wp

                                !DIR$ unroll 6
                                !$acc loop seq 
                                do q = -2, 3
                                    dvel_small = 0._wp
                                    !x-direction contributions
                                    !$acc loop seq 
                                    do i = -1, 1
                                        rho_L = 0._wp
                                        !$acc loop seq 
                                        do r = 1, num_fluids
                                            rho_L = rho_L + q_prim_vf(r)%sf(j+i,k+q,l)
                                        end do
                                        rho_sf_small(i) = rho_L
                                    end do

                                    dvel_small(1) = (1/(2._wp*dx(j))) * ( &
                                        q_prim_vf(momxb)%sf(j+1,k+q,l)/rho_sf_small(1) - &
                                        q_prim_vf(momxb)%sf(j-1,k+q,l)/rho_sf_small(-1))
                                    dvel_small(2) = (1/(2._wp*dx(j))) * ( &
                                        q_prim_vf(momxb+1)%sf(j+1,k+q,l)/rho_sf_small(1) - &
                                        q_prim_vf(momxb+1)%sf(j-1,k+q,l)/rho_sf_small(-1))

                                    if (q > -2) then
                                        vflux_L_arr(1) = vflux_L_arr(1) + coeff_L(q)*(dvel_small(2))
                                        vflux_L_arr(3) = vflux_L_arr(3) + coeff_L(q)*(-2._wp*dvel_small(1))/3._wp
                                    end if
                                    if (q < 3) then
                                        vflux_R_arr(1) = vflux_R_arr(1) + coeff_R(q)*(dvel_small(2))
                                        vflux_R_arr(3) = vflux_R_arr(3) + coeff_R(q)*(-2._wp*dvel_small(1))/3._wp
                                    end if

                                    !y-direction contributions
                                    !$acc loop seq 
                                    do i = -1, 1
                                        rho_L = 0._wp
                                        !$acc loop seq 
                                        do r = 1, num_fluids
                                            rho_L = rho_L + q_prim_vf(r)%sf(j,k+i+q,l)
                                        end do
                                        rho_sf_small(i) = rho_L
                                    end do

                                    dvel_small(1) = (1/(2._wp*dy(k))) * ( &
                                        q_prim_vf(momxb)%sf(j,k+1+q,l)/rho_sf_small(1) - &
                                        q_prim_vf(momxb)%sf(j,k-1+q,l)/rho_sf_small(-1))
                                    dvel_small(2) = (1/(2._wp*dy(k))) * ( &
                                        q_prim_vf(momxb+1)%sf(j,k+1+q,l)/rho_sf_small(1) - &
                                        q_prim_vf(momxb+1)%sf(j,k-1+q,l)/rho_sf_small(-1))

                                    if (q > -2) then
                                        vflux_L_arr(1) = vflux_L_arr(1) + coeff_L(q)*(dvel_small(1))
                                        vflux_L_arr(3) = vflux_L_arr(3) + coeff_L(q)*(4._wp*dvel_small(2))/3._wp
                                    end if
                                    if (q < 3) then
                                        vflux_R_arr(1) = vflux_R_arr(1) + coeff_R(q)*(dvel_small(1))
                                        vflux_R_arr(3) = vflux_R_arr(3) + coeff_R(q)*(4._wp*dvel_small(2))/3._wp
                                    end if

                                end do
                            end if

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
                            end do 

                            if(num_fluids > 1) then  
                                !$acc loop seq 
                                do i = 1,num_fluids - 1
                                    alpha_L(i) =  (1._wp/60._wp) * (-3._wp * q_prim_vf(E_idx+i)%sf(j, k-1, l) + &
                                                27._wp * q_prim_vf(E_idx+i)%sf(j, k, l) + &
                                                47._wp * q_prim_vf(E_idx+i)%sf(j, k+1, l) -   &
                                                13._wp * q_prim_vf(E_idx+i)%sf(j, k+2, l) + &
                                                2._wp * q_prim_vf(E_idx+i)%sf(j, k+3, l))                          
                                    alpha_R(i) = (1._wp/60._wp) * (-3._wp * q_prim_vf(E_idx+i)%sf(j, k+2, l) + &
                                                        27._wp * q_prim_vf(E_idx+i)%sf(j, k+1, l) + &
                                                        47._wp * q_prim_vf(E_idx+i)%sf(j, k, l) -   &
                                                        13._wp * q_prim_vf(E_idx+i)%sf(j, k-1, l) + &
                                                        2._wp * q_prim_vf(E_idx+i)%sf(j, k-2, l))
                                end do
                                alpha_L(num_fluids) = 1._wp - sum(alpha_L(1:num_fluids-1))
                                alpha_R(num_fluids) = 1._wp - sum(alpha_R(1:num_fluids-1))
                            else 
                                alpha_L(1) = 1._wp
                                alpha_R(1) = 1._wp 
                            end if

                            rho_L = sum(alpha_rho_L)
                            gamma_L = sum(alpha_L*gammas)
                            pi_inf_L = sum(alpha_L*pi_infs)

                            rho_R = sum(alpha_rho_R)
                            gamma_R = sum(alpha_R*gammas)
                            pi_inf_R = sum(alpha_R*pi_infs)

                            !$acc loop seq 
                            do i = 1, num_dims
                                vel_L(i) = (1._wp/60._wp) * (-3._wp * q_prim_vf(momxb+i-1)%sf(j,k-1,l) + &
                                                27._wp * q_prim_vf(momxb+i-1)%sf(j,k,l) + &
                                                47._wp * q_prim_vf(momxb+i-1)%sf(j,k+1,l) -   &
                                                13._wp * q_prim_vf(momxb+i-1)%sf(j,k+2,l) + &
                                                2._wp * q_prim_vf(momxb+i-1)%sf(j,k+3,l)) / rho_L
                            end do 

                            !$acc loop seq 
                            do i = 1, num_dims
                                vel_R(i) = (1._wp/60._wp) * (-3._wp * q_prim_vf(momxb+i-1)%sf(j,k+2,l) + &
                                                27._wp * q_prim_vf(momxb+i-1)%sf(j,k+1,l) + &
                                                47._wp * q_prim_vf(momxb+i-1)%sf(j,k,l) -   &
                                                13._wp * q_prim_vf(momxb+i-1)%sf(j,k-1,l) + &
                                                2._wp * q_prim_vf(momxb+i-1)%sf(j,k-2,l)) / rho_R
                            end do

                            if(viscous) then

                                if(num_fluids > 1) then 
                                    mu_L = 0._wp
                                    mu_R = 0._wp
                                    !$acc loop seq 
                                    do i = 1, num_fluids - 1
                                        mu_L = (1._wp/60._wp) * (-3._wp * q_prim_vf(E_idx+i)%sf(j, k-1, l) + &
                                                27._wp * q_prim_vf(E_idx+i)%sf(j, k, l) + &
                                                47._wp * q_prim_vf(E_idx+i)%sf(j, k+1, l) -   &
                                                13._wp * q_prim_vf(E_idx+i)%sf(j, k+2, l) + &
                                                2._wp * q_prim_vf(E_idx+i)%sf(j, k+3, l)) / Res(1, i) + mu_L

                                        mu_R = (1._wp/60._wp) * (-3._wp * q_prim_vf(E_idx+i)%sf(j, k+2, l) + &
                                                27._wp * q_prim_vf(E_idx+i)%sf(j, k+1, l) + &
                                                47._wp * q_prim_vf(E_idx+i)%sf(j, k, l) -   &
                                                13._wp * q_prim_vf(E_idx+i)%sf(j, k-1, l) + &
                                                2._wp * q_prim_vf(E_idx+i)%sf(j, k-2, l)) / Res(1, i) + mu_R
                                    end do

                                    mu_L = 1._wp / Res(1, num_fluids) + mu_L
                                    mu_R = 1._wp / Res(1, num_fluids) + mu_R

                                    !$acc loop seq 
                                    do i = 1, num_fluids - 1
                                        mu_L = -(1._wp/60._wp) * (-3._wp * q_prim_vf(E_idx+i)%sf(j, k-1, l) + &
                                                27._wp * q_prim_vf(E_idx+i)%sf(j, k, l) + &
                                                47._wp * q_prim_vf(E_idx+i)%sf(j, k+1, l) -   &
                                                13._wp * q_prim_vf(E_idx+i)%sf(j, k+2, l) + &
                                                2._wp * q_prim_vf(E_idx+i)%sf(j, k+3, l)) / Res(1, num_fluids) + mu_L

                                        mu_R = -(1._wp/60._wp) * (-3._wp * q_prim_vf(E_idx+i)%sf(j, k+2, l) + &
                                                27._wp * q_prim_vf(E_idx+i)%sf(j, k+1, l) + &
                                                47._wp * q_prim_vf(E_idx+i)%sf(j, k, l) -   &
                                                13._wp * q_prim_vf(E_idx+i)%sf(j, k-1, l) + &
                                                2._wp * q_prim_vf(E_idx+i)%sf(j, k-2, l)) / Res(1, num_fluids) + mu_R
                                    end do
                                else
                                    mu_L = 1._wp / Res(1, 1) 
                                    mu_R = 1._wp / Res(1, 1) 
                                end if


                                @:ATOMIC_ADD(rhs_vf(momxb)%sf(j,k+1,l),(- 0.5_wp*mu_L*vflux_L_arr(1)*(1._wp/dy(k+1))*dt))
                                !!$omp atomic update
                                !!$acc atomic
                                !rhs_vf(momxb)%sf(j,k+1,l) = rhs_vf(momxb)%sf(j,k+1,l) - 0.5_wp*mu_L*vflux_L_arr(1)*(1._wp/dy(k+1))
                                @:ATOMIC_ADD(rhs_vf(E_idx)%sf(j,k+1,l),(- 0.5_wp*mu_L*vflux_L_arr(1)*vel_L(1)*(1._wp/dy(k+1))*dt))
                                !!$omp atomic update
                                !!$acc atomic
                                !rhs_vf(E_idx)%sf(j,k+1,l) = rhs_vf(E_idx)%sf(j,k+1,l) - 0.5_wp*mu_L*vflux_L_arr(1)*vel_L(1)*(1._wp/dy(k+1)) 

                                @:ATOMIC_ADD(rhs_vf(momxb)%sf(j,k,l),(0.5_wp*mu_L*vflux_L_arr(1)*(1._wp/dy(k))*dt))
                                !!$omp atomic update
                                !!$acc atomic
                                !rhs_vf(momxb)%sf(j,k,l) = rhs_vf(momxb)%sf(j,k,l) + 0.5_wp*mu_L*vflux_L_arr(1)*(1._wp/dy(k))
                                @:ATOMIC_ADD(rhs_vf(E_idx)%sf(j,k,l),(0.5_wp*mu_L*vflux_L_arr(1)*vel_L(1)*(1._wp/dy(k))*dt))
                                !!$omp atomic update
                                !!$acc atomic
                                !rhs_vf(E_idx)%sf(j,k,l) = rhs_vf(E_idx)%sf(j,k,l) + 0.5_wp*mu_L*vflux_L_arr(1)*vel_L(1)*(1._wp/dy(k))
                            
                                @:ATOMIC_ADD(rhs_vf(momxb)%sf(j,k+1,l),(- 0.5_wp*mu_R*vflux_R_arr(1)*(1._wp/dy(k+1))*dt))
                                !!$omp atomic update
                                !!$acc atomic
                                !rhs_vf(momxb)%sf(j,k+1,l) = rhs_vf(momxb)%sf(j,k+1,l) - 0.5_wp*mu_R*vflux_R_arr(1)*(1._wp/dy(k+1))
                                @:ATOMIC_ADD(rhs_vf(E_idx)%sf(j,k+1,l),(- 0.5_wp*mu_R*vflux_R_arr(1)*vel_R(1)*(1._wp/dy(k+1))*dt))
                                !!$omp atomic update
                                !!$acc atomic
                                !rhs_vf(E_idx)%sf(j,k+1,l) = rhs_vf(E_idx)%sf(j,k+1,l) - 0.5_wp*mu_R*vflux_R_arr(1)*vel_R(1)*(1._wp/dy(k+1))
                            
                                @:ATOMIC_ADD(rhs_vf(momxb)%sf(j,k,l),(0.5_wp*mu_R*vflux_R_arr(1)*(1._wp/dy(k))*dt))
                                !!$omp atomic update
                                !!$acc atomic
                                !rhs_vf(momxb)%sf(j,k,l) = rhs_vf(momxb)%sf(j,k,l) + 0.5_wp*mu_R*vflux_R_arr(1)*(1._wp/dy(k))
                                @:ATOMIC_ADD(rhs_vf(E_idx)%sf(j,k,l),(0.5_wp*mu_R*vflux_R_arr(1)*vel_R(1)*(1._wp/dy(k))*dt))
                                !!$omp atomic update
                                !!$acc atomic
                                !rhs_vf(E_idx)%sf(j,k,l) = rhs_vf(E_idx)%sf(j,k,l) + 0.5_wp*mu_R*vflux_R_arr(1)*vel_R(1)*(1._wp/dy(k))


                                @:ATOMIC_ADD(rhs_vf(momxb+1)%sf(j,k+1,l),(- 0.5_wp*mu_L*vflux_L_arr(3)*(1._wp/dy(k+1))*dt))
                                !!$omp atomic update
                                !!$acc atomic
                                !rhs_vf(momxb+1)%sf(j,k+1,l) = rhs_vf(momxb+1)%sf(j,k+1,l) - 0.5_wp*mu_L*vflux_L_arr(3)*(1._wp/dy(k+1))
                                @:ATOMIC_ADD(rhs_vf(E_idx)%sf(j,k+1,l),(- 0.5_wp*mu_L*vflux_L_arr(3)*vel_L(2)*(1._wp/dy(k+1))*dt))
                                !!$omp atomic update
                                !!$acc atomic
                                !rhs_vf(E_idx)%sf(j,k+1,l) = rhs_vf(E_idx)%sf(j,k+1,l) - 0.5_wp*mu_L*vflux_L_arr(3)*vel_L(2)*(1._wp/dy(k+1))
                        
                                @:ATOMIC_ADD(rhs_vf(momxb+1)%sf(j,k,l),(0.5_wp*mu_L*vflux_L_arr(3)*(1._wp/dy(k))*dt))
                                !!$omp atomic update
                                !!$acc atomic
                                !rhs_vf(momxb+1)%sf(j,k,l) = rhs_vf(momxb+1)%sf(j,k,l) + 0.5_wp*mu_L*vflux_L_arr(3)*(1._wp/dy(k))
                                @:ATOMIC_ADD(rhs_vf(E_idx)%sf(j,k,l),(0.5_wp*mu_L*vflux_L_arr(3)*vel_L(2)*(1._wp/dy(k))*dt))
                                !!$omp atomic update
                                !!$acc atomic
                                !rhs_vf(E_idx)%sf(j,k,l) = rhs_vf(E_idx)%sf(j,k,l) + 0.5_wp*mu_L*vflux_L_arr(3)*vel_L(2)*(1._wp/dy(k))                    
                             
                                @:ATOMIC_ADD(rhs_vf(momxb+1)%sf(j,k+1,l),(- 0.5_wp*mu_R*vflux_R_arr(3)*(1._wp/dy(k+1))*dt))
                                !!$omp atomic update
                                !!$acc atomic
                                !rhs_vf(momxb+1)%sf(j,k+1,l) = rhs_vf(momxb+1)%sf(j,k+1,l) - 0.5_wp*mu_R*vflux_R_arr(3)*(1._wp/dy(k+1))
                                @:ATOMIC_ADD(rhs_vf(E_idx)%sf(j,k+1,l),(- 0.5_wp*mu_R*vflux_R_arr(3)*vel_R(2)*(1._wp/dy(k+1))*dt))
                                !!$omp atomic update
                                !!$acc atomic
                                !rhs_vf(E_idx)%sf(j,k+1,l) = rhs_vf(E_idx)%sf(j,k+1,l) - 0.5_wp*mu_R*vflux_R_arr(3)*vel_R(2)*(1._wp/dy(k+1))
                             
                                @:ATOMIC_ADD(rhs_vf(momxb+1)%sf(j,k,l),(0.5_wp*mu_R*vflux_R_arr(3)*(1._wp/dy(k))*dt))
                                !!$omp atomic update
                                !!$acc atomic
                                !rhs_vf(momxb+1)%sf(j,k,l) = rhs_vf(momxb+1)%sf(j,k,l) + 0.5_wp*mu_R*vflux_R_arr(3)*(1._wp/dy(k))
                                @:ATOMIC_ADD(rhs_vf(E_idx)%sf(j,k,l),(0.5_wp*mu_R*vflux_R_arr(3)*vel_R(2)*(1._wp/dy(k))*dt))
                                !!$omp atomic update
                                !!$acc atomic
                                !rhs_vf(E_idx)%sf(j,k,l) = rhs_vf(E_idx)%sf(j,k,l) + 0.5_wp*mu_R*vflux_R_arr(3)*vel_R(2)*(1._wp/dy(k))
                            
                            end if

                            F_L = (1._wp/60._wp) * (-3._wp * jac(j, k-1, l) + &
                                                27._wp * jac(j, k, l) + &
                                                47._wp * jac(j, k+1, l) -   &
                                                13._wp * jac(j, k+2, l) + &
                                                2._wp * jac(j, k+3, l))

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

                            pres_L = (E_L - pi_inf_L - 0.5_wp*rho_L*(vel_L(1)**2._wp + vel_L(2)**2._wp ))/gamma_L                    

                            pres_R = (E_R - pi_inf_R - 0.5_wp*rho_R*(vel_R(1)**2._wp + vel_R(2)**2._wp ))/gamma_R 

                            a_L = sqrt((pres_L*(1._wp/gamma_L + 1._wp) + pi_inf_L / gamma_L) / rho_L)
                            a_R = sqrt((pres_R*(1._wp/gamma_R + 1._wp) + pi_inf_R / gamma_R) / rho_R)

                            cfl = (max(sqrt(vel_L(1)**2._wp + vel_L(2)**2._wp),sqrt(vel_R(1)**2._wp + vel_R(2)**2._wp)) + max(a_L,a_R))


                            !$acc loop seq
                            do i = 1, num_fluids
                                @:ATOMIC_ADD(rhs_vf(i)%sf(j,k+1,l),(0.5_wp * (alpha_rho_L(i) * &
                                    vel_L(2))*(1._wp/dy(k+1)) - &
                                    0.5_wp*cfl * (alpha_rho_L(i))*(1._wp/dy(k+1)))*dt)
                                ! !$omp atomic update
                                ! !$acc atomic
                                ! rhs_vf(i)%sf(j,k+1,l) = rhs_vf(i)%sf(j,k+1,l) + &
                                !     (0.5_wp * (alpha_rho_L(i) * &
                                !     vel_L(2))*(1._wp/dy(k+1)) - &
                                !     0.5_wp*cfl * (alpha_rho_L(i))*(1._wp/dy(k+1)))

                                @:ATOMIC_ADD(rhs_vf(i)%sf(j,k,l),(-0.5_wp * (alpha_rho_L(i) * &
                                    vel_L(2))*(1._wp/dy(k)) + &
                                    0.5_wp*cfl * (alpha_rho_L(i))*(1._wp/dy(k)))*dt)
                                ! !$omp atomic update
                                ! !$acc atomic
                                ! rhs_vf(i)%sf(j,k,l) = rhs_vf(i)%sf(j,k,l) - &
                                !     (0.5_wp * (alpha_rho_L(i) * &
                                !     vel_L(2))*(1._wp/dy(k)) - &
                                !     0.5_wp*cfl * (alpha_rho_L(i))*(1._wp/dy(k))) 
                            end do

                            if(num_fluids > 1) then 
                                !$acc loop seq
                                do i = 1, num_fluids - 1
                                    @:ATOMIC_ADD(rhs_vf(advxb+i-1)%sf(j,k+1,l),(0.5_wp * (alpha_L(i) * &
                                        vel_L(2))*(1._wp/dy(k+1)) - &
                                        0.5_wp*cfl*(alpha_L(i))*(1._wp/dy(k+1)))*dt)
                                    ! !$omp atomic update
                                    ! !$acc atomic
                                    ! rhs_vf(advxb+i-1)%sf(j,k+1,l) = rhs_vf(advxb+i-1)%sf(j,k+1,l) + &
                                    !     (0.5_wp * (alpha_L(i) * &
                                    !     vel_L(2))*(1._wp/dy(k+1)) - &
                                    !     0.5_wp*cfl*(alpha_L(i))*(1._wp/dy(k+1)))

                                    @:ATOMIC_ADD(rhs_vf(advxb+i-1)%sf(j,k+1,l),(-0.5_wp * q_prim_vf(advxb+i-1)%sf(j,k+1,l)  * vel_L(2)*(1._wp/dy(k+1)))*dt)
                                    ! !$omp atomic update
                                    ! !$acc atomic
                                    ! rhs_vf(advxb+i-1)%sf(j,k+1,l) = rhs_vf(advxb+i-1)%sf(j,k+1,l) &
                                    ! - (0.5_wp * q_prim_vf(advxb+i-1)%sf(j,k+1,l)  * vel_L(2)*(1._wp/dy(k+1)))

                                    @:ATOMIC_ADD(rhs_vf(advxb+i-1)%sf(j,k,l),(-0.5_wp * (alpha_L(i) * &
                                        vel_L(2))*(1._wp/dy(k)) + &
                                        0.5_wp*cfl*(alpha_L(i))*(1._wp/dy(k)))*dt)
                                    ! !$omp atomic update
                                    ! !$acc atomic
                                    ! rhs_vf(advxb+i-1)%sf(j,k,l) = rhs_vf(advxb+i-1)%sf(j,k,l) - &
                                    !     (0.5_wp * (alpha_L(i) * &
                                    !     vel_L(2))*(1._wp/dy(k)) - &
                                    !     0.5_wp*cfl*(alpha_L(i))*(1._wp/dy(k)))

                                    @:ATOMIC_ADD(rhs_vf(advxb+i-1)%sf(j,k,l),(0.5_wp * q_prim_vf(advxb+i-1)%sf(j,k,l)  * vel_L(2)*(1._wp/dy(k)))*dt)
                                    ! !$omp atomic update
                                    ! !$acc atomic
                                    ! rhs_vf(advxb+i-1)%sf(j,k,l) = rhs_vf(advxb+i-1)%sf(j,k,l) &
                                    ! + (0.5_wp * q_prim_vf(advxb+i-1)%sf(j,k,l)  * vel_L(2)*(1._wp/dy(k)))
                                end do
                            end if

                            @:ATOMIC_ADD(rhs_vf(momxb+1)%sf(j,k+1,l),(0.5_wp* (rho_L * (vel_L(2))**2.0 + &
                                 pres_L+F_L)*(1._wp/dy(k+1)) - &
                                 0.5_wp*cfl * (rho_L*vel_L(2))*(1._wp/dy(k+1)))*dt)
                            ! !$omp atomic update
                            ! !$acc atomic
                            ! rhs_vf(momxb+1)%sf(j,k+1,l) = rhs_vf(momxb+1)%sf(j,k+1,l) + &
                            !      (0.5_wp* (rho_L * (vel_L(2))**2.0 + &
                            !      pres_L+F_L)*(1._wp/dy(k+1)) - &
                            !      0.5_wp*cfl * (rho_L*vel_L(2))*(1._wp/dy(k+1)))

                            @:ATOMIC_ADD(rhs_vf(momxb)%sf(j, k+1, l),(0.5_wp * rho_L * vel_L(1)*vel_L(2)*(1._wp/dy(k+1)) - &
                                0.5_wp*cfl * (rho_L*vel_L(1))*(1._wp/dy(k+1)))*dt)
                            ! !$omp atomic update
                            ! !$acc atomic
                            ! rhs_vf(momxb)%sf(j, k+1, l) =  rhs_vf(momxb)%sf(j, k+1, l) + &
                            !     (0.5_wp * rho_L * vel_L(1)*vel_L(2)*(1._wp/dy(k+1)) - &
                            !     0.5_wp*cfl * (rho_L*vel_L(1))*(1._wp/dy(k+1)))

                            @:ATOMIC_ADD(rhs_vf(E_idx)%sf(j, k+1, l),(0.5_wp * (vel_L(2) * (E_L + &
                                pres_L+F_L) )*(1._wp/dy(k+1)) - &
                                0.5_wp*cfl * (E_L)*(1._wp/dy(k+1)))*dt)
                            ! !$omp atomic update
                            ! !$acc atomic
                            ! rhs_vf(E_idx)%sf(j, k+1, l) = rhs_vf(E_idx)%sf(j, k+1, l) +  &
                            !     (0.5_wp * (vel_L(2) * (E_L + &
                            !     pres_L+F_L) )*(1._wp/dy(k+1)) - &
                            !     0.5_wp*cfl * (E_L)*(1._wp/dy(k+1)))

                            @:ATOMIC_ADD(rhs_vf(momxb+1)%sf(j,k,l),(-0.5_wp* (rho_L * (vel_L(2))**2.0 + &
                                 pres_L+F_L)*(1._wp/dy(k)) + &
                                 0.5_wp*cfl * (rho_L*vel_L(2))*(1._wp/dy(k)))*dt)
                            ! !$omp atomic update
                            ! !$acc atomic
                            ! rhs_vf(momxb+1)%sf(j,k,l) = rhs_vf(momxb+1)%sf(j,k,l) - &
                            !      (0.5_wp* (rho_L * (vel_L(2))**2.0 + &
                            !      pres_L+F_L)*(1._wp/dy(k)) - &
                            !      0.5_wp*cfl * (rho_L*vel_L(2))*(1._wp/dy(k)))

                            @:ATOMIC_ADD(rhs_vf(momxb)%sf(j, k, l),(-0.5_wp * rho_L * vel_L(1)*vel_L(2)*(1._wp/dy(k)) + &
                                0.5_wp*cfl * (rho_L*vel_L(1))*(1._wp/dy(k)))*dt)
                            ! !$omp atomic update
                            ! !$acc atomic
                            ! rhs_vf(momxb)%sf(j, k, l) =  rhs_vf(momxb)%sf(j, k, l) - &
                            !     (0.5_wp * rho_L * vel_L(1)*vel_L(2)*(1._wp/dy(k)) - &
                            !     0.5_wp*cfl * (rho_L*vel_L(1))*(1._wp/dy(k)))


                            @:ATOMIC_ADD(rhs_vf(E_idx)%sf(j, k, l),(-0.5_wp * (vel_L(2) * (E_L + &
                                pres_L+F_L) )*(1._wp/dy(k)) + &
                                0.5_wp*cfl * (E_L)*(1._wp/dy(k)))*dt)
                            ! !$omp atomic update
                            ! !$acc atomic
                            ! rhs_vf(E_idx)%sf(j, k, l) = rhs_vf(E_idx)%sf(j, k, l) -  &
                            !     (0.5_wp * (vel_L(2) * (E_L + &
                            !     pres_L+F_L) )*(1._wp/dy(k)) - &
                            !     0.5_wp*cfl * (E_L)*(1._wp/dy(k)))


                            F_L = (1._wp/60._wp) * (-3._wp * jac(j, k+2, l) + &
                                                27._wp * jac(j, k+1, l) + &
                                                47._wp * jac(j, k, l) -   &
                                                13._wp * jac(j, k-1, l) + &
                                                2._wp * jac(j, k-2, l))

                            !$acc loop seq
                            do i = 1, num_fluids
                                @:ATOMIC_ADD(rhs_vf(i)%sf(j,k+1,l),(0.5_wp * (alpha_rho_R(i) * &
                                    vel_R(2))*(1._wp/dy(k+1)) + &
                                    0.5_wp*cfl * (alpha_rho_R(i))*(1._wp/dy(k+1)))*dt)
                                ! !$omp atomic update
                                ! !$acc atomic
                                ! rhs_vf(i)%sf(j,k+1,l) = rhs_vf(i)%sf(j,k+1,l) + &
                                !     (0.5_wp * (alpha_rho_R(i) * &
                                !     vel_R(2))*(1._wp/dy(k+1)) + &
                                !     0.5_wp*cfl * (alpha_rho_R(i))*(1._wp/dy(k+1)))

                                @:ATOMIC_ADD(rhs_vf(i)%sf(j,k,l),(-0.5_wp * (alpha_rho_R(i) * &
                                    vel_R(2))*(1._wp/dy(k)) - &
                                    0.5_wp*cfl * (alpha_rho_R(i))*(1._wp/dy(k)))*dt)
                                ! !$omp atomic update
                                ! !$acc atomic
                                ! rhs_vf(i)%sf(j,k,l) = rhs_vf(i)%sf(j,k,l) - &
                                !     (0.5_wp * (alpha_rho_R(i) * &
                                !     vel_R(2))*(1._wp/dy(k)) + &
                                !     0.5_wp*cfl * (alpha_rho_R(i))*(1._wp/dy(k))) 
                            end do

                            if(num_fluids > 1) then 
                                !$acc loop seq
                                do i = 1, num_fluids - 1
                                    @:ATOMIC_ADD(rhs_vf(advxb+i-1)%sf(j,k+1,l),(0.5_wp * (alpha_R(i) * &
                                        vel_R(2))*(1._wp/dy(k+1)) + &
                                        0.5_wp*cfl*(alpha_R(i))*(1._wp/dy(k+1)))*dt)
                                    ! !$omp atomic update
                                    ! !$acc atomic
                                    ! rhs_vf(advxb+i-1)%sf(j,k+1,l) = rhs_vf(advxb+i-1)%sf(j,k+1,l) + &
                                    !     (0.5_wp * (alpha_R(i) * &
                                    !     vel_R(2))*(1._wp/dy(k+1)) + &
                                    !     0.5_wp*cfl*(alpha_R(i))*(1._wp/dy(k+1)))

                                    @:ATOMIC_ADD(rhs_vf(advxb+i-1)%sf(j,k+1,l),(-0.5_wp * q_prim_vf(advxb+i-1)%sf(j,k+1,l)  * vel_R(2)*(1._wp/dy(k+1)))*dt)
                                    ! !$omp atomic update
                                    ! !$acc atomic
                                    ! rhs_vf(advxb+i-1)%sf(j,k+1,l) = rhs_vf(advxb+i-1)%sf(j,k+1,l) &
                                    ! - (0.5_wp * q_prim_vf(advxb+i-1)%sf(j,k+1,l)  * vel_R(2)*(1._wp/dy(k+1)))

                                    @:ATOMIC_ADD(rhs_vf(advxb+i-1)%sf(j,k,l),(-0.5_wp * (alpha_R(i) * &
                                        vel_R(2))*(1._wp/dy(k)) - &
                                        0.5_wp*cfl*(alpha_R(i))*(1._wp/dy(k)))*dt)
                                    ! !$omp atomic update
                                    ! !$acc atomic
                                    ! rhs_vf(advxb+i-1)%sf(j,k,l) = rhs_vf(advxb+i-1)%sf(j,k,l) - &
                                    !     (0.5_wp * (alpha_R(i) * &
                                    !     vel_R(2))*(1._wp/dy(k)) + &
                                    !     0.5_wp*cfl*(alpha_R(i))*(1._wp/dy(k)))

                                    @:ATOMIC_ADD(rhs_vf(advxb+i-1)%sf(j,k,l),(0.5_wp * q_prim_vf(advxb+i-1)%sf(j,k,l)  * vel_R(2)*(1._wp/dy(k)))*dt)
                                    ! !$omp atomic update
                                    ! !$acc atomic
                                    ! rhs_vf(advxb+i-1)%sf(j,k,l) = rhs_vf(advxb+i-1)%sf(j,k,l) &
                                    ! + (0.5_wp * q_prim_vf(advxb+i-1)%sf(j,k,l)  * vel_R(2)*(1._wp/dy(k)))
                                end do
                            end if

                            @:ATOMIC_ADD(rhs_vf(momxb+1)%sf(j,k+1,l),(0.5_wp* (rho_R * (vel_R(2))**2.0 + &
                                 pres_R+F_L)*(1._wp/dy(k+1)) + &
                                 0.5_wp*cfl * (rho_R*vel_R(2))*(1._wp/dy(k+1)))*dt)
                            ! !$omp atomic update
                            ! !$acc atomic
                            ! rhs_vf(momxb+1)%sf(j,k+1,l) = rhs_vf(momxb+1)%sf(j,k+1,l) + &
                            !      (0.5_wp* (rho_R * (vel_R(2))**2.0 + &
                            !      pres_R+F_L)*(1._wp/dy(k+1)) + &
                            !      0.5_wp*cfl * (rho_R*vel_R(2))*(1._wp/dy(k+1)))

                            @:ATOMIC_ADD(rhs_vf(momxb)%sf(j, k+1, l),(0.5_wp * rho_R * vel_R(2)*vel_R(1)*(1._wp/dy(k+1)) + &
                                0.5_wp*cfl * (rho_R*vel_R(1))*(1._wp/dy(k+1)))*dt)
                            ! !$omp atomic update
                            ! !$acc atomic
                            ! rhs_vf(momxb)%sf(j, k+1, l) =  rhs_vf(momxb)%sf(j, k+1, l) + &
                            !     (0.5_wp * rho_R * vel_R(2)*vel_R(1)*(1._wp/dy(k+1)) + &
                            !     0.5_wp*cfl * (rho_R*vel_R(1))*(1._wp/dy(k+1)))


                            @:ATOMIC_ADD(rhs_vf(E_idx)%sf(j, k+1, l),(0.5_wp * (vel_R(2) * (E_R+ &
                                pres_R+F_L ) )*(1._wp/dy(k+1)) + &
                                0.5_wp*cfl * (E_R)*(1._wp/dy(k+1)))*dt)
                            ! !$omp atomic update
                            ! !$acc atomic
                            !  rhs_vf(E_idx)%sf(j, k+1, l) = rhs_vf(E_idx)%sf(j, k+1, l) +  &
                            !     (0.5_wp * (vel_R(2) * (E_R+ &
                            !     pres_R+F_L ) )*(1._wp/dy(k+1)) + &
                            !     0.5_wp*cfl * (E_R)*(1._wp/dy(k+1)))

                            @:ATOMIC_ADD(rhs_vf(momxb+1)%sf(j,k,l),(-0.5_wp* (rho_R * (vel_R(2))**2.0 + &
                                 pres_R+F_L)*(1._wp/dy(k)) - &
                                 0.5_wp*cfl * (rho_R*vel_R(2))*(1._wp/dy(k)))*dt)
                            ! !$omp atomic update
                            ! !$acc atomic
                            ! rhs_vf(momxb+1)%sf(j,k,l) = rhs_vf(momxb+1)%sf(j,k,l) - &
                            !      (0.5_wp* (rho_R * (vel_R(2))**2.0 + &
                            !      pres_R+F_L)*(1._wp/dy(k)) + &
                            !      0.5_wp*cfl * (rho_R*vel_R(2))*(1._wp/dy(k)))

                            @:ATOMIC_ADD(rhs_vf(momxb)%sf(j, k, l),(-0.5_wp * rho_R * vel_R(2)*vel_R(1)*(1._wp/dy(k)) - &
                                0.5_wp*cfl * (rho_R*vel_R(1))*(1._wp/dy(k)))*dt)
                            ! !$omp atomic update
                            ! !$acc atomic
                            ! rhs_vf(momxb)%sf(j, k, l) =  rhs_vf(momxb)%sf(j, k, l) - &
                            !     (0.5_wp * rho_R * vel_R(2)*vel_R(1)*(1._wp/dy(k)) + &
                            !     0.5_wp*cfl * (rho_R*vel_R(1))*(1._wp/dy(k)))


                            @:ATOMIC_ADD(rhs_vf(E_idx)%sf(j, k, l),(-0.5_wp * (vel_R(2) * (E_R + &
                                pres_R+F_L) )*(1._wp/dy(k)) - &
                                0.5_wp*cfl * (E_R)*(1._wp/dy(k)))*dt)
                            ! !$omp atomic update
                            ! !$acc atomic
                            !  rhs_vf(E_idx)%sf(j, k, l) = rhs_vf(E_idx)%sf(j, k, l) -  &
                            !     (0.5_wp * (vel_R(2) * (E_R + &
                            !     pres_R+F_L) )*(1._wp/dy(k)) + &
                            !     0.5_wp*cfl * (E_R)*(1._wp/dy(k)))                            

                        end do
                    end do
                end do
            else
                !$omp target teams loop collapse(3) private(rho_L,gamma_L,pi_inf_L,mu_L,a_L,vel_L,vel_R,pres_L,alpha_L,alpha_R,alpha_rho_L,cfl,F_L,E_L,mu_R,rho_sf_small,alpha_rho_R,vflux_L_arr,vflux_R_arr,dvel_small)
                !$acc parallel loop collapse(3) gang vector default(present) private(rho_L,gamma_L,pi_inf_L,mu_L,a_L,vel_L,vel_R,pres_L,alpha_L,alpha_R,alpha_rho_L,cfl,F_L,E_L,mu_R,rho_sf_small,alpha_rho_R,vflux_L_arr,vflux_R_arr,dvel_small)
                do l = 0, p
                    do k = -buff_size+5, n+buff_size-5
                        do j = 0, m

                            if(viscous) then 
                                vflux_L_arr = 0._wp
                                vflux_R_arr = 0._wp

                                !DIR$ unroll 6
                                !$acc loop seq 
                                do q = -2, 3
                                    dvel_small = 0._wp
                                    !x-direction contributions
                                    !$acc loop seq 
                                    do i = -1, 1
                                        rho_L = 0._wp
                                        !$acc loop seq 
                                        do r = 1, num_fluids
                                            rho_L = rho_L + q_prim_vf(r)%sf(j+i,k+q,l)
                                        end do
                                        rho_sf_small(i) = rho_L
                                    end do

                                    dvel_small(1) = (1/(2._wp*dx(j))) * ( &
                                        q_prim_vf(momxb)%sf(j+1,k+q,l)/rho_sf_small(1) - &
                                        q_prim_vf(momxb)%sf(j-1,k+q,l)/rho_sf_small(-1))
                                    dvel_small(2) = (1/(2._wp*dx(j))) * ( &
                                        q_prim_vf(momxb+1)%sf(j+1,k+q,l)/rho_sf_small(1) - &
                                        q_prim_vf(momxb+1)%sf(j-1,k+q,l)/rho_sf_small(-1))

                                    if (q > -2) then
                                        vflux_L_arr(1) = vflux_L_arr(1) + coeff_L(q)*(dvel_small(2))
                                        vflux_L_arr(3) = vflux_L_arr(3) + coeff_L(q)*(-2._wp*dvel_small(1))/3._wp
                                    end if
                                    if (q < 3) then
                                        vflux_R_arr(1) = vflux_R_arr(1) + coeff_R(q)*(dvel_small(2))
                                        vflux_R_arr(3) = vflux_R_arr(3) + coeff_R(q)*(-2._wp*dvel_small(1))/3._wp
                                    end if

                                    !y-direction contributions
                                    !$acc loop seq 
                                    do i = -1, 1
                                        rho_L = 0._wp
                                        !$acc loop seq 
                                        do r = 1, num_fluids
                                            rho_L = rho_L + q_prim_vf(r)%sf(j,k+i+q,l)
                                        end do
                                        rho_sf_small(i) = rho_L
                                    end do

                                    dvel_small(1) = (1/(2._wp*dy(k))) * ( &
                                        q_prim_vf(momxb)%sf(j,k+1+q,l)/rho_sf_small(1) - &
                                        q_prim_vf(momxb)%sf(j,k-1+q,l)/rho_sf_small(-1))
                                    dvel_small(2) = (1/(2._wp*dy(k))) * ( &
                                        q_prim_vf(momxb+1)%sf(j,k+1+q,l)/rho_sf_small(1) - &
                                        q_prim_vf(momxb+1)%sf(j,k-1+q,l)/rho_sf_small(-1))
                                    dvel_small(3) = (1/(2._wp*dy(k))) * ( &
                                        q_prim_vf(momxb+2)%sf(j,k+1+q,l)/rho_sf_small(1) - &
                                        q_prim_vf(momxb+2)%sf(j,k-1+q,l)/rho_sf_small(-1))

                                    if (q > -2) then
                                        vflux_L_arr(1) = vflux_L_arr(1) + coeff_L(q)*(dvel_small(1))
                                        vflux_L_arr(2) = vflux_L_arr(2) + coeff_L(q)*(dvel_small(3))
                                        vflux_L_arr(3) = vflux_L_arr(3) + coeff_L(q)*(4._wp*dvel_small(2))/3._wp
                                    end if
                                    if (q < 3) then
                                        vflux_R_arr(1) = vflux_R_arr(1) + coeff_R(q)*(dvel_small(1))
                                        vflux_R_arr(2) = vflux_R_arr(2) + coeff_R(q)*(dvel_small(3))
                                        vflux_R_arr(3) = vflux_R_arr(3) + coeff_R(q)*(4._wp*dvel_small(2))/3._wp
                                    end if

                                    !z-direction contributions
                                    !$acc loop seq 
                                    do i = -1, 1
                                        rho_L = 0._wp
                                        !$acc loop seq 
                                        do r = 1, num_fluids
                                            rho_L = rho_L + q_prim_vf(r)%sf(j,k+q,l+i)
                                        end do
                                        rho_sf_small(i) = rho_L
                                    end do

                                    dvel_small(2) = (1/(2._wp*dz(l))) * ( &
                                        q_prim_vf(momxb+1)%sf(j,k+q,l+1)/rho_sf_small(1) - &
                                        q_prim_vf(momxb+1)%sf(j,k+q,l-1)/rho_sf_small(-1))
                                    dvel_small(3) = (1/(2._wp*dz(l))) * ( &
                                        q_prim_vf(momxb+2)%sf(j,k+q,l+1)/rho_sf_small(1) - &
                                        q_prim_vf(momxb+2)%sf(j,k+q,l-1)/rho_sf_small(-1))
                                    if (q > -2) then
                                        vflux_L_arr(2) = vflux_L_arr(2) + coeff_L(q)*(dvel_small(2))
                                        vflux_L_arr(3) = vflux_L_arr(3) + coeff_L(q)*(-2._wp*dvel_small(3))/3._wp
                                    end if
                                    if (q < 3) then
                                        vflux_R_arr(2) = vflux_R_arr(2) + coeff_R(q)*(dvel_small(2))
                                        vflux_R_arr(3) = vflux_R_arr(3) + coeff_R(q)*(-2._wp*dvel_small(3))/3._wp
                                    end if
                                end do
                            end if

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
                            end do 

                            if(num_fluids > 1) then  
                                !$acc loop seq 
                                do i = 1,num_fluids - 1
                                    alpha_L(i) =  (1._wp/60._wp) * (-3._wp * q_prim_vf(E_idx+i)%sf(j, k-1, l) + &
                                                27._wp * q_prim_vf(E_idx+i)%sf(j, k, l) + &
                                                47._wp * q_prim_vf(E_idx+i)%sf(j, k+1, l) -   &
                                                13._wp * q_prim_vf(E_idx+i)%sf(j, k+2, l) + &
                                                2._wp * q_prim_vf(E_idx+i)%sf(j, k+3, l))                          
                                    alpha_R(i) = (1._wp/60._wp) * (-3._wp * q_prim_vf(E_idx+i)%sf(j, k+2, l) + &
                                                        27._wp * q_prim_vf(E_idx+i)%sf(j, k+1, l) + &
                                                        47._wp * q_prim_vf(E_idx+i)%sf(j, k, l) -   &
                                                        13._wp * q_prim_vf(E_idx+i)%sf(j, k-1, l) + &
                                                        2._wp * q_prim_vf(E_idx+i)%sf(j, k-2, l))
                                end do
                                alpha_L(num_fluids) = 1._wp - sum(alpha_L(1:num_fluids-1))
                                alpha_R(num_fluids) = 1._wp - sum(alpha_R(1:num_fluids-1))
                            else 
                                alpha_L(1) = 1._wp
                                alpha_R(1) = 1._wp 
                            end if

                            rho_L = sum(alpha_rho_L)
                            gamma_L = sum(alpha_L*gammas)
                            pi_inf_L = sum(alpha_L*pi_infs)

                            rho_R = sum(alpha_rho_R)
                            gamma_R = sum(alpha_R*gammas)
                            pi_inf_R = sum(alpha_R*pi_infs)

                            !$acc loop seq 
                            do i = 1, num_dims
                                vel_L(i) = (1._wp/60._wp) * (-3._wp * q_prim_vf(momxb+i-1)%sf(j,k-1,l) + &
                                                27._wp * q_prim_vf(momxb+i-1)%sf(j,k,l) + &
                                                47._wp * q_prim_vf(momxb+i-1)%sf(j,k+1,l) -   &
                                                13._wp * q_prim_vf(momxb+i-1)%sf(j,k+2,l) + &
                                                2._wp * q_prim_vf(momxb+i-1)%sf(j,k+3,l)) / rho_L
                            end do 

                            !$acc loop seq 
                            do i = 1, num_dims
                                vel_R(i) = (1._wp/60._wp) * (-3._wp * q_prim_vf(momxb+i-1)%sf(j,k+2,l) + &
                                                27._wp * q_prim_vf(momxb+i-1)%sf(j,k+1,l) + &
                                                47._wp * q_prim_vf(momxb+i-1)%sf(j,k,l) -   &
                                                13._wp * q_prim_vf(momxb+i-1)%sf(j,k-1,l) + &
                                                2._wp * q_prim_vf(momxb+i-1)%sf(j,k-2,l)) / rho_R
                            end do

                            if(viscous) then
       
                                if(num_fluids > 1) then 
                                    mu_L = 0._wp
                                    mu_R = 0._wp
                                    !$acc loop seq 
                                    do i = 1, num_fluids - 1
                                        mu_L = (1._wp/60._wp) * (-3._wp * q_prim_vf(E_idx+i)%sf(j, k-1, l) + &
                                                27._wp * q_prim_vf(E_idx+i)%sf(j, k, l) + &
                                                47._wp * q_prim_vf(E_idx+i)%sf(j, k+1, l) -   &
                                                13._wp * q_prim_vf(E_idx+i)%sf(j, k+2, l) + &
                                                2._wp * q_prim_vf(E_idx+i)%sf(j, k+3, l)) / Res(1, i) + mu_L

                                        mu_R = (1._wp/60._wp) * (-3._wp * q_prim_vf(E_idx+i)%sf(j, k+2, l) + &
                                                27._wp * q_prim_vf(E_idx+i)%sf(j, k+1, l) + &
                                                47._wp * q_prim_vf(E_idx+i)%sf(j, k, l) -   &
                                                13._wp * q_prim_vf(E_idx+i)%sf(j, k-1, l) + &
                                                2._wp * q_prim_vf(E_idx+i)%sf(j, k-2, l)) / Res(1, i) + mu_R
                                    end do

                                    mu_L = 1._wp / Res(1, num_fluids) + mu_L
                                    mu_R = 1._wp / Res(1, num_fluids) + mu_R

                                    !$acc loop seq 
                                    do i = 1, num_fluids - 1
                                        mu_L = -(1._wp/60._wp) * (-3._wp * q_prim_vf(E_idx+i)%sf(j, k-1, l) + &
                                                27._wp * q_prim_vf(E_idx+i)%sf(j, k, l) + &
                                                47._wp * q_prim_vf(E_idx+i)%sf(j, k+1, l) -   &
                                                13._wp * q_prim_vf(E_idx+i)%sf(j, k+2, l) + &
                                                2._wp * q_prim_vf(E_idx+i)%sf(j, k+3, l)) / Res(1, num_fluids) + mu_L

                                        mu_R = -(1._wp/60._wp) * (-3._wp * q_prim_vf(E_idx+i)%sf(j, k+2, l) + &
                                                27._wp * q_prim_vf(E_idx+i)%sf(j, k+1, l) + &
                                                47._wp * q_prim_vf(E_idx+i)%sf(j, k, l) -   &
                                                13._wp * q_prim_vf(E_idx+i)%sf(j, k-1, l) + &
                                                2._wp * q_prim_vf(E_idx+i)%sf(j, k-2, l)) / Res(1, num_fluids) + mu_R
                                    end do
                                else
                                    mu_L = 1._wp / Res(1, 1) 
                                    mu_R = 1._wp / Res(1, 1) 
                                end if

                                @:ATOMIC_ADD(rhs_vf(momxb)%sf(j,k+1,l),(- 0.5_wp*mu_L*vflux_L_arr(1)*(1._wp/dy(k+1))*dt))
                                !!$omp atomic update
                                !!$acc atomic
                                !rhs_vf(momxb)%sf(j,k+1,l) = rhs_vf(momxb)%sf(j,k+1,l) - 0.5_wp*mu_L*vflux_L_arr(1)*(1._wp/dy(k+1))
                                @:ATOMIC_ADD(rhs_vf(E_idx)%sf(j,k+1,l),(- 0.5_wp*mu_L*vflux_L_arr(1)*vel_L(1)*(1._wp/dy(k+1))*dt))
                                !!$omp atomic update
                                !!$acc atomic
                                !rhs_vf(E_idx)%sf(j,k+1,l) = rhs_vf(E_idx)%sf(j,k+1,l) - 0.5_wp*mu_L*vflux_L_arr(1)*vel_L(1)*(1._wp/dy(k+1)) 

                                @:ATOMIC_ADD(rhs_vf(momxb)%sf(j,k,l),(0.5_wp*mu_L*vflux_L_arr(1)*(1._wp/dy(k))*dt))
                                !!$omp atomic update
                                !!$acc atomic
                                !rhs_vf(momxb)%sf(j,k,l) = rhs_vf(momxb)%sf(j,k,l) + 0.5_wp*mu_L*vflux_L_arr(1)*(1._wp/dy(k))
                                @:ATOMIC_ADD(rhs_vf(E_idx)%sf(j,k,l),(0.5_wp*mu_L*vflux_L_arr(1)*vel_L(1)*(1._wp/dy(k))*dt))
                                !!$omp atomic update
                                !!$acc atomic
                                !rhs_vf(E_idx)%sf(j,k,l) = rhs_vf(E_idx)%sf(j,k,l) + 0.5_wp*mu_L*vflux_L_arr(1)*vel_L(1)*(1._wp/dy(k))
                            
                                @:ATOMIC_ADD(rhs_vf(momxb)%sf(j,k+1,l),(- 0.5_wp*mu_R*vflux_R_arr(1)*(1._wp/dy(k+1))*dt))
                                !!$omp atomic update
                                !!$acc atomic
                                !rhs_vf(momxb)%sf(j,k+1,l) = rhs_vf(momxb)%sf(j,k+1,l) - 0.5_wp*mu_R*vflux_R_arr(1)*(1._wp/dy(k+1))
                                @:ATOMIC_ADD(rhs_vf(E_idx)%sf(j,k+1,l),(- 0.5_wp*mu_R*vflux_R_arr(1)*vel_R(1)*(1._wp/dy(k+1))*dt))
                                !!$omp atomic update
                                !!$acc atomic
                                !rhs_vf(E_idx)%sf(j,k+1,l) = rhs_vf(E_idx)%sf(j,k+1,l) - 0.5_wp*mu_R*vflux_R_arr(1)*vel_R(1)*(1._wp/dy(k+1))
                            
                                @:ATOMIC_ADD(rhs_vf(momxb)%sf(j,k,l),(0.5_wp*mu_R*vflux_R_arr(1)*(1._wp/dy(k))*dt))
                                !!$omp atomic update
                                !!$acc atomic
                                !rhs_vf(momxb)%sf(j,k,l) = rhs_vf(momxb)%sf(j,k,l) + 0.5_wp*mu_R*vflux_R_arr(1)*(1._wp/dy(k))
                                @:ATOMIC_ADD(rhs_vf(E_idx)%sf(j,k,l),(0.5_wp*mu_R*vflux_R_arr(1)*vel_R(1)*(1._wp/dy(k))*dt))
                                !!$omp atomic update
                                !!$acc atomic
                                !rhs_vf(E_idx)%sf(j,k,l) = rhs_vf(E_idx)%sf(j,k,l) + 0.5_wp*mu_R*vflux_R_arr(1)*vel_R(1)*(1._wp/dy(k))


                                @:ATOMIC_ADD(rhs_vf(momxb+2)%sf(j,k+1,l),(- 0.5_wp*mu_L*vflux_L_arr(2)*(1._wp/dy(k+1))*dt))
                                !!$omp atomic update
                                !!$acc atomic
                                !rhs_vf(momxb+2)%sf(j,k+1,l) = rhs_vf(momxb+2)%sf(j,k+1,l) - 0.5_wp*mu_L*vflux_L_arr(2)*(1._wp/dy(k+1))
                                @:ATOMIC_ADD(rhs_vf(E_idx)%sf(j,k+1,l),(- 0.5_wp*mu_L*vflux_L_arr(2)*vel_L(3)*(1._wp/dy(k+1))*dt))
                                !!$omp atomic update
                                !!$acc atomic
                                !rhs_vf(E_idx)%sf(j,k+1,l) = rhs_vf(E_idx)%sf(j,k+1,l) - 0.5_wp*mu_L*vflux_L_arr(2)*vel_L(3)*(1._wp/dy(k+1))
                            
                                @:ATOMIC_ADD(rhs_vf(momxb+2)%sf(j,k,l),(0.5_wp*mu_L*vflux_L_arr(2)*(1._wp/dy(k))*dt))
                                !!$omp atomic update
                                !!$acc atomic
                                !rhs_vf(momxb+2)%sf(j,k,l) = rhs_vf(momxb+2)%sf(j,k,l) + 0.5_wp*mu_L*vflux_L_arr(2)*(1._wp/dy(k))
                                @:ATOMIC_ADD(rhs_vf(E_idx)%sf(j,k,l),(0.5_wp*mu_L*vflux_L_arr(2)*vel_L(3)*(1._wp/dy(k))*dt))
                                !!$omp atomic update
                                !!$acc atomic
                                !rhs_vf(E_idx)%sf(j,k,l) = rhs_vf(E_idx)%sf(j,k,l) + 0.5_wp*mu_L*vflux_L_arr(2)*vel_L(3)*(1._wp/dy(k))
                            
                                @:ATOMIC_ADD(rhs_vf(momxb+2)%sf(j,k+1,l),(- 0.5_wp*mu_R*vflux_R_arr(2)*(1._wp/dy(k+1))*dt))
                                !!$omp atomic update
                                !!$acc atomic
                                !rhs_vf(momxb+2)%sf(j,k+1,l) = rhs_vf(momxb+2)%sf(j,k+1,l) - 0.5_wp*mu_R*vflux_R_arr(2)*(1._wp/dy(k+1))
                                @:ATOMIC_ADD(rhs_vf(E_idx)%sf(j,k+1,l),(- 0.5_wp*mu_R*vflux_R_arr(2)*vel_R(3)*(1._wp/dy(k+1))*dt))
                                !!$omp atomic update
                                !!$acc atomic
                                !rhs_vf(E_idx)%sf(j,k+1,l) = rhs_vf(E_idx)%sf(j,k+1,l) - 0.5_wp*mu_R*vflux_R_arr(2)*vel_R(3)*(1._wp/dy(k+1))
                            
                                @:ATOMIC_ADD(rhs_vf(momxb+2)%sf(j,k,l),(0.5_wp*mu_R*vflux_R_arr(2)*(1._wp/dy(k))*dt))
                                !!$omp atomic update
                                !!$acc atomic
                                !rhs_vf(momxb+2)%sf(j,k,l) = rhs_vf(momxb+2)%sf(j,k,l) + 0.5_wp*mu_R*vflux_R_arr(2)*(1._wp/dy(k))
                                @:ATOMIC_ADD(rhs_vf(E_idx)%sf(j,k,l),(0.5_wp*mu_R*vflux_R_arr(2)*vel_R(3)*(1._wp/dy(k))*dt))
                                !!$omp atomic update
                                !!$acc atomic
                                !rhs_vf(E_idx)%sf(j,k,l) = rhs_vf(E_idx)%sf(j,k,l) + 0.5_wp*mu_R*vflux_R_arr(2)*vel_R(3)*(1._wp/dy(k))


                                @:ATOMIC_ADD(rhs_vf(momxb+1)%sf(j,k+1,l),(- 0.5_wp*mu_L*vflux_L_arr(3)*(1._wp/dy(k+1))*dt))
                                !!$omp atomic update
                                !!$acc atomic
                                !rhs_vf(momxb+1)%sf(j,k+1,l) = rhs_vf(momxb+1)%sf(j,k+1,l) - 0.5_wp*mu_L*vflux_L_arr(3)*(1._wp/dy(k+1))
                                @:ATOMIC_ADD(rhs_vf(E_idx)%sf(j,k+1,l),(- 0.5_wp*mu_L*vflux_L_arr(3)*vel_L(2)*(1._wp/dy(k+1))*dt))
                                !!$omp atomic update
                                !!$acc atomic
                                !rhs_vf(E_idx)%sf(j,k+1,l) = rhs_vf(E_idx)%sf(j,k+1,l) - 0.5_wp*mu_L*vflux_L_arr(3)*vel_L(2)*(1._wp/dy(k+1))
                        
                                @:ATOMIC_ADD(rhs_vf(momxb+1)%sf(j,k,l),(0.5_wp*mu_L*vflux_L_arr(3)*(1._wp/dy(k))*dt))
                                !!$omp atomic update
                                !!$acc atomic
                                !rhs_vf(momxb+1)%sf(j,k,l) = rhs_vf(momxb+1)%sf(j,k,l) + 0.5_wp*mu_L*vflux_L_arr(3)*(1._wp/dy(k))
                                @:ATOMIC_ADD(rhs_vf(E_idx)%sf(j,k,l),(0.5_wp*mu_L*vflux_L_arr(3)*vel_L(2)*(1._wp/dy(k))*dt))
                                !!$omp atomic update
                                !!$acc atomic
                                !rhs_vf(E_idx)%sf(j,k,l) = rhs_vf(E_idx)%sf(j,k,l) + 0.5_wp*mu_L*vflux_L_arr(3)*vel_L(2)*(1._wp/dy(k))                    
                             
                                @:ATOMIC_ADD(rhs_vf(momxb+1)%sf(j,k+1,l),(- 0.5_wp*mu_R*vflux_R_arr(3)*(1._wp/dy(k+1))*dt))
                                !!$omp atomic update
                                !!$acc atomic
                                !rhs_vf(momxb+1)%sf(j,k+1,l) = rhs_vf(momxb+1)%sf(j,k+1,l) - 0.5_wp*mu_R*vflux_R_arr(3)*(1._wp/dy(k+1))
                                @:ATOMIC_ADD(rhs_vf(E_idx)%sf(j,k+1,l),(- 0.5_wp*mu_R*vflux_R_arr(3)*vel_R(2)*(1._wp/dy(k+1))*dt))
                                !!$omp atomic update
                                !!$acc atomic
                                !rhs_vf(E_idx)%sf(j,k+1,l) = rhs_vf(E_idx)%sf(j,k+1,l) - 0.5_wp*mu_R*vflux_R_arr(3)*vel_R(2)*(1._wp/dy(k+1))
                             
                                @:ATOMIC_ADD(rhs_vf(momxb+1)%sf(j,k,l),(0.5_wp*mu_R*vflux_R_arr(3)*(1._wp/dy(k))*dt))
                                !!$omp atomic update
                                !!$acc atomic
                                !rhs_vf(momxb+1)%sf(j,k,l) = rhs_vf(momxb+1)%sf(j,k,l) + 0.5_wp*mu_R*vflux_R_arr(3)*(1._wp/dy(k))
                                @:ATOMIC_ADD(rhs_vf(E_idx)%sf(j,k,l),(0.5_wp*mu_R*vflux_R_arr(3)*vel_R(2)*(1._wp/dy(k))*dt))
                                !!$omp atomic update
                                !!$acc atomic
                                !rhs_vf(E_idx)%sf(j,k,l) = rhs_vf(E_idx)%sf(j,k,l) + 0.5_wp*mu_R*vflux_R_arr(3)*vel_R(2)*(1._wp/dy(k))

                                
                            endif

                            F_L = (1._wp/60._wp) * (-3._wp * jac(j, k-1, l) + &
                                                27._wp * jac(j, k, l) + &
                                                47._wp * jac(j, k+1, l) -   &
                                                13._wp * jac(j, k+2, l) + &
                                                2._wp * jac(j, k+3, l))

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

                            pres_L = (E_L - pi_inf_L - 0.5_wp*rho_L*(vel_L(1)**2._wp + vel_L(2)**2._wp + vel_L(3)**2._wp ))/gamma_L                    

                            pres_R = (E_R - pi_inf_R - 0.5_wp*rho_R*(vel_R(1)**2._wp + vel_R(2)**2._wp + vel_R(3)**2._wp ))/gamma_R 

                            a_L = sqrt((pres_L*(1._wp/gamma_L + 1._wp) + pi_inf_L / gamma_L) / rho_L)
                            a_R = sqrt((pres_R*(1._wp/gamma_R + 1._wp) + pi_inf_R / gamma_R) / rho_R)

                            cfl = (max(sqrt(vel_L(1)**2._wp + vel_L(2)**2._wp + vel_L(3)**2._wp),sqrt(vel_R(1)**2._wp + vel_R(2)**2._wp + vel_R(3)**2._wp)) + max(a_L,a_R))

                            !$acc loop seq
                            do i = 1, num_fluids
                                @:ATOMIC_ADD(rhs_vf(i)%sf(j,k+1,l),(0.5_wp * (alpha_rho_L(i) * &
                                    vel_L(2))*(1._wp/dy(k+1)) - &
                                    0.5_wp*cfl * (alpha_rho_L(i))*(1._wp/dy(k+1)))*dt)
                                ! !$omp atomic update
                                ! !$acc atomic
                                ! rhs_vf(i)%sf(j,k+1,l) = rhs_vf(i)%sf(j,k+1,l) + &
                                !     (0.5_wp * (alpha_rho_L(i) * &
                                !     vel_L(2))*(1._wp/dy(k+1)) - &
                                !     0.5_wp*cfl * (alpha_rho_L(i))*(1._wp/dy(k+1)))

                                @:ATOMIC_ADD(rhs_vf(i)%sf(j,k,l),(-0.5_wp * (alpha_rho_L(i) * &
                                    vel_L(2))*(1._wp/dy(k)) + &
                                    0.5_wp*cfl * (alpha_rho_L(i))*(1._wp/dy(k)))*dt)
                                ! !$omp atomic update
                                ! !$acc atomic
                                ! rhs_vf(i)%sf(j,k,l) = rhs_vf(i)%sf(j,k,l) - &
                                !     (0.5_wp * (alpha_rho_L(i) * &
                                !     vel_L(2))*(1._wp/dy(k)) - &
                                !     0.5_wp*cfl * (alpha_rho_L(i))*(1._wp/dy(k))) 
                            end do

                            if(num_fluids > 1) then 
                                !$acc loop seq
                                do i = 1, num_fluids - 1
                                    @:ATOMIC_ADD(rhs_vf(advxb+i-1)%sf(j,k+1,l),(0.5_wp * (alpha_L(i) * &
                                        vel_L(2))*(1._wp/dy(k+1)) - &
                                        0.5_wp*cfl*(alpha_L(i))*(1._wp/dy(k+1)))*dt)
                                    ! !$omp atomic update
                                    ! !$acc atomic
                                    ! rhs_vf(advxb+i-1)%sf(j,k+1,l) = rhs_vf(advxb+i-1)%sf(j,k+1,l) + &
                                    !     (0.5_wp * (alpha_L(i) * &
                                    !     vel_L(2))*(1._wp/dy(k+1)) - &
                                    !     0.5_wp*cfl*(alpha_L(i))*(1._wp/dy(k+1)))

                                    @:ATOMIC_ADD(rhs_vf(advxb+i-1)%sf(j,k+1,l),(-0.5_wp * q_prim_vf(advxb+i-1)%sf(j,k+1,l)  * vel_L(2)*(1._wp/dy(k+1)))*dt)
                                    ! !$omp atomic update
                                    ! !$acc atomic
                                    ! rhs_vf(advxb+i-1)%sf(j,k+1,l) = rhs_vf(advxb+i-1)%sf(j,k+1,l) &
                                    ! - (0.5_wp * q_prim_vf(advxb+i-1)%sf(j,k+1,l)  * vel_L(2)*(1._wp/dy(k+1)))

                                    @:ATOMIC_ADD(rhs_vf(advxb+i-1)%sf(j,k,l),(-0.5_wp * (alpha_L(i) * &
                                        vel_L(2))*(1._wp/dy(k)) + &
                                        0.5_wp*cfl*(alpha_L(i))*(1._wp/dy(k)))*dt)
                                    ! !$omp atomic update
                                    ! !$acc atomic
                                    ! rhs_vf(advxb+i-1)%sf(j,k,l) = rhs_vf(advxb+i-1)%sf(j,k,l) - &
                                    !     (0.5_wp * (alpha_L(i) * &
                                    !     vel_L(2))*(1._wp/dy(k)) - &
                                    !     0.5_wp*cfl*(alpha_L(i))*(1._wp/dy(k)))

                                    @:ATOMIC_ADD(rhs_vf(advxb+i-1)%sf(j,k,l),(0.5_wp * q_prim_vf(advxb+i-1)%sf(j,k,l)  * vel_L(2)*(1._wp/dy(k)))*dt)
                                    ! !$omp atomic update
                                    ! !$acc atomic
                                    ! rhs_vf(advxb+i-1)%sf(j,k,l) = rhs_vf(advxb+i-1)%sf(j,k,l) &
                                    ! + (0.5_wp * q_prim_vf(advxb+i-1)%sf(j,k,l)  * vel_L(2)*(1._wp/dy(k)))
                                end do
                            end if

                            @:ATOMIC_ADD(rhs_vf(momxb+1)%sf(j,k+1,l),(0.5_wp* (rho_L * (vel_L(2))**2.0 + &
                                 pres_L+F_L)*(1._wp/dy(k+1)) - &
                                 0.5_wp*cfl * (rho_L*vel_L(2))*(1._wp/dy(k+1)))*dt)
                            ! !$omp atomic update
                            ! !$acc atomic
                            ! rhs_vf(momxb+1)%sf(j,k+1,l) = rhs_vf(momxb+1)%sf(j,k+1,l) + &
                            !      (0.5_wp* (rho_L * (vel_L(2))**2.0 + &
                            !      pres_L+F_L)*(1._wp/dy(k+1)) - &
                            !      0.5_wp*cfl * (rho_L*vel_L(2))*(1._wp/dy(k+1)))

                            @:ATOMIC_ADD(rhs_vf(momxb)%sf(j, k+1, l),(0.5_wp * rho_L * vel_L(1)*vel_L(2)*(1._wp/dy(k+1)) - &
                                0.5_wp*cfl * (rho_L*vel_L(1))*(1._wp/dy(k+1)))*dt)
                            ! !$omp atomic update
                            ! !$acc atomic
                            ! rhs_vf(momxb)%sf(j, k+1, l) =  rhs_vf(momxb)%sf(j, k+1, l) + &
                            !     (0.5_wp * rho_L * vel_L(1)*vel_L(2)*(1._wp/dy(k+1)) - &
                            !     0.5_wp*cfl * (rho_L*vel_L(1))*(1._wp/dy(k+1)))

                            @:ATOMIC_ADD(rhs_vf(momxb+2)%sf(j, k+1, l),(0.5_wp * rho_L * vel_L(3)*vel_L(2)*(1._wp/dy(k+1)) - &
                                0.5_wp*cfl * (rho_L*vel_L(3))*(1._wp/dy(k+1)))*dt)
                            ! !$omp atomic update
                            ! !$acc atomic
                            ! rhs_vf(momxb+2)%sf(j, k+1, l) =  rhs_vf(momxb+2)%sf(j, k+1, l) + &
                            !     (0.5_wp * rho_L * vel_L(3)*vel_L(2)*(1._wp/dy(k+1)) - &
                            !     0.5_wp*cfl * (rho_L*vel_L(3))*(1._wp/dy(k+1)))

                            @:ATOMIC_ADD(rhs_vf(E_idx)%sf(j, k+1, l),(0.5_wp * (vel_L(2) * (E_L + &
                                pres_L+F_L) )*(1._wp/dy(k+1)) - &
                                0.5_wp*cfl * (E_L)*(1._wp/dy(k+1)))*dt)
                            ! !$omp atomic update
                            ! !$acc atomic
                            ! rhs_vf(E_idx)%sf(j, k+1, l) = rhs_vf(E_idx)%sf(j, k+1, l) +  &
                            !     (0.5_wp * (vel_L(2) * (E_L + &
                            !     pres_L+F_L) )*(1._wp/dy(k+1)) - &
                            !     0.5_wp*cfl * (E_L)*(1._wp/dy(k+1)))

                            @:ATOMIC_ADD(rhs_vf(momxb+1)%sf(j,k,l),(-0.5_wp* (rho_L * (vel_L(2))**2.0 + &
                                 pres_L+F_L)*(1._wp/dy(k)) + &
                                 0.5_wp*cfl * (rho_L*vel_L(2))*(1._wp/dy(k)))*dt)
                            ! !$omp atomic update
                            ! !$acc atomic
                            ! rhs_vf(momxb+1)%sf(j,k,l) = rhs_vf(momxb+1)%sf(j,k,l) - &
                            !      (0.5_wp* (rho_L * (vel_L(2))**2.0 + &
                            !      pres_L+F_L)*(1._wp/dy(k)) - &
                            !      0.5_wp*cfl * (rho_L*vel_L(2))*(1._wp/dy(k)))

                            @:ATOMIC_ADD(rhs_vf(momxb)%sf(j, k, l),(-0.5_wp * rho_L * vel_L(1)*vel_L(2)*(1._wp/dy(k)) + &
                                0.5_wp*cfl * (rho_L*vel_L(1))*(1._wp/dy(k)))*dt)
                            ! !$omp atomic update
                            ! !$acc atomic
                            ! rhs_vf(momxb)%sf(j, k, l) =  rhs_vf(momxb)%sf(j, k, l) - &
                            !     (0.5_wp * rho_L * vel_L(1)*vel_L(2)*(1._wp/dy(k)) - &
                            !     0.5_wp*cfl * (rho_L*vel_L(1))*(1._wp/dy(k)))

                            @:ATOMIC_ADD(rhs_vf(momxb+2)%sf(j, k, l),(-0.5_wp * rho_L * vel_L(3)*vel_L(2)*(1._wp/dy(k)) + &
                                0.5_wp*cfl * (rho_L*vel_L(3))*(1._wp/dy(k)))*dt)
                            ! !$omp atomic update
                            ! !$acc atomic
                            ! rhs_vf(momxb+2)%sf(j, k, l) =  rhs_vf(momxb+2)%sf(j, k, l) - &
                            !     (0.5_wp * rho_L * vel_L(3)*vel_L(2)*(1._wp/dy(k)) - &
                            !     0.5_wp*cfl * (rho_L*vel_L(3))*(1._wp/dy(k)))

                            @:ATOMIC_ADD(rhs_vf(E_idx)%sf(j, k, l),(-0.5_wp * (vel_L(2) * (E_L + &
                                pres_L+F_L) )*(1._wp/dy(k)) + &
                                0.5_wp*cfl * (E_L)*(1._wp/dy(k)))*dt)
                            ! !$omp atomic update
                            ! !$acc atomic
                            ! rhs_vf(E_idx)%sf(j, k, l) = rhs_vf(E_idx)%sf(j, k, l) -  &
                            !     (0.5_wp * (vel_L(2) * (E_L + &
                            !     pres_L+F_L) )*(1._wp/dy(k)) - &
                            !     0.5_wp*cfl * (E_L)*(1._wp/dy(k)))


                            F_L = (1._wp/60._wp) * (-3._wp * jac(j, k+2, l) + &
                                                27._wp * jac(j, k+1, l) + &
                                                47._wp * jac(j, k, l) -   &
                                                13._wp * jac(j, k-1, l) + &
                                                2._wp * jac(j, k-2, l))

                            !$acc loop seq
                            do i = 1, num_fluids
                                @:ATOMIC_ADD(rhs_vf(i)%sf(j,k+1,l),(0.5_wp * (alpha_rho_R(i) * &
                                    vel_R(2))*(1._wp/dy(k+1)) + &
                                    0.5_wp*cfl * (alpha_rho_R(i))*(1._wp/dy(k+1)))*dt)
                                ! !$omp atomic update
                                ! !$acc atomic
                                ! rhs_vf(i)%sf(j,k+1,l) = rhs_vf(i)%sf(j,k+1,l) + &
                                !     (0.5_wp * (alpha_rho_R(i) * &
                                !     vel_R(2))*(1._wp/dy(k+1)) + &
                                !     0.5_wp*cfl * (alpha_rho_R(i))*(1._wp/dy(k+1)))

                                @:ATOMIC_ADD(rhs_vf(i)%sf(j,k,l),(-0.5_wp * (alpha_rho_R(i) * &
                                    vel_R(2))*(1._wp/dy(k)) - &
                                    0.5_wp*cfl * (alpha_rho_R(i))*(1._wp/dy(k)))*dt)
                                ! !$omp atomic update
                                ! !$acc atomic
                                ! rhs_vf(i)%sf(j,k,l) = rhs_vf(i)%sf(j,k,l) - &
                                !     (0.5_wp * (alpha_rho_R(i) * &
                                !     vel_R(2))*(1._wp/dy(k)) + &
                                !     0.5_wp*cfl * (alpha_rho_R(i))*(1._wp/dy(k))) 
                            end do

                            if(num_fluids > 1) then 
                                !$acc loop seq
                                do i = 1, num_fluids - 1
                                    @:ATOMIC_ADD(rhs_vf(advxb+i-1)%sf(j,k+1,l),(0.5_wp * (alpha_R(i) * &
                                        vel_R(2))*(1._wp/dy(k+1)) + &
                                        0.5_wp*cfl*(alpha_R(i))*(1._wp/dy(k+1)))*dt)
                                    ! !$omp atomic update
                                    ! !$acc atomic
                                    ! rhs_vf(advxb+i-1)%sf(j,k+1,l) = rhs_vf(advxb+i-1)%sf(j,k+1,l) + &
                                    !     (0.5_wp * (alpha_R(i) * &
                                    !     vel_R(2))*(1._wp/dy(k+1)) + &
                                    !     0.5_wp*cfl*(alpha_R(i))*(1._wp/dy(k+1)))

                                    @:ATOMIC_ADD(rhs_vf(advxb+i-1)%sf(j,k+1,l),(-0.5_wp * q_prim_vf(advxb+i-1)%sf(j,k+1,l)  * vel_R(2)*(1._wp/dy(k+1)))*dt)
                                    ! !$omp atomic update
                                    ! !$acc atomic
                                    ! rhs_vf(advxb+i-1)%sf(j,k+1,l) = rhs_vf(advxb+i-1)%sf(j,k+1,l) &
                                    ! - (0.5_wp * q_prim_vf(advxb+i-1)%sf(j,k+1,l)  * vel_R(2)*(1._wp/dy(k+1)))

                                    @:ATOMIC_ADD(rhs_vf(advxb+i-1)%sf(j,k,l),(-0.5_wp * (alpha_R(i) * &
                                        vel_R(2))*(1._wp/dy(k)) - &
                                        0.5_wp*cfl*(alpha_R(i))*(1._wp/dy(k)))*dt)
                                    ! !$omp atomic update
                                    ! !$acc atomic
                                    ! rhs_vf(advxb+i-1)%sf(j,k,l) = rhs_vf(advxb+i-1)%sf(j,k,l) - &
                                    !     (0.5_wp * (alpha_R(i) * &
                                    !     vel_R(2))*(1._wp/dy(k)) + &
                                    !     0.5_wp*cfl*(alpha_R(i))*(1._wp/dy(k)))

                                    @:ATOMIC_ADD(rhs_vf(advxb+i-1)%sf(j,k,l),(0.5_wp * q_prim_vf(advxb+i-1)%sf(j,k,l)  * vel_R(2)*(1._wp/dy(k)))*dt)
                                    ! !$omp atomic update
                                    ! !$acc atomic
                                    ! rhs_vf(advxb+i-1)%sf(j,k,l) = rhs_vf(advxb+i-1)%sf(j,k,l) &
                                    ! + (0.5_wp * q_prim_vf(advxb+i-1)%sf(j,k,l)  * vel_R(2)*(1._wp/dy(k)))
                                end do
                            end if

                            @:ATOMIC_ADD(rhs_vf(momxb+1)%sf(j,k+1,l),(0.5_wp* (rho_R * (vel_R(2))**2.0 + &
                                 pres_R+F_L)*(1._wp/dy(k+1)) + &
                                 0.5_wp*cfl * (rho_R*vel_R(2))*(1._wp/dy(k+1)))*dt)
                            ! !$omp atomic update
                            ! !$acc atomic
                            ! rhs_vf(momxb+1)%sf(j,k+1,l) = rhs_vf(momxb+1)%sf(j,k+1,l) + &
                            !      (0.5_wp* (rho_R * (vel_R(2))**2.0 + &
                            !      pres_R+F_L)*(1._wp/dy(k+1)) + &
                            !      0.5_wp*cfl * (rho_R*vel_R(2))*(1._wp/dy(k+1)))

                            @:ATOMIC_ADD(rhs_vf(momxb)%sf(j, k+1, l),(0.5_wp * rho_R * vel_R(2)*vel_R(1)*(1._wp/dy(k+1)) + &
                                0.5_wp*cfl * (rho_R*vel_R(1))*(1._wp/dy(k+1)))*dt)
                            ! !$omp atomic update
                            ! !$acc atomic
                            ! rhs_vf(momxb)%sf(j, k+1, l) =  rhs_vf(momxb)%sf(j, k+1, l) + &
                            !     (0.5_wp * rho_R * vel_R(2)*vel_R(1)*(1._wp/dy(k+1)) + &
                            !     0.5_wp*cfl * (rho_R*vel_R(1))*(1._wp/dy(k+1)))

                            @:ATOMIC_ADD(rhs_vf(momxb+2)%sf(j, k+1, l),(0.5_wp * rho_R * vel_R(2)*vel_R(3)*(1._wp/dy(k+1)) + &
                                0.5_wp*cfl * (rho_R*vel_R(3))*(1._wp/dy(k+1)))*dt)
                            ! !$omp atomic update
                            ! !$acc atomic
                            ! rhs_vf(momxb+2)%sf(j, k+1, l) =  rhs_vf(momxb+2)%sf(j, k+1, l) + &
                            !     (0.5_wp * rho_R * vel_R(2)*vel_R(3)*(1._wp/dy(k+1)) + &
                            !     0.5_wp*cfl * (rho_R*vel_R(3))*(1._wp/dy(k+1)))

                            @:ATOMIC_ADD(rhs_vf(E_idx)%sf(j, k+1, l),(0.5_wp * (vel_R(2) * (E_R+ &
                                pres_R+F_L ) )*(1._wp/dy(k+1)) + &
                                0.5_wp*cfl * (E_R)*(1._wp/dy(k+1)))*dt)
                            ! !$omp atomic update
                            ! !$acc atomic
                            !  rhs_vf(E_idx)%sf(j, k+1, l) = rhs_vf(E_idx)%sf(j, k+1, l) +  &
                            !     (0.5_wp * (vel_R(2) * (E_R+ &
                            !     pres_R+F_L ) )*(1._wp/dy(k+1)) + &
                            !     0.5_wp*cfl * (E_R)*(1._wp/dy(k+1)))

                            @:ATOMIC_ADD(rhs_vf(momxb+1)%sf(j,k,l),(-0.5_wp* (rho_R * (vel_R(2))**2.0 + &
                                 pres_R+F_L)*(1._wp/dy(k)) - &
                                 0.5_wp*cfl * (rho_R*vel_R(2))*(1._wp/dy(k)))*dt)
                            ! !$omp atomic update
                            ! !$acc atomic
                            ! rhs_vf(momxb+1)%sf(j,k,l) = rhs_vf(momxb+1)%sf(j,k,l) - &
                            !      (0.5_wp* (rho_R * (vel_R(2))**2.0 + &
                            !      pres_R+F_L)*(1._wp/dy(k)) + &
                            !      0.5_wp*cfl * (rho_R*vel_R(2))*(1._wp/dy(k)))

                            @:ATOMIC_ADD(rhs_vf(momxb)%sf(j, k, l),(-0.5_wp * rho_R * vel_R(2)*vel_R(1)*(1._wp/dy(k)) - &
                                0.5_wp*cfl * (rho_R*vel_R(1))*(1._wp/dy(k)))*dt)
                            ! !$omp atomic update
                            ! !$acc atomic
                            ! rhs_vf(momxb)%sf(j, k, l) =  rhs_vf(momxb)%sf(j, k, l) - &
                            !     (0.5_wp * rho_R * vel_R(2)*vel_R(1)*(1._wp/dy(k)) + &
                            !     0.5_wp*cfl * (rho_R*vel_R(1))*(1._wp/dy(k)))

                            @:ATOMIC_ADD(rhs_vf(momxb+2)%sf(j, k, l),(-0.5_wp * rho_R * vel_R(2)*vel_R(3)*(1._wp/dy(k)) - &
                                0.5_wp*cfl * (rho_R*vel_R(3))*(1._wp/dy(k)))*dt)
                            ! !$omp atomic update
                            ! !$acc atomic
                            ! rhs_vf(momxb+2)%sf(j, k, l) =  rhs_vf(momxb+2)%sf(j, k, l) - &
                            !     (0.5_wp * rho_R * vel_R(2)*vel_R(3)*(1._wp/dy(k)) + &
                            !     0.5_wp*cfl * (rho_R*vel_R(3))*(1._wp/dy(k)))

                            @:ATOMIC_ADD(rhs_vf(E_idx)%sf(j, k, l),(-0.5_wp * (vel_R(2) * (E_R + &
                                pres_R+F_L) )*(1._wp/dy(k)) - &
                                0.5_wp*cfl * (E_R)*(1._wp/dy(k)))*dt)
                            ! !$omp atomic update
                            ! !$acc atomic
                            !  rhs_vf(E_idx)%sf(j, k, l) = rhs_vf(E_idx)%sf(j, k, l) -  &
                            !     (0.5_wp * (vel_R(2) * (E_R + &
                            !     pres_R+F_L) )*(1._wp/dy(k)) + &
                            !     0.5_wp*cfl * (E_R)*(1._wp/dy(k)))

                        end do
                    end do
                end do
            end if
        elseif (idir == 3) then
                !$omp target teams loop collapse(3) private(rho_L,gamma_L,pi_inf_L,mu_L,a_L,vel_L,vel_R,pres_L,alpha_L,alpha_R,alpha_rho_L,cfl,F_L,E_L,mu_R,rho_sf_small,alpha_rho_R,vflux_L_arr,vflux_R_arr,dvel_small)
                !$acc parallel loop collapse(3) gang vector default(present) private(rho_L,gamma_L,pi_inf_L,mu_L,a_L,vel_L,vel_R,pres_L,alpha_L,alpha_R,alpha_rho_L,cfl,F_L,E_L,mu_R,rho_sf_small,alpha_rho_R,vflux_L_arr,vflux_R_arr,dvel_small)
                do l = -buff_size+5, p+buff_size-5
                    do k = 0, n
                        do j = 0, m
                            if(viscous) then 
                                vflux_L_arr = 0._wp
                                vflux_R_arr = 0._wp

                                !DIR$ unroll 6
                                !$acc loop seq 
                                do q = -2, 3
                                    dvel_small = 0._wp
                                    !x-direction contributions
                                    !$acc loop seq 
                                    do i = -1, 1
                                        rho_L = 0._wp
                                        !$acc loop seq 
                                        do r = 1, num_fluids
                                            rho_L = rho_L + q_prim_vf(r)%sf(j+i,k,l+q)
                                        end do
                                        rho_sf_small(i) = rho_L
                                    end do

                                    dvel_small(1) = (1/(2._wp*dx(j))) * ( &
                                        q_prim_vf(momxb)%sf(j+1,k,l+q)/rho_sf_small(1) - &
                                        q_prim_vf(momxb)%sf(j-1,k,l+q)/rho_sf_small(-1))
                                    dvel_small(3) = (1/(2._wp*dx(j))) * ( &
                                        q_prim_vf(momxb+2)%sf(j+1,k,l+q)/rho_sf_small(1) - &
                                        q_prim_vf(momxb+2)%sf(j-1,k,l+q)/rho_sf_small(-1))

                                    if (q > -2) then
                                        vflux_L_arr(1) = vflux_L_arr(1) + coeff_L(q)*(dvel_small(3))
                                        vflux_L_arr(3) = vflux_L_arr(3) + coeff_L(q)*(-2._wp*dvel_small(1))/3._wp
                                    end if
                                    if (q < 3) then
                                        vflux_R_arr(1) = vflux_R_arr(1) + coeff_R(q)*(dvel_small(3))
                                        vflux_R_arr(3) = vflux_R_arr(3) + coeff_R(q)*(-2._wp*dvel_small(1))/3._wp
                                    end if

                                    !y-direction contributions
                                    !$acc loop seq 
                                    do i = -1, 1
                                        rho_L = 0._wp
                                        !$acc loop seq 
                                        do r = 1, num_fluids
                                            rho_L = rho_L + q_prim_vf(r)%sf(j,k+i,l+q)
                                        end do
                                        rho_sf_small(i) = rho_L
                                    end do

                                    dvel_small(2) = (1/(2._wp*dy(k))) * ( &
                                        q_prim_vf(momxb+1)%sf(j,k+1,l+q)/rho_sf_small(1) - &
                                        q_prim_vf(momxb+1)%sf(j,k-1,l+q)/rho_sf_small(-1))
                                    dvel_small(3) = (1/(2._wp*dy(k))) * ( &
                                        q_prim_vf(momxb+2)%sf(j,k+1,l+q)/rho_sf_small(1) - &
                                        q_prim_vf(momxb+2)%sf(j,k-1,l+q)/rho_sf_small(-1))

                                    if (q > -2) then
                                        vflux_L_arr(2) = vflux_L_arr(2) + coeff_L(q)*(dvel_small(3))
                                        vflux_L_arr(3) = vflux_L_arr(3) + coeff_L(q)*(-2._wp*dvel_small(2))/3._wp
                                    end if
                                    if (q < 3) then
                                        vflux_R_arr(2) = vflux_R_arr(2) + coeff_R(q)*(dvel_small(3))
                                        vflux_R_arr(3) = vflux_R_arr(3) + coeff_R(q)*(-2._wp*dvel_small(2))/3._wp
                                    end if

                                    !z-direction contributions
                                    !$acc loop seq 
                                    do i = -1, 1
                                        rho_L = 0._wp
                                        !$acc loop seq 
                                        do r = 1, num_fluids
                                            rho_L = rho_L + q_prim_vf(r)%sf(j,k,l+i+q)
                                        end do
                                        rho_sf_small(i) = rho_L
                                    end do
                                    dvel_small(1) = (1/(2._wp*dz(l))) * ( &
                                        q_prim_vf(momxb)%sf(j,k,l+1+q)/rho_sf_small(1) - &
                                        q_prim_vf(momxb)%sf(j,k,l-1+q)/rho_sf_small(-1))
                                    dvel_small(2) = (1/(2._wp*dz(l))) * ( &
                                        q_prim_vf(momxb+1)%sf(j,k,l+1+q)/rho_sf_small(1) - &
                                        q_prim_vf(momxb+1)%sf(j,k,l-1+q)/rho_sf_small(-1))
                                    dvel_small(3) = (1/(2._wp*dz(l))) * ( &
                                        q_prim_vf(momxb+2)%sf(j,k,l+1+q)/rho_sf_small(1) - &
                                        q_prim_vf(momxb+2)%sf(j,k,l-1+q)/rho_sf_small(-1))
                                    if (q > -2) then
                                        vflux_L_arr(1) = vflux_L_arr(1) + coeff_L(q)*(dvel_small(1))
                                        vflux_L_arr(2) = vflux_L_arr(2) + coeff_L(q)*(dvel_small(2))
                                        vflux_L_arr(3) = vflux_L_arr(3) + coeff_L(q)*(4._wp*dvel_small(3))/3._wp
                                    end if
                                    if (q < 3) then
                                        vflux_R_arr(1) = vflux_R_arr(1) + coeff_R(q)*(dvel_small(1))
                                        vflux_R_arr(2) = vflux_R_arr(2) + coeff_R(q)*(dvel_small(2))
                                        vflux_R_arr(3) = vflux_R_arr(3) + coeff_R(q)*(4._wp*dvel_small(3))/3._wp
                                    end if
                                end do
                            end if

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
                            end do 

                            if(num_fluids > 1) then  
                                !$acc loop seq 
                                do i = 1,num_fluids - 1
                                    alpha_L(i) =  (1._wp/60._wp) * (-3._wp * q_prim_vf(E_idx+i)%sf(j, k, l-1) + &
                                                    27._wp * q_prim_vf(E_idx+i)%sf(j, k, l) + &
                                                    47._wp * q_prim_vf(E_idx+i)%sf(j, k, l+1) -   &
                                                    13._wp * q_prim_vf(E_idx+i)%sf(j, k, l+2) + &
                                                    2._wp * q_prim_vf(E_idx+i)%sf(j, k, l+3))                          
                                    alpha_R(i) = (1._wp/60._wp) * (-3._wp * q_prim_vf(E_idx+i)%sf(j, k, l+2) + &
                                                        27._wp * q_prim_vf(E_idx+i)%sf(j, k, l+1) + &
                                                        47._wp * q_prim_vf(E_idx+i)%sf(j, k, l) -   &
                                                        13._wp * q_prim_vf(E_idx+i)%sf(j, k, l-1) + &
                                                        2._wp * q_prim_vf(E_idx+i)%sf(j, k, l-2))
                                end do 
                                alpha_L(num_fluids) = 1._wp - sum(alpha_L(1:num_fluids-1))
                                alpha_R(num_fluids) = 1._wp - sum(alpha_R(1:num_fluids-1))
                            else 
                                alpha_L(1) = 1._wp
                                alpha_R(1) = 1._wp
                            end if

                            rho_L = sum(alpha_rho_L)
                            gamma_L = sum(alpha_L*gammas)
                            pi_inf_L = sum(alpha_L*pi_infs)

                            rho_R = sum(alpha_rho_R)
                            gamma_R = sum(alpha_R*gammas)
                            pi_inf_R = sum(alpha_R*pi_infs)

                            !$acc loop seq 
                            do i = 1, num_dims
                                vel_L(i) = (1._wp/60._wp) * (-3._wp * q_prim_vf(momxb+i-1)%sf(j,k,l-1) + &
                                                27._wp * q_prim_vf(momxb+i-1)%sf(j,k,l) + &
                                                47._wp * q_prim_vf(momxb+i-1)%sf(j,k,l+1) -   &
                                                13._wp * q_prim_vf(momxb+i-1)%sf(j,k,l+2) + &
                                                2._wp * q_prim_vf(momxb+i-1)%sf(j,k,l+3)) / rho_L
                            end do 

                            !$acc loop seq 
                            do i = 1, num_dims
                                vel_R(i) = (1._wp/60._wp) * (-3._wp * q_prim_vf(momxb+i-1)%sf(j,k,l+2) + &
                                                27._wp * q_prim_vf(momxb+i-1)%sf(j,k,l+1) + &
                                                47._wp * q_prim_vf(momxb+i-1)%sf(j,k,l) -   &
                                                13._wp * q_prim_vf(momxb+i-1)%sf(j,k,l-1) + &
                                                2._wp * q_prim_vf(momxb+i-1)%sf(j,k,l-2)) / rho_R
                            end do

                            if(viscous) then

                                if(num_fluids > 1) then 
                                    mu_L = 0._wp
                                    mu_R = 0._wp
                                    !$acc loop seq 
                                    do i = 1, num_fluids - 1
                                        mu_L = (1._wp/60._wp) * (-3._wp * q_prim_vf(E_idx+i)%sf(j, k, l-1) + &
                                                27._wp * q_prim_vf(E_idx+i)%sf(j, k, l) + &
                                                47._wp * q_prim_vf(E_idx+i)%sf(j, k, l+1) -   &
                                                13._wp * q_prim_vf(E_idx+i)%sf(j, k, l+2) + &
                                                2._wp * q_prim_vf(E_idx+i)%sf(j, k, l+3)) / Res(1, i) + mu_L

                                        mu_R = (1._wp/60._wp) * (-3._wp * q_prim_vf(E_idx+i)%sf(j, k, l+2) + &
                                                27._wp * q_prim_vf(E_idx+i)%sf(j, k, l+1) + &
                                                47._wp * q_prim_vf(E_idx+i)%sf(j, k, l) -   &
                                                13._wp * q_prim_vf(E_idx+i)%sf(j, k, l-1) + &
                                                2._wp * q_prim_vf(E_idx+i)%sf(j, k, l-2)) / Res(1, i) + mu_R
                                    end do

                                    mu_L = 1._wp / Res(1, num_fluids) + mu_L
                                    mu_R = 1._wp / Res(1, num_fluids) + mu_R

                                    !$acc loop seq 
                                    do i = 1, num_fluids - 1
                                        mu_L = -(1._wp/60._wp) * (-3._wp * q_prim_vf(E_idx+i)%sf(j, k, l-1) + &
                                                27._wp * q_prim_vf(E_idx+i)%sf(j, k, l) + &
                                                47._wp * q_prim_vf(E_idx+i)%sf(j, k, l+1) -   &
                                                13._wp * q_prim_vf(E_idx+i)%sf(j, k, l+2) + &
                                                2._wp * q_prim_vf(E_idx+i)%sf(j, k, l+3)) / Res(1, num_fluids) + mu_L

                                        mu_R = -(1._wp/60._wp) * (-3._wp * q_prim_vf(E_idx+i)%sf(j, k, l+2) + &
                                                27._wp * q_prim_vf(E_idx+i)%sf(j, k, l+1) + &
                                                47._wp * q_prim_vf(E_idx+i)%sf(j, k, l) -   &
                                                13._wp * q_prim_vf(E_idx+i)%sf(j, k, l-1) + &
                                                2._wp * q_prim_vf(E_idx+i)%sf(j, k, l-2)) / Res(1, num_fluids) + mu_R
                                    end do
                                else
                                    mu_L = 1._wp / Res(1, 1) 
                                    mu_R = 1._wp / Res(1, 1) 
                                end if

                                @:ATOMIC_ADD(rhs_vf(momxb)%sf(j,k,l+1),(- 0.5_wp*mu_L*vflux_L_arr(1)*(1._wp/dz(l+1))*dt))
                                !!$omp atomic update
                                !!$acc atomic
                                !rhs_vf(momxb)%sf(j,k,l+1) = rhs_vf(momxb)%sf(j,k,l+1) - 0.5_wp*mu_L*vflux_L_arr(1)*(1._wp/dz(l+1))
                                @:ATOMIC_ADD(rhs_vf(E_idx)%sf(j,k,l+1),(- 0.5_wp*mu_L*vflux_L_arr(1)*vel_L(1)*(1._wp/dz(l+1))*dt))
                                !!$omp atomic update
                                !!$acc atomic
                                !rhs_vf(E_idx)%sf(j,k,l+1) = rhs_vf(E_idx)%sf(j,k,l+1) - 0.5_wp*mu_L*vflux_L_arr(1)*vel_L(1)*(1._wp/dz(l+1)) 

                                @:ATOMIC_ADD(rhs_vf(momxb)%sf(j,k,l),(0.5_wp*mu_L*vflux_L_arr(1)*(1._wp/dz(l))*dt))
                                !!$omp atomic update
                                !!$acc atomic
                                !rhs_vf(momxb)%sf(j,k,l) = rhs_vf(momxb)%sf(j,k,l) + 0.5_wp*mu_L*vflux_L_arr(1)*(1._wp/dz(l))
                                @:ATOMIC_ADD(rhs_vf(E_idx)%sf(j,k,l),(0.5_wp*mu_L*vflux_L_arr(1)*vel_L(1)*(1._wp/dz(l))*dt))
                                !!$omp atomic update
                                !!$acc atomic
                                !rhs_vf(E_idx)%sf(j,k,l) = rhs_vf(E_idx)%sf(j,k,l) + 0.5_wp*mu_L*vflux_L_arr(1)*vel_L(1)*(1._wp/dz(l))
                            
                                @:ATOMIC_ADD(rhs_vf(momxb)%sf(j,k,l+1),(- 0.5_wp*mu_R*vflux_R_arr(1)*(1._wp/dz(l+1))*dt))
                                !!$omp atomic update
                                !!$acc atomic
                                !rhs_vf(momxb)%sf(j,k,l+1) = rhs_vf(momxb)%sf(j,k,l+1) - 0.5_wp*mu_R*vflux_R_arr(1)*(1._wp/dz(l+1))
                                @:ATOMIC_ADD(rhs_vf(E_idx)%sf(j,k,l+1),(- 0.5_wp*mu_R*vflux_R_arr(1)*vel_R(1)*(1._wp/dz(l+1))*dt))
                                !!$omp atomic update
                                !!$acc atomic
                                !rhs_vf(E_idx)%sf(j,k,l+1) = rhs_vf(E_idx)%sf(j,k,l+1) - 0.5_wp*mu_R*vflux_R_arr(1)*vel_R(1)*(1._wp/dz(l+1))
                            
                                @:ATOMIC_ADD(rhs_vf(momxb)%sf(j,k,l),(0.5_wp*mu_R*vflux_R_arr(1)*(1._wp/dz(l))*dt))
                                !!$omp atomic update
                                !!$acc atomic
                                !rhs_vf(momxb)%sf(j,k,l) = rhs_vf(momxb)%sf(j,k,l) + 0.5_wp*mu_R*vflux_R_arr(1)*(1._wp/dz(l))
                                @:ATOMIC_ADD(rhs_vf(E_idx)%sf(j,k,l),(0.5_wp*mu_R*vflux_R_arr(1)*vel_R(1)*(1._wp/dz(l))*dt))
                                !!$omp atomic update
                                !!$acc atomic
                                !rhs_vf(E_idx)%sf(j,k,l) = rhs_vf(E_idx)%sf(j,k,l) + 0.5_wp*mu_R*vflux_R_arr(1)*vel_R(1)*(1._wp/dz(l))

                                @:ATOMIC_ADD(rhs_vf(momxb+1)%sf(j,k,l+1),(- 0.5_wp*mu_L*vflux_L_arr(2)*(1._wp/dz(l+1))*dt))
                                !!$omp atomic update
                                !!$acc atomic
                                !rhs_vf(momxb+1)%sf(j,k,l+1) = rhs_vf(momxb+1)%sf(j,k,l+1) - 0.5_wp*mu_L*vflux_L_arr(2)*(1._wp/dz(l+1))
                                @:ATOMIC_ADD(rhs_vf(E_idx)%sf(j,k,l+1),(- 0.5_wp*mu_L*vflux_L_arr(2)*vel_L(2)*(1._wp/dz(l+1))*dt))
                                !!$omp atomic update
                                !!$acc atomic
                                !rhs_vf(E_idx)%sf(j,k,l+1) = rhs_vf(E_idx)%sf(j,k,l+1) - 0.5_wp*mu_L*vflux_L_arr(2)*vel_L(2)*(1._wp/dz(l+1))
                        
                                @:ATOMIC_ADD(rhs_vf(momxb+1)%sf(j,k,l),(0.5_wp*mu_L*vflux_L_arr(2)*(1._wp/dz(l))*dt))
                                !!$omp atomic update
                                !!$acc atomic
                                !rhs_vf(momxb+1)%sf(j,k,l) = rhs_vf(momxb+1)%sf(j,k,l) + 0.5_wp*mu_L*vflux_L_arr(2)*(1._wp/dz(l))
                                @:ATOMIC_ADD(rhs_vf(E_idx)%sf(j,k,l),(0.5_wp*mu_L*vflux_L_arr(2)*vel_L(2)*(1._wp/dz(l))*dt))
                                !!$omp atomic update
                                !!$acc atomic
                                !rhs_vf(E_idx)%sf(j,k,l) = rhs_vf(E_idx)%sf(j,k,l) + 0.5_wp*mu_L*vflux_L_arr(2)*vel_L(2)*(1._wp/dz(l))                    
                             
                                @:ATOMIC_ADD(rhs_vf(momxb+1)%sf(j,k,l+1),(- 0.5_wp*mu_R*vflux_R_arr(2)*(1._wp/dz(l+1))*dt))
                                !!$omp atomic update
                                !!$acc atomic
                                !rhs_vf(momxb+1)%sf(j,k,l+1) = rhs_vf(momxb+1)%sf(j,k,l+1) - 0.5_wp*mu_R*vflux_R_arr(2)*(1._wp/dz(l+1))
                                @:ATOMIC_ADD(rhs_vf(E_idx)%sf(j,k,l+1),(- 0.5_wp*mu_R*vflux_R_arr(2)*vel_R(2)*(1._wp/dz(l+1))*dt))
                                !!$omp atomic update
                                !!$acc atomic
                                !rhs_vf(E_idx)%sf(j,k,l+1) = rhs_vf(E_idx)%sf(j,k,l+1) - 0.5_wp*mu_R*vflux_R_arr(2)*vel_R(2)*(1._wp/dz(l+1))
                             
                                @:ATOMIC_ADD(rhs_vf(momxb+1)%sf(j,k,l),(0.5_wp*mu_R*vflux_R_arr(2)*(1._wp/dz(l))*dt))
                                !!$omp atomic update
                                !!$acc atomic
                                !rhs_vf(momxb+1)%sf(j,k,l) = rhs_vf(momxb+1)%sf(j,k,l) + 0.5_wp*mu_R*vflux_R_arr(2)*(1._wp/dz(l))
                                @:ATOMIC_ADD(rhs_vf(E_idx)%sf(j,k,l),(0.5_wp*mu_R*vflux_R_arr(2)*vel_R(2)*(1._wp/dz(l))*dt))
                                !!$omp atomic update
                                !!$acc atomic
                                !rhs_vf(E_idx)%sf(j,k,l) = rhs_vf(E_idx)%sf(j,k,l) + 0.5_wp*mu_R*vflux_R_arr(2)*vel_R(2)*(1._wp/dz(l))



                                @:ATOMIC_ADD(rhs_vf(momxb+2)%sf(j,k,l+1),(- 0.5_wp*mu_L*vflux_L_arr(3)*(1._wp/dz(l+1))*dt))
                                !!$omp atomic update
                                !!$acc atomic
                                !rhs_vf(momxb+2)%sf(j,k,l+1) = rhs_vf(momxb+2)%sf(j,k,l+1) - 0.5_wp*mu_L*vflux_L_arr(3)*(1._wp/dz(l+1))
                                @:ATOMIC_ADD(rhs_vf(E_idx)%sf(j,k,l+1),(- 0.5_wp*mu_L*vflux_L_arr(3)*vel_L(3)*(1._wp/dz(l+1))*dt))
                                !!$omp atomic update
                                !!$acc atomic
                                !rhs_vf(E_idx)%sf(j,k,l+1) = rhs_vf(E_idx)%sf(j,k,l+1) - 0.5_wp*mu_L*vflux_L_arr(3)*vel_L(3)*(1._wp/dz(l+1))
                            
                                @:ATOMIC_ADD(rhs_vf(momxb+2)%sf(j,k,l),(0.5_wp*mu_L*vflux_L_arr(3)*(1._wp/dz(l))*dt))
                                !!$omp atomic update
                                !!$acc atomic
                                !rhs_vf(momxb+2)%sf(j,k,l) = rhs_vf(momxb+2)%sf(j,k,l) + 0.5_wp*mu_L*vflux_L_arr(3)*(1._wp/dz(l))
                                @:ATOMIC_ADD(rhs_vf(E_idx)%sf(j,k,l),(0.5_wp*mu_L*vflux_L_arr(3)*vel_L(3)*(1._wp/dz(l))*dt))
                                !!$omp atomic update
                                !!$acc atomic
                                !rhs_vf(E_idx)%sf(j,k,l) = rhs_vf(E_idx)%sf(j,k,l) + 0.5_wp*mu_L*vflux_L_arr(3)*vel_L(3)*(1._wp/dz(l))
                            
                                @:ATOMIC_ADD(rhs_vf(momxb+2)%sf(j,k,l+1),(- 0.5_wp*mu_R*vflux_R_arr(3)*(1._wp/dz(l+1))*dt))
                                !!$omp atomic update
                                !!$acc atomic
                                !rhs_vf(momxb+2)%sf(j,k,l+1) = rhs_vf(momxb+2)%sf(j,k,l+1) - 0.5_wp*mu_R*vflux_R_arr(3)*(1._wp/dz(l+1))
                                @:ATOMIC_ADD(rhs_vf(E_idx)%sf(j,k,l+1),(- 0.5_wp*mu_R*vflux_R_arr(3)*vel_R(3)*(1._wp/dz(l+1))*dt))
                                !!$omp atomic update
                                !!$acc atomic
                                !rhs_vf(E_idx)%sf(j,k,l+1) = rhs_vf(E_idx)%sf(j,k,l+1) - 0.5_wp*mu_R*vflux_R_arr(3)*vel_R(3)*(1._wp/dz(l+1))
                            
                                @:ATOMIC_ADD(rhs_vf(momxb+2)%sf(j,k,l),(0.5_wp*mu_R*vflux_R_arr(3)*(1._wp/dz(l))*dt))
                                !!$omp atomic update
                                !!$acc atomic
                                !rhs_vf(momxb+2)%sf(j,k,l) = rhs_vf(momxb+2)%sf(j,k,l) + 0.5_wp*mu_R*vflux_R_arr(3)*(1._wp/dz(l))
                                @:ATOMIC_ADD(rhs_vf(E_idx)%sf(j,k,l),(0.5_wp*mu_R*vflux_R_arr(3)*vel_R(3)*(1._wp/dz(l))*dt))
                                !!$omp atomic update
                                !!$acc atomic
                                !rhs_vf(E_idx)%sf(j,k,l) = rhs_vf(E_idx)%sf(j,k,l) + 0.5_wp*mu_R*vflux_R_arr(3)*vel_R(3)*(1._wp/dz(l))

                                
                            endif

                            F_L = (1._wp/60._wp) * (-3._wp * jac(j, k, l-1) + &
                                                27._wp * jac(j, k, l) + &
                                                47._wp * jac(j, k, l+1) -   &
                                                13._wp * jac(j, k, l+2) + &
                                                2._wp * jac(j, k, l+3))

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

                            pres_L = (E_L - pi_inf_L - 0.5_wp*rho_L*(vel_L(1)**2._wp + vel_L(2)**2._wp + vel_L(3)**2._wp ))/gamma_L                    

                            pres_R = (E_R - pi_inf_R - 0.5_wp*rho_R*(vel_R(1)**2._wp + vel_R(2)**2._wp + vel_R(3)**2._wp ))/gamma_R 

                            a_L = sqrt((pres_L*(1._wp/gamma_L + 1._wp) + pi_inf_L / gamma_L) / rho_L)
                            a_R = sqrt((pres_R*(1._wp/gamma_R + 1._wp) + pi_inf_R / gamma_R) / rho_R)

                            cfl = (max(sqrt(vel_L(1)**2._wp + vel_L(2)**2._wp + vel_L(3)**2._wp),sqrt(vel_R(1)**2._wp + vel_R(2)**2._wp + vel_R(3)**2._wp)) + max(a_L,a_R))

                            !$acc loop seq
                            do i = 1, num_fluids
                                @:ATOMIC_ADD(rhs_vf(i)%sf(j,k,l+1),(0.5_wp * (alpha_rho_L(i) * &
                                    vel_L(3))*(1._wp/dz(l+1)) - &
                                    0.5_wp*cfl * (alpha_rho_L(i))*(1._wp/dz(l+1)))*dt)
                                ! !$omp atomic update
                                ! !$acc atomic
                                ! rhs_vf(i)%sf(j,k,l+1) = rhs_vf(i)%sf(j,k,l+1) + &
                                !     (0.5_wp * (alpha_rho_L(i) * &
                                !     vel_L(3))*(1._wp/dz(l+1)) - &
                                !     0.5_wp*cfl * (alpha_rho_L(i))*(1._wp/dz(l+1)))

                                @:ATOMIC_ADD(rhs_vf(i)%sf(j,k,l),(-0.5_wp * (alpha_rho_L(i) * &
                                    vel_L(3))*(1._wp/dz(l)) + &
                                    0.5_wp*cfl * (alpha_rho_L(i))*(1._wp/dz(l)))*dt)
                                ! !$omp atomic update
                                ! !$acc atomic
                                ! rhs_vf(i)%sf(j,k,l) = rhs_vf(i)%sf(j,k,l) - &
                                !     (0.5_wp * (alpha_rho_L(i) * &
                                !     vel_L(3))*(1._wp/dz(l)) - &
                                !     0.5_wp*cfl * (alpha_rho_L(i))*(1._wp/dz(l))) 
                            end do

                            if(num_fluids > 1) then 
                                !$acc loop seq
                                do i = 1, num_fluids - 1
                                    @:ATOMIC_ADD(rhs_vf(advxb+i-1)%sf(j,k,l+1),(0.5_wp * (alpha_L(i) * &
                                        vel_L(3))*(1._wp/dz(l+1)) - &
                                        0.5_wp*cfl*(alpha_L(i))*(1._wp/dz(l+1)))*dt)
                                    ! !$omp atomic update
                                    ! !$acc atomic
                                    ! rhs_vf(advxb+i-1)%sf(j,k,l+1) = rhs_vf(advxb+i-1)%sf(j,k,l+1) + &
                                    !     (0.5_wp * (alpha_L(i) * &
                                    !     vel_L(3))*(1._wp/dz(l+1)) - &
                                    !     0.5_wp*cfl*(alpha_L(i))*(1._wp/dz(l+1)))

                                    @:ATOMIC_ADD(rhs_vf(advxb+i-1)%sf(j,k,l+1),(-0.5_wp * q_prim_vf(advxb+i-1)%sf(j,k,l+1)  * vel_L(3)*(1._wp/dz(l+1)))*dt)
                                    ! !$omp atomic update
                                    ! !$acc atomic
                                    ! rhs_vf(advxb+i-1)%sf(j,k,l+1) = rhs_vf(advxb+i-1)%sf(j,k,l+1) &
                                    ! - (0.5_wp * q_prim_vf(advxb+i-1)%sf(j,k,l+1)  * vel_L(3)*(1._wp/dz(l+1)))

                                    @:ATOMIC_ADD(rhs_vf(advxb+i-1)%sf(j,k,l),(-0.5_wp * (alpha_L(i) * &
                                        vel_L(3))*(1._wp/dz(l)) + &
                                        0.5_wp*cfl*(alpha_L(i))*(1._wp/dz(l)))*dt)
                                    ! !$omp atomic update
                                    ! !$acc atomic
                                    ! rhs_vf(advxb+i-1)%sf(j,k,l) = rhs_vf(advxb+i-1)%sf(j,k,l) - &
                                    !     (0.5_wp * (alpha_L(i) * &
                                    !     vel_L(3))*(1._wp/dz(l)) - &
                                    !     0.5_wp*cfl*(alpha_L(i))*(1._wp/dz(l)))

                                    @:ATOMIC_ADD(rhs_vf(advxb+i-1)%sf(j,k,l),(0.5_wp * q_prim_vf(advxb+i-1)%sf(j,k,l)  * vel_L(3)*(1._wp/dz(l)))*dt)
                                    ! !$omp atomic update
                                    ! !$acc atomic
                                    ! rhs_vf(advxb+i-1)%sf(j,k,l) = rhs_vf(advxb+i-1)%sf(j,k,l) &
                                    ! + (0.5_wp * q_prim_vf(advxb+i-1)%sf(j,k,l)  * vel_L(3)*(1._wp/dz(l)))
                                end do
                            end if

                            @:ATOMIC_ADD(rhs_vf(momxb+2)%sf(j,k,l+1),(0.5_wp* (rho_L * (vel_L(3))**2.0 + &
                                 pres_L+F_L)*(1._wp/dz(l+1)) - &
                                 0.5_wp*cfl * (rho_L*vel_L(3))*(1._wp/dz(l+1)))*dt)
                            ! !$omp atomic update
                            ! !$acc atomic
                            ! rhs_vf(momxb+2)%sf(j,k,l+1) = rhs_vf(momxb+2)%sf(j,k,l+1) + &
                            !      (0.5_wp* (rho_L * (vel_L(3))**2.0 + &
                            !      pres_L+F_L)*(1._wp/dz(l+1)) - &
                            !      0.5_wp*cfl * (rho_L*vel_L(3))*(1._wp/dz(l+1)))

                            @:ATOMIC_ADD(rhs_vf(momxb)%sf(j, k,l+1),(0.5_wp * rho_L * vel_L(1)*vel_L(3)*(1._wp/dz(l+1)) - &
                                0.5_wp*cfl * (rho_L*vel_L(1))*(1._wp/dz(l+1)))*dt)
                            ! !$omp atomic update
                            ! !$acc atomic
                            ! rhs_vf(momxb)%sf(j, k,l+1) =  rhs_vf(momxb)%sf(j, k,l+1) + &
                            !     (0.5_wp * rho_L * vel_L(1)*vel_L(3)*(1._wp/dz(l+1)) - &
                            !     0.5_wp*cfl * (rho_L*vel_L(1))*(1._wp/dz(l+1)))

                            @:ATOMIC_ADD(rhs_vf(momxb+1)%sf(j,k,l+1),(0.5_wp * rho_L * vel_L(2)*vel_L(3)*(1._wp/dz(l+1)) - &
                                0.5_wp*cfl * (rho_L*vel_L(2))*(1._wp/dz(l+1)))*dt)
                            ! !$omp atomic update
                            ! !$acc atomic
                            ! rhs_vf(momxb+1)%sf(j,k,l+1) =  rhs_vf(momxb+1)%sf(j,k,l+1) + &
                            !     (0.5_wp * rho_L * vel_L(2)*vel_L(3)*(1._wp/dz(l+1)) - &
                            !     0.5_wp*cfl * (rho_L*vel_L(2))*(1._wp/dz(l+1)))

                            @:ATOMIC_ADD(rhs_vf(E_idx)%sf(j,k,l+1),(0.5_wp * (vel_L(3) * (E_L + &
                                pres_L+F_L) )*(1._wp/dz(l+1)) - &
                                0.5_wp*cfl * (E_L)*(1._wp/dz(l+1)))*dt)
                            ! !$omp atomic update
                            ! !$acc atomic
                            ! rhs_vf(E_idx)%sf(j,k,l+1) = rhs_vf(E_idx)%sf(j,k,l+1) +  &
                            !     (0.5_wp * (vel_L(3) * (E_L + &
                            !     pres_L+F_L) )*(1._wp/dz(l+1)) - &
                            !     0.5_wp*cfl * (E_L)*(1._wp/dz(l+1)))

                            @:ATOMIC_ADD(rhs_vf(momxb+2)%sf(j,k,l),(-0.5_wp* (rho_L * (vel_L(3))**2.0 + &
                                 pres_L+F_L)*(1._wp/dz(l)) + &
                                 0.5_wp*cfl * (rho_L*vel_L(3))*(1._wp/dz(l)))*dt)
                            ! !$omp atomic update
                            ! !$acc atomic
                            ! rhs_vf(momxb+2)%sf(j,k,l) = rhs_vf(momxb+2)%sf(j,k,l) - &
                            !      (0.5_wp* (rho_L * (vel_L(3))**2.0 + &
                            !      pres_L+F_L)*(1._wp/dz(l)) - &
                            !      0.5_wp*cfl * (rho_L*vel_L(3))*(1._wp/dz(l)))

                            @:ATOMIC_ADD(rhs_vf(momxb)%sf(j, k, l),(-0.5_wp * rho_L * vel_L(1)*vel_L(3)*(1._wp/dz(l)) + &
                                0.5_wp*cfl * (rho_L*vel_L(1))*(1._wp/dz(l)))*dt)
                            ! !$omp atomic update
                            ! !$acc atomic
                            ! rhs_vf(momxb)%sf(j, k, l) =  rhs_vf(momxb)%sf(j, k, l) - &
                            !     (0.5_wp * rho_L * vel_L(1)*vel_L(3)*(1._wp/dz(l)) - &
                            !     0.5_wp*cfl * (rho_L*vel_L(1))*(1._wp/dz(l)))

                            @:ATOMIC_ADD(rhs_vf(momxb+1)%sf(j, k, l),(-0.5_wp * rho_L * vel_L(2)*vel_L(3)*(1._wp/dz(l)) + &
                                0.5_wp*cfl * (rho_L*vel_L(2))*(1._wp/dz(l)))*dt)
                            ! !$omp atomic update
                            ! !$acc atomic
                            ! rhs_vf(momxb+1)%sf(j, k, l) =  rhs_vf(momxb+1)%sf(j, k, l) - &
                            !     (0.5_wp * rho_L * vel_L(2)*vel_L(3)*(1._wp/dz(l)) - &
                            !     0.5_wp*cfl * (rho_L*vel_L(2))*(1._wp/dz(l)))

                            @:ATOMIC_ADD(rhs_vf(E_idx)%sf(j, k, l),(-0.5_wp * (vel_L(3) * (E_L + &
                                pres_L+F_L) )*(1._wp/dz(l)) + &
                                0.5_wp*cfl * (E_L)*(1._wp/dz(l)))*dt)
                            ! !$omp atomic update
                            ! !$acc atomic
                            ! rhs_vf(E_idx)%sf(j, k, l) = rhs_vf(E_idx)%sf(j, k, l) -  &
                            !     (0.5_wp * (vel_L(3) * (E_L + &
                            !     pres_L+F_L) )*(1._wp/dz(l)) - &
                            !     0.5_wp*cfl * (E_L)*(1._wp/dz(l)))

                            F_L = (1._wp/60._wp) * (-3._wp * jac(j, k, l+2) + &
                                                27._wp * jac(j, k, l+1) + &
                                                47._wp * jac(j, k, l) -   &
                                                13._wp * jac(j, k, l-1) + &
                                                2._wp * jac(j, k, l-2))
                            !$acc loop seq
                            do i = 1, num_fluids
                                @:ATOMIC_ADD(rhs_vf(i)%sf(j,k,l+1),(0.5_wp * (alpha_rho_R(i) * &
                                    vel_R(3))*(1._wp/dz(l+1)) + &
                                    0.5_wp*cfl * (alpha_rho_R(i))*(1._wp/dz(l+1)))*dt)
                                ! !$omp atomic update
                                ! !$acc atomic
                                ! rhs_vf(i)%sf(j,k,l+1) = rhs_vf(i)%sf(j,k,l+1) + &
                                !     (0.5_wp * (alpha_rho_R(i) * &
                                !     vel_R(3))*(1._wp/dz(l+1)) + &
                                !     0.5_wp*cfl * (alpha_rho_R(i))*(1._wp/dz(l+1)))

                                @:ATOMIC_ADD(rhs_vf(i)%sf(j,k,l),(-0.5_wp * (alpha_rho_R(i) * &
                                    vel_R(3))*(1._wp/dz(l)) - &
                                    0.5_wp*cfl * (alpha_rho_R(i))*(1._wp/dz(l)))*dt)
                                ! !$omp atomic update
                                ! !$acc atomic
                                ! rhs_vf(i)%sf(j,k,l) = rhs_vf(i)%sf(j,k,l) - &
                                !     (0.5_wp * (alpha_rho_R(i) * &
                                !     vel_R(3))*(1._wp/dz(l)) + &
                                !     0.5_wp*cfl * (alpha_rho_R(i))*(1._wp/dz(l)))
                            end do

                            if(num_fluids > 1) then 
                                !$acc loop seq
                                do i = 1, num_fluids - 1
                                    @:ATOMIC_ADD(rhs_vf(advxb+i-1)%sf(j,k,l+1),(0.5_wp * (alpha_R(i) * &
                                        vel_R(3))*(1._wp/dz(l+1)) + &
                                        0.5_wp*cfl*(alpha_R(i))*(1._wp/dz(l+1)))*dt)
                                    ! !$omp atomic update
                                    ! !$acc atomic
                                    ! rhs_vf(advxb+i-1)%sf(j,k,l+1) = rhs_vf(advxb+i-1)%sf(j,k,l+1) + &
                                    !     (0.5_wp * (alpha_R(i) * &
                                    !     vel_R(3))*(1._wp/dz(l+1)) + &
                                    !     0.5_wp*cfl*(alpha_R(i))*(1._wp/dz(l+1)))

                                    @:ATOMIC_ADD(rhs_vf(advxb+i-1)%sf(j,k,l+1),(-0.5_wp * q_prim_vf(advxb+i-1)%sf(j,k,l+1)  * vel_R(3)*(1._wp/dz(l+1)))*dt)
                                    ! !$omp atomic update
                                    ! !$acc atomic
                                    ! rhs_vf(advxb+i-1)%sf(j,k,l+1) = rhs_vf(advxb+i-1)%sf(j,k,l+1) &
                                    ! - (0.5_wp * q_prim_vf(advxb+i-1)%sf(j,k,l+1)  * vel_R(3)*(1._wp/dz(l+1)))

                                    @:ATOMIC_ADD(rhs_vf(advxb+i-1)%sf(j,k,l),(-0.5_wp * (alpha_R(i) * &
                                        vel_R(3))*(1._wp/dz(l)) - &
                                        0.5_wp*cfl*(alpha_R(i))*(1._wp/dz(l)))*dt)
                                    ! !$omp atomic update
                                    ! !$acc atomic
                                    ! rhs_vf(advxb+i-1)%sf(j,k,l) = rhs_vf(advxb+i-1)%sf(j,k,l) - &
                                    !     (0.5_wp * (alpha_R(i) * &
                                    !     vel_R(3))*(1._wp/dz(l)) + &
                                    !     0.5_wp*cfl*(alpha_R(i))*(1._wp/dz(l)))

                                    @:ATOMIC_ADD(rhs_vf(advxb+i-1)%sf(j,k,l),(0.5_wp * q_prim_vf(advxb+i-1)%sf(j,k,l)  * vel_R(3)*(1._wp/dz(l)))*dt)
                                    ! !$omp atomic update
                                    ! !$acc atomic
                                    ! rhs_vf(advxb+i-1)%sf(j,k,l) = rhs_vf(advxb+i-1)%sf(j,k,l) &
                                    ! + (0.5_wp * q_prim_vf(advxb+i-1)%sf(j,k,l)  * vel_R(3)*(1._wp/dz(l)))
                                end do
                            end if

                            @:ATOMIC_ADD(rhs_vf(momxb+2)%sf(j,k,l+1),(0.5_wp* (rho_R * (vel_R(3))**2.0 + &
                                 pres_R+F_L)*(1._wp/dz(l+1)) + &
                                 0.5_wp*cfl * (rho_R*vel_R(3))*(1._wp/dz(l+1)))*dt)
                            ! !$omp atomic update
                            ! !$acc atomic
                            ! rhs_vf(momxb+2)%sf(j,k,l+1) = rhs_vf(momxb+2)%sf(j,k,l+1) + &
                            !      (0.5_wp* (rho_R * (vel_R(3))**2.0 + &
                            !      pres_R+F_L)*(1._wp/dz(l+1)) + &
                            !      0.5_wp*cfl * (rho_R*vel_R(3))*(1._wp/dz(l+1)))

                            @:ATOMIC_ADD(rhs_vf(momxb)%sf(j, k,l+1),(0.5_wp * rho_R * vel_R(1)*vel_R(3)*(1._wp/dz(l+1)) + &
                                0.5_wp*cfl * (rho_R*vel_R(1))*(1._wp/dz(l+1)))*dt)
                            ! !$omp atomic update
                            ! !$acc atomic
                            ! rhs_vf(momxb)%sf(j, k,l+1) =  rhs_vf(momxb)%sf(j, k,l+1) + &
                            !     (0.5_wp * rho_R * vel_R(1)*vel_R(3)*(1._wp/dz(l+1)) + &
                            !     0.5_wp*cfl * (rho_R*vel_R(1))*(1._wp/dz(l+1)))

                            @:ATOMIC_ADD(rhs_vf(momxb+1)%sf(j,k,l+1),(0.5_wp * rho_R * vel_R(2)*vel_R(3)*(1._wp/dz(l+1)) + &
                                0.5_wp*cfl * (rho_R*vel_R(2))*(1._wp/dz(l+1)))*dt)
                            ! !$omp atomic update
                            ! !$acc atomic
                            ! rhs_vf(momxb+1)%sf(j,k,l+1) =  rhs_vf(momxb+1)%sf(j,k,l+1) + &
                            !     (0.5_wp * rho_R * vel_R(2)*vel_R(3)*(1._wp/dz(l+1)) + &
                            !     0.5_wp*cfl * (rho_R*vel_R(2))*(1._wp/dz(l+1)))

                            @:ATOMIC_ADD(rhs_vf(E_idx)%sf(j,k,l+1),(0.5_wp * (vel_R(3) * (E_R + &
                                pres_R+F_L) )*(1._wp/dz(l+1)) + &
                                0.5_wp*cfl * (E_R)*(1._wp/dz(l+1)))*dt)
                            ! !$omp atomic update
                            ! !$acc atomic
                            ! rhs_vf(E_idx)%sf(j,k,l+1) = rhs_vf(E_idx)%sf(j,k,l+1) +  &
                            !     (0.5_wp * (vel_R(3) * (E_R + &
                            !     pres_R+F_L) )*(1._wp/dz(l+1)) + &
                            !     0.5_wp*cfl * (E_R)*(1._wp/dz(l+1)))

                            @:ATOMIC_ADD(rhs_vf(momxb+2)%sf(j,k,l),(-0.5_wp* (rho_R * (vel_R(3))**2.0 + &
                                 pres_R+F_L)*(1._wp/dz(l)) - &
                                 0.5_wp*cfl * (rho_R*vel_R(3))*(1._wp/dz(l)))*dt)
                            ! !$omp atomic update
                            ! !$acc atomic
                            ! rhs_vf(momxb+2)%sf(j,k,l) = rhs_vf(momxb+2)%sf(j,k,l) - &
                            !      (0.5_wp* (rho_R * (vel_R(3))**2.0 + &
                            !      pres_R+F_L)*(1._wp/dz(l)) + &
                            !      0.5_wp*cfl * (rho_R*vel_R(3))*(1._wp/dz(l)))

                            @:ATOMIC_ADD(rhs_vf(momxb)%sf(j, k, l),(-0.5_wp * rho_R * vel_R(1)*vel_R(3)*(1._wp/dz(l)) - &
                                0.5_wp*cfl * (rho_R*vel_R(1))*(1._wp/dz(l)))*dt)
                            ! !$omp atomic update
                            ! !$acc atomic
                            ! rhs_vf(momxb)%sf(j, k, l) =  rhs_vf(momxb)%sf(j, k, l) - &
                            !     (0.5_wp * rho_R * vel_R(1)*vel_R(3)*(1._wp/dz(l)) + &
                            !     0.5_wp*cfl * (rho_R*vel_R(1))*(1._wp/dz(l)))

                            @:ATOMIC_ADD(rhs_vf(momxb+1)%sf(j, k, l),(-0.5_wp * rho_R * vel_R(2)*vel_R(3)*(1._wp/dz(l)) - &
                                0.5_wp*cfl * (rho_R*vel_R(2))*(1._wp/dz(l)))*dt)
                            ! !$omp atomic update
                            ! !$acc atomic
                            ! rhs_vf(momxb+1)%sf(j, k, l) =  rhs_vf(momxb+1)%sf(j, k, l) - &
                            !     (0.5_wp * rho_R * vel_R(2)*vel_R(3)*(1._wp/dz(l)) + &
                            !     0.5_wp*cfl * (rho_R*vel_R(2))*(1._wp/dz(l)))

                            @:ATOMIC_ADD(rhs_vf(E_idx)%sf(j, k, l),(-0.5_wp * (vel_R(3) * (E_R + &
                                pres_R+F_L) )*(1._wp/dz(l)) - &
                                0.5_wp*cfl * (E_R)*(1._wp/dz(l)))*dt)
                            ! !$omp atomic update
                            ! !$acc atomic
                            ! rhs_vf(E_idx)%sf(j, k, l) = rhs_vf(E_idx)%sf(j, k, l) -  &
                            !     (0.5_wp * (vel_R(3) * (E_R + &
                            !     pres_R+F_L) )*(1._wp/dz(l)) + &
                            !     0.5_wp*cfl * (E_R)*(1._wp/dz(l)))

                        end do
                    end do
                end do
        end if

    end subroutine s_igr_riemann_solver

    subroutine s_initialize_igr(q_prim_vf, rhs_vf)

        type(scalar_field_half), dimension(sys_size), intent(in) :: q_prim_vf, rhs_vf

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
                    do i = 1, vec_size
                        rhs_vf(i)%sf(j,k,l) = 0._2 
                    end do
                end do
            end do
        end do

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
