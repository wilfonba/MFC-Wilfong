#:include 'macros.fpp'
module m_muscl

    ! Dependencies =============================================================
    use m_derived_types        !< Definitions of the derived types

    use m_global_parameters    !< Definitions of the global parameters

    use m_variables_conversion !< State variables type conversion procedures

#ifdef _OPENACC
    use openacc
#endif

    use m_mpi_proxy
    ! ==========================================================================

    private; public :: s_initialize_muscl_module, &
        s_muscl, &
        s_finalize_muscl_module

    integer :: v_size 
    !$acc declare create(v_size)

    type(int_bounds_info) :: is1, is2, is3
    !$acc declare create(is1, is2, is3)

    !> @name The cell-average variables that will be MUSCL-reconstructed. Formerly, they
    !! are stored in v_vf. However, they are transferred to v_rs_wsL and v_rs_wsR
    !! as to be reshaped (RS) and/or characteristically decomposed. The reshaping
    !! allows the muscl procedure to be independent of the coordinate direction of
    !! the reconstruction. Lastly, notice that the left (L) and right (R) results
    !! of the characteristic decomposition are stored in custom-constructed muscl-
    !! stencils (WS) that are annexed to each position of a given scalar field.
    !> @{
    real(kind(0d0)), allocatable, dimension(:, :, :, :) :: v_rs_ws_x, v_rs_ws_y, v_rs_ws_z
    !> @}
    !$acc declare create (v_rs_ws_x, v_rs_ws_y, v_rs_ws_z)

contains

    subroutine s_initialize_muscl_module()

        ! Initializing in x-direction
        is1%beg = -buff_size; is1%end = m - is1%beg
        if (n == 0) then
            is2%beg = 0
        else
            is2%beg = -buff_size; 
        end if

        is2%end = n - is2%beg

        if (p == 0) then
            is3%beg = 0
        else
            is3%beg = -buff_size
        end if

        is3%end = p - is3%beg

        @:ALLOCATE(v_rs_ws_x(is1%beg:is1%end, &
                    is2%beg:is2%end, is3%beg:is3%end, 1:sys_size))

        if (n == 0) return

        ! initializing in y-direction
        is2%beg = -buff_size; is2%end = n - is2%beg
        is1%beg = -buff_size; is1%end = m - is1%beg

        if (p == 0) then
            is3%beg = 0
        else
            is3%beg = -buff_size
        end if

        is3%end = p - is3%beg

        @:ALLOCATE(v_rs_ws_y(is2%beg:is2%end, &
                    is1%beg:is1%end, is3%beg:is3%end, 1:sys_size))

        if (p == 0) return

        ! initializing in z-direction
        is2%beg = -buff_size; is2%end = n - is2%beg
        is1%beg = -buff_size; is1%end = m - is1%beg
        is3%beg = -buff_size; is3%end = p - is3%beg

        @:ALLOCATE(v_rs_ws_z(is3%beg:is3%end, &
                    is2%beg:is2%end, is1%beg:is1%end, 1:sys_size))

    end subroutine s_initialize_muscl_module

    subroutine s_muscl(v_vf, vL_rs_vf_x, vL_rs_vf_y, vL_rs_vf_z, vR_rs_vf_x, vR_rs_vf_y, vR_rs_vf_z, & 
                          norm_dir, muscl_dir, &
                          is1_d, is2_d, is3_d)

        type(scalar_field), dimension(1:), intent(IN) :: v_vf
        real(kind(0d0)), dimension(startx:, starty:, startz:, 1:), intent(INOUT) :: &
            vL_rs_vf_x, vL_rs_vf_y, &
            vL_rs_vf_z, vR_rs_vf_x, &
            vR_rs_vf_y, vR_rs_vf_z
        integer, intent(IN) :: norm_dir
        integer, intent(IN) :: muscl_dir
        type(int_bounds_info), intent(IN) :: is1_d, is2_d, is3_d

        integer :: j, k, l, i, q
        real(kind(0d0)) :: phir, delta, top, bottom, sign

        is1 = is1_d
        is2 = is2_d
        is3 = is3_d

!$acc update device(is1, is2, is3)

        if (muscl_order /= 1) then
            call s_initialize_muscl(v_vf, &
                                   norm_dir, muscl_dir)
        end if

        if (muscl_order == 1) then
            if (muscl_dir == 1) then
!$acc parallel loop collapse(4) default(present)
                do i = 1, ubound(v_vf, 1)
                    do l = is3%beg, is3%end
                        do k = is2%beg, is2%end
                            do j = is1%beg, is1%end
                                vL_rs_vf_x(j, k, l, i) = v_vf(i)%sf(j, k, l)
                                vR_rs_vf_x(j, k, l, i) = v_vf(i)%sf(j, k, l)
                            end do
                        end do
                    end do
                end do
!$acc end parallel loop
            else if (muscl_dir == 2) then
!$acc parallel loop collapse(4) default(present)
                do i = 1, ubound(v_vf, 1)
                    do l = is3%beg, is3%end
                        do k = is2%beg, is2%end
                            do j = is1%beg, is1%end
                                vL_rs_vf_y(j, k, l, i) = v_vf(i)%sf(k, j, l)
                                vR_rs_vf_y(j, k, l, i) = v_vf(i)%sf(k, j, l)
                            end do
                        end do
                    end do
                end do
!$acc end parallel loop
            else if (muscl_dir == 3) then
!$acc parallel loop collapse(4) default(present)
                do i = 1, ubound(v_vf, 1)
                    do l = is3%beg, is3%end
                        do k = is2%beg, is2%end
                            do j = is1%beg, is1%end
                                vL_rs_vf_z(j, k, l, i) = v_vf(i)%sf(l, k, j)
                                vR_rs_vf_z(j, k, l, i) = v_vf(i)%sf(l, k, j)
                            end do
                        end do
                    end do
                end do
!$acc end parallel loop
            end if
        
        else if (muscl_order == 2) then
            ! MUSCL Reconstruction
            #:for MUSCL_DIR, XYZ in [(1, 'x'), (2, 'y'), (3, 'z')]
            if (muscl_dir == ${MUSCL_DIR}$) then
!$acc parallel loop collapse(4) gang vector default(present) private(top, &
!$acc bottom, sign, phir, delta)
                do l = is3%beg, is3%end
                    do k = is2%beg, is2%end
                            do j = is1%beg, is1%end
                                do i = 1, v_size

                                top = v_rs_ws_${XYZ}$(j, k, l, i) - &
                                    v_rs_ws_${XYZ}$(j - 1, k, l, i)
                                bottom = v_rs_ws_${XYZ}$(j + 1, k, l, i) - &
                                    v_rs_ws_${XYZ}$(j, k, l,  i )

                                call s_compute_slope_limiter(top, bottom, phir)
                                
                                delta = (1d0/2d0)*(v_rs_ws_${XYZ}$(j + 1, k, l, i) - &
                                                v_rs_ws_${XYZ}$(j - 1, k, l, i))

                                ! reconstruct from left side
                                vL_rs_vf_${XYZ}$(j, k, l, i) = &
                                   v_rs_ws_${XYZ}$(j, k, l, i) + 5d-1*phir*delta
                                
                                ! reconstruct from the right side
                                vR_rs_vf_${XYZ}$(j, k, l, i) = &
                                   v_rs_ws_${XYZ}$(j, k, l, i) - 5d-1*phir*delta
                            end do
                        end do
                    end do
                end do
            end if
            #:endfor

        endif

        if (int_comp) then
            call s_interface_compression(vL_rs_vf_x, vL_rs_vf_y, vL_rs_vf_z, &
                vR_rs_vf_x, vR_rs_vf_y, vR_rs_vf_z, &
                norm_dir, muscl_dir, is1_d, is2_d, is3_d)
        endif

    end subroutine s_muscl

    subroutine s_compute_slope_limiter(top, bottom, phir)
!$acc routine seq
    
        real(kind(0d0)) :: top, bottom, phir
        real(kind(0d0)) :: r

        if (r <= 0) then
            phir = 0d0
        else
            if (muscl_lim == 1) then ! ULTRABEE like
                phir = min(2d0*r/(1d0 + r), 2d0/(1d0 + r))
            elseif (muscl_lim == 2) then ! SUPERBEE like 
                if (r <= 5d-1) phir = 2d0*r
                if (r < 1d0) phir = 1d0
                if (r > 1d0) phir = min(r, 2d0/(1d0 + r), 2d0)
            elseif (muscl_lim == 3) then ! van-Alabada like
                phir = min(r*(1d0 + r)/(1d0 + r**2d0), 2d0/(1d0 + r))
            elseif (muscl_lim == 4) then ! MINBEE like
                if (r < 1d0) phir = r
                if (r >=1d0) phir = min(1d0, 2d0/(1d0 + r))
            end if
        end if

    end subroutine s_compute_slope_limiter

    subroutine s_interface_compression(vL_rs_vf_x, vL_rs_vf_y, vL_rs_vf_z, vR_rs_vf_x, vR_rs_vf_y, vR_rs_vf_z, & 
            norm_dir, muscl_dir, &
            is1_d, is2_d, is3_d)

        real(kind(0d0)), dimension(startx:, starty:, startz:, 1:), intent(INOUT) :: &
            vL_rs_vf_x, vL_rs_vf_y, &
            vL_rs_vf_z, vR_rs_vf_x, &
            vR_rs_vf_y, vR_rs_vf_z
        integer, intent(IN) :: norm_dir
        integer, intent(IN) :: muscl_dir
        type(int_bounds_info), intent(IN) :: is1_d, is2_d, is3_d

        integer :: j, k, l, i, q

        real(kind(0d0)) :: iceps, aCL, aCR, aC, aTHINC, qmin, qmax, A, B, C, beta, sign, moncon

        iceps = 1d-4; beta = 1.6d0

        is1 = is1_d
        is2 = is2_d
        is3 = is3_d
!$acc update device(is1, is2, is3)

        #:for MUSCL_DIR, XYZ in [(1, 'x'), (2, 'y'), (3, 'z')]
            if (muscl_dir == ${MUSCL_DIR}$) then
                
!$acc parallel loop collapse(3) gang vector default(present) private(aCL, aC, &
!$acc aCR, aTHINC, moncon, sign, qmin, qmax)
                do l = is3%beg, is3%end
                    do k = is2%beg, is2%end
                        do j = is1%beg, is1%end
                            
                            aCL = v_rs_ws_${XYZ}$(j - 1, k, l, advxb)
                            aC = v_rs_ws_${XYZ}$(j, k, l, advxb)
                            aCR = v_rs_ws_${XYZ}$(j + 1, k, l, advxb)

                            moncon = (aCR - aC)*(aC - aCL)

                            if (aC >= iceps .and. aC <= 1d0 - iceps .and. moncon > 1d-8) then ! Interface cell
                                     
                                if (aCR - aCL > 0d0) then
                                    sign = 1d0
                                else
                                    sign = -1d0
                                endif

                                qmin = min(aCR, aCL)
                                qmax = max(aCR, aCL) - qmin

                                C = (aC - qmin + sgm_eps) / (qmax + sgm_eps)
                                B = exp(sign*beta*(2d0*C - 1d0))
                                A = (B / cosh(beta) - 1d0) / tanh(beta)

                                ! Left reconstruction
                                aTHINC = qmin + 5d-1*qmax*(1d0 + sign * A)
                                if (aTHINC < iceps) aTHINC = iceps
                                if (aTHINC > 1 - iceps) aTHINC = 1 - iceps
                                vL_rs_vf_${XYZ}$(j, k, l, contxb) = vL_rs_vf_${XYZ}$(j, k, l, contxb) / &
                                    vL_rs_vf_${XYZ}$(j, k, l, advxb) * aTHINC
                                vL_rs_vf_${XYZ}$(j, k, l, contxe) = vL_rs_vf_${XYZ}$(j, k, l, contxe) / &
                                    (1d0 - vL_rs_vf_${XYZ}$(j, k, l, advxb)) * (1d0 - aTHINC)
                                vL_rs_vf_${XYZ}$(j, k, l, advxb) = aTHINC
                                vL_rs_vf_${XYZ}$(j, k, l, advxe) = 1 - aTHINC

                                ! Right reconstruction
                                aTHINC = qmin + 5d-1*qmax*(1d0 + sign*(tanh(beta) + A) / (1d0 + A*tanh(beta)))
                                if (aTHINC < iceps) aTHINC = iceps
                                if (aTHINC > 1 - iceps) aTHINC = 1 - iceps
                                vR_rs_vf_${XYZ}$(j, k, l, contxb) = vL_rs_vf_${XYZ}$(j, k, l, contxb) / &
                                    vL_rs_vf_${XYZ}$(j, k, l, advxb) * aTHINC
                                vR_rs_vf_${XYZ}$(j, k, l, contxe) = vL_rs_vf_${XYZ}$(j, k, l, contxe) / &
                                    (1d0 - vL_rs_vf_${XYZ}$(j, k, l, advxb)) * (1d0 - aTHINC)
                                vR_rs_vf_${XYZ}$(j, k, l, advxb) = aTHINC
                                vR_rs_vf_${XYZ}$(j, k, l, advxe) = 1 - aTHINC

                            end if

                        end do
                    end do
                end do
            end if
            #:endfor

    end subroutine s_interface_compression

    subroutine s_initialize_muscl(v_vf, & ! ---------
                                 norm_dir, muscl_dir)

        type(scalar_field), dimension(:), intent(IN) :: v_vf

        integer, intent(IN) :: norm_dir
        integer, intent(IN) :: muscl_dir

        integer :: i, j, k, l, q !< Generic loop iterators

        ! Determining the number of cell-average variables which will be
        ! muscl-reconstructed and mapping their indical bounds in the x-,
        ! y- and z-directions to those in the s1-, s2- and s3-directions
        ! as to reshape the inputted data in the coordinate direction of
        ! the muscl reconstruction
        v_size = ubound(v_vf, 1)    

        !$acc update device(v_size)

        if (muscl_dir == 1) then
!$acc parallel loop collapse(4) gang vector default(present)
            do j = 1, v_size
                do q = is3%beg, is3%end
                    do l = is2%beg, is2%end
                        do k = is1%beg - muscl_polyn, is1%end + muscl_polyn
                            v_rs_ws_x(k, l, q, j) = v_vf(j)%sf(k, l, q)
                        end do
                    end do
                end do
            end do
!$acc end parallel loop
        end if
        ! ==================================================================

        ! Reshaping/Projecting onto Characteristic Fields in y-direction ===
        if (n == 0) return

        if (muscl_dir == 2) then
#if MFC_cuTENSOR
            if (cu_tensor) then
                if (p == 0) then
                    block
                        use CuTensorEx

                        !$acc host_data use_device(v_rs_ws_x, v_rs_ws_y)
     v_rs_ws_y = reshape(v_rs_ws_x, shape=[n + 1 + 2*buff_size, m + 2*buff_size + 1, p + 1, sys_size], order=[2, 1, 3, 4])
                        !$acc end host_data
                    end block
                else
                    block
                        use CuTensorEx

                        !$acc host_data use_device(v_rs_ws_x, v_rs_ws_y)
 v_rs_ws_y = reshape(v_rs_ws_x, shape = [n+1+2*buff_size, m+2*buff_size+1,p+1+2*buff_size,sys_size], order = [2, 1, 3, 4])
                        !$acc end host_data
                    end block
                end if
            else
#endif
!$acc parallel loop collapse(4) gang vector default(present)
                do j = 1, v_size
                    do q = is3%beg, is3%end
                        do l = is2%beg, is2%end
                            do k = is1%beg - muscl_polyn, is1%end + muscl_polyn
                                v_rs_ws_y(k, l, q, j) = v_vf(j)%sf(l, k, q)
                            end do
                        end do
                    end do
                end do
!$acc end parallel loop
#if MFC_cuTENSOR
            end if
#endif
        end if

        ! ==================================================================

        ! Reshaping/Projecting onto Characteristic Fields in z-direction ===
        if (p == 0) return
        if (muscl_dir == 3) then
#if MFC_cuTENSOR
            if (cu_tensor) then
                block
                    use CuTensorEx

                    !$acc host_data use_device(v_rs_ws_x, v_rs_ws_z)
 v_rs_ws_z = reshape(v_rs_ws_x, shape = [p+1+2*buff_size, n+2*buff_size+1,m+2*buff_size+1,sys_size], order = [3, 2, 1, 4])
                    !$acc end host_data
                end block
            else
#endif
!$acc parallel loop collapse(4) gang vector default(present)
                do j = 1, v_size
                    do q = is3%beg, is3%end
                        do l = is2%beg, is2%end
                            do k = is1%beg - muscl_polyn, is1%end + muscl_polyn
                                v_rs_ws_z(k, l, q, j) = v_vf(j)%sf(q, l, k)
                            end do
                        end do
                    end do
                end do
!$acc end parallel loop
#if MFC_cuTENSOR
            end if
#endif
        end if

        ! ==================================================================

    end subroutine s_initialize_muscl ! -------------------------------------

    subroutine s_finalize_muscl_module()

        @:DEALLOCATE(v_rs_ws_x)

        if (n == 0) return

        @:DEALLOCATE(v_rs_ws_y)

        if (p == 0) return

        @:DEALLOCATE(v_rs_ws_z)

    end subroutine s_finalize_muscl_module

end module m_muscl