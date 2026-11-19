!>
!!@file
!!@brief Contains module m_amr_seam

#! AMD OpenMP: assert allocatables present on every kernel (see OMP_DEFAULT_STR); a kernel may name only module arrays allocated
#! under its own launch condition.
#:set MFC_OMP_PRESENT_ALLOCATABLE = True
#:include 'macros.fpp'

!> @brief Same-level fine-fine seams and the coarse conservative halo.
module m_amr_seam

    use m_derived_types
    use m_global_parameters
    use m_mpi_proxy, only: s_mpi_abort  ! @:ASSERT expands to it
    use m_mpi_common, only: s_mpi_sendrecv_variables_buffers
    use m_amr_state
    use m_amr_wave
    use m_amr_distribution

    implicit none

    private
    public :: f_amr_seam, f_amr_seam_dim, s_amr_build_seam_pairs, s_amr_exchange_coarse_cons_halo, s_amr_fine_fine_post, &
        & s_amr_fine_fine_drain, s_amr_fine_fine_level

contains

    !> Exchange coarse conservative ghost layers at internal rank boundaries (beg then end per direction); the solver fills only
    !! prim ghosts. Collective.
    impure subroutine s_amr_exchange_coarse_cons_halo(q_cons)

        type(scalar_field), dimension(sys_size), intent(inout) :: q_cons

        #:for D, X in [(1, 'x'), (2, 'y'), (3, 'z')]
            if (amr_dim(${D}$)) then
                if (bc_${X}$%beg >= 0) call s_mpi_sendrecv_variables_buffers(q_cons, ${D}$, -1, sys_size)
                if (bc_${X}$%end >= 0) call s_mpi_sendrecv_variables_buffers(q_cons, ${D}$, 1, sys_size)
            end if
        #:endfor

    end subroutine s_amr_exchange_coarse_cons_halo

    !> True iff yb sits directly above xb on xb's high face in dim d with identical transverse coarse extents (one ordered pair per
    !! seam).
    pure logical function f_amr_seam(xb, yb, d) result(seam)

        integer, intent(in) :: xb, yb, d
        integer             :: t

        seam = amr_region_lo_all(d, yb) == amr_region_hi_all(d, xb) + 1
        do t = 1, 3
            if (t /= d) seam = seam .and. amr_region_lo_all(t, xb) == amr_region_lo_all(t, yb) .and. amr_region_hi_all(t, &
                & xb) == amr_region_hi_all(t, yb)
        end do

    end function f_amr_seam

    !> Global fine-index interior box [flo, fhi] of block k (0 in collapsed dims).
    pure subroutine s_amr_fine_box(k, flo, fhi)

        integer, intent(in)  :: k
        integer, intent(out) :: flo(3), fhi(3)
        integer              :: fmul

        fmul = amr_ref_ratio**amr_block_level(k)
        flo = merge(fmul*amr_region_lo_all(:,k), 0, amr_dim)
        fhi = merge(fmul*(amr_region_hi_all(:,k) + 1) - 1, 0, amr_dim)

    end subroutine s_amr_fine_box

    !> Part of dst's ghost layer covered by src's interior, in global fine indices (empty if any hi < lo). Same-level blocks are
    !! disjoint, so it is all ghost cells, and it spans faces, edges and corners alike.
    pure subroutine s_amr_halo_box(dst, src, lo, hi)

        integer, intent(in)  :: dst, src
        integer, intent(out) :: lo(3), hi(3)
        integer              :: dlo(3), dhi(3), slo(3), shi(3)

        call s_amr_fine_box(dst, dlo, dhi); call s_amr_fine_box(src, slo, shi)
        lo = max(dlo - merge(buff_size, 0, amr_dim), slo)
        hi = min(dhi + merge(buff_size, 0, amr_dim), shi)

    end subroutine s_amr_halo_box

    !> Same-rank halo copies in one kernel: per transfer t, src slot pt(1) box [pt(3:5), +pt(6:8)) -> dst slot pt(2) shifted by
    !! pt(9:11). Sources are interiors and targets ghosts of disjoint blocks, so transfers never overlap. Threads span the max box.
    impure subroutine s_amr_fine_seam_exchange(nt, pt)

        integer, intent(in) :: nt, pt(:,:)
        integer             :: t, i, g, gmax, n1, n2, j1, j2, j3

        gmax = 0
        do t = 1, nt
            gmax = max(gmax, pt(6, t)*pt(7, t)*pt(8, t))
        end do

        $:GPU_PARALLEL_LOOP(collapse=3, copyin='[pt]', private='[n1, n2, j1, j2, j3]')
        do t = 1, nt
            do i = 1, sys_size
                do g = 0, gmax - 1
                    n1 = pt(6, t); n2 = pt(7, t)
                    if (g < n1*n2*pt(8, t)) then
                        j1 = pt(3, t) + mod(g, n1); j2 = pt(4, t) + mod(g/n1, n2); j3 = pt(5, t) + g/(n1*n2)
                        amr_cons_st(j1 + pt(9, t), j2 + pt(10, t), j3 + pt(11, t), i, pt(2, t)) = amr_cons_st(j1, j2, j3, i, &
                                    & pt(1, t))
                    end if
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_amr_fine_seam_exchange

    !> Seam dim of the ordered pair (xb, yb): the d in which yb is xb's same-level high-face neighbour, else 0 (adjacent in at most
    !! one dim).
    pure integer function f_amr_seam_dim(xb, yb) result(d)

        integer, intent(in) :: xb, yb

        d = 0
        if (xb == yb) return
        if (amr_block_level(xb) /= amr_block_level(yb)) return
        if (f_amr_seam(xb, yb, 1)) d = 1
        if (n_glb > 0) then; if (f_amr_seam(xb, yb, 2)) d = 2; end if
        if (p_glb > 0) then; if (f_amr_seam(xb, yb, 3)) d = 3; end if

    end function f_amr_seam_dim

    !> Rebuild the same-level halo-pair list (xb < yb whose ghost layers reach each other's interior) once per regrid/restart in
    !! identical order on all ranks, and the per-block overlap-rank lists (amr_ovl_gather/scatter).
    impure subroutine s_amr_build_seam_pairs()

        integer :: xb, yb, np, k, mx, pass
        integer :: plo(3), phi(3), rlo(3), rhi(3)

        if (allocated(amr_seam_pairs)) deallocate (amr_seam_pairs)

        ! pass 1 counts, pass 2 fills: keeps amr_seam_pairs exactly sized. Emission order must match on all ranks or sends and
        ! recvs mismatch and deadlock.
        do pass = 1, 2
            np = 0
            do xb = 1, amr_num_blocks - 1
                do yb = xb + 1, amr_num_blocks
                    if (amr_block_level(xb) /= amr_block_level(yb)) cycle
                    call s_amr_halo_box(xb, yb, rlo, rhi)
                    if (any(rhi < rlo)) cycle
                    np = np + 1
                    if (pass == 2) amr_seam_pairs(:,np) = [xb, yb]
                end do
            end do
            if (pass == 1) then
                amr_num_seam_pairs = np
                allocate (amr_seam_pairs(2, max(np, 1)))
            end if
        end do
        ! per-block overlap-rank lists, rank-ascending (same MPI order as a 0..num_procs-1 scan); the first dim is the max
        ! overlap, reallocated each build. Consumers run behind a build_seam_pairs guard.
        mx = 1
        do k = 1, amr_num_blocks
            call s_amr_patch_box(k, plo, phi); call s_amr_region_box(k, rlo, rhi)
            mx = max(mx, f_amr_overlap_count(plo, phi), f_amr_overlap_count(rlo, rhi))
        end do
        if (allocated(amr_ovl_gather)) deallocate (amr_ovl_gather)
        if (allocated(amr_ovl_scatter)) deallocate (amr_ovl_scatter)
        allocate (amr_ovl_gather(mx, amr_max_blocks), amr_ovl_scatter(mx, amr_max_blocks))
        do k = 1, amr_num_blocks
            call s_amr_patch_box(k, plo, phi); call s_amr_region_box(k, rlo, rhi)
            call s_amr_ranks_overlapping(plo, phi, amr_ovl_gather(:,k), amr_ovl_gather_n(k))
            call s_amr_ranks_overlapping(rlo, rhi, amr_ovl_scatter(:,k), amr_ovl_scatter_n(k))
        end do
        amr_seam_pairs_nblk = amr_num_blocks
        amr_seam_pairs_dirty = .false.

    end subroutine s_amr_build_seam_pairs

    !> Fine-fine halo of level lev alone, posted and drained.
    impure subroutine s_amr_fine_fine_level(lev)

        integer, intent(in) :: lev

        call s_amr_fine_fine_post(lev)
        call s_amr_fine_fine_drain()

    end subroutine s_amr_fine_fine_level

    !> Fine-fine halo between same-level blocks (of level lev only, if given): overwrite every ghost cell a neighbour's interior
    !! covers (faces, edges, corners) with that stage-entry interior, cross-rank via the seam wave, same-rank by local copy. Wire is
    !! wp, cast to stp on unpack.
    impure subroutine s_amr_fine_fine_post(lev)

        integer, intent(in), optional :: lev
        integer                       :: idx, s, src, dst, rS, rD, cnt, wlo, whi, lo(3), hi(3), slo(3), dlo(3), shi(3)

        amr_sw_nsame = 0
        if (.not. amr) return
        if (amr_num_blocks < 2) return

        if (amr_seam_pairs_dirty .or. amr_seam_pairs_nblk /= amr_num_blocks) call s_amr_build_seam_pairs()
        if (allocated(amr_sw_pt)) deallocate (amr_sw_pt)
        allocate (amr_sw_pt(11, 2*amr_num_seam_pairs))
        ! Seam wave: one aggregated message per peer and direction; both owners walk the same replicated pair list, so the
        ! layout agrees. The wave stays open across the fill waves (own handle and sides).
        call s_amr_wave_open(amr_wave_seam, amr_band_seam)
        ! transfer records: blk = the block packed/unpacked, (bl, bh) = its local box
        call s_amr_wave_reset(amr_wseam_snd); call s_amr_wave_reset(amr_wseam_rcv)
        do idx = 1, amr_num_seam_pairs
            if (present(lev)) then
                if (amr_block_level(amr_seam_pairs(1, idx)) /= lev) cycle
            end if
            do s = 1, 2
                src = amr_seam_pairs(s, idx); dst = amr_seam_pairs(3 - s, idx)
                rS = amr_block_owner(src); rD = amr_block_owner(dst)
                if (proc_rank /= rS .and. proc_rank /= rD) cycle
                call s_amr_halo_box(dst, src, lo, hi)
                call s_amr_fine_box(src, slo, shi); call s_amr_fine_box(dst, dlo, shi)
                cnt = sys_size*product(hi - lo + 1)
                if (rS == rD) then  ! same rank owns both: deferred to the batched kernel at drain
                    amr_sw_nsame = amr_sw_nsame + 1
                    amr_sw_pt(:,amr_sw_nsame) = [amr_loc_of(src), amr_loc_of(dst), lo - slo, hi - lo + 1, slo - dlo]
                else if (proc_rank == rS) then
                    call s_amr_wave_add(amr_wseam_snd, rD, src, lo - slo, hi - slo, cnt)
                else
                    call s_amr_wave_add(amr_wseam_rcv, rS, dst, lo - dlo, hi - dlo, cnt)
                end if
            end do
        end do
        call s_amr_wave_close(amr_wseam_snd, amr_sw_sq, .false.)
        call s_amr_wave_close(amr_wseam_rcv, amr_sw_rq, .false.)
        call s_amr_wave_post(amr_wave_seam, amr_wseam_rcv, amr_sw_rq, .false.)
        do idx = 1, amr_wseam_snd%nx
            call s_amr_wave_slice(amr_wseam_snd, idx, wlo, whi)
            call s_amr_fine_box_xfer(amr_wseam_snd%blk(idx), amr_wseam_snd%bl(:,idx), amr_wseam_snd%bh(:,idx), &
                                     & amr_sw_sq(wlo:whi), 1)
        end do
        call s_amr_wave_send(amr_wave_seam, amr_wseam_snd, amr_sw_sq, .false.)

    end subroutine s_amr_fine_fine_post

    impure subroutine s_amr_fine_fine_drain()

        integer :: idx, lo, hi

        if (f_amr_wave_nreq(amr_wave_seam) == 0 .and. amr_sw_nsame == 0) return
        call s_amr_wave_wait(amr_wave_seam)
        do idx = 1, amr_wseam_rcv%nx
            call s_amr_wave_slice(amr_wseam_rcv, idx, lo, hi)
            call s_amr_fine_box_xfer(amr_wseam_rcv%blk(idx), amr_wseam_rcv%bl(:,idx), amr_wseam_rcv%bh(:,idx), amr_sw_rq(lo:hi), -1)
        end do
        if (amr_sw_nsame > 0) call s_amr_fine_seam_exchange(amr_sw_nsame, amr_sw_pt(:,1:amr_sw_nsame))
        call s_amr_select_slot(1)

    end subroutine s_amr_fine_fine_drain

end module m_amr_seam
