!>
!!@file
!!@brief Contains module m_amr_nest

#:include 'macros.fpp'

!> @brief Hierarchical child nesting for AMR regrid: each level's children are tagged by the density-gradient sensor on the old
!! parent-level fine solution, clustered, and replayed in rank-invariant order. Also holds the shared tag sensor.
module m_amr_nest

#ifdef MFC_MPI
    use mpi
#endif

    use m_derived_types
    use m_global_parameters
    use m_amr_cluster, only: s_amr_cluster
    use m_eos, only: s_compute_mixture_coefficients, f_pressure
    use m_amr, only: amr_slots, amr_cons_st, amr_loc_of, amr_my_blk, amr_n_my, s_amr_refresh_my_blocks, amr_maxc_fit, &
        & amr_cpat_mar, s_amr_body_bbox, s_amr_expand_box_over_bodies, s_amr_tile_box, f_amr_boxes_overlap

    implicit none

    private
    public :: s_amr_regrid_nest_children, s_amr_nest_window, f_amr_tag_val_sf, f_amr_tag_ind, f_amr_boxes_compatible, &
        & s_amr_nest_seams, f_amr_nest_abuts_ok, f_amr_nest_thr, f_amr_nest_snap

    !> Nested children were dropped for lack of block capacity (warn once per regrid)
    logical :: amr_child_dropped = .false.

contains

    !> Grow the pack buffers sidx/skb geometrically to hold nloc+extra entries, preserving the first nloc.
    impure subroutine s_amr_grow_pack(sidx, skb, nloc, extra)

        integer(8), allocatable, intent(inout) :: sidx(:)
        integer, allocatable, intent(inout)    :: skb(:)
        integer, intent(in)                    :: nloc, extra
        integer                                :: cap, newcap
        integer(8), allocatable                :: t8(:)
        integer, allocatable                   :: ti(:)

        cap = 0
        if (allocated(sidx)) cap = size(sidx)
        if (nloc + extra <= cap) return
        newcap = max(2*cap, max(nloc + extra, 1024))
        allocate (t8(newcap), ti(newcap))
        if (nloc > 0) then
            t8(1:nloc) = sidx(1:nloc)
            ti(1:nloc) = skb(1:nloc)
        end if
        call move_alloc(t8, sidx)
        call move_alloc(ti, skb)

    end subroutine s_amr_grow_pack

    !> Append this rank's tagged cells of child window [mlo:mhi] as (int8 linear index, parent kb) pairs; the encode matches the
    !! s_amr_nest_window_tags decode. gwin is not modified.
    impure subroutine s_amr_pack_gwin_pairs(gwin, mlo, mhi, mg, ng, kb, sidx, skb, nloc)

        integer, intent(in)                    :: mlo(3), mhi(3), mg, ng, kb
        logical, intent(in)                    :: gwin(mlo(1):,mlo(2):,mlo(3):)
        integer(8), allocatable, intent(inout) :: sidx(:)
        integer, allocatable, intent(inout)    :: skb(:)
        integer, intent(inout)                 :: nloc
        integer                                :: gi, gj, gk

        do gk = mlo(3), mhi(3)
            do gj = mlo(2), mhi(2)
                do gi = mlo(1), mhi(1)
                    if (.not. gwin(gi, gj, gk)) cycle
                    call s_amr_grow_pack(sidx, skb, nloc, 1)
                    nloc = nloc + 1
                    sidx(nloc) = int(gi, 8) + int(mg + 1, 8)*(int(gj, 8) + int(ng + 1, 8)*int(gk, 8))
                    skb(nloc) = kb
                end do
            end do
        end do

    end subroutine s_amr_pack_gwin_pairs

    !> Regrid phase 3b: nest child boxes inside each level-(l-1) box for l = 2..amr_max_level (parents first) and set box_level. A
    !! child's extent is the density-gradient sensor on the old parent-level fine solution, coarsened and clustered; a new region
    !! with no fine data gets a centred inset.
    !!
    !! Per level: pass 1 tags each parent's window from owned old blocks and packs (index, kb) pairs; an alltoall routes them to the
    !! round-robin owner of the parent's seam group; pass 2 rebuilds the windows and cuts each group's pooled tags into children;
    !! an allgatherv plus a stable sort by kb replays the children in serial order, so the list and its truncation are
    !! rank-invariant. IB nesting is np=1 only.
    impure subroutine s_amr_regrid_nest_children(boxes, nboxes, box_level)

        type(t_box), allocatable, intent(inout) :: boxes(:)
        integer, intent(inout)                  :: nboxes, box_level(:)
        type(t_box), allocatable                :: grown(:)
        integer                                 :: lev, plo, phi, newlo, kb, ob, jch, nloc_send, ntot_g, nmych, ntot_ch
        integer                                 :: mlo(3), mhi(3)
        integer, allocatable                    :: mlo_all(:,:), mhi_all(:,:), powner(:), mych(:,:), gch(:,:), chord(:)
        integer, allocatable                    :: skb(:), gkb(:), groot(:), cpar(:)
        integer(8), allocatable                 :: sidx(:), gidx(:)
        logical, allocatable                    :: mine(:), covered(:), gwin(:,:,:), sl_all(:,:), sh_all(:,:)

        box_level(1:nboxes) = 1
        if (amr_max_level < 2) return
        amr_child_dropped = .false.

        ! grow to amr_max_blocks for the appends (the IB path leaves boxes at the cluster count)
        if (size(boxes) < amr_max_blocks) then
            allocate (grown(amr_max_blocks)); grown(1:nboxes) = boxes(1:nboxes); call move_alloc(grown, boxes)
        end if
        ! host-refresh the old blocks: the sensor reads the flat store on the host, stale until the stash's update
        do ob = 1, amr_num_blocks
            if (.not. amr_owns_all(ob)) cycle
            $:GPU_UPDATE(host='[amr_cons_st(:, :, :, :, amr_loc_of(ob))]')
        end do
        call s_amr_refresh_my_blocks()

        allocate (cpar(amr_max_blocks))
        plo = 1; phi = nboxes  ! parents
        do lev = 2, amr_max_level
            if (phi < plo) exit
            newlo = nboxes + 1
            allocate (covered(plo:phi), mine(plo:phi), mlo_all(3,plo:phi), mhi_all(3,plo:phi), powner(plo:phi), mych(7, &
                      & amr_max_blocks), groot(plo:phi), sl_all(3,plo:phi), sh_all(3,plo:phi))
            call s_amr_nest_mark_parents(boxes, plo, phi, lev, mine, covered)
            call s_amr_nest_groups(boxes, plo, phi, covered, groot)
            sl_all = .false.; sh_all = .false.
            do kb = plo, phi
                powner(kb) = mod(groot(kb) - plo, max(num_procs, 1))
                if (covered(kb)) call s_amr_nest_seams(boxes(kb), boxes(plo:phi), covered, sl_all(:,kb), sh_all(:,kb))
                call s_amr_nest_window(boxes(kb), lev, mlo_all(:,kb), mhi_all(:,kb), sl_all(:,kb), sh_all(:,kb))
            end do

            ! pass 1: build a window only where this rank holds an overlapping old block (IB: the owner also builds one for the
            ! body)
            nloc_send = 0
            do kb = plo, phi
                mlo = mlo_all(:,kb); mhi = mhi_all(:,kb)
                if (f_amr_nest_window_empty(mlo, mhi)) cycle
                if (.not. (mine(kb) .or. (ib .and. powner(kb) == proc_rank))) cycle
                allocate (gwin(mlo(1):mhi(1),mlo(2):mhi(2),mlo(3):mhi(3)))
                call s_amr_nest_tag_parent(boxes(kb), lev, mlo, mhi, gwin, covered(kb))
                call s_amr_pack_gwin_pairs(gwin, mlo, mhi, m_glb, n_glb, kb, sidx, skb, nloc_send)
                deallocate (gwin)
            end do
            call s_amr_nest_route_pairs(powner, plo, sidx, skb, nloc_send, gidx, gkb, ntot_g)

            ! pass 2: own group roots only (an uncovered parent is its own group)
            nmych = 0
            do kb = plo, phi
                if (groot(kb) /= kb .or. powner(kb) /= proc_rank) cycle
                if (nmych + 1 > amr_max_blocks) exit
                if (covered(kb)) then
                    call s_amr_nest_group(kb, boxes, plo, groot, mlo_all, mhi_all, sl_all, sh_all, lev, gidx, gkb, ntot_g, mych, &
                                          & nmych)
                else if (.not. f_amr_nest_window_empty(mlo_all(:,kb), mhi_all(:,kb))) then
                    call s_amr_nest_inset_child(boxes(kb), lev, kb, mych, nmych)
                end if
            end do

            call s_amr_nest_gather_children(mych, nmych, plo, phi, gch, chord, ntot_ch)
            do jch = 1, ntot_ch
                if (nboxes + 1 > amr_max_blocks) exit  ! canonical order gives the same truncation
                nboxes = nboxes + 1
                boxes(nboxes)%lo = gch(1:3,chord(jch)); boxes(nboxes)%hi = gch(4:6,chord(jch)); box_level(nboxes) = lev
                cpar(nboxes) = gch(7, chord(jch))
            end do
            call s_amr_nest_drop_unmatched(boxes, newlo, nboxes, box_level, cpar, plo, sl_all, sh_all)
            deallocate (gch, chord, powner, mych, gidx, gkb, covered, mine, mlo_all, mhi_all, groot, sl_all, sh_all)
            plo = newlo; phi = nboxes
        end do
        if (nboxes >= amr_max_blocks .and. proc_rank == 0) print '(A)', &
            & ' [amr] NOTE: block pool full during multi-level nesting; some boxes were not refined further'

    end subroutine s_amr_regrid_nest_children

    !> Drop children newlo..nboxes (parents cpar) that end on a seam without an exact partner across (the block cap can truncate the
    !! partner), repeating since a drop can orphan the child across it. Order is kept, so the list stays rank-invariant.
    pure subroutine s_amr_nest_drop_unmatched(boxes, newlo, nboxes, box_level, cpar, plo, sl_all, sh_all)

        type(t_box), intent(inout) :: boxes(:)
        integer, intent(in)        :: newlo, plo, cpar(:)
        integer, intent(inout)     :: nboxes, box_level(:)
        logical, intent(in)        :: sl_all(:,plo:), sh_all(:,plo:)
        logical                    :: keep(newlo:nboxes), gone
        integer                    :: k, n, pk

        keep = .true.; gone = .true.
        do while (gone)
            gone = .false.
            do k = newlo, nboxes
                pk = cpar(k)
                if (.not. keep(k)) cycle
                if (f_amr_nest_abuts_ok(boxes(k), boxes(pk), sl_all(:,pk), sh_all(:,pk), boxes(newlo:nboxes), keep)) cycle
                keep(k) = .false.; gone = .true.
            end do
        end do
        n = newlo - 1
        do k = newlo, nboxes
            if (.not. keep(k)) cycle
            n = n + 1; boxes(n) = boxes(k); box_level(n) = box_level(k)
        end do
        nboxes = n

    end subroutine s_amr_nest_drop_unmatched

    !> mine(kb): this rank holds a level-(lev-1) block overlapping parent kb; covered(kb): any rank does (one LOR; every block has
    !! one owner).
    impure subroutine s_amr_nest_mark_parents(boxes, plo, phi, lev, mine, covered)

        type(t_box), intent(in) :: boxes(:)
        integer, intent(in)     :: plo, phi, lev
        logical, intent(out)    :: mine(plo:), covered(plo:)
        integer                 :: obi, ob, kb

#ifdef MFC_MPI
        integer :: ierr
#endif

        mine = .false.
        do obi = 1, amr_n_my
            ob = amr_my_blk(obi)
            if (amr_block_level(ob) /= lev - 1) cycle
            do kb = plo, phi
                if (f_amr_boxes_overlap(boxes(kb)%lo, boxes(kb)%hi, amr_region_lo_all(:,ob), amr_region_hi_all(:, &
                    & ob))) mine(kb) = .true.
            end do
        end do
        covered = mine
#ifdef MFC_MPI
        if (num_procs > 1) call MPI_ALLREDUCE(MPI_IN_PLACE, covered, phi - plo + 1, MPI_LOGICAL, MPI_LOR, MPI_COMM_WORLD, ierr)
#endif

    end subroutine s_amr_nest_mark_parents

    !> n cells of level lev in L0 cells, rounded up.
    pure integer function f_amr_l0(n, lev) result(c)

        integer, intent(in) :: n, lev

        c = (n + amr_ref_ratio**lev - 1)/amr_ref_ratio**lev

    end function f_amr_l0

    !> Nesting inset of a level-lev child, L0 cells: amr_cpat_mar parent cells.
    pure integer function f_amr_nest_ins(lev) result(ins)

        integer, intent(in) :: lev

        ins = f_amr_l0(amr_cpat_mar, lev - 1)

    end function f_amr_nest_ins

    !> Tag padding of level lev, L0 cells. Level 1 (and IB): amr_buf. Deeper: amr_buf parent cells, at least the level's stencil
    !! (buff_size own cells, +1 on rectangle levels for the seam pull-back), and enough that the next level's padded tags fit its
    !! nesting window unclipped.
    pure recursive integer function f_amr_nest_pad(lev) result(pad)

        integer, intent(in) :: lev

        pad = amr_buf
        if (lev == 1 .or. ib) return
        pad = max(f_amr_l0(amr_buf, lev - 1), f_amr_l0(buff_size, lev) + merge(0, 1, lev == amr_max_level))
        if (lev < amr_max_level) pad = max(pad, f_amr_nest_pad(lev + 1) + f_amr_nest_ins(lev + 1))

    end function f_amr_nest_pad

    !> Raw-gap merge threshold for clustering level lev, L0 cells: padded boxes stay buff_size parent cells (>= 2, the
    !! f_amr_boxes_compatible gap) apart. Level 1: buff_size + 2*amr_buf.
    pure integer function f_amr_nest_thr(lev) result(thr)

        integer, intent(in) :: lev

        thr = merge(buff_size, max(f_amr_l0(buff_size, lev - 1), 2), lev == 1 .or. ib) + 2*f_amr_nest_pad(lev)

    end function f_amr_nest_thr

    !> Regrid hysteresis of level lev, L0 cells: amr_snap, deeper levels capped so a snapped box keeps the level's stencil
    !! (f_amr_l0(buff_size, lev)) around its tags.
    pure integer function f_amr_nest_snap(lev) result(snap)

        integer, intent(in) :: lev

        snap = amr_snap
        if (lev > 1 .and. .not. ib) snap = max(0, min(amr_snap, f_amr_nest_pad(lev) - f_amr_l0(buff_size, lev)))

    end function f_amr_nest_snap

    !> The nesting window of a level-(lev-1) parent: its box inset by f_amr_nest_ins(lev), so a level-lev child's ghost prolongation
    !! reads valid parent interior; no inset on a seam face (sl/sh, s_amr_nest_seams).
    pure subroutine s_amr_nest_window(box, lev, mlo, mhi, sl, sh)

        type(t_box), intent(in) :: box
        integer, intent(in)     :: lev
        integer, intent(out)    :: mlo(3), mhi(3)
        logical, intent(in)     :: sl(3), sh(3)
        integer                 :: ins

        ins = f_amr_nest_ins(lev)
        mlo = box%lo + merge(ins, 0, amr_dim .and. .not. sl); mhi = box%hi - merge(ins, 0, amr_dim .and. .not. sh)

    end subroutine s_amr_nest_window

    !> +d (-d) if box b abuts box a's high (low) face in dim d with identical transverse extents, else 0.
    pure integer function f_amr_nest_face(a, b) result(f)

        type(t_box), intent(in) :: a, b
        integer                 :: d
        logical                 :: tm(3)

        f = 0
        tm = (a%lo == b%lo .and. a%hi == b%hi) .or. .not. amr_dim
        do d = 1, 3
            if (.not. amr_dim(d) .or. .not. all(tm .or. [1, 2, 3] == d)) cycle
            if (b%lo(d) == a%hi(d) + 1) f = d
            if (b%hi(d) == a%lo(d) - 1) f = -d
        end do

    end function f_amr_nest_face

    !> Seams are on unless IB (body nesting keeps every inset) or the patch reach exceeds the parent ghost layer.
    pure logical function f_amr_nest_seams_on() result(on)

        on = .not. ib .and. amr_cpat_mar <= buff_size

    end function f_amr_nest_seams_on

    !> Seams of parent box `box`: its faces shared exactly with a box of cand where use is set (sl/sh: low/high per dim). A child
    !! may end on a seam, since the parent's ghosts there are its neighbour's interior (filled before the child's patch is read) and
    !! the child meets a matching child across (f_amr_nest_abuts_ok).
    pure subroutine s_amr_nest_seams(box, cand, use, sl, sh)

        type(t_box), intent(in) :: box, cand(:)
        logical, intent(in)     :: use(:)
        logical, intent(out)    :: sl(3), sh(3)
        integer                 :: j, f

        sl = .false.; sh = .false.
        if (.not. f_amr_nest_seams_on()) return
        do j = 1, size(cand)
            if (.not. use(j)) cycle
            f = f_amr_nest_face(box, cand(j))
            if (f > 0) sh(f) = .true.
            if (f < 0) sl(-f) = .true.
        end do

    end subroutine s_amr_nest_seams

    !> Seam groups of parents plo..phi: groot(kb) = the lowest covered parent joined to kb through seams (s_amr_nest_seams), so one
    !! rank cuts a group's children and they match across its seams.
    pure subroutine s_amr_nest_groups(boxes, plo, phi, covered, groot)

        type(t_box), intent(in) :: boxes(:)
        integer, intent(in)     :: plo, phi
        logical, intent(in)     :: covered(plo:)
        integer, intent(out)    :: groot(plo:)
        integer                 :: a, b
        logical                 :: grown

        groot = [(a, a=plo, phi)]
        if (.not. f_amr_nest_seams_on()) return
        grown = .true.
        do while (grown)
            grown = .false.
            do a = plo, phi
                do b = a + 1, phi
                    if (.not. (covered(a) .and. covered(b)) .or. groot(a) == groot(b)) cycle
                    if (f_amr_nest_face(boxes(a), boxes(b)) == 0) cycle
                    groot(a) = min(groot(a), groot(b)); groot(b) = groot(a); grown = .true.
                end do
            end do
        end do

    end subroutine s_amr_nest_groups

    !> Child c of parent par (seams sl/sh) puts each face it lays on a seam against a box of kids (use set) exactly: a c/f face
    !! there would face the neighbouring parent, which nothing refluxes.
    pure logical function f_amr_nest_abuts_ok(c, par, sl, sh, kids, use) result(ok)

        type(t_box), intent(in) :: c, par, kids(:)
        logical, intent(in)     :: sl(3), sh(3), use(:)
        logical                 :: need(-3:3)
        integer                 :: j, d, f

        need = .false.
        do d = 1, 3
            if (.not. amr_dim(d)) cycle
            need(d) = sh(d) .and. c%hi(d) == par%hi(d)
            need(-d) = sl(d) .and. c%lo(d) == par%lo(d)
        end do
        do j = 1, size(kids)
            if (.not. use(j)) cycle
            f = f_amr_nest_face(c, kids(j))
            if (f /= 0) need(f) = .false.
        end do
        ok = .not. any(need)

    end function f_amr_nest_abuts_ok

    !> True iff two same-level blocks can coexist: disjoint, and either >= 2 empty coarse cells apart or face-adjacent with an exact
    !! transverse match. Mirrors s_amr_check_seam_topology.
    pure logical function f_amr_boxes_compatible(alo, ahi, blo, bhi) result(ok)

        integer, intent(in) :: alo(3), ahi(3), blo(3), bhi(3)
        integer             :: d, t, gap
        logical             :: tmatch, tover

        ok = .true.
        if (all(alo <= bhi .and. blo <= ahi)) then; ok = .false.; return; end if
        do d = 1, 3
            if (.not. amr_dim(d)) cycle
            gap = max(blo(d) - ahi(d), alo(d) - bhi(d)) - 1  ! empty cells along d (< 0: overlap)
            if (gap >= 2 .or. gap < 0) cycle
            tover = .true.; tmatch = .true.
            do t = 1, 3
                if (t == d .or. .not. amr_dim(t)) cycle
                tover = tover .and. alo(t) <= bhi(t) .and. blo(t) <= ahi(t)
                tmatch = tmatch .and. alo(t) == blo(t) .and. ahi(t) == bhi(t)
            end do
            if (tover .and. .not. (gap == 0 .and. tmatch)) then; ok = .false.; return; end if
        end do

    end function f_amr_boxes_compatible

    !> Too small to nest a child in some active dimension.
    pure logical function f_amr_nest_window_empty(lo, hi) result(e)

        integer, intent(in) :: lo(3), hi(3)

        e = any(amr_dim .and. hi < lo)

    end function f_amr_nest_window_empty

    !> Pass 1 for one parent: tag its window from owned level-(lev-1) blocks (amr_block_level still holds the old levels). IB: also
    !! mark the body bbox with margin max(amr_buf, 4) + amr_cpat_mar, so the clamp to the window keeps a full image-point stencil.
    impure subroutine s_amr_nest_tag_parent(box, lev, mlo, mhi, gwin, covered)

        type(t_box), intent(in) :: box
        integer, intent(in)     :: lev, mlo(3), mhi(3)
        logical, intent(out)    :: gwin(mlo(1):,mlo(2):,mlo(3):)
        logical, intent(inout)  :: covered
        integer                 :: obi, ob, ib_i, bb_lo(3), bb_hi(3)
        logical                 :: any_tag

        gwin = .false.; any_tag = .false.
        do obi = 1, amr_n_my
            ob = amr_my_blk(obi)
            if (amr_block_level(ob) /= lev - 1) cycle
            if (.not. f_amr_boxes_overlap(box%lo, box%hi, amr_region_lo_all(:,ob), amr_region_hi_all(:,ob))) cycle
            call s_amr_tag_child_from_fine(ob, mlo, mhi, gwin, any_tag)
        end do
        if (.not. ib) return
        do ib_i = 1, num_ibs
            call s_amr_body_bbox(ib_i, max(amr_buf, 4) + amr_cpat_mar, bb_lo, bb_hi)
            bb_lo = max(bb_lo, mlo); bb_hi = min(bb_hi, mhi)
            if (f_amr_nest_window_empty(bb_lo, bb_hi)) cycle
            covered = .true.
            gwin(bb_lo(1):bb_hi(1),bb_lo(2):bb_hi(2),bb_lo(3):bb_hi(3)) = .true.
        end do

    end subroutine s_amr_nest_tag_parent

    !> Route each parent's (index, kb) pairs to its clustering owner with a stable counting sort by owner. Serial: the send list is
    !! the gathered list.
    impure subroutine s_amr_nest_route_pairs(powner, plo, sidx, skb, nloc, gidx, gkb, ntot)

        integer, intent(in)                    :: plo, powner(plo:)
        integer(8), allocatable, intent(inout) :: sidx(:)
        integer, allocatable, intent(inout)    :: skb(:)
        integer, intent(in)                    :: nloc
        integer(8), allocatable, intent(out)   :: gidx(:)
        integer, allocatable, intent(out)      :: gkb(:)
        integer, intent(out)                   :: ntot

#ifdef MFC_MPI
        integer                 :: i, ip, ierr
        integer, allocatable    :: rcnt(:), rdsp(:), scnt(:), sdsp(:), phead(:), pord(:), tkb(:)
        integer(8), allocatable :: tidx(:)
#endif

        if (.not. allocated(sidx)) allocate (sidx(0), skb(0))  ! this rank owned no tags at this level
#ifdef MFC_MPI
        if (num_procs > 1) then
            allocate (rcnt(num_procs), rdsp(num_procs), scnt(num_procs), sdsp(num_procs), phead(num_procs), pord(max(nloc, 1)))
            allocate (tidx(max(nloc, 1)), tkb(max(nloc, 1)))
            phead = 0
            do i = 1, nloc
                ip = powner(skb(i)) + 1; phead(ip) = phead(ip) + 1
            end do
            scnt = phead
            sdsp(1) = 0
            do ip = 2, num_procs
                sdsp(ip) = sdsp(ip - 1) + scnt(ip - 1)
            end do
            phead = sdsp + 1
            do i = 1, nloc
                ip = powner(skb(i)) + 1; pord(phead(ip)) = i; phead(ip) = phead(ip) + 1
            end do
            do i = 1, nloc
                tidx(i) = sidx(pord(i)); tkb(i) = skb(pord(i))
            end do
            call MPI_ALLTOALL(scnt, 1, MPI_INTEGER, rcnt, 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)
            rdsp(1) = 0
            do ip = 2, num_procs
                rdsp(ip) = rdsp(ip - 1) + rcnt(ip - 1)
            end do
            ntot = rdsp(num_procs) + rcnt(num_procs)
            allocate (gidx(max(ntot, 1)), gkb(max(ntot, 1)))
            call MPI_ALLTOALLV(tidx, scnt, sdsp, MPI_INTEGER8, gidx, rcnt, rdsp, MPI_INTEGER8, MPI_COMM_WORLD, ierr)
            call MPI_ALLTOALLV(tkb, scnt, sdsp, MPI_INTEGER, gkb, rcnt, rdsp, MPI_INTEGER, MPI_COMM_WORLD, ierr)
            deallocate (sidx, skb)
            return
        end if
#endif
        call move_alloc(sidx, gidx); call move_alloc(skb, gkb)
        ntot = nloc

    end subroutine s_amr_nest_route_pairs

    !> Pass 2 for one parent: rebuild its dense window from the routed pairs with gkb == kb (dedups replicated tags) and extract the
    !! tags in (k,j,i) order. The int8 decode matches s_amr_pack_gwin_pairs.
    impure subroutine s_amr_nest_window_tags(kb, mlo, mhi, gidx, gkb, ntot, ctags, nct)

        integer, intent(in)               :: kb, mlo(3), mhi(3), ntot, gkb(:)
        integer(8), intent(in)            :: gidx(:)
        integer, allocatable, intent(out) :: ctags(:,:)
        integer, intent(out)              :: nct
        logical, allocatable              :: gwin(:,:,:)
        integer                           :: i, gi, gj, gk
        integer(8)                        :: jrem

        allocate (gwin(mlo(1):mhi(1),mlo(2):mhi(2),mlo(3):mhi(3)))
        gwin = .false.
        do i = 1, ntot
            if (gkb(i) /= kb) cycle
            gk = int(gidx(i)/(int(m_glb + 1, 8)*int(n_glb + 1, 8)))
            jrem = gidx(i) - int(gk, 8)*int(m_glb + 1, 8)*int(n_glb + 1, 8)
            gj = int(jrem/int(m_glb + 1, 8))
            gi = int(jrem - int(gj, 8)*int(m_glb + 1, 8))
            gwin(gi, gj, gk) = .true.
        end do
        nct = count(gwin)
        allocate (ctags(3, max(nct, 1)))
        nct = 0
        do gk = mlo(3), mhi(3); do gj = mlo(2), mhi(2); do gi = mlo(1), mhi(1)
            if (gwin(gi, gj, gk)) then
                nct = nct + 1
                ctags(1, nct) = gi; ctags(2, nct) = gj; ctags(3, nct) = gk
            end if
        end do; end do; end do

    end subroutine s_amr_nest_window_tags

    !> Pass 2 for the seam group rooted at root: cut the members' pooled tags jointly (s_amr_nest_cut) so children meet exactly
    !! across each seam. If the cut breaks a seam rule (misaligned parents), redo each member alone in its fully inset window.
    impure subroutine s_amr_nest_group(root, boxes, plo, groot, mlo_all, mhi_all, sl_all, sh_all, lev, gidx, gkb, ntot, mych, nmych)

        integer, intent(in)     :: root, plo, lev, ntot, groot(plo:), mlo_all(:,plo:), mhi_all(:,plo:), gkb(:)
        type(t_box), intent(in) :: boxes(:)
        logical, intent(in)     :: sl_all(:,plo:), sh_all(:,plo:)
        integer(8), intent(in)  :: gidx(:)
        integer, intent(inout)  :: mych(:,:), nmych
        integer, allocatable    :: mem(:), tags(:,:), ctags(:,:)
        integer                 :: kb, nm, nt, nct, n0, im, i, ilo(3), ihi(3)
        logical, parameter      :: no(3, 1) = .false.

        allocate (mem(count(groot == root)), tags(3, max(count([(groot(gkb(i)) == root, i=1, ntot)]), 1)))
        nm = 0; nt = 0
        do kb = root, ubound(groot, 1)
            if (groot(kb) /= root .or. f_amr_nest_window_empty(mlo_all(:,kb), mhi_all(:,kb))) cycle
            nm = nm + 1; mem(nm) = kb
            call s_amr_nest_window_tags(kb, mlo_all(:,kb), mhi_all(:,kb), gidx, gkb, ntot, ctags, nct)
            tags(:,nt + 1:nt + nct) = ctags(:,1:nct); nt = nt + nct
            deallocate (ctags)
        end do
        if (nm == 0) return
        n0 = nmych
        call s_amr_nest_cut(tags, nt, boxes, mem(1:nm), mlo_all(:,mem(1:nm)), mhi_all(:,mem(1:nm)), sl_all(:,mem(1:nm)), sh_all(:, &
                            & mem(1:nm)), lev, mych, nmych)
        if (.not. any(sl_all(:,mem(1:nm)) .or. sh_all(:,mem(1:nm)))) return
        if (f_amr_nest_tiles_ok(mych(:,n0 + 1:nmych), boxes, plo, sl_all, sh_all)) return
        nmych = n0
        do im = 1, nm
            kb = mem(im)
            call s_amr_nest_window(boxes(kb), lev, ilo, ihi, no(:,1), no(:,1))
            if (f_amr_nest_window_empty(ilo, ihi)) cycle
            call s_amr_nest_window_tags(kb, mlo_all(:,kb), mhi_all(:,kb), gidx, gkb, ntot, ctags, nct)
            nt = 0
            do i = 1, nct
                if (any(ctags(:,i) < ilo .or. ctags(:,i) > ihi)) cycle
                nt = nt + 1; tags(:,nt) = ctags(:,i)
            end do
            deallocate (ctags)
            call s_amr_nest_cut(tags, nt, boxes, [kb], reshape(ilo, [3, 1]), reshape(ihi, [3, 1]), no, no, lev, mych, nmych)
        end do

    end subroutine s_amr_nest_group

    !> Cluster tags into children: pad each cluster by f_amr_nest_pad, pull an end that stops exactly on a member's seam back one
    !! cell (the piece on that side would abut nothing), and cut it by each member's window (mem, wlo/whi) into that parent's
    !! children. One member without seams is the plain pad-and-clamp. IB: expand over bodies and re-clamp (the window holds the body
    !! plus margin). The finest level (no IB) is tiled on a lattice instead (s_amr_nest_lattice).
    impure subroutine s_amr_nest_cut(tags, nt, boxes, mem, wlo, whi, sl, sh, lev, mych, nmych)

        integer, intent(in)      :: tags(:,:), nt, mem(:), wlo(:,:), whi(:,:), lev
        type(t_box), intent(in)  :: boxes(:)
        logical, intent(in)      :: sl(:,:), sh(:,:)
        integer, intent(inout)   :: mych(:,:), nmych
        type(t_box), allocatable :: cboxes(:)
        integer                  :: ncb, kc, im, kb, d, pad, clo(3), chi(3), ilo(3), ihi(3)

        if (nt == 0) return
        if (lev == amr_max_level .and. .not. ib) then
            call s_amr_nest_lattice(tags, nt, boxes, mem, wlo, whi, sl, sh, lev, mych, nmych); return
        end if
        pad = f_amr_nest_pad(lev)
        call s_amr_cluster(tags, nt, cboxes, ncb, .false., f_amr_nest_thr(lev))
        do kc = 1, ncb
            clo = merge(cboxes(kc)%lo - pad, 0, amr_dim); chi = merge(cboxes(kc)%hi + pad, 0, amr_dim)
            do im = 1, size(mem)
                kb = mem(im)
                if (.not. f_amr_boxes_overlap(clo, chi, boxes(kb)%lo, boxes(kb)%hi)) cycle
                do d = 1, 3
                    if (sh(d, im) .and. chi(d) == boxes(kb)%hi(d)) chi(d) = chi(d) - 1
                    if (sl(d, im) .and. clo(d) == boxes(kb)%lo(d)) clo(d) = clo(d) + 1
                end do
            end do
            do im = 1, size(mem)
                ilo = clo; ihi = chi
                call s_amr_nest_clamp(ilo, ihi, wlo(:,im), whi(:,im), 0)
                if (f_amr_nest_window_empty(ilo, ihi)) cycle
                if (ib) then
                    call s_amr_expand_box_over_bodies(ilo, ihi)
                    call s_amr_nest_clamp(ilo, ihi, wlo(:,im), whi(:,im), 0)
                end if
                call s_amr_nest_emit(ilo, ihi, lev, mem(im), mych, nmych)
            end do
        end do
        if (allocated(cboxes)) deallocate (cboxes)

    end subroutine s_amr_nest_cut

    !> Finest level: tiles of a lattice that fills the member windows evenly (s_amr_nest_grid) and meet the tags padded by
    !! f_amr_nest_pad, plus the tile across any seam a tile ends on (its piece needs a partner), each cut by the member windows; a
    !! piece whose partner falls in an empty window pulls back a cell. Tiles abut exactly or sit apart by a whole tile or an inset
    !! gap, so the level follows the tags with no merge. Coarser levels keep rectangles: a child must not end on a parent tile face
    !! without a partner across, which tiled parents would force.
    impure subroutine s_amr_nest_lattice(tags, nt, boxes, mem, wlo, whi, sl, sh, lev, mych, nmych)

        integer, intent(in)     :: tags(:,:), nt, mem(:), wlo(:,:), whi(:,:), lev
        type(t_box), intent(in) :: boxes(:)
        logical, intent(in)     :: sl(:,:), sh(:,:)
        integer, intent(inout)  :: mych(:,:), nmych
        logical, allocatable    :: sel(:,:,:)
        integer, allocatable    :: cs(:,:), tix(:,:)
        integer                 :: glo(3), ghi(3), nc(3), lo(3), hi(3), ilo(3), ihi(3), ti(3), pad, i, j, k, im, kb, d, sd
        logical                 :: grown

        glo = huge(0); ghi = -huge(0)
        do im = 1, size(mem)
            glo = min(glo, boxes(mem(im))%lo); ghi = max(ghi, boxes(mem(im))%hi)
        end do
        allocate (cs(maxval(ghi - glo) + 2, 3), tix(0:maxval(ghi - glo),3))
        do d = 1, 3
            call s_amr_nest_grid(d, lev, glo(d), ghi(d), wlo, whi, cs(:,d), tix(:,d), nc(d))
        end do
        allocate (sel(nc(1), nc(2), nc(3)))
        sel = .false.
        pad = f_amr_nest_pad(lev)
        do i = 1, nt
            do d = 1, 3
                lo(d) = tix(max(tags(d, i) - pad, glo(d)) - glo(d), d); hi(d) = tix(min(tags(d, i) + pad, ghi(d)) - glo(d), d)
            end do
            sel(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)) = .true.
        end do
        grown = .true.
        do while (grown)
            grown = .false.
            do k = 1, nc(3); do j = 1, nc(2); do i = 1, nc(1)
                if (.not. sel(i, j, k)) cycle
                ti = [i, j, k]
                lo = [(cs(ti(d), d), d=1, 3)]; hi = [(cs(ti(d) + 1, d) - 1, d=1, 3)]
                do im = 1, size(mem)
                    kb = mem(im)
                    if (.not. f_amr_boxes_overlap(lo, hi, boxes(kb)%lo, boxes(kb)%hi)) cycle
                    do d = 1, 3
                        do sd = -1, 1, 2
                            ti = [i, j, k]
                            if (sd > 0 .and. .not. (sh(d, im) .and. hi(d) == boxes(kb)%hi(d))) cycle
                            if (sd < 0 .and. .not. (sl(d, im) .and. lo(d) == boxes(kb)%lo(d))) cycle
                            ti(d) = ti(d) + sd
                            if (any(ti < 1 .or. ti > nc)) cycle
                            if (sel(ti(1), ti(2), ti(3))) cycle
                            sel(ti(1), ti(2), ti(3)) = .true.; grown = .true.
                        end do
                    end do
                end do
            end do; end do; end do
        end do
        do k = 1, nc(3); do j = 1, nc(2); do i = 1, nc(1)
            if (.not. sel(i, j, k)) cycle
            do im = 1, size(mem)
                kb = mem(im)
                ilo = [cs(i, 1), cs(j, 2), cs(k, 3)]; ihi = [cs(i + 1, 1), cs(j + 1, 2), cs(k + 1, 3)] - 1
                call s_amr_nest_clamp(ilo, ihi, wlo(:,im), whi(:,im), 0)
                ! a piece ending on a seam whose tile across has no piece (that parent's window is empty there) pulls back a cell
                do d = 1, 3
                    ti = [i, j, k]; ti(d) = ti(d) + 1
                    if (sh(d, im) .and. ihi(d) == boxes(kb)%hi(d)) then
                        if (.not. f_tile_cut(ti)) ihi(d) = ihi(d) - 1
                    end if
                    ti(d) = ti(d) - 2
                    if (sl(d, im) .and. ilo(d) == boxes(kb)%lo(d)) then
                        if (.not. f_tile_cut(ti)) ilo(d) = ilo(d) + 1
                    end if
                end do
                if (.not. f_amr_nest_window_empty(ilo, ihi)) call s_amr_nest_emit(ilo, ihi, lev, kb, mych, nmych)
            end do
        end do; end do; end do

    contains

        !> Tile t is selected and some member window holds part of it.
        logical function f_tile_cut(t) result(cut)

            integer, intent(in) :: t(3)
            integer             :: m, plo(3), phi(3)

            cut = .false.
            if (any(t < 1 .or. t > nc)) return
            if (.not. sel(t(1), t(2), t(3))) return
            do m = 1, size(mem)
                plo = [cs(t(1), 1), cs(t(2), 2), cs(t(3), 3)]; phi = [cs(t(1) + 1, 1), cs(t(2) + 1, 2), cs(t(3) + 1, 3)] - 1
                call s_amr_nest_clamp(plo, phi, wlo(:,m), whi(:,m), 0)
                cut = cut .or. .not. f_amr_nest_window_empty(plo, phi)
            end do

        end function f_tile_cut

    end subroutine s_amr_nest_lattice

    !> Lattice lines of s_amr_nest_lattice along dim d over [glo:ghi]: the member window edges, each gap between them split evenly
    !! into tiles of at most amr_blocking_factor cells (>= 2, <= the level cap). Tile t is [cs(t):cs(t+1)-1], t = 1..nc, and
    !! tix(x-glo) is the tile holding x. A collapsed dim is one tile.
    pure subroutine s_amr_nest_grid(d, lev, glo, ghi, wlo, whi, cs, tix, nc)

        integer, intent(in)  :: d, lev, glo, ghi, wlo(:,:), whi(:,:)
        integer, intent(out) :: cs(:), tix(0:), nc
        logical              :: cut(glo:ghi + 1)
        integer              :: tw, im, q, x, a, k, n

        cut = .false.; cut(glo) = .true.; cut(ghi + 1) = .true.
        if (amr_dim(d)) then
            do im = 1, size(wlo, 2)
                cut(min(max(wlo(d, im), glo), ghi + 1)) = .true.; cut(min(max(whi(d, im) + 1, glo), ghi + 1)) = .true.
            end do
            tw = min(max(2, amr_blocking_factor), max(amr_maxc_fit(d)/amr_ref_ratio**(lev - 1), 1))
            a = glo
            do x = glo + 1, ghi + 1
                if (.not. cut(x)) cycle
                n = (x - a + tw - 1)/tw
                do k = 1, n - 1
                    cut(a + (k*(x - a))/n) = .true.
                end do
                a = x
            end do
        end if
        nc = 0
        do x = glo, ghi + 1
            if (.not. cut(x)) cycle
            nc = nc + 1; cs(nc) = x
        end do
        nc = nc - 1
        q = 1
        do x = glo, ghi
            if (x == cs(q + 1)) q = q + 1
            tix(x - glo) = q
        end do

    end subroutine s_amr_nest_grid

    !> The cut children (rows lo, hi, parent kb) keep the seam rules: pairwise compatible, and every face on a seam abuts exactly.
    pure logical function f_amr_nest_tiles_ok(ch, boxes, plo, sl_all, sh_all) result(ok)

        integer, intent(in)     :: ch(:,:), plo
        type(t_box), intent(in) :: boxes(:)
        logical, intent(in)     :: sl_all(:,plo:), sh_all(:,plo:)
        type(t_box)             :: t(size(ch, 2))
        integer                 :: i, j

        do i = 1, size(t)
            t(i)%lo = ch(1:3,i); t(i)%hi = ch(4:6,i)
        end do
        ok = .true.
        do i = 1, size(t)
            do j = i + 1, size(t)
                ok = ok .and. f_amr_boxes_compatible(t(i)%lo, t(i)%hi, t(j)%lo, t(j)%hi)
            end do
            ok = ok .and. f_amr_nest_abuts_ok(t(i), boxes(ch(7, i)), sl_all(:,ch(7, i)), sh_all(:,ch(7, i)), t, [(.true., j=1, &
                                              & size(t))])
        end do

    end function f_amr_nest_tiles_ok

    !> Brand-new region (no old fine data yet): a centred inset so the child still appears this regrid.
    impure subroutine s_amr_nest_inset_child(box, lev, kb, mych, nmych)

        type(t_box), intent(in) :: box
        integer, intent(in)     :: lev, kb
        integer, intent(inout)  :: mych(:,:), nmych
        integer                 :: ins(3), clo(3), chi(3)

        ins = merge(max((box%hi - box%lo + 1)/4, amr_cpat_mar), 0, amr_dim)
        clo = box%lo + ins; chi = box%hi - ins
        if (f_amr_nest_window_empty(clo, chi)) return
        call s_amr_nest_emit(clo, chi, lev, kb, mych, nmych)

    end subroutine s_amr_nest_inset_child

    !> Pad a child by `pad` and clamp it into the nesting window; collapsed dims are pinned to 0.
    pure subroutine s_amr_nest_clamp(clo, chi, mlo, mhi, pad)

        integer, intent(inout) :: clo(3), chi(3)
        integer, intent(in)    :: mlo(3), mhi(3), pad

        clo = merge(max(clo - pad, mlo), 0, amr_dim); chi = merge(min(chi + pad, mhi), 0, amr_dim)

    end subroutine s_amr_nest_clamp

    !> Tile a child to the level's slot cap, amr_maxc_fit/amr_ref_ratio**(lev-1) coarse cells, and append the tiles to this rank's
    !! list.
    impure subroutine s_amr_nest_emit(clo, chi, lev, kb, mych, nmych)

        integer, intent(in)    :: clo(3), chi(3), lev, kb
        integer, intent(inout) :: mych(:,:), nmych
        type(t_box)            :: tiles(amr_max_blocks)
        integer                :: nt, capped, it

        nt = 0; capped = 0
        call s_amr_tile_box(clo, chi, tiles, nt, amr_max_blocks, capped, amr_maxc_fit/amr_ref_ratio**(lev - 1))
        do it = 1, nt
            if (nmych + 1 > amr_max_blocks) then
                ! refinement is dropped: warn once
                if (.not. amr_child_dropped) then
                    write (0, '(A,I0,A)') ' [amr] WARNING: nested refinement dropped, box capacity amr_max_blocks = ', &
                           & amr_max_blocks, ' reached; raise amr_max_blocks'
                    amr_child_dropped = .true.
                end if
                exit
            end if
            nmych = nmych + 1
            mych(1:3,nmych) = tiles(it)%lo; mych(4:6,nmych) = tiles(it)%hi; mych(7, nmych) = kb
        end do

    end subroutine s_amr_nest_emit

    !> Assemble every rank's children (7 ints each: lo, hi, parent kb) and order them by parent kb with a stable counting sort:
    !! chord(1:ntot) indexes gch in the canonical (kb ascending, emission) order.
    impure subroutine s_amr_nest_gather_children(mych, nmych, plo, phi, gch, chord, ntot)

        integer, intent(in)               :: mych(:,:), nmych, plo, phi
        integer, allocatable, intent(out) :: gch(:,:), chord(:)
        integer, intent(out)              :: ntot
        integer, allocatable              :: chhead(:)
        integer                           :: ich, jch, kb

#ifdef MFC_MPI
        integer              :: ip, ierr
        integer, allocatable :: rcnt(:), rdsp(:)

        if (num_procs > 1) then
            allocate (rcnt(num_procs), rdsp(num_procs))
            call MPI_ALLGATHER(nmych, 1, MPI_INTEGER, rcnt, 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)
            rdsp(1) = 0
            do ip = 2, num_procs
                rdsp(ip) = rdsp(ip - 1) + rcnt(ip - 1)
            end do
            ntot = rdsp(num_procs) + rcnt(num_procs)
            allocate (gch(7, max(ntot, 1)))
            rcnt = rcnt*7; rdsp = rdsp*7
            call MPI_ALLGATHERV(mych, nmych*7, MPI_INTEGER, gch, rcnt, rdsp, MPI_INTEGER, MPI_COMM_WORLD, ierr)
        else
#endif
            ntot = nmych
            allocate (gch(7, max(ntot, 1))); gch(:,1:ntot) = mych(:,1:ntot)
#ifdef MFC_MPI
        end if
#endif
        allocate (chhead(plo:phi), chord(max(ntot, 1)))
        chhead = 0
        do ich = 1, ntot
            chhead(gch(7, ich)) = chhead(gch(7, ich)) + 1
        end do
        jch = 1
        do kb = plo, phi
            ich = chhead(kb); chhead(kb) = jch; jch = jch + ich
        end do
        do ich = 1, ntot
            kb = gch(7, ich)
            chord(chhead(kb)) = ich; chhead(kb) = chhead(kb) + 1
        end do

    end subroutine s_amr_nest_gather_children

    !> Sensor-on-fine child tagging: OR the density-gradient tags of an old fine block into an L0 tag grid over a parent window.
    !! Reads the flat store on the host (caller refreshes it); one-sided differences at the fine edges avoid stale ghosts.
    impure subroutine s_amr_tag_child_from_fine(ob, win_lo, win_hi, ctag, any_tag)

        integer, intent(in)    :: ob, win_lo(3), win_hi(3)
        logical, intent(inout) :: ctag(win_lo(1):,win_lo(2):,win_lo(3):)
        logical, intent(inout) :: any_tag
        integer                :: rr, ci, cj, ck, fi, fj, fk, d1, d2, d3, fm1, fm2, fm3, olo(3), lo(3), hi(3)
        real(wp)               :: g
        logical                :: tagged

        rr = amr_ref_ratio**amr_block_level(ob)  ! block-fine cells per L0 cell
        olo = amr_region_lo_all(:,ob)
        fm1 = amr_slots(ob)%m; fm2 = amr_slots(ob)%n; fm3 = amr_slots(ob)%p
        lo = merge(max(win_lo, amr_region_lo_all(:,ob)), 0, amr_dim); hi = merge(min(win_hi, amr_region_hi_all(:,ob)), 0, amr_dim)
        do ck = lo(3), hi(3)
            do cj = lo(2), hi(2)
                do ci = lo(1), hi(1)
                    tagged = .false.
                    do d3 = 0, merge(rr - 1, 0, p_glb > 0)
                        fk = (ck - olo(3))*rr + d3
                        do d2 = 0, merge(rr - 1, 0, n_glb > 0)
                            fj = (cj - olo(2))*rr + d2
                            do d1 = 0, rr - 1
                                fi = (ci - olo(1))*rr + d1
                                g = f_amr_tag_grad_st(amr_loc_of(ob), fi, fj, fk, fm1, fm2, fm3)
                                if (g > amr_tag_eps) tagged = .true.
                            end do
                        end do
                    end do
                    if (tagged) then
                        ctag(ci, cj, ck) = .true.
                        any_tag = .true.
                    end if
                end do
            end do
        end do

    end subroutine s_amr_tag_child_from_fine

    !> The tagged variable at one cell (amr_tag_var): 1 total density (sum of the continuity variables), 2 volume fraction of fluid
    !! 1, 3 velocity magnitude, 4 pressure (five-equation mixture). `_st` reads the flat store, `_sf` the level-0 field.
    #:for RSFX, RSRC in [('st', 'amr_cons_st'), ('sf', '')]
        #:set RARG = 'loc' if RSRC else 'q'
        #:set RHO = 'f_amr_tag_val_' + RSFX
        #:set QPRE = 'real(amr_cons_st(ci, cj, ck, ' if RSRC else 'real(q('
        #:set QPOST = ', loc), wp)' if RSRC else ')%sf(ci, cj, ck), wp)'
        function f_amr_tag_val_${RSFX}$(${RARG}$, ci, cj, ck) result(r)

            #:if RSRC
                integer, intent(in) :: loc  !< flat-store slot
            #:else
                type(scalar_field), dimension(:), intent(in) :: q
            #:endif
            integer, intent(in) :: ci, cj, ck
            real(wp)            :: r, u2, gam, pinf, qv
            integer             :: f
            #:if not MFC_CASE_OPTIMIZATION and USING_AMD
                real(wp), dimension(3) :: ar, al
            #:else
                real(wp), dimension(num_fluids) :: ar, al
            #:endif

            r = 0._wp
            if (amr_tag_var == 2) then
                r = ${QPRE}$eqn_idx%adv%beg${QPOST}$
                return
            end if
            do f = eqn_idx%cont%beg, eqn_idx%cont%end
                r = r + ${QPRE}$f${QPOST}$
            end do
            if (amr_tag_var < 3) return
            u2 = 0._wp
            do f = eqn_idx%mom%beg, eqn_idx%mom%end
                u2 = u2 + ${QPRE}$f${QPOST}$**2
            end do
            if (amr_tag_var == 3) then
                r = sqrt(u2)/max(r, 1.e-30_wp)
                return
            end if
            ar = 0._wp; al = 0._wp
            do f = 1, num_fluids
                ar(f) = ${QPRE}$eqn_idx%cont%beg + f - 1${QPOST}$
                al(f) = ${QPRE}$eqn_idx%adv%beg + f - 1${QPOST}$
            end do
            call s_compute_mixture_coefficients(ar, al, r, gam, pinf, qv)
            r = f_pressure(${QPRE}$eqn_idx%E${QPOST}$ - 0.5_wp*u2/max(r, 1.e-30_wp), gam, pinf, qv)

        end function f_amr_tag_val_${RSFX}$
    #:endfor

    !> Tag indicator (jump of the amr_tag_var variable) at one cell of a block in the flat store, max over the active dimensions.
    !! The stencil is clamped to the block (fm1:fm3 = block extents) since its ghosts are stale.
    function f_amr_tag_grad_st(loc, ci, cj, ck, fm1, fm2, fm3) result(g)

        integer, intent(in) :: loc, fm1, fm2, fm3
        integer, intent(in) :: ci, cj, ck
        real(wp)            :: g, r0
        integer             :: ip(3), im(3)

        ip = min([ci, cj, ck] + 1, [fm1, fm2, fm3]); im = max([ci, cj, ck] - 1, 0)
        r0 = f_amr_tag_val_st(loc, ci, cj, ck)
        g = f_amr_tag_ind(f_amr_tag_val_st(loc, ip(1), cj, ck), r0, f_amr_tag_val_st(loc, im(1), cj, ck))
        if (n_glb > 0) g = max(g, f_amr_tag_ind(f_amr_tag_val_st(loc, ci, ip(2), ck), r0, f_amr_tag_val_st(loc, ci, im(2), ck)))
        if (p_glb > 0) g = max(g, f_amr_tag_ind(f_amr_tag_val_st(loc, ci, cj, ip(3)), r0, f_amr_tag_val_st(loc, ci, cj, im(3))))

    end function f_amr_tag_grad_st

    !> Tag indicator from the tagged variable at i+1, i, i-1: central difference over 2*r_i, or with amr_tag_abs the larger absolute
    !! face jump.
    pure real(wp) function f_amr_tag_ind(rp, r0, rm) result(v)
        real(wp), intent(in) :: rp, r0, rm

        if (amr_tag_abs) then
            v = max(abs(rp - r0), abs(r0 - rm))
        else
            v = abs(rp - rm)/(2._wp*max(abs(r0), 1.e-30_wp))
        end if

    end function f_amr_tag_ind

end module m_amr_nest
