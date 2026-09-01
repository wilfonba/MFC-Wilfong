!>
!! @file m_particle_cloud.fpp
!! @brief Generates particle beds: converts particle_cloud specifications into
!!        individual sphere/circle particle_cloud_ibs entries before reduction.

#:include 'macros.fpp'

!> @brief Generates particle beds by converting particle_cloud patch specifications into individual immersed boundary patches before
!! domain reduction. Each rank runs the same deterministic placement so no MPI broadcast of particle positions is needed.
module m_particle_cloud

    use m_global_parameters
    use m_constants
    use m_mpi_common
    use m_collisions

    implicit none

    private

    public :: s_generate_particle_clouds, s_add_cloud_particle

contains

    !> Generate all particle beds and fill particle_cloud_ibs. Called on all ranks before s_reduce_ib_patch_array. Each packing
    !! method owns and allocates its own per-cloud working array (see s_particle_cloud_lattice / s_particle_cloud_rejection_pack)
    !! and hands back only the entries that fall within this rank's IB neighborhood. Only the first num_particle_cloud_ibs of them
    !! are actually written - callers must use that count, not size(particle_cloud_ibs), since the remainder of the array is left
    !! uninitialized.
    impure subroutine s_generate_particle_clouds(particle_cloud_ibs, num_particle_cloud_ibs)

        type(ib_patch_parameters), allocatable, intent(out), dimension(:) :: particle_cloud_ibs
        integer, intent(out)                                              :: num_particle_cloud_ibs
        type(ib_patch_parameters), allocatable                            :: cloud_ibs(:)
        integer                                                           :: cloud_idx, glbl_idx, num_cloud_ibs, n_total_particles
        real(wp)                                                          :: t_start, t_end

        if (num_particle_clouds == 0) then
            allocate (particle_cloud_ibs(0))
            num_particle_cloud_ibs = 0
            return
        end if

        call cpu_time(t_start)

        n_total_particles = 0
        do cloud_idx = 1, num_particle_clouds
            n_total_particles = n_total_particles + particle_cloud(cloud_idx)%num_particles
        end do
        allocate (particle_cloud_ibs(min(num_ib_patches_max_namelist, n_total_particles)))

        num_particle_cloud_ibs = 0
        glbl_idx = num_ibs

        do cloud_idx = 1, num_particle_clouds
            ! Dispatch on packing method only: rejection sampling and relaxation both handle box and
            ! hemisphere-shell geometries (the per-candidate geometry sampling lives in s_sample_cloud_candidate),
            ! while lattice packing is box-only - the hemisphere-shell + lattice combination is rejected in
            ! case_validator.py.
            select case (particle_cloud(cloud_idx)%packing_method)
            case (1)  ! rejection (random) packing method
                call s_particle_cloud_rejection_pack(cloud_idx, glbl_idx, cloud_ibs, num_cloud_ibs)
            case (2)  ! lattice packing method
                call s_particle_cloud_lattice(cloud_idx, glbl_idx, cloud_ibs, num_cloud_ibs)
            case (3)  ! relaxation (collective rearrangement) packing method
                call s_particle_cloud_relaxation_pack(cloud_idx, glbl_idx, cloud_ibs, num_cloud_ibs)
            case default
                call s_mpi_abort("Particle cloud packing method is not a known packing method of MFC. Exiting.")
            end select

            @:PROHIBIT(num_particle_cloud_ibs + num_cloud_ibs > num_ib_patches_max_namelist, &
                       & "Too many particle-cloud IBs in one rank's neighborhood. Modify case file or increase num_ib_patches_max_namelist.")
            particle_cloud_ibs(num_particle_cloud_ibs + 1:num_particle_cloud_ibs + num_cloud_ibs) = cloud_ibs(1:num_cloud_ibs)
            num_particle_cloud_ibs = num_particle_cloud_ibs + num_cloud_ibs
            deallocate (cloud_ibs)
        end do

        call cpu_time(t_end)
        if (proc_rank == 0) print '(a,i0,a,f0.3,a)', 'Particle beds placed ', glbl_idx - num_ibs, ' particles in ', &
            & t_end - t_start, ' seconds.'

    end subroutine s_generate_particle_clouds

    !> Rejection-samples particle centres into a box or hemisphere-shell region with a minimum centre-to-centre spacing. Rejection
    !! sampling needs every placed particle tracked (regardless of which rank's neighborhood it falls in) to detect overlaps
    !! deterministically, so cloud_ibs is allocated here to the cloud's full requested particle count and only pared down to this
    !! rank's neighborhood afterwards, via s_reduce_particle_cloud_ibs. Only the per-candidate geometry sampling differs between box
    !! and hemisphere shell; it is delegated to s_sample_cloud_candidate, and every other step (overlap rejection via the spatial
    !! hash, acceptance, reduction) is geometry-independent.
    subroutine s_particle_cloud_rejection_pack(cloud_idx, glbl_idx, cloud_ibs, num_cloud_ibs)

        integer, intent(in)                                               :: cloud_idx
        integer, intent(inout)                                            :: glbl_idx
        type(ib_patch_parameters), allocatable, intent(out), dimension(:) :: cloud_ibs
        integer, intent(out)                                              :: num_cloud_ibs
        integer                                                           :: ib_idx, n_placed, geom, seed, alloc_stat
        integer(8)                                                        :: n_attempts, max_attempts
        real(wp)                                                          :: min_dist, rx, ry, rz
        logical                                                           :: overlaps, reject, periodic_pack
        real(wp), allocatable                                             :: placed(:,:)
        integer                                                           :: hash_size, slot
        integer                                                           :: bx, by, bz, nx_bins, ny_bins, nz_bins
        integer, allocatable                                              :: hash_head(:), chain_next(:)
        real(wp)                                                          :: xmin, ymin, zmin, length_x, length_y, length_z

        allocate (cloud_ibs(particle_cloud(cloud_idx)%num_particles), stat=alloc_stat)
        if (alloc_stat /= 0) then
            call s_mpi_abort("Error :: Ran out of CPU memory trying to allocate particle cloud IB array. " &
                             & // "Current system resources cannot perform rejection packing with the specified number of particles.")
        end if
        ib_idx = 0

        min_dist = 2._wp*particle_cloud(cloud_idx)%radius + particle_cloud(cloud_idx)%min_spacing
        periodic_pack = particle_cloud(cloud_idx)%cloud_geometry == 1 .and. particle_cloud(cloud_idx)%periodic == 1
        length_x = particle_cloud(cloud_idx)%length_x
        length_y = particle_cloud(cloud_idx)%length_y
        length_z = particle_cloud(cloud_idx)%length_z
        xmin = particle_cloud(cloud_idx)%x_centroid - 0.5_wp*length_x
        ymin = particle_cloud(cloud_idx)%y_centroid - 0.5_wp*length_y
        zmin = particle_cloud(cloud_idx)%z_centroid - 0.5_wp*length_z
        nx_bins = max(1, ceiling(length_x/min_dist))
        ny_bins = max(1, ceiling(length_y/min_dist))
        nz_bins = max(1, ceiling(length_z/min_dist))
        if (num_dims < 3) nz_bins = 1

        if (num_dims < 3) then
            geom = 2  ! circle for 2D
        else
            geom = 8  ! sphere for 3D
        end if

        max_attempts = int(particle_cloud(cloud_idx)%num_particles, 8)*1000_8
        n_placed = 0
        n_attempts = 0
        seed = particle_cloud(cloud_idx)%seed
        if (seed == 0) seed = 1 + cloud_idx*1013904223

        allocate (placed(3, particle_cloud(cloud_idx)%num_particles))

        ! Hash table: 4x overprovisioned for ~25% load factor, minimum 16 buckets. chain_next(i) links placed particle i to the
        ! previous occupant of its bucket.
        hash_size = max(16, 4*particle_cloud(cloud_idx)%num_particles)
        allocate (hash_head(hash_size))
        allocate (chain_next(particle_cloud(cloud_idx)%num_particles))
        hash_head = -1
        chain_next = -1

        do while (n_placed < particle_cloud(cloud_idx)%num_particles .and. n_attempts < max_attempts)
            n_attempts = n_attempts + 1

            call s_sample_cloud_candidate(cloud_idx, seed, rx, ry, rz, reject)
            if (reject) cycle

            call s_check_cloud_particle_overlap(rx, ry, rz, placed, hash_head, chain_next, hash_size, min_dist, periodic_pack, &
                                                & xmin, ymin, zmin, length_x, length_y, length_z, nx_bins, ny_bins, nz_bins, &
                                                & overlaps, bx, by, bz)

            if (.not. overlaps) then
                n_placed = n_placed + 1
                placed(1, n_placed) = rx
                placed(2, n_placed) = ry
                placed(3, n_placed) = rz

                ! Insert into hash grid as head of bucket chain
                slot = f_bin_hash(bx, by, bz, hash_size)
                chain_next(n_placed) = hash_head(slot)
                hash_head(slot) = n_placed

                glbl_idx = glbl_idx + 1
                call s_add_cloud_particle(cloud_idx, ib_idx, glbl_idx, geom, rx, ry, rz, cloud_ibs)
            end if
        end do

        if (n_placed < particle_cloud(cloud_idx)%num_particles) then
            call s_mpi_abort("Error :: Failed to place all particles in particle bed")
        end if

        deallocate (placed, hash_head, chain_next)

        call s_reduce_particle_cloud_ibs(cloud_ibs, ib_idx)
        num_cloud_ibs = ib_idx

    end subroutine s_particle_cloud_rejection_pack

    !> Draws one rejection-sampling candidate centre (rx, ry, rz) for cloud_idx, advancing seed in place. For box geometry the
    !! candidate is uniform in the box and never rejected. For a hemisphere shell the candidate is uniform in the shell volume - 2D
    !! uses theta uniform on [0, pi] with the sqrt radial CDF; 3D uses uniform phi, uniform cos(polar) on [0, 1], and the cube-root
    !! radial CDF - and reject is set when it lands within one particle radius of the flat face (the plane through the centroid
    !! perpendicular to shell_axis: 1=x, 2=y, 3=z; 2D has no z-axis so any value other than 1 falls back to y, matching the
    !! pre-shell_axis default), a hard geometric cut applied after sampling that preserves uniformity over the remaining region.
    subroutine s_sample_cloud_candidate(cloud_idx, seed, rx, ry, rz, reject)

        integer, intent(in)    :: cloud_idx
        integer, intent(inout) :: seed
        real(wp), intent(out)  :: rx, ry, rz
        logical, intent(out)   :: reject
        real(wp)               :: xmin, xmax, ymin, ymax, zmin, zmax
        real(wp)               :: theta, phi, r_shell, rho, u, zdir, r_inner, r_outer

        reject = .false.

        select case (particle_cloud(cloud_idx)%cloud_geometry)
        case (1)  ! box
            xmin = particle_cloud(cloud_idx)%x_centroid - 0.5_wp*particle_cloud(cloud_idx)%length_x
            xmax = particle_cloud(cloud_idx)%x_centroid + 0.5_wp*particle_cloud(cloud_idx)%length_x
            ymin = particle_cloud(cloud_idx)%y_centroid - 0.5_wp*particle_cloud(cloud_idx)%length_y
            ymax = particle_cloud(cloud_idx)%y_centroid + 0.5_wp*particle_cloud(cloud_idx)%length_y
            zmin = particle_cloud(cloud_idx)%z_centroid - 0.5_wp*particle_cloud(cloud_idx)%length_z
            zmax = particle_cloud(cloud_idx)%z_centroid + 0.5_wp*particle_cloud(cloud_idx)%length_z

            rx = xmin + f_xorshift(seed)*(xmax - xmin)
            ry = ymin + f_xorshift(seed)*(ymax - ymin)
            if (num_dims < 3) then
                rz = particle_cloud(cloud_idx)%z_centroid
            else
                rz = zmin + f_xorshift(seed)*(zmax - zmin)
            end if
        case (2)  ! hemisphere shell
            r_inner = particle_cloud(cloud_idx)%shell_inner_radius + particle_cloud(cloud_idx)%radius
            r_outer = particle_cloud(cloud_idx)%shell_outer_radius - particle_cloud(cloud_idx)%radius

            if (num_dims < 3) then
                theta = pi*f_xorshift(seed)
                u = f_xorshift(seed)
                r_shell = sqrt((r_outer**2 - r_inner**2)*u + r_inner**2)
                rz = particle_cloud(cloud_idx)%z_centroid
                if (particle_cloud(cloud_idx)%shell_axis == 1) then
                    rx = particle_cloud(cloud_idx)%x_centroid + r_shell*sin(theta)
                    ry = particle_cloud(cloud_idx)%y_centroid + r_shell*cos(theta)
                    if (rx < particle_cloud(cloud_idx)%x_centroid + particle_cloud(cloud_idx)%radius) reject = .true.
                else
                    rx = particle_cloud(cloud_idx)%x_centroid + r_shell*cos(theta)
                    ry = particle_cloud(cloud_idx)%y_centroid + r_shell*sin(theta)
                    if (ry < particle_cloud(cloud_idx)%y_centroid + particle_cloud(cloud_idx)%radius) reject = .true.
                end if
            else
                phi = 2._wp*pi*f_xorshift(seed)
                zdir = f_xorshift(seed)
                rho = sqrt(max(0._wp, 1._wp - zdir**2))
                u = f_xorshift(seed)
                r_shell = ((r_outer**3 - r_inner**3)*u + r_inner**3)**(1._wp/3._wp)
                select case (particle_cloud(cloud_idx)%shell_axis)
                case (1)  ! opens toward +x
                    rx = particle_cloud(cloud_idx)%x_centroid + r_shell*zdir
                    ry = particle_cloud(cloud_idx)%y_centroid + r_shell*rho*cos(phi)
                    rz = particle_cloud(cloud_idx)%z_centroid + r_shell*rho*sin(phi)
                    if (rx < particle_cloud(cloud_idx)%x_centroid + particle_cloud(cloud_idx)%radius) reject = .true.
                case (2)  ! opens toward +y
                    ry = particle_cloud(cloud_idx)%y_centroid + r_shell*zdir
                    rx = particle_cloud(cloud_idx)%x_centroid + r_shell*rho*cos(phi)
                    rz = particle_cloud(cloud_idx)%z_centroid + r_shell*rho*sin(phi)
                    if (ry < particle_cloud(cloud_idx)%y_centroid + particle_cloud(cloud_idx)%radius) reject = .true.
                case default  ! 3: opens toward +z
                    rz = particle_cloud(cloud_idx)%z_centroid + r_shell*zdir
                    rx = particle_cloud(cloud_idx)%x_centroid + r_shell*rho*cos(phi)
                    ry = particle_cloud(cloud_idx)%y_centroid + r_shell*rho*sin(phi)
                    if (rz < particle_cloud(cloud_idx)%z_centroid + particle_cloud(cloud_idx)%radius) reject = .true.
                end select
            end if
        case default
            call s_mpi_abort("Particle cloud geometry is not a known cloud geometry of MFC. Exiting.")
        end select

    end subroutine s_sample_cloud_candidate

    !> Places particles on the optimally dense lattice for the cloud region: a triangular lattice in 2D, a face-centered cubic
    !! lattice in 3D. The lattice spacing is set by the particle density (num_particles over the region area/volume); if that
    !! spacing falls below the required centre-to-centre distance (2*radius + min_spacing), the region is too dense and the run is
    !! aborted. No two lattice sites can overlap, so unlike rejection packing each site's IB neighborhood membership
    !! (get_neighbor_bounds() must already have run) is checked as it is generated and only in-neighborhood sites are stored;
    !! cloud_ibs is therefore allocated to the neighborhood-sized cap rather than the cloud's full particle count.
    subroutine s_particle_cloud_lattice(cloud_idx, glbl_idx, cloud_ibs, num_cloud_ibs)

        integer, intent(in)                                               :: cloud_idx
        integer, intent(inout)                                            :: glbl_idx
        type(ib_patch_parameters), allocatable, intent(out), dimension(:) :: cloud_ibs
        integer, intent(out)                                              :: num_cloud_ibs
        integer                                                           :: ib_idx, n_placed, n_target, geom
        integer                                                           :: row, col, ncx, ncy, ix, jy, kz, b
        real(wp)                                                          :: xmin, xmax, ymin, ymax, zmin, zmax, min_dist
        real(wp)                                                          :: spacing, row_dy, cell, x0, px, py
        real(wp), dimension(4)                                            :: bx_off, by_off, bz_off
        real(wp), dimension(3)                                            :: centroid

        allocate (cloud_ibs(min(num_ib_patches_max_namelist, particle_cloud(cloud_idx)%num_particles)))
        ib_idx = 0

        xmin = particle_cloud(cloud_idx)%x_centroid - 0.5_wp*particle_cloud(cloud_idx)%length_x
        xmax = particle_cloud(cloud_idx)%x_centroid + 0.5_wp*particle_cloud(cloud_idx)%length_x
        ymin = particle_cloud(cloud_idx)%y_centroid - 0.5_wp*particle_cloud(cloud_idx)%length_y
        ymax = particle_cloud(cloud_idx)%y_centroid + 0.5_wp*particle_cloud(cloud_idx)%length_y
        zmin = particle_cloud(cloud_idx)%z_centroid - 0.5_wp*particle_cloud(cloud_idx)%length_z
        zmax = particle_cloud(cloud_idx)%z_centroid + 0.5_wp*particle_cloud(cloud_idx)%length_z

        min_dist = 2._wp*particle_cloud(cloud_idx)%radius + particle_cloud(cloud_idx)%min_spacing
        n_target = particle_cloud(cloud_idx)%num_particles
        n_placed = 0

        if (num_dims < 3) then
            geom = 2  ! circle for 2D
            ! Triangular lattice: area per particle = (sqrt(3)/2)*spacing**2.
            spacing = sqrt(2._wp*(xmax - xmin)*(ymax - ymin)/(sqrt(3._wp)*real(n_target, wp)))
        else
            geom = 8  ! sphere for 3D
            ! Face-centered cubic lattice: volume per particle = spacing**3/sqrt(2).
            spacing = (sqrt(2._wp)*(xmax - xmin)*(ymax - ymin)*(zmax - zmin)/real(n_target, wp))**(1._wp/3._wp)
        end if

        if (spacing < min_dist) then
            call s_mpi_abort("Error :: Particle cloud is too dense for lattice packing; " &
                             & // "reduce num_particles or min_spacing, or enlarge the cloud region")
        end if

        if (num_dims < 3) then
            ! Triangular lattice: rows pitched by spacing*sqrt(3)/2, odd rows shifted by half a spacing.
            row_dy = spacing*sqrt(3._wp)/2._wp
            row = 0
            do while (n_placed < n_target)
                py = ymin + real(row, wp)*row_dy
                x0 = xmin
                if (mod(row, 2) == 1) x0 = xmin + 0.5_wp*spacing
                col = 0
                px = x0
                do while (px <= xmax .and. n_placed < n_target)
                    glbl_idx = glbl_idx + 1
                    centroid = [px, py, particle_cloud(cloud_idx)%z_centroid]
                    if (f_neighborhood_ranks_own_location(centroid)) then
                        call s_add_cloud_particle(cloud_idx, ib_idx, glbl_idx, geom, centroid(1), centroid(2), centroid(3), &
                                                  & cloud_ibs)
                    end if
                    n_placed = n_placed + 1
                    col = col + 1
                    px = x0 + real(col, wp)*spacing
                end do
                row = row + 1
            end do
        else
            ! Face-centered cubic lattice via the conventional cubic cell (side = spacing*sqrt(2)) and its four basis points.
            cell = spacing*sqrt(2._wp)
            bx_off = [0._wp, 0.5_wp, 0.5_wp, 0._wp]*cell
            by_off = [0._wp, 0.5_wp, 0._wp, 0.5_wp]*cell
            bz_off = [0._wp, 0._wp, 0.5_wp, 0.5_wp]*cell
            ncx = max(1, ceiling((xmax - xmin)/cell))
            ncy = max(1, ceiling((ymax - ymin)/cell))
            kz = 0
            do while (n_placed < n_target)
                do jy = 0, ncy - 1
                    do ix = 0, ncx - 1
                        do b = 1, 4
                            if (n_placed >= n_target) exit
                            centroid = [xmin + real(ix, wp)*cell + bx_off(b), ymin + real(jy, wp)*cell + by_off(b), &
                                                    & zmin + real(kz, wp)*cell + bz_off(b)]
                            glbl_idx = glbl_idx + 1
                            if (f_neighborhood_ranks_own_location(centroid)) then
                                call s_add_cloud_particle(cloud_idx, ib_idx, glbl_idx, geom, centroid(1), centroid(2), &
                                                          & centroid(3), cloud_ibs)
                            end if
                            n_placed = n_placed + 1
                        end do
                    end do
                end do
                kz = kz + 1
            end do
        end if

        num_cloud_ibs = ib_idx

    end subroutine s_particle_cloud_lattice

    !> Places particles by collective rearrangement: n_target centres are dropped into the cloud region at random, overlaps and all,
    !! and are then relaxed apart until none remain. Rejection packing never moves a particle once placed, so it saturates at the
    !! random sequential addition limit - roughly 0.38 by volume in 3D and 0.55 by area in 2D, and well below either once the region
    !! is only a few particles thick. Relaxation reaches random close packing (~0.63 in 3D) instead, and the bed stays disordered,
    !! unlike the lattice method. Every sweep is a Jacobi step - each overlapping pair contributes a displacement along its line of
    !! centres, the contributions are summed per particle and applied together - so the result does not depend on the order pairs
    !! are visited in, and every rank relaxing the same cloud reaches the same bed.
    subroutine s_particle_cloud_relaxation_pack(cloud_idx, glbl_idx, cloud_ibs, num_cloud_ibs)

        integer, intent(in)                                               :: cloud_idx
        integer, intent(inout)                                            :: glbl_idx
        type(ib_patch_parameters), allocatable, intent(out), dimension(:) :: cloud_ibs
        integer, intent(out)                                              :: num_cloud_ibs

        ! A bed too dense to relax never clears, so the sweep cap is what turns that into an abort rather than a hang. The
        ! tolerance is relative to min_dist; a residual overlap this small is orders of magnitude below any grid spacing the
        ! particles could be resolved on.
        integer, parameter  :: max_sweeps = 10000
        real(wp), parameter :: converged_overlap = 1.e-8_wp
        ! Jacobi damping: a particle wedged between several neighbours receives one displacement per overlapping pair, and
        ! applying their sum undamped overshoots and sets the bed ringing instead of settling.
        real(wp), parameter   :: damping = 0.5_wp
        integer               :: i, ib_idx, n_target, n_placed, geom, seed, sweep, hash_size, alloc_stat
        real(wp)              :: min_dist, rx, ry, rz, max_overlap
        logical               :: reject
        real(wp), allocatable :: placed(:,:), disp(:,:)
        integer, allocatable  :: hash_head(:), chain_next(:), last_seen(:)

        n_target = particle_cloud(cloud_idx)%num_particles
        allocate (cloud_ibs(n_target), stat=alloc_stat)
        if (alloc_stat /= 0) then
            call s_mpi_abort("Error :: Ran out of CPU memory trying to allocate particle cloud IB array. " &
                             & // "Current system resources cannot perform relaxation packing with the specified number of particles.")
        end if
        ib_idx = 0

        ! Relax against a hair more than the required separation. A sweep stops with a residual overlap of up to
        ! converged_overlap of this distance, and the margin is what keeps that residual from eating into the gap the case asked
        ! for, so a relaxed bed clears 2*radius + min_spacing as strictly as a rejection-packed one does.
        min_dist = (2._wp*particle_cloud(cloud_idx)%radius + particle_cloud(cloud_idx)%min_spacing)*(1._wp + converged_overlap)
        geom = merge(2, 8, num_dims < 3)  ! circle for 2D, sphere for 3D

        seed = particle_cloud(cloud_idx)%seed
        if (seed == 0) seed = 1 + cloud_idx*1013904223

        ! Hash table sized as in rejection packing: 4x overprovisioned for a ~25% load factor, minimum 16 buckets.
        hash_size = max(16, 4*n_target)
        allocate (placed(3, n_target), disp(3, n_target))
        allocate (hash_head(hash_size), chain_next(n_target), last_seen(n_target))

        ! Seed the bed by drawing centres from the same sampler rejection packing uses, but keeping every draw: overlaps are what
        ! the relaxation is for, and starting from a uniform draw is what keeps the final bed disordered.
        n_placed = 0
        do while (n_placed < n_target)
            call s_sample_cloud_candidate(cloud_idx, seed, rx, ry, rz, reject)
            if (reject) cycle
            n_placed = n_placed + 1
            placed(1, n_placed) = rx
            placed(2, n_placed) = ry
            placed(3, n_placed) = rz
        end do

        do sweep = 1, max_sweeps
            call s_relax_cloud_overlaps(n_target, placed, disp, hash_head, chain_next, last_seen, hash_size, min_dist, damping, &
                                        & max_overlap)
            do i = 1, n_target
                call s_clamp_to_cloud_region(cloud_idx, placed(1, i), placed(2, i), placed(3, i))
            end do
            if (max_overlap < converged_overlap*min_dist) exit
        end do

        if (max_overlap >= converged_overlap*min_dist) then
            call s_mpi_abort("Error :: Particle cloud is too dense for relaxation packing; " &
                             & // "reduce num_particles or void_fraction or min_spacing, or enlarge the cloud region")
        end if

        do i = 1, n_target
            glbl_idx = glbl_idx + 1
            call s_add_cloud_particle(cloud_idx, ib_idx, glbl_idx, geom, placed(1, i), placed(2, i), placed(3, i), cloud_ibs)
        end do

        deallocate (placed, disp, hash_head, chain_next, last_seen)

        call s_reduce_particle_cloud_ibs(cloud_ibs, ib_idx)
        num_cloud_ibs = ib_idx

    end subroutine s_particle_cloud_relaxation_pack

    !> One Jacobi relaxation sweep: rebuilds the spatial hash from the particles' current positions, sums the separating
    !! displacement every overlapping pair asks for, and applies the damped total. Hands back the largest overlap it saw so the
    !! caller can stop once the bed is clear. Only neighbours with a higher index contribute, which visits each pair exactly once
    !! and keeps its two halves exactly antisymmetric. last_seen stamps the particles already paired with i this pass: distinct bins
    !! can hash to one slot, so without it a colliding slot would be walked twice and its particles counted twice.
    subroutine s_relax_cloud_overlaps(n_target, placed, disp, hash_head, chain_next, last_seen, hash_size, min_dist, damping, &
                                      & max_overlap)

        integer, intent(in)                     :: n_target, hash_size
        real(wp), intent(inout), dimension(:,:) :: placed
        real(wp), intent(out), dimension(:,:)   :: disp
        integer, intent(out), dimension(:)      :: hash_head, chain_next, last_seen
        real(wp), intent(in)                    :: min_dist, damping
        real(wp), intent(out)                   :: max_overlap
        integer                                 :: i, j, slot, bx, by, bz, dx_b, dy_b, dz_b, dz_lo, dz_hi
        real(wp)                                :: dx, dy, dz, dist, overlap, push
        real(wp), dimension(3)                  :: sep_dir

        hash_head = -1
        chain_next = -1
        last_seen = 0
        do i = 1, n_target
            call s_get_cloud_bin(placed(1, i), placed(2, i), placed(3, i), min_dist, .false., 0._wp, 0._wp, 0._wp, 1, 1, 1, bx, &
                                 & by, bz)
            slot = f_bin_hash(bx, by, bz, hash_size)
            chain_next(i) = hash_head(slot)
            hash_head(slot) = i
        end do

        dz_lo = -1
        dz_hi = 1
        if (num_dims < 3) then
            dz_lo = 0
            dz_hi = 0
        end if

        disp = 0._wp
        max_overlap = 0._wp

        do i = 1, n_target
            call s_get_cloud_bin(placed(1, i), placed(2, i), placed(3, i), min_dist, .false., 0._wp, 0._wp, 0._wp, 1, 1, 1, bx, &
                                 & by, bz)
            do dx_b = -1, 1
                do dy_b = -1, 1
                    do dz_b = dz_lo, dz_hi
                        slot = f_bin_hash(bx + dx_b, by + dy_b, bz + dz_b, hash_size)
                        j = hash_head(slot)
                        do while (j > 0)
                            if (j > i .and. last_seen(j) /= i) then
                                last_seen(j) = i
                                dx = placed(1, i) - placed(1, j)
                                dy = placed(2, i) - placed(2, j)
                                dz = 0._wp
                                if (num_dims == 3) dz = placed(3, i) - placed(3, j)
                                dist = sqrt(dx**2 + dy**2 + dz**2)
                                if (dist < min_dist) then
                                    overlap = min_dist - dist
                                    max_overlap = max(max_overlap, overlap)
                                    ! Coincident centres have no line of centres to separate along; x serves, and the pair still
                                    ! resolves because i and j take opposite signs.
                                    if (dist > sqrt(tiny(1._wp))) then
                                        sep_dir = [dx, dy, dz]/dist
                                    else
                                        sep_dir = [1._wp, 0._wp, 0._wp]
                                    end if
                                    push = 0.5_wp*overlap
                                    disp(:,i) = disp(:,i) + push*sep_dir
                                    disp(:,j) = disp(:,j) - push*sep_dir
                                end if
                            end if
                            j = chain_next(j)
                        end do
                    end do
                end do
            end do
        end do

        do i = 1, n_target
            placed(:,i) = placed(:,i) + damping*disp(:,i)
        end do

    end subroutine s_relax_cloud_overlaps

    !> Pushes a relaxed particle centre back inside its cloud region. The bounds are the ones s_sample_cloud_candidate draws from,
    !! so a relaxed bed fills exactly the region a rejection-packed one would: box centres stay within the box, and shell centres
    !! stay between shell_inner_radius + radius and shell_outer_radius - radius, one particle radius clear of the flat face. The
    !! shell is clamped in the (axial, transverse) half-plane of its opening axis - the axial component is clamped first and then
    !! held fixed while the transverse distance is clamped into the annulus band it leaves, so restoring the radius can never undo
    !! the flat-face clearance.
    subroutine s_clamp_to_cloud_region(cloud_idx, px, py, pz)

        integer, intent(in)     :: cloud_idx
        real(wp), intent(inout) :: px, py, pz
        integer                 :: axis, t1, t2
        real(wp)                :: r_inner, r_outer, axial, trans, trans_min, trans_max
        real(wp), dimension(3)  :: pos, centroid, half_length

        centroid = [particle_cloud(cloud_idx)%x_centroid, particle_cloud(cloud_idx)%y_centroid, &
                                   & particle_cloud(cloud_idx)%z_centroid]
        pos = [px, py, pz]

        select case (particle_cloud(cloud_idx)%cloud_geometry)
        case (1)  ! box
            half_length = 0.5_wp*[particle_cloud(cloud_idx)%length_x, particle_cloud(cloud_idx)%length_y, &
                                                 & particle_cloud(cloud_idx)%length_z]
            pos(1:num_dims) = min(max(pos(1:num_dims), centroid(1:num_dims) - half_length(1:num_dims)), &
                & centroid(1:num_dims) + half_length(1:num_dims))
        case (2)  ! hemisphere shell
            r_inner = particle_cloud(cloud_idx)%shell_inner_radius + particle_cloud(cloud_idx)%radius
            r_outer = particle_cloud(cloud_idx)%shell_outer_radius - particle_cloud(cloud_idx)%radius

            ! The axis the shell opens toward; 2D has no z-axis, so anything other than x opens toward y, matching the sampler.
            axis = particle_cloud(cloud_idx)%shell_axis
            if (num_dims < 3) then
                if (axis /= 1) axis = 2
                t1 = 3 - axis  ! the other in-plane axis
                t2 = 3  ! 2D has no z-extent, so this contributes nothing to the transverse distance
            else
                t1 = 1 + mod(axis, 3)
                t2 = 1 + mod(axis + 1, 3)
            end if

            axial = min(max(pos(axis) - centroid(axis), particle_cloud(cloud_idx)%radius), r_outer)
            trans = sqrt((pos(t1) - centroid(t1))**2 + (pos(t2) - centroid(t2))**2)
            trans_max = sqrt(max(0._wp, r_outer**2 - axial**2))
            trans_min = sqrt(max(0._wp, r_inner**2 - axial**2))

            pos(axis) = centroid(axis) + axial
            if (trans > sqrt(tiny(1._wp))) then
                pos(t1) = centroid(t1) + (pos(t1) - centroid(t1))*min(max(trans, trans_min), trans_max)/trans
                pos(t2) = centroid(t2) + (pos(t2) - centroid(t2))*min(max(trans, trans_min), trans_max)/trans
            else
                ! On the axis with an inner radius to clear, there is no transverse direction to scale; t1 serves.
                pos(t1) = centroid(t1) + trans_min
                pos(t2) = centroid(t2)
            end if
        end select

        px = pos(1)
        py = pos(2)
        if (num_dims == 3) pz = pos(3)

    end subroutine s_clamp_to_cloud_region

    !> Writes a single placed particle into particle_cloud_ibs at the next free slot, advancing ib_idx. The caller decides whether
    !! this particle belongs in the array (neighborhood membership, for lattice packing, or unconditionally for rejection packing -
    !! see s_particle_cloud_lattice / s_particle_cloud_rejection_pack) and supplies its already-assigned, absolute global patch id
    !! via glbl_idx - s_reduce_ib_patch_array copies gbl_patch_id as-is. Shared by all packing methods so the per-particle
    !! ib_patch_parameters setup stays in one place.
    subroutine s_add_cloud_particle(cloud_idx, ib_idx, glbl_idx, geom, px, py, pz, particle_cloud_ibs)

        integer, intent(in)                                    :: cloud_idx, glbl_idx, geom
        integer, intent(inout)                                 :: ib_idx
        real(wp), intent(in)                                   :: px, py, pz
        type(ib_patch_parameters), intent(inout), dimension(:) :: particle_cloud_ibs

        ib_idx = ib_idx + 1
        @:PROHIBIT(ib_idx > size(particle_cloud_ibs), &
                   & "Too many particle-cloud IBs in one rank's neighborhood. Modify case file or increase num_ib_patches_max_namelist.")

        particle_cloud_ibs(ib_idx)%gbl_patch_id = glbl_idx
        particle_cloud_ibs(ib_idx)%geometry = geom
        particle_cloud_ibs(ib_idx)%x_centroid = px
        particle_cloud_ibs(ib_idx)%y_centroid = py
        particle_cloud_ibs(ib_idx)%z_centroid = pz
        particle_cloud_ibs(ib_idx)%step_x_centroid = px
        particle_cloud_ibs(ib_idx)%step_y_centroid = py
        particle_cloud_ibs(ib_idx)%step_z_centroid = pz
        particle_cloud_ibs(ib_idx)%angles(:) = 0._wp
        particle_cloud_ibs(ib_idx)%step_angles(:) = 0._wp
        particle_cloud_ibs(ib_idx)%vel(:) = 0._wp
        particle_cloud_ibs(ib_idx)%step_vel(:) = 0._wp
        particle_cloud_ibs(ib_idx)%angular_vel(:) = 0._wp
        particle_cloud_ibs(ib_idx)%step_angular_vel(:) = 0._wp
        particle_cloud_ibs(ib_idx)%force(:) = 0._wp
        particle_cloud_ibs(ib_idx)%torque(:) = 0._wp
        particle_cloud_ibs(ib_idx)%centroid_offset(:) = 0._wp
        particle_cloud_ibs(ib_idx)%rotation_matrix = 0._wp
        particle_cloud_ibs(ib_idx)%rotation_matrix(1, 1) = 1._wp
        particle_cloud_ibs(ib_idx)%rotation_matrix(2, 2) = 1._wp
        particle_cloud_ibs(ib_idx)%rotation_matrix(3, 3) = 1._wp
        particle_cloud_ibs(ib_idx)%rotation_matrix_inverse = particle_cloud_ibs(ib_idx)%rotation_matrix
        particle_cloud_ibs(ib_idx)%radius = particle_cloud(cloud_idx)%radius
        particle_cloud_ibs(ib_idx)%mass = particle_cloud(cloud_idx)%mass
        particle_cloud_ibs(ib_idx)%moment = dflt_real
        particle_cloud_ibs(ib_idx)%moving_ibm = particle_cloud(cloud_idx)%moving_ibm
        particle_cloud_ibs(ib_idx)%slip = .false.

        ! Particles are inert surfaces. These must be set explicitly: particle_cloud_ibs is
        ! allocated (not default-initialized) and s_reduce_ib_patch_array copies the whole
        ! struct into patch_ib, overwriting the defaults from
        ! s_assign_default_values_to_user_inputs -- so anything left unset here reaches the
        ! solver as uninitialized memory (a nonzero v_blow injects a garbage wall-normal
        ! velocity and NaNs the field).
        particle_cloud_ibs(ib_idx)%v_blow = 0._wp
        particle_cloud_ibs(ib_idx)%inj_species = 0
        particle_cloud_ibs(ib_idx)%burn_rate_exp = 0._wp
        particle_cloud_ibs(ib_idx)%burn_rate_pref = 0._wp

    end subroutine s_add_cloud_particle

    !> Convert a candidate particle centre to spatial-hash bin coordinates.
    subroutine s_get_cloud_bin(px, py, pz, min_dist, periodic_pack, xmin, ymin, zmin, nx_bins, ny_bins, nz_bins, bx, by, bz)

        real(wp), intent(in) :: px, py, pz, min_dist
        logical, intent(in)  :: periodic_pack
        real(wp), intent(in) :: xmin, ymin, zmin
        integer, intent(in)  :: nx_bins, ny_bins, nz_bins
        integer, intent(out) :: bx, by, bz

        if (periodic_pack) then
            bx = modulo(int(floor((px - xmin)/min_dist)), nx_bins)
            by = modulo(int(floor((py - ymin)/min_dist)), ny_bins)
            if (num_dims < 3) then
                bz = 0
            else
                bz = modulo(int(floor((pz - zmin)/min_dist)), nz_bins)
            end if
        else
            bx = int(floor(px/min_dist))
            by = int(floor(py/min_dist))
            if (num_dims < 3) then
                bz = 0
            else
                bz = int(floor(pz/min_dist))
            end if
        end if

    end subroutine s_get_cloud_bin

    !> Check whether a candidate particle centre overlaps any already-placed particle in neighbouring spatial-hash bins. Scans the
    !! 3x3(x3) bin neighborhood - O(1) average via hash lookup - and also hands back the candidate's own bin so the caller can
    !! insert it without recomputing.
    subroutine s_check_cloud_particle_overlap(px, py, pz, placed, hash_head, chain_next, hash_size, min_dist, periodic_pack, &
        & xmin, ymin, zmin, length_x, length_y, length_z, nx_bins, ny_bins, nz_bins, overlaps, bx, by, bz)

        real(wp), intent(in)                 :: px, py, pz, min_dist
        real(wp), intent(in), dimension(:,:) :: placed
        integer, intent(in), dimension(:)    :: hash_head, chain_next
        integer, intent(in)                  :: hash_size
        logical, intent(in)                  :: periodic_pack
        real(wp), intent(in)                 :: xmin, ymin, zmin, length_x, length_y, length_z
        integer, intent(in)                  :: nx_bins, ny_bins, nz_bins
        logical, intent(out)                 :: overlaps
        integer, intent(out)                 :: bx, by, bz
        integer                              :: nbx, nby, nbz, slot
        integer                              :: dx_b, dy_b, dz_b, dz_lo, dz_hi, j
        real(wp)                             :: dist_sq, min_dist_sq, dx, dy, dz

        call s_get_cloud_bin(px, py, pz, min_dist, periodic_pack, xmin, ymin, zmin, nx_bins, ny_bins, nz_bins, bx, by, bz)

        dz_lo = -1
        dz_hi = 1
        if (num_dims < 3) then
            dz_lo = 0
            dz_hi = 0
        end if

        min_dist_sq = min_dist**2
        overlaps = .false.

        do dx_b = -1, 1
            do dy_b = -1, 1
                do dz_b = dz_lo, dz_hi
                    nbx = bx + dx_b
                    nby = by + dy_b
                    nbz = bz + dz_b
                    if (periodic_pack) then
                        nbx = modulo(nbx, nx_bins)
                        nby = modulo(nby, ny_bins)
                        if (num_dims == 3) nbz = modulo(nbz, nz_bins)
                    end if
                    slot = f_bin_hash(nbx, nby, nbz, hash_size)
                    j = hash_head(slot)
                    do while (j > 0)
                        dx = abs(px - placed(1, j))
                        dy = abs(py - placed(2, j))
                        if (periodic_pack) then
                            dx = min(dx, length_x - dx)
                            dy = min(dy, length_y - dy)
                        end if
                        if (num_dims < 3) then
                            dist_sq = dx**2 + dy**2
                        else
                            dz = abs(pz - placed(3, j))
                            if (periodic_pack) dz = min(dz, length_z - dz)
                            dist_sq = dx**2 + dy**2 + dz**2
                        end if
                        if (dist_sq < min_dist_sq) then
                            overlaps = .true.
                            return
                        end if
                        j = chain_next(j)
                    end do
                end do
            end do
        end do

    end subroutine s_check_cloud_particle_overlap

    !> Compacts cloud_ibs(1:num_ibs) in place, discarding entries outside this rank's IB neighborhood (get_neighbor_bounds() must
    !! already have run) and updating num_ibs to the retained count. Used by rejection packing, which cannot filter as it places
    !! particles (see s_particle_cloud_rejection_pack), to pare its full, unfiltered placement down to this rank's neighborhood.
    subroutine s_reduce_particle_cloud_ibs(cloud_ibs, num_cloud_ibs)

        type(ib_patch_parameters), intent(inout), dimension(:) :: cloud_ibs
        integer, intent(inout)                                 :: num_cloud_ibs
        integer                                                :: i, write_idx
        real(wp), dimension(3)                                 :: centroid

        write_idx = 0
        do i = 1, num_cloud_ibs
            centroid = [cloud_ibs(i)%x_centroid, cloud_ibs(i)%y_centroid, 0._wp]
            if (num_dims == 3) centroid(3) = cloud_ibs(i)%z_centroid
            if (f_neighborhood_ranks_own_location(centroid)) then
                write_idx = write_idx + 1
                if (write_idx /= i) cloud_ibs(write_idx) = cloud_ibs(i)
            end if
        end do
        num_cloud_ibs = write_idx

    end subroutine s_reduce_particle_cloud_ibs

    !> Xorshift PRNG. Advances seed in-place and returns a value in [0, 1).
    function f_xorshift(seed) result(rval)

        integer, intent(inout) :: seed
        real(wp)               :: rval

        seed = ieor(seed, ishft(seed, 13))
        seed = ieor(seed, ishft(seed, -17))
        seed = ieor(seed, ishft(seed, 5))

        ! Mask off the sign bit rather than abs(): at seed = -huge-1, abs() overflows and returns the value
        ! unchanged, yielding rval = -1 and, in the shell sampler, a NaN centre that then passes every overlap
        ! comparison (all comparisons against NaN are false) and gets placed.
        rval = real(iand(seed, huge(seed)), wp)/real(huge(seed), wp)

    end function f_xorshift

    !> Hash bin coordinates to a 1-indexed slot in [1, hash_size]. Uses large prime multipliers to spread bins across buckets. Hash
    !! collisions are benign: the distance check catches false neighbours.
    function f_bin_hash(bx, by, bz, hash_size) result(slot)

        integer, intent(in) :: bx, by, bz, hash_size
        integer             :: slot
        integer(8)          :: key

        key = ieor(ieor(int(bx, 8)*73856093_8, int(by, 8)*19349663_8), int(bz, 8)*83492791_8)
        slot = int(mod(abs(key), int(hash_size, 8))) + 1

    end function f_bin_hash

end module m_particle_cloud
