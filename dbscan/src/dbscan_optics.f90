module dbscan_optics
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_positive_inf, ieee_value
    use dbscan_kinds, only : dp
    use dbscan_types, only : optics_result, frnn_result
    use dbscan_nn, only : frnn, knn_dist
    implicit none
    private

    public :: optics, extract_dbscan, extract_xi

contains

    subroutine optics(x, min_pts, result, eps, status)
        real(dp), intent(in) :: x(:, :) !! Data matrix with observations in rows.
        integer, intent(in) :: min_pts !! OPTICS MinPts, including the point itself.
        type(optics_result), intent(out) :: result !! OPTICS ordering and reachability/core-distance arrays.
        real(dp), intent(in), optional :: eps !! Maximum neighborhood radius; defaults to max minPts-NN distance.
        integer, intent(out), optional :: status !! Zero on success; nonzero for invalid input.
        type(frnn_result) :: nn
        real(dp), allocatable :: kd(:)
        real(dp) :: eps_use
        integer :: stat

        if (min_pts < 2 .or. min_pts >= size(x, 1)) then
            if (present(status)) status = 1
            call empty_result(result)
            return
        end if
        if (present(eps)) then
            eps_use = eps
        else
            kd = knn_dist(x, min_pts, stat)
            if (stat /= 0) then
                if (present(status)) status = stat
                call empty_result(result)
                return
            end if
            eps_use = maxval(kd)
        end if
        call frnn(x, eps_use, nn, stat)
        if (stat /= 0) then
            if (present(status)) status = stat
            call empty_result(result)
            return
        end if
        call optics_from_frnn(nn, min_pts, result)
        result%eps = eps_use
        if (present(status)) status = 0
    end subroutine optics

    subroutine extract_dbscan(object, eps_cl)
        type(optics_result), intent(inout) :: object !! Existing OPTICS ordering modified with flat DBSCAN-like labels.
        real(dp), intent(in) :: eps_cl !! Reachability/core-distance cut threshold.
        integer, allocatable :: ordered_cluster(:)
        integer :: cluster_id
        integer :: i
        integer :: n
        integer :: p

        n = size(object%order)
        if (allocated(object%cluster)) deallocate(object%cluster)
        allocate(object%cluster(n), ordered_cluster(n))
        ordered_cluster = 0
        cluster_id = 0
        do i = 1, n
            p = object%order(i)
            if (object%reachdist(p) > eps_cl) then
                if (object%coredist(p) <= eps_cl) then
                    cluster_id = cluster_id + 1
                    ordered_cluster(i) = cluster_id
                else
                    ordered_cluster(i) = 0
                end if
            else
                ordered_cluster(i) = cluster_id
            end if
        end do
        object%cluster = 0
        do i = 1, n
            object%cluster(object%order(i)) = ordered_cluster(i)
        end do
        object%eps_cl = eps_cl
        object%xi = -1.0_dp
    end subroutine extract_dbscan

    subroutine extract_xi(object, xi, minimum, correct_predecessors, status)
        type(optics_result), intent(inout) :: object !! OPTICS result modified with Xi hierarchy labels and intervals.
        real(dp), intent(in) :: xi !! Xi steepness in the open interval (0,1).
        logical, intent(in), optional :: minimum !! Return only non-overlapping local clusters when true.
        logical, intent(in), optional :: correct_predecessors !! Apply predecessor correction; default true.
        integer, intent(out), optional :: status !! Zero on success; nonzero for invalid xi.
        real(dp), allocatable :: rd(:)
        real(dp), allocatable :: sda_max(:)
        real(dp), allocatable :: sda_mib(:)
        integer, allocatable :: sda_s(:)
        integer, allocatable :: sda_e(:)
        integer, allocatable :: cs(:)
        integer, allocatable :: ce(:)
        integer, allocatable :: corder(:)
        integer :: nsda
        integer :: ncl
        integer :: n
        integer :: index
        integer :: startsteep
        integer :: endsteep
        integer :: cstart
        integer :: cend
        integer :: a
        integer :: b
        integer :: i
        integer :: j
        integer :: p
        real(dp) :: mib
        real(dp) :: startval
        real(dp) :: esuccr
        real(dp) :: ixi
        logical :: use_minimum
        logical :: use_correction
        logical :: predecessor_found

        if (present(status)) status = 0
        if (xi <= 0.0_dp .or. xi >= 1.0_dp) then
            if (present(status)) status = 1
            return
        end if
        use_minimum = .false.
        if (present(minimum)) use_minimum = minimum
        use_correction = .true.
        if (present(correct_predecessors)) use_correction = correct_predecessors

        n = size(object%order)
        allocate(rd(n), sda_max(n), sda_mib(n), sda_s(n), sda_e(n))
        allocate(cs(max(1, n * n)), ce(max(1, n * n)))
        rd = object%reachdist(object%order)
        ixi = 1.0_dp - xi
        nsda = 0
        ncl = 0
        mib = 0.0_dp
        index = 1

        do while (index <= n)
            mib = max(mib, rd(index))
            if (index + 1 > n) exit

            if (steep_down(rd, index, ixi)) then
                call filter_sda(mib, ixi, nsda, sda_s, sda_e, sda_max, sda_mib)
                startval = rd(index)
                mib = 0.0_dp
                startsteep = index
                endsteep = index + 1
                do while (index + 1 <= n)
                    index = index + 1
                    if (steep_down(rd, index, ixi)) then
                        endsteep = index + 1
                        cycle
                    end if
                    if (.not. steep_down(rd, index, 1.0_dp)) exit
                    if (index - endsteep > object%min_pts) exit
                end do
                nsda = nsda + 1
                sda_s(nsda) = startsteep
                sda_e(nsda) = min(endsteep, n)
                sda_max(nsda) = startval
                sda_mib(nsda) = 0.0_dp
                cycle
            end if

            if (steep_up(rd, index, ixi)) then
                call filter_sda(mib, ixi, nsda, sda_s, sda_e, sda_max, sda_mib)
                startsteep = index
                endsteep = index + 1
                mib = rd(index)
                if (index + 1 > n) then
                    esuccr = ieee_value(1.0_dp, ieee_positive_inf)
                else
                    esuccr = rd(index + 1)
                end if
                if (ieee_is_finite(esuccr)) then
                    do while (index + 1 <= n)
                        index = index + 1
                        if (steep_up(rd, index, ixi)) then
                            endsteep = index + 1
                            mib = rd(index)
                            if (index + 1 <= n) then
                                esuccr = rd(index + 1)
                            else
                                esuccr = ieee_value(1.0_dp, ieee_positive_inf)
                            end if
                            if (.not. ieee_is_finite(esuccr)) then
                                endsteep = endsteep - 1
                                exit
                            end if
                            cycle
                        end if
                        if (.not. steep_up(rd, index, 1.0_dp)) exit
                        if (index - endsteep > object%min_pts) exit
                    end do
                else
                    endsteep = endsteep - 1
                    index = index + 1
                end if

                do a = nsda, 1, -1
                    if (mib * ixi < sda_mib(a)) cycle
                    cstart = sda_s(a)
                    cend = min(endsteep, n)
                    if (use_correction) then
                        do while (cend > cstart .and. .not. ieee_is_finite(rd(cend)))
                            cend = cend - 1
                        end do
                    end if

                    if (sda_max(a) * ixi >= esuccr) then
                        do while (cstart < cend .and. rd(cstart + 1) > esuccr)
                            cstart = cstart + 1
                        end do
                    else if (esuccr * ixi >= sda_max(a)) then
                        do while (cend > cstart .and. rd(cend - 1) > sda_max(a))
                            cend = cend - 1
                        end do
                    end if

                    if (use_correction) then
                        do while (cend > cstart)
                            p = object%predecessor(object%order(cend))
                            predecessor_found = .false.
                            if (p > 0) then
                                do b = cstart, cend - 1
                                    if (object%order(b) == p) then
                                        predecessor_found = .true.
                                        exit
                                    end if
                                end do
                            end if
                            if (predecessor_found) exit
                            cend = cend - 1
                        end do
                    end if
                    if (index > 1) then
                        if (steep_up(rd, index - 1, ixi)) cend = cend - 1
                    end if
                    if (cend - cstart + 1 < object%min_pts) cycle
                    ncl = ncl + 1
                    cs(ncl) = cstart
                    ce(ncl) = cend
                end do
            else
                index = index + 1
            end if
        end do

        object%xi = xi
        object%eps_cl = -1.0_dp
        if (allocated(object%xi_start)) deallocate(object%xi_start)
        if (allocated(object%xi_end)) deallocate(object%xi_end)
        if (allocated(object%cluster)) deallocate(object%cluster)
        allocate(object%cluster(n))
        object%cluster = 0
        if (ncl == 0) then
            allocate(object%xi_start(0), object%xi_end(0))
            return
        end if

        call sort_cluster_bounds(cs, ce, ncl)
        allocate(object%xi_start(ncl), object%xi_end(ncl))
        object%xi_start = cs(1:ncl)
        object%xi_end = ce(1:ncl)
        allocate(corder(ncl))
        do i = 1, ncl
            corder(i) = i
        end do
        call sort_cluster_sizes(cs, ce, corder, ncl, use_minimum)
        do j = 1, ncl
            i = corder(j)
            if (use_minimum) then
                if (any(object%cluster(object%order(cs(i):ce(i))) /= 0)) cycle
            end if
            object%cluster(object%order(cs(i):ce(i))) = i
        end do
    end subroutine extract_xi

    subroutine optics_from_frnn(nn, min_pts, result)
        type(frnn_result), intent(in) :: nn !! Fixed-radius neighbor lists excluding self.
        integer, intent(in) :: min_pts !! MinPts including the point itself.
        type(optics_result), intent(out) :: result !! OPTICS ordering and distance summaries.
        logical, allocatable :: visited(:)
        logical, allocatable :: in_seed(:)
        integer, allocatable :: seeds(:)
        real(dp), allocatable :: reach(:)
        real(dp), allocatable :: core(:)
        integer, allocatable :: pre(:)
        integer :: p
        integer :: q
        integer :: root
        integer :: n
        integer :: nseed
        integer :: norder
        real(dp) :: inf

        n = nn%n
        inf = ieee_value(1.0_dp, ieee_positive_inf)
        allocate(visited(n), in_seed(n), seeds(n), reach(n), core(n), pre(n))
        allocate(result%order(n), result%reachdist(n), result%coredist(n), result%predecessor(n))
        visited = .false.
        in_seed = .false.
        reach = inf
        core = inf
        pre = 0
        norder = 0
        nseed = 0

        do p = 1, n
            if (visited(p)) cycle
            root = p
            call visit_point(p, root, nn, min_pts, visited, in_seed, seeds, nseed, reach, core, pre, &
                result%order, norder)
            do while (nseed > 0)
                q = pop_seed(seeds, nseed, in_seed, reach)
                call visit_point(q, root, nn, min_pts, visited, in_seed, seeds, nseed, reach, core, pre, &
                    result%order, norder)
            end do
        end do
        result%reachdist = reach
        result%coredist = core
        result%predecessor = pre
        result%min_pts = min_pts
        result%eps_cl = -1.0_dp
        result%xi = -1.0_dp
    end subroutine optics_from_frnn

    subroutine visit_point(p, root, nn, min_pts, visited, in_seed, seeds, nseed, reach, core, pre, order, norder)
        integer, intent(in) :: p !! Point to mark visited and expand.
        integer, intent(in) :: root !! Root of the current expansion component.
        type(frnn_result), intent(in) :: nn !! Fixed-radius neighbors excluding self.
        integer, intent(in) :: min_pts !! MinPts including self.
        logical, intent(inout) :: visited(:) !! Visitation state.
        logical, intent(inout) :: in_seed(:) !! Seed-membership state.
        integer, intent(inout) :: seeds(:) !! Updateable-priority-queue storage.
        integer, intent(inout) :: nseed !! Number of active seed entries.
        real(dp), intent(inout) :: reach(:) !! Reachability distances.
        real(dp), intent(inout) :: core(:) !! Core distances.
        integer, intent(inout) :: pre(:) !! One-based predecessor indices; zero means undefined.
        integer, intent(inout) :: order(:) !! OPTICS visitation order.
        integer, intent(inout) :: norder !! Number of entries currently stored in order.
        integer :: count_n

        if (visited(p)) return
        visited(p) = .true.
        count_n = nn%offset(p + 1) - nn%offset(p)
        if (count_n + 1 >= min_pts) core(p) = kth_neighbor_distance(nn, p, min_pts - 1)
        if (pre(p) == 0 .and. p /= root) pre(p) = root
        norder = norder + 1
        order(norder) = p
        if (.not. ieee_is_finite(core(p))) return
        call update_neighbors(p, nn, visited, in_seed, seeds, nseed, reach, core, pre)
    end subroutine visit_point

    subroutine update_neighbors(p, nn, visited, in_seed, seeds, nseed, reach, core, pre)
        integer, intent(in) :: p !! Core point used to update neighboring reachability distances.
        type(frnn_result), intent(in) :: nn !! Fixed-radius neighbor lists.
        logical, intent(in) :: visited(:) !! Visitation state.
        logical, intent(inout) :: in_seed(:) !! Seed-membership state.
        integer, intent(inout) :: seeds(:) !! Active seed identifiers.
        integer, intent(inout) :: nseed !! Number of active seeds.
        real(dp), intent(inout) :: reach(:) !! Current reachability distances.
        real(dp), intent(in) :: core(:) !! Current core distances.
        integer, intent(inout) :: pre(:) !! Predecessor identifiers.
        integer :: idx
        integer :: o
        real(dp) :: new_reach

        do idx = nn%offset(p), nn%offset(p + 1) - 1
            o = nn%id(idx)
            if (visited(o)) cycle
            new_reach = max(core(p), nn%dist(idx))
            if (.not. ieee_is_finite(reach(o))) then
                reach(o) = new_reach
                if (.not. in_seed(o)) then
                    nseed = nseed + 1
                    seeds(nseed) = o
                    in_seed(o) = .true.
                end if
            else if (new_reach < reach(o)) then
                reach(o) = new_reach
                pre(o) = p
            end if
        end do
    end subroutine update_neighbors

    integer function pop_seed(seeds, nseed, in_seed, reach) result(q)
        integer, intent(inout) :: seeds(:) !! Active seed identifiers.
        integer, intent(inout) :: nseed !! Number of active seed identifiers.
        logical, intent(inout) :: in_seed(:) !! Seed-membership state.
        real(dp), intent(in) :: reach(:) !! Reachability distances used as priorities.
        integer :: best
        integer :: i

        best = 1
        do i = 2, nseed
            if (reach(seeds(i)) < reach(seeds(best))) then
                best = i
            else if (.not. (reach(seeds(i)) > reach(seeds(best)))) then
                if (seeds(i) > seeds(best)) best = i
            end if
        end do
        q = seeds(best)
        in_seed(q) = .false.
        seeds(best) = seeds(nseed)
        nseed = nseed - 1
    end function pop_seed

    pure real(dp) function kth_neighbor_distance(nn, row, k) result(value)
        type(frnn_result), intent(in) :: nn !! Sorted fixed-radius neighbors.
        integer, intent(in) :: row !! One-based row index.
        integer, intent(in) :: k !! One-based neighbor rank excluding self.
        value = nn%dist(nn%offset(row) + k - 1)
    end function kth_neighbor_distance

    pure logical function steep_up(rd, i, ixi) result(value)
        real(dp), intent(in) :: rd(:) !! Reachability distances in OPTICS order.
        integer, intent(in) :: i !! One-based position in rd.
        real(dp), intent(in) :: ixi !! Multiplicative steepness factor.

        if (i < 1 .or. i > size(rd)) then
            value = .false.
        else if (.not. ieee_is_finite(rd(i))) then
            value = .false.
        else if (i + 1 > size(rd)) then
            value = .true.
        else
            value = rd(i) <= rd(i + 1) * ixi
        end if
    end function steep_up

    pure logical function steep_down(rd, i, ixi) result(value)
        real(dp), intent(in) :: rd(:) !! Reachability distances in OPTICS order.
        integer, intent(in) :: i !! One-based position in rd.
        real(dp), intent(in) :: ixi !! Multiplicative steepness factor.

        if (i < 1 .or. i + 1 > size(rd)) then
            value = .false.
        else if (.not. ieee_is_finite(rd(i + 1))) then
            value = .false.
        else
            value = rd(i) * ixi >= rd(i + 1)
        end if
    end function steep_down

    subroutine filter_sda(mib, ixi, nsda, sda_s, sda_e, sda_max, sda_mib)
        real(dp), intent(in) :: mib !! Current maximum-in-between reachability distance.
        real(dp), intent(in) :: ixi !! Xi multiplicative factor.
        integer, intent(inout) :: nsda !! Number of retained steep-down areas.
        integer, intent(inout) :: sda_s(:) !! Start positions of steep-down areas.
        integer, intent(inout) :: sda_e(:) !! End positions of steep-down areas.
        real(dp), intent(inout) :: sda_max(:) !! Starting maxima of steep-down areas.
        real(dp), intent(inout) :: sda_mib(:) !! Maximum-in-between values for steep-down areas.
        integer :: i
        integer :: keep

        keep = 0
        do i = 1, nsda
            if (sda_max(i) * ixi <= mib) cycle
            keep = keep + 1
            sda_s(keep) = sda_s(i)
            sda_e(keep) = sda_e(i)
            sda_max(keep) = sda_max(i)
            sda_mib(keep) = max(sda_mib(i), mib)
        end do
        nsda = keep
    end subroutine filter_sda

    subroutine sort_cluster_bounds(cs, ce, ncl)
        integer, intent(inout) :: cs(:) !! Cluster start positions to sort ascending.
        integer, intent(inout) :: ce(:) !! Cluster end positions permuted with cs.
        integer, intent(in) :: ncl !! Number of active clusters.
        integer :: i
        integer :: j
        integer :: ks
        integer :: ke

        do i = 2, ncl
            ks = cs(i)
            ke = ce(i)
            j = i - 1
            do while (j >= 1)
                if (cs(j) < ks) exit
                if (cs(j) == ks .and. ce(j) <= ke) exit
                cs(j + 1) = cs(j)
                ce(j + 1) = ce(j)
                j = j - 1
            end do
            cs(j + 1) = ks
            ce(j + 1) = ke
        end do
    end subroutine sort_cluster_bounds

    subroutine sort_cluster_sizes(cs, ce, order, ncl, minimum)
        integer, intent(in) :: cs(:) !! Cluster start positions.
        integer, intent(in) :: ce(:) !! Cluster end positions.
        integer, intent(inout) :: order(:) !! Cluster indices sorted by interval size.
        integer, intent(in) :: ncl !! Number of active clusters.
        logical, intent(in) :: minimum !! Sort ascending size when true; descending otherwise.
        integer :: i
        integer :: j
        integer :: key
        integer :: s1
        integer :: s2
        logical :: move

        do i = 2, ncl
            key = order(i)
            j = i - 1
            do while (j >= 1)
                s1 = ce(order(j)) - cs(order(j))
                s2 = ce(key) - cs(key)
                if (minimum) then
                    move = s1 > s2
                else
                    move = s1 < s2
                end if
                if (.not. move) exit
                order(j + 1) = order(j)
                j = j - 1
            end do
            order(j + 1) = key
        end do
    end subroutine sort_cluster_sizes

    subroutine empty_result(result)
        type(optics_result), intent(out) :: result !! OPTICS result initialized with zero-length arrays.
        allocate(result%order(0), result%predecessor(0), result%reachdist(0), result%coredist(0))
    end subroutine empty_result

end module dbscan_optics
