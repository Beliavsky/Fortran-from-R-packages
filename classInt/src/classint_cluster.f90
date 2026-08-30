! SPDX-License-Identifier: GPL-2.0-or-later
! SPDX-FileComment: Modern Fortran translation of R package classInt computational code.
module classint_cluster
    use classint_kinds, only: dp
    use classint_types, only: hclust_1d_model, classint_options
    use classint_utils, only: sort_real, unique_sorted, lower_text
    use e1071, only: rng_state, rng_integer, bclust_model, bclust_fit, bclust_centers, bclust_clusters
    implicit none
    private

    public :: kmeans_breaks, hclust_fit_1d, hclust_cut_breaks, bclust_fit_breaks, bclust_cut_breaks
    public :: breaks_from_labels

contains

    subroutine kmeans_breaks(x, k, rng, iter_max, nstart, breaks)
        real(dp), intent(in) :: x(:) !! Finite one-dimensional observations used by Lloyd k-means classification.
        integer, intent(in) :: k !! Number of requested clusters; must not exceed the number of distinct observations.
        type(rng_state), intent(inout) :: rng !! Mutable RNG state used to choose initial distinct centers for each start.
        integer, intent(in) :: iter_max !! Maximum Lloyd iterations per random start; must be positive.
        integer, intent(in) :: nstart !! Number of independent initial-center sets; must be positive.
        real(dp), allocatable, intent(out) :: breaks(:) !! K+1 ascending breaks formed midway between adjacent cluster ranges.
        real(dp), allocatable :: unique_x(:)
        real(dp), allocatable :: centers(:)
        real(dp), allocatable :: new_centers(:)
        real(dp), allocatable :: best_centers(:)
        integer, allocatable :: labels(:)
        integer, allocatable :: best_labels(:)
        integer, allocatable :: chosen(:)
        integer, allocatable :: counts(:)
        real(dp) :: objective
        real(dp) :: best_objective
        real(dp) :: distance
        integer :: start
        integer :: iter
        integer :: i
        integer :: j
        integer :: c
        integer :: pick
        logical :: changed

        if (size(x) < 1) error stop "kmeans_breaks: empty input"
        if (iter_max < 1 .or. nstart < 1) error stop "kmeans_breaks: invalid iteration controls"
        unique_x = unique_sorted(x)
        if (k < 1 .or. k > size(unique_x)) error stop "kmeans_breaks: invalid cluster count"
        allocate(centers(k), new_centers(k), labels(size(x)), chosen(size(unique_x)), counts(k))
        best_objective = huge(1.0_dp)
        do start = 1, nstart
            chosen = 0
            do c = 1, k
                do
                    pick = rng_integer(rng, size(unique_x))
                    if (chosen(pick) == 0) exit
                end do
                chosen(pick) = 1
                centers(c) = unique_x(pick)
            end do
            labels = 0
            do iter = 1, iter_max
                changed = .false.
                do i = 1, size(x)
                    c = 1
                    distance = abs(x(i) - centers(1))
                    do j = 2, k
                        if (abs(x(i) - centers(j)) < distance) then
                            c = j
                            distance = abs(x(i) - centers(j))
                        end if
                    end do
                    if (labels(i) /= c) changed = .true.
                    labels(i) = c
                end do
                new_centers = 0.0_dp
                counts = 0
                do i = 1, size(x)
                    c = labels(i)
                    new_centers(c) = new_centers(c) + x(i)
                    counts(c) = counts(c) + 1
                end do
                do c = 1, k
                    if (counts(c) == 0) then
                        call reinitialize_empty_center(x, centers, labels, new_centers(c))
                    else
                        new_centers(c) = new_centers(c) / real(counts(c), dp)
                    end if
                end do
                if (maxval(abs(new_centers - centers)) <= 1.0e-12_dp * max(1.0_dp, maxval(abs(centers)))) then
                    centers = new_centers
                    exit
                end if
                centers = new_centers
                if (.not. changed) exit
            end do
            call assign_nearest(x, centers, labels)
            objective = 0.0_dp
            do i = 1, size(x)
                objective = objective + (x(i) - centers(labels(i)))**2
            end do
            if (objective < best_objective) then
                best_objective = objective
                best_centers = centers
                best_labels = labels
            end if
        end do
        call order_labels_by_center(best_centers, best_labels)
        call breaks_from_labels(x, best_labels, k, breaks)
    end subroutine kmeans_breaks

    subroutine reinitialize_empty_center(x, centers, labels, center)
        real(dp), intent(in) :: x(:) !! Observations searched for the point farthest from its currently assigned center.
        real(dp), intent(in) :: centers(:) !! Current center vector used to measure each observation's residual distance.
        integer, intent(in) :: labels(:) !! Current one-based cluster assignment for every observation.
        real(dp), intent(out) :: center !! Replacement center set to the most poorly represented observation.
        real(dp) :: best
        real(dp) :: d
        integer :: i

        best = -1.0_dp
        center = x(1)
        do i = 1, size(x)
            d = abs(x(i) - centers(labels(i)))
            if (d > best) then
                best = d
                center = x(i)
            end if
        end do
    end subroutine reinitialize_empty_center

    subroutine assign_nearest(x, centers, labels)
        real(dp), intent(in) :: x(:) !! Observations assigned to the nearest center by absolute one-dimensional distance.
        real(dp), intent(in) :: centers(:) !! Candidate cluster centers.
        integer, intent(out) :: labels(:) !! One-based nearest-center label for each observation; same length as x.
        real(dp) :: best
        real(dp) :: d
        integer :: i
        integer :: c

        if (size(labels) /= size(x)) error stop "assign_nearest: shape mismatch"
        do i = 1, size(x)
            labels(i) = 1
            best = abs(x(i) - centers(1))
            do c = 2, size(centers)
                d = abs(x(i) - centers(c))
                if (d < best) then
                    best = d
                    labels(i) = c
                end if
            end do
        end do
    end subroutine assign_nearest

    subroutine order_labels_by_center(centers, labels)
        real(dp), intent(in) :: centers(:) !! Unordered center values whose numeric order defines class order.
        integer, intent(inout) :: labels(:) !! Cluster labels remapped in place so class one has the smallest center.
        integer, allocatable :: order(:)
        integer, allocatable :: map(:)
        integer :: i
        integer :: j
        integer :: tmp

        allocate(order(size(centers)), map(size(centers)))
        do i = 1, size(centers)
            order(i) = i
        end do
        do i = 2, size(order)
            tmp = order(i)
            j = i - 1
            do while (j >= 1)
                if (centers(order(j)) <= centers(tmp)) exit
                order(j + 1) = order(j)
                j = j - 1
            end do
            order(j + 1) = tmp
        end do
        do i = 1, size(order)
            map(order(i)) = i
        end do
        do i = 1, size(labels)
            labels(i) = map(labels(i))
        end do
    end subroutine order_labels_by_center

    subroutine breaks_from_labels(x, labels, k, breaks)
        real(dp), intent(in) :: x(:) !! Finite observations partitioned into ordered classes.
        integer, intent(in) :: labels(:) !! One-based ordered class label for each observation; same length as x.
        integer, intent(in) :: k !! Number of classes represented in labels.
        real(dp), allocatable, intent(out) :: breaks(:) !! K+1 break vector using class extremes and adjacent midpoints.
        real(dp), allocatable :: lo(:)
        real(dp), allocatable :: hi(:)
        integer :: i
        integer :: c

        if (size(labels) /= size(x)) error stop "breaks_from_labels: shape mismatch"
        allocate(lo(k), source=huge(1.0_dp))
        allocate(hi(k), source=-huge(1.0_dp))
        do i = 1, size(x)
            c = labels(i)
            if (c < 1 .or. c > k) error stop "breaks_from_labels: invalid label"
            lo(c) = min(lo(c), x(i))
            hi(c) = max(hi(c), x(i))
        end do
        if (any(lo >= huge(1.0_dp) / 2.0_dp)) error stop "breaks_from_labels: empty class"
        allocate(breaks(k + 1))
        breaks(1) = lo(1)
        do c = 1, k - 1
            breaks(c + 1) = 0.5_dp * (hi(c) + lo(c + 1))
        end do
        breaks(k + 1) = hi(k)
    end subroutine breaks_from_labels

    subroutine hclust_fit_1d(x, method, model)
        real(dp), intent(in) :: x(:) !! Finite one-dimensional observations; duplicates are retained as distinct leaves.
        character(len=*), intent(in) :: method !! Linkage name: complete, single, average, centroid, ward.D, or ward.D2.
        type(hclust_1d_model), intent(out) :: model !! Agglomerative tree represented by parent links for all leaf/internal nodes.
        real(dp), allocatable :: srt(:)
        real(dp), allocatable :: cmean(:)
        real(dp), allocatable :: cmin(:)
        real(dp), allocatable :: cmax(:)
        integer, allocatable :: count_cluster(:)
        integer, allocatable :: node(:)
        character(len=:), allocatable :: use_method
        real(dp) :: d
        real(dp) :: best
        integer :: n
        integer :: m
        integer :: step
        integer :: j
        integer :: pair
        integer :: newnode
        integer :: total

        n = size(x)
        if (n < 2) error stop "hclust_fit_1d: at least two observations are required"
        use_method = trim(adjustl(lower_text(method)))
        select case (use_method)
        case ("complete", "single", "average", "centroid", "ward.d", "ward.d2")
        case default
            error stop "hclust_fit_1d: unsupported linkage"
        end select
        srt = x
        call sort_real(srt)
        allocate(cmean(n), cmin(n), cmax(n), count_cluster(n), node(n))
        cmean = srt
        cmin = srt
        cmax = srt
        count_cluster = 1
        do j = 1, n
            node(j) = j
        end do
        allocate(model%parent(2 * n - 1), source=0)
        allocate(model%node_mean(2 * n - 1), source=0.0_dp)
        model%node_mean(1:n) = srt
        m = n
        do step = 1, n - 1
            pair = 1
            best = linkage_distance(use_method, cmean(1), cmean(2), cmin(1), cmax(1), cmin(2), cmax(2), &
                                    count_cluster(1), count_cluster(2))
            do j = 2, m - 1
                d = linkage_distance(use_method, cmean(j), cmean(j + 1), cmin(j), cmax(j), cmin(j + 1), cmax(j + 1), &
                                     count_cluster(j), count_cluster(j + 1))
                if (d < best) then
                    best = d
                    pair = j
                end if
            end do
            newnode = n + step
            model%parent(node(pair)) = newnode
            model%parent(node(pair + 1)) = newnode
            total = count_cluster(pair) + count_cluster(pair + 1)
            cmean(pair) = (real(count_cluster(pair), dp) * cmean(pair) + &
                           real(count_cluster(pair + 1), dp) * cmean(pair + 1)) / real(total, dp)
            cmin(pair) = min(cmin(pair), cmin(pair + 1))
            cmax(pair) = max(cmax(pair), cmax(pair + 1))
            count_cluster(pair) = total
            node(pair) = newnode
            model%node_mean(newnode) = cmean(pair)
            do j = pair + 1, m - 1
                cmean(j) = cmean(j + 1)
                cmin(j) = cmin(j + 1)
                cmax(j) = cmax(j + 1)
                count_cluster(j) = count_cluster(j + 1)
                node(j) = node(j + 1)
            end do
            m = m - 1
        end do
        model%n = n
        model%method = use_method
    end subroutine hclust_fit_1d

    pure elemental function linkage_distance(method, mean1, mean2, min1, max1, min2, max2, n1, n2) result(d)
        !! Computes the distance between two adjacent clusters for the requested linkage rule.
        character(len=*), intent(in) :: method !! Normalized linkage name controlling the adjacent-cluster merge criterion.
        real(dp), intent(in) :: mean1 !! Arithmetic mean of the lower-valued active cluster.
        real(dp), intent(in) :: mean2 !! Arithmetic mean of the higher-valued active cluster.
        real(dp), intent(in) :: min1 !! Minimum observation in the lower-valued cluster.
        real(dp), intent(in) :: max1 !! Maximum observation in the lower-valued cluster.
        real(dp), intent(in) :: min2 !! Minimum observation in the higher-valued cluster.
        real(dp), intent(in) :: max2 !! Maximum observation in the higher-valued cluster.
        integer, intent(in) :: n1 !! Number of observations in the lower-valued cluster.
        integer, intent(in) :: n2 !! Number of observations in the higher-valued cluster.
        real(dp) :: d

        select case (method)
        case ("single")
            d = max(0.0_dp, min2 - max1)
        case ("complete")
            d = max2 - min1
        case ("average", "centroid")
            d = abs(mean2 - mean1)
        case ("ward.d", "ward.d2")
            d = real(n1 * n2, dp) / real(n1 + n2, dp) * (mean2 - mean1)**2
        case default
            error stop "linkage_distance: unsupported linkage"
        end select
    end function linkage_distance

    subroutine hclust_cut_breaks(x, model, k, breaks)
        real(dp), intent(in) :: x(:) !! Original finite observations used to reconstruct ordered class ranges for a tree cut.
        type(hclust_1d_model), intent(in) :: model !! Hierarchical model previously fitted to the same observation multiset.
        integer, intent(in) :: k !! Requested number of clusters between one and model%n.
        real(dp), allocatable, intent(out) :: breaks(:) !! K+1 ascending break vector for the requested hierarchical cut.
        real(dp), allocatable :: srt(:)
        integer, allocatable :: labels(:)

        if (size(x) /= model%n) error stop "hclust_cut_breaks: data size differs from fitted tree"
        srt = x
        call sort_real(srt)
        call hclust_cut_labels(model, k, labels)
        call breaks_from_labels(srt, labels, k, breaks)
    end subroutine hclust_cut_breaks

    subroutine hclust_cut_labels(model, k, labels)
        type(hclust_1d_model), intent(in) :: model !! Hierarchical parent tree whose first n nodes are sorted observations.
        integer, intent(in) :: k !! Number of clusters requested from the tree cut.
        integer, allocatable, intent(out) :: labels(:) !! Ordered cluster label for each sorted leaf observation.
        integer, allocatable :: roots(:)
        integer, allocatable :: unique_roots(:)
        integer, allocatable :: order(:)
        integer :: threshold
        integer :: i
        integer :: r
        integer :: nroot
        integer :: j
        integer :: tmp

        if (k < 1 .or. k > model%n) error stop "hclust_cut_labels: invalid k"
        threshold = model%n + (model%n - k)
        allocate(roots(model%n))
        do i = 1, model%n
            r = i
            do while (model%parent(r) > 0 .and. model%parent(r) <= threshold)
                r = model%parent(r)
            end do
            roots(i) = r
        end do
        allocate(unique_roots(k), source=0)
        nroot = 0
        do i = 1, model%n
            if (.not. any(unique_roots(1:nroot) == roots(i))) then
                nroot = nroot + 1
                if (nroot > k) error stop "hclust_cut_labels: too many roots"
                unique_roots(nroot) = roots(i)
            end if
        end do
        if (nroot /= k) error stop "hclust_cut_labels: root count mismatch"
        allocate(order(k))
        do i = 1, k
            order(i) = i
        end do
        do i = 2, k
            tmp = order(i)
            j = i - 1
            do while (j >= 1)
                if (model%node_mean(unique_roots(order(j))) <= model%node_mean(unique_roots(tmp))) exit
                order(j + 1) = order(j)
                j = j - 1
            end do
            order(j + 1) = tmp
        end do
        allocate(labels(model%n))
        do i = 1, model%n
            do j = 1, k
                if (roots(i) == unique_roots(order(j))) then
                    labels(i) = j
                    exit
                end if
            end do
        end do
    end subroutine hclust_cut_labels

    subroutine bclust_fit_breaks(x, k, rng, options, model, breaks)
        real(dp), intent(in) :: x(:) !! Finite observations supplied as a one-column matrix to translated e1071 bagged clustering.
        integer, intent(in) :: k !! Requested final bagged-clustering class count.
        type(rng_state), intent(inout) :: rng !! Mutable RNG shared with e1071 bootstrap and base-clustering initialization.
        type(classint_options), intent(in) :: options !! Bagged-clustering controls mirroring the classInt/e1071 style options.
        type(bclust_model), intent(out) :: model !! Fitted translated e1071 bagged-clustering model retained for later cuts.
        real(dp), allocatable, intent(out) :: breaks(:) !! K+1 breaks derived from the fitted bagged-cluster assignment.
        real(dp), allocatable :: matrix_x(:, :)
        integer :: base_centers
        integer :: maxcluster

        allocate(matrix_x(size(x), 1))
        matrix_x(:, 1) = x
        base_centers = min(options%bclust_base_centers, size(x))
        maxcluster = min(options%bclust_maxcluster, options%bclust_iter_base * base_centers)
        if (base_centers < 2 .or. k > maxcluster) error stop "bclust_fit_breaks: incompatible class controls"
        call bclust_fit(matrix_x, k, rng, model, iter_base=options%bclust_iter_base, &
                        minsize=options%bclust_minsize, distance_method=trim(options%bclust_distance), &
                        hclust_method=trim(options%bclust_hclust), base_method=trim(options%bclust_base_method), &
                        base_centers=base_centers, final_kmeans=options%bclust_final_kmeans, &
                        resample=options%bclust_resample, maxcluster=maxcluster)
        call bclust_cut_breaks(x, model, k, breaks)
    end subroutine bclust_fit_breaks

    subroutine bclust_cut_breaks(x, model, k, breaks)
        real(dp), intent(in) :: x(:) !! Finite observations classified through the retained bagged-clustering committee.
        type(bclust_model), intent(in) :: model !! Fitted translated e1071 bagged-clustering model.
        integer, intent(in) :: k !! Requested hierarchical cut count supported by the fitted committee.
        real(dp), allocatable, intent(out) :: breaks(:) !! K+1 ordered breaks derived from bagged-cluster assignments.
        real(dp), allocatable :: matrix_x(:, :)
        real(dp), allocatable :: centers(:, :)
        real(dp), allocatable :: center1(:)
        integer, allocatable :: labels(:)

        allocate(matrix_x(size(x), 1))
        matrix_x(:, 1) = x
        labels = bclust_clusters(model, matrix_x, k)
        centers = bclust_centers(model, k)
        center1 = centers(:, 1)
        call order_labels_by_center(center1, labels)
        call breaks_from_labels(x, labels, k, breaks)
    end subroutine bclust_cut_breaks
end module classint_cluster
