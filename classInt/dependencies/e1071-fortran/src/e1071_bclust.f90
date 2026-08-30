module e1071_bclust
    use e1071_kinds, only: dp
    use e1071_rng, only: rng_state, rng_integer
    use e1071_fuzzy, only: fuzzy_cluster_model, cmeans_fit
    use proxy, only: proxy_dist_auto
    implicit none
    private

    type, public :: bclust_model
        real(dp), allocatable :: all_centers(:, :)
        integer, allocatable :: all_cluster(:)
        integer, allocatable :: members(:, :)
        integer, allocatable :: cluster(:)
        real(dp), allocatable :: centers(:, :)
        real(dp), allocatable :: data_mean(:)
        integer :: iter_base = 0
        integer :: base_centers = 0
        integer :: max_cluster = 0
        character(len=32) :: distance_method = "euclidean"
        character(len=16) :: hclust_method = "average"
    end type bclust_model

    public :: bclust_fit, bclust_centers, bclust_clusters

contains

    subroutine bclust_fit(x, centers, rng, model, iter_base, minsize, distance_method, hclust_method, base_method, &
                          base_centers, final_kmeans, resample, maxcluster)
        real(dp), intent(in) :: x(:, :) !! Observation-by-variable data matrix used for base clusterings and final assignments.
        integer, intent(in) :: centers !! Requested final committee cluster count; must be between two and maxcluster.
        type(rng_state), intent(inout) :: rng !! Mutable RNG used for bootstrap samples and base-center initialization.
        type(bclust_model), intent(out) :: model !! Fitted bagged-clustering committee, hierarchical memberships, centers, and
        !! labels.
        integer, intent(in), optional :: iter_base !! Number of bootstrap committee members; defaults to 10 and must be positive.
        integer, intent(in), optional :: minsize !! Minimum observations assigned to a retained base center; zero disables pruning.
        character(len=*), intent(in), optional :: distance_method !! proxy distance used between base centers; defaults to
        !! Euclidean.
        character(len=*), intent(in), optional :: hclust_method !! Agglomerative linkage: average, single, or complete; default
        !! average.
        character(len=*), intent(in), optional :: base_method !! Base clustering algorithm: "kmeans" or "cmeans"; default kmeans.
        integer, intent(in), optional :: base_centers !! Number of centers in each bootstrap base fit; defaults to 20.
        logical, intent(in), optional :: final_kmeans !! If true, refine final committee centers by k-means on x; default false.
        logical, intent(in), optional :: resample !! If true, bootstrap rows for each base fit; default true.
        integer, intent(in), optional :: maxcluster !! Largest hierarchical cut retained; defaults to base_centers.
        real(dp), allocatable :: sample(:, :)
        real(dp), allocatable :: initial(:, :)
        real(dp), allocatable :: base_fit(:, :)
        real(dp), allocatable :: tmp_centers(:, :)
        integer, allocatable :: labels(:)
        integer :: use_iter
        integer :: use_base
        integer :: use_min
        integer :: use_max
        integer :: b
        integer :: i
        integer :: row
        integer :: offset
        logical :: do_resample
        logical :: do_final
        character(len=:), allocatable :: use_base_method
        character(len=:), allocatable :: use_dist
        character(len=:), allocatable :: use_link

        use_iter = 10
        if (present(iter_base)) use_iter = iter_base
        use_base = 20
        if (present(base_centers)) use_base = base_centers
        use_min = 0
        if (present(minsize)) use_min = minsize
        use_max = use_base
        if (present(maxcluster)) use_max = maxcluster
        do_resample = .true.
        if (present(resample)) do_resample = resample
        do_final = .false.
        if (present(final_kmeans)) do_final = final_kmeans
        use_base_method = "kmeans"
        if (present(base_method)) use_base_method = trim(adjustl(base_method))
        use_dist = "euclidean"
        if (present(distance_method)) use_dist = trim(adjustl(distance_method))
        use_link = "average"
        if (present(hclust_method)) use_link = trim(adjustl(hclust_method))
        if (use_iter < 1 .or. use_base < 2 .or. use_base > size(x, 1)) error stop "bclust_fit: invalid base controls"
        if (use_max < 2 .or. use_max > use_iter * use_base) error stop "bclust_fit: invalid maxcluster"
        if (centers < 2 .or. centers > use_max) error stop "bclust_fit: invalid final centers"

        allocate(model%all_centers(use_iter * use_base, size(x, 2)))
        allocate(sample(size(x, 1), size(x, 2)), initial(use_base, size(x, 2)))
        offset = 0
        do b = 1, use_iter
            if (do_resample) then
                do i = 1, size(x, 1)
                    row = rng_integer(rng, size(x, 1))
                    sample(i, :) = x(row, :)
                end do
            else
                sample = x
            end if
            call random_initial_centers(sample, use_base, rng, initial)
            if (use_base_method == "cmeans") then
                call base_cmeans(sample, initial, base_fit)
            else if (use_base_method == "kmeans") then
                call kmeans_fit(sample, initial, base_fit, labels)
            else
                error stop "bclust_fit: base_method must be kmeans or cmeans"
            end if
            model%all_centers(offset + 1:offset + use_base, :) = base_fit
            offset = offset + use_base
        end do
        call nearest_center_indices(model%all_centers, x, model%all_cluster)
        if (use_min > 0) call prune_centers(model%all_centers, x, use_min, model%all_cluster)
        use_max = min(use_max, size(model%all_centers, 1))
        if (centers > use_max) error stop "bclust_fit: pruning left fewer centers than requested"
        call hierarchical_memberships(model%all_centers, use_dist, use_link, use_max, model%members)
        call bclust_centers_from_members(model%all_centers, model%members(:, centers - 1), centers, tmp_centers)
        call nearest_center_indices(model%all_centers, x, model%all_cluster)
        allocate(model%cluster(size(x, 1)))
        do i = 1, size(x, 1)
            model%cluster(i) = model%members(model%all_cluster(i), centers - 1)
        end do
        if (do_final) then
            call kmeans_fit(x, tmp_centers, model%centers, model%cluster)
        else
            model%centers = tmp_centers
        end if
        allocate(model%data_mean(size(x, 2)))
        model%data_mean = sum(x, dim=1) / real(size(x, 1), dp)
        model%iter_base = use_iter
        model%base_centers = use_base
        model%max_cluster = use_max
        model%distance_method = use_dist
        model%hclust_method = use_link
    end subroutine bclust_fit

    function bclust_centers(model, k) result(centers)
        type(bclust_model), intent(in) :: model !! Fitted bagged-clustering model with hierarchical base-center memberships.
        integer, intent(in) :: k !! Requested hierarchical cut count between two and model%max_cluster.
        real(dp), allocatable :: centers(:, :)

        if (k < 2 .or. k > model%max_cluster) error stop "bclust_centers: invalid k"
        call bclust_centers_from_members(model%all_centers, model%members(:, k - 1), k, centers)
    end function bclust_centers

    function bclust_clusters(model, x, k) result(cluster)
        type(bclust_model), intent(in) :: model !! Fitted bagged-clustering model.
        real(dp), intent(in) :: x(:, :) !! Observations to classify through their nearest retained base center.
        integer, intent(in) :: k !! Requested hierarchical cut count between two and model%max_cluster.
        integer, allocatable :: cluster(:)
        integer, allocatable :: nearest(:)
        integer :: i

        if (k < 2 .or. k > model%max_cluster) error stop "bclust_clusters: invalid k"
        call nearest_center_indices(model%all_centers, x, nearest)
        allocate(cluster(size(x, 1)))
        do i = 1, size(x, 1)
            cluster(i) = model%members(nearest(i), k - 1)
        end do
    end function bclust_clusters

    subroutine base_cmeans(x, initial, centers)
        real(dp), intent(in) :: x(:, :) !! Bootstrap sample passed to fuzzy c-means.
        real(dp), intent(in) :: initial(:, :) !! Initial fuzzy center matrix.
        real(dp), allocatable, intent(out) :: centers(:, :) !! Final fuzzy c-means center matrix.
        type(fuzzy_cluster_model) :: fit

        call cmeans_fit(x, initial, fit)
        centers = fit%centers
    end subroutine base_cmeans

    subroutine kmeans_fit(x, initial, centers, labels)
        real(dp), intent(in) :: x(:, :) !! Observation-by-variable matrix used for Lloyd k-means.
        real(dp), intent(in) :: initial(:, :) !! Initial center matrix; its row count fixes k.
        real(dp), allocatable, intent(out) :: centers(:, :) !! Final arithmetic-mean cluster centers.
        integer, allocatable, intent(out) :: labels(:) !! One-based nearest-center assignment for every observation.
        real(dp), allocatable :: previous(:, :)
        integer, allocatable :: count_cluster(:)
        integer :: iter
        integer :: i
        integer :: c

        centers = initial
        allocate(previous(size(initial, 1), size(initial, 2)), labels(size(x, 1)), count_cluster(size(initial, 1)))
        do iter = 1, 100
            previous = centers
            do i = 1, size(x, 1)
                labels(i) = nearest_row(centers, x(i, :))
            end do
            centers = 0.0_dp
            count_cluster = 0
            do i = 1, size(x, 1)
                c = labels(i)
                centers(c, :) = centers(c, :) + x(i, :)
                count_cluster(c) = count_cluster(c) + 1
            end do
            do c = 1, size(centers, 1)
                if (count_cluster(c) > 0) then
                    centers(c, :) = centers(c, :) / real(count_cluster(c), dp)
                else
                    centers(c, :) = previous(c, :)
                end if
            end do
            if (maxval(abs(centers - previous)) < 1.0e-10_dp) exit
        end do
        do i = 1, size(x, 1)
            labels(i) = nearest_row(centers, x(i, :))
        end do
    end subroutine kmeans_fit

    subroutine random_initial_centers(x, k, rng, centers)
        real(dp), intent(in) :: x(:, :) !! Candidate observation rows from which initial centers are sampled without replacement.
        integer, intent(in) :: k !! Number of distinct row indices selected.
        type(rng_state), intent(inout) :: rng !! Mutable random generator used for sampling row indices.
        real(dp), intent(out) :: centers(:, :) !! Preallocated k-by-nvar matrix receiving sampled rows.
        integer, allocatable :: available(:)
        integer :: i
        integer :: pick
        integer :: row

        allocate(available(size(x, 1)))
        available = [(i, i = 1, size(x, 1))]
        do i = 1, k
            pick = rng_integer(rng, size(x, 1) - i + 1)
            row = available(pick)
            centers(i, :) = x(row, :)
            available(pick) = available(size(x, 1) - i + 1)
        end do
    end subroutine random_initial_centers

    subroutine nearest_center_indices(centers, x, index)
        real(dp), intent(in) :: centers(:, :) !! Candidate center matrix used in Euclidean 1-nearest-neighbor assignment.
        real(dp), intent(in) :: x(:, :) !! Observation-by-variable matrix assigned to centers.
        integer, allocatable, intent(out) :: index(:) !! One-based nearest center for every observation.
        integer :: i

        allocate(index(size(x, 1)))
        do i = 1, size(x, 1)
            index(i) = nearest_row(centers, x(i, :))
        end do
    end subroutine nearest_center_indices

    pure function nearest_row(centers, x) result(index)
        real(dp), intent(in) :: centers(:, :) !! Center rows compared by squared Euclidean distance.
        real(dp), intent(in) :: x(:) !! Observation vector with the same variable count as centers.
        integer :: index
        real(dp) :: best
        real(dp) :: d
        integer :: c

        index = 1
        best = sum((centers(1, :) - x)**2)
        do c = 2, size(centers, 1)
            d = sum((centers(c, :) - x)**2)
            if (d < best) then
                best = d
                index = c
            end if
        end do
    end function nearest_row

    subroutine prune_centers(centers, x, minsize, all_cluster)
        real(dp), allocatable, intent(inout) :: centers(:, :) !! Retained base centers compacted repeatedly until every center
        !! meets minsize.
        real(dp), intent(in) :: x(:, :) !! Original observations used to count nearest-center assignments.
        integer, intent(in) :: minsize !! Minimum number of observations assigned to every retained center; must be positive here.
        integer, allocatable, intent(inout) :: all_cluster(:) !! Nearest-center labels updated after each pruning pass.
        integer, allocatable :: counts(:)
        logical, allocatable :: keep(:)
        real(dp), allocatable :: reduced(:, :)
        integer :: i
        integer :: n

        do
            call nearest_center_indices(centers, x, all_cluster)
            allocate(counts(size(centers, 1)), keep(size(centers, 1)))
            counts = 0
            do i = 1, size(all_cluster)
                counts(all_cluster(i)) = counts(all_cluster(i)) + 1
            end do
            keep = counts >= minsize
            if (all(keep)) then
                deallocate(counts, keep)
                exit
            end if
            n = count(keep)
            if (n == 0) error stop "prune_centers: minsize removed every base center"
            allocate(reduced(n, size(centers, 2)))
            n = 0
            do i = 1, size(centers, 1)
                if (.not. keep(i)) cycle
                n = n + 1
                reduced(n, :) = centers(i, :)
            end do
            call move_alloc(reduced, centers)
            deallocate(counts, keep)
        end do
    end subroutine prune_centers


    subroutine hierarchical_memberships(centers, distance_method, linkage_method, maxcluster, members)
        real(dp), intent(in) :: centers(:, :) !! Base-center rows that are agglomerated into the bagged-clustering hierarchy.
        character(len=*), intent(in) :: distance_method !! proxy distance name used to measure dissimilarity between base centers.
        character(len=*), intent(in) :: linkage_method !! Agglomerative linkage rule: "average", "single", or "complete".
        integer, intent(in) :: maxcluster !! Largest cluster count retained in `members`; must be at least two.
        integer, allocatable, intent(out) :: members(:, :) !! Base-center labels; column k-1 contains the k-cluster cut.
        real(dp), allocatable :: distance(:, :)
        integer, allocatable :: label(:)
        integer :: n
        integer :: current
        integer :: a
        integer :: b
        integer :: best_a
        integer :: best_b
        integer :: status
        real(dp) :: d
        real(dp) :: best

        n = size(centers, 1)
        if (n < 2) error stop "hierarchical_memberships: at least two centers are required"
        if (maxcluster < 2 .or. maxcluster > n) error stop "hierarchical_memberships: invalid maxcluster"
        call proxy_dist_auto(centers, distance_method, distance, status=status)
        if (status /= 0) error stop "hierarchical_memberships: unknown or invalid proxy distance"
        allocate(label(n), members(n, maxcluster - 1))
        label = [(a, a = 1, n)]
        members = 0
        current = n
        if (current <= maxcluster) members(:, current - 1) = label
        do while (current > 2)
            best = huge(1.0_dp)
            best_a = 1
            best_b = 2
            do a = 1, current - 1
                do b = a + 1, current
                    d = cluster_linkage(distance, label, a, b, linkage_method)
                    if (d < best) then
                        best = d
                        best_a = a
                        best_b = b
                    end if
                end do
            end do
            call merge_labels(label, best_a, best_b)
            current = current - 1
            if (current <= maxcluster) members(:, current - 1) = label
        end do
    end subroutine hierarchical_memberships

    function cluster_linkage(distance, label, a, b, method) result(value)
        real(dp), intent(in) :: distance(:, :) !! Symmetric base-center distance matrix used by the linkage calculation.
        integer, intent(in) :: label(:) !! Current one-based cluster label for every base center.
        integer, intent(in) :: a !! First current cluster label whose linkage to `b` is requested.
        integer, intent(in) :: b !! Second current cluster label; it must differ from `a`.
        character(len=*), intent(in) :: method !! Linkage rule: "average", "single", or "complete".
        real(dp) :: value
        integer :: i
        integer :: j
        integer :: count_pair
        character(len=:), allocatable :: key

        key = trim(adjustl(method))
        if (key == "single") then
            value = huge(1.0_dp)
        else if (key == "complete") then
            value = -huge(1.0_dp)
        else if (key == "average") then
            value = 0.0_dp
        else
            error stop "cluster_linkage: unsupported linkage method"
        end if
        count_pair = 0
        do i = 1, size(label)
            if (label(i) /= a) cycle
            do j = 1, size(label)
                if (label(j) /= b) cycle
                count_pair = count_pair + 1
                if (key == "single") then
                    value = min(value, distance(i, j))
                else if (key == "complete") then
                    value = max(value, distance(i, j))
                else
                    value = value + distance(i, j)
                end if
            end do
        end do
        if (count_pair < 1) error stop "cluster_linkage: empty cluster encountered"
        if (key == "average") value = value / real(count_pair, dp)
    end function cluster_linkage

    subroutine merge_labels(label, keep_label, remove_label)
        integer, intent(inout) :: label(:) !! Current cluster labels, compacted in place after the requested merge.
        integer, intent(in) :: keep_label !! Lower cluster label retained for the merged cluster.
        integer, intent(in) :: remove_label !! Higher cluster label removed and shifted out of the compact label sequence.
        integer :: i

        if (keep_label >= remove_label) error stop "merge_labels: labels must satisfy keep_label < remove_label"
        do i = 1, size(label)
            if (label(i) == remove_label) label(i) = keep_label
            if (label(i) > remove_label) label(i) = label(i) - 1
        end do
    end subroutine merge_labels

    subroutine bclust_centers_from_members(all_centers, member, k, centers)
        real(dp), intent(in) :: all_centers(:, :) !! Base-center matrix whose rows are averaged within hierarchical groups.
        integer, intent(in) :: member(:) !! One-based hierarchical group label for every row of `all_centers`.
        integer, intent(in) :: k !! Expected number of groups represented by `member`.
        real(dp), allocatable, intent(out) :: centers(:, :) !! Arithmetic mean center of each hierarchical group.
        integer, allocatable :: count_cluster(:)
        integer :: i
        integer :: c

        if (size(member) /= size(all_centers, 1)) error stop "bclust_centers_from_members: shape mismatch"
        allocate(centers(k, size(all_centers, 2)), count_cluster(k))
        centers = 0.0_dp
        count_cluster = 0
        do i = 1, size(member)
            c = member(i)
            if (c < 1 .or. c > k) error stop "bclust_centers_from_members: invalid group label"
            centers(c, :) = centers(c, :) + all_centers(i, :)
            count_cluster(c) = count_cluster(c) + 1
        end do
        do c = 1, k
            if (count_cluster(c) < 1) error stop "bclust_centers_from_members: empty group"
            centers(c, :) = centers(c, :) / real(count_cluster(c), dp)
        end do
    end subroutine bclust_centers_from_members

end module e1071_bclust
