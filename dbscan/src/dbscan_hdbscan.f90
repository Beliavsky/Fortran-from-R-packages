module dbscan_hdbscan
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan
    use dbscan_kinds, only : dp
    use dbscan_types, only : hdbscan_result, dendrogram_result
    use dbscan_nn, only : knn_dist
    use dbscan_distance, only : pairwise_distance_matrix
    implicit none
    private

    public :: coredist, mrdist, mst_dense, hdbscan, glosh, extract_fosc

    type :: condensed_tree_data
        integer :: n = 0
        integer :: ncluster = 0
        integer, allocatable :: parent(:)
        integer, allocatable :: first_child(:)
        integer, allocatable :: next_sibling(:)
        integer, allocatable :: point_owner(:)
        integer, allocatable :: n_children(:)
        real(dp), allocatable :: birth_eps(:)
        real(dp), allocatable :: death_eps(:)
        real(dp), allocatable :: stability(:)
        real(dp), allocatable :: point_eps(:)
    end type condensed_tree_data

contains

    function coredist(x, min_pts, status) result(cd)
        real(dp), intent(in) :: x(:, :) !! Data matrix with observations in rows.
        integer, intent(in) :: min_pts !! HDBSCAN MinPts including self; at least two.
        integer, intent(out), optional :: status !! Zero on success; nonzero for invalid input.
        real(dp), allocatable :: cd(:)
        integer :: stat

        if (min_pts < 2 .or. min_pts > size(x, 1)) then
            if (present(status)) status = 1
            allocate(cd(0))
            return
        end if
        cd = knn_dist(x, min_pts - 1, stat)
        if (present(status)) status = stat
    end function coredist

    function mrdist(x, min_pts, status) result(mrd)
        real(dp), intent(in) :: x(:, :) !! Data matrix with observations in rows.
        integer, intent(in) :: min_pts !! HDBSCAN MinPts including self.
        integer, intent(out), optional :: status !! Zero on success; nonzero for invalid input.
        real(dp), allocatable :: mrd(:, :)
        real(dp), allocatable :: d(:, :)
        real(dp), allocatable :: cd(:)
        integer :: i
        integer :: j
        integer :: n
        integer :: stat

        cd = coredist(x, min_pts, stat)
        if (stat /= 0) then
            if (present(status)) status = stat
            allocate(mrd(0, 0))
            return
        end if
        call pairwise_distance_matrix(x, d, status = stat)
        if (stat /= 0) then
            if (present(status)) status = stat
            allocate(mrd(0, 0))
            return
        end if
        n = size(x, 1)
        allocate(mrd(n, n))
        mrd = 0.0_dp
        do i = 1, n
            do j = i + 1, n
                mrd(i, j) = max(d(i, j), max(cd(i), cd(j)))
                mrd(j, i) = mrd(i, j)
            end do
        end do
        if (present(status)) status = 0
    end function mrdist

    function mst_dense(dmat, status) result(edges)
        real(dp), intent(in) :: dmat(:, :) !! Full symmetric distance matrix.
        integer, intent(out), optional :: status !! Zero on success; nonzero for invalid input.
        real(dp), allocatable :: edges(:, :) !! MST rows: from, to, weight with one-based vertex IDs.
        integer, allocatable :: parent(:)
        real(dp), allocatable :: weight(:)
        logical, allocatable :: visited(:)
        integer :: edge_count
        integer :: i
        integer :: n
        integer :: next_node
        integer :: node
        real(dp) :: next_weight

        n = size(dmat, 1)
        if (n < 1 .or. size(dmat, 2) /= n) then
            if (present(status)) status = 1
            allocate(edges(0, 3))
            return
        end if
        allocate(edges(max(0, n - 1), 3), parent(n), weight(n), visited(n))
        parent = 0
        weight = huge(1.0_dp)
        visited = .false.
        parent(1) = -1
        weight(1) = 0.0_dp
        next_node = 1
        edge_count = 0

        do while (next_node > 0)
            node = next_node
            next_node = 0
            next_weight = huge(1.0_dp)
            visited(node) = .true.
            if (node > 1) then
                edge_count = edge_count + 1
                edges(edge_count, 1) = real(node, dp)
                edges(edge_count, 2) = real(parent(node), dp)
                edges(edge_count, 3) = weight(node)
            end if
            do i = 2, n
                if (visited(i) .or. i == node) cycle
                if (dmat(node, i) < weight(i)) then
                    weight(i) = dmat(node, i)
                    parent(i) = node
                end if
                if (weight(i) < next_weight) then
                    next_weight = weight(i)
                    next_node = i
                end if
            end do
        end do
        if (edge_count /= n - 1) then
            if (present(status)) status = 2
        else if (present(status)) then
            status = 0
        end if
    end function mst_dense

    subroutine hdbscan(x, min_pts, result, status, cluster_selection_epsilon)
        real(dp), intent(in) :: x(:, :) !! Data matrix with observations in rows.
        integer, intent(in) :: min_pts !! Minimum cluster size and core-distance MinPts, including self.
        type(hdbscan_result), intent(out) :: result !! HDBSCAN hierarchy, flat EOM clustering, and scores.
        integer, intent(out), optional :: status !! Zero on success; nonzero for invalid input or extraction failure.
        real(dp), intent(in), optional :: cluster_selection_epsilon !! HDBSCAN(e) distance threshold; default zero.
        real(dp), allocatable :: d(:, :)
        real(dp), allocatable :: mrd(:, :)
        real(dp), allocatable :: mst(:, :)
        real(dp), allocatable :: cd(:)
        real(dp) :: selection_eps
        integer :: stat

        if (min_pts < 2 .or. min_pts > size(x, 1)) then
            if (present(status)) status = 1
            allocate(result%cluster(0), result%membership_prob(0), result%outlier_scores(0))
            return
        end if
        selection_eps = 0.0_dp
        if (present(cluster_selection_epsilon)) selection_eps = cluster_selection_epsilon
        if (selection_eps < 0.0_dp .or. .not. ieee_is_finite(selection_eps)) then
            if (present(status)) status = 1
            allocate(result%cluster(0), result%membership_prob(0), result%outlier_scores(0))
            return
        end if
        cd = coredist(x, min_pts, stat)
        if (stat /= 0) then
            if (present(status)) status = stat
            return
        end if
        call pairwise_distance_matrix(x, d, status = stat)
        if (stat /= 0) then
            if (present(status)) status = stat
            return
        end if
        call build_mrd(d, cd, mrd)
        mst = mst_dense(mrd, stat)
        if (stat /= 0) then
            if (present(status)) status = stat
            return
        end if
        result%min_pts = min_pts
        result%core_dist = cd
        result%mst = mst
        call mst_to_dendrogram(mst, result%hierarchy)
        call eom_from_hierarchy(result%hierarchy, min_pts, result%cluster, result%cluster_scores, &
            result%outlier_scores, stat, cluster_selection_epsilon = selection_eps)
        if (stat /= 0) then
            if (present(status)) status = stat
            return
        end if
        allocate(result%membership_prob(size(result%cluster)))
        call membership_from_core_distance(result%cluster, result%core_dist, result%membership_prob)
        if (present(status)) status = 0
    end subroutine hdbscan

    subroutine membership_from_core_distance(cluster, core_dist, probability)
        integer, intent(in) :: cluster(:) !! Flat HDBSCAN labels; zero denotes noise.
        real(dp), intent(in) :: core_dist(:) !! Core distance for each observation.
        real(dp), intent(inout) :: probability(:) !! Membership probabilities replaced by the upstream core-distance rule.
        integer :: c
        integer :: i
        integer :: max_cluster
        real(dp) :: max_core

        probability = 0.0_dp
        if (size(cluster) == 0) return
        max_cluster = maxval(cluster)
        do c = 1, max_cluster
            max_core = 0.0_dp
            do i = 1, size(cluster)
                if (cluster(i) == c) max_core = max(max_core, core_dist(i))
            end do
            if (max_core <= 0.0_dp) then
                do i = 1, size(cluster)
                    if (cluster(i) == c) probability(i) = 1.0_dp
                end do
            else
                do i = 1, size(cluster)
                    if (cluster(i) /= c) cycle
                    probability(i) = max(0.0_dp, min(1.0_dp, (max_core - core_dist(i)) / max_core))
                end do
            end if
        end do
    end subroutine membership_from_core_distance

    function glosh(x, k, status) result(score)
        real(dp), intent(in) :: x(:, :) !! Data matrix with observations in rows.
        integer, intent(in) :: k !! Neighborhood size used by the hierarchy; at least two.
        integer, intent(out), optional :: status !! Zero on success; nonzero for invalid input.
        real(dp), allocatable :: score(:)
        type(hdbscan_result) :: fit
        integer :: stat

        call hdbscan(x, k, fit, stat)
        if (stat /= 0) then
            allocate(score(0))
        else
            score = fit%outlier_scores
        end if
        if (present(status)) status = stat
    end function glosh

    subroutine extract_fosc(hierarchy, min_pts, cluster, stability, status, constraints, alpha, &
        prune_unstable, constraint_score, total_score)
        type(dendrogram_result), intent(in) :: hierarchy !! Single-link hierarchy to condense and extract by FOSC.
        integer, intent(in) :: min_pts !! Minimum condensed-tree cluster size; at least two.
        integer, allocatable, intent(out) :: cluster(:) !! FOSC flat labels; zero denotes noise.
        real(dp), allocatable, intent(out) :: stability(:) !! Raw stability scores for selected clusters.
        integer, intent(out), optional :: status !! Zero on success; nonzero for invalid hierarchy or constraints.
        integer, intent(in), optional :: constraints(:, :) !! Symmetric n-by-n matrix: 1 should-link, -1 should-not-link, 0 none.
        real(dp), intent(in), optional :: alpha !! Mixed-objective stability weight in [0, 1]; default zero.
        logical, intent(in), optional :: prune_unstable !! Prune unstable descendant branches when true; default false.
        real(dp), allocatable, intent(out), optional :: constraint_score(:) !! Constraint score for each selected cluster.
        real(dp), allocatable, intent(out), optional :: total_score(:) !! FOSC propagated objective score for each selected cluster.
        real(dp), allocatable :: outlier(:)
        real(dp), allocatable :: selected_constraint(:)
        real(dp), allocatable :: selected_total(:)
        real(dp) :: alpha_value
        logical :: prune_value
        integer :: stat

        alpha_value = 0.0_dp
        if (present(alpha)) alpha_value = alpha
        prune_value = .false.
        if (present(prune_unstable)) prune_value = prune_unstable
        if (min_pts < 2 .or. size(hierarchy%merge, 1) < 1) then
            allocate(cluster(0), stability(0))
            if (present(constraint_score)) allocate(constraint_score(0))
            if (present(total_score)) allocate(total_score(0))
            if (present(status)) status = 1
            return
        end if
        if (alpha_value < 0.0_dp .or. alpha_value > 1.0_dp .or. .not. ieee_is_finite(alpha_value)) then
            allocate(cluster(0), stability(0))
            if (present(constraint_score)) allocate(constraint_score(0))
            if (present(total_score)) allocate(total_score(0))
            if (present(status)) status = 1
            return
        end if
        if (present(constraints)) then
            call eom_from_hierarchy(hierarchy, min_pts, cluster, stability, outlier, stat, &
                constraints = constraints, alpha = alpha_value, prune_unstable = prune_value, &
                selected_constraint = selected_constraint, selected_total = selected_total)
        else
            call eom_from_hierarchy(hierarchy, min_pts, cluster, stability, outlier, stat, &
                prune_unstable = prune_value, selected_constraint = selected_constraint, &
                selected_total = selected_total)
        end if
        if (stat /= 0) then
            if (.not. allocated(cluster)) allocate(cluster(0))
            if (.not. allocated(stability)) allocate(stability(0))
            if (.not. allocated(selected_constraint)) allocate(selected_constraint(0))
            if (.not. allocated(selected_total)) allocate(selected_total(0))
        end if
        if (present(constraint_score)) constraint_score = selected_constraint
        if (present(total_score)) total_score = selected_total
        if (present(status)) status = stat
    end subroutine extract_fosc

    subroutine build_mrd(d, cd, mrd)
        real(dp), intent(in) :: d(:, :) !! Pairwise distance matrix.
        real(dp), intent(in) :: cd(:) !! Core distance for each observation.
        real(dp), allocatable, intent(out) :: mrd(:, :) !! Mutual-reachability distance matrix.
        integer :: i
        integer :: j
        integer :: n

        n = size(d, 1)
        allocate(mrd(n, n))
        mrd = 0.0_dp
        do i = 1, n
            do j = i + 1, n
                mrd(i, j) = max(d(i, j), max(cd(i), cd(j)))
                mrd(j, i) = mrd(i, j)
            end do
        end do
    end subroutine build_mrd

    subroutine mst_to_dendrogram(mst, tree)
        real(dp), intent(in) :: mst(:, :) !! MST rows: from, to, weight.
        type(dendrogram_result), intent(out) :: tree !! Single-link dendrogram in hclust-like merge form.
        integer, allocatable :: ord(:)
        integer, allocatable :: label(:)
        integer, allocatable :: parent(:)
        integer, allocatable :: size_set(:)
        integer :: i
        integer :: j
        integer :: n
        integer :: u
        integer :: v
        integer :: ru
        integer :: rv
        integer :: labu
        integer :: labv

        n = size(mst, 1) + 1
        allocate(ord(n - 1), label(n), parent(n), size_set(n))
        do i = 1, n - 1
            ord(i) = i
        end do
        call sort_mst_order(mst, ord)
        do i = 1, n
            parent(i) = i
            size_set(i) = 1
            label(i) = -i
        end do
        allocate(tree%merge(n - 1, 2), tree%height(n - 1), tree%order(n))
        do i = 1, n - 1
            j = ord(i)
            u = nint(mst(j, 1))
            v = nint(mst(j, 2))
            ru = uf_find(parent, u)
            rv = uf_find(parent, v)
            labu = label(ru)
            labv = label(rv)
            tree%merge(i, 1) = labu
            tree%merge(i, 2) = labv
            tree%height(i) = mst(j, 3)
            call uf_union(parent, size_set, ru, rv, i, label)
        end do
        i = 0
        call dendro_visit(tree%merge, n - 1, 1, tree%order, i)
        call dendro_visit(tree%merge, n - 1, 2, tree%order, i)
    end subroutine mst_to_dendrogram

    subroutine eom_from_hierarchy(tree, min_size, labels, selected_stability, outlier, status, &
        cluster_selection_epsilon, constraints, alpha, prune_unstable, selected_constraint, selected_total)
        type(dendrogram_result), intent(in) :: tree !! Single-link dendrogram in hclust-compatible merge form.
        integer, intent(in) :: min_size !! Minimum condensed-tree cluster size; at least two.
        integer, allocatable, intent(out) :: labels(:) !! Flat selected-cluster labels; zero denotes noise.
        real(dp), allocatable, intent(out) :: selected_stability(:) !! Raw stability scores for selected clusters.
        real(dp), allocatable, intent(out) :: outlier(:) !! GLOSH scores computed from condensed-cluster death distances.
        integer, intent(out) :: status !! Zero on success; nonzero for invalid input or constraints.
        real(dp), intent(in), optional :: cluster_selection_epsilon !! Distance threshold for HDBSCAN(e); default zero.
        integer, intent(in), optional :: constraints(:, :) !! Symmetric full pairwise FOSC constraint matrix.
        real(dp), intent(in), optional :: alpha !! Mixed FOSC objective weight in [0, 1]; default zero.
        logical, intent(in), optional :: prune_unstable !! Prune weak descendant branches when true.
        real(dp), allocatable, intent(out), optional :: selected_constraint(:) !! Constraint scores of selected clusters.
        real(dp), allocatable, intent(out), optional :: selected_total(:) !! Propagated FOSC objective scores of selected clusters.
        type(condensed_tree_data) :: ct
        logical, allocatable :: selected(:)
        real(dp), allocatable :: objective(:)
        real(dp), allocatable :: constraint_value(:)
        real(dp), allocatable :: vscore(:)
        real(dp) :: alpha_value
        real(dp) :: ignored_constraint
        real(dp) :: ignored_score
        real(dp) :: selection_eps
        real(dp) :: total_stability
        logical :: prune_value
        integer :: cid
        integer :: i
        integer :: k
        integer :: n_constraints
        integer :: nsel

        status = 0
        selection_eps = 0.0_dp
        if (present(cluster_selection_epsilon)) selection_eps = cluster_selection_epsilon
        alpha_value = 0.0_dp
        if (present(alpha)) alpha_value = alpha
        prune_value = .false.
        if (present(prune_unstable)) prune_value = prune_unstable
        if (min_size < 2 .or. size(tree%merge, 1) < 1 .or. selection_eps < 0.0_dp .or. &
            .not. ieee_is_finite(selection_eps) .or. alpha_value < 0.0_dp .or. alpha_value > 1.0_dp .or. &
            .not. ieee_is_finite(alpha_value)) then
            status = 1
            allocate(labels(0), selected_stability(0), outlier(0))
            if (present(selected_constraint)) allocate(selected_constraint(0))
            if (present(selected_total)) allocate(selected_total(0))
            return
        end if
        call build_condensed_tree(tree, min_size, ct, status)
        if (status /= 0) then
            allocate(labels(0), selected_stability(0), outlier(0))
            if (present(selected_constraint)) allocate(selected_constraint(0))
            if (present(selected_total)) allocate(selected_total(0))
            return
        end if
        call compute_glosh_scores(ct, outlier)
        allocate(selected(ct%ncluster), objective(ct%ncluster), constraint_value(ct%ncluster))
        selected = .false.
        objective = 0.0_dp
        constraint_value = 0.0_dp

        if (present(constraints)) then
            if (.not. valid_constraint_matrix(constraints, ct%n)) then
                status = 2
                allocate(labels(0), selected_stability(0))
                if (present(selected_constraint)) allocate(selected_constraint(0))
                if (present(selected_total)) allocate(selected_total(0))
                return
            end if
            n_constraints = count(constraints /= 0)
            if (n_constraints == 0) then
                ignored_score = select_unsupervised(ct, 1, .true., selection_eps, prune_value, selected, objective)
                if (ignored_score < 0.0_dp) selected(1) = .false.
            else
                selected = .false.
                objective = 0.0_dp
                ignored_score = select_unsupervised(ct, 1, .true., selection_eps, .false., selected, objective)
                total_stability = 0.0_dp
                do cid = 2, ct%ncluster
                    if (selected(cid)) total_stability = total_stability + ct%stability(cid)
                end do
                if (total_stability <= tiny(1.0_dp) .or. .not. ieee_is_finite(total_stability)) then
                    total_stability = 1.0_dp
                end if
                allocate(vscore(ct%ncluster))
                vscore = 0.0_dp
                do cid = 2, ct%ncluster
                    vscore(cid) = compute_virtual_cluster(ct, cid, constraints, n_constraints)
                end do
                selected = .false.
                objective = 0.0_dp
                constraint_value = vscore
                call select_semisupervised(ct, 1, .true., selection_eps, prune_value, constraints, &
                    n_constraints, alpha_value, total_stability, vscore, selected, objective, &
                    constraint_value, ignored_score, ignored_constraint)
            end if
        else
            ignored_score = select_unsupervised(ct, 1, .true., selection_eps, prune_value, selected, objective)
            if (ignored_score < 0.0_dp) selected(1) = .false.
        end if
        selected(1) = .false.

        nsel = count(selected)
        allocate(labels(ct%n), selected_stability(nsel))
        labels = 0
        if (present(selected_constraint)) allocate(selected_constraint(nsel))
        if (present(selected_total)) allocate(selected_total(nsel))
        k = 0
        do cid = 2, ct%ncluster
            if (.not. selected(cid)) cycle
            k = k + 1
            selected_stability(k) = ct%stability(cid)
            if (present(selected_constraint)) selected_constraint(k) = constraint_value(cid)
            if (present(selected_total)) selected_total(k) = objective(cid)
            do i = 1, ct%n
                if (cluster_descends_from(ct, ct%point_owner(i), cid)) labels(i) = k
            end do
        end do
    end subroutine eom_from_hierarchy

    subroutine build_condensed_tree(tree, min_size, ct, status)
        type(dendrogram_result), intent(in) :: tree !! Single-link dendrogram to condense using upstream dbscan rules.
        integer, intent(in) :: min_size !! Minimum cluster size used by HDBSCAN/FOSC.
        type(condensed_tree_data), intent(out) :: ct !! Condensed hierarchy, stability data, and point exit information.
        integer, intent(out) :: status !! Zero on success; nonzero for malformed merge labels or dimensions.
        integer, allocatable :: cl_tracker(:)
        integer, allocatable :: member_size(:)
        logical, allocatable :: processed(:)
        integer :: cid
        integer :: cl
        integer :: cr
        integer :: global_cid
        integer :: k
        integer :: lm
        integer :: max_cluster
        integer :: n
        integer :: rm
        integer :: singleton
        real(dp) :: noise_eps

        n = size(tree%merge, 1) + 1
        status = 0
        if (size(tree%merge, 2) /= 2 .or. size(tree%height) /= n - 1 .or. n < 2) then
            status = 1
            return
        end if
        max_cluster = max(2 * n, 2)
        ct%n = n
        allocate(ct%parent(max_cluster), ct%first_child(max_cluster), ct%next_sibling(max_cluster))
        allocate(ct%n_children(max_cluster))
        allocate(ct%birth_eps(max_cluster), ct%death_eps(max_cluster), ct%stability(max_cluster))
        allocate(ct%point_owner(n), ct%point_eps(n))
        allocate(cl_tracker(n - 1), member_size(n - 1), processed(max_cluster))
        ct%parent = 0
        ct%first_child = 0
        ct%next_sibling = 0
        ct%n_children = 0
        ct%birth_eps = 0.0_dp
        ct%death_eps = -1.0_dp
        ct%stability = 0.0_dp
        ct%point_owner = 0
        ct%point_eps = 0.0_dp
        cl_tracker = 1
        member_size = 0
        processed = .false.

        do k = 1, n - 1
            lm = tree%merge(k, 1)
            rm = tree%merge(k, 2)
            if (.not. valid_merge_label(lm, k, n) .or. .not. valid_merge_label(rm, k, n)) then
                status = 1
                return
            end if
            if (lm < 0 .and. rm < 0) then
                member_size(k) = 2
            else if (lm < 0 .or. rm < 0) then
                if (lm > 0) then
                    member_size(k) = member_size(lm) + 1
                else
                    member_size(k) = member_size(rm) + 1
                end if
            else
                member_size(k) = member_size(lm) + member_size(rm)
            end if
        end do

        ct%birth_eps(1) = tree%height(n - 1)
        global_cid = 1
        do k = n - 1, 1, -1
            cid = cl_tracker(k)
            lm = tree%merge(k, 1)
            rm = tree%merge(k, 2)
            if (lm < 0 .and. rm < 0) then
                if (processed(cid)) then
                    noise_eps = ct%death_eps(cid)
                else
                    noise_eps = tree%height(k)
                end if
                call record_point_exit(-lm, cid, noise_eps, ct)
                call record_point_exit(-rm, cid, noise_eps, ct)
                if (.not. processed(cid)) then
                    if (ct%death_eps(cid) < 0.0_dp) then
                        ct%death_eps(cid) = tree%height(k)
                    else
                        ct%death_eps(cid) = min(ct%death_eps(cid), tree%height(k))
                    end if
                end if
            else if (lm < 0 .or. rm < 0) then
                if (lm < 0) then
                    singleton = -lm
                    cl_tracker(rm) = cid
                else
                    singleton = -rm
                    cl_tracker(lm) = cid
                end if
                if (processed(cid)) then
                    noise_eps = ct%death_eps(cid)
                else
                    noise_eps = tree%height(k)
                end if
                call record_point_exit(singleton, cid, noise_eps, ct)
            else if (member_size(lm) >= min_size .and. member_size(rm) >= min_size) then
                ct%death_eps(cid) = tree%height(k)
                ct%n_children(cid) = member_size(lm) + member_size(rm)
                processed(cid) = .true.
                global_cid = global_cid + 1
                cl = global_cid
                global_cid = global_cid + 1
                cr = global_cid
                if (cr > max_cluster) then
                    status = 3
                    return
                end if
                ct%first_child(cid) = cl
                ct%next_sibling(cl) = cr
                ct%parent(cl) = cid
                ct%parent(cr) = cid
                cl_tracker(lm) = cl
                cl_tracker(rm) = cr
                ct%birth_eps(cl) = tree%height(k)
                ct%birth_eps(cr) = tree%height(k)
                ct%death_eps(cl) = tree%height(lm)
                ct%death_eps(cr) = tree%height(rm)
            else
                cl_tracker(lm) = cid
                cl_tracker(rm) = cid
            end if
        end do
        ct%ncluster = global_cid
        if (any(ct%point_owner == 0)) then
            status = 4
            return
        end if
        do cid = 1, ct%ncluster
            if (ct%death_eps(cid) < 0.0_dp) ct%death_eps(cid) = ct%birth_eps(cid)
        end do
        do k = 1, n
            cid = ct%point_owner(k)
            ct%stability(cid) = ct%stability(cid) + inverse_height(ct%point_eps(k)) - &
                inverse_height(ct%birth_eps(cid))
        end do
        do cid = 1, ct%ncluster
            if (ct%n_children(cid) == 0) cycle
            ct%stability(cid) = ct%stability(cid) + real(ct%n_children(cid), dp) * &
                (inverse_height(ct%death_eps(cid)) - inverse_height(ct%birth_eps(cid)))
        end do
    end subroutine build_condensed_tree

    pure logical function valid_merge_label(label, row, n) result(valid)
        integer, intent(in) :: label !! hclust merge label to validate.
        integer, intent(in) :: row !! Current one-based merge row.
        integer, intent(in) :: n !! Number of dendrogram leaves.

        if (label < 0) then
            valid = -label >= 1 .and. -label <= n
        else
            valid = label >= 1 .and. label < row
        end if
    end function valid_merge_label

    pure subroutine record_point_exit(point, cid, eps, ct)
        integer, intent(in) :: point !! One-based observation index leaving its current condensed cluster.
        integer, intent(in) :: cid !! Condensed-cluster identifier owning the observation at exit.
        real(dp), intent(in) :: eps !! Dissimilarity distance at which the observation exits the cluster.
        type(condensed_tree_data), intent(inout) :: ct !! Condensed-tree point-exit arrays to update.

        ct%point_owner(point) = cid
        ct%point_eps(point) = eps
    end subroutine record_point_exit

    pure logical function cluster_descends_from(ct, candidate, ancestor) result(descends)
        type(condensed_tree_data), intent(in) :: ct !! Condensed hierarchy providing parent links.
        integer, intent(in) :: candidate !! Candidate cluster identifier.
        integer, intent(in) :: ancestor !! Cluster identifier tested as an ancestor or self.
        integer :: current

        descends = .false.
        current = candidate
        do while (current > 0)
            if (current == ancestor) then
                descends = .true.
                return
            end if
            current = ct%parent(current)
        end do
    end function cluster_descends_from

    pure subroutine compute_glosh_scores(ct, outlier)
        type(condensed_tree_data), intent(in) :: ct !! Condensed hierarchy and per-point exit distances.
        real(dp), allocatable, intent(out) :: outlier(:) !! Exact upstream-style GLOSH score for each observation.
        integer :: cid
        integer :: i
        real(dp) :: eps_max
        real(dp) :: value

        allocate(outlier(ct%n))
        outlier = 0.0_dp
        do i = 1, ct%n
            cid = ct%point_owner(i)
            eps_max = leaf_death_epsilon(ct, cid)
            if (ct%point_eps(i) <= tiny(1.0_dp)) then
                value = 0.0_dp
            else
                value = 1.0_dp - eps_max / ct%point_eps(i)
                if (ieee_is_nan(value)) value = 0.0_dp
            end if
            outlier(i) = value
        end do
    end subroutine compute_glosh_scores

    pure recursive real(dp) function leaf_death_epsilon(ct, cid) result(eps_min)
        type(condensed_tree_data), intent(in) :: ct !! Condensed hierarchy containing death distances.
        integer, intent(in) :: cid !! Cluster whose leaf-descendant death distance is requested.
        integer :: child
        real(dp) :: child_eps

        child = ct%first_child(cid)
        if (child == 0) then
            eps_min = ct%death_eps(cid)
            return
        end if
        eps_min = huge(1.0_dp)
        do while (child /= 0)
            child_eps = leaf_death_epsilon(ct, child)
            eps_min = min(eps_min, child_eps)
            child = ct%next_sibling(child)
        end do
    end function leaf_death_epsilon

    recursive real(dp) function select_unsupervised(ct, cid, is_root, selection_eps, prune_unstable, &
        selected, objective) result(score)
        type(condensed_tree_data), intent(in) :: ct !! Condensed hierarchy and raw stability scores.
        integer, intent(in) :: cid !! Cluster identifier currently visited.
        logical, intent(in) :: is_root !! True only for the synthetic root cluster, which is never selected.
        real(dp), intent(in) :: selection_eps !! HDBSCAN(e) distance threshold; zero disables epsilon pruning.
        logical, intent(in) :: prune_unstable !! Apply upstream unstable-descendant pruning when true.
        logical, intent(inout) :: selected(:) !! Salient-cluster flags updated during traversal.
        real(dp), intent(inout) :: objective(:) !! Propagated stability objective for each cluster.
        integer :: child
        integer :: child_id(2)
        integer :: nchild
        real(dp) :: child_score(2)
        real(dp) :: children_total
        logical :: keep_children

        child = ct%first_child(cid)
        if (child == 0) then
            if (.not. is_root) selected(cid) = .true.
            score = ct%stability(cid)
            objective(cid) = score
            return
        end if
        child_id = 0
        child_score = 0.0_dp
        children_total = 0.0_dp
        nchild = 0
        do while (child /= 0)
            nchild = nchild + 1
            child_id(nchild) = child
            child_score(nchild) = select_unsupervised(ct, child, .false., selection_eps, prune_unstable, &
                selected, objective)
            children_total = children_total + child_score(nchild)
            child = ct%next_sibling(child)
        end do
        keep_children = ct%stability(cid) < children_total
        if (.not. is_root .and. ct%death_eps(cid) < selection_eps) keep_children = .false.
        if (.not. keep_children .and. .not. is_root) then
            call clear_descendants(cid, ct%first_child, ct%next_sibling, selected)
            selected(cid) = .true.
            score = ct%stability(cid)
        else
            selected(cid) = .false.
            score = children_total
            if (keep_children .and. prune_unstable) then
                if (any(child_score(1:nchild) >= ct%stability(cid))) then
                    do child = 1, nchild
                        if (child_score(child) < ct%stability(cid)) then
                            call clear_descendants(child_id(child), ct%first_child, ct%next_sibling, selected)
                        end if
                    end do
                end if
            end if
        end if
        objective(cid) = score
    end function select_unsupervised

    recursive subroutine select_semisupervised(ct, cid, is_root, selection_eps, prune_unstable, constraints, &
        n_constraints, alpha, total_stability, vscore, selected, objective, constraint_value, score, constraint_score)
        type(condensed_tree_data), intent(in) :: ct !! Condensed hierarchy and raw stability scores.
        integer, intent(in) :: cid !! Cluster identifier currently visited.
        logical, intent(in) :: is_root !! True only for the synthetic root cluster.
        real(dp), intent(in) :: selection_eps !! Optional HDBSCAN(e) distance threshold.
        logical, intent(in) :: prune_unstable !! Apply upstream unstable-descendant pruning when true.
        integer, intent(in) :: constraints(:, :) !! Symmetric pairwise FOSC constraints.
        integer, intent(in) :: n_constraints !! Number of directed nonzero entries in constraints.
        real(dp), intent(in) :: alpha !! Stability weight in the mixed objective, between zero and one.
        real(dp), intent(in) :: total_stability !! Unsupervised optimum used to normalize stability.
        real(dp), intent(inout) :: vscore(:) !! Cluster constraint scores, updated as child solutions improve them.
        logical, intent(inout) :: selected(:) !! Salient-cluster flags updated during traversal.
        real(dp), intent(inout) :: objective(:) !! Propagated mixed-objective score for each cluster.
        real(dp), intent(inout) :: constraint_value(:) !! Propagated constraint score for each cluster.
        real(dp), intent(out) :: score !! Objective score returned to the parent cluster.
        real(dp), intent(out) :: constraint_score !! Constraint score returned to the parent cluster.
        integer :: child
        integer :: child_id(2)
        integer :: nchild
        real(dp) :: child_constraint(2)
        real(dp) :: child_score(2)
        real(dp) :: new_constraint
        real(dp) :: new_stability
        real(dp) :: old_constraint
        real(dp) :: old_stability
        logical :: keep_children

        child = ct%first_child(cid)
        if (child == 0) then
            if (.not. is_root) selected(cid) = .true.
            score = ct%stability(cid)
            constraint_score = vscore(cid)
            objective(cid) = score
            constraint_value(cid) = constraint_score
            return
        end if
        child_id = 0
        child_score = 0.0_dp
        child_constraint = 0.0_dp
        nchild = 0
        do while (child /= 0)
            nchild = nchild + 1
            child_id(nchild) = child
            call select_semisupervised(ct, child, .false., selection_eps, prune_unstable, constraints, &
                n_constraints, alpha, total_stability, vscore, selected, objective, constraint_value, &
                child_score(nchild), child_constraint(nchild))
            child = ct%next_sibling(child)
        end do

        old_stability = ct%stability(cid) / total_stability
        new_stability = sum(child_score(1:nchild)) / total_stability
        old_constraint = vscore(cid)
        new_constraint = sum(child_constraint(1:nchild)) + &
            compute_virtual_owner(ct, cid, constraints, n_constraints)
        keep_children = .true.
        if (old_constraint < new_constraint .and. .not. is_root) then
            vscore(cid) = new_constraint
            score = alpha * new_stability + (1.0_dp - alpha) * new_constraint
            constraint_score = new_constraint
        else if (old_constraint > new_constraint .and. .not. is_root) then
            score = alpha * old_stability + (1.0_dp - alpha) * old_constraint
            constraint_score = old_constraint
            keep_children = .false.
        else
            if (old_stability < new_stability) then
                score = new_stability / total_stability
            else
                score = old_stability / total_stability
                keep_children = .false.
            end if
            constraint_score = old_constraint
            vscore(cid) = old_constraint
        end if
        if (.not. is_root .and. ct%death_eps(cid) < selection_eps) keep_children = .false.

        if (.not. keep_children .and. .not. is_root) then
            call clear_descendants(cid, ct%first_child, ct%next_sibling, selected)
            selected(cid) = .true.
        else if (keep_children .and. prune_unstable) then
            if (any(child_score(1:nchild) >= old_stability)) then
                do child = 1, nchild
                    if (child_score(child) < old_stability) then
                        call clear_descendants(child_id(child), ct%first_child, ct%next_sibling, selected)
                    end if
                end do
            end if
        end if
        objective(cid) = score
        constraint_value(cid) = constraint_score
    end subroutine select_semisupervised

    pure real(dp) function compute_virtual_cluster(ct, cid, constraints, n_constraints) result(score)
        type(condensed_tree_data), intent(in) :: ct !! Condensed hierarchy defining full cluster membership.
        integer, intent(in) :: cid !! Cluster whose full subtree membership is scored.
        integer, intent(in) :: constraints(:, :) !! Symmetric pairwise constraint matrix.
        integer, intent(in) :: n_constraints !! Directed count of nonzero constraints used for normalization.
        logical :: in_cluster(ct%n)
        integer :: i

        do i = 1, ct%n
            in_cluster(i) = cluster_descends_from(ct, ct%point_owner(i), cid)
        end do
        score = compute_virtual_set(in_cluster, constraints, n_constraints)
    end function compute_virtual_cluster

    pure real(dp) function compute_virtual_owner(ct, cid, constraints, n_constraints) result(score)
        type(condensed_tree_data), intent(in) :: ct !! Condensed hierarchy defining each cluster's own noise points.
        integer, intent(in) :: cid !! Cluster whose directly contained points are scored.
        integer, intent(in) :: constraints(:, :) !! Symmetric pairwise constraint matrix.
        integer, intent(in) :: n_constraints !! Directed count of nonzero constraints used for normalization.
        logical :: in_set(ct%n)
        integer :: i

        do i = 1, ct%n
            in_set(i) = ct%point_owner(i) == cid
        end do
        score = compute_virtual_set(in_set, constraints, n_constraints)
    end function compute_virtual_owner

    pure real(dp) function compute_virtual_set(in_set, constraints, n_constraints) result(score)
        logical, intent(in) :: in_set(:) !! Membership mask for the point set used as a FOSC virtual node.
        integer, intent(in) :: constraints(:, :) !! Pairwise constraints with values -1, 0, or 1.
        integer, intent(in) :: n_constraints !! Directed number of nonzero constraints; must be positive.
        integer :: i
        integer :: j
        integer :: satisfied

        satisfied = 0
        do i = 1, size(in_set)
            if (.not. in_set(i)) cycle
            do j = 1, size(in_set)
                if (constraints(i, j) > 0) then
                    if (in_set(j)) satisfied = satisfied + 1
                else if (constraints(i, j) < 0) then
                    if (.not. in_set(j)) satisfied = satisfied + 1
                end if
            end do
        end do
        score = real(satisfied, dp) / real(n_constraints, dp)
    end function compute_virtual_set

    pure logical function valid_constraint_matrix(constraints, n) result(valid)
        integer, intent(in) :: constraints(:, :) !! Candidate full pairwise matrix of FOSC constraints.
        integer, intent(in) :: n !! Number of observations expected in each matrix dimension.
        integer :: i
        integer :: j

        valid = size(constraints, 1) == n .and. size(constraints, 2) == n
        if (.not. valid) return
        do i = 1, n
            if (constraints(i, i) /= 0) then
                valid = .false.
                return
            end if
            do j = i + 1, n
                if (abs(constraints(i, j)) > 1 .or. abs(constraints(j, i)) > 1) then
                    valid = .false.
                    return
                end if
                if (constraints(i, j) /= constraints(j, i)) then
                    valid = .false.
                    return
                end if
            end do
        end do
    end function valid_constraint_matrix

    pure recursive subroutine clear_descendants(cid, first_child, next_sibling, selected)
        integer, intent(in) :: cid !! Parent cluster whose descendants are deselected.
        integer, intent(in) :: first_child(:) !! First child of each condensed cluster.
        integer, intent(in) :: next_sibling(:) !! Next sibling in each parent's child list.
        logical, intent(inout) :: selected(:) !! Selection flags to clear below cid.
        integer :: child

        child = first_child(cid)
        do while (child /= 0)
            selected(child) = .false.
            call clear_descendants(child, first_child, next_sibling, selected)
            child = next_sibling(child)
        end do
    end subroutine clear_descendants

    pure real(dp) function inverse_height(h) result(lambda)
        real(dp), intent(in) :: h !! Nonnegative hierarchy distance; zero maps to a large finite reciprocal.

        if (h <= tiny(1.0_dp)) then
            lambda = 1.0_dp / tiny(1.0_dp)
        else
            lambda = 1.0_dp / h
        end if
    end function inverse_height

    subroutine sort_mst_order(mst, order)
        real(dp), intent(in) :: mst(:, :) !! MST rows whose weight is column three.
        integer, intent(inout) :: order(:) !! Row indices sorted by ascending weight with deterministic ties.
        integer :: i
        integer :: j
        integer :: key
        logical :: move

        do i = 2, size(order)
            key = order(i)
            j = i - 1
            do while (j >= 1)
                if (mst(order(j), 3) > mst(key, 3)) then
                    move = .true.
                else if (mst(order(j), 3) < mst(key, 3)) then
                    move = .false.
                else
                    move = order(j) > key
                end if
                if (.not. move) exit
                order(j + 1) = order(j)
                j = j - 1
            end do
            order(j + 1) = key
        end do
    end subroutine sort_mst_order

    recursive integer function uf_find(parent, x) result(root)
        integer, intent(inout) :: parent(:) !! Union-find parent vector.
        integer, intent(in) :: x !! Vertex whose set representative is requested.
        if (parent(x) == x) then
            root = x
        else
            parent(x) = uf_find(parent, parent(x))
            root = parent(x)
        end if
    end function uf_find

    subroutine uf_union(parent, size_set, ra, rb, merge_id, label)
        integer, intent(inout) :: parent(:) !! Union-find parent vector.
        integer, intent(inout) :: size_set(:) !! Union-find set sizes.
        integer, intent(in) :: ra !! Root of the first set.
        integer, intent(in) :: rb !! Root of the second set.
        integer, intent(in) :: merge_id !! One-based dendrogram merge row being created.
        integer, intent(inout) :: label(:) !! hclust-style label associated with each current union-find root.
        integer :: root
        integer :: other

        if (size_set(ra) >= size_set(rb)) then
            root = ra
            other = rb
        else
            root = rb
            other = ra
        end if
        parent(other) = root
        size_set(root) = size_set(ra) + size_set(rb)
        label(root) = merge_id
    end subroutine uf_union

    recursive subroutine dendro_visit(merge, row, col, order, pos)
        integer, intent(in) :: merge(:, :) !! hclust-style merge matrix.
        integer, intent(in) :: row !! Merge row to visit.
        integer, intent(in) :: col !! Child column, one or two.
        integer, intent(inout) :: order(:) !! Output leaf ordering.
        integer, intent(inout) :: pos !! Number of leaves already emitted.
        integer :: value

        value = merge(row, col)
        if (value < 0) then
            pos = pos + 1
            order(pos) = -value
        else
            call dendro_visit(merge, value, 1, order, pos)
            call dendro_visit(merge, value, 2, order, pos)
        end if
    end subroutine dendro_visit

end module dbscan_hdbscan
