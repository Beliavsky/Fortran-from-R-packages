module dbscan_snn
    use dbscan_kinds, only : dp
    use dbscan_types, only : knn_result, frnn_result, clustering_result
    use dbscan_nn, only : knn
    use dbscan_cluster, only : dbscan_from_frnn
    use dbscan_graph, only : relabel_components
    implicit none
    private

    public :: snn, snnclust, jpclust

contains

    subroutine snn(x, k, result, jp, status)
        real(dp), intent(in) :: x(:, :) !! Data matrix used to build the k-nearest-neighbor graph.
        integer, intent(in) :: k !! Number of nearest neighbors.
        type(knn_result), intent(out) :: result !! kNN result plus shared-neighbor counts.
        logical, intent(in), optional :: jp !! Require mutual-neighbor links for similarity when true.
        integer, intent(out), optional :: status !! Zero on success; nonzero for invalid input.
        logical :: use_jp
        integer :: i
        integer :: j
        integer :: p
        integer :: q
        integer :: shared
        integer :: stat

        call knn(x, k, result, stat)
        if (present(status)) status = stat
        if (stat /= 0) return
        use_jp = .false.
        if (present(jp)) use_jp = jp
        allocate(result%shared(size(result%id, 1), k))
        result%shared = 0

        do i = 1, size(result%id, 1)
            do p = 1, k
                j = result%id(i, p)
                if (use_jp .and. .not. any(result%id(j, :) == i)) cycle
                shared = 0
                do q = 1, k
                    if (any(result%id(j, :) == result%id(i, q))) shared = shared + 1
                end do
                if (any(result%id(j, :) == i)) shared = shared + 1
                result%shared(i, p) = shared
            end do
        end do
    end subroutine snn

    subroutine snnclust(x, k, eps, min_pts, result, border_points, status)
        real(dp), intent(in) :: x(:, :) !! Data matrix used to construct the shared-neighbor graph.
        integer, intent(in) :: k !! kNN sparsification size.
        integer, intent(in) :: eps !! Minimum shared-neighbor count for a graph edge.
        integer, intent(in) :: min_pts !! Minimum number of qualifying neighbors including self for a core point.
        type(clustering_result), intent(out) :: result !! Shared-nearest-neighbor clustering labels.
        logical, intent(in), optional :: border_points !! Assign border points when true; default true.
        integer, intent(out), optional :: status !! Zero on success; nonzero for invalid input.
        type(knn_result) :: nn
        type(frnn_result) :: graph
        integer, allocatable :: counts(:)
        integer :: i
        integer :: p
        integer :: pos
        integer :: stat

        call snn(x, k, nn, .true., stat)
        if (stat /= 0) then
            if (present(status)) status = stat
            allocate(result%cluster(0))
            return
        end if
        allocate(counts(size(x, 1)))
        counts = 0
        do i = 1, size(x, 1)
            counts(i) = count(nn%shared(i, :) >= eps)
        end do
        graph%n = size(x, 1)
        graph%eps = real(eps, dp)
        allocate(graph%offset(graph%n + 1))
        graph%offset(1) = 1
        do i = 1, graph%n
            graph%offset(i + 1) = graph%offset(i) + counts(i)
        end do
        allocate(graph%id(sum(counts)), graph%dist(sum(counts)))
        do i = 1, graph%n
            pos = graph%offset(i)
            do p = 1, k
                if (nn%shared(i, p) < eps) cycle
                graph%id(pos) = nn%id(i, p)
                graph%dist(pos) = real(k - nn%shared(i, p), dp)
                pos = pos + 1
            end do
        end do
        call dbscan_from_frnn(graph, min_pts, result, border_points, status = status)
        result%k = k
        result%eps = real(eps, dp)
    end subroutine snnclust

    subroutine jpclust(x, k, kt, result, status)
        real(dp), intent(in) :: x(:, :) !! Data matrix used to construct the nearest-neighbor graph.
        integer, intent(in) :: k !! kNN sparsification size.
        integer, intent(in) :: kt !! Minimum shared-neighbor threshold including mutual endpoints.
        type(clustering_result), intent(out) :: result !! Jarvis-Patrick connected-component labels.
        integer, intent(out), optional :: status !! Zero on success; nonzero for invalid input.
        type(knn_result) :: nn
        integer, allocatable :: label(:)
        integer :: i
        integer :: j
        integer :: p
        integer :: q
        integer :: shared
        integer :: newlabel
        integer :: oldlabel
        integer :: stat

        if (kt < 1 .or. kt > k) then
            if (present(status)) status = 1
            allocate(result%cluster(0))
            return
        end if
        call knn(x, k, nn, stat)
        if (stat /= 0) then
            if (present(status)) status = stat
            allocate(result%cluster(0))
            return
        end if
        allocate(label(size(x, 1)))
        do i = 1, size(label)
            label(i) = i
        end do
        do i = 1, size(label)
            do p = 1, k
                j = nn%id(i, p)
                if (j < i) cycle
                if (.not. any(nn%id(j, :) == i)) cycle
                shared = 1
                do q = 1, k
                    if (any(nn%id(j, :) == nn%id(i, q))) shared = shared + 1
                end do
                if (shared < kt .or. label(i) == label(j)) cycle
                newlabel = min(label(i), label(j))
                oldlabel = max(label(i), label(j))
                where (label == oldlabel) label = newlabel
            end do
        end do
        call relabel_components(label)
        result%cluster = label
        result%k = k
        result%kt = kt
        if (present(status)) status = 0
    end subroutine jpclust

end module dbscan_snn
