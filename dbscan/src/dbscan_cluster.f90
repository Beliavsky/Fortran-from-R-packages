module dbscan_cluster
    use dbscan_kinds, only : dp
    use dbscan_types, only : clustering_result, frnn_result
    use dbscan_nn, only : frnn
    implicit none
    private

    public :: dbscan, dbscan_from_frnn, is_corepoint

contains

    subroutine dbscan(x, eps, min_pts, result, border_points, weights, status)
        real(dp), intent(in) :: x(:, :) !! Numeric data matrix with observations in rows.
        real(dp), intent(in) :: eps !! Euclidean DBSCAN neighborhood radius.
        integer, intent(in) :: min_pts !! Minimum neighborhood weight/count including the point itself.
        type(clustering_result), intent(out) :: result !! DBSCAN cluster labels; zero denotes noise.
        logical, intent(in), optional :: border_points !! Assign border points when true; default true.
        real(dp), intent(in), optional :: weights(:) !! Optional nonnegative point weights used instead of counts.
        integer, intent(out), optional :: status !! Zero on success; nonzero for invalid input.
        type(frnn_result) :: nn
        integer :: stat

        call frnn(x, eps, nn, stat)
        if (stat /= 0) then
            if (present(status)) status = stat
            allocate(result%cluster(0))
            return
        end if
        call dbscan_from_frnn(nn, min_pts, result, border_points, weights, status)
        result%eps = eps
    end subroutine dbscan

    subroutine dbscan_from_frnn(nn, min_pts, result, border_points, weights, status)
        type(frnn_result), intent(in) :: nn !! Fixed-radius neighbors excluding self.
        integer, intent(in) :: min_pts !! Minimum neighborhood weight/count including the point itself.
        type(clustering_result), intent(out) :: result !! DBSCAN cluster labels; zero denotes noise.
        logical, intent(in), optional :: border_points !! Assign border points when true; default true.
        real(dp), intent(in), optional :: weights(:) !! Optional point weights; self is included in every neighborhood.
        integer, intent(out), optional :: status !! Zero on success; nonzero for invalid input.
        logical, allocatable :: visited(:)
        integer, allocatable :: queue(:)
        logical :: use_border
        integer :: cluster_id
        integer :: i
        integer :: j
        integer :: p
        integer :: qhead
        integer :: qtail
        integer :: n
        real(dp) :: neighborhood_weight

        if (present(status)) status = 0
        n = nn%n
        if (min_pts < 1 .or. n < 1) then
            if (present(status)) status = 1
            allocate(result%cluster(0))
            return
        end if
        if (present(weights)) then
            if (size(weights) /= n) then
                if (present(status)) status = 2
                allocate(result%cluster(0))
                return
            end if
        end if

        use_border = .true.
        if (present(border_points)) use_border = border_points
        allocate(result%cluster(n), visited(n), queue(max(1, n * n)))
        result%cluster = 0
        visited = .false.
        cluster_id = 0
        result%min_pts = min_pts
        result%eps = nn%eps

        do i = 1, n
            if (visited(i)) cycle
            neighborhood_weight = local_weight(nn, i, weights)
            if (neighborhood_weight < real(min_pts, dp)) cycle

            cluster_id = cluster_id + 1
            result%cluster(i) = cluster_id
            visited(i) = .true.
            qhead = 1
            qtail = 0
            do p = nn%offset(i), nn%offset(i + 1) - 1
                qtail = qtail + 1
                queue(qtail) = nn%id(p)
            end do

            do while (qhead <= qtail)
                j = queue(qhead)
                qhead = qhead + 1
                if (visited(j)) cycle
                visited(j) = .true.
                neighborhood_weight = local_weight(nn, j, weights)
                if (neighborhood_weight >= real(min_pts, dp)) then
                    do p = nn%offset(j), nn%offset(j + 1) - 1
                        qtail = qtail + 1
                        queue(qtail) = nn%id(p)
                    end do
                end if
                if (neighborhood_weight >= real(min_pts, dp) .or. use_border) then
                    result%cluster(j) = cluster_id
                end if
            end do
        end do
    end subroutine dbscan_from_frnn

    function is_corepoint(x, eps, min_pts, weights, status) result(core)
        real(dp), intent(in) :: x(:, :) !! Numeric data matrix with observations in rows.
        real(dp), intent(in) :: eps !! Euclidean neighborhood radius.
        integer, intent(in) :: min_pts !! Minimum neighborhood weight/count including self.
        real(dp), intent(in), optional :: weights(:) !! Optional point weights used instead of counts.
        integer, intent(out), optional :: status !! Zero on success; nonzero for invalid input.
        logical, allocatable :: core(:)
        type(frnn_result) :: nn
        integer :: i
        integer :: stat

        call frnn(x, eps, nn, stat)
        if (present(status)) status = stat
        if (stat /= 0) then
            allocate(core(0))
            return
        end if
        allocate(core(nn%n))
        do i = 1, nn%n
            core(i) = local_weight(nn, i, weights) >= real(min_pts, dp)
        end do
    end function is_corepoint

    pure real(dp) function local_weight(nn, i, weights) result(value)
        type(frnn_result), intent(in) :: nn !! Fixed-radius neighbors excluding self.
        integer, intent(in) :: i !! One-based point whose neighborhood is measured.
        real(dp), intent(in), optional :: weights(:) !! Optional point weights.
        integer :: p

        if (present(weights)) then
            value = weights(i)
            do p = nn%offset(i), nn%offset(i + 1) - 1
                value = value + weights(nn%id(p))
            end do
        else
            value = real(1 + nn%offset(i + 1) - nn%offset(i), dp)
        end if
    end function local_weight

end module dbscan_cluster
