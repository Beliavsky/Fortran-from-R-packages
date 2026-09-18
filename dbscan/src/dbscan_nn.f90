module dbscan_nn
    use dbscan_kinds, only : dp
    use dbscan_types, only : knn_result, frnn_result
    use dbscan_utils, only : contains_nan, sort_pairs
    implicit none
    private

    public :: knn, knn_query, frnn, frnn_query, knn_from_dist, frnn_from_dist
    public :: adjacencylist_knn, knn_dist

contains

    subroutine knn(x, k, result, status)
        real(dp), intent(in) :: x(:, :) !! Data matrix with observations in rows and features in columns.
        integer, intent(in) :: k !! Number of neighbors; must satisfy 1 <= k < nrow(x).
        type(knn_result), intent(out) :: result !! Sorted exact Euclidean k-nearest neighbors, excluding self.
        integer, intent(out), optional :: status !! Zero on success; nonzero for invalid input.
        integer :: i
        integer :: j
        integer :: m
        integer :: n
        real(dp), allocatable :: d(:)
        integer, allocatable :: id(:)

        if (present(status)) status = 0
        n = size(x, 1)
        if (n < 2 .or. size(x, 2) < 1 .or. k < 1 .or. k >= n .or. contains_nan(x)) then
            if (present(status)) status = 1
            result%k = 0
            allocate(result%id(0, 0), result%dist(0, 0))
            return
        end if

        result%k = k
        allocate(result%id(n, k), result%dist(n, k))
        allocate(d(n - 1), id(n - 1))
        do i = 1, n
            m = 0
            do j = 1, n
                if (j == i) cycle
                m = m + 1
                id(m) = j
                d(m) = sqrt(sum((x(i, :) - x(j, :)) ** 2))
            end do
            call sort_pairs(d, id, m)
            result%id(i, :) = id(1:k)
            result%dist(i, :) = d(1:k)
        end do
    end subroutine knn

    subroutine knn_query(x, query, k, result, status)
        real(dp), intent(in) :: x(:, :) !! Reference data matrix with observations in rows.
        real(dp), intent(in) :: query(:, :) !! Query matrix with the same number of columns as x.
        integer, intent(in) :: k !! Number of reference neighbors to return for each query.
        type(knn_result), intent(out) :: result !! Sorted exact Euclidean neighbors for each query row.
        integer, intent(out), optional :: status !! Zero on success; nonzero for invalid input.
        integer :: i
        integer :: j
        integer :: n
        integer :: nq
        real(dp), allocatable :: d(:)
        integer, allocatable :: id(:)

        if (present(status)) status = 0
        n = size(x, 1)
        nq = size(query, 1)
        if (size(x, 2) < 1 .or. size(query, 2) /= size(x, 2) .or. k < 1 .or. k > n) then
            if (present(status)) status = 1
            allocate(result%id(0, 0), result%dist(0, 0))
            result%k = 0
            return
        end if
        if (contains_nan(x) .or. contains_nan(query)) then
            if (present(status)) status = 2
            allocate(result%id(0, 0), result%dist(0, 0))
            result%k = 0
            return
        end if

        result%k = k
        allocate(result%id(nq, k), result%dist(nq, k))
        allocate(d(n), id(n))
        do i = 1, nq
            do j = 1, n
                id(j) = j
                d(j) = sqrt(sum((query(i, :) - x(j, :)) ** 2))
            end do
            call sort_pairs(d, id, n)
            result%id(i, :) = id(1:k)
            result%dist(i, :) = d(1:k)
        end do
    end subroutine knn_query

    subroutine knn_from_dist(dmat, k, result, status)
        real(dp), intent(in) :: dmat(:, :) !! Full symmetric distance matrix with zero diagonal.
        integer, intent(in) :: k !! Number of neighbors; must satisfy 1 <= k < size(dmat,1).
        type(knn_result), intent(out) :: result !! Sorted nearest neighbors based on dmat.
        integer, intent(out), optional :: status !! Zero on success; nonzero for invalid input.
        integer :: i
        integer :: j
        integer :: m
        integer :: n
        real(dp), allocatable :: d(:)
        integer, allocatable :: id(:)

        if (present(status)) status = 0
        n = size(dmat, 1)
        if (size(dmat, 2) /= n .or. k < 1 .or. k >= n) then
            if (present(status)) status = 1
            allocate(result%id(0, 0), result%dist(0, 0))
            result%k = 0
            return
        end if
        result%k = k
        allocate(result%id(n, k), result%dist(n, k))
        allocate(d(n - 1), id(n - 1))
        do i = 1, n
            m = 0
            do j = 1, n
                if (j == i) cycle
                m = m + 1
                d(m) = dmat(i, j)
                id(m) = j
            end do
            call sort_pairs(d, id, m)
            result%id(i, :) = id(1:k)
            result%dist(i, :) = d(1:k)
        end do
    end subroutine knn_from_dist

    subroutine frnn(x, eps, result, status)
        real(dp), intent(in) :: x(:, :) !! Data matrix with observations in rows and features in columns.
        real(dp), intent(in) :: eps !! Nonnegative Euclidean neighborhood radius.
        type(frnn_result), intent(out) :: result !! Fixed-radius neighbors excluding self, sorted by distance and id.
        integer, intent(out), optional :: status !! Zero on success; nonzero for invalid input.
        real(dp), allocatable :: drow(:)
        integer, allocatable :: irow(:)
        integer, allocatable :: counts(:)
        integer :: i
        integer :: j
        integer :: m
        integer :: n
        integer :: pos
        real(dp) :: dij

        if (present(status)) status = 0
        n = size(x, 1)
        if (n < 1 .or. size(x, 2) < 1 .or. eps < 0.0_dp .or. contains_nan(x)) then
            if (present(status)) status = 1
            result%n = 0
            allocate(result%offset(1), result%id(0), result%dist(0))
            result%offset = 1
            return
        end if

        allocate(counts(n))
        counts = 0
        do i = 1, n
            do j = 1, n
                if (j == i) cycle
                dij = sqrt(sum((x(i, :) - x(j, :)) ** 2))
                if (dij <= eps) counts(i) = counts(i) + 1
            end do
        end do

        result%n = n
        result%eps = eps
        allocate(result%offset(n + 1))
        result%offset(1) = 1
        do i = 1, n
            result%offset(i + 1) = result%offset(i) + counts(i)
        end do
        allocate(result%id(sum(counts)), result%dist(sum(counts)))
        allocate(drow(max(1, n - 1)), irow(max(1, n - 1)))
        do i = 1, n
            m = 0
            do j = 1, n
                if (j == i) cycle
                dij = sqrt(sum((x(i, :) - x(j, :)) ** 2))
                if (dij <= eps) then
                    m = m + 1
                    drow(m) = dij
                    irow(m) = j
                end if
            end do
            if (m > 1) call sort_pairs(drow, irow, m)
            pos = result%offset(i)
            if (m > 0) then
                result%id(pos:pos + m - 1) = irow(1:m)
                result%dist(pos:pos + m - 1) = drow(1:m)
            end if
        end do
    end subroutine frnn

    subroutine frnn_query(x, query, eps, result, status)
        real(dp), intent(in) :: x(:, :) !! Reference data matrix with observations in rows.
        real(dp), intent(in) :: query(:, :) !! Query points with the same number of columns as x.
        real(dp), intent(in) :: eps !! Nonnegative Euclidean neighborhood radius.
        type(frnn_result), intent(out) :: result !! Fixed-radius reference neighbors for each query row.
        integer, intent(out), optional :: status !! Zero on success; nonzero for invalid input.
        integer, allocatable :: counts(:)
        real(dp), allocatable :: drow(:)
        integer, allocatable :: irow(:)
        integer :: i
        integer :: j
        integer :: m
        integer :: n
        integer :: nq
        integer :: pos
        real(dp) :: dij

        if (present(status)) status = 0
        n = size(x, 1)
        nq = size(query, 1)
        if (size(query, 2) /= size(x, 2) .or. eps < 0.0_dp) then
            if (present(status)) status = 1
            result%n = 0
            allocate(result%offset(1), result%id(0), result%dist(0))
            result%offset = 1
            return
        end if
        if (contains_nan(x) .or. contains_nan(query)) then
            if (present(status)) status = 2
            result%n = 0
            allocate(result%offset(1), result%id(0), result%dist(0))
            result%offset = 1
            return
        end if

        allocate(counts(nq))
        counts = 0
        do i = 1, nq
            do j = 1, n
                dij = sqrt(sum((query(i, :) - x(j, :)) ** 2))
                if (dij <= eps) counts(i) = counts(i) + 1
            end do
        end do
        result%n = nq
        result%eps = eps
        allocate(result%offset(nq + 1))
        result%offset(1) = 1
        do i = 1, nq
            result%offset(i + 1) = result%offset(i) + counts(i)
        end do
        allocate(result%id(sum(counts)), result%dist(sum(counts)))
        allocate(drow(max(1, n)), irow(max(1, n)))
        do i = 1, nq
            m = 0
            do j = 1, n
                dij = sqrt(sum((query(i, :) - x(j, :)) ** 2))
                if (dij <= eps) then
                    m = m + 1
                    drow(m) = dij
                    irow(m) = j
                end if
            end do
            if (m > 1) call sort_pairs(drow, irow, m)
            pos = result%offset(i)
            if (m > 0) then
                result%id(pos:pos + m - 1) = irow(1:m)
                result%dist(pos:pos + m - 1) = drow(1:m)
            end if
        end do
    end subroutine frnn_query

    subroutine frnn_from_dist(dmat, eps, result, status)
        real(dp), intent(in) :: dmat(:, :) !! Full symmetric distance matrix.
        real(dp), intent(in) :: eps !! Nonnegative neighborhood radius.
        type(frnn_result), intent(out) :: result !! Fixed-radius neighbors excluding self.
        integer, intent(out), optional :: status !! Zero on success; nonzero for invalid input.
        integer, allocatable :: counts(:)
        real(dp), allocatable :: drow(:)
        integer, allocatable :: irow(:)
        integer :: i
        integer :: j
        integer :: m
        integer :: n
        integer :: pos

        if (present(status)) status = 0
        n = size(dmat, 1)
        if (size(dmat, 2) /= n .or. eps < 0.0_dp) then
            if (present(status)) status = 1
            result%n = 0
            allocate(result%offset(1), result%id(0), result%dist(0))
            result%offset = 1
            return
        end if
        allocate(counts(n))
        counts = 0
        do i = 1, n
            do j = 1, n
                if (j /= i .and. dmat(i, j) <= eps) counts(i) = counts(i) + 1
            end do
        end do
        result%n = n
        result%eps = eps
        allocate(result%offset(n + 1))
        result%offset(1) = 1
        do i = 1, n
            result%offset(i + 1) = result%offset(i) + counts(i)
        end do
        allocate(result%id(sum(counts)), result%dist(sum(counts)))
        allocate(drow(max(1, n - 1)), irow(max(1, n - 1)))
        do i = 1, n
            m = 0
            do j = 1, n
                if (j == i .or. dmat(i, j) > eps) cycle
                m = m + 1
                drow(m) = dmat(i, j)
                irow(m) = j
            end do
            if (m > 1) call sort_pairs(drow, irow, m)
            pos = result%offset(i)
            if (m > 0) then
                result%id(pos:pos + m - 1) = irow(1:m)
                result%dist(pos:pos + m - 1) = drow(1:m)
            end if
        end do
    end subroutine frnn_from_dist

    function adjacencylist_knn(nn) result(adj)
        type(knn_result), intent(in) :: nn !! k-nearest-neighbor object to convert.
        type(frnn_result) :: adj
        integer :: i
        integer :: k
        integer :: n
        integer :: p

        n = size(nn%id, 1)
        k = size(nn%id, 2)
        adj%n = n
        adj%eps = huge(1.0_dp)
        allocate(adj%offset(n + 1), adj%id(n * k), adj%dist(n * k))
        adj%offset(1) = 1
        do i = 1, n
            adj%offset(i + 1) = adj%offset(i) + k
            p = adj%offset(i)
            adj%id(p:p + k - 1) = nn%id(i, :)
            adj%dist(p:p + k - 1) = nn%dist(i, :)
        end do
    end function adjacencylist_knn

    function knn_dist(x, k, status) result(distance)
        real(dp), intent(in) :: x(:, :) !! Data matrix with observations in rows.
        integer, intent(in) :: k !! Neighbor rank whose distance is returned.
        integer, intent(out), optional :: status !! Zero on success; nonzero for invalid input.
        real(dp), allocatable :: distance(:)
        type(knn_result) :: nn
        integer :: stat

        call knn(x, k, nn, stat)
        if (present(status)) status = stat
        if (stat /= 0) then
            allocate(distance(0))
            return
        end if
        allocate(distance(size(x, 1)))
        distance = nn%dist(:, k)
    end function knn_dist

end module dbscan_nn
