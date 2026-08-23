module rnanoflann_search
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use rnanoflann_kinds, only: dp
    use rnanoflann_types, only: nn_result
    use rnanoflann_metrics, only: metric_code, metric_hellinger, metric_minkowski, metric_distance_prepared
!$  use omp_lib, only: omp_set_num_threads
    implicit none
    private
    public :: nn

contains

    function nn(data, points, k, method, search, eps, square, sorted, radius, trans, leafs, p, parallel, cores) result(res)
        real(dp), intent(in) :: data(:,:), points(:,:)
        integer, intent(in), optional :: k
        character(*), intent(in), optional :: method, search
        real(dp), intent(in), optional :: eps, radius, p
        logical, intent(in), optional :: square, sorted, trans, parallel
        integer, intent(in), optional :: leafs, cores
        type(nn_result) :: res

        integer :: kk, code, ndata, nq, ndim, leaf_count, ncores
        real(dp) :: eps_value, radius_value, power
        logical :: use_square, use_sorted, use_trans, use_parallel
        character(len=:), allocatable :: method_name, search_name
        integer, allocatable :: idx_work(:,:)
        real(dp), allocatable :: dist_work(:,:)
        real(dp), allocatable :: data_h(:,:), points_h(:,:)

        ndata = size(data, 1)
        ndim = size(data, 2)
        nq = size(points, 1)
        if (size(points, 2) /= ndim) error stop "nn: data and points must have the same number of columns"
        if (ndata < 1 .or. nq < 1 .or. ndim < 1) error stop "nn: matrices must be nonempty"

        kk = ndata
        if (present(k)) kk = k
        if (kk < 1 .or. kk > ndata) error stop "nn: k must satisfy 1 <= k <= size(data,1)"

        method_name = "euclidean"
        if (present(method)) method_name = trim(adjustl(method))
        code = metric_code(method_name)
        if (code == 0) error stop "nn: unsupported distance method"

        search_name = "standard"
        if (present(search)) search_name = lower_ascii(trim(adjustl(search)))
        if (search_name /= "standard" .and. search_name /= "radius") then
            error stop "nn: search must be 'standard' or 'radius'"
        end if

        use_square = .false.
        if (present(square)) use_square = square
        use_sorted = .false.
        if (present(sorted)) use_sorted = sorted
        use_trans = .true.
        if (present(trans)) use_trans = trans
        use_parallel = .false.
        if (present(parallel)) use_parallel = parallel

        eps_value = 0.0_dp
        if (present(eps)) eps_value = eps
        if (eps_value < 0.0_dp) error stop "nn: eps must be nonnegative"

        radius_value = 0.0_dp
        if (present(radius)) radius_value = radius
        if (search_name == "radius" .and. radius_value < 0.0_dp) error stop "nn: radius must be nonnegative"

        power = 0.0_dp
        if (present(p)) power = p
        if (code == metric_minkowski .and. power <= 0.0_dp) error stop "nn: Minkowski p must be positive"

        leaf_count = 10
        if (present(leafs)) leaf_count = leafs
        if (leaf_count < 1) error stop "nn: leafs must be positive"
        ncores = 0
        if (present(cores)) ncores = cores
        if (ncores < 0) error stop "nn: cores must be nonnegative"

        allocate(idx_work(kk, nq), dist_work(kk, nq), res%counts(nq))
        idx_work = 0
        dist_work = sqrt(huge(0.0_dp))
        res%counts = 0

        if (code == metric_hellinger) then
            if (any(data < 0.0_dp) .or. any(points < 0.0_dp)) then
                error stop "nn: Hellinger inputs must be nonnegative"
            end if
            allocate(data_h(ndata, ndim), points_h(nq, ndim))
            data_h = sqrt(data)
            points_h = sqrt(points)
            call search_queries(data_h, points_h, kk, code, search_name, use_square, use_sorted, &
                radius_value, power, use_parallel, ncores, idx_work, dist_work, res%counts)
        else
            call search_queries(data, points, kk, code, search_name, use_square, use_sorted, &
                radius_value, power, use_parallel, ncores, idx_work, dist_work, res%counts)
        end if

        ! eps and leafs are retained for source/API compatibility. In the upstream
        ! implementation SearchParameters is not passed to either search call, so
        ! eps does not affect results. The custom metric adaptors also return zero
        ! from accum_dist, preventing distance-bound pruning.

        if (use_trans) then
            allocate(res%indices(nq, kk), res%distances(nq, kk))
            res%indices = transpose(idx_work)
            res%distances = transpose(dist_work)
        else
            allocate(res%indices(kk, nq), res%distances(kk, nq))
            res%indices = idx_work
            res%distances = dist_work
        end if
    end function nn

    subroutine search_queries(data, points, k, code, search_name, square, sorted, radius, p, &
        parallel, cores, indices, distances, counts)
        real(dp), intent(in) :: data(:,:), points(:,:)
        integer, intent(in) :: k, code, cores
        character(*), intent(in) :: search_name
        logical, intent(in) :: square, sorted, parallel
        real(dp), intent(in) :: radius, p
        integer, intent(out) :: indices(:,:), counts(:)
        real(dp), intent(out) :: distances(:,:)
        integer :: iq

        if (parallel) then
!$          if (cores > 0) call omp_set_num_threads(cores)
!$omp parallel do default(shared) private(iq) schedule(static)
            do iq = 1, size(points, 1)
                call search_one(data, points(iq,:), k, code, search_name, square, sorted, radius, p, &
                    indices(:,iq), distances(:,iq), counts(iq))
            end do
!$omp end parallel do
        else
            do iq = 1, size(points, 1)
                call search_one(data, points(iq,:), k, code, search_name, square, sorted, radius, p, &
                    indices(:,iq), distances(:,iq), counts(iq))
            end do
        end if
        if (cores < 0) error stop "search_queries: unreachable core guard"
    end subroutine search_queries

    subroutine search_one(data, query, k, code, search_name, square, sorted, radius, p, indices, distances, count)
        real(dp), intent(in) :: data(:,:), query(:)
        integer, intent(in) :: k, code
        character(*), intent(in) :: search_name
        logical, intent(in) :: square, sorted
        real(dp), intent(in) :: radius, p
        integer, intent(out) :: indices(:), count
        real(dp), intent(out) :: distances(:)

        integer :: i, nfound, take
        real(dp) :: d
        integer, allocatable :: ri(:)
        real(dp), allocatable :: rd(:)

        indices = 0
        distances = sqrt(huge(0.0_dp))
        count = 0

        if (search_name == "standard") then
            do i = 1, size(data, 1)
                d = metric_distance_prepared(data(i,:), query, code, square, p)
                if (ieee_is_finite(d)) call insert_best(indices, distances, i, d)
            end do
            count = 0
            do i = 1, size(indices)
                if (indices(i) > 0) count = count + 1
            end do
            return
        end if

        allocate(ri(size(data, 1)), rd(size(data, 1)))
        nfound = 0
        do i = 1, size(data, 1)
            d = metric_distance_prepared(data(i,:), query, code, square, p)
            if (ieee_is_finite(d) .and. d < radius) then
                nfound = nfound + 1
                ri(nfound) = i
                rd(nfound) = d
            end if
        end do
        if (sorted .and. nfound > 1) call sort_pairs(ri, rd, nfound)
        take = min(k, nfound)
        if (take > 0) then
            indices(1:take) = ri(1:take)
            distances(1:take) = rd(1:take)
        end if
        count = nfound
    end subroutine search_one

    subroutine insert_best(indices, distances, index_value, distance_value)
        integer, intent(inout) :: indices(:)
        real(dp), intent(inout) :: distances(:)
        integer, intent(in) :: index_value
        real(dp), intent(in) :: distance_value
        integer :: j, pos

        pos = size(distances) + 1
        do j = 1, size(distances)
            if (distance_value < distances(j)) then
                pos = j
                exit
            else if (distance_value <= distances(j) .and. index_value < indices(j)) then
                pos = j
                exit
            end if
        end do
        if (pos > size(distances)) return
        do j = size(distances), pos + 1, -1
            distances(j) = distances(j - 1)
            indices(j) = indices(j - 1)
        end do
        distances(pos) = distance_value
        indices(pos) = index_value
    end subroutine insert_best

    subroutine sort_pairs(indices, distances, n)
        integer, intent(inout) :: indices(:)
        real(dp), intent(inout) :: distances(:)
        integer, intent(in) :: n
        integer :: i, j, ii
        real(dp) :: dd

        do i = 2, n
            ii = indices(i)
            dd = distances(i)
            j = i - 1
            do while (j >= 1)
                if (distances(j) < dd) exit
                if (distances(j) <= dd .and. indices(j) <= ii) exit
                distances(j + 1) = distances(j)
                indices(j + 1) = indices(j)
                j = j - 1
            end do
            distances(j + 1) = dd
            indices(j + 1) = ii
        end do
    end subroutine sort_pairs

    pure function lower_ascii(text) result(out)
        character(*), intent(in) :: text
        character(len=len(text)) :: out
        integer :: i, c
        do i = 1, len(text)
            c = iachar(text(i:i))
            if (c >= iachar('A') .and. c <= iachar('Z')) then
                out(i:i) = achar(c + 32)
            else
                out(i:i) = text(i:i)
            end if
        end do
    end function lower_ascii

end module rnanoflann_search
