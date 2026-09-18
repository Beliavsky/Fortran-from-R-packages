module dbscan_dbcv
    use dbscan_kinds, only : dp
    use dbscan_types, only : dbcv_result
    use dbscan_distance, only : pairwise_distance_matrix
    use dbscan_hdbscan, only : mst_dense
    implicit none
    private

    public :: dbcv

contains

    subroutine dbcv(x, cluster, result, metric, status)
        real(dp), intent(in) :: x(:, :) !! Data matrix with observations in rows.
        integer, intent(in) :: cluster(:) !! Cluster labels; zero denotes noise.
        type(dbcv_result), intent(out) :: result !! DBCV overall and per-cluster validity measures.
        character(len=*), intent(in), optional :: metric !! euclidean, sqeuclidean, or manhattan.
        integer, intent(out), optional :: status !! Zero on success; nonzero for invalid/undefined clustering.
        character(len=:), allocatable :: met
        real(dp), allocatable :: dmat(:, :)
        real(dp), allocatable :: core(:)
        real(dp), allocatable :: mrd(:, :)
        real(dp), allocatable :: sub(:, :)
        real(dp), allocatable :: mst(:, :)
        integer, allocatable :: ids(:)
        integer, allocatable :: valid_label(:)
        integer, allocatable :: count_label(:)
        integer, allocatable :: point_cluster(:)
        integer, allocatable :: internal(:)
        integer, allocatable :: degree(:)
        integer, allocatable :: internal_count(:)
        integer :: n
        integer :: d
        integer :: ncl
        integer :: i
        integer :: j
        integer :: c
        integer :: ci
        integer :: cj
        integer :: p
        integer :: q
        integer :: m
        integer :: stat
        integer :: edge
        integer :: u
        integer :: v
        integer :: ni
        integer :: nj
        real(dp) :: sum_inv
        real(dp) :: dij
        real(dp) :: min_sep
        real(dp) :: best

        if (present(status)) status = 0
        n = size(x, 1)
        d = size(x, 2)
        if (n < 1 .or. d < 1 .or. size(cluster) /= n) then
            if (present(status)) status = 1
            return
        end if
        met = 'euclidean'
        if (present(metric)) met = trim(adjustl(metric))
        call pairwise_distance_matrix(x, dmat, met, stat)
        if (stat /= 0) then
            if (present(status)) status = stat
            return
        end if

        allocate(valid_label(n), count_label(n))
        ncl = 0
        count_label = 0
        do i = 1, n
            if (cluster(i) <= 0) cycle
            c = 0
            do j = 1, ncl
                if (valid_label(j) == cluster(i)) then
                    c = j
                    exit
                end if
            end do
            if (c == 0) then
                ncl = ncl + 1
                valid_label(ncl) = cluster(i)
                c = ncl
            end if
            count_label(c) = count_label(c) + 1
        end do

        c = 0
        do i = 1, ncl
            if (count_label(i) >= 3) then
                c = c + 1
                valid_label(c) = valid_label(i)
                count_label(c) = count_label(i)
            end if
        end do
        ncl = c
        if (ncl <= 1) then
            result%score = -1.0_dp
            result%n = n
            result%d = d
            allocate(result%cluster_size(max(0, ncl)), result%dsc(max(0, ncl)))
            allocate(result%dspc(max(0, ncl), max(0, ncl)), result%validity(max(0, ncl)))
            if (present(status)) status = 2
            return
        end if

        allocate(point_cluster(n), result%cluster_size(ncl), internal_count(ncl))
        point_cluster = 0
        do i = 1, n
            do c = 1, ncl
                if (cluster(i) == valid_label(c)) then
                    point_cluster(i) = c
                    exit
                end if
            end do
        end do
        result%cluster_size = count_label(1:ncl)
        result%n = n
        result%d = d
        allocate(core(n), mrd(n, n), result%dsc(ncl), result%dspc(ncl, ncl), result%validity(ncl))
        core = 0.0_dp
        mrd = 0.0_dp
        result%dsc = 0.0_dp
        result%dspc = huge(1.0_dp)
        internal_count = 0
        allocate(internal(n))
        internal = 0

        do c = 1, ncl
            m = result%cluster_size(c)
            allocate(ids(m))
            p = 0
            do i = 1, n
                if (point_cluster(i) == c) then
                    p = p + 1
                    ids(p) = i
                end if
            end do
            do p = 1, m
                sum_inv = 0.0_dp
                do q = 1, m
                    if (q == p) cycle
                    dij = dmat(ids(p), ids(q))
                    if (dij <= tiny(1.0_dp)) then
                        sum_inv = huge(1.0_dp)
                        exit
                    end if
                    sum_inv = sum_inv + (1.0_dp / dij) ** d
                end do
                if (sum_inv >= huge(1.0_dp) / 2.0_dp) then
                    core(ids(p)) = 0.0_dp
                else
                    core(ids(p)) = (sum_inv / real(m - 1, dp)) ** (-1.0_dp / real(d, dp))
                end if
            end do
            deallocate(ids)
        end do

        do i = 1, n
            if (point_cluster(i) == 0) cycle
            do j = i + 1, n
                if (point_cluster(j) == 0) cycle
                mrd(i, j) = max(dmat(i, j), max(core(i), core(j)))
                mrd(j, i) = mrd(i, j)
            end do
        end do

        do c = 1, ncl
            m = result%cluster_size(c)
            allocate(ids(m), sub(m, m), degree(m))
            p = 0
            do i = 1, n
                if (point_cluster(i) == c) then
                    p = p + 1
                    ids(p) = i
                end if
            end do
            do p = 1, m
                do q = 1, m
                    sub(p, q) = mrd(ids(p), ids(q))
                end do
            end do
            mst = mst_dense(sub, stat)
            degree = 0
            do edge = 1, size(mst, 1)
                u = nint(mst(edge, 1))
                v = nint(mst(edge, 2))
                degree(u) = degree(u) + 1
                degree(v) = degree(v) + 1
            end do
            result%dsc(c) = -huge(1.0_dp)
            do edge = 1, size(mst, 1)
                u = nint(mst(edge, 1))
                v = nint(mst(edge, 2))
                if (degree(u) > 1 .and. degree(v) > 1) result%dsc(c) = max(result%dsc(c), mst(edge, 3))
            end do
            if (result%dsc(c) <= -huge(1.0_dp) / 2.0_dp) result%dsc(c) = maxval(mst(:, 3))
            do p = 1, m
                if (degree(p) > 1) then
                    internal_count(c) = internal_count(c) + 1
                    internal(ids(p)) = 1
                end if
            end do
            deallocate(ids, sub, degree)
            if (allocated(mst)) deallocate(mst)
        end do

        do ci = 1, ncl
            result%dspc(ci, ci) = 0.0_dp
            do cj = ci + 1, ncl
                best = huge(1.0_dp)
                ni = internal_count(ci)
                nj = internal_count(cj)
                if ((ni > 1 .or. nj > 1) .and. ni > 0 .and. nj > 0) then
                    do i = 1, n
                        if (point_cluster(i) /= ci .or. internal(i) == 0) cycle
                        do j = 1, n
                            if (point_cluster(j) /= cj .or. internal(j) == 0) cycle
                            best = min(best, mrd(i, j))
                        end do
                    end do
                end if
                result%dspc(ci, cj) = best
                result%dspc(cj, ci) = best
            end do
        end do

        result%score = 0.0_dp
        do ci = 1, ncl
            min_sep = huge(1.0_dp)
            do cj = 1, ncl
                if (ci == cj) cycle
                min_sep = min(min_sep, result%dspc(ci, cj))
            end do
            result%validity(ci) = (min_sep - result%dsc(ci)) / max(min_sep, result%dsc(ci))
            result%score = result%score + real(result%cluster_size(ci), dp) / real(n, dp) * result%validity(ci)
        end do
        if (present(status)) status = 0
    end subroutine dbcv

end module dbscan_dbcv
