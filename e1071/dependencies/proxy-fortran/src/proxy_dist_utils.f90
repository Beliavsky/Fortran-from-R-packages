module proxy_dist_utils
    use, intrinsic :: ieee_arithmetic, only: ieee_is_nan, ieee_is_finite
    use proxy_kinds, only: dp
    use proxy_ieee, only: proxy_nan
    implicit none
    private

    public :: proxy_pack_dist, proxy_unpack_dist, proxy_subset_dist
    public :: proxy_row_sums_dist, proxy_row_means_dist, proxy_row_dist, proxy_col_dist

contains

    subroutine proxy_pack_dist(matrix, packed)
        real(dp), intent(in) :: matrix(:, :) !! Square symmetric proximity matrix; only the strict lower triangle is stored in R
        !! `dist` column order.
        real(dp), allocatable, intent(out) :: packed(:) !! Allocated packed strict-lower-triangle vector of length `n*(n-1)/2`.
        integer :: i
        integer :: j
        integer :: k
        integer :: n

        n = size(matrix, 1)
        if (size(matrix, 2) /= n) then
            allocate(packed(0))
            return
        end if
        allocate(packed(n * (n - 1) / 2))
        k = 0
        do j = 1, n - 1
            do i = j + 1, n
                k = k + 1
                packed(k) = matrix(i, j)
            end do
        end do
    end subroutine proxy_pack_dist

    subroutine proxy_unpack_dist(packed, n, matrix, diagonal)
        real(dp), intent(in) :: packed(:) !! Packed strict-lower-triangle values in R `dist` column order.
        integer, intent(in) :: n !! Intended square-matrix dimension; requires `size(packed)=n*(n-1)/2`.
        real(dp), allocatable, intent(out) :: matrix(:, :) !! Allocated symmetric `n` by `n` matrix, or empty on an inconsistent
        !! packed length.
        real(dp), intent(in), optional :: diagonal !! Optional diagonal value; defaults to zero for distance-like matrices.
        integer :: i
        integer :: j
        integer :: k
        real(dp) :: diag_value

        if (size(packed) /= n * (n - 1) / 2) then
            allocate(matrix(0, 0))
            return
        end if
        diag_value = 0.0_dp
        if (present(diagonal)) diag_value = diagonal
        allocate(matrix(n, n))
        matrix = 0.0_dp
        do i = 1, n
            matrix(i, i) = diag_value
        end do
        k = 0
        do j = 1, n - 1
            do i = j + 1, n
                k = k + 1
                matrix(i, j) = packed(k)
                matrix(j, i) = packed(k)
            end do
        end do
    end subroutine proxy_unpack_dist

    subroutine proxy_subset_dist(packed, n, indices, subset)
        real(dp), intent(in) :: packed(:) !! Original packed `dist` vector using R lower-triangle storage order.
        integer, intent(in) :: n !! Size of the original symmetric proximity matrix.
        integer, intent(in) :: indices(:) !! One-based object indices to retain, including duplicates; duplicates produce NaN
        !! pair entries like `subset.dist`.
        real(dp), allocatable, intent(out) :: subset(:) !! Allocated packed vector for the selected objects in the order given by
        !! `indices`.
        integer :: i
        integer :: j
        integer :: k
        integer :: si
        integer :: sj

        if (size(packed) /= n * (n - 1) / 2 .or. any(indices < 1) .or. any(indices > n)) then
            allocate(subset(0))
            return
        end if
        allocate(subset(size(indices) * (size(indices) - 1) / 2))
        k = 0
        do j = 1, size(indices) - 1
            sj = indices(j)
            do i = j + 1, size(indices)
                si = indices(i)
                k = k + 1
                if (si == sj) then
                    subset(k) = proxy_nan()
                else
                    subset(k) = packed(packed_index(max(si, sj), min(si, sj), n))
                end if
            end do
        end do
    end subroutine proxy_subset_dist

    subroutine proxy_row_sums_dist(packed, n, sums, na_rm)
        real(dp), intent(in) :: packed(:) !! Packed symmetric proximity values in R `dist` storage order.
        integer, intent(in) :: n !! Number of objects represented by `packed`.
        real(dp), allocatable, intent(out) :: sums(:) !! Allocated row/column sums excluding the implicit diagonal.
        logical, intent(in), optional :: na_rm !! If true, omit NaN values; default false propagates NaN to each affected row as
        !! in proxy.
        logical :: remove_na
        integer :: i
        integer :: j
        integer :: k
        real(dp) :: value

        remove_na = .false.
        if (present(na_rm)) remove_na = na_rm
        if (size(packed) /= n * (n - 1) / 2) then
            allocate(sums(0))
            return
        end if
        allocate(sums(n))
        sums = 0.0_dp
        k = 0
        do j = 1, n - 1
            do i = j + 1, n
                k = k + 1
                value = packed(k)
                if (.not. ieee_is_finite(value)) then
                    if (ieee_is_nan(value)) then
                        if (remove_na) cycle
                        sums(i) = proxy_nan()
                        sums(j) = proxy_nan()
                    else
                        sums(i) = value
                        sums(j) = value
                    end if
                    exit
                end if
                if (.not. ieee_is_nan(sums(i))) sums(i) = sums(i) + value
                if (.not. ieee_is_nan(sums(j))) sums(j) = sums(j) + value
            end do
        end do
    end subroutine proxy_row_sums_dist

    subroutine proxy_row_means_dist(packed, n, means, na_rm, diag)
        real(dp), intent(in) :: packed(:) !! Packed symmetric proximity values in R `dist` order.
        integer, intent(in) :: n !! Number of represented objects.
        real(dp), allocatable, intent(out) :: means(:) !! Allocated row means with denominator semantics matching `rowMeans.dist`.
        logical, intent(in), optional :: na_rm !! If true, divide by each row's count of nonmissing off-diagonal values plus
        !! optional diagonal; default false.
        logical, intent(in), optional :: diag !! If true, include an implicit diagonal element in the denominator; default true,
        !! matching proxy.
        real(dp), allocatable :: sums(:)
        integer, allocatable :: counts(:)
        logical :: include_diag
        logical :: remove_na
        integer :: i
        integer :: j
        integer :: k

        include_diag = .true.
        remove_na = .false.
        if (present(diag)) include_diag = diag
        if (present(na_rm)) remove_na = na_rm
        call proxy_row_sums_dist(packed, n, sums, remove_na)
        if (size(sums) == 0) then
            allocate(means(0))
            return
        end if
        allocate(means(n))
        if (.not. remove_na) then
            means = sums / real(n - merge(0, 1, include_diag), dp)
            return
        end if
        allocate(counts(n))
        counts = merge(1, 0, include_diag)
        k = 0
        do j = 1, n - 1
            do i = j + 1, n
                k = k + 1
                if (.not. ieee_is_nan(packed(k))) then
                    counts(i) = counts(i) + 1
                    counts(j) = counts(j) + 1
                end if
            end do
        end do
        do i = 1, n
            if (counts(i) > 0) then
                means(i) = sums(i) / real(counts(i), dp)
            else
                means(i) = proxy_nan()
            end if
        end do
    end subroutine proxy_row_means_dist

    subroutine proxy_row_dist(n, indices)
        integer, intent(in) :: n !! Size of the symmetric matrix represented by an R-style packed `dist` vector.
        integer, allocatable, intent(out) :: indices(:) !! Row indices corresponding to each packed lower-triangle value,
        !! matching `row.dist`.
        integer :: i
        integer :: j
        integer :: k

        allocate(indices(n * (n - 1) / 2))
        k = 0
        do j = 1, n - 1
            do i = j + 1, n
                k = k + 1
                indices(k) = i
            end do
        end do
    end subroutine proxy_row_dist

    subroutine proxy_col_dist(n, indices)
        integer, intent(in) :: n !! Size of the symmetric matrix represented by an R-style packed `dist` vector.
        integer, allocatable, intent(out) :: indices(:) !! Column indices corresponding to each packed lower-triangle value,
        !! matching `col.dist`.
        integer :: i
        integer :: j
        integer :: k

        allocate(indices(n * (n - 1) / 2))
        k = 0
        do j = 1, n - 1
            do i = j + 1, n
                k = k + 1
                indices(k) = j
            end do
        end do
    end subroutine proxy_col_dist

    pure function packed_index(i, j, n) result(index)
        integer, intent(in) :: i !! One-based row index with `i>j` in the symmetric matrix.
        integer, intent(in) :: j !! One-based column index with `j<i`.
        integer, intent(in) :: n !! Matrix size used by the R `dist` packed-index formula.
        integer :: index

        index = i + (j - 1) * (n - 1) - j * (j - 1) / 2 - 1
    end function packed_index

end module proxy_dist_utils
