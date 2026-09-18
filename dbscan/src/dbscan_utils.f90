module dbscan_utils
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
    use dbscan_kinds, only : dp
    implicit none
    private

    public :: contains_nan, sort_pairs, compact_labels, unique_positive_count
    public :: count_zero, median_value

contains

    pure logical function contains_nan(x) result(found)
        real(dp), intent(in) :: x(:, :) !! Numeric matrix to inspect for NaN values.
        integer :: i
        integer :: j

        found = .false.
        do j = 1, size(x, 2)
            do i = 1, size(x, 1)
                if (ieee_is_nan(x(i, j))) then
                    found = .true.
                    return
                end if
            end do
        end do
    end function contains_nan

    pure subroutine sort_pairs(dist, id, n, decreasing)
        real(dp), intent(inout) :: dist(:) !! Distances to sort in place.
        integer, intent(inout) :: id(:) !! Integer identifiers permuted with distances.
        integer, intent(in) :: n !! Number of active entries in dist and id.
        logical, intent(in), optional :: decreasing !! Reverse ordering when true.
        integer :: i
        integer :: j
        integer :: key_id
        real(dp) :: key_dist
        logical :: dec
        logical :: move

        dec = .false.
        if (present(decreasing)) dec = decreasing

        do i = 2, n
            key_dist = dist(i)
            key_id = id(i)
            j = i - 1
            do while (j >= 1)
                if (dec) then
                    if (dist(j) < key_dist) then
                        move = .true.
                    else if (dist(j) > key_dist) then
                        move = .false.
                    else
                        move = id(j) < key_id
                    end if
                else
                    if (dist(j) > key_dist) then
                        move = .true.
                    else if (dist(j) < key_dist) then
                        move = .false.
                    else
                        move = id(j) > key_id
                    end if
                end if
                if (.not. move) exit
                dist(j + 1) = dist(j)
                id(j + 1) = id(j)
                j = j - 1
            end do
            dist(j + 1) = key_dist
            id(j + 1) = key_id
        end do
    end subroutine sort_pairs

    pure subroutine compact_labels(labels)
        integer, intent(inout) :: labels(:) !! Cluster labels; positive labels are compacted in first-seen order.
        integer, allocatable :: old(:)
        integer :: i
        integer :: j
        integer :: nold
        logical :: found

        allocate(old(size(labels)))
        nold = 0
        do i = 1, size(labels)
            if (labels(i) <= 0) cycle
            found = .false.
            do j = 1, nold
                if (old(j) == labels(i)) then
                    labels(i) = j
                    found = .true.
                    exit
                end if
            end do
            if (.not. found) then
                nold = nold + 1
                old(nold) = labels(i)
                labels(i) = nold
            end if
        end do
    end subroutine compact_labels

    pure integer function unique_positive_count(labels) result(n)
        integer, intent(in) :: labels(:) !! Cluster labels where zero denotes noise.
        integer, allocatable :: seen(:)
        integer :: i

        if (size(labels) == 0) then
            n = 0
            return
        end if
        allocate(seen(max(1, maxval(labels))))
        seen = 0
        do i = 1, size(labels)
            if (labels(i) > 0) seen(labels(i)) = 1
        end do
        n = sum(seen)
    end function unique_positive_count

    pure integer function count_zero(labels) result(n)
        integer, intent(in) :: labels(:) !! Cluster labels where zero denotes noise.
        n = count(labels == 0)
    end function count_zero

    pure real(dp) function median_value(x) result(med)
        real(dp), intent(in) :: x(:) !! Values whose median is requested.
        real(dp), allocatable :: y(:)
        integer, allocatable :: ids(:)
        integer :: n
        integer :: i

        n = size(x)
        if (n == 0) then
            med = 0.0_dp
            return
        end if
        allocate(y(n), ids(n))
        y = x
        do i = 1, n
            ids(i) = i
        end do
        call sort_pairs(y, ids, n)
        if (mod(n, 2) == 1) then
            med = y((n + 1) / 2)
        else
            med = 0.5_dp * (y(n / 2) + y(n / 2 + 1))
        end if
    end function median_value

end module dbscan_utils
