! SPDX-License-Identifier: GPL-2.0-only
module slam_utils
    use iso_fortran_env, only : int64
    implicit none
    private

    public :: linear_index
    public :: linear_indices
    public :: coordinates_from_linear
    public :: argsort_int64
    public :: product_dims

contains

    pure function product_dims(dims) result(n)
        integer, intent(in) :: dims(:)
        integer(int64) :: n
        integer :: k

        n = 1_int64
        do k = 1, size(dims)
            if (dims(k) < 0) then
                n = -1_int64
                return
            end if
            n = n * int(dims(k), int64)
        end do
    end function product_dims

    pure function linear_index(coords, dims) result(idx)
        integer, intent(in) :: coords(:)
        integer, intent(in) :: dims(:)
        integer(int64) :: idx
        integer(int64) :: stride
        integer :: k

        if (size(coords) /= size(dims)) then
            idx = 0_int64
            return
        end if
        if (size(dims) == 0) then
            idx = 0_int64
            return
        end if

        idx = int(coords(1), int64)
        if (coords(1) < 1 .or. coords(1) > dims(1)) then
            idx = 0_int64
            return
        end if
        stride = int(dims(1), int64)
        do k = 2, size(dims)
            if (coords(k) < 1 .or. coords(k) > dims(k)) then
                idx = 0_int64
                return
            end if
            idx = idx + stride * int(coords(k) - 1, int64)
            stride = stride * int(dims(k), int64)
        end do
    end function linear_index

    pure function linear_indices(coords, dims) result(idx)
        integer, intent(in) :: coords(:,:)
        integer, intent(in) :: dims(:)
        integer(int64), allocatable :: idx(:)
        integer :: r

        allocate(idx(size(coords, 1)))
        if (size(coords, 2) /= size(dims)) then
            idx = 0_int64
            return
        end if
        do r = 1, size(coords, 1)
            idx(r) = linear_index(coords(r, :), dims)
        end do
    end function linear_indices

    pure subroutine coordinates_from_linear(idx, dims, coords, ok)
        integer(int64), intent(in) :: idx
        integer, intent(in) :: dims(:)
        integer, intent(out) :: coords(:)
        logical, intent(out), optional :: ok
        integer(int64) :: z, n
        integer :: k
        logical :: good

        good = size(coords) == size(dims)
        n = product_dims(dims)
        good = good .and. n >= 0_int64 .and. idx >= 1_int64 .and. idx <= n
        if (.not. good) then
            if (size(coords) > 0) coords = 0
            if (present(ok)) ok = .false.
            return
        end if

        z = idx - 1_int64
        do k = 1, size(dims)
            coords(k) = int(mod(z, int(dims(k), int64))) + 1
            z = z / int(dims(k), int64)
        end do
        if (present(ok)) ok = .true.
    end subroutine coordinates_from_linear

    subroutine argsort_int64(keys, order)
        integer(int64), intent(in) :: keys(:)
        integer, allocatable, intent(out) :: order(:)
        integer(int64), allocatable :: work(:)
        integer :: k

        allocate(order(size(keys)), work(size(keys)))
        do k = 1, size(keys)
            order(k) = k
        end do
        work = keys
        if (size(keys) > 1) call quicksort(work, order, 1, size(keys))
    contains
        recursive subroutine quicksort(a, p, left, right)
            integer(int64), intent(inout) :: a(:)
            integer, intent(inout) :: p(:)
            integer, intent(in) :: left, right
            integer :: i, j, tp
            integer(int64) :: pivot, tv

            i = left
            j = right
            pivot = a((left + right) / 2)
            do
                do while (a(i) < pivot)
                    i = i + 1
                end do
                do while (a(j) > pivot)
                    j = j - 1
                end do
                if (i <= j) then
                    tv = a(i); a(i) = a(j); a(j) = tv
                    tp = p(i); p(i) = p(j); p(j) = tp
                    i = i + 1
                    j = j - 1
                end if
                if (i > j) exit
            end do
            if (left < j) call quicksort(a, p, left, j)
            if (i < right) call quicksort(a, p, i, right)
        end subroutine quicksort
    end subroutine argsort_int64

end module slam_utils
