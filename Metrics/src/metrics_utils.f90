! SPDX-License-Identifier: BSD-3-Clause
module metrics_utils
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan, ieee_positive_inf, ieee_negative_inf
    use metrics_kinds, only : dp
    implicit none
    private

    public :: quiet_nan, positive_infinity, negative_infinity
    public :: mean_value, median_value, average_ranks

contains

    pure real(dp) function quiet_nan() result(x)
        x = ieee_value(0.0_dp, ieee_quiet_nan)
    end function quiet_nan

    pure real(dp) function positive_infinity() result(x)
        x = ieee_value(0.0_dp, ieee_positive_inf)
    end function positive_infinity

    pure real(dp) function negative_infinity() result(x)
        x = ieee_value(0.0_dp, ieee_negative_inf)
    end function negative_infinity

    pure real(dp) function mean_value(x) result(value)
        real(dp), intent(in) :: x(:)

        if (size(x) == 0) then
            value = quiet_nan()
        else
            value = sum(x) / real(size(x), dp)
        end if
    end function mean_value

    pure real(dp) function median_value(x) result(value)
        real(dp), intent(in) :: x(:)
        real(dp), allocatable :: work(:)
        integer :: n

        n = size(x)
        if (n == 0) then
            value = quiet_nan()
            return
        end if

        work = x
        call insertion_sort(work)
        if (mod(n, 2) == 1) then
            value = work((n + 1) / 2)
        else
            value = 0.5_dp * (work(n / 2) + work(n / 2 + 1))
        end if
    end function median_value

    pure subroutine average_ranks(x, ranks)
        real(dp), intent(in) :: x(:)
        real(dp), intent(out) :: ranks(:)
        integer, allocatable :: order(:)
        integer :: i, j, n
        real(dp) :: rank_value

        n = size(x)
        if (size(ranks) /= n) return
        allocate(order(n))
        order = [(i, i = 1, n)]
        call sort_indices(x, order)

        i = 1
        do while (i <= n)
            j = i
            do while (j < n)
                if (x(order(j + 1)) /= x(order(i))) exit
                j = j + 1
            end do
            rank_value = 0.5_dp * real(i + j, dp)
            ranks(order(i:j)) = rank_value
            i = j + 1
        end do
    end subroutine average_ranks

    pure subroutine insertion_sort(x)
        real(dp), intent(inout) :: x(:)
        real(dp) :: key
        integer :: i, j

        do i = 2, size(x)
            key = x(i)
            j = i - 1
            do while (j >= 1)
                if (x(j) <= key) exit
                x(j + 1) = x(j)
                j = j - 1
            end do
            x(j + 1) = key
        end do
    end subroutine insertion_sort

    pure subroutine sort_indices(x, order)
        real(dp), intent(in) :: x(:)
        integer, intent(inout) :: order(:)
        integer :: i, j, key

        do i = 2, size(order)
            key = order(i)
            j = i - 1
            do while (j >= 1)
                if (x(order(j)) <= x(key)) exit
                order(j + 1) = order(j)
                j = j - 1
            end do
            order(j + 1) = key
        end do
    end subroutine sort_indices

end module metrics_utils
