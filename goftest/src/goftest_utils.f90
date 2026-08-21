! SPDX-License-Identifier: GPL-2.0-or-later
module goftest_utils
    use goftest_kinds, only : dp
    implicit none
    private

    public :: sort_real
    public :: clamp01
    public :: shuffle_int

contains

    subroutine sort_real(x)
        real(dp), intent(inout) :: x(:)
        integer :: i, j
        real(dp) :: key

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
    end subroutine sort_real

    pure real(dp) function clamp01(x) result(y)
        real(dp), intent(in) :: x
        y = min(1.0_dp, max(0.0_dp, x))
    end function clamp01

    subroutine shuffle_int(x)
        integer, intent(inout) :: x(:)
        integer :: i, j, tmp
        real(dp) :: u

        do i = size(x), 2, -1
            call random_number(u)
            j = 1 + int(u * real(i, dp))
            if (j > i) j = i
            tmp = x(i)
            x(i) = x(j)
            x(j) = tmp
        end do
    end subroutine shuffle_int

end module goftest_utils
