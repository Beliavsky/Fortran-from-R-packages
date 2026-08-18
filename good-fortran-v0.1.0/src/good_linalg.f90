! SPDX-License-Identifier: GPL-2.0-or-later
! Modern Fortran translation of computational code from R package good 1.0.2.

module good_linalg
    use good_kinds, only : dp
    implicit none
    private

    public :: invert_matrix
    public :: outer_product

contains

    pure function outer_product(x, y) result(a)
        real(dp), intent(in) :: x(:), y(:)
        real(dp) :: a(size(x), size(y))
        integer :: i, j

        do j = 1, size(y)
            do i = 1, size(x)
                a(i, j) = x(i) * y(j)
            end do
        end do
    end function outer_product

    subroutine invert_matrix(a, ainv, status)
        real(dp), intent(in) :: a(:, :)
        real(dp), intent(out) :: ainv(:, :)
        integer, intent(out) :: status

        real(dp), allocatable :: aug(:, :), rowtmp(:)
        real(dp) :: pivot, factor, scale
        integer :: n, i, k, p

        n = size(a, 1)
        status = 0
        if (size(a, 2) /= n .or. size(ainv, 1) /= n .or. size(ainv, 2) /= n) then
            status = -1
            return
        end if

        allocate(aug(n, 2 * n), rowtmp(2 * n))
        aug = 0.0_dp
        aug(:, 1:n) = a
        do i = 1, n
            aug(i, n + i) = 1.0_dp
        end do

        scale = max(1.0_dp, maxval(abs(a)))
        do k = 1, n
            p = k
            do i = k + 1, n
                if (abs(aug(i, k)) > abs(aug(p, k))) p = i
            end do
            if (abs(aug(p, k)) <= 100.0_dp * epsilon(1.0_dp) * scale) then
                status = 1
                ainv = 0.0_dp
                return
            end if
            if (p /= k) then
                rowtmp = aug(k, :)
                aug(k, :) = aug(p, :)
                aug(p, :) = rowtmp
            end if

            pivot = aug(k, k)
            aug(k, :) = aug(k, :) / pivot
            do i = 1, n
                if (i == k) cycle
                factor = aug(i, k)
                if (abs(factor) > tiny(1.0_dp)) aug(i, :) = aug(i, :) - factor * aug(k, :)
            end do
        end do

        ainv = aug(:, n + 1:2 * n)
    end subroutine invert_matrix

end module good_linalg
