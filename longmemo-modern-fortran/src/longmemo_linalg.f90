! Part of the modern Fortran translation of longmemo 1.1-4.
! This file was created or modified for the Fortran project on 2026-07-23.
! Original longmemo authors retain copyright; see ORIGINAL_PACKAGE.txt.
! SPDX-License-Identifier: GPL-2.0-or-later

module longmemo_linalg
    use longmemo_kinds, only : dp
    implicit none
    private

    public :: inverse_spd, solve_spd, cross_product

contains

    subroutine cross_product(x, a)
        real(dp), intent(in) :: x(:, :)
        real(dp), allocatable, intent(out) :: a(:, :)

        allocate(a(size(x, 2), size(x, 2)))
        a = matmul(transpose(x), x)
    end subroutine cross_product


    subroutine inverse_spd(a, ainv, info)
        real(dp), intent(in) :: a(:, :)
        real(dp), allocatable, intent(out) :: ainv(:, :)
        integer, intent(out) :: info
        real(dp), allocatable :: l(:, :), e(:), y(:), x(:)
        integer :: n, j

        n = size(a, 1)
        if (size(a, 2) /= n) then
            info = -1
            allocate(ainv(0, 0))
            return
        end if

        allocate(l(n, n), ainv(n, n), e(n), y(n), x(n))
        call cholesky_lower(a, l, info)
        if (info /= 0) then
            ainv = 0.0_dp
            return
        end if

        do j = 1, n
            e = 0.0_dp
            e(j) = 1.0_dp
            call forward_substitution(l, e, y)
            call backward_substitution(transpose(l), y, x)
            ainv(:, j) = x
        end do
        ainv = 0.5_dp*(ainv + transpose(ainv))
    end subroutine inverse_spd


    subroutine solve_spd(a, b, x, info)
        real(dp), intent(in) :: a(:, :), b(:)
        real(dp), allocatable, intent(out) :: x(:)
        integer, intent(out) :: info
        real(dp), allocatable :: l(:, :), y(:)
        integer :: n

        n = size(a, 1)
        allocate(x(n), l(n, n), y(n))
        if (size(a, 2) /= n .or. size(b) /= n) then
            info = -1
            x = 0.0_dp
            return
        end if

        call cholesky_lower(a, l, info)
        if (info /= 0) then
            x = 0.0_dp
            return
        end if
        call forward_substitution(l, b, y)
        call backward_substitution(transpose(l), y, x)
    end subroutine solve_spd


    subroutine cholesky_lower(a, l, info)
        real(dp), intent(in) :: a(:, :)
        real(dp), intent(out) :: l(:, :)
        integer, intent(out) :: info
        real(dp) :: s
        integer :: n, i, j, k

        n = size(a, 1)
        l = 0.0_dp
        info = 0
        do i = 1, n
            do j = 1, i
                s = a(i, j)
                do k = 1, j - 1
                    s = s - l(i, k)*l(j, k)
                end do
                if (i == j) then
                    if (s <= epsilon(1.0_dp)*max(1.0_dp, abs(a(i, i)))) then
                        info = i
                        return
                    end if
                    l(i, j) = sqrt(s)
                else
                    l(i, j) = s/l(j, j)
                end if
            end do
        end do
    end subroutine cholesky_lower


    subroutine forward_substitution(l, b, x)
        real(dp), intent(in) :: l(:, :), b(:)
        real(dp), intent(out) :: x(:)
        integer :: i

        do i = 1, size(b)
            x(i) = (b(i) - dot_product(l(i, 1:i - 1), x(1:i - 1)))/l(i, i)
        end do
    end subroutine forward_substitution


    subroutine backward_substitution(u, b, x)
        real(dp), intent(in) :: u(:, :), b(:)
        real(dp), intent(out) :: x(:)
        integer :: i, n

        n = size(b)
        do i = n, 1, -1
            x(i) = (b(i) - dot_product(u(i, i + 1:n), x(i + 1:n)))/u(i, i)
        end do
    end subroutine backward_substitution

end module longmemo_linalg
