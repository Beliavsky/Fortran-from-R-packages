! SPDX-License-Identifier: GPL-2.0-or-later
module mla_linalg
    use mla_kinds, only : dp
    implicit none
    private

    public :: packed_size, pack_upper, unpack_upper
    public :: solve_packed_spd, invert_packed_spd, solve_spd, invert_spd

contains

    pure integer function packed_size(n) result(np)
        integer, intent(in) :: n
        np = n * (n + 1) / 2
    end function packed_size

    pure subroutine pack_upper(a, p)
        real(dp), intent(in) :: a(:, :)
        real(dp), intent(out) :: p(:)
        integer :: i, j, k, n

        n = size(a, 1)
        k = 0
        do j = 1, n
            do i = 1, j
                k = k + 1
                p(k) = a(i, j)
            end do
        end do
    end subroutine pack_upper

    pure subroutine unpack_upper(p, a)
        real(dp), intent(in) :: p(:)
        real(dp), intent(out) :: a(:, :)
        integer :: i, j, k, n

        n = size(a, 1)
        a = 0.0_dp
        k = 0
        do j = 1, n
            do i = 1, j
                k = k + 1
                a(i, j) = p(k)
                a(j, i) = p(k)
            end do
        end do
    end subroutine unpack_upper

    subroutine cholesky_factor(a, l, ok, logdet)
        real(dp), intent(in) :: a(:, :)
        real(dp), intent(out) :: l(:, :)
        logical, intent(out) :: ok
        real(dp), intent(out), optional :: logdet
        integer :: i, j, k, n
        real(dp) :: s, ld, tol

        n = size(a, 1)
        l = 0.0_dp
        ld = 0.0_dp
        tol = epsilon(1.0_dp)
        ok = .false.

        do i = 1, n
            do j = 1, i
                s = a(i, j)
                do k = 1, j - 1
                    s = s - l(i, k) * l(j, k)
                end do
                if (i == j) then
                    if (s <= tol * max(1.0_dp, abs(a(i, i)))) return
                    l(i, j) = sqrt(s)
                    ld = ld + 2.0_dp * log(l(i, j))
                else
                    l(i, j) = s / l(j, j)
                end if
            end do
        end do

        ok = .true.
        if (present(logdet)) logdet = ld
    end subroutine cholesky_factor

    subroutine solve_spd(a, b, x, ok)
        real(dp), intent(in) :: a(:, :), b(:)
        real(dp), intent(out) :: x(:)
        logical, intent(out) :: ok
        real(dp), allocatable :: l(:, :), y(:)
        integer :: i, k, n

        n = size(a, 1)
        allocate(l(n, n), y(n))
        call cholesky_factor(a, l, ok)
        if (.not. ok) then
            x = 0.0_dp
            return
        end if

        do i = 1, n
            y(i) = b(i)
            do k = 1, i - 1
                y(i) = y(i) - l(i, k) * y(k)
            end do
            y(i) = y(i) / l(i, i)
        end do

        do i = n, 1, -1
            x(i) = y(i)
            do k = i + 1, n
                x(i) = x(i) - l(k, i) * x(k)
            end do
            x(i) = x(i) / l(i, i)
        end do
    end subroutine solve_spd

    subroutine invert_spd(a, ainv, ok, logdet)
        real(dp), intent(in) :: a(:, :)
        real(dp), intent(out) :: ainv(:, :)
        logical, intent(out) :: ok
        real(dp), intent(out), optional :: logdet
        real(dp), allocatable :: l(:, :), rhs(:), sol(:)
        real(dp) :: ld
        integer :: j, n

        n = size(a, 1)
        allocate(l(n, n), rhs(n), sol(n))
        call cholesky_factor(a, l, ok, ld)
        if (.not. ok) then
            ainv = 0.0_dp
            if (present(logdet)) logdet = 0.0_dp
            return
        end if

        ainv = 0.0_dp
        do j = 1, n
            rhs = 0.0_dp
            rhs(j) = 1.0_dp
            call solve_with_factor(l, rhs, sol)
            ainv(:, j) = sol
        end do
        ainv = 0.5_dp * (ainv + transpose(ainv))
        if (present(logdet)) logdet = ld
    contains
        subroutine solve_with_factor(ll, bb, xx)
            real(dp), intent(in) :: ll(:, :), bb(:)
            real(dp), intent(out) :: xx(:)
            real(dp) :: yy(size(bb))
            integer :: ii, kk, nn
            nn = size(bb)
            do ii = 1, nn
                yy(ii) = bb(ii)
                do kk = 1, ii - 1
                    yy(ii) = yy(ii) - ll(ii, kk) * yy(kk)
                end do
                yy(ii) = yy(ii) / ll(ii, ii)
            end do
            do ii = nn, 1, -1
                xx(ii) = yy(ii)
                do kk = ii + 1, nn
                    xx(ii) = xx(ii) - ll(kk, ii) * xx(kk)
                end do
                xx(ii) = xx(ii) / ll(ii, ii)
            end do
        end subroutine solve_with_factor
    end subroutine invert_spd

    subroutine solve_packed_spd(p, rhs, x, idpos)
        real(dp), intent(in) :: p(:), rhs(:)
        real(dp), intent(out) :: x(:)
        integer, intent(out) :: idpos
        real(dp), allocatable :: a(:, :)
        logical :: ok
        integer :: n

        n = size(rhs)
        allocate(a(n, n))
        call unpack_upper(p, a)
        call solve_spd(a, rhs, x, ok)
        if (ok) then
            idpos = 0
        else
            idpos = 1
        end if
    end subroutine solve_packed_spd

    subroutine invert_packed_spd(p, pinv, ier, logdet)
        real(dp), intent(in) :: p(:)
        real(dp), intent(out) :: pinv(:)
        integer, intent(out) :: ier
        real(dp), intent(out), optional :: logdet
        real(dp), allocatable :: a(:, :), ai(:, :)
        real(dp) :: ld
        logical :: ok
        integer :: n

        n = int((sqrt(1.0_dp + 8.0_dp * real(size(p), dp)) - 1.0_dp) / 2.0_dp)
        allocate(a(n, n), ai(n, n))
        call unpack_upper(p, a)
        call invert_spd(a, ai, ok, ld)
        if (ok) then
            ier = 0
            call pack_upper(ai, pinv)
        else
            ier = -1
            pinv = 0.0_dp
        end if
        if (present(logdet)) logdet = ld
    end subroutine invert_packed_spd

end module mla_linalg
