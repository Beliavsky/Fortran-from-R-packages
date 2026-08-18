! Translation of bivgeom 1.0 computational code.
! Upstream DESCRIPTION declares License: GPL. See LICENSE.md.
module bivgeom_math
    use, intrinsic :: iso_fortran_env, only : int64
    use bivgeom_kinds, only : dp
    implicit none
    private

    public :: qgeom_theta
    public :: seed_rng
    public :: solve_3x3
    public :: upper_string

contains

    integer function qgeom_theta(prob, theta) result(q)
        real(dp), intent(in) :: prob
        real(dp), intent(in) :: theta
        real(dp) :: z

        if (prob <= 0.0_dp) then
            q = 0
            return
        end if
        if (prob >= 1.0_dp) then
            q = huge(1)
            return
        end if
        if (theta <= 0.0_dp) then
            q = 0
            return
        end if
        if (theta >= 1.0_dp) then
            q = huge(1)
            return
        end if

        z = log1p_safe(-prob) / log(theta) - 1.0_dp
        q = max(0, ceiling(z - 16.0_dp * epsilon(1.0_dp)))
    end function qgeom_theta

    subroutine seed_rng(seed)
        integer(int64), intent(in) :: seed
        integer :: n, i
        integer, allocatable :: put(:)
        integer(int64) :: x
        integer(int64), parameter :: modulus = 2147483647_int64
        integer(int64), parameter :: multiplier = 48271_int64

        call random_seed(size=n)
        allocate(put(n))
        x = modulo(abs(seed), modulus - 1_int64) + 1_int64
        do i = 1, n
            x = modulo(multiplier * x, modulus)
            put(i) = int(x, kind(put(i)))
        end do
        call random_seed(put=put)
    end subroutine seed_rng

    subroutine solve_3x3(a, b, x, ok)
        real(dp), intent(in) :: a(3, 3), b(3)
        real(dp), intent(out) :: x(3)
        logical, intent(out) :: ok
        real(dp) :: aug(3, 4), tmp(4), factor
        integer :: i, j, k, piv

        aug(:, 1:3) = a
        aug(:, 4) = b
        ok = .true.

        do k = 1, 3
            piv = k
            do i = k + 1, 3
                if (abs(aug(i, k)) > abs(aug(piv, k))) piv = i
            end do
            if (abs(aug(piv, k)) <= 100.0_dp * epsilon(1.0_dp)) then
                ok = .false.
                x = 0.0_dp
                return
            end if
            if (piv /= k) then
                tmp = aug(k, :)
                aug(k, :) = aug(piv, :)
                aug(piv, :) = tmp
            end if
            do i = k + 1, 3
                factor = aug(i, k) / aug(k, k)
                do j = k, 4
                    aug(i, j) = aug(i, j) - factor * aug(k, j)
                end do
            end do
        end do

        do i = 3, 1, -1
            x(i) = aug(i, 4)
            do j = i + 1, 3
                x(i) = x(i) - aug(i, j) * x(j)
            end do
            x(i) = x(i) / aug(i, i)
        end do
    end subroutine solve_3x3

    function upper_string(s) result(out)
        character(len=*), intent(in) :: s
        character(len=len(s)) :: out
        integer :: i, c

        out = s
        do i = 1, len(s)
            c = iachar(out(i:i))
            if (c >= iachar('a') .and. c <= iachar('z')) then
                out(i:i) = achar(c - iachar('a') + iachar('A'))
            end if
        end do
    end function upper_string

    pure real(dp) function log1p_safe(x) result(y)
        real(dp), intent(in) :: x
        real(dp) :: term, sum
        integer :: k

        if (abs(x) > 1.0e-4_dp) then
            y = log(1.0_dp + x)
            return
        end if

        term = x
        sum = 0.0_dp
        do k = 1, 50
            if (mod(k, 2) == 1) then
                sum = sum + term / real(k, dp)
            else
                sum = sum - term / real(k, dp)
            end if
            term = term * x
            if (abs(term) <= epsilon(1.0_dp) * max(1.0_dp, abs(sum))) exit
        end do
        y = sum
    end function log1p_safe

end module bivgeom_math
