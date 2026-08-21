! SPDX-License-Identifier: GPL-3.0-or-later
module cubature_utils
    use cubature_kinds, only : dp, i8
    use cubature_types, only : errors_converged
    implicit none
    private
    public :: errors_converged, halton_point, primes_first, clenshaw_curtis_rule
    public :: safe_pow_int, volume_box

contains

    integer(i8) function safe_pow_int(base, exponent, overflow) result(v)
        integer, intent(in) :: base, exponent
        logical, intent(out) :: overflow
        integer :: k
        v = 1_i8
        overflow = .false.
        if (exponent < 0 .or. base < 0) then
            overflow = .true.
            return
        end if
        do k = 1, exponent
            if (v > huge(v) / max(1_i8, int(base, i8))) then
                overflow = .true.
                v = huge(v)
                return
            end if
            v = v * int(base, i8)
        end do
    end function safe_pow_int

    real(dp) function volume_box(lower, upper) result(v)
        real(dp), intent(in) :: lower(:), upper(:)
        if (size(lower) /= size(upper)) then
            v = 0.0_dp
        else
            v = product(upper - lower)
        end if
    end function volume_box

    subroutine primes_first(n, primes)
        integer, intent(in) :: n
        integer, intent(out) :: primes(n)
        integer :: candidate, count, j
        logical :: prime
        if (n <= 0) return
        primes(1) = 2
        count = 1
        candidate = 3
        do while (count < n)
            prime = .true.
            do j = 1, count
                if (primes(j) * primes(j) > candidate) exit
                if (mod(candidate, primes(j)) == 0) then
                    prime = .false.
                    exit
                end if
            end do
            if (prime) then
                count = count + 1
                primes(count) = candidate
            end if
            candidate = candidate + 2
        end do
    end subroutine primes_first

    subroutine halton_point(index, dim, x, shift)
        integer(i8), intent(in) :: index
        integer, intent(in) :: dim
        real(dp), intent(out) :: x(dim)
        real(dp), intent(in), optional :: shift(dim)
        integer, allocatable :: primes(:)
        integer(i8) :: n
        integer :: d, b
        real(dp) :: f, r
        allocate(primes(dim))
        call primes_first(dim, primes)
        do d = 1, dim
            b = primes(d)
            n = index
            f = 1.0_dp / real(b, dp)
            r = 0.0_dp
            do while (n > 0_i8)
                r = r + f * real(modulo(n, int(b, i8)), dp)
                n = n / int(b, i8)
                f = f / real(b, dp)
            end do
            if (present(shift)) then
                x(d) = modulo(r + shift(d), 1.0_dp)
            else
                x(d) = r
            end if
        end do
    end subroutine halton_point

    subroutine clenshaw_curtis_rule(n, x, w)
        integer, intent(in) :: n
        real(dp), allocatable, intent(out) :: x(:), w(:)
        integer :: nint, j, k
        real(dp), allocatable :: theta(:), v(:)
        real(dp), parameter :: pi = acos(-1.0_dp)

        if (n < 2) then
            allocate(x(1), w(1))
            x = 0.0_dp
            w = 2.0_dp
            return
        end if
        nint = n - 1
        allocate(x(n), w(n), theta(n))
        do j = 0, nint
            theta(j + 1) = pi * real(j, dp) / real(nint, dp)
            x(j + 1) = cos(theta(j + 1))
        end do
        w = 0.0_dp
        if (nint == 1) then
            w = 1.0_dp
            return
        end if
        allocate(v(nint - 1))
        v = 1.0_dp
        if (mod(nint, 2) == 0) then
            w(1) = 1.0_dp / real(nint * nint - 1, dp)
            w(n) = w(1)
            do k = 1, nint / 2 - 1
                do j = 1, nint - 1
                    v(j) = v(j) - 2.0_dp * cos(2.0_dp * real(k, dp) * theta(j + 1)) / &
                           real(4 * k * k - 1, dp)
                end do
            end do
            do j = 1, nint - 1
                v(j) = v(j) - cos(real(nint, dp) * theta(j + 1)) / real(nint * nint - 1, dp)
            end do
        else
            w(1) = 1.0_dp / real(nint * nint, dp)
            w(n) = w(1)
            do k = 1, (nint - 1) / 2
                do j = 1, nint - 1
                    v(j) = v(j) - 2.0_dp * cos(2.0_dp * real(k, dp) * theta(j + 1)) / &
                           real(4 * k * k - 1, dp)
                end do
            end do
        end if
        do j = 1, nint - 1
            w(j + 1) = 2.0_dp * v(j) / real(nint, dp)
        end do
    end subroutine clenshaw_curtis_rule

end module cubature_utils
