! SPDX-License-Identifier: GPL-2.0-only
!
! Numerical support for the Fortran translation of CRAN package signal.
! Upstream signal sources are GPL-2 and include work by Paul Kienzle,
! David Billinghurst, John W. Eaton, and other Octave/Octave-Forge authors.
module signal_utils
    use signal_kinds, only : dp, signal_pi
    implicit none
    private

    public :: convolve
    public :: dft_real
    public :: dft_complex
    public :: inverse_dft
    public :: polyval_real_complex
    public :: polyval_complex
    public :: polynomial_from_roots
    public :: polynomial_roots
    public :: sinc_value
    public :: solve_linear_system
    public :: bessel_i0
    public :: factorial_real
    public :: reverse_real
    public :: reverse_complex

contains

    pure function convolve(x, y) result(z)
        real(dp), intent(in) :: x(:) !! First real sequence to convolve.
        real(dp), intent(in) :: y(:) !! Second real sequence to convolve.
        real(dp), allocatable :: z(:)
        integer :: i
        integer :: j

        if (size(x) == 0 .or. size(y) == 0) then
            allocate(z(0))
            return
        end if
        allocate(z(size(x) + size(y) - 1))
        z = 0.0_dp
        do i = 1, size(x)
            do j = 1, size(y)
                z(i + j - 1) = z(i + j - 1) + x(i) * y(j)
            end do
        end do
    end function convolve

    pure function dft_real(x) result(y)
        real(dp), intent(in) :: x(:) !! Real samples whose forward discrete Fourier transform is required.
        complex(dp), allocatable :: y(:)
        complex(dp), allocatable :: xc(:)

        allocate(xc(size(x)))
        xc = cmplx(x, 0.0_dp, dp)
        y = dft_complex(xc, .false.)
    end function dft_real

    pure function dft_complex(x, inverse) result(y)
        complex(dp), intent(in) :: x(:) !! Complex samples to transform.
        logical, intent(in) :: inverse !! True requests the unnormalised inverse-sign transform.
        complex(dp), allocatable :: y(:)
        integer :: j
        integer :: k
        real(dp) :: sign_phase
        real(dp) :: phase

        allocate(y(size(x)))
        y = cmplx(0.0_dp, 0.0_dp, dp)
        if (size(x) == 0) return
        sign_phase = -1.0_dp
        if (inverse) sign_phase = 1.0_dp
        do k = 1, size(x)
            do j = 1, size(x)
                phase = sign_phase * 2.0_dp * signal_pi * real((j - 1) * (k - 1), dp) / real(size(x), dp)
                y(k) = y(k) + x(j) * exp(cmplx(0.0_dp, phase, dp))
            end do
        end do
    end function dft_complex

    pure function inverse_dft(x) result(y)
        complex(dp), intent(in) :: x(:) !! Complex frequency-domain sequence to invert.
        complex(dp), allocatable :: y(:)

        y = dft_complex(x, .true.)
        if (size(x) > 0) y = y / real(size(x), dp)
    end function inverse_dft

    pure function polyval_real_complex(coef, z) result(y)
        real(dp), intent(in) :: coef(:) !! Polynomial coefficients in descending powers.
        complex(dp), intent(in) :: z(:) !! Complex evaluation points.
        complex(dp), allocatable :: y(:)
        integer :: i

        allocate(y(size(z)))
        y = cmplx(0.0_dp, 0.0_dp, dp)
        do i = 1, size(coef)
            y = y * z + coef(i)
        end do
    end function polyval_real_complex

    pure function polyval_complex(coef, z) result(y)
        complex(dp), intent(in) :: coef(:) !! Complex polynomial coefficients in descending powers.
        complex(dp), intent(in) :: z(:) !! Complex evaluation points.
        complex(dp), allocatable :: y(:)
        integer :: i

        allocate(y(size(z)))
        y = cmplx(0.0_dp, 0.0_dp, dp)
        do i = 1, size(coef)
            y = y * z + coef(i)
        end do
    end function polyval_complex

    pure function polynomial_from_roots(r) result(coef)
        complex(dp), intent(in) :: r(:) !! Polynomial roots; output has leading coefficient one.
        complex(dp), allocatable :: coef(:)
        complex(dp), allocatable :: next(:)
        integer :: i
        integer :: j

        allocate(coef(1))
        coef(1) = cmplx(1.0_dp, 0.0_dp, dp)
        do i = 1, size(r)
            allocate(next(size(coef) + 1))
            next = cmplx(0.0_dp, 0.0_dp, dp)
            do j = 1, size(coef)
                next(j) = next(j) + coef(j)
                next(j + 1) = next(j + 1) - r(i) * coef(j)
            end do
            call move_alloc(next, coef)
        end do
    end function polynomial_from_roots

    pure function polynomial_roots(coef) result(r)
        real(dp), intent(in) :: coef(:) !! Real polynomial coefficients in descending powers; first coefficient must be nonzero.
        complex(dp), allocatable :: r(:)
        complex(dp), allocatable :: old(:)
        complex(dp) :: denom
        complex(dp) :: delta
        complex(dp) :: value
        real(dp) :: radius
        real(dp) :: tol
        integer :: degree
        integer :: i
        integer :: j
        integer :: iter

        degree = size(coef) - 1
        if (degree <= 0 .or. abs(coef(1)) <= tiny(1.0_dp)) then
            allocate(r(0))
            return
        end if
        allocate(r(degree), old(degree))
        radius = 1.0_dp + maxval(abs(coef(2:)) / abs(coef(1)))
        do i = 1, degree
            r(i) = radius * exp(cmplx(0.0_dp, 2.0_dp * signal_pi * real(i - 1, dp) / real(degree, dp), dp))
        end do
        tol = 100.0_dp * epsilon(1.0_dp)
        do iter = 1, 2000
            old = r
            do i = 1, degree
                value = cmplx(coef(1), 0.0_dp, dp)
                do j = 2, size(coef)
                    value = value * old(i) + coef(j)
                end do
                denom = cmplx(1.0_dp, 0.0_dp, dp)
                do j = 1, degree
                    if (j /= i) denom = denom * (old(i) - old(j))
                end do
                if (abs(denom) < tiny(1.0_dp)) denom = denom + cmplx(tiny(1.0_dp), tiny(1.0_dp), dp)
                delta = value / denom
                r(i) = old(i) - delta
            end do
            if (maxval(abs(r - old)) <= tol * max(1.0_dp, maxval(abs(r)))) exit
        end do
    end function polynomial_roots

    pure elemental function sinc_value(x) result(y)
        real(dp), intent(in) :: x !! Dimensionless sinc argument using sin(pi*x)/(pi*x).
        real(dp) :: y

        if (abs(x) <= tiny(1.0_dp)) then
            y = 1.0_dp
        else
            y = sin(signal_pi * x) / (signal_pi * x)
        end if
    end function sinc_value

    pure subroutine solve_linear_system(a, b, x, ok)
        real(dp), intent(in) :: a(:, :) !! Square coefficient matrix.
        real(dp), intent(in) :: b(:) !! Right-hand side vector with length equal to the matrix order.
        real(dp), allocatable, intent(out) :: x(:) !! Solution vector, or zeros when the system is singular.
        logical, intent(out) :: ok !! True when pivoting found a nonsingular system.
        real(dp), allocatable :: aug(:, :)
        real(dp), allocatable :: rowtmp(:)
        real(dp) :: factor
        real(dp) :: pivot_abs
        integer :: i
        integer :: j
        integer :: k
        integer :: n
        integer :: pivot

        n = size(b)
        allocate(x(n))
        x = 0.0_dp
        ok = size(a, 1) == n .and. size(a, 2) == n
        if (.not. ok) return
        allocate(aug(n, n + 1), rowtmp(n + 1))
        aug(:, 1:n) = a
        aug(:, n + 1) = b
        do k = 1, n
            pivot = k
            pivot_abs = abs(aug(k, k))
            do i = k + 1, n
                if (abs(aug(i, k)) > pivot_abs) then
                    pivot = i
                    pivot_abs = abs(aug(i, k))
                end if
            end do
            if (pivot_abs <= 100.0_dp * epsilon(1.0_dp)) then
                ok = .false.
                return
            end if
            if (pivot /= k) then
                rowtmp = aug(k, :)
                aug(k, :) = aug(pivot, :)
                aug(pivot, :) = rowtmp
            end if
            do i = k + 1, n
                factor = aug(i, k) / aug(k, k)
                aug(i, k:n + 1) = aug(i, k:n + 1) - factor * aug(k, k:n + 1)
            end do
        end do
        do i = n, 1, -1
            x(i) = aug(i, n + 1)
            do j = i + 1, n
                x(i) = x(i) - aug(i, j) * x(j)
            end do
            x(i) = x(i) / aug(i, i)
        end do
    end subroutine solve_linear_system

    pure elemental function bessel_i0(x) result(y)
        real(dp), intent(in) :: x !! Real argument of the modified Bessel function I0.
        real(dp) :: y
        real(dp) :: ax
        real(dp) :: t

        ax = abs(x)
        if (ax < 3.75_dp) then
            t = (ax / 3.75_dp) ** 2
            y = 1.0_dp + t * (3.5156229_dp + t * (3.0899424_dp + t * (1.2067492_dp + &
                t * (0.2659732_dp + t * (0.0360768_dp + t * 0.0045813_dp)))))
        else
            t = 3.75_dp / ax
            y = exp(ax) / sqrt(ax) * (0.39894228_dp + t * (0.01328592_dp + t * (0.00225319_dp + &
                t * (-0.00157565_dp + t * (0.00916281_dp + t * (-0.02057706_dp + &
                t * (0.02635537_dp + t * (-0.01647633_dp + t * 0.00392377_dp))))))))
        end if
    end function bessel_i0

    pure elemental function factorial_real(n) result(value)
        integer, intent(in) :: n !! Nonnegative integer whose factorial is returned as real(dp).
        real(dp) :: value
        integer :: i

        value = 1.0_dp
        do i = 2, n
            value = value * real(i, dp)
        end do
    end function factorial_real

    pure function reverse_real(x) result(y)
        real(dp), intent(in) :: x(:) !! Real vector to return in reverse order.
        real(dp), allocatable :: y(:)
        integer :: i

        allocate(y(size(x)))
        do i = 1, size(x)
            y(i) = x(size(x) - i + 1)
        end do
    end function reverse_real

    pure function reverse_complex(x) result(y)
        complex(dp), intent(in) :: x(:) !! Complex vector to return in reverse order.
        complex(dp), allocatable :: y(:)
        integer :: i

        allocate(y(size(x)))
        do i = 1, size(x)
            y(i) = x(size(x) - i + 1)
        end do
    end function reverse_complex

end module signal_utils
