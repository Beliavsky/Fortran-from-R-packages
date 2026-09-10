! SPDX-License-Identifier: GPL-2.0-only
!
! Window functions translated from CRAN signal / Octave-Forge sources.
module signal_windows
    use signal_kinds, only : dp, signal_pi
    use signal_utils, only : bessel_i0, dft_complex
    implicit none
    private

    public :: bartlett_window
    public :: blackman_window
    public :: boxcar_window
    public :: chebyshev_window
    public :: flattop_window
    public :: gaussian_window
    public :: hamming_window
    public :: hanning_window
    public :: kaiser_window
    public :: triangular_window

contains

    pure function bartlett_window(n) result(w)
        integer, intent(in) :: n !! Number of points; must be positive.
        real(dp), allocatable :: w(:)
        integer :: i
        integer :: m

        if (n <= 0) then
            allocate(w(0))
            return
        end if
        allocate(w(n))
        if (n == 1) then
            w = 1.0_dp
            return
        end if
        m = n - 1
        do i = 0, m
            if (i <= m / 2) then
                w(i + 1) = 2.0_dp * real(i, dp) / real(m, dp)
            else
                w(i + 1) = 2.0_dp - 2.0_dp * real(i, dp) / real(m, dp)
            end if
        end do
    end function bartlett_window

    pure function blackman_window(n) result(w)
        integer, intent(in) :: n !! Number of points; must be positive.
        real(dp), allocatable :: w(:)
        real(dp) :: x
        integer :: i

        if (n <= 0) then
            allocate(w(0))
            return
        end if
        allocate(w(n))
        if (n == 1) then
            w = 1.0_dp
            return
        end if
        do i = 0, n - 1
            x = real(i, dp) / real(n - 1, dp)
            w(i + 1) = 0.42_dp - 0.5_dp * cos(2.0_dp * signal_pi * x) + &
                0.08_dp * cos(4.0_dp * signal_pi * x)
        end do
    end function blackman_window

    pure function boxcar_window(n) result(w)
        integer, intent(in) :: n !! Number of points; must be positive.
        real(dp), allocatable :: w(:)

        allocate(w(max(0, n)))
        if (n > 0) w = 1.0_dp
    end function boxcar_window

    pure function flattop_window(n, periodic) result(w)
        integer, intent(in) :: n !! Number of points; must be positive.
        logical, intent(in), optional :: periodic !! True uses the periodic denominator n; default is symmetric.
        real(dp), allocatable :: w(:)
        logical :: use_periodic
        real(dp) :: divisor
        real(dp) :: x
        integer :: i

        if (n <= 0) then
            allocate(w(0))
            return
        end if
        allocate(w(n))
        use_periodic = .false.
        if (present(periodic)) use_periodic = periodic
        if (n == 1) then
            if (use_periodic) then
                w = (1.0_dp - 1.93_dp + 1.29_dp - 0.388_dp + 0.0322_dp) / 4.6402_dp
            else
                w = 1.0_dp
            end if
            return
        end if
        divisor = real(n - 1, dp)
        if (use_periodic) divisor = real(n, dp)
        do i = 0, n - 1
            x = 2.0_dp * signal_pi * real(i, dp) / divisor
            w(i + 1) = (1.0_dp - 1.93_dp * cos(x) + 1.29_dp * cos(2.0_dp * x) - &
                0.388_dp * cos(3.0_dp * x) + 0.0322_dp * cos(4.0_dp * x)) / 4.6402_dp
        end do
    end function flattop_window

    pure function gaussian_window(n, width) result(w)
        integer, intent(in) :: n !! Number of points; must be positive.
        real(dp), intent(in), optional :: width !! Gaussian width parameter; defaults to 2.5 as in signal::gausswin.
        real(dp), allocatable :: w(:)
        real(dp) :: alpha
        real(dp) :: coordinate
        integer :: i

        if (n <= 0) then
            allocate(w(0))
            return
        end if
        alpha = 2.5_dp
        if (present(width)) alpha = width
        allocate(w(n))
        do i = 1, n
            coordinate = real(2 * i - n - 1, dp)
            w(i) = exp(-0.5_dp * (alpha * coordinate / real(n, dp)) ** 2)
        end do
    end function gaussian_window

    pure function hamming_window(n) result(w)
        integer, intent(in) :: n !! Number of points; must be positive.
        real(dp), allocatable :: w(:)
        integer :: i

        if (n <= 0) then
            allocate(w(0))
            return
        end if
        allocate(w(n))
        if (n == 1) then
            w = 1.0_dp
            return
        end if
        do i = 0, n - 1
            w(i + 1) = 0.54_dp - 0.46_dp * cos(2.0_dp * signal_pi * real(i, dp) / real(n - 1, dp))
        end do
    end function hamming_window

    pure function hanning_window(n) result(w)
        integer, intent(in) :: n !! Number of points; must be positive.
        real(dp), allocatable :: w(:)
        integer :: i

        if (n <= 0) then
            allocate(w(0))
            return
        end if
        allocate(w(n))
        if (n == 1) then
            w = 1.0_dp
            return
        end if
        do i = 0, n - 1
            w(i + 1) = 0.5_dp - 0.5_dp * cos(2.0_dp * signal_pi * real(i, dp) / real(n - 1, dp))
        end do
    end function hanning_window

    pure function kaiser_window(n, beta) result(w)
        integer, intent(in) :: n !! Number of points; must be positive.
        real(dp), intent(in) :: beta !! Kaiser shape parameter.
        real(dp), allocatable :: w(:)
        real(dp) :: argument
        integer :: i

        if (n <= 0) then
            allocate(w(0))
            return
        end if
        allocate(w(n))
        if (n == 1) then
            w = 1.0_dp
            return
        end if
        do i = 0, n - 1
            argument = 2.0_dp * beta / real(n - 1, dp) * &
                sqrt(real(i * (n - 1 - i), dp))
            w(i + 1) = bessel_i0(argument) / bessel_i0(beta)
        end do
    end function kaiser_window

    pure function triangular_window(n) result(w)
        integer, intent(in) :: n !! Number of points; must be positive.
        real(dp), allocatable :: w(:)
        real(dp) :: denominator
        real(dp) :: coordinate
        integer :: i

        if (n <= 0) then
            allocate(w(0))
            return
        end if
        allocate(w(n))
        denominator = real(n + modulo(n, 2), dp)
        do i = 1, n
            coordinate = real(2 * i - n - 1, dp)
            w(i) = 1.0_dp - abs(coordinate / denominator)
        end do
    end function triangular_window

    pure function chebyshev_window(n, attenuation_db) result(w)
        integer, intent(in) :: n !! Number of points; must be positive.
        real(dp), intent(in) :: attenuation_db !! Sidelobe attenuation in decibels, matching signal::chebwin.
        real(dp), allocatable :: w(:)
        real(dp), allocatable :: p(:)
        complex(dp), allocatable :: pc(:)
        complex(dp), allocatable :: spectrum(:)
        real(dp) :: beta
        real(dp) :: gamma
        real(dp) :: x
        integer :: i
        integer :: m

        if (n <= 0) then
            allocate(w(0))
            return
        end if
        allocate(w(n))
        if (n == 1) then
            w = 1.0_dp
            return
        end if
        gamma = 10.0_dp ** (-attenuation_db / 20.0_dp)
        beta = cosh(acosh(1.0_dp / gamma) / real(n - 1, dp))
        allocate(p(n), pc(n))
        do i = 0, n - 1
            x = beta * cos(signal_pi * real(i, dp) / real(n, dp))
            p(i + 1) = chebyshev_t(n - 1, x)
            pc(i + 1) = cmplx(p(i + 1), 0.0_dp, dp)
        end do
        if (modulo(n, 2) == 0) then
            do i = 0, n - 1
                pc(i + 1) = pc(i + 1) * exp(cmplx(0.0_dp, signal_pi * real(i, dp) / real(n, dp), dp))
            end do
        end if
        spectrum = dft_complex(pc, .false.)
        if (modulo(n, 2) == 1) then
            m = (n + 1) / 2
            spectrum = spectrum / real(spectrum(1), dp)
            do i = 1, m - 1
                w(i) = real(spectrum(m - i + 1), dp)
            end do
            do i = m, n
                w(i) = real(spectrum(i - m + 1), dp)
            end do
        else
            m = n / 2 + 1
            spectrum = spectrum / real(spectrum(2), dp)
            do i = 1, m - 1
                w(i) = real(spectrum(m - i + 1), dp)
            end do
            do i = m, n
                w(i) = real(spectrum(i - m + 2), dp)
            end do
        end if
    end function chebyshev_window

    pure elemental function chebyshev_t(n, x) result(value)
        integer, intent(in) :: n !! Nonnegative Chebyshev polynomial degree.
        real(dp), intent(in) :: x !! Real evaluation point.
        real(dp) :: value

        if (abs(x) <= 1.0_dp) then
            value = cos(real(n, dp) * acos(x))
        else if (x > 1.0_dp) then
            value = cosh(real(n, dp) * acosh(x))
        else
            value = (-1.0_dp) ** n * cosh(real(n, dp) * acosh(-x))
        end if
    end function chebyshev_t

end module signal_windows
