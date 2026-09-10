! SPDX-License-Identifier: GPL-2.0-only
!
! Core filtering and filter-object operations translated from CRAN signal.
module signal_filters
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
    use signal_kinds, only : dp, signal_pi
    use signal_types, only : arma_filter, zpg_filter, filter_order, fft_filter, median_filter
    use signal_utils, only : convolve, polynomial_from_roots, polynomial_roots
    implicit none
    private

    public :: make_arma
    public :: make_ma
    public :: make_zpg
    public :: arma_from_zpg
    public :: zpg_from_arma
    public :: make_filter_order
    public :: make_fft_filter
    public :: make_median_filter
    public :: make_spencer_filter
    public :: unit_phasor_degrees
    public :: filter_signal
    public :: fft_filter_signal
    public :: zero_phase_filter
    public :: median_filter_signal
    public :: spencer_smooth
    public :: levinson_durbin
    public :: unwrap_phase
    public :: fractional_difference

contains

    pure function make_arma(b, a) result(filt)
        real(dp), intent(in) :: b(:) !! Numerator coefficients in descending powers of z for transfer-function operations.
        real(dp), intent(in) :: a(:) !! Denominator coefficients; a(1) is the filter normalisation coefficient.
        type(arma_filter) :: filt

        filt%b = b
        filt%a = a
    end function make_arma

    pure function make_ma(b) result(filt)
        real(dp), intent(in) :: b(:) !! FIR numerator coefficients.
        type(arma_filter) :: filt

        filt%b = b
        filt%a = [1.0_dp]
    end function make_ma

    pure function make_zpg(zero, pole, gain) result(filt)
        complex(dp), intent(in) :: zero(:) !! Filter zeros in the complex plane.
        complex(dp), intent(in) :: pole(:) !! Filter poles in the complex plane.
        real(dp), intent(in) :: gain !! Scalar zero-pole-gain multiplier.
        type(zpg_filter) :: filt

        filt%zero = zero
        filt%pole = pole
        filt%gain = gain
    end function make_zpg

    pure function arma_from_zpg(filt) result(arma)
        type(zpg_filter), intent(in) :: filt !! Zero-pole-gain filter to convert to polynomial coefficients.
        type(arma_filter) :: arma
        complex(dp), allocatable :: bc(:)
        complex(dp), allocatable :: ac(:)

        bc = polynomial_from_roots(filt%zero)
        ac = polynomial_from_roots(filt%pole)
        allocate(arma%b(size(bc)), arma%a(size(ac)))
        arma%b = real(filt%gain * bc, dp)
        arma%a = real(ac, dp)
    end function arma_from_zpg

    pure function zpg_from_arma(filt) result(zpg)
        type(arma_filter), intent(in) :: filt !! Polynomial transfer-function coefficients to convert to zero-pole-gain form.
        type(zpg_filter) :: zpg

        zpg%zero = polynomial_roots(filt%b)
        zpg%pole = polynomial_roots(filt%a)
        if (size(filt%b) > 0 .and. size(filt%a) > 0 .and. abs(filt%a(1)) > tiny(1.0_dp)) then
            zpg%gain = filt%b(1) / filt%a(1)
        else
            zpg%gain = 0.0_dp
        end if
    end function zpg_from_arma

    pure function make_filter_order(n, wc, filter_type, rp, rs, beta) result(order_spec)
        integer, intent(in) :: n !! Integer filter order.
        real(dp), intent(in) :: wc(:) !! One or two normalised critical frequencies.
        character(len=*), intent(in) :: filter_type !! Filter kind such as low, high, stop, pass, DC-0, or DC-1.
        real(dp), intent(in), optional :: rp !! Optional pass-band ripple in decibels.
        real(dp), intent(in), optional :: rs !! Optional stop-band attenuation in decibels.
        real(dp), intent(in), optional :: beta !! Optional Kaiser window beta parameter.
        type(filter_order) :: order_spec

        order_spec%n = n
        order_spec%wc = wc
        order_spec%filter_type = filter_type
        if (present(rp)) order_spec%rp = rp
        if (present(rs)) order_spec%rs = rs
        if (present(beta)) order_spec%beta = beta
    end function make_filter_order

    pure function make_fft_filter(b, n) result(filt)
        real(dp), intent(in) :: b(:) !! FIR coefficients.
        integer, intent(in), optional :: n !! Requested overlap-add transform length; retained as metadata in this translation.
        type(fft_filter) :: filt

        filt%b = b
        filt%n = 0
        if (present(n)) filt%n = n
    end function make_fft_filter

    pure function make_median_filter(n) result(filt)
        integer, intent(in), optional :: n !! Odd running-median width; defaults to three.
        type(median_filter) :: filt

        filt%n = 3
        if (present(n)) filt%n = n
    end function make_median_filter

    pure function make_spencer_filter() result(filt)
        type(arma_filter) :: filt

        filt%b = real([-3, -6, -5, 3, 21, 46, 67, 74, 67, 46, 21, 3, -5, -6, -3], dp) / 320.0_dp
        filt%a = [1.0_dp]
    end function make_spencer_filter

    pure elemental function unit_phasor_degrees(degrees) result(value)
        real(dp), intent(in) :: degrees !! Angle in degrees.
        complex(dp) :: value

        value = exp(cmplx(0.0_dp, degrees * signal_pi / 180.0_dp, dp))
    end function unit_phasor_degrees

    pure function filter_signal(b, a, x, init_x, init_y) result(y)
        real(dp), intent(in) :: b(:) !! Feed-forward coefficients; b(1) acts on the current input sample.
        real(dp), intent(in) :: a(:) !! Feedback coefficients; a(1) must be nonzero.
        real(dp), intent(in) :: x(:) !! Input time series.
        real(dp), intent(in), optional :: init_x(:) !! Prior inputs ordered newest to oldest; absent values are zero.
        real(dp), intent(in), optional :: init_y(:) !! Prior outputs ordered newest to oldest; absent values are zero.
        real(dp), allocatable :: y(:)
        real(dp), allocatable :: xb(:)
        real(dp), allocatable :: yb(:)
        real(dp) :: total
        integer :: i
        integer :: j
        integer :: nb
        integer :: na

        allocate(y(size(x)))
        y = 0.0_dp
        if (size(a) == 0 .or. abs(a(1)) <= tiny(1.0_dp)) return
        nb = size(b)
        na = size(a)
        allocate(xb(max(0, nb - 1)), yb(max(0, na - 1)))
        xb = 0.0_dp
        yb = 0.0_dp
        if (present(init_x)) then
            do j = 1, min(size(xb), size(init_x))
                xb(j) = init_x(j)
            end do
        end if
        if (present(init_y)) then
            do j = 1, min(size(yb), size(init_y))
                yb(j) = init_y(j)
            end do
        end if
        do i = 1, size(x)
            total = 0.0_dp
            if (nb > 0) total = b(1) * x(i)
            do j = 2, nb
                if (i - j + 1 >= 1) then
                    total = total + b(j) * x(i - j + 1)
                else
                    total = total + b(j) * xb(j - i)
                end if
            end do
            do j = 2, na
                if (i - j + 1 >= 1) then
                    total = total - a(j) * y(i - j + 1)
                else
                    total = total - a(j) * yb(j - i)
                end if
            end do
            y(i) = total / a(1)
        end do
    end function filter_signal

    pure function fft_filter_signal(b, x, n) result(y)
        real(dp), intent(in) :: b(:) !! FIR coefficients.
        real(dp), intent(in) :: x(:) !! Input time series.
        integer, intent(in), optional :: n !! Requested FFT length; accepted for API parity but direct convolution is used.
        real(dp), allocatable :: y(:)
        real(dp), allocatable :: full(:)

        if (present(n)) continue
        full = convolve(b, x)
        allocate(y(size(x)))
        if (size(x) > 0) y = full(1:size(x))
    end function fft_filter_signal

    pure function zero_phase_filter(b, a, x) result(y)
        real(dp), intent(in) :: b(:) !! Feed-forward filter coefficients.
        real(dp), intent(in) :: a(:) !! Feedback filter coefficients.
        real(dp), intent(in) :: x(:) !! Input signal to filter forward and backward.
        real(dp), allocatable :: y(:)
        real(dp), allocatable :: padded(:)
        real(dp), allocatable :: forward(:)
        real(dp), allocatable :: backward(:)
        integer :: npad
        integer :: i
        integer :: n

        n = size(x)
        npad = 2 * max(size(a), size(b))
        allocate(padded(n + npad))
        padded = 0.0_dp
        if (n > 0) padded(1:n) = x
        forward = filter_signal(b, a, padded)
        backward = filter_signal(b, a, forward(size(forward):1:-1))
        allocate(y(n))
        do i = 1, n
            y(i) = backward(size(backward) - i + 1)
        end do
    end function zero_phase_filter

    pure function median_filter_signal(x, n) result(y)
        real(dp), intent(in) :: x(:) !! Input sequence for running-median smoothing.
        integer, intent(in), optional :: n !! Odd window width; defaults to three.
        real(dp), allocatable :: y(:)
        real(dp), allocatable :: work(:)
        integer :: width
        integer :: half
        integer :: i
        integer :: lo
        integer :: hi

        width = 3
        if (present(n)) width = max(1, n)
        if (modulo(width, 2) == 0) width = width + 1
        half = width / 2
        allocate(y(size(x)))
        do i = 1, size(x)
            lo = max(1, i - half)
            hi = min(size(x), i + half)
            work = x(lo:hi)
            call insertion_sort(work)
            y(i) = work((size(work) + 1) / 2)
        end do
    end function median_filter_signal

    pure function spencer_smooth(x) result(y)
        real(dp), intent(in) :: x(:) !! Input sequence; the seven samples at each end are returned as NaN.
        real(dp), allocatable :: y(:)
        real(dp), allocatable :: raw(:)
        real(dp), allocatable :: b(:)
        real(dp) :: nan_value
        integer :: n

        b = real([-3, -6, -5, 3, 21, 46, 67, 74, 67, 46, 21, 3, -5, -6, -3], dp) / 320.0_dp
        raw = fft_filter_signal(b, x)
        n = size(x)
        allocate(y(n))
        nan_value = ieee_value(0.0_dp, ieee_quiet_nan)
        y = nan_value
        if (n > 14) y(8:n - 7) = raw(15:n)
    end function spencer_smooth

    pure subroutine levinson_durbin(acf, p, a, v, reflection)
        real(dp), intent(in) :: acf(:) !! Autocovariance sequence beginning at lag zero.
        integer, intent(in), optional :: p !! Requested recursion order; defaults to size(acf)-1.
        real(dp), allocatable, intent(out) :: a(:) !! AR polynomial coefficients with leading coefficient one.
        real(dp), intent(out) :: v !! Final prediction-error variance.
        real(dp), allocatable, intent(out) :: reflection(:) !! Reflection coefficients for lags one through p.
        real(dp), allocatable :: current(:)
        real(dp), allocatable :: next(:)
        real(dp) :: g
        integer :: i
        integer :: j
        integer :: order

        order = max(0, size(acf) - 1)
        if (present(p)) order = min(order, max(0, p))
        allocate(reflection(order), a(order + 1))
        a = 0.0_dp
        a(1) = 1.0_dp
        if (order == 0 .or. size(acf) == 0 .or. abs(acf(1)) <= tiny(1.0_dp)) then
            v = 0.0_dp
            return
        end if
        allocate(current(order))
        current = 0.0_dp
        g = -acf(2) / acf(1)
        current(1) = g
        reflection(1) = g
        v = (1.0_dp - g * g) * acf(1)
        do i = 2, order
            g = acf(i + 1)
            do j = 1, i - 1
                g = g + current(j) * acf(i - j + 1)
            end do
            g = -g / v
            allocate(next(order))
            next = current
            do j = 1, i - 1
                next(j) = current(j) + g * current(i - j)
            end do
            next(i) = g
            call move_alloc(next, current)
            v = v * (1.0_dp - g * g)
            reflection(i) = g
        end do
        a(2:) = current(1:order)
    end subroutine levinson_durbin

    pure function unwrap_phase(a, tol) result(y)
        real(dp), intent(in) :: a(:) !! Phase sequence in radians.
        real(dp), intent(in), optional :: tol !! Absolute jump threshold; defaults to pi.
        real(dp), allocatable :: y(:)
        real(dp) :: threshold
        real(dp) :: correction
        real(dp) :: delta
        integer :: i

        threshold = signal_pi
        if (present(tol)) threshold = abs(tol)
        allocate(y(size(a)))
        if (size(a) == 0) return
        y(1) = a(1)
        correction = 0.0_dp
        do i = 2, size(a)
            delta = a(i - 1) - a(i)
            if (delta > threshold) correction = correction + 2.0_dp * signal_pi
            if (delta < -threshold) correction = correction - 2.0_dp * signal_pi
            y(i) = a(i) + correction
        end do
    end function unwrap_phase

    pure function fractional_difference(x, d) result(y)
        real(dp), intent(in) :: x(:) !! Input vector to fractionally difference with the lag operator.
        real(dp), intent(in) :: d !! Differencing exponent; must be greater than -1, with integer parts applied first.
        real(dp), allocatable :: y(:)
        real(dp), allocatable :: work(:)
        real(dp), allocatable :: weights(:)
        real(dp) :: fractional_order
        integer :: integer_order
        integer :: k
        integer :: n

        if (size(x) < 2 .or. d <= -1.0_dp) then
            allocate(y(0))
            return
        end if

        work = x
        integer_order = 0
        if (d >= 1.0_dp) integer_order = int(floor(d))
        do k = 1, integer_order
            if (size(work) < 2) then
                allocate(y(0))
                return
            end if
            work = work(2:) - work(:size(work) - 1)
        end do

        fractional_order = modulo(d, 1.0_dp)
        if (d < 0.0_dp) fractional_order = fractional_order - 1.0_dp
        if (abs(fractional_order) <= epsilon(1.0_dp)) then
            y = work
            return
        end if

        allocate(weights(101))
        weights(1) = 1.0_dp
        do n = 1, 100
            weights(n + 1) = weights(n) * (real(n - 1, dp) - fractional_order) / real(n, dp)
        end do
        y = fft_filter_signal(weights, work)
    end function fractional_difference

    pure subroutine insertion_sort(x)
        real(dp), intent(inout) :: x(:) !! Real vector sorted in ascending order in place.
        real(dp) :: key
        integer :: i
        integer :: j

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

end module signal_filters
