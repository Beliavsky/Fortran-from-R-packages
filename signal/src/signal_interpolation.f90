! SPDX-License-Identifier: GPL-2.0-only
!
! Interpolation and sample-rate conversion translated from CRAN signal.
module signal_interpolation
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
    use signal_kinds, only : dp, signal_pi
    use signal_filters, only : fft_filter_signal, zero_phase_filter
    use signal_fir_design, only : fir1_filter
    use signal_iir_design, only : cheby1_filter
    use signal_types, only : arma_filter
    use signal_utils, only : sinc_value
    implicit none
    private

    public :: interp1_signal
    public :: pchip_interpolate
    public :: interpolate_signal
    public :: resample_signal
    public :: decimate_signal

contains

    pure function pchip_interpolate(x, y, xi) result(yi)
        real(dp), intent(in) :: x(:) !! Strictly monotone interpolation abscissae.
        real(dp), intent(in) :: y(:) !! Function values corresponding one-for-one with x.
        real(dp), intent(in) :: xi(:) !! Query abscissae; extrapolation uses the end cubic pieces.
        real(dp), allocatable :: yi(:)
        real(dp), allocatable :: xx(:)
        real(dp), allocatable :: yy(:)
        real(dp), allocatable :: h(:)
        real(dp), allocatable :: delta(:)
        real(dp), allocatable :: d(:)
        real(dp) :: t
        real(dp) :: hseg
        integer :: i
        integer :: idx
        integer :: n

        n = size(x)
        allocate(yi(size(xi)))
        if (n == 0 .or. size(y) /= n) then
            yi = 0.0_dp
            return
        end if
        if (n == 1) then
            yi = y(1)
            return
        end if
        allocate(xx(n), yy(n))
        if (x(n) < x(1)) then
            xx = x(n:1:-1)
            yy = y(n:1:-1)
        else
            xx = x
            yy = y
        end if
        allocate(h(n - 1), delta(n - 1), d(n))
        h = xx(2:) - xx(:n - 1)
        delta = (yy(2:) - yy(:n - 1)) / h
        call pchip_slopes(h, delta, d)
        do i = 1, size(xi)
            idx = locate_interval(xx, xi(i))
            hseg = xx(idx + 1) - xx(idx)
            t = (xi(i) - xx(idx)) / hseg
            yi(i) = (2.0_dp * t ** 3 - 3.0_dp * t ** 2 + 1.0_dp) * yy(idx) + &
                (t ** 3 - 2.0_dp * t ** 2 + t) * hseg * d(idx) + &
                (-2.0_dp * t ** 3 + 3.0_dp * t ** 2) * yy(idx + 1) + &
                (t ** 3 - t ** 2) * hseg * d(idx + 1)
        end do
    end function pchip_interpolate

    pure function interp1_signal(x, y, xi, method, extrapolate, extrap_value) result(yi)
        real(dp), intent(in) :: x(:) !! Monotone table abscissae.
        real(dp), intent(in) :: y(:) !! Table ordinates corresponding to x.
        real(dp), intent(in) :: xi(:) !! Query abscissae.
        character(len=*), intent(in), optional :: method !! linear, nearest, pchip, cubic, or spline; default is linear.
        logical, intent(in), optional :: extrapolate !! True extrapolates outside x; default false.
        real(dp), intent(in), optional :: extrap_value !! Fill value outside x when extrapolate is false; default is quiet NaN.
        real(dp), allocatable :: yi(:)
        real(dp), allocatable :: pchip_values(:)
        real(dp) :: fill_value
        real(dp) :: value
        character(len=12) :: kind
        logical :: do_extrapolate
        integer :: i
        integer :: idx

        kind = 'linear'
        if (present(method)) kind = adjustl(method)
        do_extrapolate = .false.
        if (present(extrapolate)) do_extrapolate = extrapolate
        fill_value = ieee_value(0.0_dp, ieee_quiet_nan)
        if (present(extrap_value)) fill_value = extrap_value
        allocate(yi(size(xi)))
        yi = fill_value
        if (size(x) < 2 .or. size(y) /= size(x)) return
        if (kind == 'pchip' .or. kind == 'spline') pchip_values = pchip_interpolate(x, y, xi)
        do i = 1, size(xi)
            if (.not. do_extrapolate) then
                if (xi(i) < min(x(1), x(size(x))) .or. xi(i) > max(x(1), x(size(x)))) cycle
            end if
            select case (trim(kind))
            case ('nearest')
                idx = nearest_index(x, xi(i))
                yi(i) = y(idx)
            case ('pchip', 'spline')
                yi(i) = pchip_values(i)
            case ('cubic')
                yi(i) = cubic_lagrange(x, y, xi(i))
            case default
                idx = locate_interval_monotone(x, xi(i))
                if (abs(x(idx + 1) - x(idx)) <= tiny(1.0_dp)) then
                    yi(i) = y(idx)
                else
                    value = (xi(i) - x(idx)) / (x(idx + 1) - x(idx))
                    yi(i) = y(idx) + value * (y(idx + 1) - y(idx))
                end if
            end select
        end do
    end function interp1_signal

    pure function interpolate_signal(x, q, n, wc) result(y)
        real(dp), intent(in) :: x(:) !! Input samples to interpolate by an integer factor.
        integer, intent(in) :: q !! Positive integer interpolation factor.
        integer, intent(in), optional :: n !! Half-length control for the anti-imaging FIR; default four.
        real(dp), intent(in), optional :: wc !! Relative cutoff control before division by q; default 0.5.
        real(dp), allocatable :: y(:)
        real(dp), allocatable :: padded(:)
        real(dp), allocatable :: filtered(:)
        type(arma_filter) :: fir
        integer :: order_control
        integer :: delay
        integer :: i
        real(dp) :: cutoff

        if (q <= 0) then
            allocate(y(0))
            return
        end if
        order_control = 4
        if (present(n)) order_control = max(1, n)
        cutoff = 0.5_dp
        if (present(wc)) cutoff = wc
        allocate(padded(size(x) * q + q * order_control + 1))
        padded = 0.0_dp
        do i = 1, size(x)
            padded(1 + (i - 1) * q) = x(i)
        end do
        fir = fir1_filter(2 * q * order_control + 1, [cutoff / real(q, dp)], 'low')
        filtered = real(q, dp) * fft_filter_signal(fir%b, padded)
        delay = q * order_control + 1
        allocate(y(size(x) * q))
        if (size(y) > 0) y = filtered(delay + 1:delay + size(y))
    end function interpolate_signal

    pure function resample_signal(x, p, q, d) result(y)
        real(dp), intent(in) :: x(:) !! Input samples to resample.
        real(dp), intent(in) :: p !! Positive output-rate numerator; non-integer ratios are allowed.
        real(dp), intent(in), optional :: q !! Positive output-rate denominator; defaults to one.
        integer, intent(in), optional :: d !! Number of neighbouring samples on each side; defaults to five.
        real(dp), allocatable :: y(:)
        real(dp), allocatable :: work(:)
        real(dp), allocatable :: padded(:)
        real(dp), allocatable :: tfrac(:)
        integer, allocatable :: idx(:)
        type(arma_filter) :: antialias
        real(dp) :: denominator
        integer :: order
        integer :: nout
        integer :: j
        integer :: i
        real(dp) :: ratio
        real(dp) :: t
        real(dp) :: limit
        real(dp) :: weight_value

        denominator = 1.0_dp
        if (present(q)) denominator = q
        order = 5
        if (present(d)) order = max(1, d)
        if (p <= 0 .or. denominator <= 0) then
            allocate(y(0))
            return
        end if
        ratio = p / denominator
        work = x
        if (ratio < 1.0_dp) then
            antialias = fir1_filter(2 * order + 1, [ratio], 'low')
            work = fft_filter_signal(antialias%b, work)
        end if
        limit = real(size(work), dp) + 1.0_dp - 1.0_dp / ratio
        nout = max(0, int(floor((limit - 1.0_dp) * ratio + &
            64.0_dp * epsilon(1.0_dp) * max(1.0_dp, abs((limit - 1.0_dp) * ratio)))) + 1)
        allocate(idx(nout), tfrac(nout), y(nout))
        do j = 1, nout
            t = 1.0_dp + real(j - 1, dp) / ratio
            idx(j) = int(t)
            tfrac(j) = t - real(idx(j), dp)
        end do
        allocate(padded(size(work) + 2 * order))
        padded = 0.0_dp
        padded(order + 1:order + size(work)) = work
        y = 0.0_dp
        do i = -order, order
            do j = 1, nout
                weight_value = sinc_value(tfrac(j) - real(i, dp)) * &
                    (0.5_dp + 0.5_dp * cos(signal_pi * (tfrac(j) - real(i, dp)) / &
                    (real(order, dp) + 0.5_dp)))
                y(j) = y(j) + padded(idx(j) + i + order) * weight_value
            end do
        end do
    end function resample_signal

    pure function decimate_signal(x, q, n, ftype) result(y)
        real(dp), intent(in) :: x(:) !! Input sequence to low-pass filter and downsample.
        integer, intent(in) :: q !! Positive integer decimation factor.
        integer, intent(in), optional :: n !! Filter order; defaults to eight for IIR and thirty for FIR.
        character(len=*), intent(in), optional :: ftype !! 'iir' or 'fir'; default is iir.
        real(dp), allocatable :: y(:)
        real(dp), allocatable :: filtered(:)
        type(arma_filter) :: filt
        character(len=8) :: kind
        integer :: order
        integer :: nout
        integer :: i

        if (q <= 0) then
            allocate(y(0))
            return
        end if
        kind = 'iir'
        if (present(ftype)) kind = adjustl(ftype)
        if (kind == 'fir') then
            order = 30
            if (present(n)) order = n
            filt = fir1_filter(order, [1.0_dp / real(q, dp)], 'low')
            filtered = fft_filter_signal(filt%b, x)
        else
            order = 8
            if (present(n)) order = n
            filt = cheby1_filter(order, 0.05_dp, [0.8_dp / real(q, dp)], 'low')
            filtered = zero_phase_filter(filt%b, filt%a, x)
        end if
        nout = (size(x) + q - 1) / q
        allocate(y(nout))
        do i = 1, nout
            y(i) = filtered(1 + (i - 1) * q)
        end do
    end function decimate_signal

    pure subroutine pchip_slopes(h, delta, d)
        real(dp), intent(in) :: h(:) !! Positive interval widths.
        real(dp), intent(in) :: delta(:) !! First divided differences for each interval.
        real(dp), intent(out) :: d(:) !! Shape-preserving endpoint and interior derivatives.
        real(dp) :: w1
        real(dp) :: w2
        integer :: i
        integer :: n

        n = size(d)
        if (n == 2) then
            d = delta(1)
            return
        end if
        d(1) = ((2.0_dp * h(1) + h(2)) * delta(1) - h(1) * delta(2)) / (h(1) + h(2))
        if (d(1) * delta(1) <= 0.0_dp) d(1) = 0.0_dp
        if (delta(1) * delta(2) < 0.0_dp .and. &
            abs(d(1)) > abs(3.0_dp * delta(1))) d(1) = 3.0_dp * delta(1)
        do i = 2, n - 1
            if (delta(i - 1) * delta(i) <= 0.0_dp) then
                d(i) = 0.0_dp
            else
                w1 = 2.0_dp * h(i) + h(i - 1)
                w2 = h(i) + 2.0_dp * h(i - 1)
                d(i) = (w1 + w2) / (w1 / delta(i - 1) + w2 / delta(i))
            end if
        end do
        d(n) = ((2.0_dp * h(n - 1) + h(n - 2)) * delta(n - 1) - &
            h(n - 1) * delta(n - 2)) / (h(n - 1) + h(n - 2))
        if (d(n) * delta(n - 1) <= 0.0_dp) d(n) = 0.0_dp
        if (delta(n - 1) * delta(n - 2) < 0.0_dp .and. &
            abs(d(n)) > abs(3.0_dp * delta(n - 1))) d(n) = 3.0_dp * delta(n - 1)
    end subroutine pchip_slopes

    pure integer function locate_interval(x, value) result(idx)
        real(dp), intent(in) :: x(:) !! Strictly increasing abscissae.
        real(dp), intent(in) :: value !! Query coordinate.
        integer :: lo
        integer :: hi
        integer :: mid

        if (value <= x(1)) then
            idx = 1
            return
        end if
        if (value >= x(size(x))) then
            idx = size(x) - 1
            return
        end if
        lo = 1
        hi = size(x)
        do while (hi - lo > 1)
            mid = (lo + hi) / 2
            if (value >= x(mid)) then
                lo = mid
            else
                hi = mid
            end if
        end do
        idx = lo
    end function locate_interval

    pure integer function locate_interval_monotone(x, value) result(idx)
        real(dp), intent(in) :: x(:) !! Strictly monotone abscissae, increasing or decreasing.
        real(dp), intent(in) :: value !! Query coordinate.
        real(dp), allocatable :: xr(:)
        integer :: reverse_idx

        if (x(size(x)) >= x(1)) then
            idx = locate_interval(x, value)
        else
            xr = x(size(x):1:-1)
            reverse_idx = locate_interval(xr, value)
            idx = size(x) - reverse_idx
        end if
    end function locate_interval_monotone

    pure integer function nearest_index(x, value) result(idx)
        real(dp), intent(in) :: x(:) !! Table abscissae.
        real(dp), intent(in) :: value !! Query coordinate.
        integer :: i

        idx = 1
        do i = 2, size(x)
            if (abs(x(i) - value) < abs(x(idx) - value)) idx = i
        end do
    end function nearest_index

    pure function cubic_lagrange(x, y, value) result(result_value)
        real(dp), intent(in) :: x(:) !! Table abscissae with at least four points.
        real(dp), intent(in) :: y(:) !! Table ordinates corresponding to x.
        real(dp), intent(in) :: value !! Query coordinate.
        real(dp) :: result_value
        integer :: center
        integer :: first
        integer :: last
        integer :: i
        integer :: j
        real(dp) :: basis

        if (size(x) < 4) then
            center = locate_interval_monotone(x, value)
            result_value = y(center) + (value - x(center)) * (y(center + 1) - y(center)) / &
                (x(center + 1) - x(center))
            return
        end if
        center = locate_interval_monotone(x, value)
        first = max(1, min(size(x) - 3, center - 1))
        last = first + 3
        result_value = 0.0_dp
        do i = first, last
            basis = 1.0_dp
            do j = first, last
                if (j /= i) basis = basis * (value - x(j)) / (x(i) - x(j))
            end do
            result_value = result_value + basis * y(i)
        end do
    end function cubic_lagrange

end module signal_interpolation
