! SPDX-License-Identifier: GPL-2.0-only
!
! FIR design routines translated from CRAN signal / Octave-Forge algorithms.
module signal_fir_design
    use signal_kinds, only : dp, signal_pi
    use signal_types, only : arma_filter, filter_order
    use signal_filters, only : make_filter_order, make_ma
    use signal_utils, only : inverse_dft, polyval_real_complex, solve_linear_system
    use signal_windows, only : hamming_window
    implicit none
    private

    public :: fir1_filter
    public :: fir2_filter
    public :: kaiser_order
    public :: remez_filter

contains

    pure function fir2_filter(n, f, magnitude, grid_n, ramp_n, window) result(filt)
        integer, intent(in) :: n !! FIR order, so the returned filter normally has n+1 taps.
        real(dp), intent(in) :: f(:) !! Nondecreasing normalised frequency knots from zero to one.
        real(dp), intent(in) :: magnitude(:) !! Desired magnitude at each frequency knot.
        integer, intent(in), optional :: grid_n !! Half-grid size used for frequency sampling; defaults to 512.
        real(dp), intent(in), optional :: ramp_n !! Discontinuity-ramp width in grid points; defaults to grid_n/20.
        real(dp), intent(in), optional :: window(:) !! Optional n+1 point time-domain window; defaults to Hamming.
        type(arma_filter) :: filt
        complex(dp), allocatable :: spectrum(:)
        complex(dp), allocatable :: time(:)
        real(dp), allocatable :: grid(:)
        real(dp), allocatable :: b(:)
        real(dp), allocatable :: win(:)
        integer :: g
        integer :: total
        integer :: i
        integer :: j
        integer :: left_count
        integer :: right_count
        real(dp) :: x
        real(dp) :: ramp_width

        g = 512
        if (present(grid_n)) g = max(2, grid_n)
        if (2 * g < n + 1) then
            g = 1
            do while (2 * g < n + 1)
                g = 2 * g
            end do
        end if
        ramp_width = real(g, dp) / 20.0_dp
        if (present(ramp_n)) ramp_width = max(0.0_dp, ramp_n)
        allocate(grid(g + 1))
        do i = 0, g
            x = real(i, dp) / real(g, dp)
            grid(i + 1) = piecewise_magnitude(f, magnitude, x, ramp_width / real(g, dp))
        end do
        if (modulo(n, 2) == 0) then
            total = 2 * g
            allocate(spectrum(total))
            spectrum(1:g + 1) = cmplx(grid, 0.0_dp, dp)
            do i = 1, g - 1
                spectrum(g + 1 + i) = cmplx(grid(g + 1 - i), 0.0_dp, dp)
            end do
            time = inverse_dft(spectrum)
            allocate(b(n + 1))
            left_count = floor(real(n + 1, dp) / 2.0_dp)
            right_count = ceiling(real(n + 1, dp) / 2.0_dp)
            do i = 1, left_count
                b(i) = real(time(total - left_count + i), dp)
            end do
            do i = 1, right_count
                b(left_count + i) = real(time(i), dp)
            end do
        else
            total = 4 * g
            allocate(spectrum(total))
            spectrum = cmplx(0.0_dp, 0.0_dp, dp)
            spectrum(1:g + 1) = cmplx(grid, 0.0_dp, dp)
            do i = 1, g - 1
                spectrum(3 * g + 1 + i) = cmplx(grid(g + 1 - i), 0.0_dp, dp)
            end do
            time = inverse_dft(spectrum)
            allocate(b(n + 1))
            j = 0
            do i = total - n + 1, total, 2
                j = j + 1
                b(j) = 2.0_dp * real(time(i), dp)
            end do
            do i = 2, n + 2, 2
                j = j + 1
                b(j) = 2.0_dp * real(time(i), dp)
            end do
        end if
        if (present(window)) then
            if (size(window) == size(b)) then
                win = window
            else
                win = hamming_window(size(b))
            end if
        else
            win = hamming_window(size(b))
        end if
        b = b * win
        filt = make_ma(b)
    end function fir2_filter

    pure function fir1_filter(n, w, filter_type, window, scale) result(filt)
        integer, intent(in) :: n !! Nominal FIR order; high-pass/band-stop designs may be promoted to the next even order.
        real(dp), intent(in) :: w(:) !! One or more normalised band edges in ascending order.
        character(len=*), intent(in), optional :: filter_type !! low, high, stop, pass, DC-0, or DC-1; inferred low by default.
        real(dp), intent(in), optional :: window(:) !! Optional time-domain window.
        logical, intent(in), optional :: scale !! True normalises gain in the first pass band; defaults to true.
        type(arma_filter) :: filt
        real(dp), allocatable :: f(:)
        real(dp), allocatable :: m(:)
        real(dp), allocatable :: win(:)
        complex(dp), allocatable :: response(:)
        complex(dp) :: point(1)
        real(dp) :: w0
        real(dp) :: renorm
        character(len=8) :: kind
        logical :: first_pass
        logical :: do_scale
        integer :: order
        integer :: bands
        integer :: i

        kind = 'low'
        if (present(filter_type)) kind = adjustl(filter_type)
        first_pass = kind == 'low' .or. kind == 'stop' .or. kind == 'DC-1'
        do_scale = .true.
        if (present(scale)) do_scale = scale
        order = max(0, n)
        bands = size(w) + 1
        allocate(f(2 * bands), m(2 * bands))
        f = 0.0_dp
        m = 0.0_dp
        f(1) = 0.0_dp
        f(2 * bands) = 1.0_dp
        do i = 1, size(w)
            f(2 * i) = w(i)
            f(2 * i + 1) = w(i)
        end do
        do i = 1, bands
            if (modulo(i - merge(1, 0, first_pass), 2) == 0) then
                m(2 * i - 1:2 * i) = 1.0_dp
            else
                m(2 * i - 1:2 * i) = 0.0_dp
            end if
        end do
        if (modulo(order, 2) == 1 .and. m(2 * bands) > 0.5_dp) order = order + 1
        if (present(window)) then
            if (size(window) == order + 1) then
                win = window
            else
                win = hamming_window(order + 1)
            end if
        else
            win = hamming_window(order + 1)
        end if
        filt = fir2_filter(order, f, m, 512, 25.6_dp, win)
        if (do_scale .and. size(filt%b) > 0) then
            if (m(1) > 0.5_dp) then
                w0 = 0.5_dp * (f(2) - f(1))
            else if (size(f) >= 4) then
                w0 = f(3) + 0.5_dp * (f(4) - f(3))
            else
                w0 = 0.5_dp
            end if
            point(1) = exp(cmplx(0.0_dp, -signal_pi * w0, dp))
            response = polyval_real_complex(filt%b, point)
            if (abs(response(1)) > tiny(1.0_dp)) then
                renorm = 1.0_dp / abs(response(1))
                filt%b = renorm * filt%b
            end if
        end if
    end function fir1_filter

    pure function kaiser_order(f, magnitude, deviation, fs) result(spec)
        real(dp), intent(in) :: f(:) !! Transition-edge frequencies, two entries for each transition.
        real(dp), intent(in) :: magnitude(:) !! Target magnitude for each alternating band.
        real(dp), intent(in) :: deviation(:) !! Allowed deviation for one or all bands; the minimum controls the design.
        real(dp), intent(in), optional :: fs !! Sampling frequency; defaults to two.
        type(filter_order) :: spec
        real(dp), allocatable :: wc(:)
        real(dp) :: sampling
        real(dp) :: dev
        real(dp) :: attenuation
        real(dp) :: beta
        real(dp) :: dw
        integer :: n
        integer :: i
        character(len=8) :: kind

        sampling = 2.0_dp
        if (present(fs)) sampling = fs
        allocate(wc(size(f) / 2))
        do i = 1, size(wc)
            wc(i) = 0.5_dp * (f(2 * i - 1) + f(2 * i)) / sampling
        end do
        if (size(wc) == 1) then
            if (magnitude(1) > magnitude(2)) then
                kind = 'low'
            else
                kind = 'high'
            end if
        else if (size(wc) == 2) then
            if (magnitude(1) > magnitude(2)) then
                kind = 'stop'
            else
                kind = 'pass'
            end if
        else
            if (magnitude(1) > magnitude(2)) then
                kind = 'DC-1'
            else
                kind = 'DC-0'
            end if
        end if
        dev = minval(deviation)
        attenuation = -20.0_dp * log10(dev)
        if (attenuation > 50.0_dp) then
            beta = 0.1102_dp * (attenuation - 8.7_dp)
        else if (attenuation >= 21.0_dp) then
            beta = 0.5842_dp * (attenuation - 21.0_dp) ** 0.4_dp + &
                0.07886_dp * (attenuation - 21.0_dp)
        else
            beta = 0.0_dp
        end if
        dw = huge(1.0_dp)
        do i = 1, size(f) / 2
            dw = min(dw, f(2 * i) - f(2 * i - 1))
        end do
        dw = 2.0_dp * signal_pi * dw / sampling
        n = max(1, ceiling((attenuation - 8.0_dp) / (2.285_dp * dw)))
        if ((magnitude(1) > magnitude(2)) .eqv. (modulo(size(wc), 2) == 0)) then
            if (modulo(n, 2) == 1) n = n + 1
        end if
        spec = make_filter_order(n, wc, kind, beta=beta)
    end function kaiser_order

    pure function remez_filter(n, f, amplitude, weight, ftype, density) result(filt)
        integer, intent(in) :: n !! Number of taps minus one in the signal R interface.
        real(dp), intent(in) :: f(:) !! Paired band-edge frequencies normalised from zero to one.
        real(dp), intent(in) :: amplitude(:) !! Desired amplitude values at the band edges.
        real(dp), intent(in), optional :: weight(:) !! Optional positive band weights.
        character(len=*), intent(in), optional :: ftype !! bandpass, differentiator, or hilbert.
        integer, intent(in), optional :: density !! Grid density; accepted for interface compatibility.
        type(arma_filter) :: filt
        real(dp), allocatable :: ata(:, :)
        real(dp), allocatable :: atb(:)
        real(dp), allocatable :: coeff(:)
        real(dp), allocatable :: rhs(:)
        real(dp), allocatable :: row(:)
        real(dp) :: freq
        real(dp) :: desired
        real(dp) :: wt
        real(dp) :: phase_shift
        character(len=16) :: kind
        logical :: ok
        integer :: taps
        integer :: m
        integer :: grid_points
        integer :: i
        integer :: j
        integer :: k

        taps = n + 1
        kind = 'bandpass'
        if (present(ftype)) kind = adjustl(ftype)
        grid_points = max(256, 32 * taps)
        if (present(density)) grid_points = max(grid_points, density * taps)
        m = (taps + 1) / 2
        allocate(ata(m, m), atb(m), row(m))
        ata = 0.0_dp
        atb = 0.0_dp
        phase_shift = 0.5_dp * real(taps - 1, dp)
        do i = 0, grid_points
            freq = real(i, dp) / real(grid_points, dp)
            call desired_response_at(f, amplitude, weight, freq, desired, wt)
            if (wt <= 0.0_dp) cycle
            if (kind == 'differentiator') desired = desired * freq
            do j = 1, m
                row(j) = cos(signal_pi * freq * (real(j - 1, dp) - phase_shift))
            end do
            if (kind == 'hilbert') then
                do j = 1, m
                    row(j) = sin(signal_pi * freq * (real(j - 1, dp) - phase_shift))
                end do
            end if
            do j = 1, m
                atb(j) = atb(j) + wt * wt * row(j) * desired
                do k = 1, m
                    ata(j, k) = ata(j, k) + wt * wt * row(j) * row(k)
                end do
            end do
        end do
        call solve_linear_system(ata, atb, rhs, ok)
        allocate(coeff(taps))
        coeff = 0.0_dp
        if (ok) then
            do j = 1, m
                coeff(j) = rhs(j)
                coeff(taps - j + 1) = rhs(j)
            end do
            if (kind == 'hilbert') then
                do j = 1, taps
                    coeff(j) = coeff(j) * sign(1.0_dp, real((taps + 1) / 2 - j, dp))
                end do
                if (modulo(taps, 2) == 1) coeff((taps + 1) / 2) = 0.0_dp
            end if
        end if
        filt = make_ma(coeff)
    end function remez_filter

    pure function piecewise_magnitude(f, magnitude, x, ramp_width) result(value)
        real(dp), intent(in) :: f(:) !! Frequency knots.
        real(dp), intent(in) :: magnitude(:) !! Magnitudes at frequency knots.
        real(dp), intent(in) :: x !! Normalised frequency at which to interpolate.
        real(dp), intent(in) :: ramp_width !! Half-width used to smooth duplicate-frequency jumps.
        real(dp) :: value
        real(dp) :: left_x
        real(dp) :: right_x
        integer :: i

        if (size(f) == 0 .or. size(f) /= size(magnitude)) then
            value = 0.0_dp
            return
        end if
        if (x <= f(1)) then
            value = magnitude(1)
            return
        end if
        if (x >= f(size(f))) then
            value = magnitude(size(magnitude))
            return
        end if
        do i = 1, size(f) - 1
            if (f(i + 1) > f(i)) then
                if (x >= f(i) .and. x <= f(i + 1)) then
                    value = magnitude(i) + (magnitude(i + 1) - magnitude(i)) * &
                        (x - f(i)) / (f(i + 1) - f(i))
                    return
                end if
            else if (abs(x - f(i)) <= ramp_width .and. ramp_width > 0.0_dp) then
                left_x = f(i) - 0.5_dp * ramp_width
                right_x = f(i) + 0.5_dp * ramp_width
                if (x <= left_x) then
                    value = magnitude(i)
                else if (x >= right_x) then
                    value = magnitude(i + 1)
                else
                    value = magnitude(i) + (magnitude(i + 1) - magnitude(i)) * &
                        (x - left_x) / (right_x - left_x)
                end if
                return
            end if
        end do
        value = magnitude(size(magnitude))
    end function piecewise_magnitude

    pure subroutine desired_response_at(f, amplitude, weight, x, desired, wt)
        real(dp), intent(in) :: f(:) !! Paired band-edge frequencies.
        real(dp), intent(in) :: amplitude(:) !! Desired amplitudes at each edge.
        real(dp), intent(in), optional :: weight(:) !! Optional weight for each frequency band.
        real(dp), intent(in) :: x !! Normalised frequency to evaluate.
        real(dp), intent(out) :: desired !! Interpolated desired response; zero outside specified bands.
        real(dp), intent(out) :: wt !! Band weight; zero outside specified bands.
        integer :: band
        integer :: lo
        integer :: hi

        desired = 0.0_dp
        wt = 0.0_dp
        do band = 1, size(f) / 2
            lo = 2 * band - 1
            hi = 2 * band
            if (x >= f(lo) .and. x <= f(hi)) then
                if (f(hi) > f(lo)) then
                    desired = amplitude(lo) + (amplitude(hi) - amplitude(lo)) * &
                        (x - f(lo)) / (f(hi) - f(lo))
                else
                    desired = amplitude(lo)
                end if
                wt = 1.0_dp
                if (present(weight)) then
                    if (band <= size(weight)) wt = weight(band)
                end if
                return
            end if
        end do
    end subroutine desired_response_at

end module signal_fir_design
