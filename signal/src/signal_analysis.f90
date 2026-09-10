! SPDX-License-Identifier: GPL-2.0-only
!
! Frequency-domain and response-analysis routines translated from CRAN signal.
module signal_analysis
    use signal_kinds, only : dp, signal_pi
    use signal_types, only : frequency_response, analog_response, group_delay_response, &
        impulse_response_data, spectrogram_data
    use signal_utils, only : convolve, polynomial_roots, polyval_real_complex, reverse_real, dft_real
    use signal_filters, only : filter_signal, fft_filter_signal
    use signal_windows, only : hanning_window
    implicit none
    private

    public :: chirp_signal
    public :: analog_frequency_response
    public :: digital_frequency_response
    public :: digital_frequency_response_at
    public :: group_delay
    public :: impulse_response
    public :: spectrogram

contains

    pure function chirp_signal(t, f0, t1, f1, form, phase_degrees) result(y)
        real(dp), intent(in) :: t(:) !! Time coordinates in seconds or any units consistent with the frequency inputs.
        real(dp), intent(in), optional :: f0 !! Starting frequency; defaults to zero.
        real(dp), intent(in), optional :: t1 !! Reference time at which f1 is specified; defaults to one.
        real(dp), intent(in), optional :: f1 !! Target frequency at t1; defaults to 100.
        character(len=*), intent(in), optional :: form !! linear, quadratic, or logarithmic; defaults to linear.
        real(dp), intent(in), optional :: phase_degrees !! Initial phase in degrees; defaults to zero.
        real(dp), allocatable :: y(:)
        real(dp) :: start_frequency
        real(dp) :: reference_time
        real(dp) :: target_frequency
        real(dp) :: phase
        real(dp) :: a
        real(dp) :: b
        real(dp) :: base
        character(len=16) :: kind

        start_frequency = 0.0_dp
        if (present(f0)) start_frequency = f0
        reference_time = 1.0_dp
        if (present(t1)) reference_time = t1
        target_frequency = 100.0_dp
        if (present(f1)) target_frequency = f1
        phase = 0.0_dp
        if (present(phase_degrees)) phase = 2.0_dp * signal_pi * phase_degrees / 360.0_dp
        kind = 'linear'
        if (present(form)) kind = adjustl(form)
        allocate(y(size(t)))
        select case (trim(kind))
        case ('quadratic')
            a = (2.0_dp / 3.0_dp) * signal_pi * (target_frequency - start_frequency) / &
                (reference_time * reference_time)
            b = 2.0_dp * signal_pi * start_frequency
            y = cos(a * t ** 3 + b * t + phase)
        case ('logarithmic')
            a = 2.0_dp * signal_pi * reference_time / log(target_frequency - start_frequency)
            b = 2.0_dp * signal_pi * start_frequency
            base = (target_frequency - start_frequency) ** (1.0_dp / reference_time)
            y = cos(a * base ** t + b * t + phase)
        case default
            a = signal_pi * (target_frequency - start_frequency) / reference_time
            b = 2.0_dp * signal_pi * start_frequency
            y = cos(a * t ** 2 + b * t + phase)
        end select
    end function chirp_signal

    pure function analog_frequency_response(b, a, w) result(response)
        real(dp), intent(in) :: b(:) !! Numerator polynomial coefficients in descending powers of s.
        real(dp), intent(in) :: a(:) !! Denominator polynomial coefficients in descending powers of s.
        real(dp), intent(in) :: w(:) !! Angular frequencies in radians per unit time.
        type(analog_response) :: response
        complex(dp), allocatable :: points(:)
        complex(dp), allocatable :: numerator(:)
        complex(dp), allocatable :: denominator(:)

        allocate(points(size(w)))
        points = cmplx(0.0_dp, w, dp)
        numerator = polyval_real_complex(b, points)
        denominator = polyval_real_complex(a, points)
        response%h = numerator / denominator
        response%w = w
    end function analog_frequency_response

    pure function digital_frequency_response(b, a, n, whole, fs) result(response)
        real(dp), intent(in) :: b(:) !! Feed-forward coefficients with b(1) multiplying the current input.
        real(dp), intent(in) :: a(:) !! Feedback polynomial coefficients with a(1) nonzero.
        integer, intent(in), optional :: n !! Number of returned frequency samples; defaults to 512.
        logical, intent(in), optional :: whole !! True spans [0,Fs); false spans [0,Fs/2); default false for real coefficients.
        real(dp), intent(in), optional :: fs !! Frequency scale corresponding to one full digital turn; defaults to 2*pi.
        type(frequency_response) :: response
        real(dp), allocatable :: f(:)
        integer :: count
        integer :: i
        logical :: full_circle
        real(dp) :: sampling

        count = 512
        if (present(n)) count = max(1, n)
        full_circle = .false.
        if (present(whole)) full_circle = whole
        sampling = 2.0_dp * signal_pi
        if (present(fs)) sampling = fs
        allocate(f(count))
        if (full_circle) then
            do i = 1, count
                f(i) = sampling * real(i - 1, dp) / real(count, dp)
            end do
        else
            do i = 1, count
                f(i) = 0.5_dp * sampling * real(i - 1, dp) / real(count, dp)
            end do
        end if
        response = digital_frequency_response_at(b, a, f, sampling)
    end function digital_frequency_response

    pure function digital_frequency_response_at(b, a, f, fs) result(response)
        real(dp), intent(in) :: b(:) !! Feed-forward coefficients with b(1) multiplying the current input.
        real(dp), intent(in) :: a(:) !! Feedback coefficients with a(1) nonzero.
        real(dp), intent(in) :: f(:) !! Frequencies in units determined by fs.
        real(dp), intent(in), optional :: fs !! Frequency scale for one full digital turn; defaults to 2*pi.
        type(frequency_response) :: response
        complex(dp), allocatable :: numerator(:)
        complex(dp), allocatable :: denominator(:)
        complex(dp), allocatable :: z(:)
        real(dp) :: sampling
        real(dp) :: omega
        integer :: i
        integer :: j

        sampling = 2.0_dp * signal_pi
        if (present(fs)) sampling = fs
        allocate(z(size(f)), numerator(size(f)), denominator(size(f)))
        numerator = cmplx(0.0_dp, 0.0_dp, dp)
        denominator = cmplx(0.0_dp, 0.0_dp, dp)
        do i = 1, size(f)
            omega = 2.0_dp * signal_pi * f(i) / sampling
            z(i) = exp(cmplx(0.0_dp, -omega, dp))
            do j = size(b), 1, -1
                numerator(i) = numerator(i) * z(i) + b(j)
            end do
            do j = size(a), 1, -1
                denominator(i) = denominator(i) * z(i) + a(j)
            end do
        end do
        response%h = numerator / denominator
        response%f = f
    end function digital_frequency_response_at

    pure function group_delay(b, a, n, whole, fs) result(response)
        real(dp), intent(in) :: b(:) !! Numerator coefficients of the digital filter.
        real(dp), intent(in) :: a(:) !! Denominator coefficients of the digital filter.
        integer, intent(in), optional :: n !! Number of returned points; defaults to 512.
        logical, intent(in), optional :: whole !! True returns a full-circle response; default false.
        real(dp), intent(in), optional :: fs !! Optional sampling-frequency scale; absent means radians/sample.
        type(group_delay_response) :: response
        real(dp), allocatable :: c(:)
        real(dp), allocatable :: cr(:)
        complex(dp) :: num
        complex(dp) :: den
        complex(dp) :: phase
        real(dp) :: sampling
        real(dp) :: omega
        integer :: count
        integer :: nfft
        integer :: oa
        integer :: i
        integer :: j
        logical :: full_circle
        logical :: use_hz

        count = 512
        if (present(n)) count = max(1, n)
        full_circle = .false.
        if (present(whole)) full_circle = whole
        use_hz = present(fs)
        sampling = 1.0_dp
        if (present(fs)) sampling = fs
        nfft = count
        if (.not. full_circle) nfft = 2 * count
        oa = max(0, size(a) - 1)
        c = convolve(b, reverse_real(a))
        allocate(cr(size(c)))
        do i = 1, size(c)
            cr(i) = c(i) * real(i - 1, dp)
        end do
        allocate(response%gd(count), response%w(count))
        do i = 1, count
            if (use_hz) then
                response%w(i) = sampling * real(i - 1, dp) / real(nfft, dp)
                omega = 2.0_dp * signal_pi * real(i - 1, dp) / real(nfft, dp)
            else
                response%w(i) = 2.0_dp * signal_pi * real(i - 1, dp) / real(nfft, dp)
                omega = response%w(i)
            end if
            phase = exp(cmplx(0.0_dp, -omega, dp))
            num = cmplx(0.0_dp, 0.0_dp, dp)
            den = cmplx(0.0_dp, 0.0_dp, dp)
            do j = size(cr), 1, -1
                num = num * phase + cr(j)
                den = den * phase + c(j)
            end do
            if (abs(den) < 2.0_dp * epsilon(1.0_dp)) then
                response%gd(i) = -real(oa, dp)
            else
                response%gd(i) = real(num / den, dp) - real(oa, dp)
            end if
        end do
    end function group_delay

    pure function impulse_response(b, a, n, fs) result(response)
        real(dp), intent(in) :: b(:) !! Numerator filter coefficients.
        real(dp), intent(in) :: a(:) !! Denominator filter coefficients.
        integer, intent(in), optional :: n !! Number of impulse samples; inferred from poles when absent.
        real(dp), intent(in), optional :: fs !! Samples per unit time for the returned time axis; defaults to one.
        type(impulse_response_data) :: response
        complex(dp), allocatable :: roots(:)
        real(dp), allocatable :: impulse(:)
        real(dp) :: sampling
        real(dp) :: maxpole
        real(dp) :: angle_min
        real(dp) :: damped_max
        integer :: count
        integer :: i
        integer :: candidate

        sampling = 1.0_dp
        if (present(fs)) sampling = fs
        if (present(n)) then
            count = max(1, n)
        else if (size(a) > 1) then
            roots = polynomial_roots(a)
            if (size(roots) == 0) then
                count = max(1, size(b))
            else
                maxpole = maxval(abs(roots))
                if (maxpole > 1.0_dp + 1.0e-6_dp) then
                    count = max(1, floor(6.0_dp / log10(maxpole)))
                else if (maxpole < 1.0_dp - 1.0e-6_dp) then
                    count = max(1, floor(-6.0_dp / log10(maxpole)))
                else
                    count = 30
                    angle_min = huge(1.0_dp)
                    damped_max = 0.0_dp
                    do i = 1, size(roots)
                        if (abs(roots(i)) >= 1.0_dp - 1.0e-6_dp .and. abs(atan2(aimag(roots(i)), real(roots(i), dp))) > 0.0_dp) then
                            angle_min = min(angle_min, abs(atan2(aimag(roots(i)), real(roots(i), dp))))
                        end if
                        if (abs(roots(i)) < 1.0_dp - 1.0e-6_dp) damped_max = max(damped_max, abs(roots(i)))
                    end do
                    if (angle_min < huge(1.0_dp)) then
                        candidate = ceiling(10.0_dp * signal_pi / angle_min)
                        count = max(count, candidate)
                    end if
                    if (damped_max > 0.0_dp) then
                        candidate = floor(-3.0_dp / log10(damped_max))
                        count = max(count, candidate)
                    end if
                end if
                count = count + size(b)
            end if
        else
            count = max(1, size(b))
        end if
        allocate(impulse(count))
        impulse = 0.0_dp
        impulse(1) = 1.0_dp
        if (size(a) == 1) then
            response%x = fft_filter_signal(b / a(1), impulse)
        else
            response%x = filter_signal(b, a, impulse)
        end if
        allocate(response%t(count))
        do i = 1, count
            response%t(i) = real(i - 1, dp) / sampling
        end do
    end function impulse_response

    pure function spectrogram(x, n, fs, window, overlap) result(result_data)
        real(dp), intent(in) :: x(:) !! Input signal samples.
        integer, intent(in), optional :: n !! DFT length; defaults to min(256,size(x)).
        real(dp), intent(in), optional :: fs !! Sampling frequency; defaults to two.
        real(dp), intent(in), optional :: window(:) !! Analysis window; defaults to Hanning of the DFT length.
        integer, intent(in), optional :: overlap !! Number of overlapped window samples; defaults to ceiling(window_length/2).
        type(spectrogram_data) :: result_data
        real(dp), allocatable :: win(:)
        real(dp), allocatable :: segment(:)
        complex(dp), allocatable :: transform(:)
        integer, allocatable :: offsets(:)
        real(dp) :: sampling
        integer :: fft_size
        integer :: win_size
        integer :: overlap_count
        integer :: step
        integer :: columns
        integer :: ret_n
        integer :: i
        integer :: j
        integer :: offset

        sampling = 2.0_dp
        if (present(fs)) sampling = fs
        fft_size = min(256, max(1, size(x)))
        if (present(n)) fft_size = max(1, n)
        if (present(window)) then
            win = window
        else
            win = hanning_window(fft_size)
        end if
        win_size = size(win)
        if (win_size > fft_size) fft_size = win_size
        overlap_count = (win_size + 1) / 2
        if (present(overlap)) overlap_count = overlap
        step = max(1, win_size - overlap_count)
        if (size(x) > win_size) then
            columns = 1 + max(0, (size(x) - win_size - 1) / step)
        else
            columns = 1
        end if
        allocate(offsets(columns))
        do j = 1, columns
            offsets(j) = 1 + (j - 1) * step
        end do
        if (modulo(fft_size, 2) == 1) then
            ret_n = (fft_size + 1) / 2
        else
            ret_n = fft_size / 2
        end if
        allocate(result_data%s(ret_n, columns), segment(fft_size))
        result_data%s = cmplx(0.0_dp, 0.0_dp, dp)
        do j = 1, columns
            segment = 0.0_dp
            offset = offsets(j)
            do i = 1, win_size
                if (offset + i - 1 <= size(x)) segment(i) = x(offset + i - 1) * win(i)
            end do
            transform = dft_real(segment)
            result_data%s(:, j) = transform(1:ret_n)
        end do
        allocate(result_data%f(ret_n), result_data%t(columns))
        do i = 1, ret_n
            result_data%f(i) = real(i - 1, dp) * sampling / real(fft_size, dp)
        end do
        do j = 1, columns
            result_data%t(j) = real(offsets(j), dp) / sampling
        end do
    end function spectrogram

end module signal_analysis
