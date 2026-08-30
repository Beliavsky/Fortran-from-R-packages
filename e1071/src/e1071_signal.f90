module e1071_signal
    use e1071_kinds, only: dp
    use e1071_constants, only: e1071_pi
    use e1071_utils, only: hanning_window, hamming_window, rectangle_window
    implicit none
    private

    type, public :: stft_result
        real(dp), allocatable :: values(:, :)
        integer :: window_size = 0
        integer :: increment = 0
        character(len=16) :: window_type = "hanning"
    end type stft_result

    public :: stft

contains

    subroutine stft(x, result, win, increment, coef, window_type)
        real(dp), intent(in) :: x(:) !! Input time series sampled at a regular interval.
        type(stft_result), intent(out) :: result !! Magnitude short-time Fourier transform and effective window metadata.
        integer, intent(in), optional :: win !! Window length; defaults to min(80,floor(n/10)) and is capped at 2*coef.
        integer, intent(in), optional :: increment !! Shift between windows; defaults to min(24,floor(n/30)) and must be positive.
        integer, intent(in), optional :: coef !! Number of nonnegative-frequency coefficients retained; defaults to 64.
        character(len=*), intent(in), optional :: window_type !! Window name: "hanning", "hamming", or "rectangle"; default hanning.
        integer :: n
        integer :: use_win
        integer :: use_inc
        integer :: use_coef
        integer :: nfft
        integer :: nwin
        integer :: i
        integer :: start
        real(dp), allocatable :: weights(:)
        complex(dp), allocatable :: z(:)
        complex(dp), allocatable :: y(:)
        character(len=:), allocatable :: use_type

        n = size(x)
        if (n < 1) error stop "stft: x must not be empty"
        use_coef = 64
        if (present(coef)) use_coef = coef
        if (use_coef < 1) error stop "stft: coef must be positive"
        nfft = 2 * use_coef
        use_win = min(80, n / 10)
        if (present(win)) use_win = win
        if (use_win < 1) use_win = min(n, 1)
        if (use_win > nfft) use_win = nfft
        if (use_win > n) use_win = n
        use_inc = min(24, n / 30)
        if (present(increment)) use_inc = increment
        if (use_inc < 1) use_inc = 1
        nwin = (n - use_win) / use_inc

        use_type = "hanning"
        if (present(window_type)) use_type = trim(adjustl(window_type))
        allocate(weights(use_win))
        select case (use_type)
        case ("hanning", "hanning.window")
            weights(:) = hanning_window(use_win)
            result%window_type = "hanning"
        case ("hamming", "hamming.window")
            weights(:) = hamming_window(use_win)
            result%window_type = "hamming"
        case ("rectangle", "rectangle.window")
            weights(:) = rectangle_window(use_win)
            result%window_type = "rectangle"
        case default
            error stop "stft: unsupported window type"
        end select

        allocate(result%values(nwin + 1, use_coef), z(nfft), y(nfft))
        do i = 0, nwin
            z = cmplx(0.0_dp, 0.0_dp, kind=dp)
            start = 1 + i * use_inc
            z(1:use_win) = cmplx(x(start:start + use_win - 1) * weights, 0.0_dp, kind=dp)
            call fft_any(z, y)
            result%values(i + 1, :) = abs(y(1:use_coef))
        end do
        result%window_size = use_win
        result%increment = use_inc
    end subroutine stft

    subroutine fft_any(x, y)
        complex(dp), intent(in) :: x(:) !! Complex input sequence to transform.
        complex(dp), intent(out) :: y(:) !! Discrete Fourier transform; must have the same length as x.
        integer :: n

        if (size(y) /= size(x)) error stop "fft_any: size mismatch"
        n = size(x)
        if (is_power_of_two(n)) then
            y = x
            call fft_radix2_inplace(y)
        else
            call dft_direct(x, y)
        end if
    end subroutine fft_any

    subroutine fft_radix2_inplace(a)
        complex(dp), intent(inout) :: a(:) !! Power-of-two complex sequence transformed in place with the forward FFT convention.
        integer :: n
        integer :: i
        integer :: j
        integer :: bit
        integer :: len
        integer :: k
        complex(dp) :: temp
        complex(dp) :: w
        complex(dp) :: wlen

        n = size(a)
        j = 1
        do i = 2, n
            bit = n / 2
            do while (j > bit)
                j = j - bit
                bit = bit / 2
                if (bit == 0) exit
            end do
            j = j + bit
            if (i < j) then
                temp = a(i)
                a(i) = a(j)
                a(j) = temp
            end if
        end do

        len = 2
        do while (len <= n)
            wlen = exp(cmplx(0.0_dp, -2.0_dp * e1071_pi / real(len, dp), kind=dp))
            do i = 1, n, len
                w = cmplx(1.0_dp, 0.0_dp, kind=dp)
                do k = 0, len / 2 - 1
                    temp = w * a(i + k + len / 2)
                    a(i + k + len / 2) = a(i + k) - temp
                    a(i + k) = a(i + k) + temp
                    w = w * wlen
                end do
            end do
            len = len * 2
        end do
    end subroutine fft_radix2_inplace

    subroutine dft_direct(x, y)
        complex(dp), intent(in) :: x(:) !! Complex input sequence for a direct O(n^2) Fourier transform.
        complex(dp), intent(out) :: y(:) !! Forward discrete Fourier transform with the same length as x.
        integer :: n
        integer :: k
        integer :: j
        complex(dp) :: phase

        n = size(x)
        y = cmplx(0.0_dp, 0.0_dp, kind=dp)
        do k = 1, n
            do j = 1, n
                phase = exp(cmplx(0.0_dp, -2.0_dp * e1071_pi * real((k - 1) * (j - 1), dp) / real(n, dp), kind=dp))
                y(k) = y(k) + x(j) * phase
            end do
        end do
    end subroutine dft_direct

    pure function is_power_of_two(n) result(ok)
        integer, intent(in) :: n !! Positive integer tested for exact power-of-two structure.
        logical :: ok
        integer :: value

        if (n < 1) then
            ok = .false.
            return
        end if
        value = n
        do while (modulo(value, 2) == 0)
            value = value / 2
        end do
        ok = value == 1
    end function is_power_of_two

end module e1071_signal
