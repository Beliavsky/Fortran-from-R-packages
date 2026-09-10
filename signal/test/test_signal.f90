! SPDX-License-Identifier: GPL-2.0-only
program test_signal
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan
    use signal
    implicit none

    type(arma_filter) :: filt
    type(frequency_response) :: response
    type(group_delay_response) :: gd
    type(sgolay_filter) :: sg
    real(dp), allocatable :: y(:)
    real(dp), allocatable :: w(:)
    real(dp), allocatable :: a(:)
    real(dp), allocatable :: reflection(:)
    complex(dp), allocatable :: z(:)
    real(dp) :: variance
    real(dp) :: flat_one
    integer :: i

    call assert_close_vector(bartlett_window(5), [0.0_dp, 0.5_dp, 1.0_dp, 0.5_dp, 0.0_dp], &
        2.0e-14_dp, 'bartlett')

    flat_one = (1.0_dp - 1.93_dp + 1.29_dp - 0.388_dp + 0.0322_dp) / 4.6402_dp
    call assert_close_vector(flattop_window(1, .true.), [flat_one], 2.0e-14_dp, &
        'flattop periodic n=1')

    call assert_close_vector(convolve([1.0_dp, 2.0_dp], [3.0_dp, 4.0_dp]), &
        [3.0_dp, 10.0_dp, 8.0_dp], 2.0e-14_dp, 'convolution')

    y = filter_signal([1.0_dp], [1.0_dp, -0.5_dp], [1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp])
    call assert_close_vector(y, [1.0_dp, 0.5_dp, 0.25_dp, 0.125_dp], 2.0e-14_dp, &
        'recursive filter')

    filt = butter_filter(2, [0.5_dp], 'low')
    call assert_close_vector(filt%b, [0.292893218813452_dp, 0.585786437626905_dp, &
        0.292893218813452_dp], 2.0e-13_dp, 'butter numerator')
    call assert_close_vector(filt%a, [1.0_dp, 0.0_dp, 0.171572875253810_dp], &
        2.0e-13_dp, 'butter denominator')

    response = digital_frequency_response(filt%b, filt%a, 32)
    call assert_close_scalar(abs(response%h(1)), 1.0_dp, 2.0e-13_dp, 'freqz dc gain')
    call assert_close_scalar(abs(response%h(17)) ** 2, 0.5_dp, 3.0e-13_dp, &
        'freqz half-power cutoff')

    gd = group_delay([0.0_dp, 1.0_dp], [1.0_dp], 4)
    call assert_close_vector(gd%gd, [1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp], 2.0e-13_dp, &
        'group delay one-sample shift')

    gd = group_delay([1.0_dp], [1.0_dp, 0.9_dp], 4)
    call assert_close_vector(gd%gd, [-0.473684210526316_dp, -0.469183776569300_dp, &
        -0.447513812154696_dp, -0.323159679149511_dp], 8.0e-9_dp, 'group delay iir')

    filt = fir1_filter(2, [0.5_dp], 'low', hanning_window(3), .true.)
    call assert_close_vector(filt%b, [0.0_dp, 1.0_dp, 0.0_dp], 3.0e-13_dp, 'fir1')

    y = pchip_interpolate([0.0_dp, 1.0_dp, 2.0_dp], [0.0_dp, 1.0_dp, 4.0_dp], &
        [0.0_dp, 0.5_dp, 1.0_dp, 1.5_dp, 2.0_dp])
    call assert_close_vector(y, [0.0_dp, 0.3125_dp, 1.0_dp, 2.1875_dp, 4.0_dp], &
        2.0e-13_dp, 'pchip')

    y = interp1_signal([0.0_dp, 1.0_dp, 2.0_dp], [0.0_dp, 2.0_dp, 4.0_dp], &
        [-1.0_dp, 0.5_dp, 3.0_dp], 'linear')
    call assert_true(ieee_is_nan(y(1)), 'interp1 lower out-of-range NaN')
    call assert_close_scalar(y(2), 1.0_dp, 2.0e-14_dp, 'interp1 linear')
    call assert_true(ieee_is_nan(y(3)), 'interp1 upper out-of-range NaN')

    y = resample_signal([0.0_dp, 1.0_dp, 0.0_dp], 1.0_dp, 0.05_dp, 2)
    call assert_true(size(y) == 60, 'resample accepts non-integer rate ratio')
    call assert_true(all(ieee_is_finite(y)), 'resample finite values')

    y = resample_signal(sin(2.0_dp * signal_pi * real([(i, i = 0, 10)], dp) / 5.0_dp), &
        1.0_dp, 0.05_dp)
    call assert_true(size(y) == 220, 'resample documented 20x endpoint count')

    sg = savitzky_golay(2, 5)
    y = savitzky_golay_filter([1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 5.0_dp], &
        2, 5, supplied_filter=sg)
    call assert_close_vector(y, [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 5.0_dp], &
        3.0e-12_dp, 'savitzky-golay linear reproduction')

    call levinson_durbin([1.0_dp, 0.5_dp, 0.25_dp], 2, a, variance, reflection)
    call assert_close_vector(a, [1.0_dp, -0.5_dp, 0.0_dp], 3.0e-13_dp, 'levinson coefficients')
    call assert_close_scalar(variance, 0.75_dp, 3.0e-13_dp, 'levinson variance')
    call assert_close_vector(reflection, [-0.5_dp, 0.0_dp], 3.0e-13_dp, &
        'levinson reflection coefficients')

    y = chirp_signal([0.0_dp, 0.25_dp, 0.5_dp], f0=0.0_dp, t1=1.0_dp, &
        f1=1.0_dp, form='linear')
    call assert_close_scalar(y(1), 1.0_dp, 2.0e-14_dp, 'chirp starts at cosine phase zero')

    z = inverse_dft([cmplx(1.0_dp, 0.0_dp, dp), cmplx(1.0_dp, 0.0_dp, dp), &
        cmplx(1.0_dp, 0.0_dp, dp), cmplx(1.0_dp, 0.0_dp, dp)])
    call assert_close_complex_vector(z, [cmplx(1.0_dp, 0.0_dp, dp), &
        cmplx(0.0_dp, 0.0_dp, dp), cmplx(0.0_dp, 0.0_dp, dp), &
        cmplx(0.0_dp, 0.0_dp, dp)], 3.0e-13_dp, 'inverse dft')

    filt = remez_filter(10, [0.0_dp, 0.2_dp, 0.3_dp, 0.5_dp], &
        [1.0_dp, 1.0_dp, 0.0_dp, 0.0_dp])
    call assert_true(size(filt%b) == 11, 'remez coefficient count')
    call assert_true(all(ieee_is_finite(filt%b)), 'remez coefficients finite')
    call assert_close_vector(filt%b, filt%b(size(filt%b):1:-1), 3.0e-11_dp, &
        'remez linear-phase symmetry')

    filt = ellip_filter(5, 3.0_dp, 40.0_dp, [0.1_dp], 'low')
    response = digital_frequency_response(filt%b, filt%a, 16)
    call assert_true(all(ieee_is_finite(real(response%h, dp))), 'elliptic response finite')
    call assert_close_scalar(abs(response%h(1)), 1.0_dp, 5.0e-10_dp, 'elliptic dc gain')

    w = kaiser_window(9, 4.0_dp)
    call assert_close_vector(w, w(size(w):1:-1), 3.0e-13_dp, 'kaiser symmetry')

    y = median_filter_signal([1.0_dp, 20.0_dp, 2.0_dp, 3.0_dp, 4.0_dp], 3)
    call assert_close_scalar(y(2), 2.0_dp, 2.0e-14_dp, 'median filter interior')

    y = unwrap_phase([0.0_dp, 0.9_dp * signal_pi, -0.9_dp * signal_pi])
    call assert_true(abs(y(3) - y(2)) < signal_pi, 'phase unwrap removes 2*pi jump')

    y = fractional_difference([1.0_dp, 2.0_dp, 3.0_dp], 1.0_dp)
    call assert_close_vector(y, [1.0_dp, 1.0_dp], 2.0e-14_dp, 'integer fractional difference')
    y = fractional_difference([1.0_dp, 2.0_dp, 3.0_dp], -0.5_dp)
    call assert_close_vector(y, [1.0_dp, 2.5_dp, 4.375_dp], 3.0e-13_dp, &
        'negative fractional difference')

    print '(a)', 'signal deterministic tests passed'

contains

    subroutine assert_true(condition, label)
        logical, intent(in) :: condition !! Condition that must evaluate true for the test to pass.
        character(len=*), intent(in) :: label !! Human-readable test label printed on failure.

        if (.not. condition) then
            write (*, '(a,1x,a)') 'FAIL:', trim(label)
            error stop 1
        end if
    end subroutine assert_true

    subroutine assert_close_scalar(actual, expected, tolerance, label)
        real(dp), intent(in) :: actual !! Computed scalar value under test.
        real(dp), intent(in) :: expected !! Reference scalar value.
        real(dp), intent(in) :: tolerance !! Maximum allowed absolute error.
        character(len=*), intent(in) :: label !! Human-readable test label printed on failure.

        if (abs(actual - expected) > tolerance) then
            write (*, '(a,1x,a,2(1x,es24.16))') 'FAIL:', trim(label), actual, expected
            error stop 1
        end if
    end subroutine assert_close_scalar

    subroutine assert_close_vector(actual, expected, tolerance, label)
        real(dp), intent(in) :: actual(:) !! Computed vector under test.
        real(dp), intent(in) :: expected(:) !! Reference vector with the same shape as actual.
        real(dp), intent(in) :: tolerance !! Maximum allowed elementwise absolute error.
        character(len=*), intent(in) :: label !! Human-readable test label printed on failure.

        call assert_true(size(actual) == size(expected), trim(label) // ' size')
        if (size(actual) > 0) then
            if (maxval(abs(actual - expected)) > tolerance) then
                write (*, '(a,1x,a,1x,es24.16)') 'FAIL:', trim(label), maxval(abs(actual - expected))
                error stop 1
            end if
        end if
    end subroutine assert_close_vector

    subroutine assert_close_complex_vector(actual, expected, tolerance, label)
        complex(dp), intent(in) :: actual(:) !! Computed complex vector under test.
        complex(dp), intent(in) :: expected(:) !! Reference complex vector with the same shape as actual.
        real(dp), intent(in) :: tolerance !! Maximum allowed complex magnitude error.
        character(len=*), intent(in) :: label !! Human-readable test label printed on failure.

        call assert_true(size(actual) == size(expected), trim(label) // ' size')
        if (size(actual) > 0) then
            if (maxval(abs(actual - expected)) > tolerance) then
                write (*, '(a,1x,a,1x,es24.16)') 'FAIL:', trim(label), maxval(abs(actual - expected))
                error stop 1
            end if
        end if
    end subroutine assert_close_complex_vector

end program test_signal
