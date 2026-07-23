! Part of the modern Fortran translation of longmemo 1.1-4.
! This file was created or modified for the Fortran project on 2026-07-23.
! Original longmemo authors retain copyright; see ORIGINAL_PACKAGE.txt.
! SPDX-License-Identifier: GPL-2.0-or-later

program test_longmemo
    use longmemo_kinds, only : dp, pi
    use longmemo_fft, only : fft_forward, fft_inverse
    use longmemo_stats, only : set_random_seed
    use longmemo_io, only : read_index_value_csv
    use longmemo
    implicit none

    complex(dp), allocatable :: cx(:), cy(:), cz(:)
    real(dp), allocatable :: ac(:), cov(:, :), pgram(:), x(:), data_x(:), b_fast(:), b_slow(:), freq(:)
    real(dp) :: eta1(1), eta3(3), relative_error
    type(spectrum_result) :: spec
    type(whittle_result) :: fit
    type(fexp_result) :: fexp
    integer :: i

    eta1 = [0.7_dp]
    eta3 = [0.7_dp, 0.6_dp, -0.3_dp]

    allocate(cx(7))
    do i = 1, size(cx)
        cx(i) = cmplx(real(i, dp), real(2*i - 1, dp), dp)
    end do
    call fft_forward(cx, cy)
    call fft_inverse(cy, cz)
    call assert_close(maxval(abs(cz - cx)), 0.0_dp, 1.0e-12_dp, "arbitrary-length FFT round trip")

    call ck_fgn0(4, 0.75_dp, ac)
    call assert_close(ac(1), 1.0_dp, 1.0e-14_dp, "fGn variance")
    call assert_close(ac(2), sqrt(2.0_dp) - 1.0_dp, 1.0e-14_dp, "fGn lag-one covariance")

    call ceta_fgn(eta1, cov, m=256)
    call assert_close(cov(1, 1), 0.4862165912_dp, 2.0e-8_dp, "CetaFGN regression value")

    call ceta_arima(eta1, 0, 0, cov, m=256)
    call assert_close(cov(1, 1), 0.7163356699_dp, 2.0e-8_dp, "CetaARIMA(0,0) regression value")

    call ceta_arima(eta3, 1, 1, cov, m=256)
    call assert_close(cov(1, 1), 14.48898898_dp, 2.0e-6_dp, "CetaARIMA(1,1) regression value")

    call spec_fgn(eta1, 16, spec)
    call assert_true(size(spec%spectrum) == 7, "fGn spectrum length")
    call assert_true(all(spec%spectrum > 0.0_dp), "positive fGn spectrum")

    call spec_arima(eta3, 1, 1, 64, spec)
    call assert_true(size(spec%spectrum) == 31, "fARIMA spectrum length")
    call assert_true(all(spec%spectrum > 0.0_dp), "positive fARIMA spectrum")

    allocate(freq(100))
    do i = 1, size(freq)
        freq(i) = pi*real(i, dp)/real(size(freq) + 1, dp)
    end do
    call b_spec_fgn(freq, 0.8_dp, b_fast)
    call b_spec_fgn(freq, 0.8_dp, b_slow, k_approx=0, nsum=1000)
    relative_error = sum(abs(b_fast - b_slow))/sum(abs(b_slow))
    call assert_true(relative_error < 1.0e-4_dp, "Paxson B-spectrum approximation")

    allocate(x(8))
    x = 1.0_dp
    call periodogram(x, pgram)
    call assert_close(pgram(1), 8.0_dp/(2.0_dp*pi), 1.0e-12_dp, "constant-series zero-frequency periodogram")
    call assert_true(maxval(abs(pgram(2:))) < 1.0e-12_dp, "constant-series nonzero periodogram")

    call set_random_seed(162)
    call sim_fgn0(120, 0.7_dp, x)
    call assert_true(sum(x*x) > 0.0_dp, "circulant fGn simulation")
    call sim_arma0(120, 0.7_dp, x)
    call assert_true(sum(x*x) > 0.0_dp, "circulant fARIMA(0,d,0) simulation")

    call sim_fgn_fft(128, 0.7_dp, x)
    call assert_close(sum(x)/real(size(x), dp), 0.0_dp, 1.0e-12_dp, "FFT fGn zero sample mean")
    call assert_true(sum(x*x) > 0.0_dp, "FFT fGn nonzero variance")

    call whittle_estimate(x, "fGn", fit, start_eta=[0.6_dp], covariance_m=256)
    call assert_true(fit%eta(1) > 0.1_dp .and. fit%eta(1) < 0.99_dp, "Whittle H bounds")
    call assert_true(fit%std_error(1) > 0.0_dp, "Whittle standard error")

    call fexp_estimate(x, 2, 0.5_dp, fexp)
    call assert_true(fexp%order_poly >= 0 .and. fexp%order_poly <= 2, "FEXP selected order")
    call assert_true(size(fexp%spectrum) == (size(x) - 1)/2, "FEXP spectrum length")

    call read_index_value_csv("data/NileMin.csv", data_x)
    call fexp_estimate(data_x, 3, 0.5_dp, fexp)
    call assert_close(fexp%hurst, 0.871798_dp, 3.0e-6_dp, "Nile FEXP H regression value")
    call assert_true(fexp%order_poly == 1, "Nile FEXP selected order")

    print '(a)', "All tests passed."

contains

    subroutine assert_close(actual, expected, tolerance, label)
        real(dp), intent(in) :: actual, expected, tolerance
        character(len=*), intent(in) :: label

        if (abs(actual - expected) > tolerance) then
            print '(a)', "FAILED: "//label
            print '(a,es24.16)', "  actual   = ", actual
            print '(a,es24.16)', "  expected = ", expected
            error stop 1
        end if
    end subroutine assert_close


    subroutine assert_true(condition, label)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: label

        if (.not. condition) then
            print '(a)', "FAILED: "//label
            error stop 1
        end if
    end subroutine assert_true

end program test_longmemo
