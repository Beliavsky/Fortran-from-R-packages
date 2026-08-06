program test_timeseries
  use rpese, only : dp, pi, rpese_options, periodogram_result, rpese_success, &
    frequency_all, frequency_decimate, frequency_truncate, fit_periodogram, &
    lag1_correlation, fit_ar1, polynomial_design
  implicit none
  integer, parameter :: n = 64
  real(dp) :: x(n), phi, intercept
  real(dp), allocatable :: residuals(:), design(:, :), means(:), scales(:)
  type(rpese_options) :: options
  type(periodogram_result) :: periodogram
  integer :: i, status, peak

  do i = 1, n
    x(i) = sin(2.0_dp * pi * 5.0_dp * real(i - 1, dp) / real(n, dp))
  end do
  options = rpese_options()
  options%frequency_mode = frequency_all
  call fit_periodogram(x, periodogram, options)
  call assert_true(periodogram%status == rpese_success, 'periodogram status')
  peak = maxloc(periodogram%spectrum, dim=1)
  call assert_close(periodogram%frequency(peak), 5.0_dp / real(n, dp), 1.0e-12_dp, 'periodogram peak')

  options%frequency_mode = frequency_decimate
  options%frequency_fraction = 0.5_dp
  call fit_periodogram(x, periodogram, options)
  call assert_true(size(periodogram%frequency) == 16, 'decimated frequency count')
  options%frequency_mode = frequency_truncate
  call fit_periodogram(x, periodogram, options)
  call assert_true(size(periodogram%frequency) == 15, 'truncated frequency count')

  x(1) = 0.25_dp
  do i = 2, n
    x(i) = 0.1_dp + 0.7_dp * x(i - 1) + 0.01_dp * sin(real(i, dp))
  end do
  call fit_ar1(x, phi, intercept, residuals, status)
  call assert_true(status == rpese_success, 'AR1 status')
  call assert_true(phi > 0.6_dp .and. phi < 0.8_dp, 'AR1 coefficient')
  call assert_true(abs(lag1_correlation(residuals)) < 0.5_dp, 'AR1 residual correlation')

  call polynomial_design([0.1_dp, 0.2_dp, 0.3_dp], 3, design, .true., means, scales, status)
  call assert_true(status == rpese_success, 'design status')
  call assert_true(all(shape(design) == [3, 4]), 'design shape')
  call assert_close(sum(design(:, 2)), 0.0_dp, 1.0e-12_dp, 'standardized design mean')

  print '(a)', 'test_timeseries: PASS'
contains
  subroutine assert_close(actual, expected, tolerance, label)
    real(dp), intent(in) :: actual, expected, tolerance
    character(len=*), intent(in) :: label
    if (abs(actual - expected) > tolerance) then
      print '(a,2es24.14)', trim(label) // ' failed: ', actual, expected
      error stop 1
    end if
  end subroutine assert_close

  subroutine assert_true(condition, label)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: label
    if (.not. condition) then
      print '(a)', trim(label) // ' failed'
      error stop 1
    end if
  end subroutine assert_true
end program test_timeseries
