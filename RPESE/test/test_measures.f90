program test_measures
  use rpese, only : dp, rpese_options, rpese_success, point_estimate, &
    mean_measure, sd_measure, semisd_measure, var_measure, es_measure, &
    lpm_measure, omega_ratio_measure, rachev_ratio_measure, robust_mean_measure
  implicit none
  real(dp) :: x(5), value
  type(rpese_options) :: options
  integer :: status
  character(len=160) :: message

  x = [-2.0_dp, -1.0_dp, 0.0_dp, 1.0_dp, 2.0_dp]
  call assert_close(mean_measure(x), 0.0_dp, 1.0e-12_dp, 'mean')
  call assert_close(sd_measure(x), sqrt(2.5_dp), 1.0e-12_dp, 'sample sd')
  call assert_close(semisd_measure(x), 1.0_dp, 1.0e-12_dp, 'semi sd')
  call assert_close(var_measure(x, 0.2_dp), 1.2_dp, 1.0e-12_dp, 'VaR')
  call assert_close(es_measure(x, 0.2_dp, status), 2.0_dp, 1.0e-12_dp, 'ES')
  call assert_true(status == rpese_success, 'ES status')
  call assert_close(lpm_measure(x, 0.0_dp, 1), 0.6_dp, 1.0e-12_dp, 'LPM1')
  call assert_close(lpm_measure(x, 0.0_dp, 2), 1.0_dp, 1.0e-12_dp, 'LPM2')
  call assert_close(omega_ratio_measure(x, 0.0_dp, status), 1.0_dp, 1.0e-12_dp, 'Omega')
  call assert_close(rachev_ratio_measure(x, 0.2_dp, 0.2_dp, status), 1.0_dp, 1.0e-12_dp, 'Rachev')
  call assert_close(robust_mean_measure(x, 'mopt', 0.95_dp, status), 0.0_dp, 1.0e-8_dp, 'robust mean')

  options = rpese_options()
  options%risk_free = 0.25_dp
  options%source_compatibility = .true.
  call point_estimate(x + 1.0_dp, 'DSR', value, options, status, message)
  call assert_true(status == rpese_success, 'DSR source status')
  call assert_close(value, 1.0_dp / sqrt(2.0_dp), 1.0e-12_dp, 'DSR source compatibility')
  options%source_compatibility = .false.
  call point_estimate(x + 1.0_dp, 'DSR', value, options, status, message)
  call assert_close(value, 0.75_dp / sqrt(2.0_dp), 1.0e-12_dp, 'DSR corrected')

  print '(a)', 'test_measures: PASS'
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
end program test_measures
