program test_nuisance_stats
  use rpeif, only : dp, nuisance_parameters, nuisance_parameters_fn, &
    lower_partial_moment, upper_partial_moment, rpeif_success
  use rpeif_stats, only : normal_cdf, normal_quantile, quantile_type7, sample_sd
  implicit none
  type(nuisance_parameters) :: pars
  real(dp) :: x(5), q, lpm, upm
  integer :: status

  x = [-2.0_dp, -1.0_dp, 0.0_dp, 1.0_dp, 4.0_dp]
  call nuisance_parameters_fn(pars, status=status)
  call assert_true(status == rpeif_success, 'default nuisance status')
  call assert_close(normal_cdf(normal_quantile(0.1_dp)), 0.1_dp, 2.0e-12_dp, 'normal inverse')
  call assert_close(pars%semisd, pars%sd / sqrt(2.0_dp), 1.0e-14_dp, 'semisd')
  call assert_close(pars%upm1, pars%lpm1 + pars%mu - pars%c, 1.0e-14_dp, 'partial-moment identity')
  call assert_true(pars%fq_alpha > 0.0_dp, 'positive density')

  q = quantile_type7(x, 0.25_dp)
  call assert_close(q, -1.0_dp, 1.0e-14_dp, 'type-7 quantile')
  call assert_close(sample_sd(x), sqrt(5.3_dp), 1.0e-14_dp, 'sample sd')
  lpm = lower_partial_moment(x, 0.0_dp, 1)
  upm = upper_partial_moment(x, 0.0_dp, 1)
  call assert_close(lpm, 0.6_dp, 1.0e-14_dp, 'LPM')
  call assert_close(upm, 1.0_dp, 1.0e-14_dp, 'UPM')

  print '(a)', 'test_nuisance_stats: PASS'
contains
  subroutine assert_close(actual, expected, tolerance, label)
    real(dp), intent(in) :: actual, expected, tolerance
    character(len=*), intent(in) :: label
    if (abs(actual - expected) > tolerance) then
      print '(a,2es24.14)', trim(label)//' failed: ', actual, expected
      error stop 1
    end if
  end subroutine assert_close

  subroutine assert_true(condition, label)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: label
    if (.not. condition) then
      print '(a)', trim(label)//' failed'
      error stop 1
    end if
  end subroutine assert_true
end program test_nuisance_stats
