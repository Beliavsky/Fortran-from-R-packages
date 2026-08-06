program test_dispatch_shape
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use rpeif, only : dp, rpeif_options, influence_result, influence_series, &
    evaluate_shape, supported_estimator, rpeif_success, rpeif_numerical_failure
  implicit none
  real(dp) :: returns(20)
  type(rpeif_options) :: opts
  type(influence_result) :: series, shape
  integer :: i

  do i = 1, size(returns)
    returns(i) = 0.003_dp + 0.025_dp * sin(0.7_dp * real(i, dp))
  end do
  returns(7) = -0.12_dp
  opts%prewhiten = .true.
  opts%ar_order = 1
  opts%clean_outliers = .true.
  opts%efficiency = 0.99_dp
  call influence_series('ES', returns, series, opts)
  call assert_true(series%status == rpeif_success .or. series%status == rpeif_numerical_failure, 'series status')
  call assert_true(size(series%values) == size(returns), 'series size')
  call assert_true(all(ieee_is_finite(series%values)), 'series finite')

  opts%prewhiten = .false.
  opts%clean_outliers = .false.
  call evaluate_shape('Mean', shape, opts, k=2, step=0.01_dp)
  call assert_true(shape%status == rpeif_success, 'shape status')
  call assert_true(size(shape%x) == 29, 'shape grid size')
  call assert_true(all(ieee_is_finite(shape%values)), 'shape finite')
  call assert_true(supported_estimator('Omega'), 'Omega alias')
  call assert_true(.not. supported_estimator('not-a-measure'), 'unknown estimator')

  print '(a)', 'test_dispatch_shape: PASS'
contains
  subroutine assert_true(condition, label)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: label
    if (.not. condition) then
      print '(a)', trim(label)//' failed'
      error stop 1
    end if
  end subroutine assert_true
end program test_dispatch_shape
