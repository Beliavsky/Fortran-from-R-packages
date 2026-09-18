program test_pbs
   use pbs_mod, only : dp, pbs_basis, pbs, predict_pbs, PBS_OK, PBS_OUTSIDE_BOUNDARY
   use, intrinsic :: ieee_arithmetic, only : ieee_is_nan, ieee_value, ieee_quiet_nan
   implicit none

   call test_periodic_explicit()
   call test_df_quantiles()
   call test_predict()
   call test_missing_values()
   call test_outside_periodic()
   call test_nonperiodic_extrapolation()
   print '(a)', 'All pbs tests passed.'

contains

   subroutine assert_true(condition, message)
      logical, intent(in) :: condition !! Assertion condition that must evaluate to true.
      character(len=*), intent(in) :: message !! Failure message printed before stopping the test executable.
      if (.not. condition) then
         print '(a)', 'FAIL: '//trim(message)
         error stop 1
      end if
   end subroutine assert_true

   subroutine assert_close(actual, expected, tol, message)
      real(dp), intent(in) :: actual(:) !! Computed values to compare elementwise.
      real(dp), intent(in) :: expected(:) !! Independently expected reference values.
      real(dp), intent(in) :: tol !! Maximum permitted absolute elementwise difference.
      character(len=*), intent(in) :: message !! Failure message identifying the tested quantity.
      call assert_true(size(actual) == size(expected), message//' size mismatch')
      call assert_true(maxval(abs(actual - expected)) <= tol, message)
   end subroutine assert_close

   subroutine test_periodic_explicit()
      type(pbs_basis) :: fit
      integer :: i
      real(dp), parameter :: x(7) = [0.0_dp, 0.5_dp, 1.0_dp, 2.0_dp, 3.0_dp, 5.5_dp, 6.0_dp]
      real(dp), parameter :: expected(7, 4) = reshape([ &
         0.0_dp, 0.0208333333333333_dp, 0.166666666666667_dp, 0.733333333333333_dp, &
         0.45_dp, 0.00208333333333333_dp, 0.0_dp, &
         0.05_dp, 0.00625_dp, 0.0_dp, 0.1_dp, 0.5_dp, 0.15625_dp, 0.05_dp, &
         0.5_dp, 0.3_dp, 0.1_dp, 0.0_dp, 0.05_dp, 0.58125_dp, 0.5_dp, &
         0.45_dp, 0.672916666666667_dp, 0.733333333333333_dp, 0.166666666666667_dp, &
         0.0_dp, 0.260416666666667_dp, 0.45_dp], [7, 4])

      fit = pbs(x, knots=[1.0_dp, 2.0_dp, 3.0_dp], degree=3, intercept=.true., &
         boundary_knots=[0.0_dp, 6.0_dp], periodic=.true.)
      call assert_true(fit%status == PBS_OK, 'periodic explicit-knot fit status')
      call assert_true(all(shape(fit%values) == [7, 4]), 'periodic explicit-knot dimensions')
      call assert_close(reshape(fit%values, [28]), reshape(expected, [28]), 2.0e-14_dp, &
         'periodic basis values')
      call assert_close(fit%values(1, :), fit%values(7, :), 2.0e-14_dp, 'periodic endpoint closure')
      call assert_close([(sum(fit%values(i, :)), i=1, 7)], [(1.0_dp, i=1, 7)], 2.0e-14_dp, &
         'periodic partition of unity')
   end subroutine test_periodic_explicit

   subroutine test_df_quantiles()
      type(pbs_basis) :: fit
      real(dp), parameter :: x(7) = [0.0_dp, 1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 5.0_dp, 6.0_dp]

      fit = pbs(x, df=4, degree=3, intercept=.false., boundary_knots=[0.0_dp, 6.0_dp], periodic=.true.)
      call assert_true(fit%status == PBS_OK, 'df-selected periodic fit status')
      call assert_true(all(shape(fit%values) == [7, 4]), 'df-selected periodic dimensions')
      call assert_close(fit%knots, [1.2_dp, 2.4_dp, 3.6_dp, 4.8_dp], 2.0e-14_dp, &
         'R type-7 quantile knots')
   end subroutine test_df_quantiles

   subroutine test_predict()
      type(pbs_basis) :: fit
      type(pbs_basis) :: pred
      type(pbs_basis) :: direct
      real(dp), parameter :: x(5) = [0.0_dp, 0.4_dp, 1.2_dp, 3.8_dp, 6.0_dp]
      real(dp), parameter :: newx(4) = [0.1_dp, 1.5_dp, 4.0_dp, 5.9_dp]

      fit = pbs(x, knots=[1.0_dp, 2.0_dp, 3.0_dp], degree=3, intercept=.false., &
         boundary_knots=[0.0_dp, 6.0_dp], periodic=.true.)
      pred = predict_pbs(fit, newx)
      direct = pbs(newx, knots=fit%knots, degree=fit%degree, intercept=fit%intercept, &
         boundary_knots=fit%boundary_knots, periodic=fit%periodic)
      call assert_true(pred%status == PBS_OK, 'predict status')
      call assert_close(reshape(pred%values, [size(pred%values)]), &
         reshape(direct%values, [size(direct%values)]), 1.0e-14_dp, 'predict stored specification')
   end subroutine test_predict

   subroutine test_missing_values()
      type(pbs_basis) :: fit
      real(dp) :: x(5)
      real(dp) :: nan_value

      nan_value = ieee_value(0.0_dp, ieee_quiet_nan)
      x = [0.0_dp, 1.0_dp, nan_value, 3.0_dp, 6.0_dp]
      fit = pbs(x, knots=[1.0_dp, 2.0_dp, 3.0_dp], degree=3, intercept=.true., &
         boundary_knots=[0.0_dp, 6.0_dp], periodic=.true.)
      call assert_true(fit%status == PBS_OK, 'missing-value fit status')
      call assert_true(all(ieee_is_nan(fit%values(3, :))), 'missing x produces missing basis row')
      call assert_true(all(.not. ieee_is_nan(fit%values(1, :))), 'non-missing x produces finite basis row')
   end subroutine test_missing_values

   subroutine test_outside_periodic()
      type(pbs_basis) :: fit
      fit = pbs([-0.1_dp, 1.0_dp, 2.0_dp], knots=[0.5_dp, 1.0_dp, 1.5_dp], degree=3, &
         boundary_knots=[0.0_dp, 2.0_dp], periodic=.true.)
      call assert_true(fit%status == PBS_OUTSIDE_BOUNDARY, 'periodic values outside boundaries are rejected')
   end subroutine test_outside_periodic

   subroutine test_nonperiodic_extrapolation()
      type(pbs_basis) :: fit
      real(dp), parameter :: x(6) = [-1.0_dp, 0.0_dp, 0.5_dp, 1.0_dp, 2.0_dp, 3.0_dp]
      real(dp), parameter :: expected(6, 5) = reshape([ &
         4.0_dp, 1.0_dp, 0.25_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
         -3.5_dp, 0.0_dp, 0.625_dp, 0.5_dp, 0.0_dp, 0.0_dp, &
         0.5_dp, 0.0_dp, 0.125_dp, 0.5_dp, 0.333333333333333_dp, 0.333333333333333_dp, &
         0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.666666666666667_dp, -3.33333333333333_dp, &
         0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 4.0_dp], [6, 5])

      fit = pbs(x, knots=[1.0_dp, 2.0_dp], degree=2, intercept=.true., &
         boundary_knots=[0.0_dp, 2.5_dp], periodic=.false.)
      call assert_true(fit%status == PBS_OK, 'nonperiodic fit status')
      call assert_close(reshape(fit%values, [30]), reshape(expected, [30]), 3.0e-14_dp, &
         'nonperiodic Taylor extrapolation')
   end subroutine test_nonperiodic_extrapolation

end program test_pbs
