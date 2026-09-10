! SPDX-License-Identifier: GPL-2.0-or-later
! Test/example code for the modern Fortran translation of the upstream Rssa package.
program test_forecast
   use rssa, only : dp, lrr_ssa, mssa_result, parestimate_ssa, period_estimate, rforecast_ssa
   use rssa, only : roots_lrr, ssa_1d, ssa_mssa, ssa_result, vforecast_mssa
   implicit none
   real(dp) :: x(6), channels(6, 2)
   real(dp), allocatable :: coefficients(:), prediction(:), mpred(:, :)
   complex(dp), allocatable :: roots(:)
   type(ssa_result) :: fit
   type(mssa_result) :: mfit
   type(period_estimate) :: estimate
   integer :: i, info
   x = [(2.0_dp**real(i - 1, dp), i=1, 6)]
   fit = ssa_1d(x, 3, 1)
   coefficients = lrr_ssa(fit, [1], info=info)
   prediction = rforecast_ssa(fit, [1], 2, only_new=.true., info=info)
   call assert_close(prediction(1), 64.0_dp, 2.0e-6_dp, 'first recurrent forecast')
   call assert_close(prediction(2), 128.0_dp, 4.0e-6_dp, 'second recurrent forecast')
   roots = roots_lrr(coefficients, info)
   call assert_true(size(roots) == size(coefficients), 'LRR roots')
   estimate = parestimate_ssa(fit, [1], info=info)
   call assert_true(estimate%info == 0, 'parameter estimate')
   channels(:, 1) = x
   channels(:, 2) = 3.0_dp * x
   mfit = ssa_mssa(channels, 3, 1)
   mpred = vforecast_mssa(mfit, [1], len=1, only_new=.true., info=info)
   call assert_true(info == 0 .and. all(shape(mpred) == [1, 2]), 'MSSA vector forecast')
   print *, 'test_forecast: PASS'
contains
   subroutine assert_true(condition, message)
      logical, intent(in) :: condition !! Condition required to be true.
      character(len=*), intent(in) :: message !! Failure diagnostic.
      if (.not. condition) error stop message
   end subroutine assert_true
   subroutine assert_close(value, expected, tolerance, message)
      real(dp), intent(in) :: value !! Computed scalar value.
      real(dp), intent(in) :: expected !! Expected scalar reference.
      real(dp), intent(in) :: tolerance !! Maximum absolute error.
      character(len=*), intent(in) :: message !! Failure diagnostic.
      if (abs(value - expected) > tolerance) error stop message
   end subroutine assert_close
end program test_forecast
