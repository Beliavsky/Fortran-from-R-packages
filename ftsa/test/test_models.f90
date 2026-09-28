! SPDX-License-Identifier: GPL-3.0-only
! Derived from R package ftsa 6.7 by Rob Hyndman and Han Lin Shang.
program test_models
   use ftsa, only : dp, ftsm_fit, ftsa_forecast_ftsm, ftsa_ftsm, ftsa_ftsm_iterative_forecasts
   implicit none
   type(ftsm_fit) :: fit
   real(dp) :: y(4, 6), mean_curve(4), basis(4), scores(6)
   real(dp), allocatable :: forecast(:, :), iterative(:, :)
   integer :: j, status

   mean_curve = [10.0_dp, 20.0_dp, 30.0_dp, 40.0_dp]
   basis = [1.0_dp, -1.0_dp, 2.0_dp, -2.0_dp]
   basis = basis / sqrt(sum(basis * basis))
   scores = [-3.0_dp, -2.0_dp, -1.0_dp, 1.0_dp, 2.0_dp, 3.0_dp]
   do j = 1, 6
      y(:, j) = mean_curve + basis * scores(j)
   end do
   call ftsa_ftsm(y, 1, fit, status)
   call check(status == 0, 'ftsm status')
   call check(maxval(abs(fit%fitted - y)) < 1.0e-9_dp, 'ftsm reconstruction')
   call ftsa_forecast_ftsm(fit, 3, forecast, status)
   call check(status == 0 .and. all(shape(forecast) == [4, 3]), 'forecast shape')
   call check(maxval(abs(forecast(:, 1) - y(:, 6))) < 1.0e-9_dp, 'random-walk forecast')
   call ftsa_ftsm_iterative_forecasts(y, 1, 3, iterative, status)
   call check(status == 0 .and. all(shape(iterative) == [4, 3]), 'iterative forecasts')
   call check(maxval(abs(iterative(:, 1) - y(:, 6))) < 1.0e-9_dp, 'iterative first forecast')

contains
   subroutine check(condition, message)
      logical, intent(in) :: condition !! Assertion condition expected to be true.
      character(len=*), intent(in) :: message !! Diagnostic emitted if the assertion fails.
      if (.not. condition) error stop message
   end subroutine check
end program test_models
