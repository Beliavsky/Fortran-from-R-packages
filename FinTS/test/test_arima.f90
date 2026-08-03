! SPDX-License-Identifier: GPL-2.0-or-later
program test_arima
   use fints
   implicit none

   integer, parameter :: n = 240
   real(dp) :: x(n), innovations(n), z(n), xreg(n, 1), random_walk(n)
   real(dp) :: seasonal_series(n)
   type(arima_result) :: fit, reg_fit, diff_fit, seasonal_fit
   integer :: i

   do i = 1, n
      innovations(i) = 0.22_dp * sin(0.71_dp * real(i, dp)) + &
         0.13_dp * cos(1.17_dp * real(i, dp))
   end do
   x(1) = 2.0_dp + innovations(1)
   do i = 2, n
      x(i) = 2.0_dp + 0.65_dp * (x(i - 1) - 2.0_dp) + innovations(i)
   end do
   call ARIMA(x, [1, 0, 0], fit, include_mean=.true., max_iterations=2500, tolerance=1.0e-9_dp)
   call check(fit%status == fints_ok .or. fit%status == fints_iteration_limit, 'arima status')
   call check(abs(fit%ar(1) - 0.65_dp) < 0.08_dp, 'arima ar coefficient')
   call check(abs(fit%intercept - 2.0_dp) < 0.08_dp, 'arima intercept')
   call check(fit%sigma2 > 0.0_dp .and. fit%n_used > 200, 'arima variance')
   call check(fit%box_test%status == fints_ok, 'arima box test')

   do i = 1, n
      xreg(i, 1) = sin(0.09_dp * real(i, dp))
      z(i) = 1.5_dp + 0.8_dp * xreg(i, 1) + innovations(i)
   end do
   call ARIMA(z, [0, 0, 0], reg_fit, xreg=xreg, include_mean=.true., &
      max_iterations=1200, tolerance=1.0e-10_dp)
   call check(abs(reg_fit%intercept - 1.5_dp) < 0.03_dp, 'regression intercept')
   call check(abs(reg_fit%regression(1) - 0.8_dp) < 0.03_dp, 'regression coefficient')
   call check(reg_fit%r_squared > 0.80_dp, 'regression r squared')

   random_walk(1) = innovations(1)
   do i = 2, n
      random_walk(i) = random_walk(i - 1) + innovations(i)
   end do
   call ARIMA(random_walk, [0, 1, 0], diff_fit)
   call check(diff_fit%status == fints_ok, 'differenced fit status')
   call check(diff_fit%n_used == n - 1, 'differenced n used')
   call check(diff_fit%sigma2 > 0.0_dp, 'differenced variance')

   do i = 1, 12
      seasonal_series(i) = 1.0_dp + 0.20_dp * sin(0.173_dp * real(i * i, dp)) + &
         0.12_dp * cos(0.319_dp * real(i * i, dp))
   end do
   do i = 13, n
      seasonal_series(i) = 1.0_dp + 0.55_dp * (seasonal_series(i - 12) - 1.0_dp) + &
         0.20_dp * sin(0.173_dp * real(i * i, dp)) + &
         0.12_dp * cos(0.319_dp * real(i * i, dp))
   end do
   call ARIMA(seasonal_series, [0, 0, 0], seasonal_fit, seasonal_order=[1, 0, 0], &
      seasonal_period=12, include_mean=.true., max_iterations=2500, tolerance=1.0e-9_dp)
   call check(seasonal_fit%status == fints_ok .or. &
      seasonal_fit%status == fints_iteration_limit, 'seasonal status')
   call check(abs(seasonal_fit%seasonal_ar(1) - 0.55_dp) < 0.10_dp, &
      'seasonal ar coefficient')
   call check(abs(seasonal_fit%intercept - 1.0_dp) < 0.10_dp, 'seasonal intercept')

   print '(a)', 'test_arima: PASS'

contains

   subroutine check(condition, label)
      logical, intent(in) :: condition
      character(len=*), intent(in) :: label
      if (.not. condition) then
         write(*, '(a,a)') 'FAIL: ', trim(label)
         error stop 1
      end if
   end subroutine check

end program test_arima
