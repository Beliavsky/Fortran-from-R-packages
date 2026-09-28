! SPDX-License-Identifier: GPL-3.0-only
! Derived from R package ftsa 6.7 by Rob Hyndman and Han Lin Shang.
program test_extended
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
   use ftsa, only : dp, ftsa_dynamic_flr, ftsa_dynupdate, ftsa_forecast_fplsr, ftsa_fplsr
   use ftsa, only : ftsa_ftsm_weight_select, ftsa_isfe, ftsa_maf_multivariate
   use ftsa, only : ftsa_pcscore_bootstrap_data, ftsa_stop_time_detect
   implicit none
   real(dp) :: y(4, 8), basis(4), beta_opt
   real(dp) :: multi(4, 30), expected(4)
   real(dp), allocatable :: boot(:, :, :), centers(:, :), coefficients(:, :)
   real(dp), allocatable :: fitted(:, :), forecast(:), forecasts(:, :), isfe(:, :)
   real(dp), allocatable :: loadings(:, :), maf(:, :), prediction(:), reconstruction(:, :), residuals(:, :)
   integer :: break_backward, break_forward, i, j, k_ratio, k_threshold, order_k, order_m, status
   integer, allocatable :: backward_components(:), forward_components(:)
   real(dp), allocatable :: backward_error(:), forward_error(:)
   integer, allocatable :: seed(:)

   basis = [1.0_dp, -0.5_dp, 2.0_dp, 0.75_dp]
   do j = 1, 8
      y(:, j) = [10.0_dp, 20.0_dp, 30.0_dp, 40.0_dp] + real(j, dp) * basis
   end do
   expected = [10.0_dp, 20.0_dp, 30.0_dp, 40.0_dp] + 9.0_dp * basis

   call ftsa_fplsr(y, 1, prediction, coefficients, fitted, residuals, status)
   call check(status == 0 .and. size(prediction) == 4, 'fplsr status')
   call check(maxval(abs(prediction - expected)) < 1.0e-8_dp, 'fplsr exact linear prediction')
   call check(maxval(abs(residuals)) < 1.0e-8_dp, 'fplsr exact fitted values')

   call ftsa_forecast_fplsr(y, 1, 2, forecasts, status)
   call check(status == 0 .and. all(shape(forecasts) == [4, 2]), 'forecastfplsr shape')
   call check(maxval(abs(forecasts(:, 1) - expected)) < 1.0e-8_dp, 'forecastfplsr first step')

   do j = 1, 30
      do i = 1, 4
         multi(i, j) = sin(0.17_dp * real(i * j, dp)) + 0.1_dp * real(i, dp) * cos(0.11_dp * real(j, dp))
      end do
   end do
   call ftsa_maf_multivariate(multi, 0.8_dp, maf, loadings, reconstruction, residuals, &
      k_threshold, k_ratio, status)
   call check(status == 0 .and. all(shape(maf) == [30, 4]), 'MAF status and shape')
   call check(all(ieee_is_finite(maf)) .and. all(ieee_is_finite(loadings)), 'MAF finite')
   call check(k_threshold >= 1 .and. k_ratio >= 1, 'MAF component selection')
   call check(sqrt(sum(residuals * residuals) / real(size(residuals), dp)) < 1.0e-8_dp, 'MAF reconstruction')

   call random_seed(size=i)
   allocate(seed(i))
   seed = [(37 + 11 * j, j=1, i)]
   call random_seed(put=seed)
   call ftsa_pcscore_bootstrap_data(y, 5, 'mean', boot, centers, status)
   call check(status == 0 .and. all(shape(boot) == [4, 8, 5]), 'bootstrap data shape')
   call check(all(shape(centers) == [4, 5]) .and. all(ieee_is_finite(centers)), 'bootstrap centers')

   call ftsa_ftsm_weight_select(y, 1, 2, 'mse', beta_opt, status)
   call check(status == 0 .and. beta_opt > 0.0_dp .and. beta_opt < 1.0_dp, 'weighted FTSM selection')

   call ftsa_isfe(y, 1, 4, [1, 2], isfe, status, drift=.true.)
   call check(status == 0 .and. all(shape(isfe) == [2, 2]), 'ISFE shape')
   call check(maxval(abs(isfe(2, :))) < 1.0e-8_dp, 'ISFE linear drift fit')

   call ftsa_dynupdate(y, expected(1:2), 1, 'ts', 0.01_dp, forecast, status, drift=.true.)
   call check(status == 0 .and. size(forecast) == 2, 'dynupdate time-series shape')
   call check(maxval(abs(forecast - expected(3:4))) < 1.0e-8_dp, 'dynupdate time-series forecast')
   call ftsa_dynupdate(y, expected(1:2), 1, 'pls', 0.01_dp, forecast, status, drift=.true.)
   call check(status == 0 .and. size(forecast) == 2 .and. all(ieee_is_finite(forecast)), 'dynupdate PLS')

   call ftsa_dynamic_flr(y, expected(1:2), 0.9_dp, 0.9_dp, forecast, order_k, order_m, status)
   call check(status == 0 .and. order_k == 1 .and. order_m == 1, 'dynamic FLR component selection')
   call check(maxval(abs(forecast - expected(3:4))) < 1.0e-8_dp, 'dynamic FLR point forecast')

   call ftsa_stop_time_detect(y, forward_error, backward_error, forward_components, backward_components, &
      break_forward, break_backward, status, drift=.true.)
   call check(status == 0 .and. size(forward_error) == 6 .and. size(backward_error) == 6, 'stop-time detection shapes')
   call check(break_forward >= 3 .and. break_forward <= 8, 'forward breakpoint range')
   call check(break_backward >= 1 .and. break_backward <= 6, 'backward breakpoint range')

contains
   subroutine check(condition, message)
      logical, intent(in) :: condition !! Assertion condition expected to be true.
      character(len=*), intent(in) :: message !! Diagnostic emitted if the assertion fails.
      if (.not. condition) error stop message
   end subroutine check
end program test_extended
