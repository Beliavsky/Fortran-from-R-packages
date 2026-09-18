! SPDX-License-Identifier: MIT
module r_stats_diagnostics
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_quiet_nan, ieee_value
   use r_kinds, only: dp
   use r_optional, only: optval
   use r_distributions, only: r_pchisq
   use r_stats_time_series, only: acf
   use r_stats_types, only: acf_result_t, chisq_test_result_t
   implicit none
   private
   public :: box_test

contains

   pure function box_test(x, lag, method, fitdf) result(test)
      !! Tests serial correlation using Box-Pierce or Ljung-Box sums of squared autocorrelations.
      real(dp), intent(in) :: x(:) !! Finite univariate observations or fitted residuals.
      integer, intent(in), optional :: lag !! Largest tested lag, at least one and below sample size.
      character(len=*), intent(in), optional :: method !! Box-Pierce (default) or Ljung-Box.
      integer, intent(in), optional :: fitdf !! Nonnegative fitted-parameter count, strictly smaller than lag.
      type(chisq_test_result_t) :: test !! Statistic, adjusted degrees of freedom, and upper-tail probability.
      type(acf_result_t) :: correlations
      character(len=:), allocatable :: selected
      integer :: maximum_lag, fitted_df, n, k
      real(dp) :: sample_size

      n = size(x)
      maximum_lag = optval(lag, 1)
      fitted_df = optval(fitdf, 0)
      selected = 'Box-Pierce'
      if (present(method)) selected = trim(method)
      if (selected /= 'Box-Pierce' .and. selected /= 'Ljung-Box') then
         error stop 'box_test: unsupported method'
      end if
      if (maximum_lag < 1 .or. maximum_lag >= n .or. &
          fitted_df < 0 .or. fitted_df >= maximum_lag) error stop 'box_test: invalid lag or fitdf'
      if (.not. all(ieee_is_finite(x))) error stop 'box_test: observations must be finite'
      test%parameter = maximum_lag - fitted_df
      test%method = merge(7, 8, selected == 'Box-Pierce')
      if (maxval(x) == minval(x)) then
         test%statistic = ieee_value(0.0_dp, ieee_quiet_nan)
         test%p_value = test%statistic
         return
      end if
      correlations = acf(x, lag_max=maximum_lag)
      if (correlations%status /= 0) error stop 'box_test: autocorrelation calculation failed'
      sample_size = real(n, dp)
      if (selected == 'Box-Pierce') then
         test%statistic = sample_size*sum(correlations%value(2:)**2)
      else
         test%statistic = 0.0_dp
         do k = 1, maximum_lag
            test%statistic = test%statistic + correlations%value(k + 1)**2/real(n - k, dp)
         end do
         test%statistic = sample_size*(sample_size + 2.0_dp)*test%statistic
      end if
      test%p_value = r_pchisq(test%statistic, real(test%parameter, dp), lower_tail=.false.)
   end function box_test

end module r_stats_diagnostics
