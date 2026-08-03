! SPDX-License-Identifier: GPL-2.0-or-later
program test_misc
   use jfe
   implicit none

   real(dp), parameter :: r(6) = [0.02_dp, -0.01_dp, 0.03_dp, -0.02_dp, 0.01_dp, -0.005_dp]
   real(dp), parameter :: e(4) = [1.0_dp, 0.5_dp, -0.25_dp, 0.125_dp]
   type(durbin_h_result) :: dh
   type(annualized_summary) :: summary
   real(dp) :: sr, adjusted

   dh = durbin_h(e, 10, 3, 0.01_dp)
   if (dh%status /= jfe_success) error stop 'FAIL: Durbin h status'
   call check_close(dh%durbin_watson, 0.7176470588235294_dp, 1.0e-14_dp, 'Durbin-Watson')
   call check_close(dh%statistic, 1.8907262612804454_dp, 1.0e-13_dp, 'Durbin h')
   call check_close(dh%p_value, 0.05866089495420775_dp, 1.0e-13_dp, 'Durbin h p-value')

   call check_close(bernardo_ledoit_ratio(r), 0.06_dp/0.035_dp, 1.0e-14_dp, &
      'Bernardo-Ledoit')
   call check_close(d_ratio(r), 3.0_dp*0.035_dp/(3.0_dp*0.06_dp), 1.0e-14_dp, 'D ratio')
   call check_close(kelly_ratio(r), mean_manual(r)/variance_manual(r)/2.0_dp, &
      1.0e-13_dp, 'Kelly ratio')
   call check_close(skewness_kurtosis_ratio(r), skewness(r, skew_moment)/ &
      kurtosis(r, kurt_moment), 1.0e-13_dp, 'skewness-kurtosis ratio')
   call check_close(prospect_ratio(r, 0.0_dp), &
      (sum(r, mask=r > 0.0_dp) + 2.25_dp*sum(r, mask=r < 0.0_dp))/ &
      (downside_deviation(r)*real(size(r), dp)), 1.0e-13_dp, 'prospect ratio')
   call check_close(martin_ratio(r, scale=12.0_dp), &
      return_annualized(r, 12.0_dp)/ulcer_index(r), 1.0e-13_dp, 'Martin ratio')
   call check_close(pain_ratio(r, scale=12.0_dp), &
      return_annualized(r, 12.0_dp)/pain_index(r), 1.0e-13_dp, 'Pain ratio')

   sr = sharpe_ratio_annualized(r, scale=12.0_dp)
   adjusted = adjusted_sharpe_ratio(r, scale=12.0_dp)
   call check_close(adjusted, sr*(1.0_dp + sr*skewness(r)/6.0_dp - &
      sr*sr*kurtosis(r, kurt_excess)/24.0_dp), 1.0e-13_dp, 'adjusted Sharpe')

   summary = table_annualized_returns(r, 12.0_dp)
   if (summary%status /= jfe_success) error stop 'FAIL: summary status'
   call check_close(summary%annualized_return, return_annualized(r, 12.0_dp), &
      1.0e-14_dp, 'summary return')
   call check_close(summary%annualized_sharpe, sharpe_ratio_annualized(r, scale=12.0_dp), &
      1.0e-14_dp, 'summary Sharpe')

   print '(a)', 'test_misc: PASS'

contains

   pure real(dp) function mean_manual(x) result(mu)
      real(dp), intent(in) :: x(:)
      mu = sum(x)/real(size(x), dp)
   end function mean_manual

   pure real(dp) function variance_manual(x) result(v)
      real(dp), intent(in) :: x(:)
      real(dp) :: mu
      mu = mean_manual(x)
      v = sum((x - mu)**2)/real(size(x) - 1, dp)
   end function variance_manual

   subroutine check_close(actual, expected, tolerance, label)
      real(dp), intent(in) :: actual, expected, tolerance
      character(*), intent(in) :: label
      if (abs(actual - expected) > tolerance) then
         print '(a,2(1x,es24.16))', 'FAIL: '//trim(label), actual, expected
         error stop 1
      end if
   end subroutine check_close

end program test_misc
