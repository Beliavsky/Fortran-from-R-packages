! SPDX-License-Identifier: GPL-3.0-or-later
program evt_and_risk
  use qrmtools, only : dp, fit_gpd_mle, fit_result, var_gpd, es_gpd, &
    hill_estimator, hill_result
  implicit none

  real(dp) :: excesses(12)
  type(fit_result) :: fit
  type(hill_result) :: hill

  excesses = [0.12_dp,0.20_dp,0.31_dp,0.48_dp,0.55_dp,0.72_dp, &
    0.93_dp,1.10_dp,1.35_dp,1.70_dp,2.20_dp,3.10_dp]

  fit = fit_gpd_mle(excesses,estimate_covariance=.false.)
  if(.not.fit%ok) error stop trim(fit%message)

  print '(a,2f12.6)', 'GPD MLE shape and scale: ', fit%parameters
  print '(a,f12.6)', '95% VaR: ', &
    var_gpd(0.95_dp,fit%parameters(1),fit%parameters(2))
  print '(a,f12.6)', '95% ES:  ', &
    es_gpd(0.95_dp,fit%parameters(1),fit%parameters(2))

  hill = hill_estimator(excesses,2,6)
  if(.not.hill%ok) error stop trim(hill%message)
  print '(a,5f12.6)', 'Hill estimates: ', hill%tail_index
end program evt_and_risk
