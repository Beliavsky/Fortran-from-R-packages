! SPDX-License-Identifier: GPL-3.0-or-later
program qrmtools_demo
  use qrmtools, only : dp, black_scholes, var_gpd, es_gpd, &
    fit_gpd_mom
  implicit none

  real(dp) :: parameters(2)
  real(dp) :: excesses(6)

  excesses = [0.2_dp,0.5_dp,0.8_dp,1.1_dp,1.7_dp,2.4_dp]
  parameters = fit_gpd_mom(excesses)

  print '(a,f12.6)', 'Black-Scholes call: ', &
    black_scholes(0.0_dp,100.0_dp,0.05_dp,0.2_dp,100.0_dp,1.0_dp)
  print '(a,2f12.6)', 'GPD MOM shape and scale: ', parameters
  print '(a,f12.6)', 'GPD 99% VaR: ', &
    var_gpd(0.99_dp,parameters(1),parameters(2))
  print '(a,f12.6)', 'GPD 99% ES:  ', &
    es_gpd(0.99_dp,parameters(1),parameters(2))
end program qrmtools_demo
