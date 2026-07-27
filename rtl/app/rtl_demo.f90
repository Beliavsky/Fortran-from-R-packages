! SPDX-License-Identifier: MIT
program rtl_demo
  use rtl, only: dp, option_result, ou_fit_result, frontier_result, path_result
  use rtl, only: gbs_option, sim_ou, fit_ou, efficient_frontier_statistics
  implicit none

  type(option_result) :: option
  type(ou_fit_result) :: fit
  type(frontier_result) :: frontier
  type(path_result) :: ou_path
  real(dp) :: covariance(2, 2)
  real(dp), allocatable :: series(:)

  option = gbs_option(100.0_dp, 100.0_dp, 1.0_dp, 0.05_dp, 0.02_dp, 0.20_dp, "call")
  print '(a,f12.6)', 'GBS call price: ', option%price

  ou_path = sim_ou(1, 8.0_dp, 5.0_dp, 0.5_dp, 0.2_dp, 3.0_dp, 1.0_dp / 12.0_dp, seed=42)
  series = ou_path%values(:, 1)
  fit = fit_ou(series, 1.0_dp / 12.0_dp)
  print '(a,3f12.6)', 'OU theta, mu, sigma: ', fit%theta, fit%mu, fit%sigma

  covariance = reshape([0.04_dp, 0.01_dp, 0.01_dp, 0.09_dp], shape(covariance))
  frontier = efficient_frontier_statistics(100, [0.08_dp, 0.12_dp], covariance, seed=42)
  print '(a,i0)', 'Minimum-risk portfolio row: ', frontier%minimum_risk_index
  print '(a,i0)', 'Maximum-Sharpe portfolio row: ', frontier%maximum_sharpe_index
end program rtl_demo
