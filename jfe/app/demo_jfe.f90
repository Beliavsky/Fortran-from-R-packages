! SPDX-License-Identifier: GPL-2.0-or-later
program demo_jfe
   use jfe
   implicit none

   real(dp), parameter :: asset(10) = [0.011_dp, -0.007_dp, 0.014_dp, 0.006_dp, -0.012_dp, &
      0.010_dp, 0.003_dp, -0.004_dp, 0.009_dp, 0.002_dp]
   real(dp), parameter :: market(10) = [0.008_dp, -0.006_dp, 0.010_dp, 0.004_dp, -0.009_dp, &
      0.007_dp, 0.002_dp, -0.003_dp, 0.006_dp, 0.001_dp]
   type(annualized_summary) :: summary

   summary = table_annualized_returns(asset, 12.0_dp)
   print '(a)', 'JFE performance demo'
   print '(a,f12.6)', 'Annualized return: ', summary%annualized_return
   print '(a,f12.6)', 'Annualized SD:     ', summary%annualized_sd
   print '(a,f12.6)', 'Annualized Sharpe: ', summary%annualized_sharpe
   print '(a,f12.6)', 'Information ratio:', information_ratio(asset, market, 12.0_dp)
   print '(a,f12.6)', 'Adjusted Sharpe:  ', adjusted_sharpe_ratio(asset, scale=12.0_dp)
end program demo_jfe
