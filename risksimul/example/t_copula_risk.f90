! SPDX-License-Identifier: GPL-2.0-only OR GPL-3.0-only
program t_copula_risk
   use risksimul
   implicit none
   type(portfolio_model) :: portfolio
   type(simulation_result) :: result
   real(dp) :: correlation(3,3), parameters(3,3)

   correlation = reshape([ &
      1.0_dp,0.45_dp,0.30_dp, &
      0.45_dp,1.0_dp,0.35_dp, &
      0.30_dp,0.35_dp,1.0_dp],[3,3])
   parameters = reshape([ &
      0.0002_dp,0.0001_dp,-0.0001_dp, &
      0.012_dp,0.015_dp,0.018_dp, &
      5.0_dp,7.0_dp,9.0_dp],[3,3])

   portfolio = new_portfolio(8.0_dp,correlation,'t',parameters, &
      weight=[0.40_dp,0.35_dp,0.25_dp])
   if (.not. portfolio%ok) error stop trim(portfolio%message)

   result = NVTCopula(50000,portfolio,[0.95_dp,0.97_dp],20260727_i8)
   if (.not. result%ok) error stop trim(result%message)

   print '(a,2f12.6)', 'Tail probabilities: ', &
      result%tail_probability(1)%estimate,result%tail_probability(2)%estimate
end program t_copula_risk
