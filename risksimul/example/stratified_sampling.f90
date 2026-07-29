! SPDX-License-Identifier: GPL-2.0-only OR GPL-3.0-only
program stratified_sampling
   use risksimul
   implicit none
   type(portfolio_model) :: portfolio
   type(simulation_result) :: result
   type(sis_control) :: control
   real(dp) :: correlation(2,2), parameters(2,3)

   correlation = reshape([1.0_dp,0.5_dp,0.5_dp,1.0_dp],[2,2])
   parameters = reshape([0.0_dp,0.0_dp,0.02_dp,0.025_dp,5.0_dp,7.0_dp],[2,3])
   portfolio = new_portfolio(7.5_dp,correlation,'t',parameters, &
      weight=[0.5_dp,0.5_dp])
   if (.not. portfolio%ok) error stop trim(portfolio%message)

   allocate(control%allocations(3))
   control%allocations = [1000,2000,4000]
   control%normal_strata = 5
   control%gamma_strata = 5
   control%multi_objective = objective_msre
   control%seed = 20260727_i8

   result = stratified_copula(portfolio,[0.93_dp,0.95_dp,0.97_dp],control)
   if (.not. result%ok) error stop trim(result%message)

   print '(a,i0)', 'Actual samples used: ',result%samples_used
   print '(a,3es14.5)', 'Tail probabilities: ', &
      result%tail_probability(:)%estimate
end program stratified_sampling
