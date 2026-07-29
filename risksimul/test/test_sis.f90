! SPDX-License-Identifier: GPL-2.0-only OR GPL-3.0-only
program test_sis
   use risksimul
   implicit none
   type(portfolio_model) :: model
   type(simulation_result) :: naive, sis
   type(sis_control) :: control
   real(dp) :: corr(2,2), params(2,3), difference, standard_error

   corr = reshape([1.0_dp,0.35_dp,0.35_dp,1.0_dp],[2,2])
   params = reshape([0.0_dp,0.0_dp, 0.025_dp,0.030_dp, 5.0_dp,7.0_dp],[2,3])
   model = new_portfolio(7.5_dp,corr,'t',params,weight=[0.5_dp,0.5_dp])
   call assert_true(model%ok)

   allocate(control%allocations(2))
   control%allocations = [1200,2400]
   control%normal_strata = 4
   control%gamma_strata = 4
   control%minimum_per_stratum = 5
   control%direction_iterations = 35
   control%seed = 24680_i8
   sis = stratified_copula(model,[0.94_dp,0.97_dp],control)
   naive = naive_copula(50000,model,[0.94_dp,0.97_dp],97531_i8)
   call assert_true(sis%ok .and. naive%ok)
   call assert_true(sis%tail_probability(1)%estimate <= sis%tail_probability(2)%estimate)
   difference = abs(sis%tail_probability(1)%estimate-naive%tail_probability(1)%estimate)
   standard_error = sqrt(sis%tail_probability(1)%variance+naive%tail_probability(1)%variance)
   call assert_true(difference <= 7.0_dp*standard_error+2.0e-3_dp)
   call assert_true(sis%samples_used >= sum(control%allocations))

   print '(a)', 'test_sis: PASS'
contains
   subroutine assert_true(condition)
      logical, intent(in) :: condition
      if (.not. condition) error stop 1
   end subroutine assert_true
end program test_sis
