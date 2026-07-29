! SPDX-License-Identifier: GPL-2.0-only OR GPL-3.0-only
program test_objectives
   use risksimul
   implicit none
   type(portfolio_model) :: model
   type(simulation_result) :: fit
   type(sis_control) :: control
   real(dp) :: corr(2,2), params(2,3)

   corr = reshape([1.0_dp,0.2_dp,0.2_dp,1.0_dp],[2,2])
   params = reshape([0.0_dp,0.0_dp, 0.02_dp,0.025_dp, 6.0_dp,8.0_dp],[2,3])
   model = new_portfolio(8.0_dp,corr,'t',params,weight=[0.45_dp,0.55_dp])
   call assert_true(model%ok)

   allocate(control%allocations(2))
   control%allocations = [500,800]
   control%normal_strata = 3
   control%gamma_strata = 3
   control%minimum_per_stratum = 4
   control%direction_iterations = 25
   control%optimize_conditional_excess = .true.
   control%multi_objective = objective_max_relative
   control%seed = 86420_i8
   fit = stratified_copula(model,[0.93_dp,0.95_dp,0.97_dp],control)
   call assert_true(fit%ok)
   call assert_true(all(fit%tail_probability(:)%estimate >= 0.0_dp))
   call assert_true(fit%tail_probability(1)%estimate <= fit%tail_probability(2)%estimate)
   call assert_true(fit%tail_probability(2)%estimate <= fit%tail_probability(3)%estimate)

   print '(a)', 'test_objectives: PASS'
contains
   subroutine assert_true(condition)
      logical, intent(in) :: condition
      if (.not. condition) error stop 1
   end subroutine assert_true
end program test_objectives
