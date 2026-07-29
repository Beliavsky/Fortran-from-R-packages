! SPDX-License-Identifier: GPL-2.0-only OR GPL-3.0-only
program test_gh
   use risksimul
   implicit none
   type(portfolio_model) :: model
   type(simulation_result) :: fit
   real(dp) :: corr(1,1), params(1,5), z(1), center, low, high

   corr = 1.0_dp
   params(1,:) = [-0.7_dp,5.0_dp,-0.4_dp,0.03_dp,0.001_dp]
   model = new_portfolio(8.0_dp,corr,'GH',params,gh_grid_size=257)
   call assert_true(model%ok)
   z = 0.0_dp
   center = portfolio_return_one(z,8.0_dp,model)
   z = -2.0_dp
   low = portfolio_return_one(z,8.0_dp,model)
   z = 2.0_dp
   high = portfolio_return_one(z,8.0_dp,model)
   call assert_true(low < center .and. center < high)
   fit = naive_copula(2000,model,[0.96_dp],777_i8)
   call assert_true(fit%ok)
   call assert_true(fit%tail_probability(1)%estimate >= 0.0_dp)
   call assert_true(fit%tail_probability(1)%estimate <= 1.0_dp)

   print '(a)', 'test_gh: PASS'
contains
   subroutine assert_true(condition)
      logical, intent(in) :: condition
      if (.not. condition) error stop 1
   end subroutine assert_true
end program test_gh
