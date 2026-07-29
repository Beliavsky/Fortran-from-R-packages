! SPDX-License-Identifier: GPL-2.0-only OR GPL-3.0-only
program test_portfolio
   use risksimul
   implicit none
   type(portfolio_model) :: model
   type(importance_parameters) :: imp
   real(dp) :: corr(2,2), params(2,3), z(2), value

   corr = reshape([1.0_dp,0.25_dp,0.25_dp,1.0_dp],[2,2])
   params = reshape([0.01_dp,-0.02_dp, 0.015_dp,0.020_dp, 5.0_dp,8.0_dp],[2,3])
   model = new_portfolio(7.0_dp,corr,'t',params,weight=[0.4_dp,0.6_dp])
   call assert_true(model%ok)
   z = 0.0_dp
   value = portfolio_return_one(z,model%copula_df,model)
   call assert_close(value,0.4_dp*exp(0.01_dp)+0.6_dp*exp(-0.02_dp),2.0e-12_dp)
   call assert_close(tail_loss_response(0.9_dp,0.95_dp),1.0_dp,0.0_dp)
   call assert_close(excess_response(0.9_dp,0.95_dp),0.1_dp,2.0e-16_dp)

   imp = algorithm_3(model,0.97_dp,40)
   call assert_true(imp%ok)
   call assert_true(size(imp%shift) == 2)
   call assert_true(imp%gamma_mean > 0.0_dp)
   call assert_true(portfolio_return_one(imp%shift,imp%gamma_mean,model) < 0.971_dp)

   print '(a)', 'test_portfolio: PASS'
contains
   subroutine assert_true(condition)
      logical, intent(in) :: condition
      if (.not. condition) error stop 1
   end subroutine assert_true
   subroutine assert_close(actual,expected,tolerance)
      real(dp), intent(in) :: actual,expected,tolerance
      if (abs(actual-expected) > tolerance) then
         print '(a,3es25.16)', 'mismatch: ',actual,expected,abs(actual-expected)
         error stop 1
      end if
   end subroutine assert_close
end program test_portfolio
