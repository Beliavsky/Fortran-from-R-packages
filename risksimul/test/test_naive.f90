! SPDX-License-Identifier: GPL-2.0-only OR GPL-3.0-only
program test_naive
   use risksimul
   implicit none
   type(portfolio_model) :: model
   type(simulation_result) :: a, b
   real(dp) :: corr(2,2), params(2,3)

   corr = reshape([1.0_dp,0.4_dp,0.4_dp,1.0_dp],[2,2])
   params = reshape([0.0_dp,0.0_dp, 0.02_dp,0.03_dp, 6.0_dp,9.0_dp],[2,3])
   model = new_portfolio(8.0_dp,corr,'t',params,weight=[0.5_dp,0.5_dp])
   call assert_true(model%ok)
   a = naive_copula(12000,model,[0.94_dp,0.97_dp],12345_i8)
   b = NVTCopula(12000,model,[0.94_dp,0.97_dp],12345_i8)
   call assert_true(a%ok .and. b%ok)
   call assert_close(a%tail_probability(1)%estimate,b%tail_probability(1)%estimate,0.0_dp)
   call assert_true(a%tail_probability(1)%estimate <= a%tail_probability(2)%estimate)
   call assert_true(a%tail_probability(1)%variance > 0.0_dp)
   call assert_true(a%conditional_excess(1)%estimate > 0.0_dp)
   call assert_true(a%samples_used == 12000)

   print '(a)', 'test_naive: PASS'
contains
   subroutine assert_true(condition)
      logical, intent(in) :: condition
      if (.not. condition) error stop 1
   end subroutine assert_true
   subroutine assert_close(actual,expected,tolerance)
      real(dp), intent(in) :: actual,expected,tolerance
      if (abs(actual-expected) > tolerance) error stop 1
   end subroutine assert_close
end program test_naive
