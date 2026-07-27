! SPDX-License-Identifier: GPL-2.0-only
! Based on FinancialMath 0.1.1, Copyright (C) 2016 Kameron Penn and Jack Schmidt.
program test_derivatives
   use financialmath
   use test_support
   implicit none
   type(option_order1_t) :: bs
   type(payoff_table_t) :: p, q
   type(forward_result_t) :: fw
   real(dp) :: call_value, put_value
   integer :: j

   call_value = black_scholes_call(100.0_dp, 100.0_dp, 0.05_dp, 1.0_dp, 0.20_dp)
   put_value = black_scholes_put(100.0_dp, 100.0_dp, 0.05_dp, 1.0_dp, 0.20_dp)
   call assert_close(call_value, 10.450583572185565_dp)
   call assert_close(put_value, 5.573526022256971_dp)
   call assert_close(call_value-put_value, 100.0_dp-100.0_dp*exp(-0.05_dp))

   bs = bls_order1(100.0_dp, 100.0_dp, 0.05_dp, 1.0_dp, 0.20_dp, 0.0_dp)
   call assert_true(bs%status%ok, 'Black-Scholes order-one failed')
   call assert_close(bs%call_price, call_value)
   call assert_close(bs%put_price, put_value)
   call assert_close(bs%call_delta-bs%put_delta, 1.0_dp)
   call assert_true(bs%vega > 0.0_dp, 'vega must be positive')

   p = option_call(100.0_dp, 100.0_dp, 0.05_dp, 1.0_dp, 0.20_dp, 'long')
   q = option_call(100.0_dp, 100.0_dp, 0.05_dp, 1.0_dp, 0.20_dp, 'short')
   do j = 1, size(p%stock)
      call assert_close(p%payoff(j)+q%payoff(j), 0.0_dp)
      call assert_close(p%profit(j)+q%profit(j), 0.0_dp)
   end do

   p = bull_call_bls(100.0_dp, 90.0_dp, 110.0_dp, 0.05_dp, 1.0_dp, 0.20_dp)
   q = bear_call_bls(100.0_dp, 90.0_dp, 110.0_dp, 0.05_dp, 1.0_dp, 0.20_dp)
   do j = 1, size(p%stock)
      call assert_close(p%payoff(j)+q%payoff(j), 0.0_dp)
      call assert_close(p%profit(j)+q%profit(j), 0.0_dp)
   end do

   fw = forward_contract(100.0_dp, 1.0_dp, 0.05_dp, 'long', 'none', 0.0_dp, 1.0_dp, 0.0_dp, -1.0_dp)
   call assert_close(fw%delivery_price, 100.0_dp*exp(0.05_dp))
   call assert_close(fw%prepaid_price, 100.0_dp)
   fw = forward_prepaid(100.0_dp, 1.0_dp, 0.05_dp, 'long', 'continuous', 0.0_dp, 1.0_dp, 0.02_dp, -1.0_dp)
   call assert_close(fw%prepaid_price, 100.0_dp*exp(-0.02_dp))

   print '(a)', 'test_derivatives: PASS'
end program test_derivatives
