! SPDX-License-Identifier: GPL-2.0-only
! Based on FinancialMath 0.1.1, Copyright (C) 2016 Kameron Penn and Jack Schmidt.
program financialmath_demo
   use financialmath
   implicit none
   type(tvm_result_t) :: tv
   type(annuity_result_t) :: ann
   type(bond_result_t) :: b
   type(option_order1_t) :: bs

   tv = solve_tvm(1000.0_dp, 0.0_dp, 10.0_dp, 0.05_dp, 1.0_dp, 'fv')
   ann = annuity_level(0.0_dp, 0.0_dp, 20.0_dp, 100.0_dp, 0.05_dp, &
      1.0_dp, 1.0_dp, .true., 'pv')
   b = bond(1000.0_dp, 0.05_dp, 1000.0_dp, 10.0_dp, 0.04_dp, &
      1.0_dp, 1.0_dp, 0.0_dp, .false.)
   bs = bls_order1(100.0_dp, 100.0_dp, 0.05_dp, 1.0_dp, 0.20_dp, 0.0_dp)

   print '(a,f12.4)', 'FV of 1000 after 10 years: ', tv%future_value
   print '(a,f12.4)', 'PV of 20 payments of 100:   ', ann%present_value
   print '(a,f12.4)', 'Bond price:                  ', b%price
   print '(a,f12.4)', 'Black-Scholes call:          ', bs%call_price
end program financialmath_demo
