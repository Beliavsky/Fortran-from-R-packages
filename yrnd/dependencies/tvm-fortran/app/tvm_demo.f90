! SPDX-License-Identifier: MIT
! Copyright (c) 2014 Juan Manuel Truppia
program tvm_demo
   use tvm
   implicit none
   type(loan_t) :: mortgage
   type(rate_curve_t) :: curve
   real(dp) :: periodic_irr

   mortgage = loan(0.05_dp, 10, 100.0_dp, "french")
   periodic_irr = irr([-100.0_dp, mortgage%cf])
   curve = rate_curve_from_rates([0.04_dp, 0.045_dp, 0.05_dp], "zero_eff")

   print '(a,f12.6)', "Payment: ", mortgage%cf(1)
   print '(a,f12.6)', "IRR:     ", periodic_irr
   print '(a,f12.6)', "PV:      ", curve%present_value(mortgage%cf)
end program tvm_demo
