! SPDX-License-Identifier: MIT
! Copyright (c) 2014 Juan Manuel Truppia
program basic_cashflows
   use tvm
   implicit none
   type(loan_t) :: l

   l = loan(0.05_dp, 10, 100.0_dp, "french")
   print '(a,f12.6)', "Payment = ", pmt(100.0_dp, 10, 0.05_dp)
   print '(a,f12.6)', "NPV = ", npv(0.01_dp, [-1.0_dp, 0.5_dp, 0.9_dp], [0.0_dp, 1.0_dp, 3.0_dp])
   print '(a,f12.6)', "IRR = ", irr([-1.0_dp, 0.5_dp, 0.9_dp], [0.0_dp, 1.0_dp, 3.0_dp])
   print '(a,*(f10.4,1x))', "Loan cashflows: ", l%cf
end program basic_cashflows
