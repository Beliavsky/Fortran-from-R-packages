! SPDX-License-Identifier: GPL-2.0-or-later
program tvm_example
   use fincal
   implicit none
   real(dp) :: rate
   integer :: status

   rate = discount_rate(5.0_dp, 0.0_dp, 600.0_dp, -100.0_dp, status = status)
   print '(a,f12.8)', 'rate per period: ', rate
   print '(a,a)', 'status: ', fincal_status_message(status)
   print '(a,f12.4)', 'payment: ', pmt(0.08_dp, 10.0_dp, -1000.0_dp, 0.0_dp)
end program tvm_example
