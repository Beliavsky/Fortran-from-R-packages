! SPDX-License-Identifier: GPL-2.0-or-later
program demo_fincal
   use fincal
   implicit none
   real(dp), parameter :: cash_flows(4) = [-5.0_dp, 1.6_dp, 2.4_dp, 2.8_dp]
   real(dp) :: rate
   type(inventory_result) :: inventory
   integer :: status

   print '(a,f12.6)', 'future value: ', fv(0.07_dp, 10.0_dp, present_value = 1000.0_dp, payment = 10.0_dp)
   print '(a,f12.6)', 'present value: ', pv(0.07_dp, 10.0_dp, future_value = 1000.0_dp, payment = 10.0_dp)

   rate = irr(cash_flows, status)
   print '(a,f12.6,2x,a)', 'IRR: ', rate, trim(fincal_status_message(status))

   inventory = cogs(2.0_dp, 2.0_dp, [3.0_dp, 5.0_dp], [3.0_dp, 5.0_dp], 7.0_dp, 'FIFO')
   print '(a,f10.2)', 'FIFO cost of goods: ', inventory%cost_of_goods
   print '(a,f10.2)', 'FIFO ending inventory: ', inventory%ending_inventory
end program demo_fincal
