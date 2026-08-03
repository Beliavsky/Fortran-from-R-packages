! SPDX-License-Identifier: GPL-3.0-or-later
program trade_classification
   use pinstimation
   implicit none
   real(dp) :: timestamp(6),price(6),bid(6),ask(6)
   integer,allocatable :: classification(:)
   integer :: status
   timestamp=[0.0_dp,1.0_dp,2.0_dp,3.0_dp,4.0_dp,5.0_dp]
   price=[10.0_dp,10.1_dp,10.1_dp,10.0_dp,10.2_dp,10.2_dp]
   bid=[9.9_dp,10.0_dp,10.0_dp,9.9_dp,10.1_dp,10.1_dp]
   ask=[10.1_dp,10.2_dp,10.2_dp,10.1_dp,10.3_dp,10.3_dp]
   call classify_trades(timestamp,price,bid,ask,'LR',classification,status=status)
   print '(a,*(i3))', 'LR classes (-1 sell, 0 unresolved, 1 buy): ',classification
end program trade_classification
