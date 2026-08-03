! SPDX-License-Identifier: GPL-3.0-or-later
program test_data_processing
   use pinstimation
   implicit none
   real(dp) :: timestamp(8),price(8),bid(8),ask(8)
   integer,allocatable :: tick(:),quote(:),lr(:),emo(:)
   integer :: group(8),status
   type(trade_counts) :: counts

   timestamp=[0.0_dp,1.0_dp,2.0_dp,3.0_dp,4.0_dp,5.0_dp,6.0_dp,7.0_dp]
   price=[10.0_dp,10.1_dp,10.1_dp,10.0_dp,10.2_dp,10.2_dp,10.3_dp,10.2_dp]
   bid=[9.9_dp,10.0_dp,10.0_dp,9.9_dp,10.1_dp,10.1_dp,10.2_dp,10.1_dp]
   ask=[10.1_dp,10.2_dp,10.2_dp,10.1_dp,10.3_dp,10.3_dp,10.4_dp,10.3_dp]
   group=[1,1,1,1,2,2,2,2]
   call classify_trades(timestamp,price,bid,ask,'Tick',tick,status=status)
   if(status/=0.or.any(tick/=[0,1,1,-1,1,1,1,-1])) error stop 'tick classification failed'
   call classify_trades(timestamp,price,bid,ask,'Quote',quote,status=status)
   call classify_trades(timestamp,price,bid,ask,'LR',lr,status=status)
   call classify_trades(timestamp,price,bid,ask,'EMO',emo,status=status)
   if(lr(3)/=tick(3).or.lr(6)/=tick(6)) error stop 'LR midpoint fallback failed'
   if(emo(2)/=trade_unresolved .and. emo(2)/=trade_buy .and. emo(2)/=trade_sell) error stop 'invalid EMO class'
   call aggregate_classifications(group,tick,counts,status)
   if(status/=0.or.any(counts%buys/=[2_i8,3_i8]).or.any(counts%sells/=[1_i8,1_i8])) error stop 'aggregation failed'
   print '(a)', 'test_data_processing: PASS'
end program test_data_processing
