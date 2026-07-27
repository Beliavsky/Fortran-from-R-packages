! SPDX-License-Identifier: GPL-2.0-or-later
program test_order_types
   use ob_analytics
   implicit none
   type(event_t)::e(8)
   type(trade_t)::t(2)
   integer::unidentified

   call ev(e(1),1,1_i8,1000_i8,100.0_dp,10.0_dp,action_created)
   call ev(e(2),2,1_i8,2000_i8,100.0_dp,10.0_dp,action_deleted)
   call ev(e(3),3,2_i8,1000_i8,99.0_dp,8.0_dp,action_created)
   call ev(e(4),4,3_i8,1000_i8,101.0_dp,6.0_dp,action_created)
   call ev(e(5),5,3_i8,1500_i8,102.0_dp,4.0_dp,action_changed)
   call ev(e(6),6,4_i8,1100_i8,103.0_dp,5.0_dp,action_created)
   call ev(e(7),7,5_i8,1200_i8,98.0_dp,7.0_dp,action_created)
   call ev(e(8),8,5_i8,1600_i8,98.0_dp,3.0_dp,action_changed)

   t(1)%maker_event_id=7; t(1)%taker_event_id=6
   t(2)%maker_event_id=3; t(2)%taker_event_id=8
   call set_order_types(e,t,unidentified)
   call assert_true(e(1)%order_type==type_flashed_limit,'flashed')
   call assert_true(e(3)%order_type==type_resting_limit,'resting maker')
   call assert_true(e(4)%order_type==type_pacman .and. e(5)%order_type==type_pacman,'pacman')
   call assert_true(e(6)%order_type==type_market,'market')
   call assert_true(e(7)%order_type==type_market_limit .and. e(8)%order_type==type_market_limit,'market limit')
   call assert_true(unidentified==0,'all identified')
   print '(a)', 'test_order_types: PASS'
contains
   subroutine ev(x,eid,id,ts,price,volume,action)
      type(event_t),intent(out)::x
      integer,intent(in)::eid,action
      integer(i8),intent(in)::id,ts
      real(dp),intent(in)::price,volume
      x%event_id=eid;x%id=id;x%timestamp_ms=ts;x%exchange_timestamp_ms=ts
      x%price=price;x%volume=volume;x%action=action;x%side=merge(side_bid,side_ask,mod(eid,2)==1)
   end subroutine ev
   subroutine assert_true(ok,label)
      logical,intent(in)::ok;character(len=*),intent(in)::label
      if(.not.ok) then; print '(a)','failed: '//label; error stop 1; end if
   end subroutine assert_true
end program test_order_types
