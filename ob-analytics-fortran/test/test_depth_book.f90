! SPDX-License-Identifier: GPL-2.0-or-later
program test_depth_book
   use ob_analytics
   implicit none
   type(event_t)::e(6)
   type(depth_update_t),allocatable::d(:),fd(:)
   type(depth_summary_t)::m
   type(spread_t)::sp
   type(order_book_t)::book

   call ev(e(1),1,1_i8,1000_i8,99.0_dp,10.0_dp,0.0_dp,action_created,side_bid)
   call ev(e(2),2,2_i8,2000_i8,101.0_dp,12.0_dp,0.0_dp,action_created,side_ask)
   call ev(e(3),3,3_i8,3000_i8,98.0_dp,4.0_dp,0.0_dp,action_created,side_bid)
   call ev(e(4),4,2_i8,4000_i8,101.0_dp,7.0_dp,5.0_dp,action_changed,side_ask)
   call ev(e(5),5,1_i8,5000_i8,99.0_dp,10.0_dp,0.0_dp,action_deleted,side_bid)
   call ev(e(6),6,2_i8,6000_i8,101.0_dp,7.0_dp,0.0_dp,action_deleted,side_ask)
   d=price_level_volume(e)
   call assert_true(size(d)==6,'depth size')
   call assert_close(d(1)%volume,10.0_dp,1.0e-14_dp,'first depth')
   call assert_close(d(4)%volume,7.0_dp,1.0e-14_dp,'partial fill')
   call assert_close(d(5)%volume,0.0_dp,1.0e-14_dp,'bid close')
   m=depth_metrics(d,bps=100,bins=2)
   call assert_close(m%best_bid_price(2),99.0_dp,1.0e-14_dp,'best bid')
   call assert_close(m%best_ask_price(2),101.0_dp,1.0e-14_dp,'best ask')
   call assert_close(m%best_ask_volume(4),7.0_dp,1.0e-14_dp,'best ask volume')
   sp=get_spread(m)
   call assert_true(size(sp%timestamp_ms)>=4,'spread changes')
   fd=filter_depth(d,2000_i8,5000_i8)
   call assert_true(any(fd%timestamp_ms==2000_i8),'filter start')
   call assert_true(any(fd%timestamp_ms==5000_i8),'filter end')
   book=reconstruct_order_book(e,4500_i8)
   call assert_true(size(book%bids)==2,'book bids')
   call assert_true(size(book%asks)==1,'book asks')
   call assert_close(book%asks(1)%volume,7.0_dp,1.0e-14_dp,'book latest change')
   call assert_close(book%bids(1)%price,99.0_dp,1.0e-14_dp,'book best bid')
   print '(a)', 'test_depth_book: PASS'
contains
   subroutine ev(x,eid,id,ts,price,volume,fill,action,side)
      type(event_t),intent(out)::x
      integer,intent(in)::eid,action,side
      integer(i8),intent(in)::id,ts
      real(dp),intent(in)::price,volume,fill
      x%event_id=eid;x%id=id;x%timestamp_ms=ts;x%exchange_timestamp_ms=ts
      x%price=price;x%volume=volume;x%fill=fill;x%action=action;x%side=side
      x%order_type=type_resting_limit
   end subroutine ev
   subroutine assert_true(ok,label)
      logical,intent(in)::ok;character(len=*),intent(in)::label
      if(.not.ok) then; print '(a)','failed: '//label; error stop 1; end if
   end subroutine assert_true
   subroutine assert_close(x,y,tol,label)
      real(dp),intent(in)::x,y,tol;character(len=*),intent(in)::label
      if(abs(x-y)>tol) then; print '(a,3es24.16)','mismatch '//label//': ',x,y,abs(x-y); error stop 1; end if
   end subroutine assert_close
end program test_depth_book
