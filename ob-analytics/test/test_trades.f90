! SPDX-License-Identifier: GPL-2.0-or-later
program test_trades
   use ob_analytics
   implicit none
   type(event_t) :: e(4)
   type(trade_t),allocatable::t(:)
   type(impact_t),allocatable::imp(:)

   call make_event(e(1),1,10_i8,1000_i8,900_i8,100.0_dp,5.0_dp,side_bid,2)
   call make_event(e(2),2,20_i8,1010_i8,950_i8,101.0_dp,5.0_dp,side_ask,1)
   call make_event(e(3),3,11_i8,2000_i8,1900_i8,99.0_dp,3.0_dp,side_bid,4)
   call make_event(e(4),4,20_i8,2010_i8,1950_i8,101.0_dp,3.0_dp,side_ask,3)
   t=match_trades(e,correct_large_jumps=.false.)
   call assert_true(size(t)==2,'trade count')
   call assert_close(t(1)%price,100.0_dp,1.0e-14_dp,'maker price')
   call assert_true(t(1)%direction==trade_sell,'sell direction')
   call assert_true(t(1)%maker_id==10_i8 .and. t(1)%taker_id==20_i8,'maker taker')
   call assert_true(t(2)%direction==trade_sell,'second direction')
   imp=trade_impacts(t)
   call assert_true(size(imp)==1,'impact count')
   call assert_true(imp(1)%hits==2,'impact hits')
   call assert_close(imp(1)%volume,8.0_dp,1.0e-14_dp,'impact volume')
   call assert_close(imp(1)%vwap,99.63_dp,1.0e-14_dp,'impact vwap')
   print '(a)', 'test_trades: PASS'
contains
   subroutine make_event(x,eid,id,ts,ex,price,fill,side,match)
      type(event_t),intent(out)::x
      integer,intent(in)::eid,side,match
      integer(i8),intent(in)::id,ts,ex
      real(dp),intent(in)::price,fill
      x%event_id=eid;x%id=id;x%timestamp_ms=ts;x%exchange_timestamp_ms=ex
      x%price=price;x%fill=fill;x%side=side;x%matching_event=match
   end subroutine make_event
   subroutine assert_true(ok,label)
      logical,intent(in)::ok;character(len=*),intent(in)::label
      if(.not.ok) then; print '(a)','failed: '//label; error stop 1; end if
   end subroutine assert_true
   subroutine assert_close(x,y,tol,label)
      real(dp),intent(in)::x,y,tol;character(len=*),intent(in)::label
      if(abs(x-y)>tol) then; print '(a,3es24.16)','mismatch '//label//': ',x,y,abs(x-y); error stop 1; end if
   end subroutine assert_close
end program test_trades
