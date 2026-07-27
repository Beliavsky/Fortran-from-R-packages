! SPDX-License-Identifier: GPL-2.0-or-later
program test_processing
   use ob_analytics
   implicit none
   type(processing_result_t)::r
   integer::status
   character(len=:),allocatable::message
   call process_data('test/events.csv',r,warmup_ms=0_i8,remove_zombies=.false.,status=status,message=message)
   call assert_true(status==0,'process status: '//message)
   call assert_true(size(r%events)==5,'sanitized event count')
   call assert_true(size(r%trades)==1,'inferred trade count')
   call assert_close(r%trades(1)%price,100.0_dp,1.0e-14_dp,'trade price')
   call assert_true(r%trades(1)%direction==trade_sell,'trade direction')
   call assert_true(size(r%depth)==3,'depth count excludes market order')
   call assert_close(r%depth(2)%volume,5.0_dp,1.0e-14_dp,'remaining volume')
   print '(a)', 'test_processing: PASS'
contains
   subroutine assert_true(ok,label)
      logical,intent(in)::ok;character(len=*),intent(in)::label
      if(.not.ok) then; print '(a)','failed: '//label; error stop 1; end if
   end subroutine assert_true
   subroutine assert_close(x,y,tol,label)
      real(dp),intent(in)::x,y,tol;character(len=*),intent(in)::label
      if(abs(x-y)>tol) then; print '(a,3es24.16)','mismatch '//label//': ',x,y,abs(x-y); error stop 1; end if
   end subroutine assert_close
end program test_processing
