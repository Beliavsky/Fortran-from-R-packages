program test_cv_bagging
   use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
   use forecast, only : dp, forecast_result, ts_cv, mean_forecast, bagged_ets_forecast
   implicit none
   real(dp)::y(36)
   real(dp),allocatable::err(:,:)
   type(forecast_result)::fc
   integer::i
   do i=1,36
   y(i)=10.0_dp+sin(0.3_dp*i)
   end do
   err=ts_cv(y,mean_cb,h=2,initial=5)
   call check(size(err,1)==36 .and. size(err,2)==2,'tsCV dimensions')
   call check(.not.ieee_is_nan(err(10,1)),'tsCV populated errors')
   fc=bagged_ets_forecast(y,3,num_boot=2,period=1)
   call check(size(fc%mean)==3 .and. all(fc%mean==fc%mean),'bagged ETS')
   print '(a)','test_cv_bagging: PASS'
contains
   function mean_cb(x,h) result(out)
      real(dp),intent(in)::x(:)
      integer,intent(in)::h
      type(forecast_result)::out
      out=mean_forecast(x,h)
   end function
   subroutine check(ok,msg)
      logical,intent(in)::ok
      character(*),intent(in)::msg
      if(.not.ok)then
      write(*,'(a)')'FAIL: '//trim(msg)
      error stop 1
      end if
   end subroutine
end program
