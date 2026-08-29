program test_dependencies
   use forecast, only : dp, arfima_model, nnetar_model, forecast_result, arfima_fit, arfima_forecast, &
      nnetar_fit, nnetar_forecast, ndiffs
   implicit none
   real(dp) :: y(80)
   integer :: i,d
   type(arfima_model) :: af
   type(nnetar_model) :: nn
   type(forecast_result) :: fc
   do i=1,80
   y(i)=0.6_dp*sin(0.15_dp*real(i,dp))+0.2_dp*cos(0.41_dp*real(i,dp))
   end do
   d=ndiffs(y,2)
   call check(d>=0 .and. d<=2,'urca-backed ndiffs')
   af=arfima_fit(y,max_p=1,max_q=1)
   fc=arfima_forecast(af,4)
   call check(size(fc%mean)==4 .and. all(fc%mean==fc%mean),'fracdiff-backed ARFIMA')
   nn=nnetar_fit(y,m=1,p=3,Pseason=0,size_hidden=2,repeats=1,decay=0.0_dp,maxit=40)
   fc=nnetar_forecast(nn,3)
   call check(size(fc%mean)==3 .and. all(fc%mean==fc%mean),'nnet-backed nnetar')
   print '(a)','test_dependencies: PASS'
contains
   subroutine check(ok,msg)
      logical,intent(in)::ok
      character(*),intent(in)::msg
      if(.not.ok)then
      write(*,'(a)')'FAIL: '//trim(msg)
      error stop 1
      end if
   end subroutine
end program
