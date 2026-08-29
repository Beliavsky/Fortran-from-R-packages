program test_bats
   use forecast, only : dp, bats_make_w, tbats_make_w, bats_model, forecast_result, bats_fit, tbats_fit, bats_forecast
   implicit none
   real(dp),allocatable::w(:)
   real(dp)::y(48)
   integer::i
   type(bats_model)::model
   type(forecast_result)::fc
   w=bats_make_w(.false.,1.0_dp,[12],[0.3_dp],[0.2_dp])
   call check(size(w)==15,'BATS W size')
   call check(abs(w(1)-1.0_dp)<1e-12_dp .and. abs(w(13)-1.0_dp)<1e-12_dp,'BATS seasonal observation')
   call check(abs(w(14)-0.3_dp)<1e-12_dp .and. abs(w(15)-0.2_dp)<1e-12_dp,'BATS ARMA placement')
   w=tbats_make_w(.true.,0.95_dp,[2],[0.4_dp],[0.1_dp])
   call check(size(w)==8 .and. abs(w(7)-0.4_dp)<1e-12_dp .and. abs(w(8)-0.1_dp)<1e-12_dp,'TBATS W placement')
   do i=1,48
   y(i)=20.0_dp+0.1_dp*i+2.0_dp*sin(2.0_dp*acos(-1.0_dp)*i/12.0_dp)
   end do
   model=bats_fit(y,periods=[12],use_trend=.true.,damped=.false.,optimize=.false.)
   fc=bats_forecast(model,3)
   call check(all(fc%mean==fc%mean),'BATS high-level fit')
   model=tbats_fit(y,[12],[2],use_trend=.true.,damped=.false.,optimize=.false.)
   fc=bats_forecast(model,3)
   call check(all(fc%mean==fc%mean),'TBATS high-level fit')
   print '(a)','test_bats: PASS'
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
