program test_core
   use forecast, only : dp, boxcox, inv_boxcox, fourier_terms, mean_forecast, random_walk_forecast, &
      seasonal_naive_forecast, fit_croston, forecast_croston, accuracy, dm_test, forecast_result, croston_fit
   implicit none
   real(dp) :: x(8), z(8), stat, pval
   real(dp), allocatable :: f(:,:)
   type(forecast_result) :: fc
   type(croston_fit) :: cr
   x=[1.0_dp,2.0_dp,3.0_dp,4.0_dp,5.0_dp,6.0_dp,7.0_dp,8.0_dp]
   z=inv_boxcox(boxcox(x,0.3_dp),0.3_dp)
   call check(maxval(abs(z-x))<1.0e-11_dp,'Box-Cox inversion')
   f=fourier_terms(12,12,2)
   call check(size(f,1)==12 .and. size(f,2)==4,'Fourier dimensions')
   fc=mean_forecast(x,3)
   call check(maxval(abs(fc%mean-4.5_dp))<1.0e-12_dp,'mean forecast')
   fc=random_walk_forecast(x,3)
   call check(all(abs(fc%mean-8.0_dp)<1.0e-12_dp),'naive forecast')
   fc=seasonal_naive_forecast(x,4,4)
   call check(maxval(abs(fc%mean-[5.0_dp,6.0_dp,7.0_dp,8.0_dp]))<1.0e-12_dp,'seasonal naive')
   cr=fit_croston([0.0_dp,2.0_dp,0.0_dp,0.0_dp,4.0_dp,0.0_dp],alpha=0.2_dp)
   fc=forecast_croston(cr,3)
   call check(all(fc%mean>=0.0_dp),'Croston positive forecast')
   call dm_test([1.0_dp,2.0_dp,1.0_dp,2.0_dp],[2.0_dp,3.0_dp,2.0_dp,3.0_dp],1,2,.true.,stat,pval)
   call check(pval>=0.0_dp .and. pval<=1.0_dp,'DM p-value')
   print '(a)','test_core: PASS'
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
