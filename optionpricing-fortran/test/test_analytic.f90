! SPDX-License-Identifier: GPL-2.0-only OR GPL-3.0-only
program test_analytic
   use optionpricing, only : dp, moments_result, eval_ecv, eval_lb, eval_eqcv, &
      asian_call_app_lord, covariance_conditional_log_prices
   implicit none
   type(moments_result) :: ecv,lb,eq,tmp
   real(dp), allocatable :: c(:,:)
   real(dp) :: h,p0,pp,pm,fd1,fd2,lord
   integer :: i

   ecv=eval_ecv(1.0_dp,12,100.0_dp,0.05_dp,0.2_dp,100.0_dp)
   lb=eval_lb(1.0_dp,12,100.0_dp,0.05_dp,0.2_dp,100.0_dp,.false.)
   eq=eval_eqcv(1.0_dp,12,100.0_dp,0.05_dp,0.2_dp,100.0_dp)
   call assert_close(ecv%price,6.154908593215828_dp,3.0e-13_dp)
   call assert_close(ecv%delta,0.5938761220547435_dp,3.0e-13_dp)
   call assert_close(ecv%gamma,0.03056755530266579_dp,3.0e-13_dp)
   call assert_close(lb%price,7.050800650973013e-4_dp,2.0e-14_dp)
   call assert_close(eq%price,1.894735078000753e-6_dp,2.0e-16_dp)
   lord=asian_call_app_lord(1.0_dp,12,100.0_dp,0.05_dp,0.2_dp,100.0_dp,.true.)
   call assert_close(lord,6.156040487001411_dp,2.0e-10_dp)

   h=1.0e-2_dp
   p0=ecv%price
   tmp=eval_ecv(1.0_dp,12,100.0_dp,0.05_dp,0.2_dp,100.0_dp+h)
   pp=tmp%price
   tmp=eval_ecv(1.0_dp,12,100.0_dp,0.05_dp,0.2_dp,100.0_dp-h)
   pm=tmp%price
   fd1=(pp-pm)/(2.0_dp*h)
   fd2=(pp-2.0_dp*p0+pm)/(h*h)
   call assert_close(fd1,ecv%delta,3.0e-8_dp)
   call assert_close(fd2,ecv%gamma,3.0e-8_dp)

   c=covariance_conditional_log_prices(1.0_dp,12,0.2_dp)
   if(maxval(abs(c-transpose(c)))>5.0e-15_dp) error stop 1
   if(minval([(c(i,i),i=1,12)]) < -1.0e-14_dp) error stop 1
   print '(a)', 'test_analytic: PASS'
contains
   subroutine assert_close(x,y,tol)
      real(dp), intent(in) :: x,y,tol
      if(abs(x-y)>tol) then
         print '(a,3(es24.16,1x))','mismatch: ',x,y,abs(x-y)
         error stop 1
      end if
   end subroutine assert_close
end program test_analytic
