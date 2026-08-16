program test_regression_tests
   use rfast
   implicit none
   real(dp) :: x(5,2), y(5), xb(8,2), yb(8), xp(6,2), yp(6)
   real(dp) :: tab(2,2)
   type(regression_result) :: fit
   type(test_result) :: tst

   x(:,1)=1.0_dp
   x(:,2)=[0.0_dp,1.0_dp,2.0_dp,3.0_dp,4.0_dp]
   y=2.0_dp+3.0_dp*x(:,2)
   fit=lmfit(x,y)
   call assert_close(fit%beta(1),2.0_dp,1e-11_dp,'lm intercept')
   call assert_close(fit%beta(2),3.0_dp,1e-11_dp,'lm slope')

   xb(:,1)=1.0_dp
   xb(:,2)=[-2.0_dp,-1.5_dp,-1.0_dp,-0.5_dp,0.5_dp,1.0_dp,1.5_dp,2.0_dp]
   yb=[0.0_dp,0.0_dp,0.0_dp,0.0_dp,1.0_dp,1.0_dp,1.0_dp,1.0_dp]
   fit=glm_logistic(xb,yb,maxiter=50)
   call assert_true(fit%beta(2)>0.0_dp,'logistic sign')

   xp(:,1)=1.0_dp
   xp(:,2)=[0.0_dp,1.0_dp,2.0_dp,3.0_dp,4.0_dp,5.0_dp]
   yp=[1.0_dp,1.0_dp,2.0_dp,3.0_dp,4.0_dp,6.0_dp]
   fit=glm_poisson(xp,yp,maxiter=100)
   call assert_true(fit%beta(2)>0.0_dp,'poisson sign')

   tst=ttest1([1.0_dp,2.0_dp,3.0_dp,4.0_dp,5.0_dp],3.0_dp)
   call assert_close(tst%statistic,0.0_dp,1e-14_dp,'ttest center')
   call assert_close(tst%pvalue,1.0_dp,1e-14_dp,'ttest p')
   tab=reshape([10.0_dp,20.0_dp,20.0_dp,40.0_dp],[2,2])
   tst=chi2_test(tab)
   call assert_close(tst%statistic,0.0_dp,1e-14_dp,'chi independence')

   print *, 'test_regression_tests: PASS'
contains
   subroutine assert_close(got,want,tol,msg)
      real(dp),intent(in)::got,want,tol
      character(*),intent(in)::msg
      if(abs(got-want)>tol)then
         print *, 'FAIL ',trim(msg),got,want
         error stop 1
      end if
   end subroutine
   subroutine assert_true(ok,msg)
      logical,intent(in)::ok
      character(*),intent(in)::msg
      if(.not.ok)then
         print *, 'FAIL ',trim(msg)
         error stop 1
      end if
   end subroutine
end program test_regression_tests
