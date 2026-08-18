program test_regression
   use hermite
   implicit none
   integer, parameter :: n=700
   integer(i64) :: y(n),z(1)
   real(dp) :: x(n,2),beta_true(2),mu,a,b
   type(hermite_glm_result) :: fit,auto
   integer :: i,fails

   fails=0
   beta_true=[0.3_dp,0.55_dp]
   do i=1,n
      x(i,1)=1.0_dp
      x(i,2)=-1.0_dp+2.0_dp*real(i-1,dp)/real(n-1,dp)
   end do

   call set_hermite_seed(24680)
   do i=1,n
      mu=exp(dot_product(x(i,:),beta_true))
      b=mu*(2.6_dp-1.0_dp)/real(3*2,dp)
      a=mu-3.0_dp*b
      call rhermite(z,a,b,3)
      y(i)=z(1)
   end do

   fit=fit_glm_hermite(y,x,link=HERMITE_LINK_LOG,m=3)
   if (fit%status>1) fails=fails+1
   if (abs(fit%beta(1)-beta_true(1)) > 0.18_dp) fails=fails+1
   if (abs(fit%beta(2)-beta_true(2)) > 0.18_dp) fails=fails+1
   if (abs(fit%dispersion-2.6_dp) > 0.35_dp) fails=fails+1
   if (.not. allocated(fit%fitted)) fails=fails+1

   auto=fit_glm_hermite(y,x,link=HERMITE_LINK_LOG,max_order=5)
   if (auto%order < 2 .or. auto%order > 5) fails=fails+1
   if (auto%loglik+1.0e-6_dp < fit%loglik) fails=fails+1

   if (fails/=0) then
      print *,'test_regression: FAIL',fails
      print *,'fixed:',fit%beta,fit%dispersion,fit%loglik,fit%status
      print *,'auto:',auto%order,auto%dispersion,auto%loglik
      error stop 1
   end if
   print *,'test_regression: PASS'
end program test_regression
