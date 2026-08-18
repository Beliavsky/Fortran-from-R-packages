program test_poisson_regression
   use hermite
   implicit none
   integer, parameter :: n=500
   integer(i64) :: y(n),z(1)
   real(dp) :: x(n,2),beta_true(2),mu
   type(hermite_glm_result) :: fit
   integer :: i,fails

   fails=0
   beta_true=[0.15_dp,-0.4_dp]
   call set_hermite_seed(97531)
   do i=1,n
      x(i,1)=1.0_dp
      x(i,2)=-0.8_dp+1.6_dp*real(i-1,dp)/real(n-1,dp)
      mu=exp(dot_product(x(i,:),beta_true))
      call rhermite(z,mu,0.0_dp,2)
      y(i)=z(1)
   end do

   fit=fit_glm_hermite(y,x,link=HERMITE_LINK_LOG,m=1)
   if (fit%status>1) fails=fails+1
   if (abs(fit%beta(1)-beta_true(1)) > 0.16_dp) fails=fails+1
   if (abs(fit%beta(2)-beta_true(2)) > 0.20_dp) fails=fails+1
   if (fit%order /= 1) fails=fails+1
   if (abs(fit%dispersion-1.0_dp) > 1.0e-14_dp) fails=fails+1

   if (fails/=0) then
      print *,'test_poisson_regression: FAIL',fails,fit%beta
      error stop 1
   end if
   print *,'test_poisson_regression: PASS'
end program test_poisson_regression
