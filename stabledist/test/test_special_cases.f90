program test_special_cases
   use r_compat, only: dp, normal_cdf, dnorm, dcauchy, pcauchy
   use stabledist
   implicit none
   real(dp), parameter :: tol=2.0e-8_dp
   real(dp) :: x, refd, refp

   x=0.7_dp
   call check_close(dstable(x,2.0_dp,0.8_dp),dnorm(x,mean=0.0_dp,sd=sqrt(2.0_dp)),1.0e-13_dp,'alpha=2 density')
   call check_close(pstable(x,2.0_dp,-0.9_dp),normal_cdf(x/sqrt(2.0_dp)),1.0e-13_dp,'alpha=2 cdf')
   call check_close(dstable(x,1.0_dp,0.0_dp),dcauchy(x),1.0e-13_dp,'Cauchy density')
   call check_close(pstable(x,1.0_dp,0.0_dp),pcauchy(x),1.0e-13_dp,'Cauchy cdf')

   ! Levy(mu=0,c=1) == Stable(alpha=1/2,beta=1,gamma=1,delta=0,pm=1).
   x=2.0_dp
   refd=sqrt(1.0_dp/(2.0_dp*acos(-1.0_dp)))*exp(-1.0_dp/(2.0_dp*x))/x**1.5_dp
   refp=erfc(sqrt(1.0_dp/(2.0_dp*x)))
   call check_close(dstable(x,0.5_dp,1.0_dp,gamma=1.0_dp,delta=0.0_dp,pm=1),refd,2.0e-10_dp,'Levy density')
   call check_close(pstable(x,0.5_dp,1.0_dp,gamma=1.0_dp,delta=0.0_dp,pm=1),refp,2.0e-9_dp,'Levy cdf')
   if (pstable(-1.1_dp,0.5_dp,1.0_dp) /= 0.0_dp) error stop 'finite support alpha=.5 beta=1'
   if (pstable(-2.1_dp,0.6_dp,1.0_dp) /= 0.0_dp) error stop 'finite support alpha=.6 beta=1'

   call check_close(c_stable_tail(0.0_dp),0.5_dp,0.0_dp,'tail constant alpha=0')
   call check_close(c_stable_tail(1.0_dp),1.0_dp/acos(-1.0_dp),1.0e-15_dp,'tail constant alpha=1')
   call check_close(c_stable_tail(2.0_dp),0.0_dp,0.0_dp,'tail constant alpha=2')

   print *, 'test_special_cases: PASS'
contains
   subroutine check_close(a,b,atol,msg)
      real(dp),intent(in)::a,b,atol
      character(len=*),intent(in)::msg
      if(abs(a-b)>atol*max(1.0_dp,abs(b)))then
         print *,trim(msg),a,b,abs(a-b)
         error stop 1
      end if
   end subroutine check_close
end program test_special_cases
