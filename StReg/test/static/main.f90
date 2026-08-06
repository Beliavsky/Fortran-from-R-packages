! SPDX-License-Identifier: GPL-2.0-only
program test_static
   use streg, only : dp, streg_fit, streg_options, stlm, streg_success
   use test_support_static
   implicit none
   integer, parameter :: n=80
   real(dp) :: x(n,1),y(n)
   integer :: i
   type(streg_fit) :: fit
   type(streg_options) :: opt
   do i=1,n
      x(i,1)=-2.0_dp+4.0_dp*real(i-1,dp)/real(n-1,dp)
      y(i)=1.2_dp+2.3_dp*x(i,1)+0.25_dp*sin(1.7_dp*real(i,dp))+0.1_dp*cos(0.3_dp*real(i,dp))
   end do
   opt%max_iter=200; opt%tolerance=1.0e-7_dp
   fit=stlm(y,x,v=8.0_dp,options=opt)
   call assert_true(fit%status==streg_success,'StLM status')
   call assert_true(fit%converged,'StLM convergence')
   call assert_close(fit%beta(1,1),1.2_dp,0.08_dp,'StLM intercept recovery')
   call assert_close(fit%beta(1,2),2.3_dp,0.05_dp,'StLM slope recovery')
   call assert_true(fit%r_squared(1)>0.98_dp,'StLM high R squared')
   call assert_true(fit%conditional_factor(1)>=1.0_dp,'conditional factor lower bound')
   call assert_true(fit%ad_test(1,2)>=0.0_dp .and. fit%ad_test(1,2)<=1.0_dp,'AD p value range')
   call assert_all_finite(fit%residuals(:,1),'finite StLM residuals')
   write(*,'(a)')'test_static: PASS'
end program test_static
