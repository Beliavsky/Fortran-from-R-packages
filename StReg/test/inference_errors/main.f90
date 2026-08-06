! SPDX-License-Identifier: GPL-2.0-only
program test_inference_errors
   use streg, only : dp, streg_fit, streg_options, stlm, star, streg_success, streg_invalid_input
   use test_support_inference_errors
   implicit none
   integer, parameter :: n=50
   real(dp) :: x(n,1),y(n),badx(n-1,1)
   integer :: i
   type(streg_fit) :: fit,bad
   type(streg_options) :: opt
   do i=1,n
      x(i,1)=real(i,dp)/10.0_dp
      y(i)=2.0_dp+1.5_dp*x(i,1)+0.2_dp*sin(real(i,dp))
   end do
   opt%max_iter=120; opt%compute_hessian=.true.
   fit=stlm(y,x,v=8.0_dp,options=opt)
   call assert_true(fit%status==streg_success,'inference fit status')
   call assert_true(size(fit%coef_se)==3,'coefficient standard-error count')
   call assert_true(size(fit%var_coef_se)==2,'variance standard-error count')
   call assert_all_finite(fit%coef_se,'finite coefficient standard errors')
   call assert_all_finite(fit%var_coef_se,'finite variance standard errors')
   bad=stlm(y,badx,v=8.0_dp)
   call assert_true(bad%status==streg_invalid_input,'row mismatch error')
   bad=star(y,lag=0,v=8.0_dp)
   call assert_true(bad%status==streg_invalid_input,'invalid lag error')
   write(*,'(a)')'test_inference_errors: PASS'
end program test_inference_errors
