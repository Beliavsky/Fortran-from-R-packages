! SPDX-License-Identifier: GPL-2.0-only
program test_dynamic
   use streg, only : dp, streg_fit, streg_options, star, stdlm, stvar, streg_success
   use test_support_dynamic
   implicit none
   integer, parameter :: n=100
   real(dp) :: y(n),data(n,2),x(n,1),trend(n,2)
   integer :: i
   type(streg_fit) :: arfit,dlmfit,varfit
   type(streg_options) :: opt
   y(1)=0.0_dp
   do i=2,n
      y(i)=0.4_dp+0.65_dp*y(i-1)+0.25_dp*sin(1.37_dp*real(i,dp))
   end do
   opt%max_iter=300; opt%tolerance=1.0e-7_dp
   arfit=star(y,lag=1,v=8.0_dp,options=opt)
   call assert_true(arfit%status==streg_success .and. arfit%converged,'StAR fit')
   call assert_true(all(shape(arfit%beta)==[1,2]),'StAR beta dimensions')
   call assert_true(arfit%r_squared(1)>0.15_dp,'StAR explanatory power')
   data(1,:)=[0.0_dp,0.0_dp]
   do i=2,n
      data(i,1)=0.2_dp+0.5_dp*data(i-1,1)+0.15_dp*data(i-1,2)+0.12_dp*sin(real(i,dp))
      data(i,2)=-0.1_dp+0.2_dp*data(i-1,1)+0.6_dp*data(i-1,2)+0.10_dp*cos(1.2_dp*real(i,dp))
   end do
   trend(:,1)=1.0_dp
   do i=1,n; trend(i,2)=real(i,dp)/real(n,dp); end do
   varfit=stvar(data,lag=1,trend=trend,v=10.0_dp,options=opt)
   call assert_true(varfit%status==streg_success .and. varfit%converged,'StVAR fit')
   call assert_true(all(shape(varfit%beta)==[2,4]),'StVAR beta dimensions')
   call assert_all_finite(reshape(varfit%fitted,[size(varfit%fitted)]),'finite StVAR fitted values')
   x(:,1)=data(:,2)
   dlmfit=stdlm(data(:,1),x,lag=1,v=10.0_dp,options=opt)
   call assert_true(dlmfit%status==streg_success .and. dlmfit%converged,'StDLM fit')
   call assert_true(all(shape(dlmfit%beta)==[1,4]),'StDLM beta dimensions')
   call assert_all_finite(dlmfit%conditional_factor,'finite StDLM variance factors')
   write(*,'(a)')'test_dynamic: PASS'
end program test_dynamic
