! SPDX-License-Identifier: GPL-2.0-only
program test_conditional
   use streg, only : dp, conditional_parameters, streg_success
   use test_support_conditional
   implicit none
   real(dp) :: theta(11)
   real(dp), allocatable :: s(:,:),b(:,:),sigma(:,:),q(:,:),delta0(:,:),varcoef(:)
   integer :: status
   theta=0.0_dp
   theta(7:9)=[2.0_dp,0.5_dp,1.0_dp]
   theta(10:11)=[1.0_dp,2.0_dp]
   call conditional_parameters(theta,2,0,1,8.0_dp,1,s,b,sigma,q,delta0,varcoef,status)
   call assert_true(status==streg_success,'conditional parameter status')
   call assert_close(s(1,1),4.25_dp,1.0e-12_dp,'joint scale 11')
   call assert_close(s(1,2),1.5_dp,1.0e-12_dp,'joint scale 12')
   call assert_close(s(2,2),1.25_dp,1.0e-12_dp,'joint scale 22')
   call assert_close(b(1,1),1.2_dp,1.0e-12_dp,'conditional slope')
   call assert_close(sigma(1,1),2.45_dp,1.0e-12_dp,'conditional innovation scale')
   call assert_close(q(1,1),0.1_dp,1.0e-12_dp,'scaled predictor precision')
   call assert_close(delta0(1,1),-1.4_dp,1.0e-12_dp,'conditional intercept')
   call assert_close(varcoef(1),2.8_dp,1.0e-12_dp,'variance constant')
   call assert_close(varcoef(2),0.28_dp,1.0e-12_dp,'variance quadratic coefficient')
   write(*,'(a)')'test_conditional: PASS'
end program test_conditional
