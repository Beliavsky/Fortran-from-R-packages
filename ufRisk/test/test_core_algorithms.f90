! SPDX-License-Identifier: GPL-3.0-only
program test_core_algorithms
   use ufrisk
   use ufrisk_math, only : normal_cdf, student_cdf, student_quantile
   use test_support
   implicit none
   real(dp), parameter :: expected_filter(8) = [ &
      0.7200000000000001_dp,-0.13295_dp,0.044933000000000015_dp, &
      0.00639859125_dp,0.01138649011_dp,0.008334609669612499_dp, &
      0.007398679208340751_dp,0.006397312689299338_dp]
   real(dp) :: ar(1),ma(1),loss(10),var(10),es(10),q
   real(dp), allocatable :: filter(:)
   type(ufrisk_coverage_result) :: coverage
   type(ufrisk_loss_result) :: loss_result
   integer :: i

   call assert_close(normal_cdf(0.0_dp),0.5_dp,1.0e-14_dp,'normal cdf at zero')
   q = student_quantile(0.975_dp,7.0_dp)
   call assert_close(student_cdf(q,7.0_dp),0.975_dp,2.0e-12_dp,'student inversion')

   ar = 0.35_dp; ma = -0.2_dp
   filter = arfilt_coefficients(ar,ma,0.17_dp,8)
   do i = 1,8
      call assert_close(filter(i),expected_filter(i),2.0e-13_dp,'arfilt coefficient')
   end do

   loss = [0.01_dp,0.03_dp,-0.01_dp,0.025_dp,0.06_dp,0.0_dp,0.04_dp,0.08_dp,0.015_dp,0.05_dp]
   var = 0.02_dp
   coverage = covtest(loss,var,0.2_dp)
   call assert_true(coverage%status==ufrisk_ok,'coverage status')
   call assert_true(all([coverage%n00,coverage%n01,coverage%n10,coverage%n11]==[0,4,3,2]), &
      'coverage transition counts')
   call assert_close(coverage%lr_unconditional,7.6381700195377515_dp,2.0e-12_dp,'LR uc')
   call assert_close(coverage%lr_independence,6.730116670092564_dp,2.0e-12_dp,'LR ind')
   call assert_close(coverage%lr_conditional,14.368286689630317_dp,2.0e-12_dp,'LR cc')
   call assert_close(coverage%p_unconditional,0.00571458640362744_dp,2.0e-11_dp,'p uc')
   call assert_close(coverage%p_independence,0.00947983990870854_dp,2.0e-11_dp,'p ind')
   call assert_close(coverage%p_conditional,0.000758518518518519_dp,2.0e-11_dp,'p cc')

   es = [0.025_dp,0.025_dp,0.025_dp,0.025_dp,0.05_dp,0.025_dp,0.035_dp,0.07_dp,0.025_dp,0.045_dp]
   loss_result = lossfunc(loss,es)
   call assert_close(loss_result%regulatory,2.75_dp,1.0e-13_dp,'regulatory loss')
   call assert_close(loss_result%firm,2.875_dp,1.0e-13_dp,'firm loss')
   call assert_close(loss_result%abad,2.835_dp,1.0e-13_dp,'Abad loss')
   call assert_close(loss_result%feng,2.825_dp,1.0e-13_dp,'Feng loss')
   call finish_tests('test_core_algorithms')
end program test_core_algorithms
