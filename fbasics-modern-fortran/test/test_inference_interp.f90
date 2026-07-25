! SPDX-License-Identifier: GPL-2.0-or-later
! Numerical translation derived from the GPL-2-or-later fBasics package.
program test_inference_interp
  use fbasics
  use test_support
  implicit none
  real(dp)::x(6),y(6),z(6),gridx(3),gridy(2),gridz(3,2),normal_sample(100)
  type(test_result)::r
  integer::i
  x=[0.0_dp,1.0_dp,0.0_dp,1.0_dp,0.5_dp,0.25_dp]
  y=[0.0_dp,0.0_dp,1.0_dp,1.0_dp,0.5_dp,0.75_dp]
  z=2.0_dp+3.0_dp*x-4.0_dp*y
  call assert_close(linear_interp([0.0_dp,1.0_dp],[2.0_dp,4.0_dp],0.25_dp),2.5_dp,1e-14_dp,'linear interpolation')
  gridx=[0.0_dp,1.0_dp,2.0_dp];gridy=[0.0_dp,1.0_dp]
  do i=1,3;gridz(i,1)=gridx(i);gridz(i,2)=gridx(i)+2.0_dp;end do
  call assert_close(bilinear_interp(gridx,gridy,gridz,0.5_dp,0.25_dp),1.0_dp,1e-14_dp,'bilinear interpolation')
  call assert_close(local_plane_interp(x,y,z,0.4_dp,0.3_dp),2.0_dp+1.2_dp-1.2_dp,1e-10_dp,'local plane interpolation')
  r=pearson_test([1.0_dp,2.0_dp,3.0_dp,4.0_dp],[2.0_dp,4.0_dp,6.0_dp,8.0_dp]);call assert_close(r%statistic,1.0_dp,1e-14_dp,'Pearson correlation')
  r=spearman_test([1.0_dp,3.0_dp,2.0_dp,4.0_dp],[10.0_dp,30.0_dp,20.0_dp,40.0_dp]);call assert_close(r%statistic,1.0_dp,1e-14_dp,'Spearman correlation')
  call set_lcg_seed(73_8);do i=1,100;normal_sample(i)=rnorm_lcg();end do
  r=jarque_bera_test(normal_sample);call assert_true(r%statistic>=0.0_dp.and.r%p_value>=0.0_dp.and.r%p_value<=1.0_dp,'Jarque-Bera range')
  r=anderson_darling_normal_test(normal_sample);call assert_true(r%p_value>=0.0_dp.and.r%p_value<=1.0_dp,'AD range')
  r=ks_two_sample_test(normal_sample,normal_sample);call assert_close(r%statistic,0.0_dp,1e-14_dp,'KS identical sample')
  call assert_close(maxdd_expectation(0.0_dp,1.0_dp,100.0_dp),sqrt(pi/2.0_dp)*10.0_dp,1e-12_dp,'zero-drift expected max drawdown')
  call assert_true(dmaxdd(10.0_dp,1.0_dp,100.0_dp,300)>0.0_dp,'max drawdown density')
  call assert_true(pmaxdd(10.0_dp,1.0_dp,100.0_dp,300)>=0.0_dp,'max drawdown probability')
  write(*,'(a)')'Interpolation, inference, and drawdown tests passed.'
end program
