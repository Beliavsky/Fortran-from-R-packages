! SPDX-License-Identifier: GPL-2.0-or-later
! This file is part of a modern Fortran translation of robustbase.
! It may be redistributed and/or modified under GPL version 2 or later.
program test_scales
   use robustbase
   use test_support
   implicit none
   real(dp)::x(5),mu,s,se,lf,uf,mc,q3,h,fd,b_lqq,c_lqq,s_lqq
   real(dp)::wx(3),ww(3),y(5)
   integer::it
   logical::out(10)
   real(dp)::x3(16),x4(10),mat(3,3),rm(3),cm(3),z(3,3),cen(3),sca(3)
   logical::keep(3)
   integer::rank
   x=[1.0_dp,2.0_dp,3.0_dp,4.0_dp,5.0_dp]
   call assert_close(median(x),3.0_dp,1.0e-14_dp,'median')
   call assert_close(qn_scale(x),1.8729763514_dp,1.0e-10_dp,'Qn')
   call assert_close(sn_scale(x),1.6112026_dp,1.0e-10_dp,'Sn')
   call assert_close(mad_scale(x),1.482602218505602_dp,1.0e-12_dp,'MAD')
   wx=[1.0_dp,4.0_dp,9.0_dp];ww=[2.0_dp,1.0_dp,3.0_dp]
   call assert_close(weighted_high_median(wx,ww),9.0_dp,0.0_dp,'weighted high median')
   call huber_location([1.0_dp,2.0_dp,3.0_dp,4.0_dp,100.0_dp],mu,s,it,k=1.5_dp,standard_error=se)
   call assert_true(mu<10.0_dp .and. mu>2.0_dp,'Huber location resists outlier')
   call assert_true(s>0.0_dp .and. se>0.0_dp .and. it>0,'Huber outputs')
   call assert_close(huber_psi(3.0_dp,1.5_dp),1.5_dp,0.0_dp,'Huber psi')
   call assert_close(tukey_weight(5.0_dp,4.685_dp),0.0_dp,0.0_dp,'Tukey rejection')
   h=1.0e-6_dp
   fd=(welsh_psi(1.3_dp+h,2.11_dp)-welsh_psi(1.3_dp-h,2.11_dp))/(2.0_dp*h)
   call assert_close(welsh_psi_derivative(1.3_dp,2.11_dp),fd,2.0e-7_dp,'Welsh psi derivative')
   call assert_close(welsh_weight(1.3_dp,2.11_dp),welsh_psi(1.3_dp,2.11_dp)/1.3_dp,1.0e-14_dp,'Welsh weight identity')
   fd=(optimal_psi(2.5_dp+h,1.060158_dp)-optimal_psi(2.5_dp-h,1.060158_dp))/(2.0_dp*h)
   call assert_close(optimal_psi_derivative(2.5_dp,1.060158_dp),fd,2.0e-6_dp,'optimal psi derivative')
   call assert_close(optimal_rho(4.0_dp,1.060158_dp),1.0_dp,0.0_dp,'optimal rho saturation')
   fd=(ggw_psi(2.0_dp+h,0.648_dp,1.0_dp,1.694_dp)-ggw_psi(2.0_dp-h,0.648_dp,1.0_dp,1.694_dp))/(2.0_dp*h)
   call assert_close(ggw_psi_derivative(2.0_dp,0.648_dp,1.0_dp,1.694_dp),fd,2.0e-6_dp,'GGW psi derivative')
   call assert_true(ggw_rho(100.0_dp,0.648_dp,1.0_dp,1.694_dp)>0.999_dp,'GGW rho saturation')
   b_lqq=1.4734061_dp;c_lqq=0.9822707_dp;s_lqq=1.5_dp
   fd=(lqq_psi(1.2_dp+h,b_lqq,c_lqq,s_lqq)-lqq_psi(1.2_dp-h,b_lqq,c_lqq,s_lqq))/(2.0_dp*h)
   call assert_close(lqq_psi_derivative(1.2_dp,b_lqq,c_lqq,s_lqq),fd,2.0e-6_dp,'LQQ psi derivative')
   call assert_close(lqq_weight(1.2_dp,b_lqq,c_lqq,s_lqq),lqq_psi(1.2_dp,b_lqq,c_lqq,s_lqq)/1.2_dp,1.0e-12_dp,'LQQ weight identity')
   call assert_true(lqq_rho(100.0_dp,b_lqq,c_lqq,s_lqq)>0.999_dp,'LQQ rho saturation')
   x3=[-2.0_dp,-1.0_dp,-1.0_dp,-1.0_dp,-1.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,2.0_dp,2.0_dp,2.0_dp,3.0_dp,4.0_dp]
   x4=[1.0_dp,2.0_dp,3.0_dp,4.0_dp,5.0_dp,7.0_dp,10.0_dp,15.0_dp,25.0_dp,1.0e15_dp]
   call assert_close(medcouple(x3),1.0_dp/3.0_dp,1.0e-12_dp,'medcouple ties')
   call assert_close(medcouple(x4),7.0_dp/12.0_dp,1.0e-12_dp,'medcouple skew')
   call adjusted_boxplot_stats(x4,lf,uf,mc,out)
   q3=quantile_type7(x4,0.75_dp)
   call assert_true(mc>0.0_dp .and. uf>q3,'adjusted boxplot')
   call huberize_vector([1.0_dp,2.0_dp,3.0_dp,4.0_dp,100.0_dp],y)
   call assert_true(y(5)<100.0_dp,'huberize clips outlier')
   mat=reshape([1.0_dp,2.0_dp,3.0_dp,4.0_dp,5.0_dp,6.0_dp,2.0_dp,4.0_dp,6.0_dp],[3,3])
   call row_medians(mat,rm);call column_medians(mat,cm)
   call assert_true(all(rm>0.0_dp) .and. all(cm>0.0_dp),'row and column medians')
   call robust_standardize(mat,z,cen,sca);call assert_true(all(sca>0.0_dp),'robust standardization')
   call independent_columns(mat,keep,rank);call assert_true(rank==2,'independent column rank')
   write(*,'(a)')'Scale, score, medcouple, adjusted-boxplot, and utility tests passed.'
end program test_scales
