! SPDX-License-Identifier: GPL-2.0-or-later
! This file is part of a modern Fortran translation of robustbase.
! It may be redistributed and/or modified under GPL version 2 or later.
program test_remaining
   use robustbase
   use test_support
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   implicit none
   integer,parameter::n=80
   real(dp)::x(n,2),y(n),trials(n),u,prob
   real(dp)::xn(n,1),yn(n),start(2),lower(2),upper(2)
   real(dp)::cloud(50,2),line_data(20,2),reduced_res(n),ellipse_center(2),ellipse_cov(2,2)
   real(dp),allocatable::ellipse(:,:),fullmat(:,:)
   integer,allocatable::kept(:)
   integer::fullrank
   type(lmrob_result)::lm,sfit,lar
   type(glmrob_result)::mqbin,mqpois,mt
   type(nlrob_method_result)::nmm,ntau,ncm,nmtl
   type(adjusted_outlyingness_result)::ao
   type(partitioned_mcd_result)::pmcd
   type(detmcd_result)::exact_mcd
   type(partitioned_lts_result)::plts
   type(robust_prediction_result)::pred
   type(robust_test_result)::wald,devtest
   type(robust_outlier_stats_result)::ostats
   integer::i

   call seed_rng(1111)
   do i=1,n
      x(i,1)=1.0_dp
      x(i,2)=-2.0_dp+4.0_dp*real(i-1,dp)/real(n-1,dp)
      y(i)=1.0_dp+2.0_dp*x(i,2)+0.08_dp*sin(0.7_dp*real(i,dp))
   end do
   y(n-11:n)=y(n-11:n)+25.0_dp
   call lmrob_s_fit(x,y,sfit,n_resample=250,sampling='nonsingular')
   call lmrob_fit(x,y,lm,method='SMDM',n_resample=250,sampling='nonsingular',covariance_method='sandwich')
   call assert_true(sfit%converged .and. lm%converged,'lmrob S and SMDM convergence')
   call assert_true(abs(lm%coefficients(2)-2.0_dp)<0.15_dp,'lmrob robust slope')
   call assert_true(size(lm%chain_scales)>=3,'lmrob SMDM chain scales')
   call assert_true(all(ieee_is_finite(lm%standard_errors)),'lmrob finite standard errors')
   call lmrob_lar_fit(x,y,lar,max_iter=400)
   call assert_true(lar%converged .and. abs(lar%coefficients(2)-2.0_dp)<0.15_dp,'lmrob LAR robust slope')
   call assert_true(all(ieee_is_finite(lar%standard_errors)),'lmrob LAR inference')
   call robust_linear_predict(x(1:4,:),lm%coefficients,lm%covariance,lm%scale,pred,level=0.95_dp,interval='prediction',df=n-2)
   call assert_true(all(pred%lower<pred%fit) .and. all(pred%upper>pred%fit),'robust prediction intervals')
   call robust_wald_test(lm%coefficients,lm%covariance,[2],wald)
   call assert_true(wald%statistic>10.0_dp .and. wald%p_value<0.01_dp,'robust Wald test')
   call assert_true(robust_r_squared(y,lm%residuals)>0.7_dp,'robust R squared')
   reduced_res=y-sum(y)/real(n,dp)
   call robust_deviance_test(lm%residuals,reduced_res,lm%scale,devtest,degrees_freedom=2)
   call assert_true(devtest%statistic>0.0_dp .and. devtest%p_value<0.05_dp .and. devtest%degrees_freedom==2,'robust deviance test')
   call robust_outlier_stats(x,lm%residuals,lm%scale,ostats)
   call assert_true(ostats%n_outliers>=8,'robust outlier statistics')

   call seed_rng(1444)
   call fast_lts_partitioned(x,y,plts,alpha=0.6_dp,n_partitions=4,n_starts=30,max_csteps=40)
   call assert_true(plts%converged .and. abs(plts%estimate%coefficients(2)-2.0_dp)<0.18_dp,'partitioned FAST-LTS slope')


   call seed_rng(2222)
   do i=1,n
      prob=1.0_dp/(1.0_dp+exp(-(-0.4_dp+1.3_dp*x(i,2))))
      call random_number(u)
      y(i)=merge(1.0_dp,0.0_dp,u<prob)
   end do
   y(n-7:n)=0.0_dp
   call glmrob_mqle_fit(x,y,'binomial',mqbin,max_iter=200)
   call assert_true(mqbin%converged,'glmrob Mqle binomial convergence')
   call assert_true(mqbin%coefficients(2)>0.3_dp,'glmrob Mqle binomial slope')
   call assert_true(all(ieee_is_finite(mqbin%standard_errors)),'glmrob Mqle inference')

   do i=1,n
      prob=exp(0.2_dp+0.35_dp*x(i,2))
      y(i)=real(nint(prob+0.25_dp*sin(real(i,dp))),dp)
      y(i)=max(y(i),0.0_dp)
   end do
   y(n-3:n)=y(n-3:n)+30.0_dp
   call glmrob_mqle_fit(x,y,'poisson',mqpois,max_iter=200)
   call assert_true(mqpois%coefficients(2)>0.05_dp,'glmrob Mqle poisson slope')
   call assert_true(all(ieee_is_finite(mqpois%standard_errors)),'glmrob poisson inference')

   trials=20.0_dp
   do i=1,n
      prob=1.0_dp/(1.0_dp+exp(-(-0.25_dp+0.9_dp*x(i,2))))
      y(i)=real(nint(trials(i)*prob),dp)
   end do
   y(n-5:n)=0.0_dp
   call glmrob_mt_fit(x,y,mt,trials=trials,max_iter=200)
   call assert_true(mt%converged,'glmrob MT convergence')
   call assert_true(mt%coefficients(2)>0.2_dp,'glmrob MT slope')

   do i=1,n
      xn(i,1)=real(i-1,dp)/real(n-1,dp)
      yn(i)=2.0_dp*exp(-0.7_dp*xn(i,1))+0.02_dp*sin(real(i,dp))
   end do
   yn(n-7:n)=yn(n-7:n)+4.0_dp
   start=[1.5_dp,-0.3_dp];lower=[0.1_dp,-3.0_dp];upper=[5.0_dp,1.0_dp]
   call nlrob_mm_fit(exp_model,xn,yn,start,nmm,lower=lower,upper=upper,n_starts=25,max_iter=120)
   call nlrob_tau_fit(exp_model,xn,yn,start,lower,upper,ntau,max_iter=300)
   call nlrob_cm_fit(exp_model,xn,yn,start,lower,upper,ncm,max_iter=300)
   call nlrob_mtl_fit(exp_model,xn,yn,start,lower,upper,nmtl,alpha=0.75_dp,max_iter=300)
   call assert_true(nmm%converged,'nlrob MM convergence')
   call assert_true(abs(nmm%parameters(1)-2.0_dp)<0.25_dp .and. abs(nmm%parameters(2)+0.7_dp)<0.25_dp,'nlrob MM parameters')
   call assert_true(abs(ntau%parameters(1)-2.0_dp)<0.4_dp,'nlrob tau parameters')
   call assert_true(abs(ncm%parameters(1)-2.0_dp)<0.5_dp,'nlrob CM parameters')
   call assert_true(abs(nmtl%parameters(1)-2.0_dp)<0.4_dp,'nlrob MTL parameters')

   do i=1,49
      cloud(i,1)=cos(0.3_dp*real(i,dp))+0.02_dp*real(mod(i,5),dp)
      cloud(i,2)=sin(0.3_dp*real(i,dp))+0.02_dp*real(mod(i,7),dp)
   end do
   cloud(50,:)=[12.0_dp,-11.0_dp]
   call seed_rng(3222)
   call fast_mcd_partitioned(cloud,pmcd,alpha=0.75_dp,n_partitions=4,max_csteps=40)
   call assert_true(pmcd%converged .and. pmcd%candidates>=2,'partitioned FAST-MCD convergence')
   call assert_true(abs(pmcd%estimate%center(1))<1.0_dp .and. abs(pmcd%estimate%center(2))<1.0_dp,'partitioned FAST-MCD robust center')
   call assert_true(pmcd%raw_consistency_factor>0.0_dp .and. pmcd%reweight_consistency_factor>0.0_dp,'MCD correction factors')
   call assert_true(mcd_consistency_factor(2,0.75_dp)>1.0_dp,'MCD consistency formula')
   call assert_true(mcd_finite_sample_factor(2,50,0.75_dp)>1.0_dp,'MCD finite-sample formula')
   do i=1,20
      line_data(i,1)=real(i,dp);line_data(i,2)=1.0_dp+2.0_dp*line_data(i,1)
   end do
   call cov_detmcd(line_data,exact_mcd,alpha=0.75_dp)
   call assert_true(exact_mcd%exact_fit .and. exact_mcd%points_on_hyperplane==20,'detMCD exact hyperplane fit')
   call assert_true(size(exact_mcd%hyperplane_coefficients)==2,'detMCD hyperplane coefficients')
   call full_rank_matrix(reshape([1.0_dp,2.0_dp,3.0_dp,2.0_dp,4.0_dp,6.0_dp],[3,2]),fullmat,kept,fullrank)
   call assert_true(fullrank==1 .and. size(kept)==1,'fullRank column reduction')
   call seed_rng(3333)
   call adjusted_outlyingness_full(cloud,ao,n_directions=80,max_iterations=20000)
   call assert_true(ao%converged .and. ao%directions_used==80,'adjusted outlyingness directions')
   call assert_true(abs(ao%outlyingness(50)-maxval(ao%outlyingness))<1.0e-12_dp,'adjusted outlyingness extreme point')
   call assert_true(.not.ao%non_outlier(50),'adjusted outlyingness cutoff')
   ellipse_center=[0.0_dp,0.0_dp];ellipse_cov=reshape([4.0_dp,1.0_dp,1.0_dp,2.0_dp],[2,2])
   call tolerance_ellipse_points(ellipse_center,ellipse_cov,ellipse,n_points=101)
   call assert_true(all(shape(ellipse)==[101,2]),'tolerance ellipse shape')
   call assert_true(maxval(abs(ellipse(1,:)-ellipse(101,:)))<1.0e-10_dp,'tolerance ellipse closed')
   write(*,'(a)')'lmrob, glmrob, nlrob, partitioned FAST algorithms, and full adjusted-outlyingness tests passed.'
end program test_remaining
