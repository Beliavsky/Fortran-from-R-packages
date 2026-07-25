! SPDX-License-Identifier: GPL-2.0-or-later
! Copyright (C) 2013-2026 the original stochvol and factorstochvol authors
! and the modern Fortran translation contributors.
! This file is free software: you may redistribute it and/or modify it
! under the terms of the GNU General Public License as published by the
! Free Software Foundation, either version 2, or (at your option) any later version.
! This file is distributed without any warranty; see LICENSE for details.
program test_factor
  use sv_kinds,only:dp
  use sv_rng,only:seed_rng
  use fsv_types
  use fsv_core
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  implicit none
  real(dp)::b(3,1),ip(3,3),fp(1,2),fp0(0,2),b0(3,0),covm(3,3),corm(3,3),ew(3,3),dens(2)
  real(dp),allocatable::cov(:,:,:),cor(:,:,:),l0(:,:),f0(:,:),u0(:),ll(:,:),llwb(:,:),scores(:),means_pred(:,:,:),vols_pred(:,:,:),evals(:,:),bdense(:,:),bsparse(:,:),ipdef(:,:),fpdef(:,:)
  real(dp)::x(2,2),means(2,2),vars(2,2,2),precwb(3,3),logdetwb
  integer::info,ord(3),info2
  logical::restrict(3,1),custom_restrict(3,1),het(4)
  type(fsv_sim_result)::sim,sim_alt
  type(fsv_options)::opt
  type(fsv_draws)::dr
  type(fsv_prediction)::pred
  call seed_rng(2468)
  allocate(bdense(5,2),bsparse(5,2),ipdef(5,3),fpdef(2,2));call make_dense_loadings(5,2,bdense);call make_sparse_loadings(5,2,.2_dp,bsparse);call default_fsv_parameters(5,2,ipdef,fpdef);if(bdense(1,1)<=0.0_dp.or.bdense(2,2)<=0.0_dp.or.any(count(abs(bsparse)>.2_dp,dim=1)<3).or.any(ipdef(:,3)<=0.0_dp))error stop 'default factor simulation helpers'
  b(:,1)=[1.0_dp,.7_dp,-.4_dp]
  ip(:,1)=[-1.5_dp,-1.2_dp,-1.0_dp];ip(:,2)=[.90_dp,.94_dp,.96_dp];ip(:,3)=[.25_dp,.20_dp,.15_dp]
  fp(1,:)=[.97_dp,.12_dp]
  het=[.false.,.true.,.true.,.false.]
  call simulate_fsv(20,b,ip,fp,sim_alt,heteroskedastic=het,df=7.0_dp)
  if(maxval(abs(sim_alt%h_idio(:,1)-ip(1,1)))>1.0e-12_dp.or.maxval(abs(sim_alt%h_factor(:,1)))>1.0e-12_dp.or.any(.not.ieee_is_finite(sim_alt%factors)))error stop 'homoskedastic and t-factor simulation'
  call simulate_fsv(100,b,ip,fp,sim)
  call fsv_covariance_path(b,sim%h_idio,sim%h_factor,cov);call fsv_correlation_path(b,sim%h_idio,sim%h_factor,cor)
  if(maxval(abs([(cor(1,1,info),info=1,100)]-1.0_dp))>1.0e-12_dp)error stop 'cor diagonal'
  call expweightcov(sim%y,.05_dp,50,ew);if(any(.not.ieee_is_finite(ew)))error stop 'ew covariance'
  if(ledermann(3)/=1)error stop 'Ledermann'
  call static_factor_initialize(sim%y,1,l0,f0,u0,info);if(info/=0.or.any(u0<=0.0_dp))error stop 'static initialize'
  call preorder(sim%y,1,ord,info);call findrestrict(sim%y,1,restrict,info);if(info/=0)error stop 'identification helpers'
  opt%draws=12;opt%burnin=8;opt%thin=1;opt%sv_sweeps=1
  call fit_fsv(sim%y,1,opt,dr)
  if(any(dr%para(3,:,:)<=0.0_dp))error stop 'factor SV sigma'
  if(dr%loadings(1,1,dr%ndraws)<0.0_dp)error stop 'loading sign'
  custom_restrict(:,1)=[.true.,.false.,.true.];opt%draws=4;opt%burnin=2;call fit_fsv(sim%y,1,opt,dr,custom_restrict);if(any(abs(dr%loadings(2,1,:))>1.0e-12_dp))error stop 'custom restriction'
  call running_covariance(dr,100,covm);call running_correlation(dr,100,corm)
  if(maxval(abs([(corm(info,info),info=1,3)]-1.0_dp))>1.0e-10_dp)error stop 'running cor'
  call predict_fsv(dr,2,2,pred);if(any(.not.ieee_is_finite(pred%cov)))error stop 'factor prediction'
  call woodbury_precision(dr%loadings(:,:,1),pred%h(1,1:3,1),pred%h(1,4:4,1),precwb,logdetwb,info2);if(info2/=0.or.maxval(abs(precwb-pred%precision(:,:,1,1)))>1.0e-8_dp)error stop 'Woodbury precision'
  call predict_conditional_fsv(dr,2,2,means_pred,vols_pred);if(any(vols_pred<=0.0_dp).or.any(.not.ieee_is_finite(means_pred)))error stop 'conditional prediction'
  call eigen_loading_diagnostics(dr,evals,info2);if(info2/=0.or.any(evals<0.0_dp))error stop 'loading eigen diagnostics'
  allocate(ll(2,size(pred%cov,4)));call predloglik_fsv(pred,sim%y(1:2,:),ll);if(any(.not.ieee_is_finite(ll)))error stop 'predictive loglik'
  allocate(llwb(2,size(pred%cov,4)),scores(2));call predloglik_fsv_woodbury(pred,sim%y(1:2,:),llwb);if(maxval(abs(ll-llwb))>1.0e-8_dp)error stop 'Woodbury predictive loglik';call aggregate_loglik_draws(ll,scores);if(any(.not.ieee_is_finite(scores)))error stop 'aggregate loglik'
  x=reshape([0.0_dp,0.0_dp,1.0_dp,-1.0_dp],[2,2]);means=0.0_dp;vars=0.0_dp;vars(1,1,:)=1.0_dp;vars(2,2,:)=1.0_dp
  call dmvnorm_columns(x,means,vars,.true.,dens);if(abs(dens(1)+log(2.0_dp*acos(-1.0_dp)))>1.0e-10_dp)error stop 'dmvnorm'
  call sign_identify(dr%loadings,dr%factors);call order_identify(dr%loadings,dr%factors)
  opt%draws=5;opt%burnin=3;opt%normal_gamma=.true.;call fit_fsv(sim%y,1,opt,dr);if(.not.allocated(dr%local_scale).or.any(dr%local_scale<=0.0_dp).or.any(dr%global_shrinkage<=0.0_dp))error stop 'Normal-Gamma shrinkage'
  call simulate_fsv(30,b0,ip,fp0,sim);if(size(sim%factors,2)/=0)error stop 'zero-factor simulation';opt%normal_gamma=.false.;opt%draws=3;opt%burnin=2;call fit_fsv(sim%y,0,opt,dr);if(dr%nfactors/=0.or.size(dr%loadings,2)/=0)error stop 'zero-factor fit'
  print '(a)','Factor-SV simulation, covariance, initialization, fitting, identification, prediction, and density tests passed.'
end program test_factor
