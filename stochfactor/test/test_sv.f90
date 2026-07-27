! SPDX-License-Identifier: GPL-2.0-or-later
! Copyright (C) 2013-2026 the original stochvol and factorstochvol authors
! and the modern Fortran translation contributors.
! This file is free software: you may redistribute it and/or modify it
! under the terms of the GNU General Public License as published by the
! Free Software Foundation, either version 2, or (at your option) any later version.
! This file is distributed without any warranty; see LICENSE for details.
program test_sv
  use sv_kinds,only:dp
  use sv_rng,only:seed_rng
  use sv_types
  use sv_core
  use sv_stats,only:mean1,log_returns
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  implicit none
  type(sv_params)::p
  type(sv_sim_result)::sim
  type(sv_prior)::pr
  type(sv_mcmc_options)::opt
  type(sv_draws)::dr
  type(sv_prediction)::pred
  real(dp)::prob(10),mm(10),vv(10),ll
  real(dp),allocatable::res(:),x(:,:),prices(:,:),rets(:,:),roll(:),prior_mean(:),prior_prec(:,:),breg(:)
  integer::n,i
  call seed_rng(12345);call get_omori_constants(prob,mm,vv)
  call assert_close(sum(prob),1.0_dp,1.0e-12_dp,'Omori probabilities')
  allocate(prices(3,1));prices(:,1)=[100.0_dp,101.0_dp,103.0_dp];call log_returns(prices,rets);call assert_close(rets(1,1),log(1.01_dp),1.0e-12_dp,'log returns')
  allocate(prior_mean(2),prior_prec(2,2),breg(2));prior_mean=0.0_dp;prior_prec=0.0_dp;prior_prec(1,1)=1.0_dp;prior_prec(2,2)=1.0_dp
  n=120;p%mu=-1.0_dp;p%phi=.92_dp;p%sigma=.25_dp
  call simulate_sv(n,p,sim)
  if(any(.not.(sim%vol>0.0_dp)))error stop 'simulation volatility'
  ll=sv_complete_loglik(sim%y,sim%latent,sim%latent0,p,sim%tau)
  if(.not.(ll<0.0_dp.and.ll>-huge(1.0_dp)))error stop 'complete likelihood'
  opt%draws=30;opt%burnin=20;opt%thin=1;opt%use_mixture=.true.;opt%store_latent=.true.
  call svsample(sim%y,p,pr,opt,dr)
  if(any(.not.ieee_is_finite(dr%para)))error stop 'SV draws finite'
  if(maxval(dr%accept(1:3))<=0.0_dp)error stop 'parameters did not move'
  allocate(res(n));call sv_residuals(sim%y,dr,res);if(abs(mean1(res))>1.5_dp)error stop 'residual mean'
  call predict_sv(dr,3,2,pred);if(any(pred%vola<=0.0_dp))error stop 'prediction volatility'
  ! Student-t plus regression and leverage general sampler.
  allocate(x(n,2));x(:,1)=1.0_dp;do i=1,n;x(i,2)=real(i,dp)/120.0_dp;end do
  call bayesian_regression_update(sim%y,x,prior_mean,prior_prec,breg);if(any(.not.ieee_is_finite(breg)))error stop 'Bayesian regression update'
  p%nu=8.0_dp;p%rho=-.35_dp;allocate(p%beta(2));p%beta=[.1_dp,-.2_dp]
  call simulate_sv(120,p,sim,x)
  opt%draws=20;opt%burnin=15;opt%use_mixture=.false.;opt%sample_rho=.true.;opt%sample_nu=.true.;opt%store_tau=.true.;opt%latent_sweeps=1
  call svsample(sim%y,p,pr,opt,dr,x)
  if(any(dr%para(4,:)<=2.0_dp).or.any(abs(dr%para(5,:))>=1.0_dp))error stop 't/leverage constraints'
  if(any(dr%tau<=0.0_dp))error stop 'tau positivity'
  opt%draws=2;opt%burnin=1;opt%latent_sweeps=1;call svsample_t_leverage(sim%y,p,pr,opt,dr,x);if(dr%ndraws/=2)error stop 't leverage wrapper'
  p%rho=0.0_dp;p%nu=huge(1.0_dp);if(allocated(p%beta))deallocate(p%beta);opt%sample_rho=.false.;opt%sample_nu=.false.;opt%use_mixture=.true.;opt%draws=2;opt%burnin=1;call rolling_sv_forecast(sim%y(1:41),40,p,pr,opt,roll);if(size(roll)/=1.or.any(roll<=0.0_dp))error stop 'rolling forecast'
  print '(a)','Univariate SV simulation, likelihood, MCMC, regression, t-error, leverage, residual, and prediction tests passed.'
contains
  subroutine assert_close(a,b,tol,msg)
    real(dp),intent(in)::a,b,tol;character(*),intent(in)::msg
    if(abs(a-b)>tol)then;print *,msg,a,b;error stop;end if
  end subroutine
end program test_sv
