! SPDX-License-Identifier: GPL-2.0-or-later
!
! Computational translation of MSGARCH, copyright (C) MSGARCH authors.
! Fortran translation copyright (C) 2026 translation contributors.
! Distributed under the GNU General Public License, version 2 or later.
module msgarch_estimation
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
   use msgarch_kinds, only : dp, pi
   use msgarch_types, only : msgarch_spec, fit_result, mcmc_result, filter_result
   use msgarch_models, only : parameter_count, default_parameters, parameter_bounds, unpack_parameters, spec_valid, &
      sort_parameters_by_variance
   use msgarch_filter, only : hamilton_filter
   use msgarch_optimizer, only : optimizer_result, nelder_mead
   use msgarch_linalg, only : numerical_hessian, invert_matrix, sample_mean_sd
   use msgarch_rng, only : random_normal, random_uniform
   implicit none
   private
   public :: fit_ml, fit_mcmc, log_posterior, dic_from_draws
contains
   function fit_ml(template,y,start,fixed_mask,fixed_values,tie_group,max_iterations,tolerance) result(fit)
      type(msgarch_spec),intent(in)::template
      real(dp),intent(in)::y(:)
      real(dp),intent(in),optional::start(:),fixed_values(:),tolerance
      logical,intent(in),optional::fixed_mask(:)
      integer,intent(in),optional::tie_group(:),max_iterations
      type(fit_result)::fit
      real(dp),allocatable::theta0(:),lower(:),upper(:),x0(:),xlower(:),xupper(:),full(:)
      real(dp),allocatable::hfree(:,:),covfree(:,:)
      logical,allocatable::fixed(:)
      integer,allocatable::tie(:),var_id(:)
      type(optimizer_result)::opt
      type(filter_result)::filtered
      logical::valid,success
      integer::nfull,nvar,i,j,g,maxit
      real(dp)::tol

      nfull=parameter_count(template);theta0=default_parameters(template)
      if(present(start))then
         if(size(start)/=nfull)error stop 'fit_ml: start has wrong size'
         theta0=start
      end if
      allocate(fixed(nfull),tie(nfull),var_id(nfull));fixed=.false.;tie=0
      if(present(fixed_mask))then
         if(size(fixed_mask)/=nfull)error stop 'fit_ml: fixed_mask has wrong size'
         fixed=fixed_mask
      end if
      if(present(fixed_values))then
         if(size(fixed_values)/=nfull)error stop 'fit_ml: fixed_values has wrong size'
         where(fixed)theta0=fixed_values
      end if
      if(present(tie_group))then
         if(size(tie_group)/=nfull)error stop 'fit_ml: tie_group has wrong size'
         tie=tie_group
      end if
      var_id=0;nvar=0
      do i=1,nfull
         if(fixed(i))cycle
         if(tie(i)>0)then
            g=0
            do j=1,i-1
               if(.not.fixed(j).and.tie(j)==tie(i))then;g=var_id(j);exit;end if
            end do
            if(g>0)then;var_id(i)=g;cycle;end if
         end if
         nvar=nvar+1;var_id(i)=nvar
      end do
      if(nvar==0)then
         call unpack_parameters(template,theta0,fit%spec,valid)
         if(valid)then;filtered=hamilton_filter(fit%spec,y);fit%loglik=filtered%loglik;fit%converged=.true.;end if
         fit%parameters=theta0;fit%aic=-2.0_dp*fit%loglik;fit%bic=-2.0_dp*fit%loglik
         allocate(fit%hessian(nfull,nfull),fit%covariance(nfull,nfull),fit%standard_error(nfull))
         fit%hessian=0.0_dp;fit%covariance=0.0_dp;fit%standard_error=0.0_dp
         return
      end if
      call parameter_bounds(template,lower,upper)
      allocate(x0(nvar),xlower(nvar),xupper(nvar));xlower=-huge(1.0_dp);xupper=huge(1.0_dp)
      do g=1,nvar
         do i=1,nfull
            if(var_id(i)==g)then;x0(g)=theta0(i);xlower(g)=max(xlower(g),lower(i));xupper(g)=min(xupper(g),upper(i));exit;end if
         end do
         do i=1,nfull
            if(var_id(i)==g)then;xlower(g)=max(xlower(g),lower(i));xupper(g)=min(xupper(g),upper(i));end if
         end do
      end do
      maxit=1000;if(present(max_iterations))maxit=max_iterations
      tol=1.0e-7_dp;if(present(tolerance))tol=tolerance
      opt=nelder_mead(objective,x0,tolerance=tol,max_iterations=maxit,lower=xlower,upper=xupper)
      full=expand(opt%x);call unpack_parameters(template,full,fit%spec,valid)
      fit%parameters=full;fit%iterations=opt%iterations;fit%evaluations=opt%evaluations;fit%converged=opt%converged.and.valid
      if(valid)then
         filtered=hamilton_filter(fit%spec,y);fit%loglik=filtered%loglik
      end if
      fit%aic=-2.0_dp*fit%loglik+2.0_dp*real(nvar,dp)
      fit%bic=-2.0_dp*fit%loglik+log(real(size(y)-1,dp))*real(nvar,dp)
      hfree=numerical_hessian(objective,opt%x,1.0e-4_dp);call invert_matrix(hfree,covfree,success)
      allocate(fit%hessian(nfull,nfull),fit%covariance(nfull,nfull),fit%standard_error(nfull))
      fit%hessian=0.0_dp;fit%covariance=0.0_dp
      do i=1,nfull
         if(var_id(i)>0)then
            do j=1,nfull
               if(var_id(j)>0)then
                  fit%hessian(i,j)=hfree(var_id(i),var_id(j))
                  if(success)fit%covariance(i,j)=covfree(var_id(i),var_id(j))
               end if
            end do
         end if
      end do
      do i=1,nfull;fit%standard_error(i)=sqrt(max(fit%covariance(i,i),0.0_dp));end do
   contains
      function expand(x) result(theta)
         real(dp),intent(in)::x(:)
         real(dp),allocatable::theta(:)
         integer::ii
         theta=theta0
         do ii=1,nfull;if(var_id(ii)>0)theta(ii)=x(var_id(ii));end do
      end function expand
      function objective(x) result(value)
         real(dp),intent(in)::x(:)
         real(dp)::value
         real(dp),allocatable::theta(:)
         type(msgarch_spec)::trial
         type(filter_result)::filt
         logical::ok
         theta=expand(x);call unpack_parameters(template,theta,trial,ok)
         if(.not.ok)then;value=1.0e100_dp;return;end if
         filt=hamilton_filter(trial,y)
         if(.not.ieee_is_finite(filt%loglik))then;value=1.0e100_dp;else;value=-filt%loglik;end if
      end function objective
   end function fit_ml

   function log_posterior(template,y,theta,prior_mean,prior_sd) result(value)
      type(msgarch_spec),intent(in)::template
      real(dp),intent(in)::y(:),theta(:)
      real(dp),intent(in),optional::prior_mean(:),prior_sd(:)
      real(dp)::value
      real(dp),allocatable::mean(:),sd(:)
      type(msgarch_spec)::spec
      type(filter_result)::filtered
      logical::valid
      integer::n
      n=size(theta);allocate(mean(n),sd(n));mean=default_parameters(template);sd=1000.0_dp
      if(present(prior_mean))mean=prior_mean
      if(present(prior_sd))sd=prior_sd
      call unpack_parameters(template,theta,spec,valid)
      if(.not.valid.or.any(sd<=0.0_dp))then;value=-huge(1.0_dp);return;end if
      filtered=hamilton_filter(spec,y)
      if(.not.ieee_is_finite(filtered%loglik))then;value=-huge(1.0_dp);return;end if
      value=filtered%loglik-0.5_dp*sum(((theta-mean)/sd)**2)-sum(log(sd))-0.5_dp*real(n,dp)*log(2.0_dp*pi)
   end function log_posterior

   function fit_mcmc(template,y,n_iterations,burn,thin,start,proposal_sd,prior_mean,prior_sd,adapt,sort_regimes) result(result)
      type(msgarch_spec),intent(in)::template
      real(dp),intent(in)::y(:)
      integer,intent(in)::n_iterations,burn,thin
      real(dp),intent(in),optional::start(:),proposal_sd(:),prior_mean(:),prior_sd(:)
      logical,intent(in),optional::adapt,sort_regimes
      type(mcmc_result)::result
      real(dp),allocatable::current(:),proposal(:),scale(:),kept_lp(:),mean(:),sd(:)
      real(dp)::lp_current,lp_proposal,rate
      integer::n,i,j,nkeep,accepted,batch_accepted,kept
      logical::do_adapt,do_sort
      if(n_iterations<=burn.or.thin<1)error stop 'fit_mcmc: invalid iteration controls'
      n=parameter_count(template);current=default_parameters(template)
      if(present(start))current=start
      allocate(scale(n),proposal(n));scale=0.03_dp*max(1.0_dp,abs(current));if(present(proposal_sd))scale=proposal_sd
      allocate(mean(n),sd(n));mean=default_parameters(template);sd=1000.0_dp
      if(present(prior_mean))mean=prior_mean
      if(present(prior_sd))sd=prior_sd
      do_adapt=.true.;if(present(adapt))do_adapt=adapt
      do_sort=.true.;if(present(sort_regimes))do_sort=sort_regimes
      nkeep=(n_iterations-burn)/thin;allocate(result%draws(nkeep,n),result%log_posterior(nkeep),kept_lp(nkeep))
      lp_current=log_posterior(template,y,current,mean,sd);accepted=0;batch_accepted=0;kept=0
      do i=1,n_iterations
         do j=1,n;proposal(j)=current(j)+scale(j)*random_normal();end do
         lp_proposal=log_posterior(template,y,proposal,mean,sd)
         if(log(random_uniform())<lp_proposal-lp_current)then
            current=proposal;lp_current=lp_proposal;accepted=accepted+1;batch_accepted=batch_accepted+1
         end if
         if(do_adapt.and.i<=burn.and.mod(i,50)==0)then
            rate=real(batch_accepted,dp)/50.0_dp
            if(rate<0.15_dp)scale=0.8_dp*scale
            if(rate>0.40_dp)scale=1.2_dp*scale
            batch_accepted=0
         end if
         if(i>burn.and.mod(i-burn,thin)==0)then
            kept=kept+1
            if(do_sort)then;result%draws(kept,:)=sort_parameters_by_variance(template,current)
            else;result%draws(kept,:)=current;end if
            result%log_posterior(kept)=lp_current
         end if
      end do
      result%acceptance_rate=real(accepted,dp)/real(n_iterations,dp)
      call sample_mean_sd(result%draws,result%posterior_mean,result%posterior_sd)
      result%dic=dic_from_draws(template,y,result%draws)
   end function fit_mcmc

   function dic_from_draws(template,y,draws) result(dic)
      type(msgarch_spec),intent(in)::template
      real(dp),intent(in)::y(:),draws(:,:)
      real(dp)::dic,mean_deviance,deviance_mean
      real(dp),allocatable::theta_mean(:)
      type(msgarch_spec)::spec
      type(filter_result)::filtered
      logical::valid
      integer::i
      mean_deviance=0.0_dp
      do i=1,size(draws,1)
         call unpack_parameters(template,draws(i,:),spec,valid)
         if(valid)then;filtered=hamilton_filter(spec,y);mean_deviance=mean_deviance-2.0_dp*filtered%loglik;end if
      end do
      mean_deviance=mean_deviance/real(size(draws,1),dp)
      theta_mean=sum(draws,dim=1)/real(size(draws,1),dp);call unpack_parameters(template,theta_mean,spec,valid)
      if(valid)then;filtered=hamilton_filter(spec,y);deviance_mean=-2.0_dp*filtered%loglik;else;deviance_mean=mean_deviance;end if
      dic=2.0_dp*mean_deviance-deviance_mean
   end function dic_from_draws
end module msgarch_estimation
