! Part of the experimental modern Fortran translation of rugarch 1.5-6.
! This file was created or modified for the Fortran project on 2026-07-23.
! Original rugarch authors retain copyright; see NOTICE and ORIGIN.md.
! SPDX-License-Identifier: GPL-3.0-only

module rugarch_fit
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
   use rugarch_kinds, only : dp
   use rugarch_distributions, only : dist_norm,dist_snorm,dist_std,dist_sstd,dist_ged,dist_sged,dist_jsu, &
      dist_nig,dist_ghyp,dist_ghst,distribution_pdf
   use rugarch_optimizer, only : optimizer_result,nelder_mead
   use rugarch_types, only : garch_spec,garch_fit_result,distribution_fit_result,make_garch_spec, &
      model_sgarch,model_gjrgarch,model_egarch,model_aparch,model_igarch,model_figarch,model_csgarch, &
      model_realgarch,model_fgarch,configure_fgarch_submodel,fgarch_garch,fgarch_tgarch,fgarch_avgarch, &
      fgarch_ngarch,fgarch_nagarch,fgarch_aparch,fgarch_allgarch,fgarch_gjrgarch
   use rugarch_models, only : garch_log_likelihood,realgarch_log_likelihood,true_persistence,figarch_weights
   implicit none
   private

   type :: dist_fit_context
      real(dp),allocatable::y(:)
      integer::kind=dist_norm
      logical::fit_shape=.false.,fit_skew=.false.,fit_lambda=.false.
   end type dist_fit_context

   type :: garch_fit_context
      real(dp),allocatable::y(:),realized(:)
      integer::model=model_sgarch,cond_dist=dist_norm,p=1,q=1,figarch_truncation=500
      integer::fgarch_submodel=fgarch_allgarch
      logical::fit_mean=.true.,fit_shape=.false.,fit_skew=.false.,fit_lambda=.false.
   end type garch_fit_context

   public :: fit_distribution,fit_model,fit_model11
   public :: fit_garch11,fit_gjrgarch11,fit_egarch11,fit_aparch11,fit_igarch11
   public :: fit_figarch,fit_figarch11,fit_csgarch,fit_csgarch11
   public :: fit_realgarch,fit_realgarch11,fit_fgarch,fit_fgarch11

contains

   function fit_distribution(y,kind,fit_shape,fit_skew,fit_lambda,max_iterations) result(fit)
      real(dp),intent(in)::y(:)
      integer,intent(in)::kind
      logical,intent(in),optional::fit_shape,fit_skew,fit_lambda
      integer,intent(in),optional::max_iterations
      type(distribution_fit_result)::fit
      type(dist_fit_context)::context
      type(optimizer_result)::opt
      real(dp),allocatable::x0(:)
      real(dp)::mean0,sd0
      logical::use_shape,use_skew,use_lambda
      integer::npar,idx,maxit

      use_shape=has_shape(kind);if(present(fit_shape))use_shape=fit_shape.and.use_shape
      use_skew=has_skew(kind);if(present(fit_skew))use_skew=fit_skew.and.use_skew
      use_lambda=has_lambda(kind);if(present(fit_lambda))use_lambda=fit_lambda.and.use_lambda
      maxit=1800;if(present(max_iterations))maxit=max_iterations
      context%y=y;context%kind=kind;context%fit_shape=use_shape;context%fit_skew=use_skew;context%fit_lambda=use_lambda
      mean0=sum(y)/real(size(y),dp)
      sd0=sqrt(max(sum((y-mean0)**2)/real(max(1,size(y)-1),dp),1.0e-12_dp))
      npar=2+merge(1,0,use_shape)+merge(1,0,use_skew)+merge(1,0,use_lambda)
      allocate(x0(npar));x0(1)=mean0;x0(2)=log(sd0);idx=3
      if(use_shape)then;x0(idx)=initial_shape_parameter(kind);idx=idx+1;end if
      if(use_skew)then;x0(idx)=initial_skew_parameter(kind);idx=idx+1;end if
      if(use_lambda)x0(idx)=initial_lambda_parameter()
      opt=nelder_mead(distribution_objective,x0,context,step=0.12_dp,tolerance=1.0e-8_dp,max_iterations=maxit)
      call decode_distribution_parameters(opt%x,context,fit%mean,fit%sd,fit%shape,fit%skew,fit%lambda)
      fit%log_likelihood=-opt%objective;fit%iterations=opt%iterations;fit%status=opt%status
   end function fit_distribution

   function fit_garch11(y,cond_dist,fit_mean,fit_shape,fit_skew,max_iterations) result(fit)
      real(dp),intent(in)::y(:);integer,intent(in),optional::cond_dist,max_iterations
      logical,intent(in),optional::fit_mean,fit_shape,fit_skew;type(garch_fit_result)::fit
      fit=fit_model(y,model_sgarch,1,1,cond_dist,fit_mean,fit_shape,fit_skew,max_iterations=max_iterations)
   end function fit_garch11

   function fit_gjrgarch11(y,cond_dist,fit_mean,fit_shape,fit_skew,max_iterations) result(fit)
      real(dp),intent(in)::y(:);integer,intent(in),optional::cond_dist,max_iterations
      logical,intent(in),optional::fit_mean,fit_shape,fit_skew;type(garch_fit_result)::fit
      fit=fit_model(y,model_gjrgarch,1,1,cond_dist,fit_mean,fit_shape,fit_skew,max_iterations=max_iterations)
   end function fit_gjrgarch11

   function fit_egarch11(y,cond_dist,fit_mean,fit_shape,fit_skew,max_iterations) result(fit)
      real(dp),intent(in)::y(:);integer,intent(in),optional::cond_dist,max_iterations
      logical,intent(in),optional::fit_mean,fit_shape,fit_skew;type(garch_fit_result)::fit
      fit=fit_model(y,model_egarch,1,1,cond_dist,fit_mean,fit_shape,fit_skew,max_iterations=max_iterations)
   end function fit_egarch11

   function fit_aparch11(y,cond_dist,fit_mean,fit_shape,fit_skew,max_iterations) result(fit)
      real(dp),intent(in)::y(:);integer,intent(in),optional::cond_dist,max_iterations
      logical,intent(in),optional::fit_mean,fit_shape,fit_skew;type(garch_fit_result)::fit
      fit=fit_model(y,model_aparch,1,1,cond_dist,fit_mean,fit_shape,fit_skew,max_iterations=max_iterations)
   end function fit_aparch11

   function fit_igarch11(y,cond_dist,fit_mean,fit_shape,fit_skew,max_iterations) result(fit)
      real(dp),intent(in)::y(:);integer,intent(in),optional::cond_dist,max_iterations
      logical,intent(in),optional::fit_mean,fit_shape,fit_skew;type(garch_fit_result)::fit
      fit=fit_model(y,model_igarch,1,1,cond_dist,fit_mean,fit_shape,fit_skew,max_iterations=max_iterations)
   end function fit_igarch11

   function fit_model11(y,model,cond_dist,fit_mean,fit_shape,fit_skew,max_iterations) result(fit)
      real(dp),intent(in)::y(:);integer,intent(in)::model
      integer,intent(in),optional::cond_dist,max_iterations
      logical,intent(in),optional::fit_mean,fit_shape,fit_skew
      type(garch_fit_result)::fit
      fit=fit_model(y,model,1,1,cond_dist,fit_mean,fit_shape,fit_skew,max_iterations=max_iterations)
   end function fit_model11

   function fit_figarch(y,p,q,cond_dist,fit_mean,fit_shape,fit_skew,fit_lambda,truncation,max_iterations) result(fit)
      real(dp),intent(in)::y(:);integer,intent(in),optional::p,q,cond_dist,truncation,max_iterations
      logical,intent(in),optional::fit_mean,fit_shape,fit_skew,fit_lambda
      type(garch_fit_result)::fit;integer::pp,qq
      pp=1;if(present(p))pp=p;qq=1;if(present(q))qq=q
      fit=fit_model(y,model_figarch,pp,qq,cond_dist,fit_mean,fit_shape,fit_skew,fit_lambda=fit_lambda, &
         figarch_truncation=truncation,max_iterations=max_iterations)
   end function fit_figarch

   function fit_figarch11(y,cond_dist,fit_mean,fit_shape,fit_skew,fit_lambda,truncation,max_iterations) result(fit)
      real(dp),intent(in)::y(:);integer,intent(in),optional::cond_dist,truncation,max_iterations
      logical,intent(in),optional::fit_mean,fit_shape,fit_skew,fit_lambda;type(garch_fit_result)::fit
      fit=fit_figarch(y,1,1,cond_dist,fit_mean,fit_shape,fit_skew,fit_lambda,truncation,max_iterations)
   end function fit_figarch11

   function fit_csgarch(y,p,q,cond_dist,fit_mean,fit_shape,fit_skew,fit_lambda,max_iterations) result(fit)
      real(dp),intent(in)::y(:);integer,intent(in),optional::p,q,cond_dist,max_iterations
      logical,intent(in),optional::fit_mean,fit_shape,fit_skew,fit_lambda
      type(garch_fit_result)::fit;integer::pp,qq
      pp=1;if(present(p))pp=p;qq=1;if(present(q))qq=q
      fit=fit_model(y,model_csgarch,pp,qq,cond_dist,fit_mean,fit_shape,fit_skew,fit_lambda=fit_lambda, &
         max_iterations=max_iterations)
   end function fit_csgarch

   function fit_csgarch11(y,cond_dist,fit_mean,fit_shape,fit_skew,fit_lambda,max_iterations) result(fit)
      real(dp),intent(in)::y(:);integer,intent(in),optional::cond_dist,max_iterations
      logical,intent(in),optional::fit_mean,fit_shape,fit_skew,fit_lambda;type(garch_fit_result)::fit
      fit=fit_csgarch(y,1,1,cond_dist,fit_mean,fit_shape,fit_skew,fit_lambda,max_iterations)
   end function fit_csgarch11

   function fit_fgarch(y,p,q,cond_dist,fit_mean,fit_shape,fit_skew,fit_lambda,max_iterations,submodel) result(fit)
      real(dp),intent(in)::y(:);integer,intent(in),optional::p,q,cond_dist,max_iterations,submodel
      logical,intent(in),optional::fit_mean,fit_shape,fit_skew,fit_lambda
      type(garch_fit_result)::fit;integer::pp,qq,sm
      pp=1;if(present(p))pp=p;qq=1;if(present(q))qq=q
      sm=fgarch_allgarch;if(present(submodel))sm=submodel
      fit=fit_model(y,model_fgarch,pp,qq,cond_dist,fit_mean,fit_shape,fit_skew,fit_lambda=fit_lambda, &
         max_iterations=max_iterations,fgarch_submodel=sm)
   end function fit_fgarch

   function fit_fgarch11(y,cond_dist,fit_mean,fit_shape,fit_skew,fit_lambda,max_iterations,submodel) result(fit)
      real(dp),intent(in)::y(:);integer,intent(in),optional::cond_dist,max_iterations,submodel
      logical,intent(in),optional::fit_mean,fit_shape,fit_skew,fit_lambda;type(garch_fit_result)::fit
      fit=fit_fgarch(y,1,1,cond_dist,fit_mean,fit_shape,fit_skew,fit_lambda,max_iterations,submodel)
   end function fit_fgarch11

   function fit_realgarch(returns,realized,p,q,cond_dist,fit_mean,fit_shape,fit_skew,fit_lambda,max_iterations) result(fit)
      real(dp),intent(in)::returns(:),realized(:);integer,intent(in),optional::p,q,cond_dist,max_iterations
      logical,intent(in),optional::fit_mean,fit_shape,fit_skew,fit_lambda
      type(garch_fit_result)::fit;integer::pp,qq
      pp=1;if(present(p))pp=p;qq=1;if(present(q))qq=q
      fit=fit_model(returns,model_realgarch,pp,qq,cond_dist,fit_mean,fit_shape,fit_skew,fit_lambda, &
         realized=realized,max_iterations=max_iterations)
   end function fit_realgarch

   function fit_realgarch11(returns,realized,cond_dist,fit_mean,fit_shape,fit_skew,fit_lambda,max_iterations) result(fit)
      real(dp),intent(in)::returns(:),realized(:);integer,intent(in),optional::cond_dist,max_iterations
      logical,intent(in),optional::fit_mean,fit_shape,fit_skew,fit_lambda;type(garch_fit_result)::fit
      fit=fit_realgarch(returns,realized,1,1,cond_dist,fit_mean,fit_shape,fit_skew,fit_lambda,max_iterations)
   end function fit_realgarch11

   function fit_model(y,model,p,q,cond_dist,fit_mean,fit_shape,fit_skew,fit_lambda,realized, &
      figarch_truncation,max_iterations,fgarch_submodel) result(fit)
      real(dp),intent(in)::y(:);integer,intent(in)::model,p,q
      integer,intent(in),optional::cond_dist,figarch_truncation,max_iterations,fgarch_submodel
      logical,intent(in),optional::fit_mean,fit_shape,fit_skew,fit_lambda
      real(dp),intent(in),optional::realized(:)
      type(garch_fit_result)::fit
      type(garch_fit_context)::context
      type(optimizer_result)::opt
      real(dp),allocatable::x0(:)
      real(dp)::mean0,var0
      integer::kind,npar,maxit,k
      logical::use_mean,use_shape,use_skew,use_lambda

      kind=dist_norm;if(present(cond_dist))kind=cond_dist
      use_mean=.true.;if(present(fit_mean))use_mean=fit_mean
      use_shape=has_shape(kind);if(present(fit_shape))use_shape=fit_shape.and.use_shape
      use_skew=has_skew(kind);if(present(fit_skew))use_skew=fit_skew.and.use_skew
      use_lambda=has_lambda(kind);if(present(fit_lambda))use_lambda=fit_lambda.and.use_lambda
      maxit=3500;if(present(max_iterations))maxit=max_iterations

      if(p<0 .or. q<0 .or. p+q==0 .or. size(y)<max(30,5*(p+q+1)))then
         fit%status=2;fit%message='invalid model order or insufficient observations';return
      end if
      if(model==model_realgarch)then
         if(.not.present(realized) .or. size(realized)/=size(y) .or. any(realized<=0.0_dp))then
            fit%status=2;fit%message='realGARCH requires a positive realized series of matching length';return
         end if
      end if

      context%y=y;context%model=model;context%cond_dist=kind;context%p=p;context%q=q
      context%fit_mean=use_mean;context%fit_shape=use_shape;context%fit_skew=use_skew;context%fit_lambda=use_lambda
      if(present(realized))context%realized=realized
      if(present(figarch_truncation))context%figarch_truncation=max(20,figarch_truncation)
      if(present(fgarch_submodel))context%fgarch_submodel=fgarch_submodel
      mean0=sum(y)/real(size(y),dp);var0=max(sum((y-mean0)**2)/real(max(1,size(y)-1),dp),1.0e-10_dp)
      npar=model_parameter_count(context)+merge(1,0,use_mean)+merge(1,0,use_shape)+merge(1,0,use_skew)+merge(1,0,use_lambda)
      allocate(x0(npar));call initial_garch_parameters(x0,context,mean0,var0)
      opt=nelder_mead(garch_objective,x0,context,step=0.10_dp,tolerance=3.0e-7_dp,max_iterations=maxit)
      call decode_garch_parameters(opt%x,context,fit%spec)
      k=size(x0)
      if(model==model_realgarch)then
         allocate(fit%residuals(size(y)),fit%sigma(size(y)),fit%measurement_residuals(size(y)))
         fit%log_likelihood=realgarch_log_likelihood(y,context%realized,fit%spec,fit%residuals,fit%sigma, &
            fit%measurement_residuals)
      else
         allocate(fit%residuals(size(y)),fit%sigma(size(y)))
         fit%log_likelihood=garch_log_likelihood(y,fit%spec,fit%residuals,fit%sigma)
      end if
      fit%aic=-2.0_dp*fit%log_likelihood+2.0_dp*real(k,dp)
      fit%bic=-2.0_dp*fit%log_likelihood+log(real(size(y),dp))*real(k,dp)
      fit%iterations=opt%iterations;fit%evaluations=opt%evaluations;fit%status=opt%status
      if(opt%status==0)then;fit%message='converged';else;fit%message='iteration limit reached; inspect estimates';end if
   end function fit_model

   function distribution_objective(x,generic_context) result(value)
      real(dp),intent(in)::x(:);class(*),intent(in)::generic_context
      real(dp)::value,mean,sd,shape,skew,lambda,density
      integer::i
      select type(context=>generic_context)
      type is(dist_fit_context)
         call decode_distribution_parameters(x,context,mean,sd,shape,skew,lambda)
         value=0.0_dp
         do i=1,size(context%y)
            density=distribution_pdf((context%y(i)-mean)/sd,context%kind,shape,skew,lambda)/sd
            if(density<=tiny(1.0_dp) .or. .not.finite_value(density))then;value=huge(1.0_dp)/100.0_dp;return;end if
            value=value-log(density)
         end do
      class default
         value=huge(1.0_dp)/100.0_dp
      end select
   end function distribution_objective

   function garch_objective(x,generic_context) result(value)
      real(dp),intent(in)::x(:);class(*),intent(in)::generic_context
      real(dp)::value,llh,persist
      type(garch_spec)::spec
      real(dp),allocatable::psi(:)
      select type(context=>generic_context)
      type is(garch_fit_context)
         call decode_garch_parameters(x,context,spec)
         if(context%model==model_figarch)then
            allocate(psi(spec%figarch_truncation));call figarch_weights(spec%frac_d,spec%alpha,spec%beta,psi)
            if(any(psi< -1.0e-10_dp))then;value=huge(1.0_dp)/100.0_dp;return;end if
         else if(context%model==model_csgarch)then
            if(spec%rho<=sum(spec%alpha)+sum(spec%beta) .or. spec%phi>=sum(spec%beta))then
               value=huge(1.0_dp)/100.0_dp;return
            end if
         else if(context%model==model_realgarch)then
            persist=sum(spec%beta)+spec%phi*sum(spec%alpha)
            if(persist>=0.999_dp)then;value=huge(1.0_dp)/100.0_dp;return;end if
         end if
         if(context%model==model_realgarch)then
            llh=realgarch_log_likelihood(context%y,context%realized,spec)
         else
            llh=garch_log_likelihood(context%y,spec)
         end if
         if(.not.finite_value(llh) .or. llh<=-huge(1.0_dp)/10.0_dp)then
            value=huge(1.0_dp)/100.0_dp
         else
            value=-llh
         end if
      class default
         value=huge(1.0_dp)/100.0_dp
      end select
   end function garch_objective

   subroutine decode_distribution_parameters(x,context,mean,sd,shape,skew,lambda)
      real(dp),intent(in)::x(:);type(dist_fit_context),intent(in)::context
      real(dp),intent(out)::mean,sd,shape,skew,lambda
      integer::idx
      mean=x(1);sd=exp(max(-30.0_dp,min(30.0_dp,x(2))));shape=default_shape(context%kind)
      skew=default_skew(context%kind);lambda=1.0_dp;idx=3
      if(context%fit_shape)then;shape=decode_shape(x(idx),context%kind);idx=idx+1;end if
      if(context%fit_skew)then;skew=decode_skew(x(idx),context%kind);idx=idx+1;end if
      if(context%fit_lambda)lambda=4.0_dp*tanh(x(idx))
   end subroutine decode_distribution_parameters

   subroutine initial_garch_parameters(x,context,mean0,var0)
      real(dp),intent(out)::x(:);type(garch_fit_context),intent(in)::context
      real(dp),intent(in)::mean0,var0
      integer::idx,j
      idx=1;if(context%fit_mean)then;x(idx)=mean0;idx=idx+1;end if
      if(context%model==model_egarch .or. context%model==model_realgarch)then
         x(idx)=log(var0)*0.05_dp
      else
         x(idx)=log(max(0.05_dp*var0,1.0e-12_dp))
      end if
      idx=idx+1
      select case(context%model)
      case(model_egarch)
         do j=1,context%p;x(idx)=atanh(0.10_dp);idx=idx+1;end do
         do j=1,context%q;x(idx)=log(0.85_dp/0.10_dp);idx=idx+1;end do
         do j=1,context%p;x(idx)=0.0_dp;idx=idx+1;end do
      case(model_gjrgarch)
         do j=1,context%p;x(idx)=log(0.06_dp/0.05_dp);idx=idx+1;end do
         do j=1,context%q;x(idx)=log(0.82_dp/0.05_dp);idx=idx+1;end do
         do j=1,context%p;x(idx)=log(0.05_dp/0.05_dp);idx=idx+1;end do
      case(model_aparch)
         do j=1,context%p;x(idx)=log(0.08_dp/0.07_dp);idx=idx+1;end do
         do j=1,context%q;x(idx)=log(0.82_dp/0.07_dp);idx=idx+1;end do
         do j=1,context%p;x(idx)=0.0_dp;idx=idx+1;end do
         x(idx)=0.0_dp;idx=idx+1
      case(model_igarch)
         do j=1,context%p;x(idx)=log(0.08_dp);idx=idx+1;end do
         do j=1,context%q;x(idx)=log(0.92_dp);idx=idx+1;end do
      case(model_figarch)
         do j=1,context%p;x(idx)=logit(0.15_dp/real(max(1,context%p),dp));idx=idx+1;end do
         do j=1,context%q;x(idx)=logit(0.40_dp/real(max(1,context%q),dp));idx=idx+1;end do
         x(idx)=logit(0.40_dp);idx=idx+1
      case(model_csgarch)
         do j=1,context%p;x(idx)=log(0.04_dp/0.10_dp);idx=idx+1;end do
         do j=1,context%q;x(idx)=log(0.70_dp/0.10_dp);idx=idx+1;end do
         x(idx)=logit((0.95_dp-0.05_dp)/0.949_dp);idx=idx+1;x(idx)=logit(0.05_dp/0.70_dp);idx=idx+1
      case(model_realgarch)
         do j=1,context%p;x(idx)=log(0.15_dp/0.10_dp);idx=idx+1;end do
         do j=1,context%q;x(idx)=log(0.65_dp/0.10_dp);idx=idx+1;end do
         x(idx)=0.0_dp;idx=idx+1;x(idx)=logit((0.90_dp-0.05_dp)/1.95_dp);idx=idx+1
         x(idx)=0.0_dp;idx=idx+1;x(idx)=0.0_dp;idx=idx+1;x(idx)=log(0.15_dp);idx=idx+1
      case(model_fgarch)
         do j=1,context%p;x(idx)=log(0.08_dp/0.07_dp);idx=idx+1;end do
         do j=1,context%q;x(idx)=log(0.82_dp/0.07_dp);idx=idx+1;end do
         select case(context%fgarch_submodel)
         case(fgarch_tgarch,fgarch_gjrgarch)
            do j=1,context%p;x(idx)=0.0_dp;idx=idx+1;end do
         case(fgarch_avgarch)
            do j=1,context%p;x(idx)=0.0_dp;idx=idx+1;end do
            do j=1,context%p;x(idx)=0.0_dp;idx=idx+1;end do
         case(fgarch_ngarch)
            x(idx)=logit((2.0_dp-0.01_dp)/3.99_dp);idx=idx+1
         case(fgarch_nagarch)
            do j=1,context%p;x(idx)=0.0_dp;idx=idx+1;end do
         case(fgarch_aparch)
            do j=1,context%p;x(idx)=0.0_dp;idx=idx+1;end do
            x(idx)=logit((1.0_dp-0.01_dp)/3.99_dp);idx=idx+1
         case(fgarch_allgarch)
            do j=1,context%p;x(idx)=0.0_dp;idx=idx+1;end do
            do j=1,context%p;x(idx)=0.0_dp;idx=idx+1;end do
            x(idx)=logit((2.0_dp-0.01_dp)/3.99_dp);idx=idx+1
         end select
      case default
         do j=1,context%p;x(idx)=log(0.08_dp/0.07_dp);idx=idx+1;end do
         do j=1,context%q;x(idx)=log(0.82_dp/0.07_dp);idx=idx+1;end do
      end select
      if(context%fit_shape)then;x(idx)=initial_shape_parameter(context%cond_dist);idx=idx+1;end if
      if(context%fit_skew)then;x(idx)=initial_skew_parameter(context%cond_dist);idx=idx+1;end if
      if(context%fit_lambda)x(idx)=initial_lambda_parameter()
   end subroutine initial_garch_parameters

   subroutine decode_garch_parameters(x,context,spec)
      real(dp),intent(in)::x(:);type(garch_fit_context),intent(in)::context;type(garch_spec),intent(out)::spec
      integer::idx,j,ncoef
      real(dp),allocatable::raw(:),coef(:)
      real(dp)::rho,beta_sum
      spec=make_garch_spec(context%p,context%q,context%model,context%cond_dist)
      spec%figarch_truncation=context%figarch_truncation
      if(context%model==model_fgarch)call configure_fgarch_submodel(spec,context%fgarch_submodel)
      idx=1;if(context%fit_mean)then;spec%mean=x(idx);idx=idx+1;end if
      if(context%model==model_egarch .or. context%model==model_realgarch)then
         spec%omega=max(-20.0_dp,min(20.0_dp,x(idx)))
      else
         spec%omega=exp(max(-35.0_dp,min(20.0_dp,x(idx))))
      end if
      idx=idx+1
      select case(context%model)
      case(model_egarch)
         do j=1,context%p;spec%alpha(j)=0.75_dp*tanh(x(idx));idx=idx+1;end do
         allocate(raw(context%q),coef(context%q));raw=x(idx:idx+context%q-1);call positive_coefficients(raw,0.995_dp,coef)
         spec%beta=coef;idx=idx+context%q
         do j=1,context%p;spec%gamma(j)=tanh(x(idx));idx=idx+1;end do
      case(model_gjrgarch)
         ncoef=2*context%p+context%q;allocate(raw(ncoef),coef(ncoef));raw=x(idx:idx+ncoef-1)
         call positive_coefficients(raw,0.985_dp,coef)
         spec%alpha=coef(1:context%p);spec%beta=coef(context%p+1:context%p+context%q)
         spec%gamma=2.0_dp*coef(context%p+context%q+1:ncoef);idx=idx+ncoef
      case(model_aparch)
         ncoef=context%p+context%q;allocate(raw(ncoef),coef(ncoef));raw=x(idx:idx+ncoef-1)
         call positive_coefficients(raw,0.97_dp,coef);spec%alpha=coef(1:context%p);spec%beta=coef(context%p+1:ncoef);idx=idx+ncoef
         do j=1,context%p;spec%gamma(j)=0.999_dp*tanh(x(idx));idx=idx+1;end do
         spec%delta=0.2_dp+3.8_dp*logistic(x(idx));idx=idx+1
      case(model_igarch)
         ncoef=context%p+context%q;allocate(raw(ncoef),coef(ncoef));raw=x(idx:idx+ncoef-1)
         call unit_coefficients(raw,coef);spec%alpha=coef(1:context%p);spec%beta=coef(context%p+1:ncoef);idx=idx+ncoef
      case(model_figarch)
         do j=1,context%p;spec%alpha(j)=0.999_dp*logistic(x(idx));idx=idx+1;end do
         do j=1,context%q;spec%beta(j)=0.999_dp*logistic(x(idx));idx=idx+1;end do
         spec%frac_d=0.001_dp+0.998_dp*logistic(x(idx));idx=idx+1
      case(model_csgarch)
         ncoef=context%p+context%q;allocate(raw(ncoef),coef(ncoef));raw=x(idx:idx+ncoef-1);idx=idx+ncoef
         rho=0.05_dp+0.949_dp*logistic(x(idx));idx=idx+1
         call positive_coefficients(raw,0.95_dp*rho,coef)
         spec%alpha=coef(1:context%p);spec%beta=coef(context%p+1:ncoef);spec%rho=rho
         beta_sum=sum(spec%beta);spec%phi=max(1.0e-8_dp,beta_sum)*0.99_dp*logistic(x(idx));idx=idx+1
      case(model_realgarch)
         ncoef=context%p+context%q;allocate(raw(ncoef),coef(ncoef));raw=x(idx:idx+ncoef-1);idx=idx+ncoef
         call positive_coefficients(raw,0.92_dp,coef);spec%alpha=coef(1:context%p);spec%beta=coef(context%p+1:ncoef)
         spec%xi=x(idx);idx=idx+1;spec%phi=0.05_dp+1.95_dp*logistic(x(idx));idx=idx+1
         spec%tau1=x(idx);idx=idx+1;spec%tau2=x(idx);idx=idx+1
         spec%measurement_sd=exp(max(-8.0_dp,min(3.0_dp,x(idx))));idx=idx+1
      case(model_fgarch)
         ncoef=context%p+context%q;allocate(raw(ncoef),coef(ncoef));raw=x(idx:idx+ncoef-1);idx=idx+ncoef
         call positive_coefficients(raw,0.95_dp,coef);spec%alpha=coef(1:context%p);spec%beta=coef(context%p+1:ncoef)
         select case(context%fgarch_submodel)
         case(fgarch_tgarch,fgarch_gjrgarch)
            do j=1,context%p;spec%eta1(j)=0.999_dp*tanh(x(idx));idx=idx+1;end do
         case(fgarch_avgarch)
            do j=1,context%p;spec%eta1(j)=0.999_dp*tanh(x(idx));idx=idx+1;end do
            do j=1,context%p;spec%eta2(j)=10.0_dp*tanh(x(idx));idx=idx+1;end do
         case(fgarch_ngarch)
            spec%fgarch_lambda=0.01_dp+3.99_dp*logistic(x(idx));idx=idx+1
         case(fgarch_nagarch)
            do j=1,context%p;spec%eta2(j)=10.0_dp*tanh(x(idx));idx=idx+1;end do
         case(fgarch_aparch)
            do j=1,context%p;spec%eta1(j)=0.999_dp*tanh(x(idx));idx=idx+1;end do
            spec%fgarch_lambda=0.01_dp+3.99_dp*logistic(x(idx));idx=idx+1
         case(fgarch_allgarch)
            do j=1,context%p;spec%eta1(j)=0.999_dp*tanh(x(idx));idx=idx+1;end do
            do j=1,context%p;spec%eta2(j)=10.0_dp*tanh(x(idx));idx=idx+1;end do
            spec%fgarch_lambda=0.01_dp+3.99_dp*logistic(x(idx));idx=idx+1
         end select
      case default
         ncoef=context%p+context%q;allocate(raw(ncoef),coef(ncoef));raw=x(idx:idx+ncoef-1)
         call positive_coefficients(raw,0.985_dp,coef);spec%alpha=coef(1:context%p);spec%beta=coef(context%p+1:ncoef);idx=idx+ncoef
      end select
      spec%shape=default_shape(context%cond_dist);spec%skew=default_skew(context%cond_dist);spec%lambda=1.0_dp
      if(context%fit_shape)then;spec%shape=decode_shape(x(idx),context%cond_dist);idx=idx+1;end if
      if(context%fit_skew)then;spec%skew=decode_skew(x(idx),context%cond_dist);idx=idx+1;end if
      if(context%fit_lambda)spec%lambda=4.0_dp*tanh(x(idx))
   end subroutine decode_garch_parameters

   pure integer function model_parameter_count(context) result(n)
      type(garch_fit_context),intent(in)::context
      select case(context%model)
      case(model_egarch);n=1+2*context%p+context%q
      case(model_gjrgarch);n=1+2*context%p+context%q
      case(model_aparch);n=2+2*context%p+context%q
      case(model_igarch);n=1+context%p+context%q
      case(model_figarch);n=2+context%p+context%q
      case(model_csgarch);n=3+context%p+context%q
      case(model_realgarch);n=6+context%p+context%q
      case(model_fgarch)
         n=1+context%p+context%q
         select case(context%fgarch_submodel)
         case(fgarch_tgarch,fgarch_gjrgarch);n=n+context%p
         case(fgarch_avgarch);n=n+2*context%p
         case(fgarch_ngarch);n=n+1
         case(fgarch_nagarch);n=n+context%p
         case(fgarch_aparch);n=n+context%p+1
         case(fgarch_allgarch);n=n+2*context%p+1
         end select
      case default;n=1+context%p+context%q
      end select
   end function model_parameter_count

   subroutine positive_coefficients(raw,total,coef)
      real(dp),intent(in)::raw(:),total;real(dp),intent(out)::coef(size(raw))
      real(dp)::m,den
      m=max(0.0_dp,maxval(raw));coef=exp(max(-50.0_dp,min(50.0_dp,raw-m)));den=exp(-m)+sum(coef)
      coef=total*coef/den
   end subroutine positive_coefficients

   subroutine unit_coefficients(raw,coef)
      real(dp),intent(in)::raw(:);real(dp),intent(out)::coef(size(raw));real(dp)::m
      m=maxval(raw);coef=exp(max(-50.0_dp,min(50.0_dp,raw-m)));coef=coef/sum(coef)
   end subroutine unit_coefficients

   pure elemental function logistic(x) result(value)
      real(dp),intent(in)::x;real(dp)::value
      if(x>=0.0_dp)then;value=1.0_dp/(1.0_dp+exp(-min(x,50.0_dp)))
      else;value=exp(max(x,-50.0_dp))/(1.0_dp+exp(max(x,-50.0_dp)));end if
   end function logistic

   pure elemental function logit(p) result(value)
      real(dp),intent(in)::p;real(dp)::value,q
      q=max(1.0e-12_dp,min(1.0_dp-1.0e-12_dp,p));value=log(q/(1.0_dp-q))
   end function logit

   pure logical function has_shape(kind)
      integer,intent(in)::kind
      has_shape=kind==dist_std.or.kind==dist_sstd.or.kind==dist_ged.or.kind==dist_sged.or.kind==dist_jsu.or. &
         kind==dist_nig.or.kind==dist_ghyp.or.kind==dist_ghst
   end function has_shape

   pure logical function has_skew(kind)
      integer,intent(in)::kind
      has_skew=kind==dist_snorm.or.kind==dist_sstd.or.kind==dist_sged.or.kind==dist_jsu.or.kind==dist_nig.or. &
         kind==dist_ghyp.or.kind==dist_ghst
   end function has_skew

   pure logical function has_lambda(kind)
      integer,intent(in)::kind;has_lambda=kind==dist_ghyp
   end function has_lambda

   pure function default_shape(kind) result(value)
      integer,intent(in)::kind;real(dp)::value
      select case(kind)
      case(dist_std,dist_sstd);value=6.0_dp
      case(dist_ghst);value=8.0_dp
      case(dist_nig,dist_ghyp);value=1.0_dp
      case(dist_jsu);value=2.0_dp
      case default;value=2.0_dp
      end select
   end function default_shape

   pure function default_skew(kind) result(value)
      integer,intent(in)::kind;real(dp)::value
      select case(kind)
      case(dist_snorm,dist_sstd,dist_sged);value=1.0_dp
      case default;value=0.0_dp
      end select
   end function default_skew

   pure function initial_shape_parameter(kind) result(value)
      integer,intent(in)::kind;real(dp)::value
      select case(kind)
      case(dist_std,dist_sstd);value=log(6.0_dp-2.05_dp)
      case(dist_ghst);value=log(8.0_dp-4.001_dp)
      case(dist_nig,dist_ghyp);value=log(1.0_dp-0.05_dp)
      case(dist_jsu);value=log(2.0_dp-0.05_dp)
      case default;value=log(2.0_dp-0.2_dp)
      end select
   end function initial_shape_parameter

   pure function initial_lambda_parameter() result(value)
      real(dp)::value
      value=atanh(0.25_dp)
   end function initial_lambda_parameter

   pure function initial_skew_parameter(kind) result(value)
      integer,intent(in)::kind;real(dp)::value
      select case(kind)
      case(dist_snorm,dist_sstd,dist_sged);value=0.0_dp
      case default;value=0.0_dp
      end select
   end function initial_skew_parameter

   pure function decode_shape(x,kind) result(value)
      real(dp),intent(in)::x;integer,intent(in)::kind;real(dp)::value
      select case(kind)
      case(dist_std,dist_sstd);value=2.05_dp+exp(max(-20.0_dp,min(7.0_dp,x)))
      case(dist_ghst);value=4.001_dp+exp(max(-20.0_dp,min(5.0_dp,x)))
      case(dist_nig,dist_ghyp);value=0.05_dp+exp(max(-20.0_dp,min(5.0_dp,x)))
      case(dist_jsu);value=0.05_dp+exp(max(-20.0_dp,min(5.0_dp,x)))
      case default;value=0.2_dp+exp(max(-20.0_dp,min(5.0_dp,x)))
      end select
   end function decode_shape

   pure function decode_skew(x,kind) result(value)
      real(dp),intent(in)::x;integer,intent(in)::kind;real(dp)::value
      select case(kind)
      case(dist_snorm,dist_sstd,dist_sged);value=exp(max(-5.0_dp,min(5.0_dp,x)))
      case(dist_nig,dist_ghyp);value=0.999_dp*tanh(x)
      case(dist_ghst);value=15.0_dp*tanh(x)
      case default;value=x
      end select
   end function decode_skew

   pure elemental logical function finite_value(x) result(value)
      real(dp),intent(in)::x;value=ieee_is_finite(x)
   end function finite_value

end module rugarch_fit
