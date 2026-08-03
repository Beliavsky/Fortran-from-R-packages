! Remaining standalone numerical workflows from rugarch.
! R object infrastructure, plotting, parallel rolling, and resumable rolling
! are intentionally outside the scope of this module.
module rugarch_complete
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
   use rugarch_kinds, only : dp
   use rugarch_types, only : garch_spec, garch_fit_result, model_sgarch, model_gjrgarch, &
      model_egarch, model_aparch, model_igarch, model_figarch, model_csgarch, &
      model_realgarch, model_fgarch
   use rugarch_fit, only : fit_model
   use rugarch_models, only : garch_filter, garch_log_likelihood, realgarch_filter, &
      realgarch_log_likelihood, forecast_volatility, simulate_garch, simulate_realgarch, &
      true_persistence, unconditional_variance
   use rugarch_distributions, only : distribution_pdf, distribution_cdf, &
      random_innovation, dist_norm, dist_snorm, dist_std, dist_sstd, dist_ged, &
      dist_sged, dist_jsu, dist_nig, dist_ghyp, dist_ghst
   use rugarch_gh, only : gh_parameters
   use rugarch_arfima, only : arfima_spec, arfima_fit_result, arfima_forecast_result, &
      make_arfima_spec, fit_arfima, forecast_arfima, simulate_arfima, arma_residuals, &
      fractional_weights
   use rugarch_backtests, only : berkowitz_test, berkowitz_result
   use rugarch_linalg, only : invert_matrix
   use rugarch_rng, only : random_normal
   implicit none
   private

   integer, parameter, public :: sampling_raw = 1
   integer, parameter, public :: sampling_kernel = 2
   integer, parameter, public :: sampling_spd = 3
   integer, parameter, public :: bootstrap_partial = 1
   integer, parameter, public :: bootstrap_full = 2

   type, public :: extended_garch_fit_result
      type(garch_fit_result) :: fit
      real(dp) :: intercept = 0.0_dp
      real(dp) :: archm = 0.0_dp
      real(dp) :: archm_power = 1.0_dp
      logical :: variance_targeting = .false.
      real(dp) :: target_variance = 0.0_dp
      real(dp), allocatable :: mean_beta(:)
      real(dp), allocatable :: variance_beta(:)
      real(dp), allocatable :: covariance(:,:)
      real(dp), allocatable :: robust_covariance(:,:)
      real(dp), allocatable :: standard_error(:)
      real(dp), allocatable :: robust_standard_error(:)
      character(len=16), allocatable :: parameter_name(:)
      integer :: covariance_status = 1
   end type extended_garch_fit_result

   type, public :: bootstrap_forecast_result
      real(dp), allocatable :: return_path(:,:)
      real(dp), allocatable :: sigma_path(:,:)
      real(dp), allocatable :: mean(:)
      real(dp), allocatable :: sigma(:)
      real(dp), allocatable :: lower(:)
      real(dp), allocatable :: upper(:)
      integer, allocatable :: status(:)
   end type bootstrap_forecast_result

   type, public :: arfima_rolling_result
      real(dp), allocatable :: actual(:), mean(:), sigma(:)
      integer, allocatable :: status(:), fit_end(:)
   end type arfima_rolling_result

   type, public :: arfima_distribution_result
      real(dp), allocatable :: mean(:), d(:), innovation_sd(:), css(:)
      real(dp), allocatable :: ar(:,:), ma(:,:)
      integer, allocatable :: status(:)
   end type arfima_distribution_result

   type, public :: arfima_bootstrap_result
      real(dp), allocatable :: path(:,:), mean(:), sigma(:), lower(:), upper(:)
   end type arfima_bootstrap_result

   type, public :: arfima_multi_fit_result
      type(arfima_fit_result), allocatable :: fit(:)
      integer, allocatable :: status(:)
   end type arfima_multi_fit_result

   type, public :: arfima_multi_forecast_result
      real(dp), allocatable :: mean(:,:), sigma(:,:)
      integer, allocatable :: status(:)
   end type arfima_multi_forecast_result

   type, public :: arfima_cv_result
      integer, allocatable :: p(:), q(:)
      real(dp), allocatable :: rmse(:), mae(:), berkowitz_p(:), criterion(:)
      integer :: best = 0
   end type arfima_cv_result

   type, public :: forecast_performance_result
      real(dp) :: me = 0.0_dp
      real(dp) :: mae = 0.0_dp
      real(dp) :: mse = 0.0_dp
      real(dp) :: rmse = 0.0_dp
      real(dp) :: mape = 0.0_dp
      real(dp) :: mase = 0.0_dp
      real(dp) :: directional_accuracy = 0.0_dp
   end type forecast_performance_result

   type, public :: distribution_moment_result
      real(dp) :: mean = 0.0_dp
      real(dp) :: variance = 1.0_dp
      real(dp) :: skewness = 0.0_dp
      real(dp) :: excess_kurtosis = 0.0_dp
      integer :: status = 0
   end type distribution_moment_result

   type, public :: gh_parameter_result
      real(dp) :: alpha = 0.0_dp
      real(dp) :: beta = 0.0_dp
      real(dp) :: delta = 0.0_dp
      real(dp) :: mu = 0.0_dp
      real(dp) :: lambda = 0.0_dp
      integer :: status = 0
   end type gh_parameter_result

   type, public :: numeric_index_result
      real(dp), allocatable :: value(:)
      integer :: status = 0
   end type numeric_index_result

   public :: fit_garch_extended, filter_garch_extended, attach_garch_covariance
   public :: garch_bootstrap_forecast
   public :: rolling_arfima_forecast, arfima_parametric_distribution
   public :: arfima_bootstrap_forecast, multifit_arfima, multiforecast_arfima
   public :: arfima_cross_validation
   public :: forecast_performance, unconditional_mean, volatility_half_life
   public :: distribution_moments, distribution_skewness, distribution_excess_kurtosis
   public :: ghyp_transform
   public :: move_numeric_index, generate_forward_numeric, filtered_numeric_sequence

contains

   function fit_garch_extended(y,model,p,q,cond_dist,mean_regressors,variance_regressors, &
      variance_targeting,target_variance,arch_in_mean,archm_power,fit_shape,fit_skew, &
      fit_lambda,max_iterations,realized,attach_covariance,robust_lags) result(ans)
      real(dp), intent(in) :: y(:)
      integer, intent(in) :: model,p,q
      integer, intent(in), optional :: cond_dist,max_iterations,robust_lags
      real(dp), intent(in), optional :: mean_regressors(:,:),variance_regressors(:,:),realized(:)
      real(dp), intent(in), optional :: target_variance,archm_power
      logical, intent(in), optional :: variance_targeting,arch_in_mean,fit_shape,fit_skew, &
         fit_lambda,attach_covariance
      type(extended_garch_fit_result) :: ans
      real(dp), allocatable :: x(:,:),beta(:),adjusted(:),base_sigma(:),target(:),vbeta(:)
      real(dp) :: tv,power,persist,mean_vreg
      integer :: n,k,info,nlag
      logical :: vt,aim,do_cov

      n=size(y); vt=.false.; if(present(variance_targeting))vt=variance_targeting
      aim=.false.; if(present(arch_in_mean))aim=arch_in_mean
      power=1.0_dp; if(present(archm_power))power=archm_power
      do_cov=.true.; if(present(attach_covariance))do_cov=attach_covariance
      nlag=0; if(present(robust_lags))nlag=max(0,robust_lags)
      allocate(adjusted(n));adjusted=y

      if(present(mean_regressors))then
         if(size(mean_regressors,1)/=n)then
            ans%fit%status=2;ans%fit%message='mean regressors have wrong row count';return
         end if
         k=size(mean_regressors,2);allocate(x(n,k+1),beta(k+1))
         x(:,1)=1.0_dp;if(k>0)x(:,2:)=mean_regressors
         call least_squares(x,y,beta,info)
         if(info/=0)then;ans%fit%status=2;ans%fit%message='singular mean-regressor design';return;end if
         ans%intercept=beta(1);allocate(ans%mean_beta(k));if(k>0)ans%mean_beta=beta(2:)
         adjusted=y-matmul(x,beta)
      else
         ans%intercept=sum(y)/real(max(1,n),dp)
         allocate(ans%mean_beta(0));adjusted=y-ans%intercept
      end if

      ! A two-step ARCH-in-mean estimate: preliminary volatility, followed by OLS
      ! with sigma**power as an additional regressor, then final volatility fit.
      if(aim)then
         if(present(realized))then
            ans%fit=fit_model(adjusted,model,p,q,cond_dist,.false.,fit_shape,fit_skew,fit_lambda, &
               realized=realized,max_iterations=max_iterations)
         else
            ans%fit=fit_model(adjusted,model,p,q,cond_dist,.false.,fit_shape,fit_skew,fit_lambda, &
               max_iterations=max_iterations)
         end if
         if(allocated(ans%fit%sigma))then
            if(allocated(x))deallocate(x,beta)
            k=0;if(present(mean_regressors))k=size(mean_regressors,2)
            allocate(x(n,k+2),beta(k+2));x(:,1)=1.0_dp
            if(k>0)x(:,2:k+1)=mean_regressors
            x(:,k+2)=ans%fit%sigma**power
            call least_squares(x,y,beta,info)
            if(info==0)then
               ans%intercept=beta(1)
               if(k>0)ans%mean_beta=beta(2:k+1)
               ans%archm=beta(k+2);ans%archm_power=power
               adjusted=y-matmul(x,beta)
            end if
         end if
      end if

      if(present(realized))then
         ans%fit=fit_model(adjusted,model,p,q,cond_dist,.false.,fit_shape,fit_skew,fit_lambda, &
            realized=realized,max_iterations=max_iterations)
      else
         ans%fit=fit_model(adjusted,model,p,q,cond_dist,.false.,fit_shape,fit_skew,fit_lambda, &
            max_iterations=max_iterations)
      end if

      ans%variance_targeting=vt
      persist=true_persistence(ans%fit%spec)
      if(vt .and. ans%fit%status<=1)then
         tv=sum(adjusted*adjusted)/real(max(1,n),dp);if(present(target_variance))tv=max(target_variance,1.0e-12_dp)
         ans%target_variance=tv
         select case(model)
         case(model_egarch)
            ans%fit%spec%omega=(1.0_dp-sum_safe(ans%fit%spec%beta))*log(tv)
         case(model_aparch)
            ans%fit%spec%omega=max(1.0e-12_dp,(tv**(0.5_dp*ans%fit%spec%delta))*(1.0_dp-persist))
         case(model_fgarch)
            ans%fit%spec%omega=max(1.0e-12_dp,(tv**(0.5_dp*ans%fit%spec%fgarch_lambda))*(1.0_dp-persist))
         case(model_realgarch)
            ans%fit%spec%omega=(1.0_dp-sum_safe(ans%fit%spec%beta)- &
               ans%fit%spec%phi*sum_safe(ans%fit%spec%alpha))*log(tv)- &
               ans%fit%spec%xi*sum_safe(ans%fit%spec%alpha)
         case default
            ans%fit%spec%omega=max(1.0e-12_dp,tv*(1.0_dp-persist))
         end select
         call refilter_fit(adjusted,realized,ans%fit)
      end if

      allocate(ans%variance_beta(0))
      if(present(variance_regressors) .and. allocated(ans%fit%sigma))then
         if(size(variance_regressors,1)/=n)then
            ans%fit%status=2;ans%fit%message='variance regressors have wrong row count';return
         end if
         k=size(variance_regressors,2)
         if(k>0)then
            allocate(target(n),vbeta(k));target=max(adjusted**2-ans%fit%sigma**2,0.0_dp)
            call least_squares(variance_regressors,target,vbeta,info)
            if(info==0)then
               vbeta=max(vbeta,0.0_dp);deallocate(ans%variance_beta);allocate(ans%variance_beta(k));ans%variance_beta=vbeta
               mean_vreg=sum(matmul(variance_regressors,vbeta))/real(n,dp)
               if(vt .and. model/=model_egarch .and. model/=model_realgarch) &
                  ans%fit%spec%omega=max(1.0e-12_dp,ans%fit%spec%omega-mean_vreg*(1.0_dp-persist))
               allocate(base_sigma(n));base_sigma=sqrt(max(ans%fit%sigma**2+matmul(variance_regressors,vbeta),1.0e-20_dp))
               ans%fit%sigma=base_sigma
               ans%fit%log_likelihood=conditional_loglik(adjusted,ans%fit%sigma,ans%fit%spec)
            end if
         end if
      end if

      if(do_cov .and. ans%fit%status<=1)call attach_garch_covariance(adjusted,ans%fit,ans,realized,nlag)
   end function fit_garch_extended

   subroutine filter_garch_extended(y,fit,mean_regressors,variance_regressors,residuals,sigma,valid)
      real(dp),intent(in)::y(:)
      type(extended_garch_fit_result),intent(in)::fit
      real(dp),intent(in),optional::mean_regressors(:,:),variance_regressors(:,:)
      real(dp),intent(out)::residuals(size(y)),sigma(size(y))
      logical,intent(out),optional::valid
      real(dp),allocatable::adj(:)
      logical::ok
      allocate(adj(size(y)));adj=y-fit%intercept
      if(present(mean_regressors).and.allocated(fit%mean_beta))then
         if(size(mean_regressors,1)==size(y).and.size(mean_regressors,2)==size(fit%mean_beta)) &
            adj=adj-matmul(mean_regressors,fit%mean_beta)
      end if
      call garch_filter(adj,fit%fit%spec,residuals,sigma,ok)
      if(abs(fit%archm)>tiny(1.0_dp))then
         adj=adj-fit%archm*sigma**fit%archm_power
         call garch_filter(adj,fit%fit%spec,residuals,sigma,ok)
      end if
      if(present(variance_regressors).and.allocated(fit%variance_beta))then
         if(size(variance_regressors,1)==size(y).and.size(variance_regressors,2)==size(fit%variance_beta)) &
            sigma=sqrt(max(sigma**2+matmul(variance_regressors,fit%variance_beta),1.0e-20_dp))
      end if
      if(present(valid))valid=ok
   end subroutine filter_garch_extended

   subroutine attach_garch_covariance(y,fit,ans,realized,newey_west_lags)
      real(dp),intent(in)::y(:)
      type(garch_fit_result),intent(in)::fit
      type(extended_garch_fit_result),intent(inout)::ans
      real(dp),intent(in),optional::realized(:)
      integer,intent(in),optional::newey_west_lags
      real(dp),allocatable::theta(:),hessian(:,:),hinv(:,:),scores(:,:),meat(:,:),tmp(:,:)
      character(len=16),allocatable::names(:)
      integer::m,n,i,j,info,lags
      real(dp)::hi,hj,f0,fip,fim,fpp,fpm,fmp,fmm
      type(garch_spec)::sp,sm,spp,spm,smp,smm

      call pack_spec(fit%spec,theta,names);m=size(theta);n=size(y)
      allocate(hessian(m,m),scores(n,m));hessian=0.0_dp;scores=0.0_dp
      f0=spec_loglik(y,fit%spec,realized)
      do i=1,m
         hi=max(1.0e-5_dp,1.0e-4_dp*(1.0_dp+abs(theta(i))))
         call perturbed_spec(fit%spec,i,hi,sp);call perturbed_spec(fit%spec,i,-hi,sm)
         fip=spec_loglik(y,sp,realized);fim=spec_loglik(y,sm,realized)
         if(finite_loglik(fip).and.finite_loglik(fim).and.finite_loglik(f0)) &
            hessian(i,i)=-(fip-2.0_dp*f0+fim)/(hi*hi)
         call score_column(y,fit%spec,i,hi,scores(:,i),realized)
         do j=i+1,m
            hj=max(1.0e-5_dp,1.0e-4_dp*(1.0_dp+abs(theta(j))))
            call double_perturbed_spec(fit%spec,i,hi,j,hj,spp)
            call double_perturbed_spec(fit%spec,i,hi,j,-hj,spm)
            call double_perturbed_spec(fit%spec,i,-hi,j,hj,smp)
            call double_perturbed_spec(fit%spec,i,-hi,j,-hj,smm)
            fpp=spec_loglik(y,spp,realized);fpm=spec_loglik(y,spm,realized)
            fmp=spec_loglik(y,smp,realized);fmm=spec_loglik(y,smm,realized)
            if(finite_loglik(fpp).and.finite_loglik(fpm).and.finite_loglik(fmp).and.finite_loglik(fmm))then
               hessian(i,j)=-(fpp-fpm-fmp+fmm)/(4.0_dp*hi*hj);hessian(j,i)=hessian(i,j)
            end if
         end do
      end do
      allocate(hinv(m,m))
      call invert_matrix(hessian,hinv,info)
      allocate(ans%parameter_name(m));ans%parameter_name=names
      if(info/=0)then;ans%covariance_status=2;return;end if
      allocate(ans%covariance(m,m),ans%standard_error(m));ans%covariance=hinv
      ans%standard_error=sqrt(max([(hinv(i,i),i=1,m)],0.0_dp))
      lags=0;if(present(newey_west_lags))lags=max(0,newey_west_lags)
      call score_long_run_covariance(scores,lags,meat)
      allocate(tmp(m,m),ans%robust_covariance(m,m),ans%robust_standard_error(m))
      tmp=matmul(hinv,meat);ans%robust_covariance=matmul(tmp,hinv)
      ans%robust_standard_error=sqrt(max([(ans%robust_covariance(i,i),i=1,m)],0.0_dp))
      ans%covariance_status=0
   end subroutine attach_garch_covariance

   function garch_bootstrap_forecast(fit,n_ahead,n_paths,sampling,method,n_bootfit, &
      confidence,max_iterations) result(ans)
      type(garch_fit_result),intent(in)::fit
      integer,intent(in)::n_ahead,n_paths
      integer,intent(in),optional::sampling,method,n_bootfit,max_iterations
      real(dp),intent(in),optional::confidence
      type(bootstrap_forecast_result)::ans
      type(garch_fit_result),allocatable::pool(:)
      type(garch_fit_result)::selected
      real(dp),allocatable::z(:),eps(:),vol(:),one(:),ysim(:),ssim(:),esim(:),realized(:)
      real(dp)::u,level
      integer::np,h,b,n,sm,bt,nfit,idx,p,q

      np=max(1,n_paths);sm=sampling_raw;if(present(sampling))sm=sampling
      bt=bootstrap_partial;if(present(method))bt=method
      nfit=max(10,min(100,np));if(present(n_bootfit))nfit=max(1,n_bootfit)
      level=0.95_dp;if(present(confidence))level=max(0.5_dp,min(0.999_dp,confidence))
      n=size(fit%residuals);p=size(fit%spec%alpha);q=size(fit%spec%beta)
      allocate(ans%return_path(n_ahead,np),ans%sigma_path(n_ahead,np),ans%mean(n_ahead), &
         ans%sigma(n_ahead),ans%lower(n_ahead),ans%upper(n_ahead),ans%status(np))
      allocate(z(n));z=fit%residuals/max(fit%sigma,1.0e-20_dp);z=(z-sum(z)/real(n,dp))/sample_sd(z)

      if(bt==bootstrap_full)then
         allocate(pool(nfit),ysim(n),ssim(n),esim(n))
         if(fit%spec%model==model_realgarch)allocate(realized(n))
         do b=1,nfit
            if(fit%spec%model==model_realgarch)then
               call simulate_realgarch(fit%spec,n,ysim,realized,ssim,esim,burn_in=200)
               pool(b)=fit_model(ysim,fit%spec%model,p,q,fit%spec%cond_dist,.true.,.true.,.true.,.true., &
                  realized=realized,max_iterations=max_iterations)
            else
               call simulate_garch(fit%spec,n,ysim,ssim,esim,burn_in=200)
               pool(b)=fit_model(ysim,fit%spec%model,p,q,fit%spec%cond_dist,.true.,.true.,.true.,.true., &
                  max_iterations=max_iterations,fgarch_submodel=fit%spec%fgarch_submodel, &
                  figarch_truncation=fit%spec%figarch_truncation)
            end if
         end do
      end if

      do b=1,np
         selected=fit
         if(bt==bootstrap_full)then
            call random_number(u);idx=1+int(u*real(nfit,dp));idx=min(nfit,max(1,idx))
            if(pool(idx)%status<=1.and.allocated(pool(idx)%residuals))selected=pool(idx)
         end if
         allocate(eps(size(selected%residuals)+n_ahead),vol(size(selected%sigma)+n_ahead),one(1))
         eps(1:size(selected%residuals))=selected%residuals;vol(1:size(selected%sigma))=selected%sigma
         do h=1,n_ahead
            call forecast_volatility(selected%spec,eps(1:size(selected%residuals)+h-1), &
               vol(1:size(selected%sigma)+h-1),1,one)
            vol(size(selected%sigma)+h)=one(1)
            eps(size(selected%residuals)+h)=one(1)*sample_standardized(z,sm)
            ans%sigma_path(h,b)=one(1)
            ans%return_path(h,b)=selected%spec%mean+eps(size(selected%residuals)+h)
         end do
         ans%status(b)=selected%status;deallocate(eps,vol,one)
      end do
      do h=1,n_ahead
         ans%mean(h)=sum(ans%return_path(h,:))/real(np,dp)
         ans%sigma(h)=sample_sd(ans%return_path(h,:))
         ans%lower(h)=empirical_quantile(ans%return_path(h,:),(1.0_dp-level)/2.0_dp)
         ans%upper(h)=empirical_quantile(ans%return_path(h,:),1.0_dp-(1.0_dp-level)/2.0_dp)
      end do
   end function garch_bootstrap_forecast

   function rolling_arfima_forecast(data,p,q,forecast_length,refit_every,window_size, &
      moving_window,estimate_d) result(ans)
      real(dp),intent(in)::data(:)
      integer,intent(in)::p,q,forecast_length
      integer,intent(in),optional::refit_every,window_size
      logical,intent(in),optional::moving_window,estimate_d
      type(arfima_rolling_result)::ans
      type(arfima_fit_result)::fit
      type(arfima_forecast_result)::fc
      integer::n,nf,first,i,end_fit,start_fit,refit,win
      logical::moving,need_fit,use_d
      n=size(data);nf=min(max(1,forecast_length),n-20);first=n-nf+1
      refit=1;if(present(refit_every))refit=max(1,refit_every)
      win=n;if(present(window_size))win=max(20,window_size)
      moving=.false.;if(present(moving_window))moving=moving_window
      use_d=.true.;if(present(estimate_d))use_d=estimate_d
      allocate(ans%actual(nf),ans%mean(nf),ans%sigma(nf),ans%status(nf),ans%fit_end(nf))
      ans%actual=data(first:n);ans%mean=0.0_dp;ans%sigma=0.0_dp;ans%status=1;ans%fit_end=0
      need_fit=.true.
      do i=1,nf
         end_fit=first+i-2;if(i==1.or.mod(i-1,refit)==0)need_fit=.true.
         if(need_fit)then
            start_fit=1;if(moving)start_fit=max(1,end_fit-win+1)
            fit=fit_arfima(data(start_fit:end_fit),p,q,use_d);need_fit=.false.
         end if
         ans%fit_end(i)=end_fit;ans%status(i)=fit%status
         if(fit%status==0)then
            fc=forecast_arfima(data(max(1,end_fit-win+1):end_fit),fit%spec,1)
            ans%mean(i)=fc%mean(1);ans%sigma(i)=fc%sigma(1)
         end if
      end do
   end function rolling_arfima_forecast

   function arfima_parametric_distribution(fit,nboot,nobs) result(ans)
      type(arfima_fit_result),intent(in)::fit
      integer,intent(in)::nboot
      integer,intent(in),optional::nobs
      type(arfima_distribution_result)::ans
      type(arfima_fit_result)::bf
      real(dp),allocatable::x(:)
      integer::b,n,p,q
      n=1000;if(present(nobs))n=max(30,nobs);p=size(fit%spec%ar);q=size(fit%spec%ma)
      allocate(ans%mean(nboot),ans%d(nboot),ans%innovation_sd(nboot),ans%css(nboot), &
         ans%ar(p,nboot),ans%ma(q,nboot),ans%status(nboot),x(n))
      if(p>0)ans%ar=0.0_dp;if(q>0)ans%ma=0.0_dp
      do b=1,nboot
         call simulate_arfima(fit%spec,n,x,burn_in=300);bf=fit_arfima(x,p,q,.true.)
         ans%mean(b)=bf%spec%mean;ans%d(b)=bf%spec%d;ans%innovation_sd(b)=bf%spec%innovation_sd
         ans%css(b)=bf%css;ans%status(b)=bf%status
         if(p>0)ans%ar(:,b)=bf%spec%ar;if(q>0)ans%ma(:,b)=bf%spec%ma
      end do
   end function arfima_parametric_distribution

   function arfima_bootstrap_forecast(history,fit,horizon,n_paths,sampling,confidence) result(ans)
      real(dp),intent(in)::history(:)
      type(arfima_fit_result),intent(in)::fit
      integer,intent(in)::horizon,n_paths
      integer,intent(in),optional::sampling
      real(dp),intent(in),optional::confidence
      type(arfima_bootstrap_result)::ans
      type(arfima_forecast_result)::fc
      real(dp),allocatable::w(:),innov(:)
      real(dp)::level
      integer::b,h,sm
      sm=sampling_raw;if(present(sampling))sm=sampling
      level=0.95_dp;if(present(confidence))level=max(0.5_dp,min(0.999_dp,confidence))
      fc=forecast_arfima(history,fit%spec,horizon);allocate(w(horizon),innov(horizon))
      call fractional_weights(fit%spec%d,w,inverse=.true.)
      allocate(ans%path(horizon,n_paths),ans%mean(horizon),ans%sigma(horizon), &
         ans%lower(horizon),ans%upper(horizon))
      do b=1,n_paths
         do h=1,horizon
            innov(h)=fit%spec%innovation_sd*sample_standardized(fit%residuals/safe_sd(fit%residuals),sm)
            ans%path(h,b)=fc%mean(h)+sum(w(1:h)*innov(h:1:-1))
         end do
      end do
      do h=1,horizon
         ans%mean(h)=sum(ans%path(h,:))/real(n_paths,dp);ans%sigma(h)=sample_sd(ans%path(h,:))
         ans%lower(h)=empirical_quantile(ans%path(h,:),(1.0_dp-level)/2.0_dp)
         ans%upper(h)=empirical_quantile(ans%path(h,:),1.0_dp-(1.0_dp-level)/2.0_dp)
      end do
   end function arfima_bootstrap_forecast

   function multifit_arfima(data,p,q,estimate_d) result(ans)
      real(dp),intent(in)::data(:,:)
      integer,intent(in)::p,q
      logical,intent(in),optional::estimate_d
      type(arfima_multi_fit_result)::ans
      integer::j
      allocate(ans%fit(size(data,2)),ans%status(size(data,2)))
      do j=1,size(data,2);ans%fit(j)=fit_arfima(data(:,j),p,q,estimate_d);ans%status(j)=ans%fit(j)%status;end do
   end function multifit_arfima

   function multiforecast_arfima(data,fits,horizon) result(ans)
      real(dp),intent(in)::data(:,:)
      type(arfima_multi_fit_result),intent(in)::fits
      integer,intent(in)::horizon
      type(arfima_multi_forecast_result)::ans
      type(arfima_forecast_result)::fc
      integer::j,m
      m=size(fits%fit);allocate(ans%mean(horizon,m),ans%sigma(horizon,m),ans%status(m))
      do j=1,m
         ans%status(j)=fits%fit(j)%status
         if(ans%status(j)==0)then
            fc=forecast_arfima(data(:,j),fits%fit(j)%spec,horizon)
            ans%mean(:,j)=fc%mean;ans%sigma(:,j)=fc%sigma
         else;ans%mean(:,j)=0.0_dp;ans%sigma(:,j)=0.0_dp;end if
      end do
   end function multiforecast_arfima

   function arfima_cross_validation(data,train_end,test_end,pmax,qmax,estimate_d) result(ans)
      real(dp),intent(in)::data(:)
      integer,intent(in)::train_end(:),test_end(:),pmax,qmax
      logical,intent(in),optional::estimate_d
      type(arfima_cv_result)::ans
      type(arfima_fit_result)::fit
      type(arfima_forecast_result)::fc
      type(berkowitz_result)::bk
      real(dp),allocatable::errors(:),z(:)
      integer::nc,p,q,k,idx,h,nerr
      logical::use_d
      use_d=.true.;if(present(estimate_d))use_d=estimate_d
      if(size(train_end)/=size(test_end))return
      nc=(max(0,pmax)+1)*(max(0,qmax)+1);allocate(ans%p(nc),ans%q(nc),ans%rmse(nc), &
         ans%mae(nc),ans%berkowitz_p(nc),ans%criterion(nc))
      idx=0;ans%best=1
      do p=0,max(0,pmax);do q=0,max(0,qmax)
         idx=idx+1;ans%p(idx)=p;ans%q(idx)=q;allocate(errors(sum(max(0,test_end-train_end))),z(sum(max(0,test_end-train_end))))
         nerr=0
         do k=1,size(train_end)
            if(train_end(k)<20.or.test_end(k)<=train_end(k).or.test_end(k)>size(data))cycle
            fit=fit_arfima(data(1:train_end(k)),p,q,use_d)
            fc=forecast_arfima(data(1:train_end(k)),fit%spec,test_end(k)-train_end(k))
            do h=1,size(fc%mean)
               nerr=nerr+1;errors(nerr)=data(train_end(k)+h)-fc%mean(h)
               z(nerr)=errors(nerr)/max(fc%sigma(h),1.0e-12_dp)
            end do
         end do
         if(nerr>0)then
            ans%rmse(idx)=sqrt(sum(errors(1:nerr)**2)/real(nerr,dp));ans%mae(idx)=sum(abs(errors(1:nerr)))/real(nerr,dp)
            if(nerr>5)then;bk=berkowitz_test(z(1:nerr));ans%berkowitz_p(idx)=bk%p_value;else;ans%berkowitz_p(idx)=0.0_dp;end if
         else;ans%rmse(idx)=huge(1.0_dp);ans%mae(idx)=huge(1.0_dp);ans%berkowitz_p(idx)=0.0_dp;end if
         ans%criterion(idx)=ans%rmse(idx)
         if(ans%criterion(idx)<ans%criterion(ans%best))ans%best=idx
         deallocate(errors,z)
      end do;end do
   end function arfima_cross_validation

   pure function forecast_performance(actual,forecast,insample) result(ans)
      real(dp),intent(in)::actual(:),forecast(:)
      real(dp),intent(in),optional::insample(:)
      type(forecast_performance_result)::ans
      real(dp)::e(min(size(actual),size(forecast))),den
      integer::n,i,correct
      n=size(e);if(n==0)return;e=forecast(1:n)-actual(1:n)
      ans%me=sum(e)/real(n,dp);ans%mae=sum(abs(e))/real(n,dp);ans%mse=sum(e*e)/real(n,dp);ans%rmse=sqrt(ans%mse)
      den=sum(abs(actual(1:n)))/real(n,dp);if(den>0.0_dp)ans%mape=100.0_dp*sum(abs(e)/max(abs(actual(1:n)),1.0e-12_dp))/real(n,dp)
      if(present(insample).and.size(insample)>1)then
         den=sum(abs(insample(2:)-insample(:size(insample)-1)))/real(size(insample)-1,dp)
         if(den>0.0_dp)ans%mase=ans%mae/den
      end if
      correct=0;do i=2,n;if((forecast(i)-actual(i-1))*(actual(i)-actual(i-1))>=0.0_dp)correct=correct+1;end do
      if(n>1)ans%directional_accuracy=real(correct,dp)/real(n-1,dp)
   end function forecast_performance

   function unconditional_mean(spec,mean_regressor_average,arch_in_mean,archm_power) result(value)
      type(garch_spec),intent(in)::spec
      real(dp),intent(in),optional::mean_regressor_average(:),arch_in_mean,archm_power
      real(dp)::value,den
      value=spec%mean;den=1.0_dp-sum_safe(spec%ar)
      if(present(mean_regressor_average))value=value+sum(mean_regressor_average)
      if(present(arch_in_mean))then
         if(present(archm_power))value=value+arch_in_mean*unconditional_variance(spec)**(0.5_dp*archm_power)
      end if
      value=value/max(den,1.0e-8_dp)
   end function unconditional_mean

   function volatility_half_life(spec) result(value)
      type(garch_spec),intent(in)::spec
      real(dp)::value,p
      p=true_persistence(spec)
      if(p>0.0_dp.and.p<1.0_dp)then;value=-log(2.0_dp)/log(p);else;value=huge(1.0_dp);end if
   end function volatility_half_life

   function distribution_moments(kind,shape,skew,lambda) result(ans)
      integer,intent(in)::kind
      real(dp),intent(in),optional::shape,skew,lambda
      type(distribution_moment_result)::ans
      real(dp)::sh,sk,la,a,b,h,x,pdf,m0,m1,m2,m3,m4
      integer::i,n,w
      sh=5.0_dp;if(present(shape))sh=shape;sk=1.0_dp;if(present(skew))sk=skew
      la=-0.5_dp;if(present(lambda))la=lambda
      if(kind==dist_norm)then;return;end if
      if(kind==dist_std.and.sh>4.0_dp)then;ans%excess_kurtosis=6.0_dp/(sh-4.0_dp);return;end if
      n=12000;a=-40.0_dp;b=40.0_dp;h=(b-a)/real(n,dp);m0=0.0_dp;m1=0.0_dp;m2=0.0_dp;m3=0.0_dp;m4=0.0_dp
      do i=0,n
         x=a+real(i,dp)*h;pdf=distribution_pdf(x,kind,sh,sk,la)
         if(i==0.or.i==n)then;w=1;else if(mod(i,2)==0)then;w=2;else;w=4;end if
         m0=m0+real(w,dp)*pdf;m1=m1+real(w,dp)*x*pdf;m2=m2+real(w,dp)*x*x*pdf
         m3=m3+real(w,dp)*x**3*pdf;m4=m4+real(w,dp)*x**4*pdf
      end do
      m0=m0*h/3.0_dp;m1=m1*h/3.0_dp;m2=m2*h/3.0_dp;m3=m3*h/3.0_dp;m4=m4*h/3.0_dp
      if(m0<=0.0_dp.or..not.ieee_is_finite(m0))then;ans%status=1;return;end if
      m1=m1/m0;m2=m2/m0-m1*m1
      ans%mean=m1;ans%variance=m2
      if(m2>0.0_dp)then
         ans%skewness=(m3/m0-3.0_dp*m1*(m2+m1*m1)+2.0_dp*m1**3)/m2**1.5_dp
         ans%excess_kurtosis=(m4/m0-4.0_dp*m1*(m3/m0)+6.0_dp*m1*m1*(m2+m1*m1)-3.0_dp*m1**4)/(m2*m2)-3.0_dp
      else;ans%status=2;end if
   end function distribution_moments

   function distribution_skewness(kind,shape,skew,lambda) result(value)
      integer,intent(in)::kind;real(dp),intent(in),optional::shape,skew,lambda
      real(dp)::value;type(distribution_moment_result)::m
      m=distribution_moments(kind,shape,skew,lambda);value=m%skewness
   end function distribution_skewness

   function distribution_excess_kurtosis(kind,shape,skew,lambda) result(value)
      integer,intent(in)::kind;real(dp),intent(in),optional::shape,skew,lambda
      real(dp)::value;type(distribution_moment_result)::m
      m=distribution_moments(kind,shape,skew,lambda);value=m%excess_kurtosis
   end function distribution_excess_kurtosis

   function ghyp_transform(mean,sd,skew,shape,lambda) result(ans)
      real(dp),intent(in)::mean,sd,skew,shape,lambda
      type(gh_parameter_result)::ans
      real(dp)::a,b,d,m
      logical::valid
      call gh_parameters(skew,shape,lambda,a,b,d,m,valid)
      if(.not.valid.or.sd<=0.0_dp)then;ans%status=1;return;end if
      ans%alpha=a/sd;ans%beta=b/sd;ans%delta=d*sd;ans%mu=m*sd+mean;ans%lambda=lambda
   end function ghyp_transform

   function move_numeric_index(index,by) result(ans)
      real(dp),intent(in)::index(:);integer,intent(in)::by
      type(numeric_index_result)::ans
      integer::n,k;n=size(index);k=max(0,by);allocate(ans%value(n))
      if(k>=n)then;ans%status=1;ans%value=index;return;end if
      if(k==0)then;ans%value=index;return;end if
      ans%value(1:n-k)=index(k+1:n);ans%value(n-k+1:n)=generate_linear(index(n),median_difference(index),k)
   end function move_numeric_index

   function generate_forward_numeric(t0,length_out,step) result(ans)
      real(dp),intent(in)::t0,step;integer,intent(in)::length_out
      type(numeric_index_result)::ans;integer::i
      allocate(ans%value(max(0,length_out)));do i=1,size(ans%value);ans%value(i)=t0+real(i,dp)*step;end do
   end function generate_forward_numeric

   function filtered_numeric_sequence(t0,length_out,step,period,interval_start,interval_end) result(ans)
      real(dp),intent(in)::t0,step,period,interval_start,interval_end;integer,intent(in)::length_out
      type(numeric_index_result)::ans
      real(dp)::x,phase;integer::n
      allocate(ans%value(max(0,length_out)));x=t0;n=0
      do while(n<length_out)
         x=x+step;phase=modulo(x,period)
         if(phase>=interval_start.and.phase<=interval_end)then;n=n+1;ans%value(n)=x;end if
      end do
   end function filtered_numeric_sequence

   subroutine refilter_fit(y,realized,fit)
      real(dp),intent(in)::y(:);real(dp),intent(in),optional::realized(:)
      type(garch_fit_result),intent(inout)::fit
      if(fit%spec%model==model_realgarch.and.present(realized))then
         if(.not.allocated(fit%measurement_residuals))allocate(fit%measurement_residuals(size(y)))
         fit%log_likelihood=realgarch_log_likelihood(y,realized,fit%spec,fit%residuals,fit%sigma,fit%measurement_residuals)
      else
         fit%log_likelihood=garch_log_likelihood(y,fit%spec,fit%residuals,fit%sigma)
      end if
   end subroutine refilter_fit

   function conditional_loglik(residuals,sigma,spec) result(value)
      real(dp),intent(in)::residuals(:),sigma(:);type(garch_spec),intent(in)::spec
      real(dp)::value,d;integer::i
      value=0.0_dp
      do i=1,size(residuals)
         d=distribution_pdf(residuals(i)/sigma(i),spec%cond_dist, &
            spec%shape,spec%skew,spec%lambda)/sigma(i)
         if(d<=tiny(1.0_dp))then
            value=-huge(1.0_dp)
            return
         end if
         value=value+log(d)
      end do
   end function conditional_loglik

   subroutine least_squares(x,y,beta,info)
      real(dp),intent(in)::x(:,:),y(:);real(dp),intent(out)::beta(size(x,2));integer,intent(out)::info
      real(dp),allocatable::xtx(:,:),inv(:,:)
      allocate(xtx(size(x,2),size(x,2)),inv(size(x,2),size(x,2)));xtx=matmul(transpose(x),x)
      call invert_matrix(xtx,inv,info);if(info==0)beta=matmul(inv,matmul(transpose(x),y));if(info/=0)beta=0.0_dp
   end subroutine least_squares

   subroutine pack_spec(spec,theta,names)
      type(garch_spec),intent(in)::spec
      real(dp),allocatable,intent(out)::theta(:)
      character(len=16),allocatable,intent(out)::names(:)
      integer::m,i,k

      m=2+size(spec%alpha)+size(spec%beta)
      select case(spec%model)
      case(model_gjrgarch,model_egarch)
         m=m+size(spec%gamma)
      case(model_aparch)
         m=m+size(spec%gamma)+1
      case(model_figarch)
         m=m+1
      case(model_csgarch)
         m=m+2
      case(model_realgarch)
         m=m+5
      case(model_fgarch)
         m=m+1+size(spec%eta1)+size(spec%eta2)
      end select
      if(has_shape_kind(spec%cond_dist))m=m+1
      if(has_skew_kind(spec%cond_dist))m=m+1
      if(spec%cond_dist==dist_ghyp)m=m+1

      allocate(theta(m),names(m));k=0
      call add(spec%mean,'mean')
      call add(spec%omega,'omega')
      do i=1,size(spec%alpha)
         call add(spec%alpha(i),'alpha')
      end do
      do i=1,size(spec%beta)
         call add(spec%beta(i),'beta')
      end do
      select case(spec%model)
      case(model_gjrgarch,model_egarch)
         do i=1,size(spec%gamma)
            call add(spec%gamma(i),'gamma')
         end do
      case(model_aparch)
         do i=1,size(spec%gamma)
            call add(spec%gamma(i),'gamma')
         end do
         call add(spec%delta,'delta')
      case(model_figarch)
         call add(spec%frac_d,'frac_d')
      case(model_csgarch)
         call add(spec%rho,'rho')
         call add(spec%phi,'phi')
      case(model_realgarch)
         call add(spec%xi,'xi')
         call add(spec%phi,'phi')
         call add(spec%tau1,'tau1')
         call add(spec%tau2,'tau2')
         call add(spec%measurement_sd,'meas_sd')
      case(model_fgarch)
         call add(spec%fgarch_lambda,'fg_lambda')
         do i=1,size(spec%eta1)
            call add(spec%eta1(i),'eta1')
         end do
         do i=1,size(spec%eta2)
            call add(spec%eta2(i),'eta2')
         end do
      end select
      if(has_shape_kind(spec%cond_dist))call add(spec%shape,'shape')
      if(has_skew_kind(spec%cond_dist))call add(spec%skew,'skew')
      if(spec%cond_dist==dist_ghyp)call add(spec%lambda,'lambda')
   contains
      subroutine add(value,name)
         real(dp),intent(in)::value
         character(len=*),intent(in)::name
         k=k+1
         theta(k)=value
         write(names(k),'(a,i0)')trim(name),k
      end subroutine add
   end subroutine pack_spec

   subroutine perturbed_spec(base,index,change,out)
      type(garch_spec),intent(in)::base
      integer,intent(in)::index
      real(dp),intent(in)::change
      type(garch_spec),intent(out)::out
      real(dp),allocatable::theta(:)
      character(len=16),allocatable::names(:)
      call pack_spec(base,theta,names)
      theta(index)=theta(index)+change
      call unpack_spec(base,theta,out)
   end subroutine perturbed_spec

   subroutine double_perturbed_spec(base,i,di,j,dj,out)
      type(garch_spec),intent(in)::base
      integer,intent(in)::i,j
      real(dp),intent(in)::di,dj
      type(garch_spec),intent(out)::out
      real(dp),allocatable::theta(:)
      character(len=16),allocatable::names(:)
      call pack_spec(base,theta,names)
      theta(i)=theta(i)+di
      theta(j)=theta(j)+dj
      call unpack_spec(base,theta,out)
   end subroutine double_perturbed_spec

   subroutine unpack_spec(base,theta,spec)
      type(garch_spec),intent(in)::base
      real(dp),intent(in)::theta(:)
      type(garch_spec),intent(out)::spec
      integer::i,k
      spec=base;k=0
      call get(spec%mean)
      call get(spec%omega)
      do i=1,size(spec%alpha)
         call get(spec%alpha(i))
      end do
      do i=1,size(spec%beta)
         call get(spec%beta(i))
      end do
      select case(spec%model)
      case(model_gjrgarch,model_egarch)
         do i=1,size(spec%gamma)
            call get(spec%gamma(i))
         end do
      case(model_aparch)
         do i=1,size(spec%gamma)
            call get(spec%gamma(i))
         end do
         call get(spec%delta)
      case(model_figarch)
         call get(spec%frac_d)
      case(model_csgarch)
         call get(spec%rho)
         call get(spec%phi)
      case(model_realgarch)
         call get(spec%xi)
         call get(spec%phi)
         call get(spec%tau1)
         call get(spec%tau2)
         call get(spec%measurement_sd)
      case(model_fgarch)
         call get(spec%fgarch_lambda)
         do i=1,size(spec%eta1)
            call get(spec%eta1(i))
         end do
         do i=1,size(spec%eta2)
            call get(spec%eta2(i))
         end do
      end select
      if(has_shape_kind(spec%cond_dist))call get(spec%shape)
      if(has_skew_kind(spec%cond_dist))call get(spec%skew)
      if(spec%cond_dist==dist_ghyp)call get(spec%lambda)
      spec%omega=max(spec%omega,1.0e-12_dp)
      spec%delta=max(spec%delta,0.05_dp)
      spec%shape=max(spec%shape,0.1_dp)
      spec%measurement_sd=max(spec%measurement_sd,1.0e-5_dp)
      spec%fgarch_lambda=max(spec%fgarch_lambda,0.05_dp)
   contains
      subroutine get(value)
         real(dp),intent(out)::value
         k=k+1
         value=theta(k)
      end subroutine get
   end subroutine unpack_spec

   function spec_loglik(y,spec,realized) result(value)
      real(dp),intent(in)::y(:);type(garch_spec),intent(in)::spec;real(dp),intent(in),optional::realized(:)
      real(dp)::value
      if(spec%model==model_realgarch.and.present(realized))then;value=realgarch_log_likelihood(y,realized,spec)
      else;value=garch_log_likelihood(y,spec);end if
   end function spec_loglik

   subroutine score_column(y,base,index,h,score,realized)
      real(dp),intent(in)::y(:),h;type(garch_spec),intent(in)::base;integer,intent(in)::index
      real(dp),intent(out)::score(size(y));real(dp),intent(in),optional::realized(:)
      type(garch_spec)::sp,sm;real(dp),allocatable::cp(:),cm(:)
      call perturbed_spec(base,index,h,sp);call perturbed_spec(base,index,-h,sm)
      allocate(cp(size(y)),cm(size(y)));call loglik_contributions(y,sp,cp,realized);call loglik_contributions(y,sm,cm,realized)
      score=(cp-cm)/(2.0_dp*h)
   end subroutine score_column

   subroutine loglik_contributions(y,spec,c,realized)
      real(dp),intent(in)::y(:);type(garch_spec),intent(in)::spec;real(dp),intent(out)::c(size(y))
      real(dp),intent(in),optional::realized(:)
      real(dp),allocatable::e(:),s(:),m(:);real(dp)::d,z;integer::i;logical::ok
      c=0.0_dp;allocate(e(size(y)),s(size(y)))
      if(spec%model==model_realgarch.and.present(realized))then
         allocate(m(size(y)));call realgarch_filter(y,realized,spec,e,s,m,ok)
         if(.not.ok)return
         do i=1,size(y);d=distribution_pdf(e(i)/s(i),spec%cond_dist,spec%shape,spec%skew,spec%lambda)/s(i)
            z=m(i)/spec%measurement_sd
            c(i)=log(max(d,tiny(1.0_dp)))-log(spec%measurement_sd) &
               -0.5_dp*log(2.0_dp*acos(-1.0_dp))-0.5_dp*z*z
         end do
      else
         call garch_filter(y,spec,e,s,ok);if(.not.ok)return
         do i=1,size(y)
            d=distribution_pdf(e(i)/s(i),spec%cond_dist,spec%shape, &
               spec%skew,spec%lambda)/s(i)
            c(i)=log(max(d,tiny(1.0_dp)))
         end do
      end if
   end subroutine loglik_contributions

   subroutine score_long_run_covariance(scores,lags,cov)
      real(dp),intent(in)::scores(:,:);integer,intent(in)::lags;real(dp),allocatable,intent(out)::cov(:,:)
      real(dp)::weight;integer::n,k,l
      n=size(scores,1);k=size(scores,2);allocate(cov(k,k));cov=matmul(transpose(scores),scores)
      do l=1,min(lags,n-1);weight=1.0_dp-real(l,dp)/real(lags+1,dp);cov=cov+weight*( &
         matmul(transpose(scores(l+1:n,:)),scores(1:n-l,:))+matmul(transpose(scores(1:n-l,:)),scores(l+1:n,:)));end do
   end subroutine score_long_run_covariance

   function sample_standardized(z,method) result(value)
      real(dp),intent(in)::z(:);integer,intent(in)::method;real(dp)::value,u,h,qlo,qhi,scale,excess
      integer::i,n;n=size(z);call random_number(u);i=min(n,max(1,1+int(u*real(n,dp))))
      select case(method)
      case(sampling_kernel);h=1.06_dp*sample_sd(z)*real(n,dp)**(-0.2_dp);value=z(i)+h*random_normal()
      case(sampling_spd)
         qlo=empirical_quantile(z,0.10_dp);qhi=empirical_quantile(z,0.90_dp);call random_number(u)
         if(u<0.10_dp)then
            scale=max(mean_tail(qlo-z,z<qlo),0.05_dp)
            call random_number(u)
            excess=-scale*log(max(u,1.0e-12_dp))
            value=qlo-excess
         else if(u>0.90_dp)then
            scale=max(mean_tail(z-qhi,z>qhi),0.05_dp)
            call random_number(u)
            excess=-scale*log(max(u,1.0e-12_dp))
            value=qhi+excess
         else;call random_number(u);i=min(n,max(1,1+int(u*real(n,dp))));value=z(i);end if
      case default;value=z(i)
      end select
   end function sample_standardized

   function empirical_quantile(x,p) result(value)
      real(dp),intent(in)::x(:),p;real(dp)::value,pos,w;real(dp),allocatable::s(:);integer::i,n
      n=size(x)
      allocate(s(n))
      s=x
      call sort_real(s)
      pos=1.0_dp+max(0.0_dp,min(1.0_dp,p))*real(n-1,dp)
      i=min(n-1,max(1,int(pos)))
      w=pos-real(i,dp)
      if(n==1)then;value=s(1);else;value=(1.0_dp-w)*s(i)+w*s(i+1);end if
   end function empirical_quantile

   subroutine sort_real(x)
      real(dp),intent(inout)::x(:)
      real(dp)::key
      integer::i,j
      do i=2,size(x)
         key=x(i)
         j=i-1
         do while(j>=1)
            if(x(j)<=key)exit
            x(j+1)=x(j)
            j=j-1
         end do
         x(j+1)=key
      end do
   end subroutine sort_real

   pure function sample_sd(x) result(value)
      real(dp),intent(in)::x(:)
      real(dp)::value,m
      if(size(x)<2)then
         value=1.0_dp
      else
         m=sum(x)/real(size(x),dp)
         value=sqrt(max(sum((x-m)**2)/real(size(x)-1,dp),1.0e-20_dp))
      end if
   end function sample_sd

   pure function safe_sd(x) result(value)
      real(dp),intent(in)::x(:)
      real(dp)::value
      value=max(sample_sd(x),1.0e-12_dp)
   end function safe_sd

   pure function sum_safe(x) result(value)
      real(dp),allocatable,intent(in)::x(:)
      real(dp)::value
      if(allocated(x))then
         value=sum(x)
      else
         value=0.0_dp
      end if
   end function sum_safe

   pure logical function finite_loglik(x)
      real(dp),intent(in)::x
      finite_loglik=ieee_is_finite(x).and.abs(x)<huge(1.0_dp)/10.0_dp
   end function finite_loglik

   pure function mean_tail(x,mask) result(value)
      real(dp),intent(in)::x(:)
      logical,intent(in)::mask(:)
      real(dp)::value
      if(count(mask)>0)then
         value=sum(x,mask)/real(count(mask),dp)
      else
         value=0.1_dp
      end if
   end function mean_tail

   pure function median_difference(x) result(value)
      real(dp),intent(in)::x(:)
      real(dp)::value
      real(dp),allocatable::d(:)
      if(size(x)<2)then
         value=1.0_dp
         return
      end if
      allocate(d(size(x)-1))
      d=x(2:)-x(:size(x)-1)
      call sort_real_pure(d)
      value=d((size(d)+1)/2)
   end function median_difference

   pure subroutine sort_real_pure(x)
      real(dp),intent(inout)::x(:)
      real(dp)::key
      integer::i,j
      do i=2,size(x)
         key=x(i)
         j=i-1
         do while(j>=1)
            if(x(j)<=key)exit
            x(j+1)=x(j)
            j=j-1
         end do
         x(j+1)=key
      end do
   end subroutine sort_real_pure

   pure function generate_linear(t0,step,n) result(x)
      real(dp),intent(in)::t0,step
      integer,intent(in)::n
      real(dp)::x(n)
      integer::i
      do i=1,n
         x(i)=t0+real(i,dp)*step
      end do
   end function generate_linear

   pure logical function has_shape_kind(kind)
      integer,intent(in)::kind
      has_shape_kind=kind==dist_std.or.kind==dist_sstd.or.kind==dist_ged.or. &
         kind==dist_sged.or.kind==dist_jsu.or.kind==dist_nig.or. &
         kind==dist_ghyp.or.kind==dist_ghst
   end function has_shape_kind

   pure logical function has_skew_kind(kind)
      integer,intent(in)::kind
      has_skew_kind=kind==dist_snorm.or.kind==dist_sstd.or.kind==dist_sged.or. &
         kind==dist_jsu.or.kind==dist_nig.or.kind==dist_ghyp.or.kind==dist_ghst
   end function has_skew_kind

end module rugarch_complete
