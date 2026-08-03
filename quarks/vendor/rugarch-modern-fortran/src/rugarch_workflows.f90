module rugarch_workflows
   use rugarch_kinds,only:dp
   use rugarch_types,only:garch_spec,garch_fit_result,forecast_result,model_realgarch
   use rugarch_fit,only:fit_model
   use rugarch_models,only:forecast_volatility,simulate_garch,simulate_realgarch,true_persistence
   use rugarch_risk,only:value_at_risk,expected_shortfall
   implicit none
   private

   type,public::multi_fit_result
      type(garch_fit_result),allocatable::fit(:)
      integer,allocatable::status(:)
   end type multi_fit_result

   type,public::multi_forecast_result
      real(dp),allocatable::mean(:,:),sigma(:,:)
      integer,allocatable::status(:)
   end type multi_forecast_result

   type,public::rolling_forecast_result
      real(dp),allocatable::actual(:),mean(:),sigma(:),var(:),es(:)
      integer,allocatable::status(:),fit_end(:)
   end type rolling_forecast_result

   type,public::bootstrap_distribution_result
      real(dp),allocatable::mean(:),omega(:),alpha1(:),beta1(:),persistence(:),log_likelihood(:)
      real(dp),allocatable::shape(:),skew(:),lambda(:)
      integer,allocatable::status(:)
   end type bootstrap_distribution_result

   public::multifit_garch,multiforecast_garch,rolling_garch_forecast
   public::garch_parametric_distribution

contains

   function multifit_garch(data,model,p,q,cond_dist,fit_mean,fit_shape,fit_skew,fit_lambda, &
      max_iterations) result(ans)
      real(dp),intent(in)::data(:,:)
      integer,intent(in)::model,p,q
      integer,intent(in),optional::cond_dist,max_iterations
      logical,intent(in),optional::fit_mean,fit_shape,fit_skew,fit_lambda
      type(multi_fit_result)::ans
      integer::j
      allocate(ans%fit(size(data,2)),ans%status(size(data,2)))
      do j=1,size(data,2)
         ans%fit(j)=fit_model(data(:,j),model,p,q,cond_dist,fit_mean,fit_shape,fit_skew, &
            fit_lambda,max_iterations=max_iterations)
         ans%status(j)=ans%fit(j)%status
      end do
   end function multifit_garch

   function multiforecast_garch(fits,horizon) result(ans)
      type(multi_fit_result),intent(in)::fits
      integer,intent(in)::horizon
      type(multi_forecast_result)::ans
      real(dp),allocatable::vol(:)
      integer::j,m
      m=size(fits%fit);allocate(ans%mean(horizon,m),ans%sigma(horizon,m),ans%status(m),vol(horizon))
      ans%mean=0.0_dp;ans%sigma=0.0_dp
      do j=1,m
         ans%status(j)=fits%fit(j)%status
         if(.not.allocated(fits%fit(j)%residuals) .or. .not.allocated(fits%fit(j)%sigma))then
            ans%status(j)=max(2,ans%status(j));cycle
         end if
         call forecast_volatility(fits%fit(j)%spec,fits%fit(j)%residuals, &
            fits%fit(j)%sigma,horizon,vol)
         ans%sigma(:,j)=vol;ans%mean(:,j)=fits%fit(j)%spec%mean
      end do
   end function multiforecast_garch

   function rolling_garch_forecast(data,model,p,q,forecast_length,refit_every,window_size, &
      moving_window,alpha,cond_dist,fit_mean,fit_shape,fit_skew,fit_lambda,max_iterations,realized) result(ans)
      real(dp),intent(in)::data(:)
      integer,intent(in)::model,p,q,forecast_length
      integer,intent(in),optional::refit_every,window_size,cond_dist,max_iterations
      logical,intent(in),optional::moving_window,fit_mean,fit_shape,fit_skew,fit_lambda
      real(dp),intent(in),optional::alpha,realized(:)
      type(rolling_forecast_result)::ans
      type(garch_fit_result)::fit
      real(dp)::prob
      real(dp),allocatable::vol(:)
      integer::n,nf,first,i,end_fit,start_fit,refit,win
      logical::moving,need_fit

      n=size(data);nf=min(max(1,forecast_length),n-30);first=n-nf+1
      refit=1;if(present(refit_every))refit=max(1,refit_every)
      win=n;if(present(window_size))win=max(30,window_size)
      moving=.false.;if(present(moving_window))moving=moving_window
      prob=0.01_dp;if(present(alpha))prob=alpha
      allocate(ans%actual(nf),ans%mean(nf),ans%sigma(nf),ans%var(nf),ans%es(nf), &
         ans%status(nf),ans%fit_end(nf),vol(1))
      ans%actual=data(first:n);ans%mean=0.0_dp;ans%sigma=0.0_dp;ans%var=0.0_dp;ans%es=0.0_dp
      ans%status=1;ans%fit_end=0;need_fit=.true.
      do i=1,nf
         end_fit=first+i-2
         if(i==1 .or. mod(i-1,refit)==0)need_fit=.true.
         if(need_fit)then
            start_fit=1;if(moving)start_fit=max(1,end_fit-win+1)
            if(model==model_realgarch)then
               if(.not.present(realized))then;ans%status(i:nf)=3;return;end if
               fit=fit_model(data(start_fit:end_fit),model,p,q,cond_dist,fit_mean,fit_shape,fit_skew, &
                  fit_lambda,realized=realized(start_fit:end_fit),max_iterations=max_iterations)
            else
               fit=fit_model(data(start_fit:end_fit),model,p,q,cond_dist,fit_mean,fit_shape,fit_skew, &
                  fit_lambda,max_iterations=max_iterations)
            end if
            need_fit=.false.
         end if
         ans%fit_end(i)=end_fit;ans%status(i)=fit%status
         if(.not.allocated(fit%residuals) .or. .not.allocated(fit%sigma))cycle
         call forecast_volatility(fit%spec,fit%residuals,fit%sigma,1,vol)
         ans%mean(i)=fit%spec%mean;ans%sigma(i)=vol(1)
         ans%var(i)=value_at_risk(prob,ans%mean(i),ans%sigma(i),fit%spec%cond_dist, &
            fit%spec%shape,fit%spec%skew,fit%spec%lambda)
         ans%es(i)=expected_shortfall(prob,ans%mean(i),ans%sigma(i),fit%spec%cond_dist, &
            fit%spec%shape,fit%spec%skew,fit%spec%lambda)
      end do
   end function rolling_garch_forecast

   function garch_parametric_distribution(fit,nboot,nobs,max_iterations) result(ans)
      type(garch_fit_result),intent(in)::fit
      integer,intent(in)::nboot
      integer,intent(in),optional::nobs,max_iterations
      type(bootstrap_distribution_result)::ans
      type(garch_fit_result)::bootfit
      real(dp),allocatable::y(:),sigma(:),eps(:),realized(:)
      integer::b,n,p,q

      n=1000;if(present(nobs))n=max(50,nobs)
      p=size(fit%spec%alpha);q=size(fit%spec%beta)
      allocate(ans%mean(nboot),ans%omega(nboot),ans%alpha1(nboot),ans%beta1(nboot), &
         ans%persistence(nboot),ans%log_likelihood(nboot),ans%shape(nboot),ans%skew(nboot), &
         ans%lambda(nboot),ans%status(nboot),y(n),sigma(n),eps(n))
      ans%mean=0.0_dp;ans%omega=0.0_dp;ans%alpha1=0.0_dp;ans%beta1=0.0_dp
      ans%persistence=0.0_dp;ans%log_likelihood=-huge(1.0_dp);ans%shape=0.0_dp
      ans%skew=0.0_dp;ans%lambda=0.0_dp;ans%status=1
      if(fit%spec%model==model_realgarch)allocate(realized(n))
      do b=1,nboot
         if(fit%spec%model==model_realgarch)then
            call simulate_realgarch(fit%spec,n,y,realized,sigma,eps,burn_in=200)
            bootfit=fit_model(y,fit%spec%model,p,q,fit%spec%cond_dist,.true.,.true.,.true.,.true., &
               realized=realized,max_iterations=max_iterations)
         else
            call simulate_garch(fit%spec,n,y,sigma,eps,burn_in=200)
            bootfit=fit_model(y,fit%spec%model,p,q,fit%spec%cond_dist,.true.,.true.,.true.,.true., &
               max_iterations=max_iterations,fgarch_submodel=fit%spec%fgarch_submodel, &
               figarch_truncation=fit%spec%figarch_truncation)
         end if
         ans%status(b)=bootfit%status
         ans%mean(b)=bootfit%spec%mean;ans%omega(b)=bootfit%spec%omega
         if(size(bootfit%spec%alpha)>0)ans%alpha1(b)=bootfit%spec%alpha(1)
         if(size(bootfit%spec%beta)>0)ans%beta1(b)=bootfit%spec%beta(1)
         ans%persistence(b)=true_persistence(bootfit%spec)
         ans%log_likelihood(b)=bootfit%log_likelihood;ans%shape(b)=bootfit%spec%shape
         ans%skew(b)=bootfit%spec%skew;ans%lambda(b)=bootfit%spec%lambda
      end do
   end function garch_parametric_distribution

end module rugarch_workflows
