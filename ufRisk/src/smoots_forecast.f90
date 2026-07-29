! SPDX-License-Identifier: GPL-3.0-only
module smoots_forecast
   use smoots_kinds, only : dp
   use smoots_status, only : sm_ok, sm_invalid_input, sm_iteration_limit
   use smoots_types
   use smoots_stats, only : normal_quantile, empirical_quantile, sample_with_replacement, seed_rng, mean_value
   use smoots_arma
   use smoots_estimation, only : msmooth
   implicit none
   private
   public :: trend_forecast, normal_forecast, bootstrap_forecast
   public :: model_forecast, rolling_backtest
contains
   subroutine trend_forecast(model,h,mode,forecast,status)
      type(smooth_result),intent(in)::model
      integer,intent(in)::h,mode
      real(dp),allocatable,intent(out)::forecast(:)
      integer,intent(out)::status
      real(dp)::slope
      integer::k,n
      allocate(forecast(max(0,h)));forecast=0.0_dp
      if(h<1.or..not.allocated(model%estimate).or.size(model%estimate)<2)then;status=sm_invalid_input;return;end if
      n=size(model%estimate);slope=0.0_dp
      if(mode==sm_trend_linear)slope=model%estimate(n)-model%estimate(n-1)
      do k=1,h;forecast(k)=model%estimate(n)+real(k,dp)*slope;end do
      status=sm_ok
   end subroutine trend_forecast

   subroutine normal_forecast(x,p,q,include_mean,h,confidence,result)
      real(dp),intent(in)::x(:),confidence
      integer,intent(in)::p,q,h
      logical,intent(in)::include_mean
      type(forecast_result),intent(out)::result
      real(dp),allocatable::psi(:)
      real(dp)::z
      integer::status,k
      result%horizon=h;result%confidence_level=confidence
      if(h<1.or.confidence<=0.0_dp.or.confidence>=1.0_dp)then;result%status=sm_invalid_input;return;end if
      call fit_arma(x,p,q,include_mean,result%model,status)
      if(status/=sm_ok.and.status/=sm_iteration_limit)then;result%status=status;return;end if
      call arma_point_forecast(x,result%model%residuals,result%model%ar,result%model%ma,result%model%mean_value,h,result%point,status)
      call ma_infinity(result%model%ar,result%model%ma,h-1,psi,status)
      allocate(result%lower(h),result%upper(h));z=normal_quantile(1.0_dp-(1.0_dp-confidence)/2.0_dp)
      do k=1,h
         result%lower(k)=result%point(k)-z*sqrt(result%model%sigma2*sum(psi(1:k)**2))
         result%upper(k)=result%point(k)+z*sqrt(result%model%sigma2*sum(psi(1:k)**2))
      end do
      result%status=sm_ok
   end subroutine normal_forecast

   subroutine bootstrap_forecast(x,p,q,include_mean,h,simulations,burn,confidence,result,seed)
      real(dp),intent(in)::x(:),confidence
      integer,intent(in)::p,q,h,simulations,burn
      logical,intent(in)::include_mean
      type(forecast_result),intent(out)::result
      integer(kind=8),intent(in),optional::seed
      type(arma_model)::main_model,star_model
      real(dp),allocatable::fi(:),eps(:),xstar(:),proxy_res(:),proxy_fit(:),star_point(:),true_point(:),future_eps(:)
      real(dp)::plo,phi
      integer::status,s,k,n
      n=size(x);result%horizon=h;result%simulations=simulations;result%confidence_level=confidence
      if(h<1.or.simulations<1.or.burn<0.or.confidence<=0.0_dp.or.confidence>=1.0_dp)then;result%status=sm_invalid_input;return;end if
      if(present(seed))call seed_rng(seed)
      call fit_arma(x,p,q,include_mean,main_model,status)
      if(status/=sm_ok.and.status/=sm_iteration_limit)then;result%status=status;return;end if
      result%model=main_model
      call arma_point_forecast(x,main_model%residuals,main_model%ar,main_model%ma,main_model%mean_value,h,result%point,status)
      allocate(fi(n));fi=main_model%residuals-mean_value(main_model%residuals)
      allocate(result%errors(simulations,h),eps(burn+n+h),future_eps(h),proxy_res(n),proxy_fit(n))
      do s=1,simulations
         call sample_with_replacement(fi,eps)
         call simulate_arma(main_model%ar,main_model%ma,main_model%mean_value,n,burn,xstar,status, &
              innovations=eps(burn+1:burn+n),start_innovations=eps(1:burn))
         call fit_arma(xstar,p,q,include_mean,star_model,status,max_iterations=100)
         if(status/=sm_ok.and.status/=sm_iteration_limit)then
            star_model=main_model
         end if
         call arma_residuals(x,star_model%ar,star_model%ma,star_model%mean_value,proxy_res,proxy_fit)
         call arma_point_forecast(x,proxy_res,star_model%ar,star_model%ma,star_model%mean_value,h,star_point,status)
         future_eps=eps(burn+n+1:burn+n+h)
         call conditional_true_forecast(x,main_model%residuals,future_eps,main_model%ar,main_model%ma,main_model%mean_value,true_point)
         result%errors(s,:)=true_point-star_point
      end do
      allocate(result%lower(h),result%upper(h));plo=(1.0_dp-confidence)/2.0_dp;phi=1.0_dp-plo
      do k=1,h
         result%lower(k)=result%point(k)+empirical_quantile(result%errors(:,k),plo)
         result%upper(k)=result%point(k)+empirical_quantile(result%errors(:,k),phi)
      end do
      result%status=sm_ok
   end subroutine bootstrap_forecast

   subroutine conditional_true_forecast(x,residuals,future_eps,ar,ma,mean_x,forecast)
      real(dp),intent(in)::x(:),residuals(:),future_eps(:),ar(:),ma(:),mean_x
      real(dp),allocatable,intent(out)::forecast(:)
      real(dp),allocatable::xx(:),ee(:)
      integer::n,h,p,q,t,j
      n=size(x);h=size(future_eps);p=size(ar);q=size(ma)
      allocate(xx(n+h),ee(n+h),forecast(h));xx(1:n)=x-mean_x;ee(1:n)=residuals;xx(n+1:n+h)=0.0_dp;ee(n+1:n+h)=future_eps
      do t=n+1,n+h
         xx(t)=ee(t)
         do j=1,p;xx(t)=xx(t)+ar(j)*xx(t-j);end do
         do j=1,q;xx(t)=xx(t)+ma(j)*ee(t-j);end do
      end do
      forecast=xx(n+1:n+h)+mean_x
   end subroutine conditional_true_forecast

   subroutine model_forecast(smooth_model,p,q,h,use_bootstrap,confidence,result,simulations,burn,trend_mode,seed)
      type(smooth_result),intent(in)::smooth_model
      integer,intent(in)::p,q,h
      logical,intent(in)::use_bootstrap
      real(dp),intent(in)::confidence
      type(forecast_result),intent(out)::result
      integer,intent(in),optional::simulations,burn,trend_mode
      integer(kind=8),intent(in),optional::seed
      real(dp),allocatable::tf(:)
      integer::sims,bb,tm,status
      sims=10000;if(present(simulations))sims=simulations;bb=1000;if(present(burn))bb=burn
      tm=sm_trend_linear;if(present(trend_mode))tm=trend_mode
      call trend_forecast(smooth_model,h,tm,tf,status)
      if(use_bootstrap)then
         if(present(seed))then
            call bootstrap_forecast(smooth_model%residuals,p,q,.false.,h,sims,bb,confidence,result,seed)
         else
            call bootstrap_forecast(smooth_model%residuals,p,q,.false.,h,sims,bb,confidence,result)
         end if
      else
         call normal_forecast(smooth_model%residuals,p,q,.false.,h,confidence,result)
      end if
      if(result%status==sm_ok)then;result%point=result%point+tf;result%lower=result%lower+tf;result%upper=result%upper+tf;end if
   end subroutine model_forecast

   subroutine rolling_backtest(y,k,p,q,use_bootstrap,confidence,result,simulations,burn,trend_mode,seed)
      real(dp),intent(in)::y(:),confidence
      integer,intent(in)::k,p,q
      logical,intent(in)::use_bootstrap
      type(rolling_result),intent(out)::result
      integer,intent(in),optional::simulations,burn,trend_mode
      integer(kind=8),intent(in),optional::seed
      real(dp),allocatable::yin(:),yout(:),trendfc(:),xwork(:),uwork(:),residual_fc(:)
      type(forecast_result)::one
      real(dp)::qlow,qhigh,den1,den2
      integer::n,nin,i,j,status,sims,bb,tm,l
      n=size(y);nin=n-k;result%horizon=k
      if(k<1.or.nin<10)then;result%status=sm_invalid_input;return;end if
      sims=10000;if(present(simulations))sims=simulations;bb=1000;if(present(burn))bb=burn
      tm=sm_trend_linear;if(present(trend_mode))tm=trend_mode
      allocate(yin(nin),yout(k));yin=y(1:nin);yout=y(nin+1:n)
      call msmooth(yin,result%smooth_model)
      call trend_forecast(result%smooth_model,k,tm,trendfc,status)
      call fit_arma(result%smooth_model%residuals,p,q,.false.,result%arma,status)
      l=max(p,q);allocate(xwork(l+k),uwork(l+k),residual_fc(k));xwork=0.0_dp;uwork=0.0_dp
      if(l>0)then;xwork(1:l)=result%smooth_model%residuals(nin-l+1:nin);uwork(1:l)=result%arma%residuals(nin-l+1:nin);end if
      do i=1,k
         residual_fc(i)=0.0_dp
         do j=1,p;residual_fc(i)=residual_fc(i)+result%arma%ar(j)*xwork(l+i-j);end do
         do j=1,q;residual_fc(i)=residual_fc(i)+result%arma%ma(j)*uwork(l+i-j);end do
         xwork(l+i)=yout(i)-trendfc(i);uwork(l+i)=xwork(l+i)-residual_fc(i)
      end do
      if(use_bootstrap)then
         if(present(seed))then;call bootstrap_forecast(result%smooth_model%residuals,p,q,.false.,1,sims,bb,confidence,one,seed)
         else;call bootstrap_forecast(result%smooth_model%residuals,p,q,.false.,1,sims,bb,confidence,one);end if
         qlow=empirical_quantile(one%errors(:,1),(1.0_dp-confidence)/2.0_dp)
         qhigh=empirical_quantile(one%errors(:,1),1.0_dp-(1.0_dp-confidence)/2.0_dp)
      else
         qlow=normal_quantile((1.0_dp-confidence)/2.0_dp)*sqrt(result%arma%sigma2)
         qhigh=normal_quantile(1.0_dp-(1.0_dp-confidence)/2.0_dp)*sqrt(result%arma%sigma2)
      end if
      allocate(result%observed(k),result%trend_forecast(k),result%residual_forecast(k),result%point(k),result%lower(k),result%upper(k),result%breach(k),result%breach_value(k))
      result%observed=yout;result%trend_forecast=trendfc;result%residual_forecast=residual_fc
      result%point=trendfc+residual_fc;result%lower=result%point+qlow;result%upper=result%point+qhigh
      result%breach=(yout<result%lower).or.(yout>result%upper);result%breach_value=0.0_dp
      where(yout<result%lower)result%breach_value=yout-result%lower
      where(yout>result%upper)result%breach_value=yout-result%upper
      den1=sum(abs(yin(2:nin)-yin(1:nin-1)))/real(nin-1,dp)
      den2=sum((yin(2:nin)-yin(1:nin-1))**2)/real(nin-1,dp)
      result%mase=sum(abs(yout-result%point))/real(k,dp)/max(den1,epsilon(1.0_dp))
      result%rmsse=sqrt((sum((yout-result%point)**2)/real(k,dp))/max(den2,epsilon(1.0_dp)))
      result%status=sm_ok
   end subroutine rolling_backtest
end module smoots_forecast
