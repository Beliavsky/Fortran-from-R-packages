! Part of the experimental modern Fortran translation of rugarch 1.5-6.
! This file was created or modified for the Fortran project on 2026-07-23.
! Original rugarch authors retain copyright; see NOTICE and ORIGIN.md.
! SPDX-License-Identifier: GPL-3.0-only

program test_complete_features
   use rugarch
   implicit none
   integer, parameter :: n=180
   type(garch_spec) :: spec
   type(garch_fit_result) :: fit
   type(extended_garch_fit_result) :: ext
   type(bootstrap_forecast_result) :: boot
   type(arfima_spec) :: aspec
   type(arfima_fit_result) :: afit
   type(arfima_rolling_result) :: aroll
   type(arfima_distribution_result) :: adist
   type(arfima_bootstrap_result) :: aboot
   type(arfima_multi_fit_result) :: amfit
   type(arfima_multi_forecast_result) :: amfor
   type(arfima_cv_result) :: cv
   type(forecast_performance_result) :: perf
   type(distribution_moment_result) :: moments
   type(gh_parameter_result) :: gh
   type(numeric_index_result) :: idx
   real(dp) :: y(n),sigma(n),eps(n),xreg(n,1),vreg(n,1),data2(n,2)
   real(dp) :: ax(n)
   integer :: i
   integer :: train_end(2),test_end(2)

   call seed_rng(24680)
   spec=make_garch_spec(1,1,model_sgarch,dist_norm)
   spec%omega=0.04_dp
   spec%alpha(1)=0.08_dp
   spec%beta(1)=0.87_dp
   call simulate_garch(spec,n,y,sigma,eps,burn_in=200)
   do i=1,n
      xreg(i,1)=sin(0.03_dp*real(i,dp))
      vreg(i,1)=0.01_dp+0.01_dp*cos(0.02_dp*real(i,dp))**2
   end do
   y=y+0.25_dp+0.40_dp*xreg(:,1)

   ext=fit_garch_extended(y,model_sgarch,1,1,dist_norm, &
      mean_regressors=xreg,variance_regressors=vreg,variance_targeting=.true., &
      attach_covariance=.false.,max_iterations=350)
   call assert_true(ext%fit%status<=1,'extended GARCH fit')
   call assert_true(size(ext%mean_beta)==1,'mean regressor coefficient')
   call assert_true(size(ext%variance_beta)==1,'variance regressor coefficient')
   call assert_true(ext%variance_targeting,'variance targeting flag')

   fit=fit_model(y,model_sgarch,1,1,dist_norm,max_iterations=350)
   call assert_true(fit%status<=1,'base fit for covariance')
   call attach_garch_covariance(y,fit,ext,newey_west_lags=2)
   call assert_true(ext%covariance_status==0,'automatic covariance')
   call assert_true(all(ext%standard_error>=0.0_dp),'standard errors')

   boot=garch_bootstrap_forecast(fit,3,12,sampling_kernel,bootstrap_partial,confidence=0.90_dp)
   call assert_true(size(boot%return_path,1)==3,'GARCH bootstrap horizon')
   call assert_true(all(boot%sigma_path>0.0_dp),'GARCH bootstrap volatility')
   boot=garch_bootstrap_forecast(fit,2,8,sampling_spd,bootstrap_partial)
   call assert_true(all(boot%upper>=boot%lower),'SPD intervals')
   boot=garch_bootstrap_forecast(fit,1,3,sampling_raw,bootstrap_full,n_bootfit=2,max_iterations=100)
   call assert_true(size(boot%return_path,2)==3,'full bootstrap forecast')

   aspec=make_arfima_spec(1,0)
   aspec%mean=0.1_dp
   aspec%d=0.20_dp
   aspec%ar(1)=0.25_dp
   aspec%innovation_sd=0.8_dp
   call simulate_arfima(aspec,n,ax,burn_in=300)
   afit=fit_arfima(ax,1,0,.true.)
   call assert_true(afit%status==0,'ARFIMA fit')
   aroll=rolling_arfima_forecast(ax,1,0,20,refit_every=5,window_size=100,moving_window=.true.)
   call assert_true(size(aroll%mean)==20,'ARFIMA rolling')
   adist=arfima_parametric_distribution(afit,3,nobs=100)
   call assert_true(size(adist%d)==3,'ARFIMA distribution')
   aboot=arfima_bootstrap_forecast(ax,afit,4,10,sampling_raw)
   call assert_true(size(aboot%path,2)==10,'ARFIMA bootstrap')
   data2(:,1)=ax
   data2(:,2)=0.5_dp*ax+0.5_dp*y
   amfit=multifit_arfima(data2,1,0,.true.)
   amfor=multiforecast_arfima(data2,amfit,3)
   call assert_true(size(amfor%mean,2)==2,'ARFIMA multi forecast')
   train_end=[100,130]
   test_end=[110,140]
   cv=arfima_cross_validation(ax,train_end,test_end,1,1,.true.)
   call assert_true(cv%best>=1,'ARFIMA cross validation')

   perf=forecast_performance(ax(2:20),ax(1:19),ax)
   call assert_true(perf%rmse>=0.0_dp,'forecast performance')
   moments=distribution_moments(dist_std,shape=8.0_dp)
   call assert_close(moments%excess_kurtosis,1.5_dp,1.0e-12_dp,'Student t kurtosis')
   gh=ghyp_transform(0.0_dp,1.0_dp,0.1_dp,2.0_dp,-0.5_dp)
   call assert_true(gh%status==0.and.gh%alpha>abs(gh%beta),'GH transform')
   idx=generate_forward_numeric(10.0_dp,3,0.5_dp)
   call assert_close(idx%value(3),11.5_dp,1.0e-12_dp,'forward numeric index')

   print '(a)','complete numerical feature tests passed'
contains
   subroutine assert_true(condition,message)
      logical,intent(in)::condition
      character(len=*),intent(in)::message
      if(.not.condition)then
         print '(a)',trim(message)
         error stop 1
      end if
   end subroutine assert_true

   subroutine assert_close(x,expected,tolerance,message)
      real(dp),intent(in)::x,expected,tolerance
      character(len=*),intent(in)::message
      if(abs(x-expected)>tolerance)then
         print '(a,2es16.7)',trim(message)//': ',x,expected
         error stop 1
      end if
   end subroutine assert_close
end program test_complete_features
