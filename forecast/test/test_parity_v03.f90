program test_parity_v03
   use forecast, only : dp, arima_model, forecast_result, ets_model, bats_model, &
      arima_fit, arima_refit, arima_forecast, auto_arima, arima_simulate, arima_diffuse_loglik, &
      arima_innovations_loglik, &
      ets_forecast, ets_forecast_simulated, &
      ETS_ADD, ETS_MULT, tapered_correlation_ci, bats_auto, bats_refit, tbats_fit, tbats_refit, quantile_type8
   use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan, ieee_is_finite
   implicit none
   real(dp) :: y(40), ymiss(40), xr(40,1), fx(3,1), ys(30), q, yrw(5), zmiss(4), ll, s2, llref
   real(dp),allocatable :: est(:),lo(:),hi(:),sim(:),dres(:),vf(:)
   integer,allocatable :: seed(:)
   integer :: i,nseed,info
   type(arima_model) :: am,aa,rw,ar
   type(forecast_result) :: fc,fs
   type(ets_model) :: em
   type(bats_model) :: bm,br,tm,tr

   do i=1,40
      xr(i,1)=sin(0.37_dp*real(i,dp))+0.2_dp*cos(0.11_dp*real(i,dp))
      y(i)=2.0_dp+3.0_dp*xr(i,1)+0.05_dp*sin(1.7_dp*real(i,dp))
   end do
   am=arima_fit(y,0,0,0,include_mean=.true.,xreg=xr)
   call check(size(am%xreg_coef)==1,'ARIMAX coefficient storage')
   call check(abs(am%xreg_coef(1)-3.0_dp)<0.03_dp,'ARIMAX slope')
   call check(abs(am%intercept-2.0_dp)<0.03_dp,'ARIMAX intercept')
   ar=arima_refit(y+0.1_dp,am,xr)
   call check(maxval(abs(ar%xreg_coef-am%xreg_coef))<1.0e-14_dp,'ARIMA refit fixed coefficients')
   call check(size(ar%residuals)==size(y).and.size(ar%xreg,1)==size(y),'ARIMA refit dimensions')
   fx(:,1)=[-0.5_dp,0.25_dp,1.0_dp]
   fc=arima_forecast(am,y,3,future_xreg=fx)
   call check(maxval(abs(fc%mean-(am%intercept+am%xreg_coef(1)*fx(:,1))))<1.0e-10_dp,'ARIMAX forecast regressors')

   aa=auto_arima(y,m=1,max_p=1,max_q=1,stepwise=.true.,xreg=xr,ic='bic',nmodels=20,d_fixed=0,sd_fixed=0)
   call check(allocated(aa%xreg_coef).and.size(aa%xreg_coef)==1,'auto_arima xreg propagation')
   call check(abs(aa%xreg_coef(1)-3.0_dp)<0.08_dp,'auto_arima ARIMAX estimate')
   aa=auto_arima(y,m=1,max_p=1,max_q=1,stepwise=.true.,approximation=.true.,xreg=xr, &
      nmodels=8,d_fixed=0,sd_fixed=0,truncate=20)
   call check(size(aa%xreg,1)==size(y),'auto_arima truncate full-sample refit')


   rw%d=1
   rw%sd=0
   rw%m=1
   rw%include_drift=.true.
   rw%drift=1.0_dp
   rw%sigma2=0.0_dp
   rw%ar=[real(dp)::]
   rw%ma=[real(dp)::]
   rw%sar=[real(dp)::]
   rw%sma=[real(dp)::]
   sim=arima_simulate(rw,4,burnin=0,initial=[10.0_dp])
   call check(maxval(abs(sim-[11.0_dp,12.0_dp,13.0_dp,14.0_dp]))<1.0e-12_dp,'integrated ARIMA simulation')

   yrw=[3.0_dp,4.0_dp,2.0_dp,5.0_dp,6.0_dp]
   call arima_diffuse_loglik(yrw,[real(dp)::],[real(dp)::],1,1,0,0.0_dp,ll,s2,dres,vf,info)
   llref=-0.5_dp*4.0_dp*(log(2.0_dp*acos(-1.0_dp)*3.75_dp)+1.0_dp)
   call check(info==0,'diffuse random-walk likelihood status')
   call check(abs(s2-3.75_dp)<1.0e-8_dp,'diffuse random-walk sigma2')
   call check(abs(ll-llref)<1.0e-7_dp,'diffuse random-walk likelihood')
   zmiss=[1.0_dp,ieee_value(0.0_dp,ieee_quiet_nan),-1.0_dp,2.0_dp]
   call arima_innovations_loglik(zmiss,[real(dp)::],[real(dp)::],0.0_dp,ll,s2,dres,vf,info)
   llref=-0.5_dp*3.0_dp*(log(2.0_dp*acos(-1.0_dp)*2.0_dp)+1.0_dp)
   call check(info==0.and.ieee_is_finite(ll),'stationary ARIMA missing-value likelihood status')
   call check(abs(s2-2.0_dp)<1.0e-10_dp.and.abs(ll-llref)<1.0e-9_dp,'stationary ARIMA missing-value likelihood')

   ymiss=y
   ymiss(17)=ieee_value(0.0_dp,ieee_quiet_nan)
   ar=arima_fit(ymiss,0,0,0,include_mean=.true.,method='ml')
   call check(ieee_is_finite(ar%loglik).and.ieee_is_finite(ar%intercept),'ARIMA ML fit with missing observation')


   em%error_type=ETS_MULT
   em%trend_type=ETS_ADD
   em%season_type=ETS_MULT
   em%m=3
   em%damped=.true.
   em%alpha=0.2_dp
   em%beta=0.05_dp
   em%gamma=0.1_dp
   em%phi=0.9_dp
   em%sigma2=0.04_dp
   em%state=[10.0_dp,0.5_dp,0.8_dp,1.0_dp,1.2_dp]
   em%residuals=[0.1_dp,-0.1_dp,0.05_dp,-0.05_dp]
   fc=ets_forecast(em,4)
   call check(maxval(abs(fc%mean-[12.6_dp,10.95_dp,9.084_dp,14.0803092_dp]))<1.0e-8_dp,'ETS class-3 mean')
   call check(maxval(abs(fc%se**2-[6.3504_dp,5.08275_dp,3.7561827744_dp,10.36785571289279_dp]))<2.0e-8_dp, &
      'ETS class-3 variance')
   call random_seed(size=nseed)
   allocate(seed(nseed))
   seed=[(1777+31*i,i=1,nseed)]
   call random_seed(put=seed)
   fs=ets_forecast_simulated(em,3,npaths=300)
   call check(all(fs%upper>fs%lower).and.all(fs%se>0.0_dp),'ETS simulated intervals')

   call random_seed(put=seed)
   call tapered_correlation_ci(y,6,level=90.0_dp,nsim=30,estimate=est,lower=lo,upper=hi)
   call check(size(est)==6.and.all(lo<=hi),'tapered ACF bootstrap intervals')
   call tapered_correlation_ci(y,5,partial=.true.,level=90.0_dp,nsim=25,estimate=est,lower=lo,upper=hi)
   call check(size(est)==5.and.all(lo<=hi),'tapered PACF bootstrap intervals')

   q=quantile_type8([1.0_dp,2.0_dp,3.0_dp,4.0_dp],0.25_dp)
   call check(abs(q-1.4166666666666667_dp)<1.0e-12_dp,'R type-8 quantile')

   do i=1,30
      ys(i)=20.0_dp+0.15_dp*real(i,dp)+1.5_dp*sin(2.0_dp*acos(-1.0_dp)*real(i,dp)/6.0_dp)
   end do
   bm=bats_auto(ys,[6],use_arma=.false.,consider_boxcox=.true.,consider_trend=.true.,consider_damped=.true., &
      boxcox_choice=-1,trend_choice=1,damped_choice=-1)
   call check(bm%has_trend.and..not.bm%use_boxcox.and.abs(bm%phi-1.0_dp)<1.0e-12_dp,'BATS forced controls')
   br=bats_refit(ys+0.2_dp,bm)
   call check(size(br%residuals)==size(ys).and.size(br%y)==size(ys),'BATS structural refit lengths')
   call check(abs(br%alpha-bm%alpha)<1.0e-14_dp.and.abs(br%beta-bm%beta)<1.0e-14_dp,'BATS refit fixed parameters')
   tm=tbats_fit(ys,[6],[1],use_trend=.true.,damped=.false.,optimize=.false.)
   tr=tbats_refit(ys+0.1_dp,tm)
   call check(tr%trigonometric.and.size(tr%residuals)==size(ys),'TBATS structural refit')

   print '(a)','test_parity_v03: PASS'
contains
   subroutine check(ok,msg)
      logical,intent(in)::ok
      character(*),intent(in)::msg
      if(.not.ok)then
         write(*,'(a)')'FAIL: '//trim(msg)
         error stop 1
      end if
   end subroutine check
end program test_parity_v03
