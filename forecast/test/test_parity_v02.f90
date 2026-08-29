program test_parity_v02
   use forecast, only : dp, forecast_result, cvar_result, arima_impulse_weights, ocsb_critical_value, &
      bagged_forecast, cvar, stlf_forecast, mean_forecast, seasonal_strength, arima_innovations_loglik, tbats_auto, bats_model, &
      month_days_sequence, easter_gregorian, easter_effect, business_days_sequence, fourier_terms_multi, tapered_acf, &
      ets_model, ses_fit, ets_forecast
   implicit none
   real(dp)::y(36),z(12),crit,ll,sig
   real(dp),allocatable::w(:),boots(:,:),med(:),res(:),vf(:),eff(:),fm(:,:),tacf(:)
   integer,allocatable::idays(:)
   type(forecast_result)::fc
   type(cvar_result)::cv
   type(bats_model)::bm
   type(ets_model)::emod
   integer::info,em,ed
   integer::i

   call arima_impulse_weights([0.5_dp],[0.2_dp],0,1,0,4,w)
   call check(maxval(abs(w-[1.0_dp,0.7_dp,0.35_dp,0.175_dp]))<1.0e-12_dp,'ARMA impulse weights')
   call arima_impulse_weights([real(dp)::],[real(dp)::],1,1,0,5,w)
   call check(maxval(abs(w-1.0_dp))<1.0e-12_dp,'integrated impulse weights')

   call arima_innovations_loglik([0.2_dp,-0.1_dp,0.4_dp,0.3_dp,-0.2_dp,0.5_dp], &
      [0.4_dp],[-0.25_dp],0.1_dp,ll,sig,res,vf,info)
   call check(info==0 .and. abs(ll+1.0324619779361517_dp)<2.0e-9_dp,'exact ARMA Gaussian likelihood')
   call check(abs(sig-0.08221494180551549_dp)<2.0e-9_dp,'profiled ARMA variance')

   crit=ocsb_critical_value(12)
   call check(abs(crit+1.8029627911815704_dp)<1.0e-12_dp,'OCSB critical value')

   emod=ses_fit([(10.0_dp+0.3_dp*real(i,dp),i=1,30)],alpha=0.35_dp)
   call check(abs(emod%alpha-0.35_dp)<1.0e-14_dp,'fixed ETS alpha')
   fc=ets_forecast(emod,4)
   call check(abs(fc%se(2)**2-emod%sigma2*(1.0_dp+0.35_dp**2))<1.0e-9_dp, &
      'ETS class-1 variance h2')
   call check(abs(fc%se(4)**2-emod%sigma2*(1.0_dp+3.0_dp*0.35_dp**2))<1.0e-9_dp, &
      'ETS class-1 variance h4')

   z=[(real(i,dp),i=1,12)]
   fc=bagged_forecast(z,3,mean_cb,num_boot=1,period=1,forecasts_boot=boots,median_forecast=med)
   call check(maxval(abs(fc%mean-6.5_dp))<1.0e-12_dp,'generic bagged mean')
   call check(maxval(abs(fc%lower(:,1)-fc%upper(:,1)))<1.0e-12_dp,'bagged ensemble range')
   call check(maxval(abs(med-fc%mean))<1.0e-12_dp,'bagged median')
   call check(size(boots,2)==1,'bagged component forecasts')

   cv=cvar(z,subset_mean_cb,k=4,blocked=.true.,lb_lags=3)
   call check(cv%k==4 .and. size(cv%fold_accuracy,1)==4,'CVar dimensions')
   call check(abs(cv%testfit(1)-8.0_dp)<1.0e-12_dp,'CVar blocked first fold')
   call check(all(cv%testfit==cv%testfit) .and. cv%lb_p_value>=0.0_dp .and. cv%lb_p_value<=1.0_dp,'CVar diagnostics')

   do i=1,36
      y(i)=20.0_dp+0.05_dp*real(i,dp)+2.0_dp*sin(2.0_dp*acos(-1.0_dp)*real(i,dp)/12.0_dp)
   end do
   call check(seasonal_strength(y,12)>0.7_dp,'MSTL seasonal strength')
   fc=stlf_forecast(y,[12],12,mean_cb)
   call check(size(fc%mean)==12 .and. all(fc%mean==fc%mean),'STLF forecast')
   call check(maxval(fc%mean)-minval(fc%mean)>2.0_dp,'STLF reseasonalization')

   bm=tbats_auto(y,[12],use_arma=.false.,consider_boxcox=.false.,consider_trend=.false.,consider_damped=.false.)
   call check(bm%aic<huge(1.0_dp)/10.0_dp .and. allocated(bm%k),'TBATS automatic harmonic search')

   idays=month_days_sequence(2024,1,12,3)
   call check(all(idays==[31,29,31]),'monthdays leap year')
   idays=month_days_sequence(2024,1,4,1)
   call check(idays(1)==91,'quarter days leap year')
   call easter_gregorian(2024,em,ed)
   call check(em==3 .and. ed==31,'Gregorian Easter date')
   eff=easter_effect(2024,3,12,2,.true.)
   call check(maxval(abs(eff-[0.75_dp,0.25_dp]))<1.0e-12_dp,'Easter monthly fractions')
   idays=business_days_sequence(2024,1,12,1,reshape([2024,1,1],[1,3]))
   call check(idays(1)==22,'business days with explicit holiday calendar')

   fm=fourier_terms_multi(5,[12.0_dp,6.0_dp],[2,1])
   call check(size(fm,2)==4,'multiple-season Fourier duplicate removal')
   fm=fourier_terms_multi(5,[4.0_dp],[2])
   call check(size(fm,2)==3,'Fourier Nyquist sine removal')
   tacf=tapered_acf(z,4)
   call check(size(tacf)==5 .and. abs(tacf(1)-1.0_dp)<1.0e-10_dp,'tapered ACF shrinkage')

   print '(a)','test_parity_v02: PASS'
contains
   function mean_cb(x,h) result(out)
      real(dp),intent(in)::x(:)
      integer,intent(in)::h
      type(forecast_result)::out
      out=mean_forecast(x,h)
   end function mean_cb

   function subset_mean_cb(x,train_mask) result(fitted)
      real(dp),intent(in)::x(:)
      logical,intent(in)::train_mask(:)
      real(dp),allocatable::fitted(:)
      real(dp)::mu
      mu=sum(pack(x,train_mask))/real(count(train_mask),dp)
      allocate(fitted(size(x)))
      fitted=mu
   end function subset_mean_cb

   subroutine check(ok,msg)
      logical,intent(in)::ok
      character(*),intent(in)::msg
      if(.not.ok)then
         write(*,'(a)')'FAIL: '//trim(msg)
         error stop 1
      end if
   end subroutine check
end program test_parity_v02
