program test_remaining_features
   use rugarch
   implicit none
   integer,parameter::n=240
   real(dp)::x(n),sigma(n),resid(n),pit(n),scores(n,3),losses(100,3),data2(n,2)
   real(dp)::hess(2,2),vcv(2,2),nw(3,3)
   integer::i,info,idx(n,4)
   type(test_result)::wb,alm
   type(sign_bias_result)::sb
   type(nyblom_result)::ny
   type(gof_result)::gof
   type(var_duration_result)::dur
   type(gmm_result)::gm
   type(hong_li_result)::hl
   type(mcs_result)::mcs
   type(garch_spec)::spec
   type(garch_fit_result)::fit
   type(multi_fit_result)::mf
   type(multi_forecast_result)::mfc
   type(rolling_forecast_result)::roll
   type(bootstrap_distribution_result)::bd
   type(arfima_spec)::aspec
   type(arfima_forecast_result)::afc
   type(arfima_order_result)::aorder

   call seed_rng(91234)
   do i=1,n
      x(i)=random_normal()
      sigma(i)=1.0_dp
      resid(i)=x(i)
      pit(i)=normal_cdf(x(i))
      scores(i,1)=x(i)
      scores(i,2)=x(i)*x(i)-1.0_dp
      scores(i,3)=x(i)**3
   end do

   wb=weighted_box_test(x,lag=8,fitdf=1)
   if(.not.(wb%p_value>=0.0_dp .and. wb%p_value<=1.0_dp))error stop 'weighted box'
   alm=arch_lm_test(x,5)
   if(.not.(alm%p_value>=0.0_dp .and. alm%p_value<=1.0_dp))error stop 'arch lm'
   sb=sign_bias_test(resid,sigma)
   if(sb%status/=0)error stop 'sign bias'
   ny=nyblom_test(scores)
   if(ny%status/=0 .or. size(ny%individual)/=3)error stop 'nyblom'
   gof=adjusted_pearson_gof(x,dist_norm,2.0_dp,1.0_dp,1.0_dp,[10,20])
   if(size(gof%p_value)/=2 .or. any(gof%p_value<0.0_dp) .or. any(gof%p_value>1.0_dp))error stop 'gof'

   nw=newey_west_covariance(scores,3,.true.)
   if(any(.not.(abs(nw)<huge(1.0_dp))))error stop 'newey west'
   hess=reshape([2.0_dp,0.0_dp,0.0_dp,4.0_dp],[2,2])
   call classical_covariance(hess,n,vcv,info)
   if(info/=0 .or. vcv(1,1)<=0.0_dp)error stop 'classical covariance'

   call bootstrap_indices(n,4,6,bootstrap_stationary,idx)
   if(any(idx<1) .or. any(idx>n))error stop 'stationary bootstrap'
   call bootstrap_indices(n,4,6,bootstrap_block,idx)
   if(any(idx<1) .or. any(idx>n))error stop 'block bootstrap'

   dur=var_duration_test(0.05_dp,x,-1.5_dp+0.0_dp*x)
   if(dur%status/=0 .or. dur%shape<=0.0_dp)error stop 'duration test'
   gm=gmm_test(x,2)
   if(gm%status/=0 .or. gm%joint_p_value<0.0_dp .or. gm%joint_p_value>1.0_dp)error stop 'gmm'
   hl=hong_li_test(pit,2)
   if(any(hl%p_value<0.0_dp) .or. any(hl%p_value>1.0_dp))error stop 'hong li'

   do i=1,100
      losses(i,1)=0.5_dp+0.1_dp*random_normal()
      losses(i,2)=0.7_dp+0.1_dp*random_normal()
      losses(i,3)=1.0_dp+0.1_dp*random_normal()
   end do
   mcs=mcs_test(losses,0.10_dp,40,5,bootstrap_stationary)
   if(.not.any(mcs%included_range) .or. .not.any(mcs%included_sq))error stop 'mcs'

   spec=make_garch_spec(1,1,model_sgarch,dist_norm)
   spec%mean=0.0_dp;spec%omega=0.05_dp;spec%alpha=0.08_dp;spec%beta=0.88_dp
   call simulate_garch(spec,n,data2(:,1),sigma,resid,burn_in=100)
   call simulate_garch(spec,n,data2(:,2),sigma,resid,burn_in=100)
   fit=fit_garch11(data2(:,1),max_iterations=120)
   if(.not.allocated(fit%sigma))error stop 'fit for workflows'
   mf=multifit_garch(data2,model_sgarch,1,1,max_iterations=80)
   if(size(mf%fit)/=2)error stop 'multifit'
   mfc=multiforecast_garch(mf,3)
   if(any(mfc%sigma<=0.0_dp))error stop 'multiforecast'
   roll=rolling_garch_forecast(data2(:,1),model_sgarch,1,1,forecast_length=5, &
      refit_every=5,max_iterations=80)
   if(size(roll%sigma)/=5 .or. any(roll%sigma<=0.0_dp))error stop 'rolling'
   bd=garch_parametric_distribution(fit,2,nobs=80,max_iterations=40)
   if(size(bd%status)/=2)error stop 'bootstrap distribution'
   aspec=make_arfima_spec(1,0);aspec%ar(1)=0.3_dp;aspec%d=0.2_dp;aspec%innovation_sd=1.0_dp
   afc=forecast_arfima(data2(:,1),aspec,4)
   if(size(afc%mean)/=4 .or. any(afc%sigma<=0.0_dp))error stop 'arfima forecast'
   aorder=auto_arfima(data2(:,1),2,1,'bic')
   if(aorder%best<1 .or. aorder%best>size(aorder%candidates))error stop 'auto arfima'

   if(abs(loss_squared_error(1.0_dp,2.0_dp)-1.0_dp)>1.0e-14_dp)error stop 'loss functions'
   print '(a)','remaining numerical feature tests passed'
end program test_remaining_features
