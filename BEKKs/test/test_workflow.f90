! SPDX-License-Identifier: MIT
program test_workflow
  use iso_fortran_env, only: int64
  use bekks
  implicit none
  type(rng_state) :: rng
  type(bekk_spec_type) :: spec
  type(bekk_fit_result) :: fit
  type(bekk_forecast_result) :: fc
  type(bekk_var_result) :: vr
  type(bekk_backtest_result) :: bt,bt_roll
  type(bekk_virf_result) :: vi
  type(bekk_portmanteau_result) :: pm
  type(bekk_parameters) :: par
  real(dp), allocatable :: theta(:),data(:,:),h(:,:,:),var(:,:),ret(:,:),mse(:)
  real(dp) :: signs(2),weights(2),shock(2)
  integer :: st,i

  signs=-1.0_dp;weights=[0.4_dp,0.6_dp];shock=[-1.5_dp,0.0_dp]
  par%model_type=bekk_scalar;par%asymmetric=.true.
  allocate(par%c(2,2),par%a(2,2),par%b(2,2),par%g(2,2))
  par%c=reshape([0.15_dp,0.03_dp,0.0_dp,0.12_dp],[2,2])
  par%a=0.0_dp;par%b=0.0_dp;par%g=0.0_dp
  par%a_scalar=0.08_dp;par%b_scalar=0.04_dp;par%g_scalar=0.82_dp
  theta=pack_parameters(par)
  call rng_seed(rng,1234567_int64)
  call simulate_sbekk_asymm(theta,140,2,rng,signs,0.25_dp,data,h,st)
  if(st/=bekk_ok)error stop 'simulate workflow'

  spec=bekk_spec(bekk_scalar,.true.,signs,theta)
  call bekk_fit(spec,data,fit,max_iter=2,criterion=1.0e-8_dp)
  if(fit%status/=bekk_ok .and. fit%status/=bekk_no_convergence)error stop 'fit workflow'
  if(.not.fit%stationary)error stop 'fit stationarity'
  if(.not.allocated(fit%h))error stop 'fit covariance'
  if(abs(sum(fit%score))>=huge(1.0_dp))error stop 'fit score'

  call forecast_bekk(fit,5,fc,0.90_dp)
  if(fc%status/=bekk_ok)error stop 'forecast workflow'
  if(any(fc%h(1,1,:)<=0.0_dp))error stop 'forecast positivity'

  call var_bekk_forecast(fit,fc,0.99_dp,vr,weights,'normal')
  if(vr%status/=bekk_ok)error stop 'var workflow'
  if(any(vr%value>=0.0_dp))error stop 'var sign'

  allocate(var(50,1),ret(50,1))
  do i=1,50
    ret(i,1)=data(90+i,1)*weights(1)+data(90+i,2)*weights(2)
    var(i,1)=-2.2_dp*sqrt(max(dot_product(weights,matmul(h(:,:,90+i),weights)),0.0_dp))
  end do
  call backtest_forecasts(ret,var,0.99_dp,bt)
  if(bt%status/=bekk_ok)error stop 'backtest workflow'
  if(bt%hit_rate(1)<0.0_dp .or. bt%hit_rate(1)>1.0_dp)error stop 'hit rate'

  call virf_bekk(fit,fit%h(:,:,size(fit%h,3)),shock,8,vi)
  if(vi%status/=bekk_ok)error stop 'virf workflow'
  if(size(vi%response,1)/=8)error stop 'virf dimensions'

  call portmanteau_test(fit,5,pm)
  if(pm%status/=bekk_ok)error stop 'portmanteau workflow'
  if(pm%p_value<0.0_dp .or. pm%p_value>1.0_dp)error stop 'portmanteau pvalue'

  call rolling_backtest(data,spec,100,0.99_dp,20,bt_roll,weights,'normal',max_iter=1)
  if(bt_roll%status/=bekk_ok)error stop 'rolling backtest workflow'
  if(size(bt_roll%var,1)/=40)error stop 'rolling backtest dimensions'

  call rng_seed(rng,8675309_int64)
  call bekk_mc_eval(theta,spec,[45],1,rng,mse,st,max_fit_iter=1)
  if(st/=bekk_ok)error stop 'Monte Carlo evaluation workflow'
  if(size(mse)/=1 .or. mse(1)<0.0_dp)error stop 'Monte Carlo evaluation result'

  print '(a)','test_workflow: PASS'
end program test_workflow
