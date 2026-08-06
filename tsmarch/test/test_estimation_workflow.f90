program test_estimation_workflow
  use ghyp_kinds, only : dp, i8
  use tsgarch, only : garch_spec, fit_options
  use tsmarch
  use test_support
  implicit none
  real(dp), allocatable :: z(:, :), data(:, :)
  real(dp) :: corr(2,2), mean_abs
  type(dcc_spec) :: ds
  type(dcc_parameters) :: dp0
  type(dcc_simulation) :: sim
  type(dcc_fit) :: df
  type(dcc_forecast) :: forecast
  type(garch_spec) :: gs
  type(fit_options) :: opt
  type(gogarch_spec) :: gos
  type(gogarch_fit) :: gofit
  type(gogarch_forecast) :: goforecast
  integer :: i,n

  n=180
  corr=reshape([1.0_dp,0.35_dp,0.35_dp,1.0_dp],[2,2])
  ds%constant_correlation=.false.
  ds%distribution='mvn'
  allocate(dp0%alpha(1),dp0%gamma(0),dp0%beta(1))
  dp0%alpha=0.04_dp
  dp0%beta=0.93_dp
  dp0%shape=8.0_dp
  sim=simulate_dcc_innovations(ds,dp0,corr,n,seed=2026_i8)
  call assert_true(sim%status==tsm_success,'workflow simulation')
  z=sim%innovations(:,:,1)
  allocate(data(n,2))
  do i=1,n
    data(i,1)=0.012_dp*z(i,1)
    data(i,2)=0.018_dp*z(i,2)
  end do

  gs%model='ewma'
  gs%distribution='norm'
  gs%constant=.false.
  gs%p=1
  gs%q=1
  opt%max_iterations=350
  opt%compute_inference=.false.
  ds%constant_correlation=.true.
  df=estimate_dcc(data,gs,ds,opt)
  call assert_true(df%status==tsm_success,'constant DCC estimation')
  call assert_true(size(df%marginals)==2,'marginal fit count')
  forecast=forecast_dcc(data,df,3,paths=60,seed=11_i8)
  call assert_true(forecast%status==tsm_success,'DCC forecast')
  call assert_true(all(forecast%sigma>0.0_dp),'DCC forecast sigma positive')
  call assert_true(size(forecast%simulated,3)==60,'DCC forecast paths')

  gos%ica_method='fastica'
  gos%factor_spec=gs
  gofit=estimate_gogarch(data,gos,opt,seed=5_i8)
  call assert_true(gofit%status==tsm_success,'GO-GARCH estimation')
  call assert_true(size(gofit%factors)==2,'GO-GARCH factor count')
  goforecast=forecast_gogarch(data,gofit,2,paths=40,seed=7_i8)
  call assert_true(goforecast%status==tsm_success,'GO-GARCH forecast')
  call assert_true(all(goforecast%factor_variance>0.0_dp),'GO-GARCH factor variance positive')
  mean_abs=sum(abs(goforecast%simulated))/real(size(goforecast%simulated),dp)
  call assert_true(mean_abs>0.0_dp,'GO-GARCH simulated paths nonzero')

  call finish_test('test_estimation_workflow')
end program test_estimation_workflow
