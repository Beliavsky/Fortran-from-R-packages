! SPDX-License-Identifier: MIT
program bekks_demo
  use iso_fortran_env, only: int64
  use bekks
  implicit none

  type(rng_state) :: rng
  type(bekk_parameters) :: parameters
  type(bekk_spec_type) :: specification
  type(bekk_fit_result) :: fit
  type(bekk_forecast_result) :: forecast
  type(bekk_var_result) :: risk
  real(dp), allocatable :: theta(:), returns(:,:), h(:,:,:)
  real(dp) :: signs(2), weights(2)
  integer :: status

  signs=[-1.0_dp,-1.0_dp]
  weights=[0.55_dp,0.45_dp]
  parameters%model_type=bekk_scalar
  parameters%asymmetric=.true.
  allocate(parameters%c(2,2),parameters%a(2,2),parameters%b(2,2),parameters%g(2,2))
  parameters%c=reshape([0.11_dp,0.025_dp,0.0_dp,0.09_dp],[2,2])
  parameters%a=0.0_dp
  parameters%b=0.0_dp
  parameters%g=0.0_dp
  parameters%a_scalar=0.08_dp
  parameters%b_scalar=0.035_dp
  parameters%g_scalar=0.84_dp
  theta=pack_parameters(parameters)

  call rng_seed(rng,20260728_int64)
  call simulate_sbekk_asymm(theta,500,2,rng,signs,0.25_dp,returns,h,status)
  if(status/=bekk_ok)error stop 'simulation failed'

  specification=bekk_spec(bekk_scalar,.true.,signs,theta)
  call bekk_fit(specification,returns,fit,max_iter=20,use_qml=.true.)
  if(fit%status/=bekk_ok .and. fit%status/=bekk_no_convergence)error stop 'fit failed'

  call forecast_bekk(fit,5,forecast,confidence_level=0.95_dp)
  if(forecast%status/=bekk_ok)error stop 'forecast failed'
  call var_bekk_forecast(fit,forecast,0.99_dp,risk,weights,'normal')
  if(risk%status/=bekk_ok)error stop 'VaR failed'

  print '(a,l1)', 'stationary: ',fit%stationary
  print '(a,l1)', 'converged:  ',fit%converged
  print '(a,f14.6)', 'log likelihood: ',fit%log_likelihood
  print '(a,f14.6)', 'AIC:            ',fit%aic
  print '(a,f14.6)', 'BIC:            ',fit%bic
  print '(a,*(f11.6,1x))', 'estimated theta: ',fit%theta
  print '(a,*(f11.6,1x))', 'five-step portfolio VaR: ',risk%value(:,1)
end program bekks_demo
