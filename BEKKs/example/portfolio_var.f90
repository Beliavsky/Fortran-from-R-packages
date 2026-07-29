! SPDX-License-Identifier: MIT
program portfolio_var
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
  real(dp) :: weights(3)
  integer :: status

  weights=[0.50_dp,0.30_dp,0.20_dp]
  parameters%model_type=bekk_scalar
  parameters%asymmetric=.false.
  allocate(parameters%c(3,3),parameters%a(3,3),parameters%b(3,3),parameters%g(3,3))
  parameters%c=reshape([0.10_dp,0.015_dp,-0.005_dp,0.0_dp,0.09_dp,0.012_dp, &
    0.0_dp,0.0_dp,0.075_dp],[3,3])
  parameters%a=0.0_dp
  parameters%b=0.0_dp
  parameters%g=0.0_dp
  parameters%a_scalar=0.10_dp
  parameters%g_scalar=0.85_dp
  theta=pack_parameters(parameters)

  call rng_seed(rng,442211_int64)
  call simulate_sbekk(theta,600,3,rng,returns,h,status)
  if(status/=bekk_ok)error stop 'simulation failed'
  specification=bekk_spec(bekk_scalar,.false.,initial_theta=theta)
  call bekk_fit(specification,returns,fit,max_iter=12)
  if(fit%status/=bekk_ok .and. fit%status/=bekk_no_convergence)error stop 'fit failed'
  call forecast_bekk(fit,10,forecast)
  call var_bekk_forecast(fit,forecast,0.99_dp,risk,weights,'t')
  if(risk%status/=bekk_ok)error stop 'VaR failed'

  print '(a,*(f11.6,1x))','10-step portfolio VaR: ',risk%value(:,1)
end program portfolio_var
