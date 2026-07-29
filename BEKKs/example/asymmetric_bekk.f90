! SPDX-License-Identifier: MIT
program asymmetric_bekk
  use iso_fortran_env, only: int64
  use bekks
  implicit none

  type(rng_state) :: rng
  type(bekk_parameters) :: parameters
  type(bekk_spec_type) :: specification
  type(bekk_fit_result) :: fit
  real(dp), allocatable :: theta(:), returns(:,:), h(:,:,:)
  real(dp) :: signs(2)
  integer :: status

  signs=[-1.0_dp,1.0_dp]
  parameters%model_type=bekk_diagonal
  parameters%asymmetric=.true.
  allocate(parameters%c(2,2),parameters%a(2,2),parameters%b(2,2),parameters%g(2,2))
  parameters%c=reshape([0.12_dp,0.02_dp,0.0_dp,0.10_dp],[2,2])
  parameters%a=0.0_dp
  parameters%b=0.0_dp
  parameters%g=0.0_dp
  parameters%a(1,1)=0.17_dp
  parameters%a(2,2)=0.12_dp
  parameters%b(1,1)=0.05_dp
  parameters%b(2,2)=0.04_dp
  parameters%g(1,1)=0.78_dp
  parameters%g(2,2)=0.82_dp
  theta=pack_parameters(parameters)

  call rng_seed(rng,77123_int64)
  call simulate_dbekk_asymm(theta,350,2,rng,signs,0.25_dp,returns,h,status)
  if(status/=bekk_ok)error stop 'simulation failed'
  specification=bekk_spec(bekk_diagonal,.true.,signs,theta)
  call bekk_fit(specification,returns,fit,max_iter=15)
  if(fit%status/=bekk_ok .and. fit%status/=bekk_no_convergence)error stop 'fit failed'

  print '(a,*(f10.5,1x))','true theta:      ',theta
  print '(a,*(f10.5,1x))','estimated theta: ',fit%theta
  print '(a,f12.5)','log likelihood:  ',fit%log_likelihood
end program asymmetric_bekk
