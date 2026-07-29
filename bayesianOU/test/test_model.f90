! SPDX-License-Identifier: MIT
program test_model
  use bayesianou
  implicit none
  integer, parameter :: t=70,s=2
  type(ou_input) :: inp
  type(ou_options) :: opt
  type(ou_summary) :: truth
  type(ou_fit_result) :: fit
  real(dp) :: phi(t,s),h(t,s)
  integer :: i,j
  allocate(inp%y(t,s),inp%x(t,s),inp%tmg(t),inp%com(t,s),inp%capital(t,s))
  allocate(truth%theta(s),truth%kappa(s),truth%a3(s),truth%beta0(s),truth%alpha(s),truth%rho(s),truth%sigma_eta(s))
  allocate(truth%kappa_p(s),truth%mu_const(s),truth%sigma_p(s),truth%a3_p(s))
  truth%theta=[0.1_dp,-0.1_dp];truth%kappa=[0.25_dp,0.18_dp];truth%a3=-0.04_dp;truth%beta0=[0.2_dp,0.1_dp]
  truth%alpha=-2.0_dp;truth%rho=0.75_dp;truth%sigma_eta=0.15_dp;truth%beta1=0.12_dp;truth%gamma=0.04_dp;truth%nu=8.0_dp
  truth%kappa_p=0.08_dp;truth%mu_const=0.0_dp;truth%sigma_p=0.2_dp;truth%a3_p=0.0_dp
  do i=1,t
    inp%tmg(i)=sin(0.13_dp*i)
    do j=1,s
      inp%x(i,j)=0.5_dp*sin(0.08_dp*i+0.2_dp*j)
      inp%com(i,j)=0.4_dp+0.02_dp*j+0.03_dp*cos(0.05_dp*i)
      inp%capital(i,j)=100.0_dp+5.0_dp*j+i
    end do
  end do
  opt%n_levels=1;opt%chains=2;opt%iterations=80;opt%warmup=40;opt%train_frac=0.7_dp;opt%seed=321
  call simulate_ou_nested(truth,opt,t,s,inp%x,inp%tmg,inp%com,inp%capital,seed=88,y=inp%y,phi=phi,h=h)
  call fit_ou_nonlinear_tmg(inp,opt,fit)
  if(fit%status/=status_ok) error stop trim(fit%message)
  if(.not.validate_ou_fit(fit)) error stop 'fit validation failed'
  if(.not.allocated(fit%diagnostics%loo%pointwise)) error stop 'LOO missing'
  if(any(.not.(fit%summary%kappa>0.0_dp))) error stop 'invalid kappa'
  print '(a,2f10.4)', 'estimated kappa: ',fit%summary%kappa
  print '(a,f10.4)', 'elpd_loo: ',fit%diagnostics%loo%elpd_loo
  print *, 'test_model: PASS'
end program test_model
