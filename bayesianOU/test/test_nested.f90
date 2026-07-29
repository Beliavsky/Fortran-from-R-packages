! SPDX-License-Identifier: MIT
program test_nested
  use bayesianou
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  implicit none
  integer,parameter::t=55,s=2
  type(ou_input)::inp
  type(ou_options)::opt
  type(ou_summary)::truth
  type(ou_fit_result)::fit
  real(dp)::phi(t,s),h(t,s),mu(t,s)
  integer::i,j
  allocate(inp%y(t,s),inp%x(t,s),inp%tmg(t),inp%com(t,s),inp%capital(t,s),inp%gprime(t),inp%value_anchor(t,s))
  allocate(truth%theta(s),truth%kappa(s),truth%a3(s),truth%beta0(s),truth%alpha(s),truth%rho(s),truth%sigma_eta(s))
  allocate(truth%kappa_p(s),truth%mu_const(s),truth%sigma_p(s),truth%a3_p(s))
  truth%theta=0;truth%kappa=[0.35_dp,0.25_dp];truth%a3=-0.03_dp;truth%beta0=0;truth%alpha=-2.5_dp;truth%rho=0.6_dp;truth%sigma_eta=0.12_dp
  truth%kappa_p=[0.12_dp,0.08_dp];truth%mu_const=[0.1_dp,-0.1_dp];truth%sigma_p=0.08_dp;truth%a3_p=0.0_dp
  truth%beta1=0.15_dp;truth%gamma=0.0_dp;truth%nu=20.0_dp;truth%m1=0.25_dp;truth%m_v=0.1_dp;truth%sigma_phi_meas=0.05_dp
  do i=1,t
    inp%tmg(i)=sin(0.11_dp*i);inp%gprime(i)=cos(0.07_dp*i)
    do j=1,s
      inp%x(i,j)=0.3_dp*sin(0.05_dp*i+0.3_dp*j);inp%value_anchor(i,j)=0.2_dp*cos(0.04_dp*i+0.2_dp*j)
      inp%com(i,j)=0.5_dp+0.02_dp*j;inp%capital(i,j)=100.0_dp+i+j
    end do
  end do
  opt%n_levels=3;opt%chains=2;opt%iterations=70;opt%warmup=35;opt%seed=66;opt%train_frac=0.72_dp
  opt%level_spec=ou_level_spec(level_both_full)
  call simulate_ou_nested(truth,opt,t,s,inp%x,inp%tmg,inp%com,inp%capital,inp%gprime,inp%value_anchor,77,inp%y,phi,h)
  call fit_ou_nested(inp,opt,fit)
  if(fit%status/=status_ok)error stop trim(fit%message)
  if(any(fit%summary%kappa_p<=0))error stop 'invalid kappa_p'
  call extract_mu_trajectory(fit,inp%gprime,inp%value_anchor,mu)
  if(.not.all(ieee_is_finite(mu)))error stop 'nonfinite mu'
  if(size(fit%diagnostics%oos)/=3)error stop 'missing oos'
  print '(a,2f9.4)','nested kappa_p: ',fit%summary%kappa_p
  print *, 'test_nested: PASS'
end program test_nested
