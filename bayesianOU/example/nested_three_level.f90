! SPDX-License-Identifier: MIT
program nested_three_level
  use bayesianou
  implicit none
  integer,parameter::t=60,s=2
  type(ou_input)::data
  type(ou_options)::options
  type(ou_summary)::truth
  type(ou_fit_result)::fit
  real(dp)::phi(t,s),h(t,s),mu(t,s)
  integer::i,j
  allocate(data%y(t,s),data%x(t,s),data%tmg(t),data%com(t,s),data%capital(t,s),data%gprime(t),data%value_anchor(t,s))
  allocate(truth%theta(s),truth%kappa(s),truth%a3(s),truth%beta0(s),truth%alpha(s),truth%rho(s),truth%sigma_eta(s))
  allocate(truth%kappa_p(s),truth%mu_const(s),truth%sigma_p(s),truth%a3_p(s))
  truth%theta=0;truth%kappa=[0.30_dp,0.22_dp];truth%a3=-0.03_dp;truth%beta0=0;truth%alpha=-2.3_dp;truth%rho=0.65_dp;truth%sigma_eta=0.12_dp
  truth%kappa_p=[0.10_dp,0.07_dp];truth%mu_const=[0.1_dp,-0.1_dp];truth%sigma_p=0.08_dp;truth%a3_p=-0.02_dp
  truth%beta1=0.12_dp;truth%gamma=0.02_dp;truth%nu=10;truth%m1=0.20_dp;truth%m_v=0.08_dp;truth%sigma_phi_meas=0.05_dp
  do i=1,t
    data%tmg(i)=sin(0.09_dp*i);data%gprime(i)=cos(0.07_dp*i)
    do j=1,s
      data%x(i,j)=0.25_dp*sin(0.05_dp*i+0.2_dp*j);data%value_anchor(i,j)=0.15_dp*cos(0.04_dp*i+0.1_dp*j)
      data%com(i,j)=0.4_dp+0.03_dp*j;data%capital(i,j)=100+i+j
    end do
  end do
  options%n_levels=3;options%level_spec=ou_level_spec(level_both_full);options%chains=2;options%iterations=250;options%warmup=120
  call simulate_ou_nested(truth,options,t,s,data%x,data%tmg,data%com,data%capital,data%gprime,data%value_anchor,333,data%y,phi,h)
  call fit_ou_nested(data,options,fit)
  call extract_mu_trajectory(fit,data%gprime,data%value_anchor,mu)
  print '(a,2f10.4)','level-2 kappa: ',fit%summary%kappa_p
  print '(a,2f10.4)','level-2 mean at t=1: ',mu(1,:)
end program nested_three_level
