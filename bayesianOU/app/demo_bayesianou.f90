! SPDX-License-Identifier: MIT
program demo_bayesianou
  use bayesianou
  implicit none
  integer,parameter::t=80,s=2
  type(ou_input)::data
  type(ou_options)::options
  type(ou_summary)::truth
  type(ou_fit_result)::fit
  type(stability_result)::stability
  real(dp)::phi(t,s),h(t,s)
  integer::i,j
  call allocate_demo(data,truth)
  do i=1,t
    data%tmg(i)=sin(0.10_dp*i)
    do j=1,s
      data%x(i,j)=0.4_dp*sin(0.06_dp*i+0.3_dp*j)
      data%com(i,j)=0.45_dp+0.03_dp*j+0.02_dp*cos(0.05_dp*i)
      data%capital(i,j)=100.0_dp+real(i+5*j,dp)
    end do
  end do
  options%n_levels=1;options%chains=2;options%iterations=300;options%warmup=150;options%seed=1234
  call simulate_ou_nested(truth,options,t,s,data%x,data%tmg,data%com,data%capital,seed=99,y=data%y,phi=phi,h=h)
  call fit_ou_nonlinear_tmg(data,options,fit)
  call kappa_stability_evidence(fit,stability)
  print '(a,2f10.4)','kappa posterior medians: ',fit%summary%kappa
  print '(a,f10.4)','beta1 posterior median: ',fit%summary%beta1
  print '(a,f10.4)','PSIS-style elpd_loo: ',fit%diagnostics%loo%elpd_loo
  print '(a,f8.3)','P(all kappa in (0,1)): ',stability%probability
contains
  subroutine allocate_demo(data,truth)
    type(ou_input),intent(out)::data
    type(ou_summary),intent(out)::truth
    allocate(data%y(t,s),data%x(t,s),data%tmg(t),data%com(t,s),data%capital(t,s))
    allocate(truth%theta(s),truth%kappa(s),truth%a3(s),truth%beta0(s),truth%alpha(s),truth%rho(s),truth%sigma_eta(s))
    allocate(truth%kappa_p(s),truth%mu_const(s),truth%sigma_p(s),truth%a3_p(s))
    truth%theta=[0.1_dp,-0.1_dp];truth%kappa=[0.25_dp,0.18_dp];truth%a3=-0.04_dp;truth%beta0=[0.18_dp,0.10_dp]
    truth%alpha=-2.0_dp;truth%rho=0.75_dp;truth%sigma_eta=0.15_dp;truth%beta1=0.10_dp;truth%gamma=0.03_dp;truth%nu=8.0_dp
    truth%kappa_p=0.08_dp;truth%mu_const=0;truth%sigma_p=0.2_dp;truth%a3_p=0
  end subroutine allocate_demo
end program demo_bayesianou
