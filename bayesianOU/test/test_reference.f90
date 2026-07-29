! SPDX-License-Identifier: MIT
program test_reference
  use bayesianou
  implicit none
  integer,parameter::t=5,s=1
  type(ou_input)::inp
  type(ou_fit_result)::fit
  type(ou_summary)::summ
  real(dp)::ll,pw(4)
  allocate(inp%y(t,s),inp%x(t,s),inp%tmg(t),inp%com(t,s),inp%capital(t,s))
  inp%y(:,1)=[0.0_dp,0.1_dp,0.15_dp,0.05_dp,0.0_dp]
  inp%x(:,1)=[0.2_dp,0.1_dp,0.0_dp,-0.1_dp,-0.2_dp]
  inp%tmg=[0.0_dp,0.5_dp,-0.5_dp,1.0_dp,0.2_dp]
  inp%com=0.4_dp;inp%capital(:,1)=[100.0_dp,101.0_dp,102.0_dp,103.0_dp,104.0_dp]
  allocate(fit%zy%mz(t,s),fit%zx%mz(t,s),fit%ztmg(t),fit%com_wmean(s),fit%com_wsd(s),fit%h_median(t,s),fit%phi_median(t,s))
  fit%zy%mz=inp%y;fit%zx%mz=inp%x;fit%ztmg=inp%tmg;fit%com_wmean=0.4_dp;fit%com_wsd=1.0_dp;fit%h_median=-2.0_dp;fit%phi_median=inp%x
  fit%t_lik=t;fit%options%n_levels=1;fit%options%com_in_mean=.false.
  allocate(summ%theta(s),summ%kappa(s),summ%a3(s),summ%beta0(s),summ%alpha(s),summ%rho(s),summ%sigma_eta(s))
  allocate(summ%kappa_p(s),summ%mu_const(s),summ%sigma_p(s),summ%a3_p(s))
  summ%theta=0.02_dp;summ%kappa=0.3_dp;summ%a3=-0.04_dp;summ%beta0=0.15_dp;summ%beta1=0.1_dp;summ%gamma=0;summ%nu=8.0_dp
  summ%alpha=-2.0_dp;summ%rho=0.7_dp;summ%sigma_eta=0.2_dp;summ%kappa_p=0.1_dp;summ%mu_const=0;summ%sigma_p=0.2_dp;summ%a3_p=0
  ll=ou_log_likelihood(inp,fit,summ,pw)
  if(abs(ll-0.1526471942737592_dp)>1e-12_dp)error stop 'fixed likelihood reference failed'
  if(maxval(abs(pw-[0.037787823642837066_dp,0.032895868388259045_dp,0.03446546735724757_dp,0.047498034885415506_dp]))>1e-12_dp) &
    error stop 'pointwise likelihood reference failed'
  print *, 'test_reference: PASS'
end program test_reference
