! SPDX-License-Identifier: GPL-2.0-or-later
! This file is part of a modern Fortran translation of fExtremes.
! Copyright (C) 2026. This program is free software under GPL-2.0-or-later.
program test_distributions_fit
  use fextremes_kinds, only: dp
  use fextremes_rng, only: rng_state,seed_rng
  use fextremes_distributions
  use fextremes_fit
  use fextremes_risk
  implicit none
  type(rng_state)::rng
  type(gev_fit_result)::gf,gpwm,gum
  type(gpd_fit_result)::pf
  type(risk_result)::risk
  type(return_level_result)::rl
  real(dp),allocatable::x(:),y(:)
  real(dp)::p,q,m,v
  integer::i
  do i=1,9
    p=0.05_dp+0.1_dp*real(i-1,dp)
    q=gev_quantile(p,0.2_dp,1.0_dp,2.0_dp)
    call assert_close(gev_cdf(q,0.2_dp,1.0_dp,2.0_dp),p,2.0e-12_dp,'GEV inversion')
    q=gpd_quantile(p,-0.15_dp,0.5_dp,1.3_dp)
    call assert_close(gpd_cdf(q,-0.15_dp,0.5_dp,1.3_dp),p,2.0e-12_dp,'GPD inversion')
  end do
  call gev_moments(0.0_dp,0.0_dp,1.0_dp,m,v)
  call assert_close(m,0.5772156649015329_dp,1.0e-14_dp,'Gumbel mean')
  call assert_close(v,acos(-1.0_dp)**2/6.0_dp,1.0e-13_dp,'Gumbel variance')
  call gpd_moments(0.2_dp,0.0_dp,1.0_dp,m,v)
  call assert_close(m,1.25_dp,1.0e-14_dp,'GPD mean')
  call seed_rng(rng,4711); allocate(x(2500),y(2500))
  call gev_sample(rng,0.18_dp,0.4_dp,1.1_dp,x)
  call fit_gev(x,gpwm,'pwm')
  call assert_true(gpwm%converged,'GEV PWM convergence')
  call fit_gumbel(x,gum,'mle')
  call assert_true(gum%converged,'Gumbel MLE convergence')
  call fit_gev(x,gf,'mle')
  call assert_true(gf%converged,'GEV fit convergence')
  call assert_true(abs(gf%xi-0.18_dp)<0.12_dp,'GEV xi recovery')
  call assert_true(abs(gf%mu-0.4_dp)<0.15_dp,'GEV mu recovery')
  call assert_true(abs(gf%beta-1.1_dp)<0.15_dp,'GEV beta recovery')
  call gev_return_level(gf,20.0_dp,0.95_dp,rl)
  call assert_true(rl%estimate>maxval(x)*0.25_dp,'GEV return level')
  call gpd_sample(rng,0.22_dp,0.0_dp,1.3_dp,y)
  call fit_gpd(y,0.0_dp,pf,'mle','expected')
  call assert_true(pf%converged,'GPD fit convergence')
  call assert_true(abs(pf%xi-0.22_dp)<0.10_dp,'GPD xi recovery')
  call assert_true(abs(pf%beta-1.3_dp)<0.15_dp,'GPD beta recovery')
  call gpd_risk_measures(pf,[0.99_dp,0.995_dp],risk)
  call assert_true(all(risk%value_at_risk(2:)>=risk%value_at_risk(:1)),'GPD risk monotonicity')
  call assert_true(all(risk%expected_shortfall>risk%value_at_risk),'GPD ES above VaR')
  call assert_true(sample_var(y,0.05_dp,.true.) > sample_var(y,0.05_dp), &
    'sample VaR tails')
  call assert_true(sample_cvar(y,0.05_dp,.true.) >= sample_var(y,0.05_dp,.true.), &
    'sample CVaR')
  print '(a)','Distribution, fitting, and risk tests passed.'
contains
  subroutine assert_close(a,b,tol,msg)
    real(dp),intent(in)::a,b,tol; character(len=*),intent(in)::msg
    if(abs(a-b)>tol) then; print *,trim(msg),a,b; error stop 1; end if
  end subroutine
  subroutine assert_true(cond,msg)
    logical,intent(in)::cond; character(len=*),intent(in)::msg
    if(.not.cond) then; print *,trim(msg); error stop 1; end if
  end subroutine
end program test_distributions_fit
