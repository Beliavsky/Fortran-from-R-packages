! SPDX-License-Identifier: GPL-2.0-or-later
! This file is part of a modern Fortran translation of fExtremes.
! Copyright (C) 2026. This program is free software under GPL-2.0-or-later.
program test_advanced
  use fextremes_kinds, only: dp
  use fextremes_rng, only: rng_state,seed_rng
  use fextremes_distributions, only: gev_random,gpd_random
  use fextremes_fit, only: gev_fit_result,gpd_fit_result,fit_gev,fit_gumbel,fit_gpd
  use fextremes_risk, only: threshold_stability_result,return_level_result,tail_profile_result, &
    gpd_threshold_stability,gev_return_level_profile,gpd_profile_risk,gpd_tail_curve
  use fextremes_extremal_index, only: theta_result,block_theta,run_theta
  implicit none
  type(rng_state)::rng
  type(gev_fit_result)::gf,gum
  type(gpd_fit_result)::pf
  type(threshold_stability_result)::stable
  type(return_level_result)::rl
  type(tail_profile_result)::tp
  type(theta_result)::tr
  real(dp),allocatable::x(:),y(:),surv(:)
  real(dp)::thr(4),probs(2)
  integer::i
  call seed_rng(rng,77); allocate(x(350),y(600),surv(3))
  do i=1,size(x); x(i)=gev_random(rng,0.0_dp,1.0_dp,0.8_dp); end do
  call fit_gumbel(x,gum,'pwm'); call assert_true(gum%converged,'Gumbel PWM')
  call fit_gev(x,gf,'mle'); call gev_return_level_profile(x,gf,10.0_dp,0.90_dp,rl,31)
  call assert_true(rl%lower<=rl%estimate .and. rl%upper>=rl%estimate,'GEV profile interval')
  do i=1,size(y); y(i)=gpd_random(rng,0.15_dp,0.0_dp,1.0_dp); end do
  call fit_gpd(y,0.0_dp,pf,'pwm'); call assert_true(pf%converged,'GPD PWM')
  call fit_gpd(y,0.0_dp,pf,'mle','observed')
  call gpd_profile_risk(pf,0.99_dp,0.90_dp,tp,31)
  call assert_true(tp%var_lower<=tp%var_estimate .and. tp%var_upper>=tp%var_estimate,'GPD VaR profile')
  call assert_true(tp%es_lower<=tp%es_estimate .and. tp%es_upper>=tp%es_estimate,'GPD ES profile')
  call gpd_tail_curve(pf,[0.5_dp,1.0_dp,2.0_dp],surv)
  call assert_true(all(surv(2:)<=surv(:2)),'tail curve monotonicity')
  thr=[0.2_dp,0.5_dp,0.8_dp,1.1_dp]
  call gpd_threshold_stability(y,thr,0.99_dp,stable)
  call assert_true(count(stable%converged)>=3,'threshold stability fits')
  probs=[0.90_dp,0.95_dp]
  call block_theta(y,20,probs,tr); call assert_true(all(tr%theta>=0.0_dp),'block theta')
  call run_theta(y,5,probs,tr); call assert_true(all(tr%theta>=0.0_dp),'run theta')
  print '(a)','Profile-likelihood and stability tests passed.'
contains
  subroutine assert_true(cond,msg)
    logical,intent(in)::cond; character(len=*),intent(in)::msg
    if(.not.cond) then; print *,trim(msg); error stop 1; end if
  end subroutine
end program test_advanced
