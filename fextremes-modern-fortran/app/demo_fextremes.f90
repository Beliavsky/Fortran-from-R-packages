! SPDX-License-Identifier: GPL-2.0-or-later
! This file is part of a modern Fortran translation of fExtremes.
! Copyright (C) 2026. This program is free software under GPL-2.0-or-later.
program demo_fextremes
  use fextremes_kinds, only: dp
  use fextremes_rng, only: rng_state,seed_rng
  use fextremes_distributions, only: gev_random,gpd_random
  use fextremes_fit, only: gev_fit_result,gpd_fit_result,fit_gev,fit_gpd
  use fextremes_risk, only: risk_result,gpd_risk_measures
  implicit none
  type(rng_state)::rng
  type(gev_fit_result)::gevfit
  type(gpd_fit_result)::gpdfit
  type(risk_result)::risk
  real(dp),allocatable::blockmax(:),loss(:)
  integer::i
  allocate(blockmax(1500),loss(2500)); call seed_rng(rng,4711)
  do i=1,size(blockmax); blockmax(i)=gev_random(rng,0.15_dp,0.0_dp,1.0_dp); end do
  call fit_gev(blockmax,gevfit,'mle')
  print '(a,3(1x,f10.5))','GEV xi mu beta:',gevfit%xi,gevfit%mu,gevfit%beta
  do i=1,size(loss); loss(i)=gpd_random(rng,0.20_dp,0.0_dp,1.2_dp); end do
  call fit_gpd(loss,0.0_dp,gpdfit,'mle','expected')
  call gpd_risk_measures(gpdfit,[0.99_dp,0.995_dp,0.999_dp],risk)
  print '(a,2(1x,f10.5))','GPD xi beta:',gpdfit%xi,gpdfit%beta
  do i=1,size(risk%probability)
    print '(f8.4,2(1x,f12.5))',risk%probability(i),risk%value_at_risk(i),risk%expected_shortfall(i)
  end do
end program demo_fextremes
