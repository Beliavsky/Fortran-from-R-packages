! SPDX-License-Identifier: GPL-3.0-or-later
program test_profiles
  use acdm
  implicit none
  integer,parameter::n=300
  real(dp)::x(n),par(3),seq(9)
  type(acd_order)::o
  type(rng_state)::rng
  type(hazard_result)::hr
  type(likelihood_profile_result)::lr
  integer::i,st,fail
  fail=0;o=acd_order(1,0,1);par=[0.2_dp,0.15_dp,0.7_dp]
  call seed_rng(rng,71234)
  call simulate_acd(n,MODEL_ACD,o,par,DIST_EXPONENTIAL,[real(dp)::],x,st,rng,burn=200)
  call hazard_diagnostics(x,DIST_EXPONENTIAL,[real(dp)::],.true.,hr,breaks=15)
  call ok(hr%status==ACDM_SUCCESS.and.size(hr%empirical_hazard)==14,'hazard')
  call ok(all(hr%implied_hazard>0._dp),'implied hazard')
  do i=1,9;seq(i)=0.1_dp+0.025_dp*real(i-1,dp);end do
  call likelihood_profile(x,MODEL_ACD,o,par,DIST_EXPONENTIAL,[real(dp)::],.true.,2,seq,lr)
  call ok(lr%status==ACDM_SUCCESS.and.size(lr%loglik,1)==9,'profile')
  call ok(maxloc(lr%loglik(:,1),dim=1)>=2.and.maxloc(lr%loglik(:,1),dim=1)<=8,'profile interior')
  if(fail>0)error stop 'test_profiles failed'
  print '(a)','test_profiles: PASS'
contains
  subroutine ok(cond,label)
    logical,intent(in)::cond;character(*),intent(in)::label
    if(.not.cond)then;fail=fail+1;print *,'FAIL ',trim(label);end if
  end subroutine
end program
