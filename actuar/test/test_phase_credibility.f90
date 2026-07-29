! SPDX-License-Identifier: GPL-2.0-or-later
program test_phase_credibility
  use actuar, only : dp, dphtype, pphtype, mphtype, mgfphtype, &
    credibility_result, buhlmann_straub, bayes_poisson_gamma
  implicit none
  real(dp) :: pi0(1),tmat(1,1),ratios(3,3),weights(3,3)
  type(credibility_result) :: fit,bayes

  pi0=1.0_dp; tmat(1,1)=-2.0_dp
  call assert_close(dphtype(1.0_dp,pi0,tmat),2.0_dp*exp(-2.0_dp),1.0e-12_dp)
  call assert_close(pphtype(1.0_dp,pi0,tmat),1.0_dp-exp(-2.0_dp),1.0e-12_dp)
  call assert_close(mphtype(2,pi0,tmat),0.5_dp,1.0e-12_dp)
  call assert_close(mgfphtype(0.5_dp,pi0,tmat),4.0_dp/3.0_dp,1.0e-12_dp)

  ratios=reshape([1.0_dp,1.3_dp,0.8_dp,1.1_dp,1.4_dp,0.9_dp, &
    0.9_dp,1.2_dp,1.0_dp],[3,3])
  weights=reshape([1.0_dp,2.0_dp,1.5_dp,1.0_dp,1.0_dp,2.0_dp, &
    2.0_dp,1.0_dp,1.0_dp],[3,3])
  fit=buhlmann_straub(ratios,weights)
  if(.not.fit%ok) error stop 1
  if(any(fit%weights<0.0_dp) .or. any(fit%weights>1.0_dp)) error stop 1
  bayes=bayes_poisson_gamma([2.0_dp,3.0_dp,1.0_dp],2.0_dp,0.5_dp)
  if(.not.bayes%ok) error stop 1
  call assert_close(bayes%collective_mean,1.0_dp,1.0e-14_dp)
  call assert_close(bayes%weights(1),0.6_dp,1.0e-14_dp)

  print '(a)', 'test_phase_credibility: PASS'
contains
  subroutine assert_close(actual,expected,tol)
    real(dp), intent(in) :: actual,expected,tol
    if(abs(actual-expected)>tol*max(1.0_dp,abs(expected))) then
      print '(a,3es24.15)', 'mismatch: ',actual,expected,abs(actual-expected)
      error stop 1
    end if
  end subroutine assert_close
end program test_phase_credibility
