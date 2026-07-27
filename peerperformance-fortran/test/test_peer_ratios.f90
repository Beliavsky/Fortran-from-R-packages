! SPDX-License-Identifier: GPL-2.0-or-later
program test_peer_ratios
  use peerperformance, only: dp, adjust_pi, compute_pizero, optimal_lambda, compute_peer_ratios
  implicit none
  real(dp) :: p(2,6), d(2,6), t(2,6), pz(2), pp(2), pn(2), lu(2), lam(1)
  real(dp) :: pz2(2), pp2(2), pn2(2), lu2(2), a, z, l1, l2

  p(1,:)=[0.01_dp,0.04_dp,0.12_dp,0.35_dp,0.72_dp,0.91_dp]
  p(2,:)=[0.02_dp,0.08_dp,0.18_dp,0.45_dp,0.62_dp,0.84_dp]
  d(1,:)=[0.5_dp,0.4_dp,0.3_dp,0.1_dp,-0.1_dp,-0.2_dp]
  d(2,:)=[-0.5_dp,-0.4_dp,-0.3_dp,-0.2_dp,0.1_dp,0.2_dp]
  t(1,:)=[2.5_dp,2.0_dp,1.4_dp,0.5_dp,-0.4_dp,-1.0_dp]
  t(2,:)=[-2.4_dp,-1.8_dp,-1.2_dp,-0.7_dp,0.3_dp,0.8_dp]
  lam=0.5_dp
  call compute_peer_ratios(p,d,t,pz,pp,pn,lu,lambda=lam,n_boot=39, &
       gamma_pos=0.4_dp,gamma_neg=0.6_dp,seed=11,fast=.true.)
  call assert_true(all(abs(pz+pp+pn-1.0_dp)<2.0e-14_dp),'ratios sum to one')
  call assert_true(all(abs(lu-0.5_dp)<1.0e-15_dp),'fixed lambda')
  call assert_true(all(pz>=0.0_dp .and. pz<=1.0_dp),'pizero range')

  call compute_peer_ratios(p,d,t,pz2,pp2,pn2,lu2,lambda=lam,n_boot=39, &
       gamma_pos=0.2_dp,gamma_neg=0.8_dp,seed=11,fast=.true.)
  call assert_true(all(abs(pz-pz2)<1.0e-14_dp),'pizero threshold invariance')
  call assert_true(any(abs(pp-pp2)>1.0e-12_dp) .or. any(abs(pn-pn2)>1.0e-12_dp), &
       'positive/negative split changes')

  z=compute_pizero(p(1,:),0.5_dp,adjust=.false.)
  call assert_close(z,2.0_dp/3.0_dp,1.0e-14_dp,'unadjusted pizero')
  a=adjust_pi(z,6,0.5_dp,.true.)
  call assert_true(a>=0.0_dp .and. a<=1.0_dp,'adjusted pi range')
  l1=optimal_lambda(p(1,:),n_boot=31,seed=19)
  l2=optimal_lambda(p(1,:),n_boot=31,seed=19)
  call assert_close(l1,l2,0.0_dp,'lambda reproducibility')
  call assert_true(any(abs(l1-[0.3_dp,0.4_dp,0.5_dp,0.6_dp,0.7_dp])<1.0e-14_dp), &
       'lambda grid')
  print '(a)', 'test_peer_ratios: PASS'
contains
  subroutine assert_close(actual,expected,tol,label)
    real(dp), intent(in) :: actual,expected,tol
    character(len=*), intent(in) :: label
    if (abs(actual-expected)>tol*(1.0_dp+abs(expected))) then
      print '(a,3(1x,es24.16))',trim(label),actual,expected,abs(actual-expected); error stop 1
    end if
  end subroutine assert_close
  subroutine assert_true(value,label)
    logical, intent(in) :: value
    character(len=*), intent(in) :: label
    if (.not.value) then
      print '(a)', 'failed: '//trim(label); error stop 1
    end if
  end subroutine assert_true
end program test_peer_ratios
