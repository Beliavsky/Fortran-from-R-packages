! SPDX-License-Identifier: GPL-2.0-or-later
! Modern Fortran translation of RMKdiscrete 0.1 by Robert M. Kirkpatrick.
! See LICENSE, COPYING, and TRANSLATION_NOTES.md for provenance.
program test_rng
  use rmkdiscrete, only : dp, rlgp, rnegbin, rbilgp, rbinegbin, rmanaclash, lgp_findmax
  implicit none
  real(dp)::theta(3),lambda(3),nu(3),p(3)
  integer::i,x,xy(2),xyz(3),failures
  failures=0
  do i=1,1000
    x=rlgp(2.0_dp,-0.2_dp)
    if(x<0 .or. real(x,dp)>lgp_findmax(2.0_dp,-0.2_dp)) failures=failures+1
    x=rnegbin(2.5_dp,0.4_dp)
    if(x<0) failures=failures+1
  end do
  theta=[1.0_dp,2.0_dp,1.5_dp]
  lambda=[0.2_dp,0.1_dp,0.0_dp]
  call rbilgp(theta,lambda,xy)
  if(any(xy<0)) failures=failures+1
  nu=[1.2_dp,2.0_dp,1.5_dp]
  p=[0.55_dp,0.6_dp,0.7_dp]
  call rbinegbin(nu,p,xy)
  if(any(xy<0)) failures=failures+1
  call rmanaclash(out=xyz)
  if(any(xyz<0) .or. xyz(1)>xyz(3) .or. xyz(2)>xyz(3)) failures=failures+1
  call rmanaclash(out=xyz,n=5)
  if(xyz(3)/=5 .or. xyz(1)>5 .or. xyz(2)>5) failures=failures+1
  if(failures==0) then
    print '(a)','test_rng: PASS'
  else
    print '(a,i0)','test_rng: FAIL ',failures
    error stop 1
  end if
end program test_rng
