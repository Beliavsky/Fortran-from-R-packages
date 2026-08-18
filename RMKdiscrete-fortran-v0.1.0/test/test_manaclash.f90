! SPDX-License-Identifier: GPL-2.0-or-later
! Modern Fortran translation of RMKdiscrete 0.1 by Robert M. Kirkpatrick.
! See LICENSE, COPYING, and TRANSLATION_NOTES.md for provenance.
program test_manaclash
  use rmkdiscrete, only : dp, dmanaclash_xyn, dmanaclash_dmg, dmanaclash_net
  implicit none
  real(dp)::s
  integer::x,y,failures
  failures=0
  call check_close(dmanaclash_xyn(1,1,1),0.0625_dp,1.0e-15_dp,'joint x,y,N')
  call check_close(dmanaclash_dmg(1,1),0.09375_dp,1.0e-15_dp,'marginal damage')
  call check_close(dmanaclash_dmg(1,1,n=1),1.0_dp/3.0_dp,2.0e-15_dp,'conditional damage')
  call check_close(dmanaclash_net(0,rel_eps=1.0e-14_dp),0.447213595499957939_dp,2.0e-13_dp,'net damage zero')
  s=0.0_dp
  do x=0,2
    do y=0,2
      s=s+dmanaclash_xyn(x,y,2)
    end do
  end do
  call check_close(s,0.140625_dp,1.0e-15_dp,'joint mass for N=2')
  if(failures==0) then
    print '(a)','test_manaclash: PASS'
  else
    print '(a,i0)','test_manaclash: FAIL ',failures
    error stop 1
  end if
contains
  subroutine check_close(actual,expected,tol,label)
    real(dp),intent(in)::actual,expected,tol
    character(*),intent(in)::label
    if(abs(actual-expected)>tol) then
      failures=failures+1
      print '(a,2es24.15)','FAIL '//trim(label)//': ',actual,expected
    end if
  end subroutine check_close
end program test_manaclash
