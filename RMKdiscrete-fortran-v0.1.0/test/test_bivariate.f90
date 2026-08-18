! SPDX-License-Identifier: GPL-2.0-or-later
! Modern Fortran translation of RMKdiscrete 0.1 by Robert M. Kirkpatrick.
! See LICENSE, COPYING, and TRANSLATION_NOTES.md for provenance.
program test_bivariate
  use rmkdiscrete, only : dp, dbilgp, dbinegbin, dlgp, dnegbin, log_moments2, &
    bilgp_logmv, binegbin_logmv
  implicit none
  real(dp)::theta(3),lambda(3),nu(3),p(3),v
  type(log_moments2)::m
  integer::failures
  failures=0
  theta=[1.0_dp,2.0_dp,1.5_dp]
  lambda=[0.2_dp,0.1_dp,0.0_dp]
  call check_close(dbilgp(2,1,theta,lambda),0.0464739085336362504_dp,3.0e-15_dp,'bivariate LGP PMF')
  theta(1)=0.0_dp
  v=dlgp(2,theta(2),lambda(2))*dlgp(1,theta(3),lambda(3))
  call check_close(dbilgp(2,1,theta,lambda),v,2.0e-15_dp,'bivariate LGP independence')
  m=bilgp_logmv(theta,lambda,tol=1.0e-12_dp)
  call check_close(m%cov,0.0_dp,2.0e-10_dp,'bivariate LGP independent log covariance')

  nu=[1.2_dp,2.0_dp,1.5_dp]
  p=[0.55_dp,0.6_dp,0.7_dp]
  call check_close(dbinegbin(3,2,nu,p),0.0384489640131341086_dp,3.0e-15_dp,'bivariate negative-binomial PMF')
  nu(1)=0.0_dp
  v=dnegbin(3,nu(2),p=p(2))*dnegbin(2,nu(3),p=p(3))
  call check_close(dbinegbin(3,2,nu,p),v,2.0e-15_dp,'bivariate negative-binomial independence')
  m=binegbin_logmv(nu,p,tol=1.0e-12_dp)
  call check_close(m%cov,0.0_dp,2.0e-10_dp,'bivariate NB independent log covariance')

  if(failures==0) then
    print '(a)','test_bivariate: PASS'
  else
    print '(a,i0)','test_bivariate: FAIL ',failures
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
end program test_bivariate
