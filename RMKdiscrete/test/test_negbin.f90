! SPDX-License-Identifier: GPL-2.0-or-later
! Modern Fortran translation of RMKdiscrete 0.1 by Robert M. Kirkpatrick.
! See LICENSE, COPYING, and TRANSLATION_NOTES.md for provenance.
program test_negbin
  use rmkdiscrete, only : dp, dnegbin, negbin_from_nu_p, negbin_from_mu_sigma2, &
    negbin_from_mu_nu, negbin_from_mu_p, negbin_from_sigma2_p
  implicit none
  real(dp)::mu,sigma2,nu,p
  integer::failures
  failures=0
  call check_close(dnegbin(4,2.5_dp,p=0.4_dp),0.118338754598821104_dp,2.0e-15_dp,'negative-binomial PMF')
  call check_close(dnegbin(4,2.5_dp,mu=3.75_dp),0.118338754598821104_dp,2.0e-15_dp,'mu parameterization')
  call check_close(dnegbin(4,p=0.4_dp,mu=3.75_dp),0.118338754598821104_dp,2.0e-15_dp,'missing-nu parameterization')
  call negbin_from_nu_p(2.5_dp,0.4_dp,mu,sigma2)
  call check_close(mu,3.75_dp,1.0e-14_dp,'mean')
  call check_close(sigma2,9.375_dp,1.0e-14_dp,'variance')
  call negbin_from_mu_sigma2(mu,sigma2,nu,p)
  call check_close(nu,2.5_dp,2.0e-14_dp,'nu round trip')
  call check_close(p,0.4_dp,2.0e-14_dp,'p round trip')
  call negbin_from_mu_nu(mu,nu,sigma2,p)
  call check_close(p,0.4_dp,2.0e-14_dp,'mu-nu conversion')
  call negbin_from_mu_p(mu,p,sigma2,nu)
  call check_close(nu,2.5_dp,2.0e-14_dp,'mu-p conversion')
  call negbin_from_sigma2_p(sigma2,p,mu,nu)
  call check_close(mu,3.75_dp,2.0e-14_dp,'sigma2-p conversion')
  if(failures==0) then
    print '(a)','test_negbin: PASS'
  else
    print '(a,i0)','test_negbin: FAIL ',failures
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
end program test_negbin
