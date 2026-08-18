! SPDX-License-Identifier: GPL-2.0-or-later
! Modern Fortran translation of RMKdiscrete 0.1 by Robert M. Kirkpatrick.
! See LICENSE, COPYING, and TRANSLATION_NOTES.md for provenance.
program test_lgp
  use rmkdiscrete, only : dp, dlgp, plgp, qlgp, lgp_get_nc, lgp_findmax, lgp_summary, slgp, &
    lgp_from_mu_sigma2, lgp_from_theta_lambda, lgp_from_mu_theta, lgp_from_sigma2_lambda, &
    lgp_from_sigma2_theta, lgp_from_mu_lambda
  implicit none
  real(dp) :: nc, s, theta, lambda, mu, sigma2
  type(lgp_summary) :: z
  integer :: x, failures
  failures=0

  call check_close(dlgp(3,2.0_dp,0.3_dp),0.154248426891461615_dp,2.0e-14_dp,'positive-lambda PMF')
  call check_close(plgp(3.0_dp,2.0_dp,0.3_dp),0.683212700930949889_dp,3.0e-14_dp,'positive-lambda CDF')
  call check_close(qlgp(0.5_dp,2.0_dp,0.3_dp),2.0_dp,0.0_dp,'quantile')

  nc=lgp_get_nc(2.0_dp,-0.2_dp)
  call check_close(nc,0.999999999999982145_dp,2.0e-14_dp,'negative-lambda normalizer')
  call check_close(lgp_findmax(2.0_dp,-0.2_dp),9.0_dp,0.0_dp,'negative-lambda support')
  s=0.0_dp
  do x=0,9
    s=s+dlgp(x,2.0_dp,-0.2_dp,nc)
  end do
  call check_close(s,1.0_dp,2.0e-14_dp,'finite-support normalization')
  call check_close(dlgp(4,2.0_dp,-0.2_dp,nc),0.0433719665153578722_dp,2.0e-14_dp,'negative-lambda PMF')

  call lgp_from_mu_sigma2(4.0_dp,8.0_dp,theta,lambda)
  call lgp_from_theta_lambda(theta,lambda,mu,sigma2)
  call check_close(mu,4.0_dp,2.0e-14_dp,'MVP round trip mean')
  call check_close(sigma2,8.0_dp,2.0e-14_dp,'MVP round trip variance')
  call lgp_from_mu_theta(mu,theta,sigma2,lambda)
  call check_close(sigma2,8.0_dp,2.0e-14_dp,'mu-theta conversion')
  call lgp_from_sigma2_lambda(sigma2,lambda,mu,theta)
  call check_close(mu,4.0_dp,2.0e-14_dp,'sigma2-lambda conversion')
  call lgp_from_sigma2_theta(sigma2,theta,mu,lambda)
  call check_close(mu,4.0_dp,2.0e-14_dp,'sigma2-theta conversion')
  call lgp_from_mu_lambda(mu,lambda,sigma2,theta)
  call check_close(sigma2,8.0_dp,2.0e-14_dp,'mu-lambda conversion')

  z=slgp(2.0_dp,0.0_dp)
  call check_close(z%mean,2.0_dp,1.0e-14_dp,'Poisson mean')
  call check_close(z%variance,2.0_dp,1.0e-14_dp,'Poisson variance')
  call check_close(z%third_central_moment,2.0_dp,1.0e-14_dp,'Poisson third moment')
  call check_close(z%fourth_central_moment,14.0_dp,1.0e-13_dp,'Poisson fourth moment')

  z=slgp(2.0_dp,-0.2_dp,do_numerically=.true.)
  call check_close(z%mean,sum([(real(x,dp)*dlgp(x,2.0_dp,-0.2_dp,nc),x=0,9)]),2.0e-14_dp, &
    'numeric negative-lambda mean')

  if(failures==0) then
    print '(a)','test_lgp: PASS'
  else
    print '(a,i0)','test_lgp: FAIL ',failures
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
end program test_lgp
