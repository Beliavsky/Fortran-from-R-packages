! SPDX-License-Identifier: GPL-2.0-or-later
program test_supplements
  use actuar, only : dp, mexp, levexp, mgfexp, mnorm, mgfnorm, &
    mbeta, mgamma, mgfgamma, mweibull, mlnorm, munif, mgfunif, &
    mchisq, mgfchisq
  implicit none

  call assert_close(mexp(2.0_dp,3.0_dp),18.0_dp,1.0e-14_dp)
  call assert_close(levexp(2.0_dp,3.0_dp,1.0_dp), &
    3.0_dp*(1.0_dp-exp(-2.0_dp/3.0_dp)),1.0e-13_dp)
  call assert_close(mgfexp(0.2_dp,3.0_dp),2.5_dp,1.0e-14_dp)
  call assert_close(mnorm(4,1.0_dp,2.0_dp),73.0_dp,1.0e-13_dp)
  call assert_close(mgfnorm(0.3_dp,1.0_dp,2.0_dp),exp(0.48_dp),1.0e-13_dp)
  call assert_close(mbeta(2.0_dp,2.0_dp,3.0_dp),0.2_dp,1.0e-13_dp)
  call assert_close(mgamma(2.0_dp,3.0_dp,2.0_dp),48.0_dp,1.0e-13_dp)
  call assert_close(mgfgamma(0.1_dp,3.0_dp,2.0_dp),1.953125_dp,1.0e-13_dp)
  call assert_close(mweibull(2.0_dp,2.0_dp,3.0_dp),9.0_dp,1.0e-13_dp)
  call assert_close(mlnorm(2.0_dp,0.5_dp,0.7_dp),exp(1.98_dp),1.0e-13_dp)
  call assert_close(munif(1.0_dp,2.0_dp,6.0_dp),4.0_dp,1.0e-14_dp)
  call assert_close(mgfunif(0.0_dp,2.0_dp,6.0_dp),1.0_dp,1.0e-14_dp)
  call assert_close(mchisq(2,4.0_dp,0.0_dp),24.0_dp,1.0e-13_dp)
  call assert_close(mgfchisq(0.1_dp,4.0_dp,0.0_dp),1.5625_dp,1.0e-13_dp)

  print '(a)', 'test_supplements: PASS'
contains
  subroutine assert_close(actual,expected,tol)
    real(dp), intent(in) :: actual,expected,tol
    if(abs(actual-expected)>tol*max(1.0_dp,abs(expected))) then
      print '(a,3es24.15)', 'mismatch: ',actual,expected,abs(actual-expected)
      error stop 1
    end if
  end subroutine assert_close
end program test_supplements
