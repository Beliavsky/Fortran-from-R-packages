! SPDX-License-Identifier: GPL-2.0-or-later
program test_densities
   use rnd, only : dp, regularized_beta, dgb, pgb, dmln, dmln_am, dew, dshimko
   use rnd, only : lognormal_pdf
   implicit none
   real(dp) :: x, sigma, te, s0, r, y, meanlog, value, skew, kurt, v

   call assert_close(regularized_beta(0.5_dp,2.0_dp,2.0_dp),0.5_dp,2.0e-13_dp,"beta cdf")
   x = 80.0_dp
   call assert_true(dgb(x,2.0_dp,100.0_dp,3.0_dp,4.0_dp) > 0.0_dp,"dgb positive")
   call assert_true(pgb(50.0_dp,2.0_dp,100.0_dp,3.0_dp,4.0_dp) < &
      pgb(150.0_dp,2.0_dp,100.0_dp,3.0_dp,4.0_dp),"pgb monotone")

   meanlog = log(100.0_dp)
   value = lognormal_pdf(x,meanlog,0.2_dp)
   call assert_close(dmln(x,1.0_dp,meanlog,meanlog+0.1_dp,0.2_dp,0.3_dp), &
      value,1.0e-14_dp,"dmln one component")
   call assert_close(dmln_am(x,meanlog,meanlog,meanlog,0.2_dp,0.2_dp,0.2_dp, &
      0.2_dp,0.3_dp),value,1.0e-14_dp,"dmln_am common components")

   sigma = 0.25_dp
   te = 0.5_dp
   s0 = 100.0_dp
   r = 0.04_dp
   y = 0.01_dp
   meanlog = log(s0)+(r-y-0.5_dp*sigma*sigma)*te
   v = sqrt(exp(sigma*sigma*te)-1.0_dp)
   skew = 3.0_dp*v+v**3
   kurt = 16.0_dp*v**2+15.0_dp*v**4+6.0_dp*v**6+v**8
   call assert_close(dew(x,r,y,te,s0,sigma,skew,kurt), &
      lognormal_pdf(x,meanlog,sigma*sqrt(te)),2.0e-13_dp,"dew baseline")
   call assert_close(dshimko(r,te,s0,x,y,sigma,0.0_dp,0.0_dp), &
      lognormal_pdf(x,meanlog,sigma*sqrt(te)),2.0e-13_dp,"dshimko constant vol")
   print '(a)', 'test_densities: PASS'
contains
   subroutine assert_close(actual,expected,tolerance,message)
      real(dp), intent(in) :: actual, expected, tolerance
      character(*), intent(in) :: message
      if (abs(actual-expected) > tolerance) then
         print '(a,2es20.10)', trim(message)//': ',actual,expected
         error stop 1
      end if
   end subroutine assert_close
   subroutine assert_true(condition,message)
      logical, intent(in) :: condition
      character(*), intent(in) :: message
      if (.not. condition) then
         print '(a)', trim(message)
         error stop 1
      end if
   end subroutine assert_true
end program test_densities
