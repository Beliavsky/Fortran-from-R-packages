! SPDX-License-Identifier: GPL-2.0-or-later
program test_identities
   use lindley_power_series
   implicit none

   real(dp), parameter :: lambda = 0.9_dp, theta = 0.55_dp, x = 1.7_dp
   integer, parameter :: m = 5
   integer :: failures

   failures = 0
   call hazard_check(dlindleygeometric(x,lambda,theta), &
      plindleygeometric(x,lambda,theta), hlindleygeometric(x,lambda,theta), &
      'geometric')
   call hazard_check(dlindleylogarithmic(x,lambda,theta), &
      plindleylogarithmic(x,lambda,theta), hlindleylogarithmic(x,lambda,theta), &
      'logarithmic')
   call hazard_check(dlindleynb(x,lambda,theta,m), plindleynb(x,lambda,theta,m), &
      hlindleynb(x,lambda,theta,m), 'negative-binomial')
   call hazard_check(dlindleybinomial(x,lambda,theta,m), &
      plindleybinomial(x,lambda,theta,m), hlindleybinomial(x,lambda,theta,m), &
      'binomial')
   call hazard_check(dlindleypoisson(x,lambda,theta), &
      plindleypoisson(x,lambda,theta), hlindleypoisson(x,lambda,theta), 'poisson')

   call check(plindleygeometric(x,lambda,theta,.true.), &
      log(plindleygeometric(x,lambda,theta)), 5.0e-14_dp, 'geometric log.p')
   call check(plindleypoisson(x,lambda,theta,.true.), &
      log(plindleypoisson(x,lambda,theta)), 5.0e-14_dp, 'poisson log.p')

   if (failures /= 0) error stop 1
   print '(a)', 'test_identities: PASS'

contains

   subroutine hazard_check(pdf, cdf, hazard, name)
      real(dp), intent(in) :: pdf, cdf, hazard
      character(len=*), intent(in) :: name

      call check(hazard, pdf/(1.0_dp-cdf), 2.0e-12_dp, trim(name)//' hazard')
   end subroutine hazard_check

   subroutine check(got, expected, tol, name)
      real(dp), intent(in) :: got, expected, tol
      character(len=*), intent(in) :: name

      if (abs(got - expected) > tol) then
         print '(a,2es24.15)', trim(name)//': FAIL ', got, expected
         failures = failures + 1
      end if
   end subroutine check

end program test_identities
