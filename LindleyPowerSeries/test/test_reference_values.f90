! SPDX-License-Identifier: GPL-2.0-or-later
program test_reference_values
   use lindley_power_series
   implicit none

   real(dp), parameter :: lambda = 1.3_dp, theta = 0.4_dp, x = 1.2_dp
   integer, parameter :: m = 3
   integer :: failures

   failures = 0
   call check(plindleygeometric(x,lambda,theta), &
      0.52411323879667105763_dp, 3.0e-13_dp, 'geometric cdf')
   call check(dlindleygeometric(x,lambda,theta), &
      0.37112524792811908726_dp, 3.0e-13_dp, 'geometric pdf')
   call check(hlindleygeometric(x,lambda,theta), &
      0.77986041677160860024_dp, 5.0e-13_dp, 'geometric hazard')

   call check(plindleylogarithmic(x,lambda,theta), &
      0.58663186776897182213_dp, 3.0e-13_dp, 'logarithmic cdf')
   call check(dlindleylogarithmic(x,lambda,theta), &
      0.35893269626886535464_dp, 3.0e-13_dp, 'logarithmic pdf')
   call check(hlindleylogarithmic(x,lambda,theta), &
      0.86831245149846411273_dp, 5.0e-13_dp, 'logarithmic hazard')

   call check(plindleynb(x,lambda,theta,m), &
      0.14397112212675173330_dp, 3.0e-13_dp, 'negative-binomial cdf')
   call check(dlindleynb(x,lambda,theta,m), &
      0.30583840154346214691_dp, 3.0e-13_dp, 'negative-binomial pdf')
   call check(hlindleynb(x,lambda,theta,m), &
      0.35727579927361688111_dp, 5.0e-13_dp, 'negative-binomial hazard')

   call check(plindleybinomial(x,lambda,theta,m), &
      0.57070336945693574938_dp, 3.0e-13_dp, 'binomial cdf')
   call check(dlindleybinomial(x,lambda,theta,m), &
      0.37044467832307869118_dp, 3.0e-13_dp, 'binomial pdf')
   call check(hlindleybinomial(x,lambda,theta,m), &
      0.86291075207010758565_dp, 5.0e-13_dp, 'binomial hazard')

   call check(plindleypoisson(x,lambda,theta), &
      0.60092398416058718766_dp, 3.0e-13_dp, 'poisson cdf')
   call check(dlindleypoisson(x,lambda,theta), &
      0.35791981713409199179_dp, 3.0e-13_dp, 'poisson pdf')
   call check(hlindleypoisson(x,lambda,theta), &
      0.89687127997719118337_dp, 5.0e-13_dp, 'poisson hazard')

   if (failures /= 0) error stop 1
   print '(a)', 'test_reference_values: PASS'

contains

   subroutine check(got, expected, tol, name)
      real(dp), intent(in) :: got, expected, tol
      character(len=*), intent(in) :: name

      if (abs(got - expected) > tol) then
         print '(a,2es24.15)', trim(name)//': FAIL ', got, expected
         failures = failures + 1
      end if
   end subroutine check

end program test_reference_values
