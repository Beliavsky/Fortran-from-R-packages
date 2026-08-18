! SPDX-License-Identifier: GPL-2.0-or-later
program test_quantiles
   use lindley_power_series
   implicit none

   real(dp), parameter :: lambda = 1.3_dp, theta = 0.4_dp, p = 0.37_dp
   integer, parameter :: m = 3
   integer :: failures
   real(dp) :: q

   failures = 0

   q = qlindleygeometric(p,lambda,theta)
   call check(q, 0.81707164692440651957_dp, 5.0e-13_dp, 'geometric q')
   call check(plindleygeometric(q,lambda,theta), p, 5.0e-13_dp, 'geometric inversion')

   q = qlindleylogarithmic(p,lambda,theta)
   call check(q, 0.68482540557271014137_dp, 5.0e-13_dp, 'logarithmic q')
   call check(plindleylogarithmic(q,lambda,theta), p, 5.0e-13_dp, &
      'logarithmic inversion')

   q = qlindleynb(p,lambda,theta,m)
   call check(q, 1.81804577374051083112_dp, 8.0e-13_dp, 'negative-binomial q')
   call check(plindleynb(q,lambda,theta,m), p, 8.0e-13_dp, &
      'negative-binomial inversion')

   q = qlindleybinomial(p,lambda,theta,m)
   call check(q, 0.72635742429819548869_dp, 5.0e-13_dp, 'binomial q')
   call check(plindleybinomial(q,lambda,theta,m), p, 5.0e-13_dp, &
      'binomial inversion')

   q = qlindleypoisson(p,lambda,theta)
   call check(q, 0.66029074536608947558_dp, 5.0e-13_dp, 'poisson q')
   call check(plindleypoisson(q,lambda,theta), p, 5.0e-13_dp, 'poisson inversion')

   ! The binomial power-series parameter theta is positive, not restricted to < 1.
   q = qlindleybinomial(0.6_dp, 0.8_dp, 2.5_dp, 4)
   call check(plindleybinomial(q,0.8_dp,2.5_dp,4), 0.6_dp, 2.0e-12_dp, &
      'binomial theta > 1')

   if (failures /= 0) error stop 1
   print '(a)', 'test_quantiles: PASS'

contains

   subroutine check(got, expected, tol, name)
      real(dp), intent(in) :: got, expected, tol
      character(len=*), intent(in) :: name

      if (abs(got - expected) > tol) then
         print '(a,2es24.15)', trim(name)//': FAIL ', got, expected
         failures = failures + 1
      end if
   end subroutine check

end program test_quantiles
