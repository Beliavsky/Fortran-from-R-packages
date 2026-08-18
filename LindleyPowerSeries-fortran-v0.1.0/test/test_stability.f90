! SPDX-License-Identifier: GPL-2.0-or-later
program test_stability
   use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
   use lindley_power_series
   implicit none

   real(dp) :: x, q, p
   integer :: failures

   failures = 0

   ! Far-tail hazards should remain finite even after the CDF rounds to one.
   x = 600.0_dp
   call finite_near(hlindleygeometric(x,1.3_dp,0.4_dp), 1.3_dp, 0.01_dp, &
      'geometric tail hazard')
   call finite_near(hlindleylogarithmic(x,1.3_dp,0.4_dp), 1.3_dp, 0.01_dp, &
      'logarithmic tail hazard')
   call finite_near(hlindleynb(x,1.3_dp,0.4_dp,3), 1.3_dp, 0.01_dp, &
      'negative-binomial tail hazard')
   call finite_near(hlindleybinomial(x,1.3_dp,2.0_dp,4), 1.3_dp, 0.01_dp, &
      'binomial tail hazard')
   call finite_near(hlindleypoisson(x,1.3_dp,25.0_dp), 1.3_dp, 0.01_dp, &
      'poisson tail hazard')

   ! Large power-series parameters must not overflow in quantile inversion.
   p = 0.8_dp
   q = qlindleypoisson(p, 1.1_dp, 1000.0_dp)
   call finite_near(plindleypoisson(q,1.1_dp,1000.0_dp), p, 2.0e-11_dp, &
      'large-theta poisson inversion')
   q = qlindleybinomial(p, 1.1_dp, 20.0_dp, 80)
   call finite_near(plindleybinomial(q,1.1_dp,20.0_dp,80), p, 2.0e-11_dp, &
      'large-parameter binomial inversion')

   if (failures /= 0) error stop 1
   print '(a)', 'test_stability: PASS'

contains

   subroutine finite_near(got, expected, tol, name)
      real(dp), intent(in) :: got, expected, tol
      character(len=*), intent(in) :: name

      if (ieee_is_nan(got) .or. abs(got-expected) > tol) then
         print '(a,2es24.15)', trim(name)//': FAIL ', got, expected
         failures = failures + 1
      end if
   end subroutine finite_near

end program test_stability
