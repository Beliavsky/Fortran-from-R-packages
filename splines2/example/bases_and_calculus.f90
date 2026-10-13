! SPDX-License-Identifier: GPL-3.0-or-later
program bases_and_calculus
   use splines2, only : dp, spline_basis, b_spline, i_spline, c_spline, nsk, derivative_spline
   implicit none
   type(spline_basis) :: b, slope, monotone, convex, cardinal
   real(dp) :: x(5), knots(3)
   x = [0.0_dp, 0.25_dp, 0.5_dp, 0.75_dp, 1.0_dp]
   knots = [0.25_dp, 0.5_dp, 0.75_dp]
   b = b_spline(x, knots=knots, intercept=.true.)
   if (b%status /= 0) error stop b%message
   slope = derivative_spline(b)
   monotone = i_spline(x, knots=knots)
   convex = c_spline(x, knots=knots)
   cardinal = nsk(x, knots=knots, intercept=.true.)
   print '(a,5f8.4)', 'B-spline row sums: ', sum(b%values, dim=2)
   print '(a,5f8.4)', 'Derivative row sums: ', sum(slope%values, dim=2)
   print '(a,7f8.4)', 'I-splines at right boundary: ', monotone%values(5, :)
   print '(a,7f8.4)', 'C-splines at right boundary: ', convex%values(5, :)
   print '(a,5f8.4)', 'nsk at middle knot: ', cardinal%values(3, :)
end program bases_and_calculus
