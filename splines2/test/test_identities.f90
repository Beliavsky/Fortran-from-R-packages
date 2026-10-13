! SPDX-License-Identifier: GPL-3.0-or-later
program test_identities
   use splines2
   use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan, ieee_positive_inf, ieee_is_nan
   implicit none
   type(spline_basis) :: b, m, ii, c, ns, z, other, deriv, plus, minus
   real(dp), parameter :: knots(3) = [0.2_dp, 0.5_dp, 0.9_dp], bounds(2) = [0.0_dp, 1.0_dp]
   real(dp) :: grid(101), interior(3), h, nan, inf, ident(5, 5), xmissing(7), qx(9)
   integer :: j, degree, checks
   checks = 0
   h = 1.0e-5_dp
   nan = ieee_value(0.0_dp, ieee_quiet_nan)
   inf = ieee_value(0.0_dp, ieee_positive_inf)
   grid = [(real(j, dp) / 100.0_dp, j = 0, 100)]
   interior = [0.12_dp, 0.38_dp, 0.76_dp]
   do degree = 0, 6
      b = b_spline(grid, knots=knots, degree=degree, intercept=.true., boundary_knots=bounds)
      call check(b%status == SPLINES2_OK, 'B basis status')
      call check(maxval(abs(sum(b%values, dim=2) - 1.0_dp)) < 1.0e-13_dp, 'Partition of unity')
      call check(minval(b%values) >= -1.0e-14_dp, 'Nonnegative B basis')
      other = b_spline(grid, knots=knots, degree=degree, intercept=.false., boundary_knots=bounds)
      call close_matrix(b%values(:, 2:), other%values, 0.0_dp, 'Intercept column removal')
      deriv = derivative_spline(b)
      call check(maxval(abs(sum(deriv%values, dim=2))) < 1.0e-10_dp, 'Derivative partition')
      m = m_spline(grid, knots=knots, degree=degree, intercept=.true., boundary_knots=bounds, integral=.true.)
      call check(maxval(abs(m%values(1, :))) < 1.0e-13_dp, 'M integral at left boundary')
      call check(maxval(abs(m%values(101, :) - 1.0_dp)) < 1.0e-13_dp, 'M unit integral')
      call check(minval(m%values(2:, :) - m%values(:100, :)) >= -1.0e-13_dp, 'I monotonicity')
      ii = i_spline(grid, knots=knots, degree=degree, boundary_knots=bounds)
      call close_matrix(m%values, ii%values, 0.0_dp, 'I equals integrated M')
      c = c_spline(grid, knots=knots, degree=degree, boundary_knots=bounds)
      call check(maxval(abs(c%values(101, :) - 1.0_dp)) < 1.0e-13_dp, 'C right normalization')
      call check(maxval(abs(c%values(1, :))) < 1.0e-13_dp, 'C left zero')
      call check(minval(c%values(3:, :) - 2.0_dp * c%values(2:100, :) + c%values(:99, :)) > &
         -1.0e-12_dp, 'C convexity')
      b = b_spline(interior, knots=knots, degree=degree, intercept=.true., boundary_knots=bounds)
      plus = ibs(interior + h, knots=knots, degree=degree, intercept=.true., boundary_knots=bounds)
      minus = ibs(interior - h, knots=knots, degree=degree, intercept=.true., boundary_knots=bounds)
      call close_matrix((plus%values - minus%values) / (2.0_dp * h), b%values, 2.0e-8_dp, 'B calculus identity')
      ii = i_spline(interior, knots=knots, degree=degree, boundary_knots=bounds)
      deriv = derivative_spline(ii)
      m = m_spline(interior, knots=knots, degree=degree, intercept=.true., boundary_knots=bounds)
      call close_matrix(deriv%values, m%values, 0.0_dp, 'dI = M')
      c = c_spline(interior, knots=knots, degree=degree, boundary_knots=bounds, scale=.false.)
      deriv = derivative_spline(c)
      call close_matrix(deriv%values, ii%values, 0.0_dp, 'Unscaled dC = I')
   end do
   b = b_spline(grid, knots=knots, intercept=.true., periodic=.true., boundary_knots=bounds)
   other = predict_spline(b, grid + 7.0_dp)
   call close_matrix(b%values, other%values, 3.0e-14_dp, 'B periodicity')
   m = m_spline([0.0_dp, 1.0_dp, 2.0_dp], knots=knots, intercept=.true., periodic=.true., &
      integral=.true., boundary_knots=bounds)
   call check(maxval(abs(m%values(2, :) - 1.0_dp)) < 1.0e-13_dp, 'Periodic M unit integral')
   call check(maxval(abs(m%values(3, :) - 2.0_dp)) < 1.0e-13_dp, 'Periodic M accumulated integral')
   b = b_spline([-1.0_dp, -0.1_dp], knots=knots, periodic=.true., integral=.true., boundary_knots=bounds)
   call check(maxval(abs(b%values)) == 0.0_dp, 'Periodic integrals zero below left boundary')
   ns = natural_spline(bounds, knots=knots, intercept=.true., derivs=2, boundary_knots=bounds)
   call check(maxval(abs(ns%values)) < 1.0e-12_dp, 'Natural second derivative boundary constraint')
   ns = natural_spline([-1.0_dp, 2.0_dp], knots=knots, intercept=.true., derivs=2, boundary_knots=bounds)
   call check(maxval(abs(ns%values)) == 0.0_dp, 'Natural linear continuation')
   ns = nsk([bounds(1), knots, bounds(2)], knots=knots, intercept=.true., boundary_knots=bounds)
   ident = 0.0_dp
   do j = 1, 5
      ident(j, j) = 1.0_dp
   end do
   call close_matrix(ns%values, ident, 1.0e-13_dp, 'nsk cardinal property')
   ns = nsk(interior, knots=knots, intercept=.true., boundary_knots=bounds)
   call check(maxval(abs(sum(ns%values, dim=2) - 1.0_dp)) < 1.0e-13_dp, 'nsk partition of unity')
   plus = nsk(interior + h, knots=knots, integral=.true., boundary_knots=bounds, intercept=.true.)
   minus = nsk(interior - h, knots=knots, integral=.true., boundary_knots=bounds, intercept=.true.)
   call close_matrix((plus%values - minus%values) / (2.0_dp * h), ns%values, 3.0e-9_dp, 'nsk integral derivative')
   do degree = 0, 10
      z = bernstein_poly(grid, degree=degree, intercept=.true.)
      call check(maxval(abs(sum(z%values, dim=2) - 1.0_dp)) < 1.0e-14_dp, 'Bernstein partition')
      b = b_spline(grid, degree=degree, intercept=.true.)
      call close_matrix(z%values, b%values, 1.0e-14_dp, 'Bernstein equals knot-free B')
   end do
   xmissing = [0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp, 1.0_dp, nan]
   b = b_spline(xmissing, df=6)
   call check(b%equidistant_fallback, 'Duplicate quantiles trigger fallback')
   call check(maxval(abs(b%knots - [0.25_dp, 0.5_dp, 0.75_dp])) == 0.0_dp, 'Fallback knot locations')
   call check(all(ieee_is_nan(b%values(7, :))), 'NaN coordinate propagation')
   ii = i_spline(xmissing, df=6)
   call check(maxval(abs(ii%knots - [1.0_dp / 3.0_dp, 2.0_dp / 3.0_dp])) == 0.0_dp, 'I default intercept')
   qx = [-2.0_dp, 0.0_dp, 0.1_dp, 0.3_dp, 0.4_dp, 0.8_dp, 1.0_dp, 3.0_dp, nan]
   b = b_spline(qx, df=6, boundary_knots=bounds)
   call check(maxval(abs(b%knots - [0.15_dp, 0.35_dp, 0.7_dp])) < 1.0e-14_dp, 'Inside-boundary quantiles')
   ns = natural_spline(qx, df=4, trim=0.2_dp)
   call check(maxval(abs(ns%boundary_knots - [0.04_dp, 0.92_dp])) < 1.0e-14_dp, 'Natural trimmed boundaries')
   call check(maxval(abs(ns%knots - [0.25_dp, 0.35_dp, 0.5_dp])) < 1.0e-14_dp, 'Natural trimmed knots')
   b = bsp(grid, knots=knots)
   other = b_spline(grid, knots=knots)
   call close_matrix(b%values, other%values, 0.0_dp, 'bsp alias')
   b = msp(grid, knots=knots)
   other = m_spline(grid, knots=knots)
   call close_matrix(b%values, other%values, 0.0_dp, 'msp alias')
   b = isp(grid, knots=knots)
   other = i_spline(grid, knots=knots)
   call close_matrix(b%values, other%values, 0.0_dp, 'isp alias')
   b = csp(grid, knots=knots)
   other = c_spline(grid, knots=knots)
   call close_matrix(b%values, other%values, 0.0_dp, 'csp alias')
   b = nsp(grid, knots=knots)
   other = natural_spline(grid, knots=knots)
   call close_matrix(b%values, other%values, 0.0_dp, 'nsp alias')
   b = bpoly(grid)
   other = bernstein_poly(grid)
   call close_matrix(b%values, other%values, 0.0_dp, 'bpoly alias')
   b = dbs(grid, knots=knots)
   other = b_spline(grid, knots=knots, derivs=1)
   call close_matrix(b%values, other%values, 0.0_dp, 'dbs default derivative')
   b = ibs(grid, knots=knots)
   other = b_spline(grid, knots=knots, integral=.true.)
   call close_matrix(b%values, other%values, 0.0_dp, 'ibs integration')
   b = derivative_spline(other)
   z = b_spline(grid, knots=knots)
   call close_matrix(b%values, z%values, 0.0_dp, 'Derivative of integrated object')
   b = predict_spline(z, interior)
   other = b_spline(interior, knots=knots, boundary_knots=bounds)
   call close_matrix(b%values, other%values, 0.0_dp, 'Predict retains original boundaries')
   b = b_spline([nan, nan], boundary_knots=bounds)
   call check(b%status /= 0, 'All missing x rejected')
   b = b_spline([0.0_dp, inf], boundary_knots=bounds)
   call check(b%status /= 0, 'Infinite x rejected')
   b = b_spline(grid, knots=[0.0_dp])
   call check(b%status /= 0, 'Boundary interior knot rejected')
   b = b_spline(grid, degree=-1)
   call check(b%status /= 0, 'Negative degree rejected')
   b = b_spline(grid, derivs=-1)
   call check(b%status /= 0, 'Negative derivative rejected')
   b = b_spline(grid, df=1)
   call check(b%status /= 0, 'Too-small df rejected')
   b = b_spline(grid, periodic=.true.)
   call check(b%status /= 0, 'Too-few periodic knots rejected')
   b = b_spline(grid, degree=0)
   call check(b%status /= 0, 'Zero remaining columns rejected')
   b = natural_spline(grid, knots=[0.2_dp, 0.2_dp])
   call check(b%status /= 0, 'Repeated natural knots rejected')
   b = natural_spline(grid, trim=0.5_dp)
   call check(b%status /= 0, 'Degenerate trim rejected')
   b = bernstein_poly([-0.1_dp, 1.0_dp], boundary_knots=bounds)
   call check(b%status /= 0, 'Bernstein outside domain rejected')
   b = b_spline(grid, knots=[nan])
   call check(b%status /= 0, 'Missing knot rejected')
   b = b_spline(grid, boundary_knots=[1.0_dp, 1.0_dp])
   call check(b%status /= 0, 'Equal boundaries rejected')
   b = b_spline(grid, knots=[0.9_dp, 0.5_dp, 0.2_dp], boundary_knots=[1.0_dp, 0.0_dp])
   other = b_spline(grid, knots=knots, boundary_knots=bounds)
   call close_matrix(b%values, other%values, 0.0_dp, 'Knots and boundaries sorted')
   b = b_spline(grid, degree=0, intercept=.true., periodic=.true.)
   call check(maxval(abs(b%values - 1.0_dp)) == 0.0_dp, 'Knot-free constant periodic basis')
   b = b_spline(grid, degree=1, intercept=.true., periodic=.true.)
   call check(maxval(abs(b%values - 1.0_dp)) < 1.0e-14_dp, 'Knot-free linear periodic basis')
   m = m_spline(grid, degree=1, intercept=.true., periodic=.true., integral=.true.)
   call check(maxval(abs(m%values(:, 1) - grid)) < 1.0e-14_dp, 'Knot-free periodic M integral')
   m = m_spline(grid, knots=knots, degree=3, periodic=.true., df=2)
   call check(m%status /= 0, 'Explicit periodic M df guard retained')
   print *, checks, 'independent identity and API checks passed'
contains
   subroutine check(condition, label)
      logical, intent(in) :: condition !! True when the mathematical identity or validation condition holds.
      character(*), intent(in) :: label !! Description printed if this test fails.
      checks = checks + 1
      if (.not. condition) then
         print *, 'FAIL: ', label
         error stop 1
      end if
   end subroutine check
   subroutine close_matrix(actual, expected, tolerance, label)
      real(dp), intent(in) :: actual(:, :) !! Numerically computed matrix to be checked.
      real(dp), intent(in) :: expected(:, :) !! Independent mathematical or API-equivalent reference matrix.
      real(dp), intent(in) :: tolerance !! Maximum allowed absolute discrepancy per element.
      character(*), intent(in) :: label !! Identity being verified, printed on failure.
      call check(all(shape(actual) == shape(expected)), label // ' shape')
      call check(maxval(abs(actual - expected)) <= tolerance, label)
   end subroutine close_matrix
end program test_identities
