! SPDX-License-Identifier: GPL-3.0-or-later
program test_reference
   use splines2
   implicit none
   type(spline_basis) :: ans
   real(dp), allocatable :: x(:), knots(:), expected(:, :)
   real(dp) :: bounds(2), err, worst
   integer :: family, degree, derivs, integ, intercept, periodic, scale, nr, nk, nc
   integer :: unit, ios, i, count
   character(256) :: header
   open(newunit=unit, file='test/reference.dat', status='old', action='read')
   read(unit, '(a)') header
   read(unit, '(a)') header
   count = 0
   worst = 0.0_dp
   do
      read(unit, *, iostat=ios) family, degree, derivs, integ, intercept, periodic, scale
      if (ios < 0) exit
      if (ios /= 0) error stop 'Malformed reference header'
      read(unit, *) nr, nk, nc
      read(unit, *) bounds
      allocate(x(nr), knots(nk), expected(nr, nc))
      read(unit, *) x
      if (nk > 0) read(unit, *) knots
      do i = 1, nc
         read(unit, *) expected(:, i)
      end do
      select case (family)
      case (1)
         ans = b_spline(x, knots=knots, degree=degree, derivs=derivs, integral=integ /= 0, &
            intercept=intercept /= 0, boundary_knots=bounds, periodic=periodic /= 0)
      case (2)
         ans = m_spline(x, knots=knots, degree=degree, derivs=derivs, integral=integ /= 0, &
            intercept=intercept /= 0, boundary_knots=bounds, periodic=periodic /= 0)
      case (3)
         ans = i_spline(x, knots=knots, degree=degree, derivs=derivs, &
            intercept=intercept /= 0, boundary_knots=bounds)
      case (4)
         ans = c_spline(x, knots=knots, degree=degree, derivs=derivs, scale=scale /= 0, &
            intercept=intercept /= 0, boundary_knots=bounds)
      case (5)
         ans = natural_spline(x, knots=knots, derivs=derivs, integral=integ /= 0, &
            intercept=intercept /= 0, boundary_knots=bounds)
      case (6)
         ans = nsk(x, knots=knots, derivs=derivs, integral=integ /= 0, &
            intercept=intercept /= 0, boundary_knots=bounds)
      case (7)
         ans = bernstein_poly(x, degree=degree, derivs=derivs, integral=integ /= 0, &
            intercept=intercept /= 0, boundary_knots=bounds)
      end select
      count = count + 1
      if (ans%status /= SPLINES2_OK) then
         print *, 'Failed reference case', count, ans%message
         error stop 'Unexpected error status'
      end if
      if (any(shape(ans%values) /= shape(expected))) error stop 'Incorrect reference shape'
      err = maxval(abs(ans%values - expected) / max(1.0_dp, abs(expected)))
      worst = max(worst, err)
      if (err > 2.0e-10_dp) then
         print *, 'Case, family, degree, derivative, integral, intercept, periodic, scale:'
         print *, count, family, degree, derivs, integ, intercept, periodic, scale
         print *, 'Normalized error:', err
         print *, 'Worst location:', maxloc(abs(ans%values - expected) / max(1.0_dp, abs(expected)))
         error stop 'Upstream reference mismatch'
      end if
      deallocate(x, knots, expected)
   end do
   close(unit)
   print *, count, 'upstream cases passed; worst normalized error:', worst
end program test_reference
