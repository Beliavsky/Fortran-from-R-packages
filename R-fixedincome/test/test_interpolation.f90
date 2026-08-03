! SPDX-License-Identifier: MIT
program test_interpolation
   use fixedincome
   implicit none
   type(spot_rate_curve_t) :: curve
   type(term_t) :: terms
   real(dp), allocatable :: values(:)
   real(dp) :: midpoint
   integer :: status, method

   terms = term([1.0_dp, 11.0_dp, 26.0_dp, 40.0_dp, 60.0_dp], 'days')
   curve = spotratecurve([0.0719_dp,0.056_dp,0.0674_dp,0.0687_dp,0.07_dp], terms, &
      'discrete', 'actual/365', 'actual', status=status)

   call set_interpolation(curve, interp_linear(), status)
   values = interpolate(curve, [6.0_dp, 11.0_dp], status)
   midpoint = 0.0719_dp + 0.5_dp*(0.056_dp-0.0719_dp)
   call check_close(values(1), midpoint, 1.0e-14_dp, 'linear midpoint')
   call check_close(values(2), 0.056_dp, 1.0e-14_dp, 'linear knot')

   call set_interpolation(curve, interp_loglinear(), status)
   values = interpolate(curve, [6.0_dp], status)
   call check_close(values(1), sqrt(0.0719_dp*0.056_dp), 1.0e-14_dp, 'loglinear midpoint')

   do method = INTERP_NATURAL_SPLINE, INTERP_MONOTONE_SPLINE
      select case (method)
      case (INTERP_NATURAL_SPLINE)
         call set_interpolation(curve, interp_naturalspline(), status)
      case (INTERP_HERMITE_SPLINE)
         call set_interpolation(curve, interp_hermitespline(), status)
      case (INTERP_MONOTONE_SPLINE)
         call set_interpolation(curve, interp_monotonespline(), status)
      end select
      values = interpolate(curve, curve%term_days, status)
      call check(maxval(abs(values-curve%rate)) < 5.0e-13_dp, 'spline knots')
   end do

   call set_interpolation(curve, interp_flatforward(), status)
   values = interpolate(curve, [6.0_dp], status)
   call check(values(1) > 0.0_dp .and. values(1) < 0.2_dp, 'flat forward finite')

   print '(a)', 'test_interpolation: PASS'
contains
   subroutine check(condition, label)
      logical, intent(in) :: condition
      character(len=*), intent(in) :: label
      if (.not. condition) then
         write(*,'(a)') 'FAIL: '//trim(label)
         error stop 1
      end if
   end subroutine check
   subroutine check_close(actual, expected, tolerance, label)
      real(dp), intent(in) :: actual, expected, tolerance
      character(len=*), intent(in) :: label
      call check(abs(actual-expected) <= tolerance, label)
   end subroutine check_close
end program test_interpolation
