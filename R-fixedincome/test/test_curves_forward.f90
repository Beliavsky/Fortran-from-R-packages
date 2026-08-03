! SPDX-License-Identifier: MIT
program test_curves_forward
   use fixedincome
   implicit none
   type(spot_rate_t) :: spot
   type(spot_rate_curve_t) :: curve, recovered, selected
   type(forward_rate_t) :: fwd, between
   type(term_t) :: terms
   real(dp), allocatable :: values(:), factors(:)
   integer :: status

   spot = spotrate(0.06_dp, 'simple', 'actual/365', 'actual', status)
   values = spot_compound(spot, [10.0_dp], status)
   call check_close(values(1), 1.0016438356164384_dp, 1.0e-13_dp, 'spot compound')
   values = spot_discount(spot, [10.0_dp], status)
   call check_close(values(1), 0.9983588621444204_dp, 1.0e-13_dp, 'spot discount')

   terms = term([1.0_dp, 11.0_dp, 26.0_dp, 27.0_dp, 28.0_dp], 'days')
   curve = spotratecurve([0.0719_dp,0.056_dp,0.0674_dp,0.0687_dp,0.07_dp], terms, &
      'discrete', 'actual/365', 'actual', gregorian_to_ordinal(2026,7,30), status)
   call check(status == FI_OK .and. curve%size() == 5, 'curve construction')
   factors = curve_compound(curve, status)
   call check(all(factors > 1.0_dp), 'curve factors')

   fwd = forwardrate_from_curve(curve, status)
   call check(status == FI_OK .and. fwd%size() == curve%size(), 'forward curve')
   recovered = as_spotratecurve(fwd, status=status)
   call check(maxval(abs(recovered%rate-curve%rate)) < 5.0e-13_dp, 'forward roundtrip')
   call check(maxval(abs(recovered%term_days-curve%term_days)) < 1.0e-14_dp, 'term roundtrip')

   between = forwardrate_between(curve, 11.0_dp, 26.0_dp, status)
   call check(status == FI_OK .and. between%size() == 1, 'between forward')
   call check(abs(between%interval_days(1)-15.0_dp) < 1.0e-14_dp, 'between term')

   selected = first(curve, term(11.0_dp, 'days'), status)
   call check(selected%size() == 2, 'first helper')
   selected = last(curve, term(2.0_dp, 'days'), status)
   call check(selected%size() == 3, 'last helper')
   selected = closest(curve, term(25.0_dp, 'days'), status)
   call check(selected%size() == 1 .and. abs(selected%term_days(1)-26.0_dp) < 1.0e-14_dp, 'closest helper')

   call insert_curve_points(curve, [10.0_dp, 11.0_dp], [0.05_dp, 0.051_dp], status)
   call check(curve%size() == 6, 'insert new and replace old')
   selected = curve_at_terms(curve, [10.0_dp, 11.0_dp], status)
   call check(maxval(abs(selected%rate-[0.05_dp,0.051_dp])) < 1.0e-14_dp, 'insert values')

   print '(a)', 'test_curves_forward: PASS'
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
end program test_curves_forward
