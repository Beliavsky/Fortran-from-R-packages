! SPDX-License-Identifier: MIT
program test_terms_compounding
   use fixedincome
   implicit none
   type(term_t) :: t, y, m, d, parsed
   type(daycount_t) :: dc
   real(dp), allocatable :: shifted(:)
   integer :: status

   dc = daycount('actual/360', status)
   call check(status == FI_OK, 'daycount status')
   call check(dib(dc) == 360, 'days in base')

   t = term([1.0_dp, 1.0_dp, 1.0_dp], 'day', status)
   y = toyears(dc, t, status)
   call check_close(y%value(1), 1.0_dp/360.0_dp, 1.0e-14_dp, 'days to years')

   t = term(1.0_dp, 'month', status)
   d = todays(dc, t, status)
   y = toyears(dc, t, status)
   call check_close(d%value(1), 30.0_dp, 1.0e-14_dp, 'month to days')
   call check_close(y%value(1), 1.0_dp/12.0_dp, 1.0e-14_dp, 'month to years')

   t = term(1.0_dp, 'year', status)
   m = tomonths(dc, t, status)
   call check_close(m%value(1), 12.0_dp, 1.0e-14_dp, 'year to months')

   parsed = parse_term('6 months', status)
   call check(status == FI_OK .and. parsed%unit(1) == UNIT_MONTH, 'parse term')
   call check_close(parsed%value(1), 6.0_dp, 1.0e-14_dp, 'parsed value')

   shifted = shift([1.0_dp, 2.0_dp, 3.0_dp], fill=0.0_dp, status=status)
   call check(all(abs(shifted-[0.0_dp,1.0_dp,2.0_dp]) < 1.0e-14_dp), 'shift')

   call check_close(compound(COMPOUND_SIMPLE, 2.0_dp, 0.05_dp), 1.1_dp, 1.0e-14_dp, 'simple compound')
   call check_close(compound(COMPOUND_DISCRETE, 2.0_dp, 0.05_dp), 1.1025_dp, 1.0e-14_dp, 'discrete compound')
   call check_close(compound(COMPOUND_CONTINUOUS, 2.0_dp, 0.05_dp), exp(0.1_dp), 1.0e-14_dp, 'continuous compound')
   call check_close(implied_rate(COMPOUND_SIMPLE, 2.0_dp, 1.1_dp), 0.05_dp, 1.0e-14_dp, 'simple implied')
   call check_close(implied_rate(COMPOUND_DISCRETE, 2.0_dp, 1.1025_dp), 0.05_dp, 1.0e-14_dp, 'discrete implied')
   call check_close(implied_rate(COMPOUND_CONTINUOUS, 2.0_dp, exp(0.1_dp)), 0.05_dp, 1.0e-14_dp, 'continuous implied')

   call check(gregorian_to_ordinal(2022,2,18)-gregorian_to_ordinal(2022,2,14) == 4, 'date ordinal')
   call check(offset_date(gregorian_to_ordinal(2022,2,18), 1, 'weekdays', status) == &
      gregorian_to_ordinal(2022,2,21), 'weekend offset')

   print '(a)', 'test_terms_compounding: PASS'
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
end program test_terms_compounding
