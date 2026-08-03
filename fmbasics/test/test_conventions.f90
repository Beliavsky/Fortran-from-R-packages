! SPDX-License-Identifier: GPL-2.0-only
program test_conventions
   use fmbasics
   implicit none
   type(currency_pair_t) :: pair
   type(index_t) :: index
   integer :: dates(2), values(2)

   pair = audusd()
   call check(pair_iso(pair) == 'AUDUSD', 'pair iso')
   call check(.not. is_t1(pair), 'AUDUSD is T+2')
   call check(pair_iso(invert(pair)) == 'USDAUD', 'pair invert')

   dates = [make_date(2014,4,16), make_date(2014,4,19)]
   values = to_spot(dates, pair)
   call check(all(values == [make_date(2014,4,22), make_date(2014,4,23)]), 'spot dates')
   values = to_fx_value(dates, 'tomorrow', pair)
   call check(values(1) == make_date(2014,4,17), 'tomorrow weekday')
   call check(values(2) == make_date(2014,4,22), 'tomorrow holiday')
   values = to_fx_value(dates, months_period(1), pair)
   call check(all(values == [make_date(2014,5,22), make_date(2014,5,23)]), 'one month forward')

   index = audbbsw(months_period(3))
   call check(to_reset(make_date(2017,1,3), index) == make_date(2017,1,3), 'index reset')
   call check(to_value(make_date(2017,1,3), index) == make_date(2017,1,3), 'index value')
   call check(to_maturity(make_date(2017,1,3), index) == make_date(2017,4,3), 'index maturity')

   index = usdlibor(months_period(3))
   call check(index%spot_lag%days == 2, 'LIBOR spot lag')
   call check(trim(index%day_basis) == 'act/360', 'LIBOR day basis')

   print '(a)', 'test_conventions: PASS'

contains

   subroutine check(condition, label)
      logical, intent(in) :: condition
      character(len=*), intent(in) :: label
      if (.not. condition) then
         write(*,'(a,1x,a)') 'FAIL:', trim(label)
         error stop 1
      end if
   end subroutine check

end program test_conventions
