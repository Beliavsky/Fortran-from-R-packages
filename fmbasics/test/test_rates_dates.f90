! SPDX-License-Identifier: GPL-2.0-only
program test_rates_dates
   use fmbasics
   implicit none
   type(interest_rate_t) :: rate, converted
   type(discount_factor_t) :: df, df1, df2, df3
   integer :: d1, d2, status

   d1 = make_date(2010, 1, 1)
   d2 = make_date(2015, 1, 1)
   call check(abs(year_frac(d1, d2, 'act/360') - 1826.0_dp/360.0_dp) < 1.0e-14_dp, 'act/360')

   rate = interest_rate(0.04_dp, 0.0_dp, 'act/360', status)
   call check(status == FM_OK, 'rate constructor')
   df = as_discount_factor(rate, d1, d2, status)
   call check(abs(df%value(1) - 0.831331978570109_dp) < 1.0e-14_dp, 'simple discount')

   rate = interest_rate(0.075_dp, COMPOUND_CONTINUOUS, 'act/365')
   df = as_discount_factor(rate, d1, d2)
   call check(abs(df%value(1) - 0.687148069474866_dp) < 1.0e-14_dp, 'continuous discount')

   rate = interest_rate(0.04_dp, 0.0_dp, 'act/360')
   converted = convert_interest_rate(rate, COMPOUND_CONTINUOUS, 'act/360')
   converted = convert_interest_rate(converted, 0.0_dp, 'act/360')
   call check(abs(converted%value(1) - 0.04_dp) < 1.0e-13_dp, 'rate conversion round trip')

   df = discount_factor(0.95_dp, d1, d2)
   rate = as_interest_rate(df, 2.0_dp, 'act/365')
   call check(abs(rate%value(1) - 0.0102793669522563_dp) < 1.0e-13_dp, 'df to rate')

   df1 = discount_factor(0.99_dp, make_date(2013,1,1), make_date(2014,1,1))
   df2 = discount_factor(0.98_dp, make_date(2014,1,1), make_date(2015,1,1))
   df3 = discount_multiply(df1, df2, status)
   call check(status == FM_OK, 'df multiplication status')
   call check(abs(df3%value(1) - 0.9702_dp) < 1.0e-15_dp, 'df multiplication')
   call check(df3%start_date(1) == make_date(2013,1,1), 'df multiplication start')
   call check(df3%end_date(1) == make_date(2015,1,1), 'df multiplication end')

   df3 = discount_divide(df3, df1, status)
   call check(status == FM_OK .and. abs(df3%value(1)-0.98_dp) < 1.0e-14_dp, 'df division')

   print '(a)', 'test_rates_dates: PASS'

contains

   subroutine check(condition, label)
      logical, intent(in) :: condition
      character(len=*), intent(in) :: label
      if (.not. condition) then
         write(*,'(a,1x,a)') 'FAIL:', trim(label)
         error stop 1
      end if
   end subroutine check

end program test_rates_dates
