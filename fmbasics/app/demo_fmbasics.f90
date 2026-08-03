! SPDX-License-Identifier: GPL-2.0-only
program demo_fmbasics
   use fmbasics
   implicit none
   type(currency_pair_t) :: pair
   type(interest_rate_t) :: rate
   type(discount_factor_t) :: df
   integer :: trade_date, status

   pair = eurusd()
   trade_date = make_date(2026, 7, 30)
   write(*,'(a,1x,a)') 'EURUSD spot date:', date_string(to_spot(trade_date, pair))

   rate = interest_rate(0.05_dp, COMPOUND_CONTINUOUS, 'act/365', status)
   df = as_discount_factor(rate, trade_date, add_months(trade_date, 12, .false.), status)
   write(*,'(a,f12.8)') 'One-year discount factor at 5%: ', df%value(1)
end program demo_fmbasics
