! SPDX-License-Identifier: GPL-2.0-only
program conventions_money_example
   use fmbasics
   implicit none
   type(currency_pair_t) :: pair
   type(index_t) :: index
   type(single_currency_money_t) :: items(3)
   type(multi_currency_money_t) :: money, totals
   integer :: trade_date, value_date, status, i

   pair = audusd()
   trade_date = make_date(2014, 4, 16)
   value_date = to_spot(trade_date, pair)
   write(*,'(a,1x,a)') pair_iso(pair), date_string(value_date)

   index = usdlibor(months_period(3))
   write(*,'(a,1x,a)') 'Three-month USD LIBOR maturity:', &
      date_string(to_maturity(make_date(2017, 1, 3), index))

   items = [single_currency_money(10.0_dp, aud()), &
      single_currency_money(7.5_dp, usd()), &
      single_currency_money(2.0_dp, aud())]
   money = multi_currency_money(items, status)
   totals = aggregate_by_currency(money)
   do i = 1, totals%size()
      write(*,'(a,1x,f10.4)') totals%currency(i)%iso, totals%value(i)
   end do
end program conventions_money_example
