! SPDX-License-Identifier: GPL-2.0-or-later
module fincal_accounting
   use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
   use fincal_kinds, only : dp
   use fincal_status, only : fincal_ok, fincal_invalid_input, fincal_size_mismatch, &
      fincal_insufficient_inventory
   use fincal_types, only : inventory_result
   implicit none
   private

   public :: eps, diluted_eps, weighted_average_shares, issuable_shares
   public :: straight_line_depreciation, double_declining_balance, cogs
   public :: was, iss, slde, ddb

   interface was
      module procedure weighted_average_shares
   end interface was
   interface iss
      module procedure issuable_shares
   end interface iss
   interface slde
      module procedure straight_line_depreciation
   end interface slde
   interface ddb
      module procedure double_declining_balance
   end interface ddb
contains
   elemental pure function eps(net_income, preferred_dividends, weighted_average_common_shares) result(value)
      real(dp), intent(in) :: net_income, preferred_dividends, weighted_average_common_shares
      real(dp) :: value
      value = (net_income - preferred_dividends) / weighted_average_common_shares
   end function eps

   elemental pure function diluted_eps(net_income, preferred_dividends, weighted_average_common_shares, &
      convertible_preferred_dividends, convertible_debt_interest, tax_rate, &
      convertible_preferred_shares, convertible_debt_shares, option_shares) result(value)
      real(dp), intent(in) :: net_income, preferred_dividends, weighted_average_common_shares
      real(dp), intent(in), optional :: convertible_preferred_dividends, convertible_debt_interest, tax_rate
      real(dp), intent(in), optional :: convertible_preferred_shares, convertible_debt_shares, option_shares
      real(dp) :: value
      real(dp) :: cpd, cdi, tax, cps, cds, options, basic

      cpd = 0.0_dp
      cdi = 0.0_dp
      tax = 0.0_dp
      cps = 0.0_dp
      cds = 0.0_dp
      options = 0.0_dp
      if (present(convertible_preferred_dividends)) cpd = convertible_preferred_dividends
      if (present(convertible_debt_interest)) cdi = convertible_debt_interest
      if (present(tax_rate)) tax = tax_rate
      if (present(convertible_preferred_shares)) cps = convertible_preferred_shares
      if (present(convertible_debt_shares)) cds = convertible_debt_shares
      if (present(option_shares)) options = option_shares

      basic = eps(net_income, preferred_dividends, weighted_average_common_shares)
      value = (net_income - preferred_dividends + cpd + cdi * (1.0_dp - tax)) / &
         (weighted_average_common_shares + cps + cds + options)
      if (value > basic) then
         value = (net_income - preferred_dividends + cpd) / &
            (weighted_average_common_shares + cps + options)
      end if
   end function diluted_eps

   function weighted_average_shares(number_of_shares, months, status) result(value)
      real(dp), intent(in) :: number_of_shares(:), months(:)
      integer, intent(out), optional :: status
      real(dp) :: value

      if (size(number_of_shares) /= size(months)) then
         value = ieee_value(0.0_dp, ieee_quiet_nan)
         if (present(status)) status = fincal_size_mismatch
         return
      end if
      value = sum(number_of_shares * months) / 12.0_dp
      if (present(status)) status = fincal_ok
   end function weighted_average_shares

   function issuable_shares(average_market_price, exercise_price, option_count, status) result(value)
      real(dp), intent(in) :: average_market_price, exercise_price, option_count
      integer, intent(out), optional :: status
      real(dp) :: value

      if (average_market_price <= exercise_price) then
         value = ieee_value(0.0_dp, ieee_quiet_nan)
         if (present(status)) status = fincal_invalid_input
         return
      end if
      value = (average_market_price - exercise_price) * option_count / average_market_price
      if (present(status)) status = fincal_ok
   end function issuable_shares

   elemental pure function straight_line_depreciation(cost, residual_value, life) result(value)
      real(dp), intent(in) :: cost, residual_value
      integer, intent(in) :: life
      real(dp) :: value
      value = (cost - residual_value) / real(life, dp)
   end function straight_line_depreciation

   function double_declining_balance(cost, residual_value, life, status) result(expense)
      real(dp), intent(in) :: cost, residual_value
      integer, intent(in) :: life
      integer, intent(out), optional :: status
      real(dp), allocatable :: expense(:)
      real(dp) :: book_value, proposed
      integer :: i

      if (life < 2 .or. residual_value > cost) then
         allocate(expense(0))
         if (present(status)) status = fincal_invalid_input
         return
      end if

      allocate(expense(life), source = 0.0_dp)
      book_value = cost
      do i = 1, life
         proposed = 2.0_dp * book_value / real(life, dp)
         expense(i) = min(proposed, max(0.0_dp, book_value - residual_value))
         book_value = book_value - expense(i)
         if (book_value <= residual_value) exit
      end do
      if (present(status)) status = fincal_ok
   end function double_declining_balance

   function cogs(beginning_units, beginning_price, units, prices, sold_units, method) result(result_value)
      real(dp), intent(in) :: beginning_units, beginning_price
      real(dp), intent(in) :: units(:), prices(:), sold_units
      character(len=*), intent(in), optional :: method
      type(inventory_result) :: result_value
      real(dp), allocatable :: lot_units(:), lot_prices(:)
      real(dp) :: total_units, total_cost, remaining_sale, take
      character(len=4) :: selected
      integer :: i, n

      if (size(units) /= size(prices)) then
         result_value%status = fincal_size_mismatch
         return
      end if
      if (beginning_units < 0.0_dp .or. sold_units < 0.0_dp .or. any(units < 0.0_dp)) then
         result_value%status = fincal_invalid_input
         return
      end if

      selected = 'FIFO'
      if (present(method)) selected = upper_method(method)
      if (selected /= 'FIFO' .and. selected /= 'LIFO' .and. selected /= 'WAC ') then
         result_value%status = fincal_invalid_input
         return
      end if

      n = size(units) + 1
      allocate(lot_units(n), lot_prices(n))
      lot_units(1) = beginning_units
      lot_prices(1) = beginning_price
      if (n > 1) then
         lot_units(2:n) = units
         lot_prices(2:n) = prices
      end if

      total_units = sum(lot_units)
      total_cost = sum(lot_units * lot_prices)
      if (sold_units > total_units + 100.0_dp * epsilon(max(1.0_dp, total_units))) then
         result_value%status = fincal_insufficient_inventory
         return
      end if

      if (selected == 'WAC ') then
         if (abs(total_units) <= tiny(1.0_dp)) then
            result_value%cost_of_goods = 0.0_dp
            result_value%ending_inventory = 0.0_dp
         else
            result_value%cost_of_goods = total_cost * sold_units / total_units
            result_value%ending_inventory = total_cost - result_value%cost_of_goods
         end if
      else
         remaining_sale = sold_units
         if (selected == 'FIFO') then
            do i = 1, n
               take = min(remaining_sale, lot_units(i))
               result_value%cost_of_goods = result_value%cost_of_goods + take * lot_prices(i)
               lot_units(i) = lot_units(i) - take
               remaining_sale = remaining_sale - take
               if (remaining_sale <= 0.0_dp) exit
            end do
         else
            do i = n, 1, -1
               take = min(remaining_sale, lot_units(i))
               result_value%cost_of_goods = result_value%cost_of_goods + take * lot_prices(i)
               lot_units(i) = lot_units(i) - take
               remaining_sale = remaining_sale - take
               if (remaining_sale <= 0.0_dp) exit
            end do
         end if
         result_value%ending_inventory = sum(lot_units * lot_prices)
      end if
      result_value%status = fincal_ok
   end function cogs

   pure function upper_method(text) result(method)
      character(len=*), intent(in) :: text
      character(len=4) :: method
      character(len=:), allocatable :: cleaned
      integer :: i, code

      cleaned = adjustl(trim(text))
      method = '    '
      do i = 1, min(4, len(cleaned))
         code = iachar(cleaned(i:i))
         if (code >= iachar('a') .and. code <= iachar('z')) then
            method(i:i) = achar(code - 32)
         else
            method(i:i) = cleaned(i:i)
         end if
      end do
   end function upper_method
end module fincal_accounting
