! SPDX-License-Identifier: GPL-2.0-or-later
module ragtop_instruments
   use ragtop_kinds, only : dp
   use ragtop_constants
   use ragtop_types, only : instrument_spec, cashflow_schedule, exercise_schedule, market_spec
   use ragtop_term_structures, only : discount_factor
   use ragtop_cashflows, only : value_from_prior_coupons, accelerated_coupon_value
   implicit none
   private
   public :: EuropeanOption, AmericanOption, ZeroCouponBond, CouponBond
   public :: CallableBond, ConvertibleBond
   public :: terminal_values, recovery_values, apply_optionality
   public :: instrument_cashflow_between, critical_times

contains

   function EuropeanOption(maturity, strike, callput, name) result(inst)
      real(dp), intent(in) :: maturity, strike
      integer, intent(in) :: callput
      character(len=*), intent(in), optional :: name
      type(instrument_spec) :: inst
      inst%kind = instrument_european_option
      inst%maturity = maturity
      inst%strike = strike
      inst%callput = callput
      if (present(name)) inst%name = name
   end function EuropeanOption

   function AmericanOption(maturity, strike, callput, name) result(inst)
      real(dp), intent(in) :: maturity, strike
      integer, intent(in) :: callput
      character(len=*), intent(in), optional :: name
      type(instrument_spec) :: inst
      inst = EuropeanOption(maturity,strike,callput,name)
      inst%kind = instrument_american_option
   end function AmericanOption

   function ZeroCouponBond(maturity, notional, recovery_rate, name) result(inst)
      real(dp), intent(in) :: maturity, notional
      real(dp), intent(in), optional :: recovery_rate
      character(len=*), intent(in), optional :: name
      type(instrument_spec) :: inst
      inst%kind = instrument_zero_coupon_bond
      inst%maturity = maturity
      inst%notional = notional
      if (present(recovery_rate)) inst%recovery_rate = recovery_rate
      if (present(name)) inst%name = name
   end function ZeroCouponBond

   function CouponBond(maturity, notional, coupons, recovery_rate, name) result(inst)
      real(dp), intent(in) :: maturity, notional
      type(cashflow_schedule), intent(in) :: coupons
      real(dp), intent(in), optional :: recovery_rate
      character(len=*), intent(in), optional :: name
      type(instrument_spec) :: inst
      inst%kind = instrument_coupon_bond
      inst%maturity = maturity
      inst%notional = notional
      inst%coupons = coupons
      if (present(recovery_rate)) inst%recovery_rate = recovery_rate
      if (present(name)) inst%name = name
   end function CouponBond

   function CallableBond(maturity, notional, coupons, calls, puts, recovery_rate, name) result(inst)
      real(dp), intent(in) :: maturity, notional
      type(cashflow_schedule), intent(in) :: coupons
      type(exercise_schedule), intent(in), optional :: calls, puts
      real(dp), intent(in), optional :: recovery_rate
      character(len=*), intent(in), optional :: name
      type(instrument_spec) :: inst
      inst%kind = instrument_callable_bond
      inst%maturity = maturity
      inst%notional = notional
      inst%coupons = coupons
      if (present(calls)) inst%calls = calls
      if (present(puts)) inst%puts = puts
      if (present(recovery_rate)) inst%recovery_rate = recovery_rate
      if (present(name)) inst%name = name
   end function CallableBond

   function ConvertibleBond(maturity, notional, conversion_ratio, coupons, recovery_rate, calls, puts, name) result(inst)
      real(dp), intent(in) :: maturity, notional, conversion_ratio
      type(cashflow_schedule), intent(in) :: coupons
      real(dp), intent(in), optional :: recovery_rate
      type(exercise_schedule), intent(in), optional :: calls, puts
      character(len=*), intent(in), optional :: name
      type(instrument_spec) :: inst
      inst%kind = instrument_convertible_bond
      inst%maturity = maturity
      inst%notional = notional
      inst%conversion_ratio = conversion_ratio
      inst%coupons = coupons
      if (present(calls)) inst%calls = calls
      if (present(puts)) inst%puts = puts
      if (present(recovery_rate)) inst%recovery_rate = recovery_rate
      if (present(name)) inst%name = name
   end function ConvertibleBond

   subroutine terminal_values(inst, stock, market, values)
      type(instrument_spec), intent(in) :: inst
      real(dp), intent(in) :: stock(:)
      type(market_spec), intent(in) :: market
      real(dp), intent(out) :: values(:)
      real(dp) :: cp
      select case(inst%kind)
      case(instrument_european_option,instrument_american_option)
         values = max(0.0_dp,real(inst%callput,dp)*(stock-inst%strike))
      case(instrument_zero_coupon_bond,instrument_coupon_bond,instrument_callable_bond)
         values = inst%notional
      case(instrument_convertible_bond)
         values = max(inst%notional,inst%conversion_ratio*stock)
      case default
         values = 0.0_dp
      end select
      cp = value_from_prior_coupons(inst%maturity,inst%coupons,market)
      if (inst%kind >= instrument_coupon_bond .and. inst%kind <= instrument_convertible_bond) then
         ! Coupons are inserted explicitly during backward stepping.  Do not add them twice.
         cp = cp*0.0_dp
      end if
   end subroutine terminal_values

   subroutine recovery_values(inst, hold_values, stock, t, market, values)
      type(instrument_spec), intent(in) :: inst
      real(dp), intent(in) :: hold_values(:), stock(:), t
      type(market_spec), intent(in) :: market
      real(dp), intent(out) :: values(:)
      select case(inst%kind)
      case(instrument_european_option)
         if (inst%callput == put_option .and. t < inst%maturity) then
            values = inst%strike*discount_factor(market,inst%maturity,t)
         else
            values = 0.0_dp
         end if
      case(instrument_american_option)
         if (inst%callput == put_option .and. t < inst%maturity) then
            values = inst%strike
         else
            values = 0.0_dp
         end if
      case(instrument_zero_coupon_bond,instrument_coupon_bond,instrument_callable_bond,instrument_convertible_bond)
         if (t < inst%maturity) then
            values = min(max(hold_values,0.0_dp),inst%recovery_rate*inst%notional)
         else
            values = 0.0_dp
         end if
      case default
         values = 0.0_dp
      end select
      values = values+0.0_dp*stock
   end subroutine recovery_values

   subroutine apply_optionality(inst, hold_values, stock, t, market, values)
      type(instrument_spec), intent(in) :: inst
      real(dp), intent(in) :: hold_values(:), stock(:), t
      type(market_spec), intent(in) :: market
      real(dp), intent(out) :: values(:)
      real(dp) :: exercise_value, accrued
      integer :: i
      logical :: at_event
      values = max(0.0_dp,hold_values)
      select case(inst%kind)
      case(instrument_european_option)
         if (t >= inst%maturity-100.0_dp*epsilon(1.0_dp)) then
            values = max(0.0_dp,real(inst%callput,dp)*(stock-inst%strike))
         end if
      case(instrument_american_option)
         values = max(values,max(0.0_dp,real(inst%callput,dp)*(stock-inst%strike)))
      case(instrument_zero_coupon_bond)
         if (t >= inst%maturity-100.0_dp*epsilon(1.0_dp)) values = inst%notional
      case(instrument_coupon_bond)
         if (t >= inst%maturity-100.0_dp*epsilon(1.0_dp)) values = inst%notional
      case(instrument_callable_bond)
         if (t >= inst%maturity-100.0_dp*epsilon(1.0_dp)) values = inst%notional
         if (allocated(inst%calls%time)) then
            do i = 1, size(inst%calls%time)
               at_event = abs(t-inst%calls%time(i)) <= 1.0e-10_dp*max(1.0_dp,inst%maturity)
               if (at_event) values = min(values,inst%calls%value(i))
            end do
         end if
         if (allocated(inst%puts%time)) then
            do i = 1, size(inst%puts%time)
               at_event = abs(t-inst%puts%time(i)) <= 1.0e-10_dp*max(1.0_dp,inst%maturity)
               if (at_event) values = max(values,inst%puts%value(i))
            end do
         end if
      case(instrument_convertible_bond)
         accrued = value_from_prior_coupons(t,inst%coupons,market)
         exercise_value = 0.0_dp
         if (inst%accelerate_future_coupons) then
            exercise_value = accelerated_coupon_value(t,inst%coupons,market,inst%acceleration_time)
         end if
         values = max(values,inst%conversion_ratio*stock+accrued+exercise_value)
         if (t >= inst%maturity-100.0_dp*epsilon(1.0_dp)) then
            values = max(inst%notional,inst%conversion_ratio*stock)
         end if
         if (allocated(inst%calls%time)) then
            do i = 1, size(inst%calls%time)
               at_event = abs(t-inst%calls%time(i)) <= 1.0e-10_dp*max(1.0_dp,inst%maturity)
               if (at_event) values = min(values,max(inst%calls%value(i),inst%conversion_ratio*stock))
            end do
         end if
         if (allocated(inst%puts%time)) then
            do i = 1, size(inst%puts%time)
               at_event = abs(t-inst%puts%time(i)) <= 1.0e-10_dp*max(1.0_dp,inst%maturity)
               if (at_event) values = max(values,inst%puts%value(i))
            end do
         end if
      end select
   end subroutine apply_optionality

   pure real(dp) function instrument_cashflow_between(inst, t0, t1, market) result(cf)
      type(instrument_spec), intent(in) :: inst
      real(dp), intent(in) :: t0, t1
      type(market_spec), intent(in) :: market
      integer :: i
      cf = 0.0_dp
      if (allocated(inst%coupons%time)) then
         do i = 1, size(inst%coupons%time)
            if (inst%coupons%time(i) > t0 .and. inst%coupons%time(i) <= t1+100.0_dp*epsilon(1.0_dp)) then
               cf = cf+inst%coupons%amount(i)*discount_factor(market,inst%coupons%time(i),t1)
            end if
         end do
      end if
   end function instrument_cashflow_between

   subroutine critical_times(inst, times)
      type(instrument_spec), intent(in) :: inst
      real(dp), allocatable, intent(out) :: times(:)
      integer :: n, k
      n = 1
      if (allocated(inst%coupons%time)) n = n+size(inst%coupons%time)
      if (allocated(inst%calls%time)) n = n+size(inst%calls%time)
      if (allocated(inst%puts%time)) n = n+size(inst%puts%time)
      allocate(times(n))
      k = 1
      times(k) = inst%maturity
      if (allocated(inst%coupons%time)) then
         times(k+1:k+size(inst%coupons%time)) = inst%coupons%time
         k = k+size(inst%coupons%time)
      end if
      if (allocated(inst%calls%time)) then
         times(k+1:k+size(inst%calls%time)) = inst%calls%time
         k = k+size(inst%calls%time)
      end if
      if (allocated(inst%puts%time)) then
         times(k+1:k+size(inst%puts%time)) = inst%puts%time
      end if
   end subroutine critical_times

end module ragtop_instruments
