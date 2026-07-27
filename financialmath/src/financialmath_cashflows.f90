! SPDX-License-Identifier: GPL-2.0-only
! Based on FinancialMath 0.1.1, Copyright (C) 2016 Kameron Penn and Jack Schmidt.
module financialmath_cashflows
   use financialmath_kinds, only : dp
   use financialmath_types, only : rate_conversion_t, tvm_result_t, cashflow_analysis_t, status_t
   use financialmath_math, only : solve_root, finite_real
   implicit none
   private
   public :: rate_conv, solve_tvm, npv, irr, cf_analysis
   public :: swap_rate, swap_commodity, yield_dollar, yield_time
   public :: annual_effective_rate, payment_period_rate

contains

   pure real(dp) function annual_effective_rate(rate, frequency) result(eff)
      real(dp), intent(in) :: rate, frequency
      if (frequency <= 0.0_dp) then
         eff = -1.0_dp
      else
         eff = (1.0_dp + rate/frequency)**frequency - 1.0_dp
      end if
   end function annual_effective_rate

   pure real(dp) function payment_period_rate(rate, compounding_frequency, payment_frequency) result(period_rate)
      real(dp), intent(in) :: rate, compounding_frequency, payment_frequency
      real(dp) :: eff
      eff = annual_effective_rate(rate, compounding_frequency)
      if (payment_frequency <= 0.0_dp .or. eff <= -1.0_dp) then
         period_rate = -1.0_dp
      else
         period_rate = (1.0_dp + eff)**(1.0_dp/payment_frequency) - 1.0_dp
      end if
   end function payment_period_rate

   function rate_conv(rate, conv, rate_type, nom) result(out)
      real(dp), intent(in) :: rate, conv, nom
      character(len=*), intent(in) :: rate_type
      type(rate_conversion_t) :: out
      real(dp) :: eff_i

      if (rate <= 0.0_dp .or. conv <= 0.0_dp .or. nom <= 0.0_dp) return
      select case (trim(rate_type))
      case ('interest')
         eff_i = (1.0_dp + rate/conv)**conv - 1.0_dp
      case ('discount')
         if (rate/conv >= 1.0_dp) return
         eff_i = (1.0_dp - rate/conv)**(-conv) - 1.0_dp
      case ('force')
         eff_i = exp(rate) - 1.0_dp
      case default
         return
      end select
      out%effective_interest = eff_i
      out%effective_discount = eff_i/(1.0_dp + eff_i)
      out%force = log(1.0_dp + eff_i)
      out%input_nominal_interest = conv*((1.0_dp + eff_i)**(1.0_dp/conv) - 1.0_dp)
      out%input_nominal_discount = conv*(1.0_dp - (1.0_dp - out%effective_discount)**(1.0_dp/conv))
      out%target_nominal_interest = nom*((1.0_dp + eff_i)**(1.0_dp/nom) - 1.0_dp)
      out%target_nominal_discount = nom*(1.0_dp - (1.0_dp - out%effective_discount)**(1.0_dp/nom))
   end function rate_conv

   function solve_tvm(pv, fv, periods, rate, compounding_frequency, unknown) result(out)
      real(dp), intent(in) :: pv, fv, periods, rate, compounding_frequency
      character(len=*), intent(in) :: unknown
      type(tvm_result_t) :: out
      real(dp) :: eff

      out%status%ok = .false.
      if (compounding_frequency <= 0.0_dp) then
         out%status%message = 'compounding frequency must be positive'
         return
      end if
      select case (trim(unknown))
      case ('pv')
         eff = annual_effective_rate(rate, compounding_frequency)
         if (fv <= 0.0_dp .or. periods <= 0.0_dp .or. eff <= -1.0_dp) then
            out%status%message = 'invalid TVM inputs'
            return
         end if
         out%present_value = fv/(1.0_dp + eff)**periods
         out%future_value = fv
         out%periods = periods
         out%effective_rate = eff
         out%nominal_rate = rate
      case ('fv')
         eff = annual_effective_rate(rate, compounding_frequency)
         if (pv <= 0.0_dp .or. periods <= 0.0_dp .or. eff <= -1.0_dp) then
            out%status%message = 'invalid TVM inputs'
            return
         end if
         out%present_value = pv
         out%future_value = pv*(1.0_dp + eff)**periods
         out%periods = periods
         out%effective_rate = eff
         out%nominal_rate = rate
      case ('periods', 'n')
         eff = annual_effective_rate(rate, compounding_frequency)
         if (pv <= 0.0_dp .or. fv <= 0.0_dp .or. eff <= 0.0_dp) then
            out%status%message = 'invalid TVM inputs'
            return
         end if
         out%present_value = pv
         out%future_value = fv
         out%periods = log(fv/pv)/log(1.0_dp + eff)
         out%effective_rate = eff
         out%nominal_rate = rate
      case ('rate', 'i')
         if (pv <= 0.0_dp .or. fv <= 0.0_dp .or. periods <= 0.0_dp) then
            out%status%message = 'invalid TVM inputs'
            return
         end if
         eff = (fv/pv)**(1.0_dp/periods) - 1.0_dp
         out%present_value = pv
         out%future_value = fv
         out%periods = periods
         out%effective_rate = eff
         out%nominal_rate = compounding_frequency*((1.0_dp + eff)**(1.0_dp/compounding_frequency) - 1.0_dp)
      case default
         out%status%message = 'unknown must be pv, fv, periods, or rate'
         return
      end select
      out%status%ok = .true.
   end function solve_tvm

   pure real(dp) function npv(cf0, cf, times, rate) result(value)
      real(dp), intent(in) :: cf0, cf(:), times(:), rate
      if (size(cf) /= size(times) .or. rate <= -1.0_dp) then
         value = huge(1.0_dp)
      else
         value = sum(cf/(1.0_dp + rate)**times) - abs(cf0)
      end if
   end function npv

   function irr(cf0, cf, times, ok, lower, upper) result(rate)
      real(dp), intent(in) :: cf0, cf(:), times(:)
      logical, intent(out) :: ok
      real(dp), intent(in), optional :: lower, upper
      real(dp) :: rate, lo, hi, f_lo, f_hi
      integer :: j

      ok = .false.
      rate = 0.0_dp
      if (size(cf) /= size(times) .or. size(cf) == 0) return
      if (any(times <= 0.0_dp)) return
      lo = -0.999999_dp
      hi = 1.0_dp
      if (present(lower)) lo = lower
      if (present(upper)) hi = upper
      f_lo = objective(lo)
      f_hi = objective(hi)
      if (f_lo*f_hi > 0.0_dp .and. .not. present(upper)) then
         do j = 1, 20
            hi = hi*2.0_dp + 1.0_dp
            f_hi = objective(hi)
            if (f_lo*f_hi <= 0.0_dp) exit
         end do
      end if
      rate = solve_root(objective, lo, hi, ok)
   contains
      function objective(x) result(y)
         real(dp), intent(in) :: x
         real(dp) :: y
         y = npv(cf0, cf, times, x)
      end function objective
   end function irr

   function cf_analysis(cf, times, rate) result(out)
      real(dp), intent(in) :: cf(:), times(:), rate
      type(cashflow_analysis_t) :: out
      real(dp) :: pv

      out%status%ok = .false.
      if (size(cf) /= size(times) .or. size(cf) == 0 .or. rate < 0.0_dp) then
         out%status%message = 'invalid cash-flow inputs'
         return
      end if
      if (any(times < 0.0_dp)) then
         out%status%message = 'times must be nonnegative'
         return
      end if
      pv = sum(cf/(1.0_dp + rate)**times)
      out%present_value = pv
      if (abs(pv) > tiny(1.0_dp)) then
         out%macaulay_duration = sum(times*cf/(1.0_dp + rate)**times)/pv
         out%modified_duration = sum(times*cf/(1.0_dp + rate)**(times+1.0_dp))/pv
         out%macaulay_convexity = sum(times*times*cf/(1.0_dp + rate)**times)/pv
         out%modified_convexity = sum(times*(times+1.0_dp)*cf/(1.0_dp + rate)**(times+2.0_dp))/pv
      end if
      out%status%ok = finite_real(pv)
      if (.not. out%status%ok) out%status%message = 'non-finite cash-flow result'
   end function cf_analysis

   function swap_rate(rates, rate_type, ok) result(value)
      real(dp), intent(in) :: rates(:)
      character(len=*), intent(in) :: rate_type
      logical, intent(out), optional :: ok
      real(dp) :: value, den
      integer :: j, n
      logical :: good
      n = size(rates)
      good = n > 0 .and. all(rates >= 0.0_dp)
      value = 0.0_dp
      if (good) then
         select case (trim(rate_type))
         case ('spot_rate')
            den = 0.0_dp
            do j = 1, n
               den = den + (1.0_dp + rates(j))**(-real(j, dp))
            end do
            if (den > 0.0_dp) then
               value = (1.0_dp - (1.0_dp + rates(n))**(-real(n, dp)))/den
            else
               good = .false.
            end if
         case ('zcb_price')
            den = sum(rates)
            if (den > 0.0_dp) then
               value = (1.0_dp-rates(n))/den
            else
               good = .false.
            end if
         case default
            good = .false.
         end select
      end if
      if (present(ok)) ok = good
   end function swap_rate

   function swap_commodity(prices, rates, rate_type, ok) result(value)
      real(dp), intent(in) :: prices(:), rates(:)
      character(len=*), intent(in) :: rate_type
      logical, intent(out), optional :: ok
      real(dp) :: value, num, den, disc
      integer :: j
      logical :: good
      good = size(prices) == size(rates) .and. size(prices) > 0
      value = 0.0_dp
      num = 0.0_dp
      den = 0.0_dp
      if (good) then
         select case (trim(rate_type))
         case ('spot_rate')
            do j = 1, size(prices)
               disc = (1.0_dp + rates(j))**(-real(j, dp))
               num = num + prices(j)*disc
               den = den + disc
            end do
         case ('zcb_price')
            num = sum(prices*rates)
            den = sum(rates)
         case default
            good = .false.
         end select
      end if
      if (good .and. abs(den) > tiny(1.0_dp)) then
         value = num/den
      else
         good = .false.
      end if
      if (present(ok)) ok = good
   end function swap_commodity

   pure real(dp) function yield_dollar(cf, times, start_value, end_value, end_time) result(value)
      real(dp), intent(in) :: cf(:), times(:), start_value, end_value, end_time
      real(dp) :: interest, denominator
      if (size(cf) /= size(times) .or. end_time <= 0.0_dp .or. any(times < 0.0_dp) .or. any(times > end_time)) then
         value = huge(1.0_dp)
         return
      end if
      interest = end_value-start_value-sum(cf)
      denominator = start_value*end_time + sum(cf*(end_time-times))
      if (abs(denominator) <= tiny(1.0_dp)) then
         value = huge(1.0_dp)
      else
         value = interest/denominator
      end if
   end function yield_dollar

   pure real(dp) function yield_time(cf, balances) result(value)
      real(dp), intent(in) :: cf(:), balances(:)
      integer :: j, n
      value = huge(1.0_dp)
      n = size(balances)-1
      if (n <= 0) return
      if (size(cf) /= n .and. size(cf) /= n+1) return
      value = 1.0_dp
      do j = 1, n
         if (abs(balances(j)+cf(j)) <= tiny(1.0_dp)) then
            value = huge(1.0_dp)
            return
         end if
         value = value*balances(j+1)/(balances(j)+cf(j))
      end do
      value = value-1.0_dp
   end function yield_time

end module financialmath_cashflows
