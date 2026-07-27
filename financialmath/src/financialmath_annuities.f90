! SPDX-License-Identifier: GPL-2.0-only
! Based on FinancialMath 0.1.1, Copyright (C) 2016 Kameron Penn and Jack Schmidt.
module financialmath_annuities
   use financialmath_kinds, only : dp
   use financialmath_types, only : annuity_result_t
   use financialmath_cashflows, only : annual_effective_rate, payment_period_rate
   use financialmath_math, only : solve_root, nearly_equal
   implicit none
   private
   public :: annuity_level, annuity_arith, annuity_geo
   public :: perpetuity_level, perpetuity_arith, perpetuity_geo
   public :: annuity_level_pv, annuity_level_fv
   public :: annuity_arith_pv, annuity_arith_fv
   public :: annuity_geo_pv, annuity_geo_fv

contains

   pure real(dp) function due_factor(period_rate, immediate) result(factor)
      real(dp), intent(in) :: period_rate
      logical, intent(in) :: immediate
      factor = 1.0_dp
      if (.not. immediate) factor = 1.0_dp + period_rate
   end function due_factor

   pure real(dp) function annuity_level_pv(payment, periods, period_rate, immediate) result(value)
      real(dp), intent(in) :: payment, periods, period_rate
      logical, intent(in) :: immediate
      real(dp) :: factor
      factor = due_factor(period_rate, immediate)
      if (abs(period_rate) < 1.0e-14_dp) then
         value = factor*payment*periods
      else
         value = factor*payment*(1.0_dp-(1.0_dp+period_rate)**(-periods))/period_rate
      end if
   end function annuity_level_pv

   pure real(dp) function annuity_level_fv(payment, periods, period_rate, immediate) result(value)
      real(dp), intent(in) :: payment, periods, period_rate
      logical, intent(in) :: immediate
      real(dp) :: factor
      factor = due_factor(period_rate, immediate)
      if (abs(period_rate) < 1.0e-14_dp) then
         value = factor*payment*periods
      else
         value = factor*payment*((1.0_dp+period_rate)**periods-1.0_dp)/period_rate
      end if
   end function annuity_level_fv

   pure real(dp) function annuity_arith_pv(first_payment, increment, periods, period_rate, immediate) result(value)
      real(dp), intent(in) :: first_payment, increment, periods, period_rate
      logical, intent(in) :: immediate
      real(dp) :: factor, an
      factor = due_factor(period_rate, immediate)
      if (abs(period_rate) < 1.0e-10_dp) then
         value = factor*(first_payment*periods + 0.5_dp*increment*periods*(periods-1.0_dp))
      else
         an = (1.0_dp-(1.0_dp+period_rate)**(-periods))/period_rate
         value = factor*(first_payment*an + increment*(an-periods*(1.0_dp+period_rate)**(-periods))/period_rate)
      end if
   end function annuity_arith_pv

   pure real(dp) function annuity_arith_fv(first_payment, increment, periods, period_rate, immediate) result(value)
      real(dp), intent(in) :: first_payment, increment, periods, period_rate
      logical, intent(in) :: immediate
      value = annuity_arith_pv(first_payment, increment, periods, period_rate, immediate)* &
         (1.0_dp+period_rate)**periods
   end function annuity_arith_fv

   pure real(dp) function annuity_geo_pv(first_payment, growth_rate, periods, period_rate, immediate) result(value)
      real(dp), intent(in) :: first_payment, growth_rate, periods, period_rate
      logical, intent(in) :: immediate
      real(dp) :: factor
      factor = due_factor(period_rate, immediate)
      if (nearly_equal(period_rate, growth_rate, 1.0e-12_dp, 1.0e-10_dp)) then
         value = factor*first_payment*periods/(1.0_dp+period_rate)
      else
         value = factor*first_payment*(1.0_dp-((1.0_dp+growth_rate)/(1.0_dp+period_rate))**periods)/ &
            (period_rate-growth_rate)
      end if
   end function annuity_geo_pv

   pure real(dp) function annuity_geo_fv(first_payment, growth_rate, periods, period_rate, immediate) result(value)
      real(dp), intent(in) :: first_payment, growth_rate, periods, period_rate
      logical, intent(in) :: immediate
      value = annuity_geo_pv(first_payment, growth_rate, periods, period_rate, immediate)* &
         (1.0_dp+period_rate)**periods
   end function annuity_geo_fv

   function annuity_level(pv, fv, periods, payment, rate, ic, pf, immediate, unknown) result(out)
      real(dp), intent(in) :: pv, fv, periods, payment, rate, ic, pf
      logical, intent(in) :: immediate
      character(len=*), intent(in) :: unknown
      type(annuity_result_t) :: out
      real(dp) :: j, nominal
      logical :: ok

      out%status%ok = .false.
      if (ic <= 0.0_dp .or. pf <= 0.0_dp) then
         out%status%message = 'frequencies must be positive'
         return
      end if
      j = payment_period_rate(rate, ic, pf)
      nominal = rate
      select case (trim(unknown))
      case ('pv')
         if (payment <= 0.0_dp .or. periods <= 0.0_dp) goto 900
         out%present_value = annuity_level_pv(payment, periods, j, immediate)
         out%future_value = annuity_level_fv(payment, periods, j, immediate)
         out%first_payment = payment
         out%periods = periods
      case ('fv')
         if (payment <= 0.0_dp .or. periods <= 0.0_dp) goto 900
         out%present_value = annuity_level_pv(payment, periods, j, immediate)
         out%future_value = annuity_level_fv(payment, periods, j, immediate)
         out%first_payment = payment
         out%periods = periods
      case ('payment', 'pmt')
         if (periods <= 0.0_dp) goto 900
         if (pv > 0.0_dp) then
            out%first_payment = pv/annuity_level_pv(1.0_dp, periods, j, immediate)
            out%present_value = pv
            out%future_value = pv*(1.0_dp+j)**periods
         else if (fv > 0.0_dp) then
            out%first_payment = fv/annuity_level_fv(1.0_dp, periods, j, immediate)
            out%future_value = fv
            out%present_value = fv/(1.0_dp+j)**periods
         else
            goto 900
         end if
         out%periods = periods
      case ('periods', 'n')
         if (payment <= 0.0_dp .or. j <= 0.0_dp) goto 900
         if (pv > 0.0_dp) then
            out%periods = solve_root(pv_objective, 1.0e-10_dp, 1.0e6_dp, ok)
            if (.not. ok) goto 900
            out%present_value = pv
            out%future_value = pv*(1.0_dp+j)**out%periods
         else if (fv > 0.0_dp) then
            out%periods = solve_root(fv_objective, 1.0e-10_dp, 1.0e6_dp, ok)
            if (.not. ok) goto 900
            out%future_value = fv
            out%present_value = fv/(1.0_dp+j)**out%periods
         else
            goto 900
         end if
         out%first_payment = payment
      case ('rate', 'i')
         if (payment <= 0.0_dp .or. periods <= 0.0_dp) goto 900
         j = solve_root(rate_objective, 1.0e-12_dp, 10.0_dp, ok)
         if (.not. ok) goto 900
         nominal = ic*((1.0_dp+j)**(pf/ic)-1.0_dp)
         out%first_payment = payment
         out%periods = periods
         if (pv > 0.0_dp) then
            out%present_value = pv
            out%future_value = pv*(1.0_dp+j)**periods
         else
            out%future_value = fv
            out%present_value = fv/(1.0_dp+j)**periods
         end if
      case default
         out%status%message = 'unknown must be pv, fv, payment, periods, or rate'
         return
      end select
      out%effective_rate = (1.0_dp+j)**pf-1.0_dp
      out%nominal_rate = nominal
      out%payment_frequency_rate = j*pf
      out%status%ok = .true.
      return
900   out%status%message = 'invalid level-annuity inputs or root not bracketed'
   contains
      function pv_objective(x) result(y)
         real(dp), intent(in) :: x
         real(dp) :: y
         y = annuity_level_pv(payment, x, j, immediate)-pv
      end function pv_objective
      function fv_objective(x) result(y)
         real(dp), intent(in) :: x
         real(dp) :: y
         y = annuity_level_fv(payment, x, j, immediate)-fv
      end function fv_objective
      function rate_objective(x) result(y)
         real(dp), intent(in) :: x
         real(dp) :: y
         if (pv > 0.0_dp) then
            y = annuity_level_pv(payment, periods, x, immediate)-pv
         else
            y = annuity_level_fv(payment, periods, x, immediate)-fv
         end if
      end function rate_objective
   end function annuity_level

   function annuity_arith(pv, fv, periods, first_payment, increment, rate, ic, pf, immediate, unknown) result(out)
      real(dp), intent(in) :: pv, fv, periods, first_payment, increment, rate, ic, pf
      logical, intent(in) :: immediate
      character(len=*), intent(in) :: unknown
      type(annuity_result_t) :: out
      real(dp) :: j, nominal, base, slope
      logical :: ok

      out%status%ok = .false.
      if (ic <= 0.0_dp .or. pf <= 0.0_dp) goto 900
      j = payment_period_rate(rate, ic, pf)
      nominal = rate
      select case (trim(unknown))
      case ('pv')
         out%present_value = annuity_arith_pv(first_payment, increment, periods, j, immediate)
         out%future_value = out%present_value*(1.0_dp+j)**periods
      case ('fv')
         out%future_value = annuity_arith_fv(first_payment, increment, periods, j, immediate)
         out%present_value = out%future_value/(1.0_dp+j)**periods
      case ('first_payment', 'p')
         if (pv > 0.0_dp) then
            base = annuity_arith_pv(0.0_dp, increment, periods, j, immediate)
            slope = annuity_arith_pv(1.0_dp, 0.0_dp, periods, j, immediate)
            out%present_value = pv
            out%first_payment = (pv-base)/slope
            out%future_value = pv*(1.0_dp+j)**periods
         else if (fv > 0.0_dp) then
            base = annuity_arith_fv(0.0_dp, increment, periods, j, immediate)
            slope = annuity_arith_fv(1.0_dp, 0.0_dp, periods, j, immediate)
            out%future_value = fv
            out%first_payment = (fv-base)/slope
            out%present_value = fv/(1.0_dp+j)**periods
         else
            goto 900
         end if
      case ('increment', 'q')
         if (pv > 0.0_dp) then
            base = annuity_arith_pv(first_payment, 0.0_dp, periods, j, immediate)
            slope = annuity_arith_pv(0.0_dp, 1.0_dp, periods, j, immediate)
            out%present_value = pv
            out%increment = (pv-base)/slope
            out%future_value = pv*(1.0_dp+j)**periods
         else if (fv > 0.0_dp) then
            base = annuity_arith_fv(first_payment, 0.0_dp, periods, j, immediate)
            slope = annuity_arith_fv(0.0_dp, 1.0_dp, periods, j, immediate)
            out%future_value = fv
            out%increment = (fv-base)/slope
            out%present_value = fv/(1.0_dp+j)**periods
         else
            goto 900
         end if
      case ('periods', 'n')
         out%periods = solve_root(n_objective, 1.0e-8_dp, 1.0e5_dp, ok)
         if (.not. ok) goto 900
         if (pv > 0.0_dp) then
            out%present_value = pv
            out%future_value = pv*(1.0_dp+j)**out%periods
         else
            out%future_value = fv
            out%present_value = fv/(1.0_dp+j)**out%periods
         end if
      case ('rate', 'i')
         j = solve_root(rate_objective, 1.0e-12_dp, 10.0_dp, ok)
         if (.not. ok) goto 900
         nominal = ic*((1.0_dp+j)**(pf/ic)-1.0_dp)
         if (pv > 0.0_dp) then
            out%present_value = pv
            out%future_value = pv*(1.0_dp+j)**periods
         else
            out%future_value = fv
            out%present_value = fv/(1.0_dp+j)**periods
         end if
      case default
         out%status%message = 'unknown must be pv, fv, first_payment, increment, periods, or rate'
         return
      end select
      if (trim(unknown) /= 'first_payment' .and. trim(unknown) /= 'p') out%first_payment = first_payment
      if (trim(unknown) /= 'increment' .and. trim(unknown) /= 'q') out%increment = increment
      if (trim(unknown) /= 'periods' .and. trim(unknown) /= 'n') out%periods = periods
      out%effective_rate = (1.0_dp+j)**pf-1.0_dp
      out%nominal_rate = nominal
      out%payment_frequency_rate = j*pf
      out%status%ok = .true.
      return
900   out%status%message = 'invalid arithmetic-annuity inputs or root not bracketed'
   contains
      function n_objective(x) result(y)
         real(dp), intent(in) :: x
         real(dp) :: y
         if (pv > 0.0_dp) then
            y = annuity_arith_pv(first_payment, increment, x, j, immediate)-pv
         else
            y = annuity_arith_fv(first_payment, increment, x, j, immediate)-fv
         end if
      end function n_objective
      function rate_objective(x) result(y)
         real(dp), intent(in) :: x
         real(dp) :: y
         if (pv > 0.0_dp) then
            y = annuity_arith_pv(first_payment, increment, periods, x, immediate)-pv
         else
            y = annuity_arith_fv(first_payment, increment, periods, x, immediate)-fv
         end if
      end function rate_objective
   end function annuity_arith

   function annuity_geo(pv, fv, periods, first_payment, growth_rate, rate, ic, pf, immediate, unknown) result(out)
      real(dp), intent(in) :: pv, fv, periods, first_payment, growth_rate, rate, ic, pf
      logical, intent(in) :: immediate
      character(len=*), intent(in) :: unknown
      type(annuity_result_t) :: out
      real(dp) :: j, nominal, base
      logical :: ok

      out%status%ok = .false.
      if (ic <= 0.0_dp .or. pf <= 0.0_dp) goto 900
      j = payment_period_rate(rate, ic, pf)
      nominal = rate
      select case (trim(unknown))
      case ('pv')
         out%present_value = annuity_geo_pv(first_payment, growth_rate, periods, j, immediate)
         out%future_value = out%present_value*(1.0_dp+j)**periods
      case ('fv')
         out%future_value = annuity_geo_fv(first_payment, growth_rate, periods, j, immediate)
         out%present_value = out%future_value/(1.0_dp+j)**periods
      case ('first_payment', 'p')
         if (pv > 0.0_dp) then
            base = annuity_geo_pv(1.0_dp, growth_rate, periods, j, immediate)
            out%first_payment = pv/base
            out%present_value = pv
            out%future_value = pv*(1.0_dp+j)**periods
         else if (fv > 0.0_dp) then
            base = annuity_geo_fv(1.0_dp, growth_rate, periods, j, immediate)
            out%first_payment = fv/base
            out%future_value = fv
            out%present_value = fv/(1.0_dp+j)**periods
         else
            goto 900
         end if
      case ('growth_rate', 'k')
         out%growth_rate = solve_root(growth_objective, -0.999999_dp, max(10.0_dp, 2.0_dp+j), ok)
         if (.not. ok) goto 900
         if (pv > 0.0_dp) then
            out%present_value = pv
            out%future_value = pv*(1.0_dp+j)**periods
         else
            out%future_value = fv
            out%present_value = fv/(1.0_dp+j)**periods
         end if
      case ('periods', 'n')
         out%periods = solve_root(n_objective, 1.0e-8_dp, 1.0e5_dp, ok)
         if (.not. ok) goto 900
         if (pv > 0.0_dp) then
            out%present_value = pv
            out%future_value = pv*(1.0_dp+j)**out%periods
         else
            out%future_value = fv
            out%present_value = fv/(1.0_dp+j)**out%periods
         end if
      case ('rate', 'i')
         j = solve_root(rate_objective, 1.0e-12_dp, 10.0_dp, ok)
         if (.not. ok) goto 900
         nominal = ic*((1.0_dp+j)**(pf/ic)-1.0_dp)
         if (pv > 0.0_dp) then
            out%present_value = pv
            out%future_value = pv*(1.0_dp+j)**periods
         else
            out%future_value = fv
            out%present_value = fv/(1.0_dp+j)**periods
         end if
      case default
         out%status%message = 'unknown must be pv, fv, first_payment, growth_rate, periods, or rate'
         return
      end select
      if (trim(unknown) /= 'first_payment' .and. trim(unknown) /= 'p') out%first_payment = first_payment
      if (trim(unknown) /= 'growth_rate' .and. trim(unknown) /= 'k') out%growth_rate = growth_rate
      if (trim(unknown) /= 'periods' .and. trim(unknown) /= 'n') out%periods = periods
      out%effective_rate = (1.0_dp+j)**pf-1.0_dp
      out%nominal_rate = nominal
      out%payment_frequency_rate = j*pf
      out%status%ok = .true.
      return
900   out%status%message = 'invalid geometric-annuity inputs or root not bracketed'
   contains
      function growth_objective(x) result(y)
         real(dp), intent(in) :: x
         real(dp) :: y
         if (pv > 0.0_dp) then
            y = annuity_geo_pv(first_payment, x, periods, j, immediate)-pv
         else
            y = annuity_geo_fv(first_payment, x, periods, j, immediate)-fv
         end if
      end function growth_objective
      function n_objective(x) result(y)
         real(dp), intent(in) :: x
         real(dp) :: y
         if (pv > 0.0_dp) then
            y = annuity_geo_pv(first_payment, growth_rate, x, j, immediate)-pv
         else
            y = annuity_geo_fv(first_payment, growth_rate, x, j, immediate)-fv
         end if
      end function n_objective
      function rate_objective(x) result(y)
         real(dp), intent(in) :: x
         real(dp) :: y
         if (pv > 0.0_dp) then
            y = annuity_geo_pv(first_payment, growth_rate, periods, x, immediate)-pv
         else
            y = annuity_geo_fv(first_payment, growth_rate, periods, x, immediate)-fv
         end if
      end function rate_objective
   end function annuity_geo

   function perpetuity_level(pv, payment, rate, ic, pf, immediate, unknown) result(out)
      real(dp), intent(in) :: pv, payment, rate, ic, pf
      logical, intent(in) :: immediate
      character(len=*), intent(in) :: unknown
      type(annuity_result_t) :: out
      real(dp) :: j, factor, nominal

      out%status%ok = .false.
      j = payment_period_rate(rate, ic, pf)
      nominal = rate
      factor = due_factor(j, immediate)
      select case (trim(unknown))
      case ('pv')
         if (j <= 0.0_dp) goto 900
         out%present_value = payment*factor/j
         out%first_payment = payment
      case ('payment', 'pmt')
         if (pv <= 0.0_dp) goto 900
         out%present_value = pv
         out%first_payment = pv*j/factor
      case ('rate', 'i')
         if (pv <= 0.0_dp .or. payment <= 0.0_dp) goto 900
         if (immediate) then
            j = payment/pv
         else
            if (pv <= payment) goto 900
            j = payment/(pv-payment)
         end if
         nominal = ic*((1.0_dp+j)**(pf/ic)-1.0_dp)
         out%present_value = pv
         out%first_payment = payment
      case default
         out%status%message = 'unknown must be pv, payment, or rate'
         return
      end select
      out%periods = huge(1.0_dp)
      out%effective_rate = (1.0_dp+j)**pf-1.0_dp
      out%nominal_rate = nominal
      out%payment_frequency_rate = j*pf
      out%status%ok = .true.
      return
900   out%status%message = 'invalid level-perpetuity inputs'
   end function perpetuity_level

   function perpetuity_arith(pv, first_payment, increment, rate, ic, pf, immediate, unknown) result(out)
      real(dp), intent(in) :: pv, first_payment, increment, rate, ic, pf
      logical, intent(in) :: immediate
      character(len=*), intent(in) :: unknown
      type(annuity_result_t) :: out
      real(dp) :: j, factor, nominal
      logical :: ok

      out%status%ok = .false.
      j = payment_period_rate(rate, ic, pf)
      nominal = rate
      factor = due_factor(j, immediate)
      select case (trim(unknown))
      case ('pv')
         if (j <= 0.0_dp) goto 900
         out%present_value = factor*(first_payment/j + increment/(j*j))
      case ('first_payment', 'p')
         if (j <= 0.0_dp) goto 900
         out%present_value = pv
         out%first_payment = j*(pv/factor-increment/(j*j))
      case ('increment', 'q')
         if (j <= 0.0_dp) goto 900
         out%present_value = pv
         out%increment = j*j*(pv/factor-first_payment/j)
      case ('rate', 'i')
         j = solve_root(rate_objective, 1.0e-12_dp, 10.0_dp, ok)
         if (.not. ok) goto 900
         nominal = ic*((1.0_dp+j)**(pf/ic)-1.0_dp)
         out%present_value = pv
      case default
         out%status%message = 'unknown must be pv, first_payment, increment, or rate'
         return
      end select
      if (trim(unknown) /= 'first_payment' .and. trim(unknown) /= 'p') out%first_payment = first_payment
      if (trim(unknown) /= 'increment' .and. trim(unknown) /= 'q') out%increment = increment
      out%periods = huge(1.0_dp)
      out%effective_rate = (1.0_dp+j)**pf-1.0_dp
      out%nominal_rate = nominal
      out%payment_frequency_rate = j*pf
      out%status%ok = .true.
      return
900   out%status%message = 'invalid arithmetic-perpetuity inputs or root not bracketed'
   contains
      function rate_objective(x) result(y)
         real(dp), intent(in) :: x
         real(dp) :: y, f
         f = due_factor(x, immediate)
         y = f*(first_payment/x+increment/(x*x))-pv
      end function rate_objective
   end function perpetuity_arith

   function perpetuity_geo(pv, first_payment, growth_rate, rate, ic, pf, immediate, unknown) result(out)
      real(dp), intent(in) :: pv, first_payment, growth_rate, rate, ic, pf
      logical, intent(in) :: immediate
      character(len=*), intent(in) :: unknown
      type(annuity_result_t) :: out
      real(dp) :: j, factor, nominal

      out%status%ok = .false.
      j = payment_period_rate(rate, ic, pf)
      nominal = rate
      factor = due_factor(j, immediate)
      select case (trim(unknown))
      case ('pv')
         if (j <= growth_rate) goto 900
         out%present_value = factor*first_payment/(j-growth_rate)
      case ('first_payment', 'p')
         if (j <= growth_rate) goto 900
         out%present_value = pv
         out%first_payment = pv*(j-growth_rate)/factor
      case ('growth_rate', 'k')
         out%present_value = pv
         out%growth_rate = j-factor*first_payment/pv
         if (out%growth_rate >= j) goto 900
      case ('rate', 'i')
         if (pv <= 0.0_dp .or. first_payment <= 0.0_dp) goto 900
         if (immediate) then
            j = growth_rate+first_payment/pv
         else
            ! pv = (1+j) p / (j-k); solve directly.
            if (abs(pv-first_payment) < 1.0e-14_dp) goto 900
            j = (pv*growth_rate+first_payment)/(pv-first_payment)
         end if
         nominal = ic*((1.0_dp+j)**(pf/ic)-1.0_dp)
         out%present_value = pv
      case default
         out%status%message = 'unknown must be pv, first_payment, growth_rate, or rate'
         return
      end select
      if (trim(unknown) /= 'first_payment' .and. trim(unknown) /= 'p') out%first_payment = first_payment
      if (trim(unknown) /= 'growth_rate' .and. trim(unknown) /= 'k') out%growth_rate = growth_rate
      out%periods = huge(1.0_dp)
      out%effective_rate = (1.0_dp+j)**pf-1.0_dp
      out%nominal_rate = nominal
      out%payment_frequency_rate = j*pf
      out%status%ok = .true.
      return
900   out%status%message = 'invalid geometric-perpetuity inputs'
   end function perpetuity_geo

end module financialmath_annuities
