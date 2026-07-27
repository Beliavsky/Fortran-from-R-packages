! SPDX-License-Identifier: GPL-2.0-only
! Based on FinancialMath 0.1.1, Copyright (C) 2016 Kameron Penn and Jack Schmidt.
module financialmath_loans
   use financialmath_kinds, only : dp
   use financialmath_types, only : amortization_result_t, amort_period_result_t, bond_result_t
   use financialmath_cashflows, only : annual_effective_rate, payment_period_rate
   implicit none
   private
   public :: amort_table, amort_period, bond

contains

   function amort_table(loan, periods, payment, rate, ic, pf, unknown) result(out)
      real(dp), intent(in) :: loan, periods, payment, rate, ic, pf
      character(len=*), intent(in) :: unknown
      type(amortization_result_t) :: out
      real(dp) :: j, n_exact, balance, interest_part, principal_part
      real(dp) :: frac, final_payment
      integer :: n_full, n_rows, row

      out%status%ok = .false.
      if (ic <= 0.0_dp .or. pf <= 0.0_dp .or. rate <= 0.0_dp) then
         out%status%message = 'rate and frequencies must be positive'
         return
      end if
      j = payment_period_rate(rate, ic, pf)
      select case (trim(unknown))
      case ('loan')
         if (periods <= 0.0_dp .or. payment <= 0.0_dp) goto 900
         out%loan = payment*(1.0_dp-(1.0_dp+j)**(-periods))/j
         out%payment = payment
         n_exact = periods
      case ('payment', 'pmt')
         if (periods <= 0.0_dp .or. loan <= 0.0_dp) goto 900
         out%loan = loan
         out%payment = loan*j/(1.0_dp-(1.0_dp+j)**(-periods))
         n_exact = periods
      case ('periods', 'n')
         if (loan <= 0.0_dp .or. payment <= loan*j) goto 900
         out%loan = loan
         out%payment = payment
         n_exact = -log(1.0_dp-loan*j/payment)/log(1.0_dp+j)
      case default
         out%status%message = 'unknown must be loan, payment, or periods'
         return
      end select

      out%periods = n_exact
      n_full = floor(n_exact+1.0e-12_dp)
      frac = n_exact-real(n_full, dp)
      n_rows = n_full
      if (frac > 1.0e-10_dp) n_rows = n_rows+1
      if (n_rows <= 0) goto 900
      allocate(out%time(n_rows), out%payments(n_rows), out%interest(n_rows), &
         out%principal(n_rows), out%balance(n_rows))
      balance = out%loan
      do row = 1, n_full
         interest_part = balance*j
         principal_part = out%payment-interest_part
         if (principal_part > balance .or. row == n_rows) then
            principal_part = balance
            out%payments(row) = interest_part+principal_part
         else
            out%payments(row) = out%payment
         end if
         balance = max(0.0_dp, balance-principal_part)
         out%time(row) = real(row, dp)/pf
         out%interest(row) = interest_part
         out%principal(row) = principal_part
         out%balance(row) = balance
      end do
      if (frac > 1.0e-10_dp) then
         row = n_rows
         interest_part = balance*((1.0_dp+j)**frac-1.0_dp)
         final_payment = balance+interest_part
         out%time(row) = n_exact/pf
         out%payments(row) = final_payment
         out%interest(row) = interest_part
         out%principal(row) = balance
         out%balance(row) = 0.0_dp
         balance = 0.0_dp
      end if
      out%total_paid = sum(out%payments)
      out%total_interest = sum(out%interest)
      if (trim(unknown) == 'periods' .or. trim(unknown) == 'n') then
         out%balloon_payment = loan*(1.0_dp+j)**real(n_full, dp) - &
            payment*((1.0_dp+j)**real(n_full, dp)-1.0_dp)/j + payment
         out%drop_payment = loan*(1.0_dp+j)**real(n_full+1, dp) - &
            payment*((1.0_dp+j)**real(n_full, dp)-1.0_dp)*(1.0_dp+j)/j
      end if
      out%status%ok = .true.
      return
900   out%status%message = 'invalid amortization inputs'
   end function amort_table

   function amort_period(loan, periods, payment, rate, ic, pf, time, unknown) result(out)
      real(dp), intent(in) :: loan, periods, payment, rate, ic, pf, time
      character(len=*), intent(in) :: unknown
      type(amort_period_result_t) :: out
      real(dp) :: j, n, pmt, principal

      out%status%ok = .false.
      if (ic <= 0.0_dp .or. pf <= 0.0_dp .or. time <= 0.0_dp) goto 900
      j = payment_period_rate(rate, ic, pf)
      select case (trim(unknown))
      case ('loan')
         if (periods <= 0.0_dp .or. payment <= 0.0_dp) goto 900
         principal = payment*(1.0_dp-(1.0_dp+j)**(-periods))/j
         n = periods
         pmt = payment
      case ('payment', 'pmt')
         if (loan <= 0.0_dp .or. periods <= 0.0_dp) goto 900
         principal = loan
         n = periods
         pmt = loan*j/(1.0_dp-(1.0_dp+j)**(-periods))
      case ('periods', 'n')
         if (loan <= 0.0_dp .or. payment <= loan*j) goto 900
         principal = loan
         pmt = payment
         n = -log(1.0_dp-loan*j/payment)/log(1.0_dp+j)
      case default
         out%status%message = 'unknown must be loan, payment, or periods'
         return
      end select
      if (time > n) goto 900
      out%loan = principal
      out%payment = pmt
      out%periods = n
      out%interest_paid = pmt*(1.0_dp-(1.0_dp+j)**(-(n-time+1.0_dp)))
      out%principal_paid = pmt*(1.0_dp+j)**(-(n-time+1.0_dp))
      out%balance = pmt*(1.0_dp-(1.0_dp+j)**(-(n-time)))/j
      out%effective_rate = annual_effective_rate(rate, ic)
      out%status%ok = .true.
      return
900   out%status%message = 'invalid amortization-period inputs'
   end function amort_period

   function bond(face, coupon_rate, redemption, periods, yield_rate, yield_frequency, coupon_frequency, &
      time, has_time) result(out)
      real(dp), intent(in) :: face, coupon_rate, redemption, periods, yield_rate
      real(dp), intent(in) :: yield_frequency, coupon_frequency, time
      logical, intent(in) :: has_time
      type(bond_result_t) :: out
      real(dp) :: j, coupon, rem, d, frac, price_d
      real(dp), allocatable :: idx(:)
      integer :: n, k

      out%status%ok = .false.
      if (face <= 0.0_dp .or. redemption <= 0.0_dp .or. periods <= 0.0_dp .or. &
         yield_rate <= 0.0_dp .or. yield_frequency <= 0.0_dp .or. coupon_frequency <= 0.0_dp) then
         out%status%message = 'bond inputs must be positive'
         return
      end if
      n = nint(periods)
      if (abs(periods-real(n, dp)) > 1.0e-10_dp) then
         out%status%message = 'periods must be an integer'
         return
      end if
      if (has_time .and. (time < 0.0_dp .or. time > periods)) then
         out%status%message = 'time must lie between zero and maturity'
         return
      end if
      j = payment_period_rate(yield_rate, yield_frequency, coupon_frequency)
      coupon = face*coupon_rate/coupon_frequency
      if (abs(j) < 1.0e-14_dp) then
         out%price = coupon*periods+redemption
      else
         out%price = coupon*(1.0_dp-(1.0_dp+j)**(-periods))/j + redemption*(1.0_dp+j)**(-periods)
      end if
      out%coupon = coupon
      out%effective_rate = annual_effective_rate(yield_rate, yield_frequency)
      out%premium = max(0.0_dp, out%price-redemption)
      out%discount = max(0.0_dp, redemption-out%price)
      allocate(idx(n))
      idx = [(real(k, dp), k=1,n)]
      out%macaulay_duration = (sum(coupon*idx*(1.0_dp+j)**(-idx)) + &
         redemption*periods*(1.0_dp+j)**(-periods))/out%price
      out%modified_duration = out%macaulay_duration/(1.0_dp+j)
      out%macaulay_convexity = (sum(coupon*idx*idx*(1.0_dp+j)**(-idx)) + &
         redemption*periods*periods*(1.0_dp+j)**(-periods))/out%price
      out%modified_convexity = (sum(coupon*idx*(idx+1.0_dp)*(1.0_dp+j)**(-idx-2.0_dp)) + &
         redemption*periods*(periods+1.0_dp)*(1.0_dp+j)**(-periods-2.0_dp))/out%price

      if (has_time .and. time > 0.0_dp) then
         d = real(floor(time), dp)
         frac = time-d
         rem = periods-time
         if (abs(j) < 1.0e-14_dp) then
            out%price_at_time = coupon*rem+redemption
         else
            out%price_at_time = coupon*(1.0_dp-(1.0_dp+j)**(-rem))/j + redemption*(1.0_dp+j)**(-rem)
         end if
         if (frac > 1.0e-12_dp) then
            rem = periods-d
            price_d = coupon*(1.0_dp-(1.0_dp+j)**(-rem))/j + redemption*(1.0_dp+j)**(-rem)
            out%full_price = price_d*(1.0_dp+j)**frac
            out%clean_price = out%full_price-frac*coupon
         end if
         out%write_change = abs(face*coupon_rate/coupon_frequency-redemption*j)/ &
            (1.0_dp+j)**(periods-time+1.0_dp)
      end if
      out%status%ok = .true.
   end function bond

end module financialmath_loans
