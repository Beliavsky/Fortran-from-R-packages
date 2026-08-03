! SPDX-License-Identifier: GPL-3.0-or-later
module qbc_bonds
   use qbc_kinds, only : dp
   use qbc_status, only : qbc_success, qbc_invalid_argument, qbc_size_mismatch, qbc_no_convergence
   use qbc_dates, only : qbc_date, add_months, discount_time, year_fraction, operator(>)
   use qbc_types, only : qbc_bond, qbc_curve, qbc_coupon_schedule, qbc_bond_sensitivity, &
      qbc_asset_tes, qbc_rate_continuous
   use qbc_cashflows, only : coupon_dates, coupon_cashflows, accrued_interest
   use qbc_curves, only : discount_factors, discount_factor, curve_rate
   implicit none
   private
   public :: valuation_bonds, dirty_to_clean, price_dirty2clean, bond_price_to_rate, bond_price2rate
   public :: bond_sensitivity, sens_bonds, average_life, accrued_interests

   interface valuation_bonds
      module procedure valuation_bonds_rates
      module procedure valuation_bonds_curve
   end interface valuation_bonds

contains

   subroutine make_bond_cashflows(bond, coupon_rates, schedule, cashflows, previous_coupon, status)
      type(qbc_bond), intent(in) :: bond
      real(dp), intent(in) :: coupon_rates(:)
      type(qbc_coupon_schedule), intent(out) :: schedule
      real(dp), allocatable, intent(out) :: cashflows(:)
      type(qbc_date), intent(out) :: previous_coupon
      integer, intent(out) :: status
      integer :: months_per_coupon
      if (bond%has_trade_date) then
         call coupon_dates(bond%maturity, bond%analysis_date, bond%frequency, bond%coupon_schedule, &
                           bond%business_convention, schedule, trade_date=bond%trade_date, status=status)
      else
         call coupon_dates(bond%maturity, bond%analysis_date, bond%frequency, bond%coupon_schedule, &
                           bond%business_convention, schedule, status=status)
      end if
      if (status /= qbc_success) then
         allocate(cashflows(0)); previous_coupon = bond%analysis_date; return
      end if
      if (bond%frequency > 0 .and. size(schedule%dates) > 0) then
         months_per_coupon = 12 / bond%frequency
         previous_coupon = add_months(schedule%dates(1), -months_per_coupon)
      else
         previous_coupon = bond%analysis_date
      end if
      call coupon_cashflows(schedule%dates, coupon_rates + bond%spread, bond%principal, bond%daycount, &
                            cashflows, previous_coupon, bond%asset_type == qbc_asset_tes, status)
   end subroutine make_bond_cashflows

   real(dp) function valuation_bonds_rates(bond, coupon_rates, rates, dirty, spread_only, status) result(price)
      type(qbc_bond), intent(in) :: bond
      real(dp), intent(in) :: coupon_rates(:), rates(:)
      logical, intent(in), optional :: dirty, spread_only
      integer, intent(out), optional :: status
      type(qbc_coupon_schedule) :: schedule
      type(qbc_date) :: previous_coupon
      real(dp), allocatable :: cashflows(:), factors(:), spread_cf(:)
      logical :: is_dirty, only_spread
      integer :: st
      real(dp) :: accrued
      price = 0.0_dp
      st = qbc_success
      is_dirty = .true.; only_spread = .false.
      if (present(dirty)) is_dirty = dirty
      if (present(spread_only)) only_spread = spread_only
      call make_bond_cashflows(bond, coupon_rates, schedule, cashflows, previous_coupon, st)
      if (st /= qbc_success .or. size(cashflows) == 0) then
         if (present(status)) status = merge(st, qbc_invalid_argument, st /= qbc_success)
         return
      end if
      call discount_factors(schedule%effective_dates, rates, bond%analysis_date, bond%rate_type, bond%frequency, factors, st)
      if (st /= qbc_success) then
         if (present(status)) status = st
         return
      end if
      if (only_spread) then
         allocate(spread_cf(size(cashflows)))
         spread_cf = 0.0_dp
         call spread_cashflows(bond, schedule%dates, previous_coupon, spread_cf)
         price = sum(spread_cf * factors)
      else
         price = sum(cashflows * factors)
         if (.not. is_dirty .and. bond%analysis_date > previous_coupon) then
            accrued = accrued_interest(previous_coupon, schedule%dates(1), bond%analysis_date, &
                                       coupon_rates(1) + bond%spread, bond%principal, bond%daycount)
            price = price - accrued
         end if
      end if
      if (present(status)) status = st
   end function valuation_bonds_rates

   subroutine spread_cashflows(bond, dates, previous_coupon, cashflows)
      type(qbc_bond), intent(in) :: bond
      type(qbc_date), intent(in) :: dates(:), previous_coupon
      real(dp), intent(out) :: cashflows(:)
      type(qbc_date) :: prev
      integer :: i
      prev = previous_coupon
      do i = 1, size(dates)
         cashflows(i) = bond%spread * bond%principal * year_fraction(prev, dates(i), bond%daycount)
         prev = dates(i)
      end do
   end subroutine spread_cashflows

   real(dp) function valuation_bonds_curve(bond, coupon_rates, curve, dirty, spread_only, status) result(price)
      type(qbc_bond), intent(in) :: bond
      real(dp), intent(in) :: coupon_rates(:)
      type(qbc_curve), intent(in) :: curve
      logical, intent(in), optional :: dirty, spread_only
      integer, intent(out), optional :: status
      type(qbc_coupon_schedule) :: schedule
      type(qbc_date) :: previous_coupon
      real(dp), allocatable :: cashflows(:), factors(:), spread_cf(:)
      real(dp) :: term, accrued
      logical :: is_dirty, only_spread
      integer :: i, st
      price = 0.0_dp; st = qbc_success
      is_dirty = .true.; only_spread = .false.
      if (present(dirty)) is_dirty = dirty
      if (present(spread_only)) only_spread = spread_only
      call make_bond_cashflows(bond, coupon_rates, schedule, cashflows, previous_coupon, st)
      if (st /= qbc_success .or. size(cashflows) == 0) then
         if (present(status)) status = st
         return
      end if
      allocate(factors(size(cashflows)))
      do i = 1, size(cashflows)
         term = discount_time(bond%analysis_date, schedule%effective_dates(i))
         factors(i) = discount_factor(curve_rate(curve, term), term, curve%rate_type, curve%frequency)
      end do
      if (only_spread) then
         allocate(spread_cf(size(cashflows)))
         call spread_cashflows(bond, schedule%dates, previous_coupon, spread_cf)
         price = sum(spread_cf * factors)
      else
         price = sum(cashflows * factors)
         if (.not. is_dirty .and. bond%analysis_date > previous_coupon) then
            accrued = accrued_interest(previous_coupon, schedule%dates(1), bond%analysis_date, &
                                       coupon_rates(1) + bond%spread, bond%principal, bond%daycount)
            price = price - accrued
         end if
      end if
      if (present(status)) status = st
   end function valuation_bonds_curve

   real(dp) function dirty_to_clean(bond, dirty_price, coupon_rate, status) result(clean_price)
      type(qbc_bond), intent(in) :: bond
      real(dp), intent(in) :: dirty_price, coupon_rate
      integer, intent(out), optional :: status
      type(qbc_coupon_schedule) :: schedule
      type(qbc_date) :: previous_coupon
      integer :: st, months_per_coupon
      clean_price = dirty_price
      if (bond%has_trade_date) then
         call coupon_dates(bond%maturity, bond%analysis_date, bond%frequency, bond%coupon_schedule, &
                           bond%business_convention, schedule, trade_date=bond%trade_date, status=st)
      else
         call coupon_dates(bond%maturity, bond%analysis_date, bond%frequency, bond%coupon_schedule, &
                           bond%business_convention, schedule, status=st)
      end if
      if (st == qbc_success .and. size(schedule%dates) > 0 .and. bond%frequency > 0) then
         months_per_coupon = 12 / bond%frequency
         previous_coupon = add_months(schedule%dates(1), -months_per_coupon)
         clean_price = dirty_price - accrued_interest(previous_coupon, schedule%dates(1), bond%analysis_date, &
                                                       coupon_rate + bond%spread, bond%principal, bond%daycount)
      end if
      if (present(status)) status = st
   end function dirty_to_clean

   real(dp) function bond_price_to_rate(bond, coupon_rates, target_price, dirty, lower, upper, &
                                        tolerance, max_iterations, status) result(rate)
      type(qbc_bond), intent(in) :: bond
      real(dp), intent(in) :: coupon_rates(:), target_price
      logical, intent(in), optional :: dirty
      real(dp), intent(in), optional :: lower, upper, tolerance
      integer, intent(in), optional :: max_iterations
      integer, intent(out), optional :: status
      real(dp) :: lo, hi, mid, flo, fhi, fmid, tol
      real(dp) :: rvec(1)
      integer :: iter, maxit, st
      logical :: is_dirty
      lo = -0.95_dp; hi = 5.0_dp; tol = 1.0e-10_dp; maxit = 200; is_dirty = .true.
      if (present(lower)) lo = lower
      if (present(upper)) hi = upper
      if (present(tolerance)) tol = tolerance
      if (present(max_iterations)) maxit = max_iterations
      if (present(dirty)) is_dirty = dirty
      rvec(1) = lo
      flo = valuation_bonds_rates(bond, coupon_rates, rvec, is_dirty) - target_price
      rvec(1) = hi
      fhi = valuation_bonds_rates(bond, coupon_rates, rvec, is_dirty) - target_price
      st = qbc_success
      if (flo * fhi > 0.0_dp) then
         st = qbc_invalid_argument
         rate = 0.0_dp
      else
         mid = 0.5_dp * (lo + hi)
         do iter = 1, maxit
            mid = 0.5_dp * (lo + hi)
            rvec(1) = mid
            fmid = valuation_bonds_rates(bond, coupon_rates, rvec, is_dirty) - target_price
            if (abs(fmid) <= tol .or. abs(hi - lo) <= tol * max(1.0_dp, abs(mid))) exit
            if (flo * fmid <= 0.0_dp) then
               hi = mid; fhi = fmid
            else
               lo = mid; flo = fmid
            end if
         end do
         rate = mid
         if (iter > maxit) st = qbc_no_convergence
      end if
      if (present(status)) status = st
   end function bond_price_to_rate

   function bond_sensitivity(bond, coupon_rates, rate, bump, dirty, status) result(sensitivity)
      type(qbc_bond), intent(in) :: bond
      real(dp), intent(in) :: coupon_rates(:), rate
      real(dp), intent(in), optional :: bump
      logical, intent(in), optional :: dirty
      integer, intent(out), optional :: status
      type(qbc_bond_sensitivity) :: sensitivity
      real(dp) :: h, p0, pup, pdown
      real(dp) :: r(1)
      logical :: is_dirty
      integer :: st
      h = 1.0e-4_dp; is_dirty = .true.; st = qbc_success
      if (present(bump)) h = bump
      if (present(dirty)) is_dirty = dirty
      r(1) = rate; p0 = valuation_bonds_rates(bond, coupon_rates, r, is_dirty, status=st)
      r(1) = rate + h; pup = valuation_bonds_rates(bond, coupon_rates, r, is_dirty)
      r(1) = rate - h; pdown = valuation_bonds_rates(bond, coupon_rates, r, is_dirty)
      sensitivity%price = p0
      if (abs(p0) > tiny(1.0_dp) .and. h > 0.0_dp) then
         sensitivity%modified_duration = (pdown - pup) / (2.0_dp * p0 * h)
         sensitivity%convexity = (pdown + pup - 2.0_dp * p0) / (p0 * h * h)
         sensitivity%dv01 = sensitivity%modified_duration * p0 * 1.0e-4_dp
      end if
      if (present(status)) status = st
   end function bond_sensitivity

   real(dp) function average_life(bond, coupon_rates, rates, discounted, status) result(wal)
      type(qbc_bond), intent(in) :: bond
      real(dp), intent(in) :: coupon_rates(:), rates(:)
      logical, intent(in), optional :: discounted
      integer, intent(out), optional :: status
      type(qbc_coupon_schedule) :: schedule
      type(qbc_date) :: previous_coupon
      real(dp), allocatable :: cashflows(:), weights(:), factors(:)
      logical :: use_discount
      integer :: i, st
      use_discount = .false.
      if (present(discounted)) use_discount = discounted
      call make_bond_cashflows(bond, coupon_rates, schedule, cashflows, previous_coupon, st)
      if (st /= qbc_success .or. size(cashflows) == 0) then
         wal = 0.0_dp
      else
         allocate(weights(size(cashflows)))
         weights = cashflows
         if (use_discount) then
            call discount_factors(schedule%effective_dates, rates, bond%analysis_date, bond%rate_type, bond%frequency, factors, st)
            weights = weights * factors
         end if
         wal = 0.0_dp
         do i = 1, size(cashflows)
            wal = wal + discount_time(bond%analysis_date, schedule%effective_dates(i)) * weights(i)
         end do
         if (abs(sum(weights)) > tiny(1.0_dp)) wal = wal / sum(weights)
      end if
      if (present(status)) status = st
   end function average_life

   real(dp) function price_dirty2clean(bond, price, coupon_rate, status) result(clean_price)
      type(qbc_bond), intent(in) :: bond
      real(dp), intent(in) :: price, coupon_rate
      integer, intent(out), optional :: status
      clean_price = dirty_to_clean(bond, price, coupon_rate, status)
   end function price_dirty2clean

   real(dp) function bond_price2rate(bond, coupon_rates, price, dirty, status) result(rate)
      type(qbc_bond), intent(in) :: bond
      real(dp), intent(in) :: coupon_rates(:), price
      logical, intent(in), optional :: dirty
      integer, intent(out), optional :: status
      rate = bond_price_to_rate(bond, coupon_rates, price, dirty=dirty, status=status)
   end function bond_price2rate

   real(dp) function sens_bonds(bond, coupon_rates, rate, status) result(duration)
      type(qbc_bond), intent(in) :: bond
      real(dp), intent(in) :: coupon_rates(:), rate
      integer, intent(out), optional :: status
      type(qbc_bond_sensitivity) :: result
      result = bond_sensitivity(bond, coupon_rates, rate, status=status)
      duration = result%modified_duration
   end function sens_bonds

   real(dp) function accrued_interests(bond, coupon_rate, status) result(value)
      type(qbc_bond), intent(in) :: bond
      real(dp), intent(in) :: coupon_rate
      integer, intent(out), optional :: status
      type(qbc_coupon_schedule) :: schedule
      type(qbc_date) :: previous_coupon
      integer :: st
      value = 0.0_dp
      if (bond%has_trade_date) then
         call coupon_dates(bond%maturity, bond%analysis_date, bond%frequency, bond%coupon_schedule, &
                           bond%business_convention, schedule, trade_date=bond%trade_date, status=st)
      else
         call coupon_dates(bond%maturity, bond%analysis_date, bond%frequency, bond%coupon_schedule, &
                           bond%business_convention, schedule, status=st)
      end if
      if (st == qbc_success .and. size(schedule%dates) > 0 .and. bond%frequency > 0) then
         previous_coupon = add_months(schedule%dates(1), -12 / bond%frequency)
         value = accrued_interest(previous_coupon, schedule%dates(1), bond%analysis_date, &
                                  coupon_rate + bond%spread, bond%principal, bond%daycount)
      end if
      if (present(status)) status = st
   end function accrued_interests

end module qbc_bonds
