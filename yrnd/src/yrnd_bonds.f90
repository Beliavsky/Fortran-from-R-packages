! SPDX-License-Identifier: GPL-3.0-only
module yrnd_bonds
   use yrnd_kinds, only : dp
   use yrnd_dates, only : date_t, add_days, add_months, year_fraction, &
      dc_act_act, operator(<), operator(<=), operator(>), operator(>=), operator(==)
   use tvm, only : xirr
   implicit none
   private

   type, public :: bond_t
      character(len=64) :: id = ""
      real(dp) :: coupon = 0.0_dp
      integer :: coupon_frequency = 2
      type(date_t) :: maturity
      real(dp) :: conversion_factor = 1.0_dp
      real(dp) :: current_yield = 0.0_dp
      real(dp) :: nominal = 100.0_dp
   end type bond_t

   type, public :: bond_context_t
      type(date_t) :: calibration_date
      type(date_t) :: option_date
      type(date_t) :: futures_date
      integer :: day_count = dc_act_act
      integer :: settlement_days = 0
      real(dp) :: option_zero_rate = 0.0_dp
      real(dp) :: futures_zero_rate = 0.0_dp
   end type bond_context_t

   type, public :: ctd_probability_result_t
      character(len=64), allocatable :: bond_id(:)
      real(dp), allocatable :: probability(:)
      integer, allocatable :: ctd_index(:)
   end type ctd_probability_result_t

   public :: accrued_interest, dirty_price_from_yield, yield_from_future_price
   public :: net_basis_from_yield, determine_ctd, ctd_probabilities
   public :: forward_rate_between

contains

   pure real(dp) function forward_rate_between(context) result(rate)
      type(bond_context_t), intent(in) :: context
      real(dp) :: t1, t2
      t1 = year_fraction(context%calibration_date, context%option_date, context%day_count)
      t2 = year_fraction(context%calibration_date, context%futures_date, context%day_count)
      if (t2 <= t1 + epsilon(1.0_dp)) then
         rate = context%futures_zero_rate
      else
         rate = (context%futures_zero_rate * t2 - context%option_zero_rate * t1) / (t2 - t1)
      end if
   end function forward_rate_between

   real(dp) function accrued_interest(bond, date, settlement_days, day_count) result(value)
      type(bond_t), intent(in) :: bond
      type(date_t), intent(in) :: date
      integer, intent(in) :: settlement_days, day_count
      type(date_t) :: previous_coupon, next_coupon, settled_previous
      real(dp) :: fraction
      call surrounding_coupon_dates(bond, date, previous_coupon, next_coupon)
      settled_previous = add_days(previous_coupon, settlement_days)
      if (day_count == dc_act_act) then
         fraction = real(max(0, serial_difference(settled_previous, date)), dp) / &
            real(max(1, serial_difference(previous_coupon, next_coupon)), dp)
      else
         fraction = max(0.0_dp, year_fraction(settled_previous, date, day_count)) * &
            real(bond%coupon_frequency, dp)
      end if
      value = bond%nominal * bond%coupon / real(bond%coupon_frequency, dp) * fraction
   end function accrued_interest

   real(dp) function dirty_price_from_yield(bond, yield, valuation_date, settlement_days, day_count) result(value)
      type(bond_t), intent(in) :: bond
      real(dp), intent(in) :: yield
      type(date_t), intent(in) :: valuation_date
      integer, intent(in) :: settlement_days, day_count
      type(date_t), allocatable :: dates(:)
      real(dp), allocatable :: tau(:), cf(:)
      integer :: i
      call future_coupon_schedule(bond, valuation_date, settlement_days, dates)
      allocate(tau(size(dates)), cf(size(dates)))
      do i = 1, size(dates)
         tau(i) = max(0.0_dp, year_fraction(valuation_date, dates(i), day_count))
         cf(i) = bond%nominal * bond%coupon / real(bond%coupon_frequency, dp)
      end do
      cf(size(cf)) = cf(size(cf)) + bond%nominal
      if (1.0_dp + yield / real(bond%coupon_frequency, dp) <= 0.0_dp) then
         value = huge(1.0_dp)
      else
         value = sum(cf / (1.0_dp + yield / real(bond%coupon_frequency, dp)) ** &
            (tau * real(bond%coupon_frequency, dp)))
      end if
   end function dirty_price_from_yield

   real(dp) function yield_from_future_price(bond, futures_price, context, carry_to_futures, status) result(value)
      type(bond_t), intent(in) :: bond
      real(dp), intent(in) :: futures_price
      type(bond_context_t), intent(in) :: context
      logical, intent(in), optional :: carry_to_futures
      integer, intent(out), optional :: status
      logical :: use_carry
      type(date_t), allocatable :: dates(:)
      real(dp), allocatable :: tau(:), cf(:)
      real(dp) :: dirty, invoice, fwd, t, coupon_fv
      integer :: i, istat

      use_carry = .false.
      if (present(carry_to_futures)) use_carry = carry_to_futures
      invoice = futures_price * bond%conversion_factor + &
         accrued_interest(bond, merge(context%futures_date, context%option_date, use_carry), &
                          context%settlement_days, context%day_count)
      if (use_carry .and. context%futures_date > context%option_date) then
         fwd = forward_rate_between(context)
         t = year_fraction(context%option_date, context%futures_date, context%day_count)
         coupon_fv = coupons_future_value(bond, context%option_date, context%futures_date, &
            context%settlement_days, context%day_count, fwd)
         dirty = (invoice + coupon_fv) * exp(-fwd * t)
      else
         dirty = invoice
      end if

      call future_coupon_schedule(bond, context%option_date, context%settlement_days, dates)
      allocate(tau(size(dates) + 1), cf(size(dates) + 1))
      tau(1) = 0.0_dp
      cf(1) = -dirty
      do i = 1, size(dates)
         tau(i + 1) = max(0.0_dp, year_fraction(context%option_date, dates(i), context%day_count))
         cf(i + 1) = bond%nominal * bond%coupon / real(bond%coupon_frequency, dp)
      end do
      cf(size(cf)) = cf(size(cf)) + bond%nominal
      value = xirr(cf, tau, comp_freq=real(bond%coupon_frequency, dp), status=istat)
      if (present(status)) status = istat
   end function yield_from_future_price

   real(dp) function net_basis_from_yield(bond, yield, futures_price, context, carry_to_futures) result(value)
      type(bond_t), intent(in) :: bond
      real(dp), intent(in) :: yield, futures_price
      type(bond_context_t), intent(in) :: context
      logical, intent(in), optional :: carry_to_futures
      logical :: use_carry
      real(dp) :: dirty, invoice, fwd, t, coupon_fv
      use_carry = .false.
      if (present(carry_to_futures)) use_carry = carry_to_futures
      dirty = dirty_price_from_yield(bond, yield, context%option_date, &
         context%settlement_days, context%day_count)
      if (use_carry .and. context%futures_date > context%option_date) then
         fwd = forward_rate_between(context)
         t = year_fraction(context%option_date, context%futures_date, context%day_count)
         coupon_fv = coupons_future_value(bond, context%option_date, context%futures_date, &
            context%settlement_days, context%day_count, fwd)
         dirty = dirty * exp(fwd * t) - coupon_fv
         invoice = futures_price * bond%conversion_factor + accrued_interest(bond, &
            context%futures_date, context%settlement_days, context%day_count)
      else
         invoice = futures_price * bond%conversion_factor + accrued_interest(bond, &
            context%option_date, context%settlement_days, context%day_count)
      end if
      value = dirty - invoice
   end function net_basis_from_yield

   integer function determine_ctd(bonds, futures_price, context, carry_to_futures) result(index)
      type(bond_t), intent(in) :: bonds(:)
      real(dp), intent(in) :: futures_price
      type(bond_context_t), intent(in) :: context
      logical, intent(in), optional :: carry_to_futures
      logical :: use_carry
      real(dp), allocatable :: basis(:)
      real(dp) :: implied, shift
      integer :: i, k, candidate, status
      use_carry = .false.
      if (present(carry_to_futures)) use_carry = carry_to_futures
      allocate(basis(size(bonds)))
      index = 0
      do k = 1, size(bonds)
         implied = yield_from_future_price(bonds(k), futures_price, context, use_carry, status)
         if (status /= 0) cycle
         shift = implied - bonds(k)%current_yield
         do i = 1, size(bonds)
            basis(i) = net_basis_from_yield(bonds(i), bonds(i)%current_yield + shift, &
               futures_price, context, use_carry)
         end do
         candidate = minloc(basis, dim=1)
         if (candidate == k) then
            index = k
            return
         end if
      end do
      ! Match the R fallback: use the last candidate's parallel yield shift.
      k = size(bonds)
      implied = yield_from_future_price(bonds(k), futures_price, context, use_carry, status)
      shift = implied - bonds(k)%current_yield
      do i = 1, size(bonds)
         basis(i) = net_basis_from_yield(bonds(i), bonds(i)%current_yield + shift, &
            futures_price, context, use_carry)
      end do
      index = minloc(basis, dim=1)
   end function determine_ctd

   subroutine ctd_probabilities(bonds, domain, density, context, result, carry_to_futures)
      type(bond_t), intent(in) :: bonds(:)
      real(dp), intent(in) :: domain(:), density(:)
      type(bond_context_t), intent(in) :: context
      type(ctd_probability_result_t), intent(out) :: result
      logical, intent(in), optional :: carry_to_futures
      logical :: use_carry
      real(dp) :: mass
      integer :: i, k
      if (size(domain) /= size(density) .or. size(domain) < 2) then
         error stop "ctd_probabilities: invalid density grid"
      end if
      use_carry = .false.
      if (present(carry_to_futures)) use_carry = carry_to_futures
      allocate(result%bond_id(size(bonds)), result%probability(size(bonds)))
      allocate(result%ctd_index(size(domain)))
      result%probability = 0.0_dp
      do i = 1, size(bonds)
         result%bond_id(i) = bonds(i)%id
      end do
      do i = 1, size(domain)
         result%ctd_index(i) = determine_ctd(bonds, domain(i), context, use_carry)
      end do
      do i = 1, size(domain) - 1
         mass = 0.5_dp * (density(i) + density(i + 1)) * (domain(i + 1) - domain(i))
         k = determine_ctd(bonds, 0.5_dp * (domain(i) + domain(i + 1)), context, use_carry)
         if (k >= 1 .and. k <= size(bonds)) result%probability(k) = result%probability(k) + mass
      end do
      if (sum(result%probability) > 0.0_dp) result%probability = result%probability / sum(result%probability)
   end subroutine ctd_probabilities

   subroutine surrounding_coupon_dates(bond, date, previous_coupon, next_coupon)
      type(bond_t), intent(in) :: bond
      type(date_t), intent(in) :: date
      type(date_t), intent(out) :: previous_coupon, next_coupon
      integer :: step_months
      step_months = 12 / bond%coupon_frequency
      next_coupon = bond%maturity
      previous_coupon = add_months(next_coupon, -step_months)
      do while (previous_coupon > date)
         next_coupon = previous_coupon
         previous_coupon = add_months(previous_coupon, -step_months)
      end do
      if (next_coupon <= date) then
         previous_coupon = next_coupon
         next_coupon = add_months(next_coupon, step_months)
      end if
   end subroutine surrounding_coupon_dates

   subroutine future_coupon_schedule(bond, valuation_date, settlement_days, dates)
      type(bond_t), intent(in) :: bond
      type(date_t), intent(in) :: valuation_date
      integer, intent(in) :: settlement_days
      type(date_t), allocatable, intent(out) :: dates(:)
      type(date_t) :: previous_coupon, next_coupon, current
      integer :: step_months, n, i
      call surrounding_coupon_dates(bond, valuation_date, previous_coupon, next_coupon)
      step_months = 12 / bond%coupon_frequency
      n = 0
      current = next_coupon
      do while (current <= bond%maturity)
         n = n + 1
         current = add_months(current, step_months)
      end do
      allocate(dates(n))
      current = next_coupon
      do i = 1, n
         if (current == bond%maturity) then
            dates(i) = current
         else
            dates(i) = add_days(current, settlement_days)
         end if
         current = add_months(current, step_months)
      end do
   end subroutine future_coupon_schedule

   real(dp) function coupons_future_value(bond, from_date, to_date, settlement_days, day_count, rate) result(value)
      type(bond_t), intent(in) :: bond
      type(date_t), intent(in) :: from_date, to_date
      integer, intent(in) :: settlement_days, day_count
      real(dp), intent(in) :: rate
      type(date_t), allocatable :: dates(:)
      real(dp) :: t
      integer :: i
      value = 0.0_dp
      call future_coupon_schedule(bond, from_date, settlement_days, dates)
      do i = 1, size(dates)
         if (dates(i) > to_date) cycle
         if (dates(i) == bond%maturity) cycle
         t = max(0.0_dp, year_fraction(dates(i), to_date, day_count))
         value = value + bond%nominal * bond%coupon / real(bond%coupon_frequency, dp) * exp(rate * t)
      end do
   end function coupons_future_value

   pure integer function serial_difference(a, b) result(n)
      use yrnd_dates, only : serial_day
      type(date_t), intent(in) :: a, b
      n = serial_day(b) - serial_day(a)
   end function serial_difference

end module yrnd_bonds
