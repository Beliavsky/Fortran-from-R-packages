! SPDX-License-Identifier: MIT
! Computational translation of credule 0.1.4 by Bertrand Le Nezet (2015).
module credule
   use r_kinds, only : dp
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_value, ieee_quiet_nan
   use credule_brent, only : brent_state, brent_initialize, brent_request, brent_update
   implicit none
   private
   public :: dp, price_cds, bootstrap_cds, discount_factor, survival_probability, cds_legs
   integer, parameter, public :: credule_success = 0
   integer, parameter, public :: credule_invalid_input = 1
   integer, parameter, public :: credule_no_bracket = 2
   integer, parameter, public :: credule_no_convergence = 3
   integer, parameter, public :: credule_numerical_failure = 4

contains

   pure logical function valid_curve(tenors, values, is_credit) result(valid)
      real(dp), intent(in) :: tenors(:) !! Nonempty, positive, strictly increasing curve times in years.
      real(dp), intent(in) :: values(:) !! Matching zero rates or survival probabilities at the knots.
      logical, intent(in) :: is_credit !! True to require nonincreasing survival probabilities in (0,1].
      integer :: n
      valid = .false.
      n = size(tenors)
      if (n == 0 .or. size(values) /= n) return
      if (.not. all(ieee_is_finite(tenors)) .or. .not. all(ieee_is_finite(values))) return
      if (any(tenors <= 0.0_dp)) return
      if (any(tenors(2:n) <= tenors(1:n - 1))) return
      if (is_credit) then
         if (any(values <= 0.0_dp) .or. any(values > 1.0_dp)) return
         if (any(values(2:n) > values(1:n - 1))) return
      end if
      valid = .true.
   end function valid_curve

   pure logical function valid_contract(maturities, recovery, premium_frequency, default_frequency) result(valid)
      real(dp), intent(in) :: maturities(:) !! Positive CDS maturities in years, each spanning at least one interval per leg.
      real(dp), intent(in) :: recovery !! Recovery fraction in [0,1].
      integer, intent(in) :: premium_frequency !! Positive number of premium dates per year.
      integer, intent(in) :: default_frequency !! Positive number of default integration intervals per year.
      real(dp) :: largest_count
      valid = .false.
      if (premium_frequency <= 0 .or. default_frequency <= 0) return
      if (.not. ieee_is_finite(recovery)) return
      if (recovery < 0.0_dp .or. recovery > 1.0_dp) return
      if (.not. all(ieee_is_finite(maturities))) return
      if (any(maturities < 1.0_dp / real(premium_frequency, dp))) return
      if (any(maturities < 1.0_dp / real(default_frequency, dp))) return
      if (size(maturities) > 0) then
         largest_count = maxval(maturities) * real(max(premium_frequency, default_frequency), dp)
         if (.not. ieee_is_finite(largest_count)) return
         if (largest_count >= real(huge(0) - 1, dp)) return
      end if
      valid = .true.
   end function valid_contract

   pure function discount_factor(tenors, rates, time) result(discount)
      real(dp), intent(in) :: tenors(:) !! Positive, increasing yield-curve knot times in years.
      real(dp), intent(in) :: rates(:) !! Continuously compounded annual zero rates, matching tenors; negative rates are allowed.
      real(dp), intent(in) :: time !! Nonnegative evaluation time in years.
      real(dp) :: discount
      discount = ieee_value(0.0_dp, ieee_quiet_nan)
      if (.not. valid_curve(tenors, rates, .false.)) return
      if (.not. ieee_is_finite(time)) return
      if (time < 0.0_dp) return
      discount = discount_raw(tenors, rates, time)
      if (.not. ieee_is_finite(discount)) discount = ieee_value(0.0_dp, ieee_quiet_nan)
   end function discount_factor

   pure function discount_raw(tenors, rates, time) result(discount)
      real(dp), intent(in) :: tenors(:) !! Validated yield-curve times in years.
      real(dp), intent(in) :: rates(:) !! Validated continuously compounded annual rates.
      real(dp), intent(in) :: time !! Nonnegative evaluation time in years.
      real(dp) :: discount, rate
      integer :: i, n
      n = size(tenors)
      rate = rates(n)
      if (time < tenors(1)) then
         rate = rates(1)
      else if (time < tenors(n)) then
         do i = 1, n - 1
            if (time < tenors(i + 1)) then
               rate = rates(i) + (rates(i + 1) - rates(i)) * (time - tenors(i)) / (tenors(i + 1) - tenors(i))
               exit
            end if
         end do
      end if
      discount = exp(-time * rate)
   end function discount_raw

   pure function survival_probability(tenors, probabilities, time) result(survival)
      real(dp), intent(in) :: tenors(:) !! Positive, increasing credit-curve knot times in years.
      real(dp), intent(in) :: probabilities(:) !! Nonincreasing knot survival probabilities in (0,1], matching tenors.
      real(dp), intent(in) :: time !! Nonnegative evaluation time; terminal hazard is extrapolated after the last knot.
      real(dp) :: survival
      survival = ieee_value(0.0_dp, ieee_quiet_nan)
      if (.not. valid_curve(tenors, probabilities, .true.)) return
      if (.not. ieee_is_finite(time)) return
      if (time < 0.0_dp) return
      survival = survival_raw(tenors, probabilities, time)
   end function survival_probability

   pure function survival_raw(tenors, probabilities, time) result(survival)
      real(dp), intent(in) :: tenors(:) !! Validated credit-curve times; may be empty only when time is zero.
      real(dp), intent(in) :: probabilities(:) !! Matching knot survival probabilities.
      real(dp), intent(in) :: time !! Nonnegative time in years.
      real(dp) :: survival, hazard
      integer :: i, n
      if (time == 0.0_dp) then
         survival = 1.0_dp
         return
      end if
      n = size(tenors)
      if (time < tenors(1)) then
         hazard = -log(probabilities(1)) / tenors(1)
         survival = exp(-time * hazard)
      else if (time == tenors(n)) then
         survival = probabilities(n)
      else if (time > tenors(n)) then
         if (n == 1) then
            hazard = -log(probabilities(n)) / tenors(n)
         else
            hazard = -log(probabilities(n) / probabilities(n - 1)) / (tenors(n) - tenors(n - 1))
         end if
         survival = probabilities(n) * exp(-(time - tenors(n)) * hazard)
      else
         survival = probabilities(n)
         do i = 1, n - 1
            if (time < tenors(i + 1)) then
               hazard = -log(probabilities(i + 1) / probabilities(i)) / (tenors(i + 1) - tenors(i))
               survival = probabilities(i) * exp(-(time - tenors(i)) * hazard)
               exit
            end if
         end do
      end if
   end function survival_raw

   pure subroutine leg_values(yield_tenors, yield_rates, credit_tenors, credit_sp, maturity, &
                             premium_frequency, default_frequency, accrued, recovery, tail_hazard, annuity, protection)
      real(dp), intent(in) :: yield_tenors(:) !! Validated yield-curve times in years.
      real(dp), intent(in) :: yield_rates(:) !! Matching continuously compounded zero rates.
      real(dp), intent(in) :: credit_tenors(:) !! Previously known credit times; empty at the first bootstrap knot.
      real(dp), intent(in) :: credit_sp(:) !! Matching survival probabilities at known knots.
      real(dp), intent(in) :: maturity !! CDS maturity in years, with truncated cash-flow grids as in upstream.
      integer, intent(in) :: premium_frequency !! Premium payments per year.
      integer, intent(in) :: default_frequency !! Default integration intervals per year.
      logical, intent(in) :: accrued !! Include the half-period default accrual approximation when true.
      real(dp), intent(in) :: recovery !! Recovery fraction in [0,1].
      real(dp), intent(in) :: tail_hazard !! Constant hazard after each leg's last known grid point.
      real(dp), intent(out) :: annuity !! Unit-spread premium present value, including requested accrual.
      real(dp), intent(out) :: protection !! Protection present value, including loss given default.
      integer :: n, known, count, frequency, leg, j, nc
      real(dp) :: time, previous, delta, anchor, sp, previous_sp, discount, anchor_sp, accrual
      logical :: full_curve
      nc = size(credit_tenors)
      full_curve = .false.
      if (nc > 0) full_curve = maturity <= credit_tenors(nc)
      annuity = 0.0_dp
      accrual = 0.0_dp
      protection = 0.0_dp
      do leg = 1, 2
         if (leg == 1) then
            frequency = premium_frequency
         else
            frequency = default_frequency
         end if
         count = int(maturity * real(frequency, dp))
         known = 0
         if (nc > 0) then
            ! Avoid converting a distant credit knot when the complete curve applies.
            if (full_curve) then
               known = count
            else
               known = int(credit_tenors(nc) * real(frequency, dp))
            end if
         end if
         delta = 1.0_dp / real(frequency, dp)
         anchor = real(known, dp) / real(frequency, dp)
         anchor_sp = 1.0_dp
         if (.not. full_curve .and. nc > 0) anchor_sp = survival_raw(credit_tenors, credit_sp, anchor)
         do j = 1, count
            n = j
            time = real(n, dp) / real(frequency, dp)
            previous = real(n - 1, dp) / real(frequency, dp)
            discount = discount_raw(yield_tenors, yield_rates, time)
            if (n <= known) then
               sp = survival_raw(credit_tenors, credit_sp, time)
               previous_sp = survival_raw(credit_tenors, credit_sp, previous)
            else
               sp = anchor_sp * exp(-tail_hazard * (time - anchor))
               if (previous <= anchor) then
                  previous_sp = survival_raw(credit_tenors, credit_sp, previous)
               else
                  previous_sp = anchor_sp * exp(-tail_hazard * (previous - anchor))
               end if
            end if
            if (leg == 1) then
               annuity = annuity + delta * discount * sp
               if (accrued) accrual = accrual + 0.5_dp * delta * discount * (previous_sp - sp)
            else
               protection = protection + discount * (previous_sp - sp)
            end if
         end do
      end do
      annuity = annuity + accrual
      protection = (1.0_dp - recovery) * protection
   end subroutine leg_values

   pure subroutine cds_legs(yield_tenors, yield_rates, credit_tenors, credit_sp, maturity, recovery, &
                           premium_annuity, default_leg, status, premium_frequency, default_frequency, accrued_premium)
      real(dp), intent(in) :: yield_tenors(:) !! Positive increasing yield-curve knot times in years.
      real(dp), intent(in) :: yield_rates(:) !! Matching continuously compounded annual zero rates.
      real(dp), intent(in) :: credit_tenors(:) !! Positive increasing credit-curve knot times in years.
      real(dp), intent(in) :: credit_sp(:) !! Matching nonincreasing survival probabilities in (0,1].
      real(dp), intent(in) :: maturity !! Positive CDS maturity in years; incomplete intervals are omitted.
      real(dp), intent(in) :: recovery !! Recovery fraction in [0,1].
      real(dp), intent(out) :: premium_annuity !! Unit-spread premium PV; NaN on failure.
      real(dp), intent(out) :: default_leg !! Protection PV per unit notional; NaN on failure.
      integer, intent(out) :: status !! One of the public credule status constants.
      integer, intent(in), optional :: premium_frequency !! Premium payments per year, default 4; must be positive.
      integer, intent(in), optional :: default_frequency !! Default intervals per year, default 12; must be positive.
      logical, intent(in), optional :: accrued_premium !! Include half-period accrual, default true.
      integer :: pf, df
      logical :: accrued
      pf = 4
      df = 12
      accrued = .true.
      if (present(premium_frequency)) pf = premium_frequency
      if (present(default_frequency)) df = default_frequency
      if (present(accrued_premium)) accrued = accrued_premium
      premium_annuity = ieee_value(0.0_dp, ieee_quiet_nan)
      default_leg = premium_annuity
      status = credule_invalid_input
      if (.not. valid_curve(yield_tenors, yield_rates, .false.)) return
      if (.not. valid_curve(credit_tenors, credit_sp, .true.)) return
      if (.not. valid_contract([maturity], recovery, pf, df)) return
      call leg_values(yield_tenors, yield_rates, credit_tenors, credit_sp, maturity, pf, df, accrued, &
                      recovery, 0.0_dp, premium_annuity, default_leg)
      status = credule_success
      if (.not. ieee_is_finite(premium_annuity) .or. .not. ieee_is_finite(default_leg)) then
         status = credule_numerical_failure
      else if (premium_annuity <= 0.0_dp) then
         status = credule_numerical_failure
      end if
      if (status /= credule_success) then
         premium_annuity = ieee_value(0.0_dp, ieee_quiet_nan)
         default_leg = premium_annuity
      end if
   end subroutine cds_legs

   pure subroutine price_cds(yield_tenors, yield_rates, credit_tenors, credit_sp, maturities, recovery, &
                            spreads, status, premium_frequency, default_frequency, accrued_premium)
      real(dp), intent(in) :: yield_tenors(:) !! Positive increasing yield-curve knot times in years.
      real(dp), intent(in) :: yield_rates(:) !! Matching continuously compounded annual zero rates.
      real(dp), intent(in) :: credit_tenors(:) !! Positive increasing credit-curve knot times in years.
      real(dp), intent(in) :: credit_sp(:) !! Matching nonincreasing knot survival probabilities in (0,1].
      real(dp), intent(in) :: maturities(:) !! CDS maturities in years, in arbitrary order; an empty vector is allowed.
      real(dp), intent(in) :: recovery !! Recovery fraction in [0,1].
      real(dp), intent(out) :: spreads(:) !! Annual decimal par spreads, same length as maturities; all NaN on failure.
      integer, intent(out) :: status !! One of the public credule status constants.
      integer, intent(in), optional :: premium_frequency !! Premium payments per year, default 4; must be positive.
      integer, intent(in), optional :: default_frequency !! Default intervals per year, default 12; must be positive.
      logical, intent(in), optional :: accrued_premium !! Include half-period accrual, default true.
      real(dp) :: annuity, protection
      integer :: i, pf, df
      logical :: accrued
      spreads = ieee_value(0.0_dp, ieee_quiet_nan)
      status = credule_invalid_input
      pf = 4
      df = 12
      accrued = .true.
      if (present(premium_frequency)) pf = premium_frequency
      if (present(default_frequency)) df = default_frequency
      if (present(accrued_premium)) accrued = accrued_premium
      if (size(spreads) /= size(maturities)) return
      if (.not. valid_curve(yield_tenors, yield_rates, .false.)) return
      if (.not. valid_curve(credit_tenors, credit_sp, .true.)) return
      if (.not. valid_contract(maturities, recovery, pf, df)) return
      status = credule_success
      do i = 1, size(maturities)
         call leg_values(yield_tenors, yield_rates, credit_tenors, credit_sp, maturities(i), pf, df, &
                         accrued, recovery, 0.0_dp, annuity, protection)
         if (.not. ieee_is_finite(annuity) .or. .not. ieee_is_finite(protection)) then
            status = credule_numerical_failure
            exit
         end if
         if (annuity <= 0.0_dp) then
            status = credule_numerical_failure
            exit
         end if
         spreads(i) = protection / annuity
         if (.not. ieee_is_finite(spreads(i))) then
            status = credule_numerical_failure
            exit
         end if
      end do
      if (status /= credule_success) spreads = ieee_value(0.0_dp, ieee_quiet_nan)
   end subroutine price_cds

   pure subroutine bootstrap_cds(yield_tenors, yield_rates, maturities, spreads, recovery, &
                                survival, hazards, status, premium_frequency, default_frequency, accrued_premium)
      real(dp), intent(in) :: yield_tenors(:) !! Positive increasing yield-curve knot times in years.
      real(dp), intent(in) :: yield_rates(:) !! Matching continuously compounded annual zero rates.
      real(dp), intent(in) :: maturities(:) !! Positive strictly increasing CDS maturities in years; may be empty.
      real(dp), intent(in) :: spreads(:) !! Matching nonnegative annual decimal CDS spreads.
      real(dp), intent(in) :: recovery !! Recovery fraction in [0,1); full recovery cannot identify hazards.
      real(dp), intent(out) :: survival(:) !! Survival probabilities at maturities; all NaN on failure.
      real(dp), intent(out) :: hazards(:) !! Noncumulative interval hazards per year; all NaN on failure.
      integer, intent(out) :: status !! Success, invalid input, no bracket, nonconvergence, or numerical failure.
      integer, intent(in), optional :: premium_frequency !! Premium payments per year, default 4; must be positive.
      integer, intent(in), optional :: default_frequency !! Default intervals per year, default 12; must be positive.
      logical, intent(in), optional :: accrued_premium !! Include half-period accrual, default true.
      integer :: pf, df, i, n, iteration
      logical :: accrued, finished
      type(brent_state) :: solver
      real(dp) :: lower_value, upper_value, value, hazard, annuity, protection, previous_time, previous_sp
      survival = ieee_value(0.0_dp, ieee_quiet_nan)
      hazards = ieee_value(0.0_dp, ieee_quiet_nan)
      status = credule_invalid_input
      n = size(maturities)
      ! Check independent output extents before any conformable array assignment.
      if (size(spreads) /= n .or. size(survival) /= n .or. size(hazards) /= n) return
      pf = 4
      df = 12
      accrued = .true.
      if (present(premium_frequency)) pf = premium_frequency
      if (present(default_frequency)) df = default_frequency
      if (present(accrued_premium)) accrued = accrued_premium
      if (.not. valid_curve(yield_tenors, yield_rates, .false.)) return
      if (.not. valid_contract(maturities, recovery, pf, df)) return
      if (recovery == 1.0_dp) return
      if (any(maturities(2:n) <= maturities(1:n - 1))) return
      if (.not. all(ieee_is_finite(spreads))) return
      if (any(spreads < 0.0_dp)) return
      status = credule_success
      previous_time = 0.0_dp
      previous_sp = 1.0_dp
      do i = 1, n
         call leg_values(yield_tenors, yield_rates, maturities(:i - 1), survival(:i - 1), maturities(i), &
                         pf, df, accrued, recovery, 0.0_dp, annuity, protection)
         lower_value = spreads(i) * annuity - protection
         call leg_values(yield_tenors, yield_rates, maturities(:i - 1), survival(:i - 1), maturities(i), &
                         pf, df, accrued, recovery, 30.0_dp, annuity, protection)
         upper_value = spreads(i) * annuity - protection
         if (.not. ieee_is_finite(lower_value) .or. .not. ieee_is_finite(upper_value)) then
            status = credule_numerical_failure
            exit
         end if
         if (lower_value == 0.0_dp) then
            hazard = 0.0_dp
         else if (upper_value == 0.0_dp) then
            hazard = 30.0_dp
         else if ((lower_value > 0.0_dp .and. upper_value > 0.0_dp) .or. &
                  (lower_value < 0.0_dp .and. upper_value < 0.0_dp)) then
            status = credule_no_bracket
            exit
         else
            call brent_initialize(solver, 0.0_dp, 30.0_dp, lower_value, upper_value)
            finished = .false.
            do iteration = 1, 256
               call brent_request(solver, hazard, finished)
               if (finished) exit
               call leg_values(yield_tenors, yield_rates, maturities(:i - 1), survival(:i - 1), maturities(i), &
                               pf, df, accrued, recovery, hazard, annuity, protection)
               value = spreads(i) * annuity - protection
               if (.not. ieee_is_finite(value)) then
                  status = credule_numerical_failure
                  exit
               end if
               call brent_update(solver, value)
            end do
            if (status /= credule_success) exit
            if (.not. finished) then
               status = credule_no_convergence
               exit
            end if
         end if
         hazards(i) = hazard
         survival(i) = previous_sp * exp(-hazard * (maturities(i) - previous_time))
         if (survival(i) <= 0.0_dp) then
            status = credule_numerical_failure
            exit
         end if
         previous_time = maturities(i)
         previous_sp = survival(i)
      end do
      if (status /= credule_success) then
         survival = ieee_value(0.0_dp, ieee_quiet_nan)
         hazards = ieee_value(0.0_dp, ieee_quiet_nan)
      end if
   end subroutine bootstrap_cds
end module credule
