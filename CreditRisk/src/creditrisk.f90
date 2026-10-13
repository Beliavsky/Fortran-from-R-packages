! CreditRisk 0.1.7 computational translation.
! Copyright (c) 2017 Alessandro Cimarelli & Nicolò Manca. MIT license.
module creditrisk
  use r_kinds, only: dp, i64, r_pi
  use r_distributions, only: r_pnorm
  use r_transforms, only: r_expm1
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_value, ieee_quiet_nan
  implicit none
  private
  public :: dp, credit_curve_t, credit_ok, credit_invalid
  public :: cum_normal_density, generalized_black_scholes, merton, merton_sim, merton_sim_from_shocks
  public :: black_cox, at1p, sbtv, cds, cds2
  public :: linear_constant, valid_times

  integer, parameter :: credit_ok = 0
  integer, parameter :: credit_invalid = 1

  type :: credit_curve_t
    !! Rows are maturities; columns follow the documented R data-frame order.
    real(dp), allocatable :: values(:,:)
    integer :: status = credit_invalid
  end type credit_curve_t

contains

  pure elemental real(dp) function cum_normal_density(x) result(value)
    !! Hastings CDF used by the upstream option-price formula, not exact pnorm.
    real(dp), intent(in) :: x !! Standard-normal variate; infinite tails are supported.
    real(dp) :: k, tail
    if (x > 40.0_dp) then
      value = 1.0_dp
      return
    else if (x < -40.0_dp) then
      value = 0.0_dp
      return
    end if
    k = 1.0_dp / (1.0_dp + 0.2316419_dp * abs(x))
    tail = exp(-x ** 2 / 2.0_dp) / sqrt(2.0_dp * r_pi) * &
      (0.319381530_dp * k - 0.356563782_dp * k ** 2 + 1.781477937_dp * k ** 3 &
      - 1.821255978_dp * k ** 4 + 1.330274429_dp * k ** 5)
    if (x >= 0.0_dp) then
      value = 1.0_dp - tail
    else
      value = tail
    end if
  end function cum_normal_density

  pure elemental real(dp) function generalized_black_scholes(option, spot, strike, time, rate, carry, sigma) result(value)
    character(len=1), intent(in) :: option !! Lowercase c for a call or p for a put.
    real(dp), intent(in) :: spot !! Positive current underlying value.
    real(dp), intent(in) :: strike !! Positive exercise price in the same currency as spot.
    real(dp), intent(in) :: time !! Nonnegative time to expiry in years.
    real(dp), intent(in) :: rate !! Finite annual continuously compounded discount rate.
    real(dp), intent(in) :: carry !! Finite annual cost-of-carry rate.
    real(dp), intent(in) :: sigma !! Nonnegative annual volatility.
    real(dp) :: d1, d2, asset, bond, logasset, logbond
    value = ieee_value(0.0_dp, ieee_quiet_nan)
    if (.not. all(ieee_is_finite([spot, strike, time, rate, carry, sigma]))) return
    if (spot <= 0.0_dp .or. strike <= 0.0_dp .or. time < 0.0_dp .or. sigma < 0.0_dp) return
    if (time == 0.0_dp) then
      if (option == 'c') value = max(spot - strike, 0.0_dp)
      if (option == 'p') value = max(strike - spot, 0.0_dp)
      return
    end if
    logasset = log(spot) + (carry - rate) * time
    logbond = log(strike) - rate * time
    if (.not. all(ieee_is_finite([logasset, logbond]))) return
    if (max(logasset, logbond) > log(huge(1.0_dp))) return
    if (abs((carry - rate) * time) < 500.0_dp) then
      asset = spot * exp((carry - rate) * time)
    else
      asset = exp(logasset)
    end if
    if (abs(rate * time) < 500.0_dp) then
      bond = strike * exp(-rate * time)
    else
      bond = exp(logbond)
    end if
    if (time == 0.0_dp .or. sigma == 0.0_dp) then
      if (option == 'c') value = max(asset - bond, 0.0_dp)
      if (option == 'p') value = max(bond - asset, 0.0_dp)
      return
    end if
    d1 = (log(spot) - log(strike) + (carry + sigma ** 2 / 2.0_dp) * time) / (sigma * sqrt(time))
    d2 = d1 - sigma * sqrt(time)
    if (option == 'c') value = asset * cum_normal_density(d1) - bond * cum_normal_density(d2)
    if (option == 'p') value = bond * cum_normal_density(-d2) - asset * cum_normal_density(-d1)
    if (.not. ieee_is_finite(value)) value = ieee_value(0.0_dp, ieee_quiet_nan)
  end function generalized_black_scholes

  pure function merton(debt, v0, sigma, rate, time) result(curve)
    real(dp), intent(in) :: debt !! Positive debt face value at the last maturity.
    real(dp), intent(in) :: v0 !! Positive firm value at time zero.
    real(dp), intent(in) :: sigma !! Positive annual asset volatility.
    real(dp), intent(in) :: rate !! Finite constant continuously compounded annual rate.
    real(dp), intent(in) :: time(:) !! Nonempty, strictly increasing positive maturities in years.
    type(credit_curve_t) :: curve
    real(dp), allocatable :: firm(:), d2(:), remaining(:)
    call initialize_curve(curve, time, 5)
    if (.not. valid_times(time)) return
    if (.not. all(ieee_is_finite([debt, v0, sigma, rate]))) return
    if (min(debt, v0, sigma) <= 0.0_dp) return
    firm = v0 * exp(rate * time)
    remaining = time(size(time)) - time
    ! Preserve upstream's rate term, which is not multiplied by time.
    d2 = ((log(v0) - log(debt)) + rate - sigma ** 2 * time / 2.0_dp) / (sigma * sqrt(time))
    curve%values(:,2) = firm
    curve%values(:,3) = generalized_black_scholes('c', firm, debt, remaining, rate, rate, sigma)
    curve%values(:,4) = debt * exp(-rate * remaining) &
      - generalized_black_scholes('p', firm, debt, remaining, rate, rate, sigma)
    curve%values(:,5) = r_pnorm(d2)
    call finalize_curve(curve)
  end function merton

  pure function black_cox(debt, v0, sigma, rate, gamma, time, barrier) result(curve)
    real(dp), intent(in) :: debt !! Positive debt face value, used as the final reported barrier.
    real(dp), intent(in) :: v0 !! Positive initial firm value above the discounted barrier.
    real(dp), intent(in) :: sigma !! Positive constant annual asset volatility.
    real(dp), intent(in) :: rate !! Finite constant annual risk-free rate.
    real(dp), intent(in) :: gamma !! Finite annual barrier-discount rate.
    real(dp), intent(in) :: time(:) !! Nonempty, strictly increasing positive maturities in years.
    real(dp), intent(in), optional :: barrier !! Positive K no greater than debt; defaults to debt.
    type(credit_curve_t) :: curve
    real(dp) :: k, h0, drift, exponent, distance
    real(dp), allocatable :: a1(:), logterm(:)
    integer :: n
    call initialize_curve(curve, time, 5)
    if (.not. valid_times(time)) return
    k = debt
    if (present(barrier)) k = barrier
    if (.not. all(ieee_is_finite([debt, v0, sigma, rate, gamma, k]))) return
    if (min(debt, v0, sigma, k) <= 0.0_dp .or. k > debt) return
    n = size(time)
    h0 = k * exp(-gamma * time(n))
    if (h0 <= 0.0_dp .or. h0 >= v0) return
    drift = rate - gamma - 0.5_dp * sigma ** 2
    exponent = 2.0_dp * drift / sigma ** 2
    distance = log(v0) - log(h0)
    a1 = r_pnorm((distance + drift * time) / (sigma * sqrt(time)))
    logterm = -exponent * distance + &
      r_pnorm((-distance + drift * time) / (sigma * sqrt(time)), log_probability=.true.)
    curve%values(:,2) = v0 * exp(rate * time)
    curve%values(:,3) = k * exp(-gamma * (time(n) - time))
    curve%values(n,3) = debt
    curve%values(:,4) = a1 - exp(logterm)
    curve%values(:,5) = survival_intensity(time, curve%values(:,4))
    call finalize_curve(curve)
  end function black_cox

  pure function at1p(v0, h0, shape, sigma, rate, time) result(curve)
    real(dp), intent(in) :: v0 !! Positive initial firm value.
    real(dp), intent(in) :: h0 !! Initial barrier strictly between zero and v0.
    real(dp), intent(in) :: shape !! Nonnegative B controlling barrier shape.
    real(dp), intent(in) :: sigma(:) !! Nonnegative interval volatilities, one per maturity.
    real(dp), intent(in) :: rate(:) !! Finite spot rates, one per maturity; used as exp(rate*time).
    real(dp), intent(in) :: time(:) !! Nonempty, strictly increasing positive maturities in years.
    type(credit_curve_t) :: curve
    real(dp) :: variance, dt, previous, distance, a1, logterm
    integer :: i
    call initialize_curve(curve, time, 5)
    if (.not. valid_times(time)) return
    if (size(sigma) /= size(time) .or. size(rate) /= size(time)) return
    if (.not. all(ieee_is_finite([v0, h0, shape]))) return
    if (.not. all(ieee_is_finite(sigma)) .or. .not. all(ieee_is_finite(rate))) return
    if (v0 <= 0.0_dp .or. h0 <= 0.0_dp .or. h0 >= v0 .or. shape < 0.0_dp .or. any(sigma < 0.0_dp)) return
    variance = 0.0_dp
    previous = 0.0_dp
    distance = log(v0) - log(h0)
    curve%values(:,2) = v0 * exp(rate * time)
    do i = 1, size(time)
      dt = time(i) - previous
      previous = time(i)
      variance = variance + sigma(i) ** 2 * dt
      curve%values(i,3) = (h0 / v0) * curve%values(i,2) * exp(-shape * variance)
      if (variance == 0.0_dp) then
        curve%values(i,4) = 1.0_dp
      else
        a1 = r_pnorm((distance + (shape - 0.5_dp) * variance) / sqrt(variance))
        logterm = -(2.0_dp * shape - 1.0_dp) * distance &
          + r_pnorm((-distance + (shape - 0.5_dp) * variance) / sqrt(variance), log_probability=.true.)
        curve%values(i,4) = a1 - exp(logterm)
      end if
    end do
    curve%values(:,5) = survival_intensity(time, curve%values(:,4))
    call finalize_curve(curve)
  end function at1p

  pure function sbtv(v0, barrier, probability, shape, sigma, rate, time) result(curve)
    real(dp), intent(in) :: v0 !! Positive initial firm value.
    real(dp), intent(in) :: barrier(:) !! Scenario barriers, each strictly between zero and v0.
    real(dp), intent(in) :: probability(:) !! Nonnegative scenario weights summing to one within 1e-12.
    real(dp), intent(in) :: shape !! Nonnegative AT1P barrier-shape parameter B.
    real(dp), intent(in) :: sigma(:) !! Nonnegative interval volatilities, one per maturity.
    real(dp), intent(in) :: rate(:) !! Finite spot interest rates, one per maturity.
    real(dp), intent(in) :: time(:) !! Nonempty, strictly increasing positive maturities in years.
    type(credit_curve_t) :: curve, scenario
    integer :: j
    call initialize_curve(curve, time, 3)
    if (.not. valid_times(time)) return
    if (size(barrier) == 0 .or. size(barrier) /= size(probability)) return
    if (.not. all(ieee_is_finite(probability))) return
    if (any(probability < 0.0_dp) .or. abs(sum(probability) - 1.0_dp) > 1.0e-12_dp) return
    curve%values(:,2) = 0.0_dp
    do j = 1, size(barrier)
      scenario = at1p(v0, barrier(j), shape, sigma, rate, time)
      if (scenario%status /= credit_ok) then
        curve%values(:,2:) = ieee_value(0.0_dp, ieee_quiet_nan)
        return
      end if
      curve%values(:,2) = curve%values(:,2) + probability(j) * scenario%values(:,4)
    end do
    curve%values(:,3) = survival_intensity(time, curve%values(:,2))
    call finalize_curve(curve)
  end function sbtv

  pure function cds(time, intensity, rate, premium, recovery, simplified) result(curve)
    real(dp), intent(in) :: time(:) !! Nonempty, strictly increasing positive payment dates in years.
    real(dp), intent(in) :: intensity(:) !! Nonnegative piecewise-constant hazards, one per payment interval.
    real(dp), intent(in) :: rate(:) !! Finite spot interest rates, length one or number of payment dates.
    real(dp), intent(in), optional :: premium !! Positive annual contractual premium; defaults to 0.005.
    real(dp), intent(in), optional :: recovery !! Bond recovery fraction in [0,1); defaults to 0.4.
    logical, intent(in), optional :: simplified !! True uses discrete legs; false uses the upstream continuous-leg formula.
    type(credit_curve_t) :: curve
    real(dp) :: coupon, rr, dt, prevtime, prevrt, rt, cumulative, q, oldq, discounted, olddiscounted
    real(dp) :: forward, combined, i0, i1, annuity, protection, increment, protection_increment
    logical :: simple
    integer :: i, j
    call initialize_curve(curve, time, 6)
    if (.not. valid_times(time)) return
    if (size(intensity) /= size(time)) return
    if (size(rate) /= 1 .and. size(rate) /= size(time)) return
    if (.not. all(ieee_is_finite(intensity)) .or. .not. all(ieee_is_finite(rate))) return
    if (any(intensity < 0.0_dp)) return
    coupon = 0.005_dp
    rr = 0.4_dp
    simple = .false.
    if (present(premium)) coupon = premium
    if (present(recovery)) rr = recovery
    if (present(simplified)) simple = simplified
    if (.not. all(ieee_is_finite([coupon, rr]))) return
    if (coupon <= 0.0_dp .or. rr < 0.0_dp .or. rr >= 1.0_dp) return
    prevtime = 0.0_dp
    prevrt = 0.0_dp
    cumulative = 0.0_dp
    oldq = 1.0_dp
    olddiscounted = 1.0_dp
    annuity = 0.0_dp
    protection = 0.0_dp
    do i = 1, size(time)
      j = min(i, size(rate))
      dt = time(i) - prevtime
      rt = rate(j) * time(i)
      cumulative = cumulative + intensity(i) * dt
      q = exp(-cumulative)
      discounted = exp(-rt - cumulative)
      if (simple) then
        increment = discounted * dt
        protection_increment = exp(-rt) * (oldq - q)
      else
        forward = (rt - prevrt) / dt
        combined = intensity(i) + forward
        call exponential_integrals(combined, dt, i0, i1)
        ! Upstream subtracts, rather than adds, accrued premium on default.
        increment = discounted * dt - olddiscounted * intensity(i) * i1
        protection_increment = olddiscounted * intensity(i) * i0
      end if
      annuity = annuity + increment
      protection = protection + (1.0_dp - rr) * protection_increment
      curve%values(i,2) = q
      curve%values(i,3) = coupon * annuity
      curve%values(i,4) = protection
      curve%values(i,5) = protection / annuity
      curve%values(i,6) = protection - coupon * annuity
      prevtime = time(i)
      prevrt = rt
      oldq = q
      olddiscounted = discounted
    end do
    call finalize_curve(curve)
  end function cds

  pure function cds2(time, maturity, rate_time, rate, intensity_time, intensity, premium, recovery, simplified) result(curve)
    real(dp), intent(in) :: time(:) !! Strictly increasing positive premium payment grid in years.
    real(dp), intent(in) :: maturity(:) !! Finite output dates; dates outside the payment grid clamp to its endpoints.
    real(dp), intent(in) :: rate_time(:) !! Strictly increasing finite spot-rate knot dates.
    real(dp), intent(in) :: rate(:) !! Finite spot rates with the same length as rate_time.
    real(dp), intent(in) :: intensity_time(:) !! Strictly increasing finite hazard knot dates.
    real(dp), intent(in) :: intensity(:) !! Nonnegative hazard values with the same length as intensity_time.
    real(dp), intent(in), optional :: premium !! Positive contractual premium; defaults to 0.005.
    real(dp), intent(in), optional :: recovery !! Recovery fraction in [0,1); defaults to 0.4.
    logical, intent(in), optional :: simplified !! True selects discrete CDS legs; defaults to false.
    type(credit_curve_t) :: curve, grid
    real(dp), allocatable :: newrate(:), newint(:)
    real(dp) :: coupon
    integer :: j
    call initialize_curve(curve, maturity, 6)
    if (.not. valid_knots(rate_time, rate) .or. .not. valid_knots(intensity_time, intensity)) return
    if (size(maturity) == 0 .or. .not. all(ieee_is_finite(maturity))) return
    newrate = linear_constant(rate_time, rate, time)
    newint = linear_constant(intensity_time, intensity, time)
    grid = cds(time, newint, newrate, premium, recovery, simplified)
    if (grid%status /= credit_ok) return
    do j = 2, 6
      curve%values(:,j) = linear_constant(time, grid%values(:,j), maturity)
    end do
    coupon = 0.005_dp
    if (present(premium)) coupon = premium
    ! Recompute the rate from interpolated legs, as cds2 does in R.
    curve%values(:,5) = coupon * curve%values(:,4) / curve%values(:,3)
    call finalize_curve(curve)
  end function cds2

  pure function merton_sim_from_shocks(v0, rate, sigma, time, shocks) result(curve)
    real(dp), intent(in) :: v0 !! Positive initial firm value.
    real(dp), intent(in) :: rate !! Finite constant annual drift rate.
    real(dp), intent(in) :: sigma !! Nonnegative constant annual volatility.
    real(dp), intent(in) :: time(:) !! Increasing nonnegative dates; an initial zero is allowed.
    real(dp), intent(in) :: shocks(:,:) !! Finite standard-normal draws, shaped (number of dates, number of paths).
    type(credit_curve_t) :: curve
    real(dp) :: previous, dt
    real(dp), allocatable :: logpath(:)
    integer :: i
    allocate(curve%values(size(time),size(shocks,2)))
    curve%values = ieee_value(0.0_dp, ieee_quiet_nan)
    if (.not. valid_times(time, allow_zero=.true.)) return
    if (size(shocks,1) /= size(time) .or. size(shocks,2) == 0) return
    if (.not. all(ieee_is_finite([v0, rate, sigma])) .or. .not. all(ieee_is_finite(shocks))) return
    if (v0 <= 0.0_dp .or. sigma < 0.0_dp) return
    allocate(logpath(size(shocks,2)))
    logpath = 0.0_dp
    previous = 0.0_dp
    do i = 1, size(time)
      dt = time(i) - previous
      logpath = logpath + (rate - sigma ** 2 / 2.0_dp) * dt + sigma * sqrt(dt) * shocks(i,:)
      curve%values(i,:) = v0 * exp(logpath)
      previous = time(i)
    end do
    call finalize_curve(curve)
  end function merton_sim_from_shocks

  pure function merton_sim(v0, rate, sigma, time, paths, seed) result(curve)
    real(dp), intent(in) :: v0 !! Positive initial firm value.
    real(dp), intent(in) :: rate !! Finite constant annual drift rate.
    real(dp), intent(in) :: sigma !! Nonnegative constant annual volatility.
    real(dp), intent(in) :: time(:) !! Increasing nonnegative dates; an initial zero is allowed.
    integer, intent(in) :: paths !! Positive number of simulated trajectories.
    integer, intent(in), optional :: seed !! Deterministic local RNG seed; default 104729, not R's RNG stream.
    type(credit_curve_t) :: curve
    integer(i64) :: state
    real(dp), allocatable :: shocks(:,:)
    real(dp) :: u1, u2
    integer :: i, j
    state = 104729_i64
    if (present(seed)) state = modulo(int(seed, i64), 2147483646_i64) + 1_i64
    allocate(shocks(size(time),max(0,paths)))
    do j = 1, paths
      do i = 1, size(time)
        state = modulo(16807_i64 * state, 2147483647_i64)
        u1 = real(state, dp) / 2147483647.0_dp
        state = modulo(16807_i64 * state, 2147483647_i64)
        u2 = real(state, dp) / 2147483647.0_dp
        shocks(i,j) = sqrt(-2.0_dp * log(u1)) * cos(2.0_dp * r_pi * u2)
      end do
    end do
    curve = merton_sim_from_shocks(v0, rate, sigma, time, shocks)
  end function merton_sim

  pure subroutine finalize_curve(curve)
    type(credit_curve_t), intent(inout) :: curve !! Output table whose status becomes invalid if any value is nonfinite.
    curve%status = credit_invalid
    if (all(ieee_is_finite(curve%values))) curve%status = credit_ok
  end subroutine finalize_curve

  pure subroutine initialize_curve(curve, time, columns)
    type(credit_curve_t), intent(out) :: curve !! Allocated NaN-filled output with invalid status until validated.
    real(dp), intent(in) :: time(:) !! Dates copied into the output's first column.
    integer, intent(in) :: columns !! Positive number of table columns required by the model.
    allocate(curve%values(size(time),columns))
    curve%values = ieee_value(0.0_dp, ieee_quiet_nan)
    curve%values(:,1) = time
    curve%status = credit_invalid
  end subroutine initialize_curve

  pure logical function valid_times(time, allow_zero) result(valid)
    real(dp), intent(in) :: time(:) !! Candidate increasing maturity grid in years.
    logical, intent(in), optional :: allow_zero !! Permit a first date of zero; false by default.
    logical :: zero
    valid = .false.
    zero = .false.
    if (present(allow_zero)) zero = allow_zero
    if (size(time) == 0) return
    if (.not. all(ieee_is_finite(time))) return
    if (time(1) < 0.0_dp) return
    if (.not. zero .and. time(1) == 0.0_dp) return
    if (any(time(2:) <= time(:size(time) - 1))) return
    valid = .true.
  end function valid_times

  pure function survival_intensity(time, survival) result(intensity)
    real(dp), intent(in) :: time(:) !! Strictly increasing positive interval endpoints.
    real(dp), intent(in) :: survival(:) !! Survival at endpoints, matching time in length.
    real(dp) :: intensity(size(time)), previous_time, previous_log, current_log
    integer :: i
    previous_time = 0.0_dp
    previous_log = 0.0_dp
    do i = 1, size(time)
      if (survival(i) <= 0.0_dp .or. .not. ieee_is_finite(survival(i))) then
        intensity = ieee_value(0.0_dp, ieee_quiet_nan)
        return
      end if
      current_log = log(survival(i))
      intensity(i) = (previous_log - current_log) / (time(i) - previous_time)
      previous_time = time(i)
      previous_log = current_log
    end do
  end function survival_intensity

  pure subroutine exponential_integrals(rate, dt, zeroth, first)
    real(dp), intent(in) :: rate !! Combined interval hazard and forward discount rate, possibly zero.
    real(dp), intent(in) :: dt !! Positive interval length in years.
    real(dp), intent(out) :: zeroth !! Integral of exp(-rate*s) over [0,dt].
    real(dp), intent(out) :: first !! Integral of s*exp(-rate*s) over [0,dt].
    real(dp) :: x, term
    integer :: k
    x = rate * dt
    if (abs(x) < 0.01_dp) then
      zeroth = 1.0_dp
      first = 0.5_dp
      term = 1.0_dp
      do k = 1, 10
        term = term * (-x) / real(k, dp)
        zeroth = zeroth + term / real(k + 1, dp)
        first = first + term / real(k + 2, dp)
      end do
      zeroth = dt * zeroth
      first = dt ** 2 * first
    else
      zeroth = -r_expm1(-x) / rate
      first = (-r_expm1(-x) - x * exp(-x)) / rate ** 2
    end if
  end subroutine exponential_integrals

  pure logical function valid_knots(knots, values) result(valid)
    real(dp), intent(in) :: knots(:) !! Strictly increasing finite interpolation knot dates.
    real(dp), intent(in) :: values(:) !! Finite ordinates, one per knot.
    valid = .false.
    if (size(knots) == 0 .or. size(knots) /= size(values)) return
    if (.not. all(ieee_is_finite(knots)) .or. .not. all(ieee_is_finite(values))) return
    if (any(knots(2:) <= knots(:size(knots) - 1))) return
    valid = .true.
  end function valid_knots

  pure function linear_constant(knots, values, targets) result(interpolated)
    real(dp), intent(in) :: knots(:) !! Strictly increasing finite knots; singleton grids are supported.
    real(dp), intent(in) :: values(:) !! Finite ordinates of the same length as knots.
    real(dp), intent(in) :: targets(:) !! Finite evaluation points; extrapolation is constant at endpoints.
    real(dp) :: interpolated(size(targets)), weight
    integer :: i, j, n
    interpolated = ieee_value(0.0_dp, ieee_quiet_nan)
    if (.not. valid_knots(knots, values) .or. .not. all(ieee_is_finite(targets))) return
    n = size(knots)
    do i = 1, size(targets)
      if (targets(i) <= knots(1)) then
        interpolated(i) = values(1)
      else if (targets(i) >= knots(n)) then
        interpolated(i) = values(n)
      else
        do j = 2, n
          if (targets(i) <= knots(j)) exit
        end do
        weight = (targets(i) - knots(j - 1)) / (knots(j) - knots(j - 1))
        interpolated(i) = values(j - 1) + weight * (values(j) - values(j - 1))
      end if
    end do
  end function linear_constant

end module creditrisk
