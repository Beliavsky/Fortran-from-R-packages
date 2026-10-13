! SPDX-License-Identifier: GPL-2.0-or-later
! Literal ESG 1.3 pricing conventions with stable LMN coefficient evaluation.
module esg_pricing
   use esg_models
   use r_distributions, only: r_pnorm
   use r_transforms, only: r_expm1, r_log1p
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_value, ieee_quiet_nan
   implicit none
   private
   public :: esg_zc_bond, esg_bond, esg_stock_call, esg_stock_put, esg_zc_call, esg_zc_put
   public :: esg_corporate_bond, esg_convertible_bond, esg_cds, esg_asset_price, esg_simulate_asset
contains
   pure real(dp) function initial_discount(curve, maturity) result(value)
      real(dp), intent(in) :: curve(:) !! Monthly rates using upstream indices 12*T for a grid maturity T.
      real(dp), intent(in) :: maturity !! Nonnegative maturity in annual model units, within the curve.
      integer :: lower, upper, i
      real(dp) :: grid(0:12), t1, t2, p1, p2
      value = ieee_value(0.0_dp, ieee_quiet_nan)
      if (.not. ieee_is_finite(maturity) .or. maturity < 0.0_dp) return
      if (12.0_dp * maturity > real(size(curve), dp)) return
      if (maturity == 0.0_dp) then
         value = 1.0_dp
         return
      end if
      do i = 0, 12
         grid(i) = real(floor(maturity), dp) + real(i, dp) * (1.0_dp / 12.0_dp)
      end do
      if (any(maturity == grid)) then
         ! R numeric vector subscripts truncate rather than round. Retain that quirk.
         upper = int(12.0_dp * maturity)
         if (upper < 1 .or. upper > size(curve)) return
         value = esg_exp(-curve(upper) * maturity)
      else
         t1 = maxval(grid, mask=grid <= maturity)
         t2 = minval(grid, mask=grid >= maturity)
         lower = int(12.0_dp * t1)
         upper = int(12.0_dp * t2)
         if (upper > size(curve) .or. upper < 1) return
         p1 = 1.0_dp
         if (lower > 0) p1 = esg_exp(-curve(lower) * t1)
         p2 = esg_exp(-curve(upper) * t2)
         value = p1 + (p2 - p1) / (t2 - t1) * (maturity - t1)
      end if
   end function

   pure subroutine price_start(scenario, t, maturity, value, status)
      type(esg_scenario), intent(in) :: scenario !! Model curves and paths needed for valuation.
      integer, intent(in) :: t !! Integer valuation date; source mixes annual path times with monthly curve indices.
      real(dp), intent(in) :: maturity !! Asset maturity at or after t, in the source model's units.
      real(dp), allocatable, intent(out) :: value(:) !! One value at t=0, otherwise one per scenario, initialized to NaN.
      integer, intent(out) :: status !! 0 valid data, 1 invalid domain, 2 missing curve/path or out-of-range date.
      integer :: n
      status = 1
      allocate(value(0))
      if (.not. ieee_is_finite(maturity)) return
      if (t < 0 .or. maturity < real(t, dp)) return
      status = 2
      if (.not. allocated(scenario%zero) .or. .not. allocated(scenario%forward)) return
      if (size(scenario%zero) < 1 .or. size(scenario%forward) < 1) return
      if (t >= size(scenario%forward) .or. maturity >= real(size(scenario%forward), dp)) return
      n = 1
      if (t > 0) then
         if (.not. allocated(scenario%rate)) return
         if (size(scenario%rate, 2) < t + 1) return
         n = size(scenario%rate, 1)
      end if
      if (n < 1) return
      status = 1
      if (.not. all(ieee_is_finite(scenario%zero)) .or. .not. all(ieee_is_finite(scenario%forward))) return
      if (.not. all(ieee_is_finite([scenario%p%k, scenario%p%rate_vol]))) return
      if (scenario%p%k <= 0.0_dp .or. scenario%p%rate_vol < 0.0_dp) return
      if (t > 0) then
         if (.not. all(ieee_is_finite(scenario%rate(:,t + 1)))) return
      end if
      deallocate(value)
      allocate(value(n))
      value = ieee_value(0.0_dp, ieee_quiet_nan)
      status = 0
   end subroutine

   pure subroutine esg_zc_bond(scenario, t, maturity, value, status)
      type(esg_scenario), intent(in) :: scenario !! Hull-White parameters, monthly zero curve, forwards and rate paths.
      integer, intent(in) :: t !! Nonnegative integer valuation date, column t+1 of short rates.
      real(dp), intent(in) :: maturity !! Zero-coupon maturity, possibly off the monthly grid, at least t.
      real(dp), allocatable, intent(out) :: value(:) !! Unit-principal bond distribution; scalar-sized at time zero.
      integer, intent(out) :: status !! 0 success, 1 invalid domain, 2 shape/range error, 4 numerical failure.
      real(dp) :: p0t, p0m, kk
      call price_start(scenario, t, maturity, value, status)
      if (status /= 0) return
      if (real(t, dp) == maturity) then
         value = 1.0_dp
         return
      end if
      p0m = initial_discount(scenario%zero, maturity)
      if (t == 0) then
         value = p0m
      else
         p0t = esg_exp(-scenario%zero(t + 1) * real(t, dp))
         if (p0t <= 0.0_dp .or. .not. ieee_is_finite(p0m)) then
            status = 4
            return
         end if
         kk = esg_k(maturity - real(t, dp), scenario%p%k)
         value = p0m / p0t * esg_exp(-kk ** 2 * esg_l(real(t, dp), scenario%p%rate_vol, scenario%p%k) / 2.0_dp &
            + kk * (scenario%forward(t + 1) - scenario%rate(:,t + 1)))
      end if
      if (.not. all(ieee_is_finite(value))) status = 4
   end subroutine

   pure real(dp) function accrued(coupons, coupon_rate, t, maturity) result(value)
      integer, intent(in) :: coupons !! Positive number of coupons per model year.
      real(dp), intent(in) :: coupon_rate !! Coupon amount per payment, not divided by coupon frequency.
      integer, intent(in) :: t !! Integer valuation date.
      real(dp), intent(in) :: maturity !! Final payment date greater than or equal to t.
      integer :: count
      count = floor(real(coupons, dp) * (maturity - real(t, dp)) + 1.0_dp)
      value = coupon_rate * (real(t, dp) - maturity + real(count, dp) / real(coupons, dp) + 3.0_dp / 360.0_dp)
   end function

   pure subroutine esg_bond(scenario, t, maturity, coupons, coupon_rate, value, status)
      type(esg_scenario), intent(in) :: scenario !! Curves and short-rate paths for unit-principal coupon bonds.
      integer, intent(in) :: t !! Nonnegative integer valuation date.
      real(dp), intent(in) :: maturity !! Final maturity greater than or equal to t.
      integer, intent(in) :: coupons !! Positive coupon frequency per model year, at most 10000.
      real(dp), intent(in) :: coupon_rate !! Payment amount per coupon, nonnegative.
      real(dp), allocatable, intent(out) :: value(:) !! Bond price plus upstream accrued-coupon convention.
      integer, intent(out) :: status !! 0 success, 1 invalid domain, 2 missing paths/range, 4 nonfinite result.
      real(dp), allocatable :: payment(:)
      integer :: i, count
      real(dp) :: date
      call esg_zc_bond(scenario, t, maturity, value, status)
      if (status /= 0) return
      if (coupons < 1 .or. coupons > 10000 .or. .not. ieee_is_finite(coupon_rate)) then
         status = 1
         return
      end if
      if (coupon_rate < 0.0_dp .or. real(coupons, dp) * (maturity - real(t, dp)) > 1.0e6_dp) then
         status = 1
         return
      end if
      if (real(t, dp) == maturity) return
      count = floor(real(coupons, dp) * (maturity - real(t, dp)) + 1.0_dp)
      do i = 0, count - 1
         date = maturity - real(i, dp) / real(coupons, dp)
         call esg_zc_bond(scenario, t, date, payment, status)
         if (status /= 0) return
         value = value + coupon_rate * payment
      end do
      value = value + accrued(coupons, coupon_rate, t, maturity)
      if (.not. all(ieee_is_finite(value))) status = 4
   end subroutine

   pure subroutine stock_option(scenario, t, maturity, strike, is_call, value, status)
      type(esg_scenario), intent(in) :: scenario !! Stock paths and rates; preserves upstream one-column stock lag.
      integer, intent(in) :: t !! Integer valuation date; positive t reads stock column t, not t+1.
      real(dp), intent(in) :: maturity !! Absolute option maturity, also used as Tau argument.
      real(dp), intent(in) :: strike !! Strictly positive exercise price.
      logical, intent(in) :: is_call !! True for call, false for put.
      real(dp), allocatable, intent(out) :: value(:) !! Upstream option distribution; equals one when t=maturity.
      integer, intent(out) :: status !! 0 success, 1 invalid domain, 2 missing paths/range, 4 numerical failure.
      real(dp), allocatable :: spot(:), d1(:), d2(:), bond(:)
      real(dp) :: variance, sd
      call esg_zc_bond(scenario, t, maturity, bond, status)
      allocate(value(size(bond)))
      value = ieee_value(0.0_dp, ieee_quiet_nan)
      if (status /= 0) return
      status = 1
      if (.not. all(ieee_is_finite([strike, scenario%p%stock_vol, scenario%p%rho, scenario%p%stock0]))) return
      if (strike <= 0.0_dp .or. scenario%p%stock_vol < 0.0_dp .or. abs(scenario%p%rho) > 1.0_dp) return
      if (real(t, dp) == maturity) then
         value = 1.0_dp
         status = 0
         return
      end if
      allocate(spot(size(bond)))
      spot = scenario%p%stock0
      if (t > 0) then
         status = 2
         if (.not. allocated(scenario%stock)) return
         if (size(scenario%stock, 1) /= size(spot) .or. size(scenario%stock, 2) < t) return
         spot = scenario%stock(:,t)
      end if
      status = 1
      if (.not. all(ieee_is_finite(spot)) .or. any(spot <= 0.0_dp) .or. any(bond <= 0.0_dp)) return
      variance = esg_tau(maturity, scenario%p%rate_vol, scenario%p%stock_vol, scenario%p%k, scenario%p%rho)
      if (.not. ieee_is_finite(variance) .or. variance < 0.0_dp) return
      if (variance == 0.0_dp) then
         value = max(spot - strike * bond, 0.0_dp)
         if (.not. is_call) value = max(strike * bond - spot, 0.0_dp)
      else
         sd = sqrt(variance)
         d1 = (log(spot) - log(strike) - log(bond) + variance / 2.0_dp) / sd
         d2 = d1 - sd
         if (is_call) then
            value = spot * r_pnorm(d1) - strike * bond * r_pnorm(d2)
         else
            value = -spot * r_pnorm(-d1) + strike * bond * r_pnorm(-d2)
         end if
      end if
      status = 0
      if (.not. all(ieee_is_finite(value))) status = 4
   end subroutine

   pure subroutine zc_option(scenario, t, maturity, underlying, strike, is_call, value, status)
      type(esg_scenario), intent(in) :: scenario !! Hull-White curves and rate paths.
      integer, intent(in) :: t !! Nonnegative integer valuation date, less than option maturity.
      real(dp), intent(in) :: maturity !! Option maturity T, less than or equal to underlying maturity.
      real(dp), intent(in) :: underlying !! Bond maturity s, at least T, within the monthly zero curve.
      real(dp), intent(in) :: strike !! Strictly positive bond strike.
      logical, intent(in) :: is_call !! True for a call, false for a put.
      real(dp), allocatable, intent(out) :: value(:) !! Upstream ZC-option distribution with log(ratio+H/2)/H.
      integer, intent(out) :: status !! 0 success, 1 domain error, 2 missing paths/range, 4 numerical failure.
      real(dp), allocatable :: pshort(:), plong(:), d1(:), d2(:)
      real(dp) :: h, k, variance
      call esg_zc_bond(scenario, t, maturity, pshort, status)
      allocate(value(size(pshort)))
      value = ieee_value(0.0_dp, ieee_quiet_nan)
      if (status /= 0) return
      status = 1
      if (.not. all(ieee_is_finite([underlying, strike]))) return
      if (strike <= 0.0_dp .or. underlying < maturity) return
      if (t == 0 .and. maturity == 0.0_dp) then
         value = 1.0_dp
         status = 0
         return
      end if
      if (maturity <= real(t, dp)) return
      status = 2
      if (underlying >= real(size(scenario%zero), dp) / 12.0_dp) return
      call esg_zc_bond(scenario, t, underlying, plong, status)
      if (status /= 0) return
      k = scenario%p%k
      variance = esg_l(maturity, scenario%p%rate_vol, k) * esg_k(underlying - maturity, k) ** 2
      status = 1
      if (.not. ieee_is_finite(variance) .or. variance < 0.0_dp) return
      if (any(pshort <= 0.0_dp)) return
      if (variance == 0.0_dp) then
         value = max(plong - strike * pshort, 0.0_dp)
         if (.not. is_call) value = max(strike * pshort - plong, 0.0_dp)
      else
         h = sqrt(variance)
         d1 = log(plong / (pshort * strike) + h / 2.0_dp) / h
         d2 = d1 - h
         if (is_call) then
            value = plong * r_pnorm(d1) - strike * pshort * r_pnorm(d2)
         else
            value = -plong * r_pnorm(-d1) + strike * pshort * r_pnorm(-d2)
         end if
      end if
      status = 0
      if (.not. all(ieee_is_finite(value))) status = 4
   end subroutine

   pure subroutine lmn(p, t, a, b, g, h)
      type(esg_parameters), intent(in) :: p !! LMN alpha, beta and strictly positive default volatility.
      real(dp), intent(in) :: t !! Nonnegative elapsed duration used by upstream A/B/G/H.
      real(dp), intent(out) :: a !! LMN coefficient A, evaluated without growing intermediate exponentials.
      real(dp), intent(out) :: b !! LMN coefficient B multiplying current default spread.
      real(dp), intent(out) :: g !! LMN coefficient G for default density.
      real(dp), intent(out) :: h !! LMN coefficient H for default density.
      real(dp) :: phi, decay, q, denominator, one_minus, d, variance, log_a, kk, x, v
      variance = p%default_vol ** 2
      if (variance < tiny(1.0_dp)) then
         ! Continuous zero-diffusion limit when squaring a positive volatility underflows.
         kk = esg_k(t, p%beta)
         if (p%beta * t < 1.0e-5_dp) then
            log_a = -p%alpha * t ** 2 * (0.5_dp - p%beta * t / 6.0_dp + (p%beta * t) ** 2 / 24.0_dp)
         else
            log_a = -p%alpha * (t - kk) / p%beta
         end if
         a = esg_exp(log_a)
         b = -kk
         g = p%alpha * kk * a
         h = a * exp(-p%beta * t)
         return
      end if
      phi = sqrt(2.0_dp * variance + p%beta ** 2)
      decay = exp(-phi * t)
      one_minus = -r_expm1(-phi * t)
      ! Rationalization prevents phi-beta and log(q) from rounding to zero at tiny sigma.
      d = variance / (phi * (phi + p%beta)) * one_minus
      q = 1.0_dp / (1.0_dp - d)
      denominator = 2.0_dp * phi * (1.0_dp - d)
      if (phi * t < 1.0e-3_dp) then
         ! Integrate B'=-1-beta*B+sigma^2*B^2/2 termwise through t^6.
         ! This avoids subtraction of O(1/sigma) terms when both beta and sigma are tiny.
         x = p%beta * t
         v = variance * t ** 2
         log_a = p%alpha * t ** 2 * (-0.5_dp + x / 6.0_dp + (v - x ** 2) / 24.0_dp &
            + x * (x ** 2 - 4.0_dp * v) / 120.0_dp &
            + (-x ** 4 + 11.0_dp * x ** 2 * v - 4.0_dp * v ** 2) / 720.0_dp)
      else
         log_a = -2.0_dp * p%alpha * t / (p%beta + phi) &
            - 2.0_dp * p%alpha * (r_log1p(-d) / variance)
      end if
      a = esg_exp(log_a)
      b = -2.0_dp * one_minus / denominator
      g = p%alpha / phi * one_minus * a * q
      h = a * decay * q ** 2
   end subroutine

   pure subroutine spread_start(scenario, t, maturity, value, lambda, gamma, status)
      type(esg_scenario), intent(in) :: scenario !! Default/liquidity paths and pricing model parameters.
      integer, intent(in) :: t !! Nonnegative integer valuation date.
      real(dp), intent(in) :: maturity !! Asset maturity, at least t.
      real(dp), allocatable, intent(out) :: value(:) !! Result buffer initialized by price_start.
      real(dp), allocatable, intent(out) :: lambda(:) !! Current nonnegative default spreads.
      real(dp), allocatable, intent(out) :: gamma(:) !! Current liquidity spreads, which may be negative.
      integer, intent(out) :: status !! 0 success, 1 invalid parameters, 2 absent/inconsistent paths.
      integer :: n
      call price_start(scenario, t, maturity, value, status)
      if (status /= 0) return
      status = 1
      if (.not. all(ieee_is_finite([scenario%p%alpha, scenario%p%beta, scenario%p%default_vol, scenario%p%eta]))) return
      if (min(scenario%p%alpha, scenario%p%beta, scenario%p%eta) < 0.0_dp) return
      if (scenario%p%default_vol <= 0.0_dp) return
      status = 2
      if (.not. allocated(scenario%default) .or. .not. allocated(scenario%liquidity)) return
      n = size(value)
      if (size(scenario%default, 1) < n .or. size(scenario%default, 2) < t + 1) return
      if (size(scenario%liquidity, 1) < n .or. size(scenario%liquidity, 2) < t + 1) return
      lambda = scenario%default(:n,t + 1)
      gamma = scenario%liquidity(:n,t + 1)
      status = 1
      if (.not. all(ieee_is_finite(lambda)) .or. .not. all(ieee_is_finite(gamma))) return
      if (any(lambda < 0.0_dp)) return
      status = 0
   end subroutine

   pure subroutine esg_cds(scenario, t, maturity, recovery, value, status)
      type(esg_scenario), intent(in) :: scenario !! Curves, rates and LMN default/liquidity state.
      integer, intent(in) :: t !! Integer valuation date, nonnegative.
      real(dp), intent(in) :: maturity !! CDS maturity at or after t.
      real(dp), intent(in) :: recovery !! Upstream omega multiplier in [0,1]; formula uses omega, not 1-omega.
      real(dp), allocatable, intent(out) :: value(:) !! Premium distribution using 100 right-endpoint integration steps.
      integer, intent(out) :: status !! 0 success, 1 domain error, 2 missing paths, 4 numerical failure.
      real(dp), allocatable :: lambda(:), gamma(:), num(:), den(:), bond(:), weight(:)
      real(dp) :: delta, u, a, b, g, h
      integer :: i
      call spread_start(scenario, t, maturity, value, lambda, gamma, status)
      if (status /= 0) return
      status = 1
      if (.not. ieee_is_finite(recovery) .or. recovery < 0.0_dp .or. recovery > 1.0_dp) return
      if (real(t, dp) == maturity) then
         value = 0.0_dp
         status = 0
         return
      end if
      num = 0.0_dp * lambda
      den = num
      delta = (maturity - real(t, dp)) / 100.0_dp
      do i = 1, 100
         u = real(i, dp) * delta
         call esg_zc_bond(scenario, t, real(t, dp) + u, bond, status)
         if (status /= 0) return
         call lmn(scenario%p, u, a, b, g, h)
         weight = delta * esg_exp(b * lambda) * bond
         num = num + weight * (g + lambda * h)
         den = den + weight * a
      end do
      if (any(den <= 0.0_dp)) then
         status = 4
         return
      end if
      value = recovery * num / den
      if (.not. all(ieee_is_finite(value))) status = 4
   end subroutine

   pure subroutine esg_corporate_bond(scenario, t, maturity, coupons, coupon_rate, recovery, value, status)
      type(esg_scenario), intent(in) :: scenario !! Curves, rates and LMN default/liquidity state.
      integer, intent(in) :: t !! Integer valuation date, nonnegative.
      real(dp), intent(in) :: maturity !! Corporate bond maturity at or after t.
      integer, intent(in) :: coupons !! Positive coupon frequency per model year, at most 10000.
      real(dp), intent(in) :: coupon_rate !! Nonnegative payment per coupon.
      real(dp), intent(in) :: recovery !! Upstream omega in [0,1], multiplied as 1-omega in the recovery term.
      real(dp), allocatable, intent(out) :: value(:) !! Upstream corporate bond value, retaining its bare A(T) principal term.
      integer, intent(out) :: status !! 0 success, 1 invalid domain, 2 missing paths/range, 4 numerical failure.
      real(dp), allocatable :: lambda(:), gamma(:), coupon_sum(:), integral(:), bond(:)
      real(dp) :: a, b, g, h, u, delta
      integer :: i, count
      call spread_start(scenario, t, maturity, value, lambda, gamma, status)
      if (status /= 0) return
      status = 1
      if (.not. all(ieee_is_finite([coupon_rate, recovery]))) return
      if (coupons < 1 .or. coupons > 10000 .or. coupon_rate < 0.0_dp) return
      if (recovery < 0.0_dp .or. recovery > 1.0_dp) return
      if (real(coupons, dp) * (maturity - real(t, dp)) > 1.0e6_dp) return
      if (real(t, dp) == maturity) then
         value = 1.0_dp
         status = 0
         return
      end if
      coupon_sum = 0.0_dp * lambda
      integral = coupon_sum
      count = floor(real(coupons, dp) * (maturity - real(t, dp)) + 1.0_dp)
      do i = 0, count - 1
         u = maturity - real(i, dp) / real(coupons, dp)
         call esg_zc_bond(scenario, t, u, bond, status)
         if (status /= 0) return
         call lmn(scenario%p, u, a, b, g, h)
         coupon_sum = coupon_sum + a * bond * esg_exp(b * lambda + scenario%p%eta ** 2 * u ** 3 / 6.0_dp - gamma * u)
      end do
      delta = (maturity - real(t, dp)) / 100.0_dp
      do i = 1, 100
         u = real(i, dp) * delta
         call esg_zc_bond(scenario, t, real(t, dp) + u, bond, status)
         if (status /= 0) return
         call lmn(scenario%p, u, a, b, g, h)
         integral = integral + delta * bond * (g + lambda * h) &
            * esg_exp(b * lambda + scenario%p%eta ** 2 * u ** 3 / 6.0_dp - gamma * u)
      end do
      call lmn(scenario%p, maturity, a, b, g, h)
      value = coupon_rate * coupon_sum + a + (1.0_dp - recovery) * integral &
         + accrued(coupons, coupon_rate, t, maturity)
      if (.not. all(ieee_is_finite(value))) status = 4
   end subroutine

   pure subroutine esg_convertible_bond(scenario, t, maturity, coupons, coupon_rate, value, status)
      type(esg_scenario), intent(in) :: scenario !! Rate and stock paths with model parameters.
      integer, intent(in) :: t !! Integer valuation date, nonnegative.
      real(dp), intent(in) :: maturity !! Maturity at or after t.
      integer, intent(in) :: coupons !! Positive coupon frequency per model year.
      real(dp), intent(in) :: coupon_rate !! Nonnegative coupon payment.
      real(dp), allocatable, intent(out) :: value(:) !! Coupon bond plus stock call with strike one, as in upstream.
      integer, intent(out) :: status !! 0 success; component domain, shape or numerical error otherwise.
      real(dp), allocatable :: option(:)
      call esg_bond(scenario, t, maturity, coupons, coupon_rate, value, status)
      if (status /= 0 .or. real(t, dp) == maturity) return
      call esg_stock_call(scenario, t, maturity, 1.0_dp, option, status)
      if (status == 0) value = value + option
   end subroutine

   pure subroutine esg_stock_call(scenario, t, maturity, strike, value, status)
      type(esg_scenario), intent(in) :: scenario !! Model curves and risk paths required by the option.
      integer, intent(in) :: t !! Nonnegative integer valuation date.
      real(dp), intent(in) :: maturity !! Absolute option maturity at or after the valuation date.
      real(dp), intent(in) :: strike !! Strictly positive option exercise price.
      real(dp), allocatable, intent(out) :: value(:) !! Stock call distribution under upstream conventions.
      integer, intent(out) :: status !! 0 success, 1 invalid domain, 2 missing data, 4 numerical failure.
      call stock_option(scenario, t, maturity, strike, .true., value, status)
   end subroutine

   pure subroutine esg_stock_put(scenario, t, maturity, strike, value, status)
      type(esg_scenario), intent(in) :: scenario !! Model curves and risk paths required by the option.
      integer, intent(in) :: t !! Nonnegative integer valuation date.
      real(dp), intent(in) :: maturity !! Absolute option maturity at or after the valuation date.
      real(dp), intent(in) :: strike !! Strictly positive option exercise price.
      real(dp), allocatable, intent(out) :: value(:) !! Stock put distribution under upstream conventions.
      integer, intent(out) :: status !! 0 success, 1 invalid domain, 2 missing data, 4 numerical failure.
      call stock_option(scenario, t, maturity, strike, .false., value, status)
   end subroutine

   pure subroutine esg_zc_call(scenario, t, maturity, underlying, strike, value, status)
      type(esg_scenario), intent(in) :: scenario !! Model curves and risk paths required by the option.
      integer, intent(in) :: t !! Nonnegative integer valuation date.
      real(dp), intent(in) :: maturity !! Absolute option maturity at or after the valuation date.
      real(dp), intent(in) :: underlying !! Underlying bond maturity, at least the option maturity.
      real(dp), intent(in) :: strike !! Strictly positive option exercise price.
      real(dp), allocatable, intent(out) :: value(:) !! Zc call distribution under upstream conventions.
      integer, intent(out) :: status !! 0 success, 1 invalid domain, 2 missing data, 4 numerical failure.
      call zc_option(scenario, t, maturity, underlying, strike, .true., value, status)
   end subroutine

   pure subroutine esg_zc_put(scenario, t, maturity, underlying, strike, value, status)
      type(esg_scenario), intent(in) :: scenario !! Model curves and risk paths required by the option.
      integer, intent(in) :: t !! Nonnegative integer valuation date.
      real(dp), intent(in) :: maturity !! Absolute option maturity at or after the valuation date.
      real(dp), intent(in) :: underlying !! Underlying bond maturity, at least the option maturity.
      real(dp), intent(in) :: strike !! Strictly positive option exercise price.
      real(dp), allocatable, intent(out) :: value(:) !! Zc put distribution under upstream conventions.
      integer, intent(out) :: status !! 0 success, 1 invalid domain, 2 missing data, 4 numerical failure.
      call zc_option(scenario, t, maturity, underlying, strike, .false., value, status)
   end subroutine

   pure subroutine esg_asset_price(scenario, kind, t, maturity, value, status, &
      coupons, coupon_rate, recovery, underlying, strike)
      type(esg_scenario), intent(in) :: scenario !! Typed curves, model parameters and previously generated paths.
      character(*), intent(in) :: kind !! One of the nine upstream asset type strings; EuroCall_S is also accepted.
      integer, intent(in) :: t !! Nonnegative integer valuation date.
      real(dp), intent(in) :: maturity !! Absolute final or option maturity, at least t.
      real(dp), allocatable, intent(out) :: value(:) !! Requested asset distribution; one element at t=0.
      integer, intent(out) :: status !! 0 success, 1 missing/invalid argument or type, 2 missing data, 4 numerical error.
      integer, intent(in), optional :: coupons !! Required coupon frequency for Bond, CBond and ConvBond.
      real(dp), intent(in), optional :: coupon_rate !! Required payment per coupon for coupon-bearing bonds.
      real(dp), intent(in), optional :: recovery !! Required omega in [0,1] for CBond and CDS.
      real(dp), intent(in), optional :: underlying !! Required underlying bond maturity for ZC options.
      real(dp), intent(in), optional :: strike !! Required positive exercise price for stock and ZC options.
      allocate(value(0))
      status = 1
      select case (kind)
      case ('Zero-Coupon')
         call esg_zc_bond(scenario, t, maturity, value, status)
      case ('Bond', 'CBond', 'ConvBond')
         if (.not. present(coupons) .or. .not. present(coupon_rate)) return
         select case (kind)
         case ('Bond')
            call esg_bond(scenario, t, maturity, coupons, coupon_rate, value, status)
         case ('ConvBond')
            call esg_convertible_bond(scenario, t, maturity, coupons, coupon_rate, value, status)
         case ('CBond')
            if (.not. present(recovery)) return
            call esg_corporate_bond(scenario, t, maturity, coupons, coupon_rate, recovery, value, status)
         end select
      case ('EuroCall_Stock', 'EuroCall_S', 'EuroPut_Stock')
         if (.not. present(strike)) return
         call stock_option(scenario, t, maturity, strike, kind /= 'EuroPut_Stock', value, status)
      case ('EuroCall_ZC', 'EuroPut_ZC')
         if (.not. present(strike) .or. .not. present(underlying)) return
         call zc_option(scenario, t, maturity, underlying, strike, kind == 'EuroCall_ZC', value, status)
      case ('CDS')
         if (.not. present(recovery)) return
         call esg_cds(scenario, t, maturity, recovery, value, status)
      end select
   end subroutine

   pure subroutine esg_simulate_asset(p, zero, innovations, kind, t, maturity, value, status, &
      coupons, coupon_rate, recovery, underlying, strike)
      type(esg_parameters), intent(in) :: p !! Risk parameters used by the requested asset; supplied in a typed record.
      real(dp), intent(in) :: zero(:) !! Monthly zero rates covering the supplied innovation horizon plus one year.
      real(dp), intent(in) :: innovations(:,:,:) !! Explicit draws in (scenario,horizon,risk) order, with five risk columns.
      character(*), intent(in) :: kind !! Upstream asset type string, with EuroCall_S accepted as call alias.
      integer, intent(in) :: t !! Nonnegative integer valuation date.
      real(dp), intent(in) :: maturity !! Asset maturity within the supplied scenario horizon.
      real(dp), allocatable, intent(out) :: value(:) !! Asset price distribution from freshly generated required risk paths.
      integer, intent(out) :: status !! 0 success; generation or valuation status otherwise.
      integer, intent(in), optional :: coupons !! Coupon frequency required for coupon-bearing bonds.
      real(dp), intent(in), optional :: coupon_rate !! Payment per coupon required for coupon-bearing bonds.
      real(dp), intent(in), optional :: recovery !! Omega required for corporate bonds and CDS.
      real(dp), intent(in), optional :: underlying !! Underlying maturity required for bond options.
      real(dp), intent(in), optional :: strike !! Positive exercise price required for options.
      type(esg_scenario) :: scenario
      character(:), allocatable :: risks
      allocate(value(0))
      risks = 'shortRate'
      select case (kind)
      case ('EuroCall_Stock', 'EuroCall_S', 'EuroPut_Stock', 'ConvBond')
         risks = 'stock'
      end select
      call esg_generate(p, zero, innovations, scenario, status, risks)
      if (status /= 0) return
      if (kind == 'CBond' .or. kind == 'CDS') then
         call esg_liquidity(p, innovations(:,:,4), scenario%liquidity, status)
         if (status /= 0) return
         call esg_default(p, innovations(:,:,5), scenario%default, status)
         if (status /= 0) return
      end if
      call esg_asset_price(scenario, kind, t, maturity, value, status, &
         coupons, coupon_rate, recovery, underlying, strike)
   end subroutine
end module
