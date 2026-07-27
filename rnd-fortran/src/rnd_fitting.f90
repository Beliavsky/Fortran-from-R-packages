! SPDX-License-Identifier: GPL-2.0-or-later
! Derived from RND 1.2, Copyright (C) 2017 Kam Hamidieh.
module rnd_fitting
   use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
   use rnd_kinds, only : dp
   use rnd_types, only : bsm_fit, gb_fit, mln_fit, ew_fit, am_fit, shimko_fit, rate_result, moe_result
   use rnd_types, only : option_prices
   use rnd_special, only : normal_cdf
   use rnd_optimize, only : nelder_mead
   use rnd_objectives, only : bsm_objective, gb_objective, mln_objective, ew_objective, mln_am_objective
   use rnd_pricing, only : price_bsm_option
   use rnd_densities, only : dshimko
   use rnd_linalg, only : quadratic_least_squares, simple_linear_regression
   implicit none
   private
   public :: extract_bsm_density, extract_gb_density, extract_mln_density
   public :: extract_ew_density, extract_am_density, extract_shimko_density
   public :: compute_implied_volatility, fit_implied_volatility_curve
   public :: extract_rates, get_point_estimate, fit_all_densities, moe

contains

   function extract_bsm_density(r, dividend_yield, te, s0, market_calls, call_strikes, &
         market_puts, put_strikes, initial_values, call_weights, put_weights, lambda, &
         hessian_flag, max_iter) result(fit)
      real(dp), intent(in) :: r, dividend_yield, te, s0
      real(dp), intent(in) :: market_calls(:), call_strikes(:), market_puts(:), put_strikes(:)
      real(dp), intent(in), optional :: initial_values(:), call_weights(:), put_weights(:), lambda
      logical, intent(in), optional :: hessian_flag
      integer, intent(in), optional :: max_iter
      type(bsm_fit) :: fit
      real(dp), allocatable :: cw(:), pw(:), initial(:)
      real(dp) :: penalty, band, candidate(2), best_value, value, mu, zeta, t
      logical :: hessian
      integer :: i, j, iterations

      call prepare_weights(size(market_calls),call_weights,cw)
      call prepare_weights(size(market_puts),put_weights,pw)
      penalty = 1.0_dp
      if (present(lambda)) penalty = lambda
      hessian = .false.
      if (present(hessian_flag)) hessian = hessian_flag
      iterations = 10000
      if (present(max_iter)) iterations = max_iter
      allocate(initial(2))
      if (present(initial_values)) then
         if (size(initial_values) == 2) initial = initial_values
         if (size(initial_values) /= 2) error stop "extract_bsm_density: initial_values must have size 2"
      else
         band = (r-dividend_yield-0.5_dp*0.3_dp**2)*te
         best_value = huge(1.0_dp)
         do i = 1, 10
            t = real(i-1,dp)/9.0_dp
            mu = (1.0_dp-t)*(log(s0)-band)+t*(log(s0)+band)
            do j = 1, 10
               zeta = sqrt(te)*(0.05_dp+real(j-1,dp)*(0.9_dp-0.05_dp)/9.0_dp)
               candidate = [mu,zeta]
               value = objective(candidate)
               if (value < best_value) then
                  best_value = value
                  initial = candidate
               end if
            end do
         end do
      end if
      fit%optimizer = nelder_mead(objective,initial,max_iter=iterations,compute_hessian=hessian)
      fit%mu = fit%optimizer%par(1)
      fit%zeta = fit%optimizer%par(2)
   contains
      real(dp) function objective(theta) result(value_out)
         real(dp), intent(in) :: theta(:)
         value_out = bsm_objective(theta,s0,r,te,dividend_yield,market_calls,call_strikes,cw, &
            market_puts,put_strikes,pw,penalty)
      end function objective
   end function extract_bsm_density

   function extract_gb_density(r, te, dividend_yield, s0, market_calls, call_strikes, &
         market_puts, put_strikes, initial_values, call_weights, put_weights, lambda, &
         hessian_flag, max_iter) result(fit)
      real(dp), intent(in) :: r, te, dividend_yield, s0
      real(dp), intent(in) :: market_calls(:), call_strikes(:), market_puts(:), put_strikes(:)
      real(dp), intent(in), optional :: initial_values(:), call_weights(:), put_weights(:), lambda
      logical, intent(in), optional :: hessian_flag
      integer, intent(in), optional :: max_iter
      type(gb_fit) :: fit
      real(dp), allocatable :: cw(:), pw(:), initial(:)
      real(dp) :: penalty, candidate(4), best_value, value
      real(dp) :: a, b, v, w
      logical :: hessian
      integer :: ia, ib, iv, iw, iterations

      call prepare_weights(size(market_calls),call_weights,cw)
      call prepare_weights(size(market_puts),put_weights,pw)
      penalty = 1.0_dp
      if (present(lambda)) penalty = lambda
      hessian = .false.
      if (present(hessian_flag)) hessian = hessian_flag
      iterations = 10000
      if (present(max_iter)) iterations = max_iter
      allocate(initial(4))
      if (present(initial_values)) then
         if (size(initial_values) == 4) initial = initial_values
         if (size(initial_values) /= 4) error stop "extract_gb_density: initial_values must have size 4"
      else
         best_value = huge(1.0_dp)
         initial = [2.0_dp,s0,2.0_dp,3.0_dp]
         do ia = 1, 10
            a = real(ia-1,dp)*10.0_dp/9.0_dp
            if (a <= 0.0_dp) cycle
            do ib = 1, 10
               b = (s0-5.0_dp)+real(ib-1,dp)*10.0_dp/9.0_dp
               if (b <= 0.0_dp) cycle
               do iv = 1, 10
                  v = 0.1_dp+real(iv-1,dp)*4.9_dp/9.0_dp
                  do iw = 1, 10
                     w = 0.1_dp+real(iw-1,dp)*4.9_dp/9.0_dp
                     if (a*w <= 4.0_dp) cycle
                     candidate = [a,b,v,w]
                     value = objective(candidate)
                     if (value < best_value) then
                        best_value = value
                        initial = candidate
                     end if
                  end do
               end do
            end do
         end do
      end if
      fit%optimizer = nelder_mead(objective,initial,max_iter=iterations,compute_hessian=hessian)
      fit%a = fit%optimizer%par(1)
      fit%b = fit%optimizer%par(2)
      fit%v = fit%optimizer%par(3)
      fit%w = fit%optimizer%par(4)
   contains
      real(dp) function objective(theta) result(value_out)
         real(dp), intent(in) :: theta(:)
         value_out = gb_objective(theta,r,te,dividend_yield,s0,market_calls,call_strikes,cw, &
            market_puts,put_strikes,pw,penalty)
      end function objective
   end function extract_gb_density

   function extract_mln_density(r, dividend_yield, te, s0, market_calls, call_strikes, &
         market_puts, put_strikes, initial_values, call_weights, put_weights, lambda, &
         hessian_flag, max_iter) result(fit)
      real(dp), intent(in) :: r, dividend_yield, te, s0
      real(dp), intent(in) :: market_calls(:), call_strikes(:), market_puts(:), put_strikes(:)
      real(dp), intent(in), optional :: initial_values(:), call_weights(:), put_weights(:), lambda
      logical, intent(in), optional :: hessian_flag
      integer, intent(in), optional :: max_iter
      type(mln_fit) :: fit
      real(dp), allocatable :: cw(:), pw(:), initial(:)
      real(dp) :: penalty, band, candidate(5), best_value, value
      real(dp) :: alpha1, meanlog1, meanlog2, sdlog1, sdlog2
      logical :: hessian
      integer :: ia, im1, im2, is1, is2, iterations

      call prepare_weights(size(market_calls),call_weights,cw)
      call prepare_weights(size(market_puts),put_weights,pw)
      penalty = 1.0_dp
      if (present(lambda)) penalty = lambda
      hessian = .false.
      if (present(hessian_flag)) hessian = hessian_flag
      iterations = 10000
      if (present(max_iter)) iterations = max_iter
      allocate(initial(5))
      if (present(initial_values)) then
         if (size(initial_values) == 5) initial = initial_values
         if (size(initial_values) /= 5) error stop "extract_mln_density: initial_values must have size 5"
      else
         band = (r-dividend_yield-0.5_dp*0.3_dp**2)*te
         best_value = huge(1.0_dp)
         initial = [0.5_dp,log(s0),log(s0),0.2_dp*sqrt(te),0.4_dp*sqrt(te)]
         do ia = 1, 17
            alpha1 = 0.1_dp+real(ia-1,dp)*0.05_dp
            do im1 = 1, 4
               meanlog1 = log(s0)-band+real(im1-1,dp)*2.0_dp*band/3.0_dp
               do im2 = 1, 4
                  meanlog2 = log(s0)-band+real(im2-1,dp)*2.0_dp*band/3.0_dp
                  do is1 = 1, 7
                     sdlog1 = sqrt(te)*(0.05_dp+real(is1-1,dp)*0.85_dp/6.0_dp)
                     do is2 = 1, 7
                        sdlog2 = sqrt(te)*(0.05_dp+real(is2-1,dp)*0.85_dp/6.0_dp)
                        candidate = [alpha1,meanlog1,meanlog2,sdlog1,sdlog2]
                        value = objective(candidate)
                        if (value < best_value) then
                           best_value = value
                           initial = candidate
                        end if
                     end do
                  end do
               end do
            end do
         end do
      end if
      fit%optimizer = nelder_mead(objective,initial,max_iter=iterations,compute_hessian=hessian)
      fit%alpha1 = fit%optimizer%par(1)
      fit%meanlog1 = fit%optimizer%par(2)
      fit%meanlog2 = fit%optimizer%par(3)
      fit%sdlog1 = fit%optimizer%par(4)
      fit%sdlog2 = fit%optimizer%par(5)
   contains
      real(dp) function objective(theta) result(value_out)
         real(dp), intent(in) :: theta(:)
         value_out = mln_objective(theta,r,dividend_yield,te,s0,market_calls,call_strikes,cw, &
            market_puts,put_strikes,pw,penalty)
      end function objective
   end function extract_mln_density

   function extract_ew_density(r, dividend_yield, te, s0, market_calls, call_strikes, &
         initial_values, call_weights, lambda, hessian_flag, max_iter) result(fit)
      real(dp), intent(in) :: r, dividend_yield, te, s0
      real(dp), intent(in) :: market_calls(:), call_strikes(:)
      real(dp), intent(in), optional :: initial_values(:), call_weights(:), lambda
      logical, intent(in), optional :: hessian_flag
      integer, intent(in), optional :: max_iter
      type(ew_fit) :: fit
      real(dp), allocatable :: cw(:), initial(:)
      real(dp) :: penalty, candidate(3), best_value, value, sigma, v, skew_lognorm, kurt_lognorm
      logical :: hessian
      integer :: i, iterations

      call prepare_weights(size(market_calls),call_weights,cw)
      penalty = 1.0_dp
      if (present(lambda)) penalty = lambda
      hessian = .false.
      if (present(hessian_flag)) hessian = hessian_flag
      iterations = 10000
      if (present(max_iter)) iterations = max_iter
      allocate(initial(3))
      if (present(initial_values)) then
         if (size(initial_values) == 3) initial = initial_values
         if (size(initial_values) /= 3) error stop "extract_ew_density: initial_values must have size 3"
      else
         best_value = huge(1.0_dp)
         do i = 1, 9
            sigma = 0.1_dp*real(i,dp)
            v = sqrt(exp(sigma*sigma*te)-1.0_dp)
            skew_lognorm = 3.0_dp*v+v**3
            kurt_lognorm = 16.0_dp*v**2+15.0_dp*v**4+6.0_dp*v**6+v**8
            candidate = [sigma,skew_lognorm,kurt_lognorm]
            value = objective(candidate)
            if (value < best_value) then
               best_value = value
               initial = candidate
            end if
         end do
      end if
      fit%optimizer = nelder_mead(objective,initial,max_iter=iterations,compute_hessian=hessian)
      fit%sigma = fit%optimizer%par(1)
      fit%skew = fit%optimizer%par(2)
      fit%kurt = fit%optimizer%par(3)
   contains
      real(dp) function objective(theta) result(value_out)
         real(dp), intent(in) :: theta(:)
         value_out = ew_objective(theta,r,dividend_yield,te,s0,market_calls,call_strikes,cw,penalty)
      end function objective
   end function extract_ew_density

   function extract_am_density(r, te, s0, market_calls, market_puts, strikes, initial_values, &
         call_weights, put_weights, lambda, hessian_flag, max_iter) result(fit)
      real(dp), intent(in) :: r, te, s0
      real(dp), intent(in) :: market_calls(:), market_puts(:), strikes(:)
      real(dp), intent(in), optional :: initial_values(:), call_weights(:), put_weights(:), lambda
      logical, intent(in), optional :: hessian_flag
      integer, intent(in), optional :: max_iter
      type(am_fit) :: fit
      real(dp), allocatable :: cw(:), pw(:), initial(:)
      real(dp) :: penalty, band, candidate(10), best_value, value
      real(dp) :: low(10), high(10)
      logical :: hessian
      integer :: mask, j, iterations

      call prepare_weights(size(market_calls),call_weights,cw)
      call prepare_weights(size(market_puts),put_weights,pw)
      penalty = 1.0_dp
      if (present(lambda)) penalty = lambda
      hessian = .false.
      if (present(hessian_flag)) hessian = hessian_flag
      iterations = 10000
      if (present(max_iter)) iterations = max_iter
      allocate(initial(10))
      if (present(initial_values)) then
         if (size(initial_values) == 10) initial = initial_values
         if (size(initial_values) /= 10) error stop "extract_am_density: initial_values must have size 10"
      else
         band = abs((r-0.5_dp*0.3_dp**2)*te)
         low = [0.1_dp,0.1_dp,log(s0)-5.0_dp*band,log(s0)-band,log(s0)+4.0_dp*band, &
            0.1_dp,0.1_dp,0.1_dp,0.1_dp,0.2_dp]
         high = [0.9_dp,0.9_dp,log(s0)-4.0_dp*band,log(s0)+band,log(s0)+5.0_dp*band, &
            0.4_dp,0.4_dp,0.4_dp,0.5_dp,0.6_dp]
         best_value = huge(1.0_dp)
         initial = 0.5_dp*(low+high)
         do mask = 0, 2**10-1
            do j = 1, 10
               if (btest(mask,j-1)) then
                  candidate(j) = high(j)
               else
                  candidate(j) = low(j)
               end if
            end do
            value = objective(candidate)
            if (value < best_value) then
               best_value = value
               initial = candidate
            end if
         end do
      end if
      fit%optimizer = nelder_mead(objective,initial,max_iter=iterations,compute_hessian=hessian)
      fit%w1 = fit%optimizer%par(1)
      fit%w2 = fit%optimizer%par(2)
      fit%u1 = fit%optimizer%par(3)
      fit%u2 = fit%optimizer%par(4)
      fit%u3 = fit%optimizer%par(5)
      fit%sigma1 = fit%optimizer%par(6)
      fit%sigma2 = fit%optimizer%par(7)
      fit%sigma3 = fit%optimizer%par(8)
      fit%p1 = fit%optimizer%par(9)
      fit%p2 = fit%optimizer%par(10)
   contains
      real(dp) function objective(theta) result(value_out)
         real(dp), intent(in) :: theta(:)
         value_out = mln_am_objective(theta,s0,r,te,market_calls,cw,market_puts,pw,strikes,penalty)
      end function objective
   end function extract_am_density

   function compute_implied_volatility(r, te, s0, strike, dividend_yield, call_price, &
         lower, upper, converged) result(sigma)
      real(dp), intent(in) :: r, te, s0, strike(:), dividend_yield, call_price(:), lower, upper
      logical, intent(out), optional :: converged(:)
      real(dp), allocatable :: sigma(:)
      type(option_prices) :: price
      real(dp) :: lo, hi, mid, flo, fhi, fmid, nan_value
      integer :: i, iteration
      logical :: ok
      if (size(strike) /= size(call_price)) error stop "compute_implied_volatility: size mismatch"
      if (present(converged)) then
         if (size(converged) /= size(strike)) error stop "compute_implied_volatility: status size mismatch"
      end if
      allocate(sigma(size(strike)))
      nan_value = ieee_value(0.0_dp,ieee_quiet_nan)
      do i = 1, size(strike)
         lo = max(lower,sqrt(epsilon(1.0_dp)))
         hi = max(upper,lo*2.0_dp)
         price = price_bsm_option(s0,strike(i:i),r,te,lo,dividend_yield)
         flo = price%call(1)-call_price(i)
         price = price_bsm_option(s0,strike(i:i),r,te,hi,dividend_yield)
         fhi = price%call(1)-call_price(i)
         do while (flo*fhi > 0.0_dp .and. hi < 20.0_dp)
            hi = 2.0_dp*hi
            price = price_bsm_option(s0,strike(i:i),r,te,hi,dividend_yield)
            fhi = price%call(1)-call_price(i)
         end do
         ok = flo*fhi <= 0.0_dp
         if (ok) then
            do iteration = 1, 5000
               mid = 0.5_dp*(lo+hi)
               price = price_bsm_option(s0,strike(i:i),r,te,mid,dividend_yield)
               fmid = price%call(1)-call_price(i)
               if (abs(fmid) <= 1.0e-12_dp .or. hi-lo <= 1.0e-12_dp*max(1.0_dp,mid)) exit
               if (flo*fmid <= 0.0_dp) then
                  hi = mid
                  fhi = fmid
               else
                  lo = mid
                  flo = fmid
               end if
            end do
            sigma(i) = 0.5_dp*(lo+hi)
         else
            sigma(i) = nan_value
         end if
         if (present(converged)) converged(i) = ok
      end do
   end function compute_implied_volatility

   subroutine fit_implied_volatility_curve(volatility, strike, a0, a1, a2, info)
      real(dp), intent(in) :: volatility(:), strike(:)
      real(dp), intent(out) :: a0, a1, a2
      integer, intent(out), optional :: info
      real(dp) :: coefficients(3)
      integer :: status
      call quadratic_least_squares(strike,volatility,coefficients,status)
      a0 = coefficients(1)
      a1 = coefficients(2)
      a2 = coefficients(3)
      if (present(info)) info = status
   end subroutine fit_implied_volatility_curve

   function extract_shimko_density(market_calls, call_strikes, r, dividend_yield, te, s0, &
         lower, upper) result(fit)
      real(dp), intent(in) :: market_calls(:), call_strikes(:), r, dividend_yield, te, s0
      real(dp), intent(in), optional :: lower, upper
      type(shimko_fit) :: fit
      real(dp) :: lo, hi
      integer :: info
      lo = 1.0e-6_dp
      if (present(lower)) lo = lower
      hi = 5.0_dp
      if (present(upper)) hi = upper
      allocate(fit%implied_volatility(size(call_strikes)))
      fit%implied_volatility = compute_implied_volatility(r,te,s0,call_strikes,dividend_yield, &
         market_calls,lo,hi)
      call fit_implied_volatility_curve(fit%implied_volatility,call_strikes,fit%a0,fit%a1,fit%a2,info)
      if (info /= 0) error stop "extract_shimko_density: singular implied-volatility regression"
      allocate(fit%density(size(call_strikes)))
      fit%density = dshimko(r,te,s0,call_strikes,dividend_yield,fit%a0,fit%a1,fit%a2)
   end function extract_shimko_density

   function extract_rates(calls, puts, s0, strike, te) result(rates)
      real(dp), intent(in) :: calls(:), puts(:), s0, strike(:), te
      type(rate_result) :: rates
      real(dp) :: intercept, slope
      integer :: info
      call simple_linear_regression(strike,puts-calls,intercept,slope,info)
      if (info /= 0 .or. slope <= 0.0_dp .or. intercept >= 0.0_dp) then
         rates%risk_free_rate = ieee_value(0.0_dp,ieee_quiet_nan)
         rates%dividend_yield = ieee_value(0.0_dp,ieee_quiet_nan)
      else
         rates%risk_free_rate = -log(slope)/te
         rates%dividend_yield = -log(-intercept/s0)/te
      end if
   end function extract_rates

   function get_point_estimate(market_calls, call_strikes, r, te) result(point_estimates)
      real(dp), intent(in) :: market_calls(:), call_strikes(:), r, te
      real(dp), allocatable :: point_estimates(:)
      integer :: i, n
      real(dp) :: factor, right_slope, left_slope
      n = size(call_strikes)
      if (size(market_calls) /= n .or. n < 3) then
         allocate(point_estimates(0))
         return
      end if
      allocate(point_estimates(n-2))
      do i = 2, n-1
         factor = 2.0_dp*exp(r*te)/(call_strikes(i+1)-call_strikes(i-1))
         right_slope = (market_calls(i+1)-market_calls(i))/(call_strikes(i+1)-call_strikes(i))
         left_slope = (market_calls(i-1)-market_calls(i))/(call_strikes(i-1)-call_strikes(i))
         point_estimates(i-1) = factor*(right_slope-left_slope)
      end do
   end function get_point_estimate

   function fit_all_densities(market_calls, call_strikes, market_puts, put_strikes, &
         s0, r, te, dividend_yield, call_weights, put_weights, lambda, max_iter) result(result)
      real(dp), intent(in) :: market_calls(:), call_strikes(:), market_puts(:), put_strikes(:)
      real(dp), intent(in) :: s0, r, te, dividend_yield
      real(dp), intent(in), optional :: call_weights(:), put_weights(:), lambda
      integer, intent(in), optional :: max_iter
      type(moe_result) :: result
      allocate(result%point_density(max(0,size(call_strikes)-2)))
      result%point_density = get_point_estimate(market_calls,call_strikes,r,te)
      result%bsm = extract_bsm_density(r,dividend_yield,te,s0,market_calls,call_strikes, &
         market_puts,put_strikes,call_weights=call_weights,put_weights=put_weights, &
         lambda=lambda,max_iter=max_iter)
      result%gb = extract_gb_density(r,te,dividend_yield,s0,market_calls,call_strikes, &
         market_puts,put_strikes,call_weights=call_weights,put_weights=put_weights, &
         lambda=lambda,max_iter=max_iter)
      result%mln = extract_mln_density(r,dividend_yield,te,s0,market_calls,call_strikes, &
         market_puts,put_strikes,call_weights=call_weights,put_weights=put_weights, &
         lambda=lambda,max_iter=max_iter)
      result%ew = extract_ew_density(r,dividend_yield,te,s0,market_calls,call_strikes, &
         call_weights=call_weights,lambda=lambda,max_iter=max_iter)
      result%shimko = extract_shimko_density(market_calls,call_strikes,r,dividend_yield,te,s0)
   end function fit_all_densities

   function moe(market_calls, call_strikes, market_puts, put_strikes, s0, r, te, &
         dividend_yield, call_weights, put_weights, lambda, max_iter) result(result)
      real(dp), intent(in) :: market_calls(:), call_strikes(:), market_puts(:), put_strikes(:)
      real(dp), intent(in) :: s0, r, te, dividend_yield
      real(dp), intent(in), optional :: call_weights(:), put_weights(:), lambda
      integer, intent(in), optional :: max_iter
      type(moe_result) :: result
      result = fit_all_densities(market_calls,call_strikes,market_puts,put_strikes,s0,r,te, &
         dividend_yield,call_weights,put_weights,lambda,max_iter)
   end function moe

   subroutine prepare_weights(n, provided, weights)
      integer, intent(in) :: n
      real(dp), intent(in), optional :: provided(:)
      real(dp), allocatable, intent(out) :: weights(:)
      allocate(weights(n))
      weights = 1.0_dp
      if (present(provided)) then
         if (size(provided) /= n) error stop "weight vector size mismatch"
         weights = provided
      end if
   end subroutine prepare_weights

end module rnd_fitting
