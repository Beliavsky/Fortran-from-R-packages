! SPDX-License-Identifier: GPL-2.0-or-later
module ragtop_calibration
   use ragtop_kinds, only : dp
   use ragtop_constants
   use ragtop_types, only : market_spec, instrument_spec, dividend_schedule, &
                            volatility_curve, option_value
   use ragtop_black_scholes, only : black_scholes, &
                                     black_scholes_on_term_structures
   use ragtop_pde, only : american, find_present_value
   use ragtop_term_structures, only : initialize_volatility_curve
   implicit none
   private
   public :: implied_volatility, implied_volatilities
   public :: implied_volatility_with_term_struct
   public :: american_implied_volatility
   public :: equivalent_jump_vola_to_bs, equivalent_bs_vola_to_jump
   public :: implied_jump_process_volatility, fit_variance_cumulation
   public :: price_with_intensity_link, penalty_with_intensity_link
   public :: fit_to_option_market

contains

   function eval_black_scholes(callput, spot, strike, rate, maturity, &
                               volatility, hazard, div_rate, borrow, &
                               dividends) result(ans)
      integer, intent(in) :: callput
      real(dp), intent(in) :: spot, strike, rate, maturity
      real(dp), intent(in) :: volatility, hazard, div_rate, borrow
      type(dividend_schedule), intent(in), optional :: dividends
      type(option_value) :: ans
      if (present(dividends)) then
         ans = black_scholes(callput, spot, strike, rate, maturity, &
                             volatility, hazard, div_rate, borrow, dividends)
      else
         ans = black_scholes(callput, spot, strike, rate, maturity, &
                             volatility, hazard, div_rate, borrow)
      end if
   end function eval_black_scholes

   function implied_volatility(option_price, callput, spot, strike, rate, &
                               maturity, default_intensity, dividend_rate, &
                               borrow_cost, dividends, relative_tolerance, &
                               max_iter, max_volatility, status) result(vol)
      real(dp), intent(in) :: option_price, spot, strike, rate, maturity
      integer, intent(in) :: callput
      real(dp), intent(in), optional :: default_intensity
      real(dp), intent(in), optional :: dividend_rate, borrow_cost
      type(dividend_schedule), intent(in), optional :: dividends
      real(dp), intent(in), optional :: relative_tolerance, max_volatility
      integer, intent(in), optional :: max_iter
      integer, intent(out), optional :: status
      real(dp) :: vol, lo, hi, tol, h, q, b, trial
      integer :: i, nmax, st
      type(option_value) :: low_value, high_value, current

      h = 0.0_dp
      q = 0.0_dp
      b = 0.0_dp
      if (present(default_intensity)) h = default_intensity
      if (present(dividend_rate)) q = dividend_rate
      if (present(borrow_cost)) b = borrow_cost
      tol = 1.0e-6_dp
      if (present(relative_tolerance)) tol = relative_tolerance
      nmax = 100
      if (present(max_iter)) nmax = max_iter
      lo = 1.0e-12_dp
      hi = 4.0_dp
      if (present(max_volatility)) hi = max_volatility
      st = ragtop_ok

      low_value = eval_black_scholes(callput, spot, strike, rate, maturity, &
                                     lo, h, q, b, dividends)
      high_value = eval_black_scholes(callput, spot, strike, rate, maturity, &
                                      hi, h, q, b, dividends)
      if (option_price < low_value%price .or. &
          option_price > high_value%price) then
         vol = 0.0_dp
         st = ragtop_no_solution
         if (present(status)) status = st
         return
      end if

      vol = min(max(0.5_dp, lo), hi)
      do i = 1, nmax
         current = eval_black_scholes(callput, spot, strike, rate, maturity, &
                                      vol, h, q, b, dividends)
         if (abs(current%price-option_price) <= &
             tol*max(1.0_dp, abs(option_price))) exit
         if (current%price > option_price) then
            hi = vol
         else
            lo = vol
         end if
         trial = vol
         if (current%vega > sqrt(epsilon(1.0_dp))) then
            trial = vol-(current%price-option_price)/current%vega
         end if
         if (trial <= lo .or. trial >= hi) trial = 0.5_dp*(lo+hi)
         vol = trial
      end do
      if (i > nmax) st = ragtop_max_iterations
      if (present(status)) status = st
   end function implied_volatility

   subroutine implied_volatilities(option_price, callput, spot, strike, &
                                   rate, maturity, vol, status)
      real(dp), intent(in) :: option_price(:), spot(:), strike(:)
      real(dp), intent(in) :: rate(:), maturity(:)
      integer, intent(in) :: callput(:)
      real(dp), intent(out) :: vol(:)
      integer, intent(out), optional :: status(:)
      integer :: i, st
      do i = 1, size(vol)
         vol(i) = implied_volatility(option_price(i), callput(i), spot(i), &
                                     strike(i), rate(i), maturity(i), &
                                     status=st)
         if (present(status)) status(i) = st
      end do
   end subroutine implied_volatilities

   function implied_volatility_with_term_struct(option_price, callput, spot, &
                                                strike, maturity, market, &
                                                dividends, relative_tolerance, &
                                                max_iter, max_volatility, &
                                                status) result(vol)
      real(dp), intent(in) :: option_price, spot, strike, maturity
      integer, intent(in) :: callput
      type(market_spec), intent(in) :: market
      type(dividend_schedule), intent(in), optional :: dividends
      real(dp), intent(in), optional :: relative_tolerance, max_volatility
      integer, intent(in), optional :: max_iter
      integer, intent(out), optional :: status
      real(dp) :: vol, lo, hi, mid, tol
      type(market_spec) :: trial_market
      type(option_value) :: current
      integer :: i, nmax, st

      tol = 1.0e-6_dp
      if (present(relative_tolerance)) tol = relative_tolerance
      nmax = 100
      if (present(max_iter)) nmax = max_iter
      lo = 1.0e-12_dp
      hi = 4.0_dp
      if (present(max_volatility)) hi = max_volatility
      st = ragtop_ok
      do i = 1, nmax
         mid = 0.5_dp*(lo+hi)
         trial_market = market
         trial_market%volatility = mid
         trial_market%use_vol_curve = .false.
         if (present(dividends)) then
            current = black_scholes_on_term_structures(callput, spot, strike, &
                                                        maturity, trial_market, &
                                                        dividends)
         else
            current = black_scholes_on_term_structures(callput, spot, strike, &
                                                        maturity, trial_market)
         end if
         if (abs(current%price-option_price) <= &
             tol*max(1.0_dp, abs(option_price))) exit
         if (current%price > option_price) then
            hi = mid
         else
            lo = mid
         end if
      end do
      vol = 0.5_dp*(lo+hi)
      if (i > nmax) st = ragtop_max_iterations
      if (present(status)) status = st
   end function implied_volatility_with_term_struct

   function american_implied_volatility(option_price, callput, spot, strike, &
                                        maturity, market, min_steps, &
                                        relative_tolerance, max_iter, &
                                        max_volatility, status) result(vol)
      real(dp), intent(in) :: option_price, spot, strike, maturity
      integer, intent(in) :: callput
      type(market_spec), intent(in) :: market
      integer, intent(in), optional :: min_steps, max_iter
      real(dp), intent(in), optional :: relative_tolerance, max_volatility
      integer, intent(out), optional :: status
      type(market_spec) :: trial_market
      real(dp) :: vol, lo, hi, mid, px, tol
      integer :: i, ns, nmax, st

      ns = 50
      if (present(min_steps)) ns = min_steps
      nmax = 100
      if (present(max_iter)) nmax = max_iter
      tol = 1.0e-4_dp
      if (present(relative_tolerance)) tol = relative_tolerance
      lo = 1.0e-6_dp
      hi = 4.0_dp
      if (present(max_volatility)) hi = max_volatility
      st = ragtop_ok
      do i = 1, nmax
         mid = 0.5_dp*(lo+hi)
         trial_market = market
         trial_market%volatility = mid
         trial_market%use_vol_curve = .false.
         px = american(callput, spot, strike, maturity, trial_market, ns)
         if (abs(px-option_price) <= &
             tol*max(1.0_dp, abs(option_price))) exit
         if (px > option_price) then
            hi = mid
         else
            lo = mid
         end if
      end do
      vol = 0.5_dp*(lo+hi)
      if (i > nmax) st = ragtop_max_iterations
      if (present(status)) status = st
   end function american_implied_volatility

   function equivalent_jump_vola_to_bs(bs_volatility, maturity, market, &
                                       status) result(jump_volatility)
      real(dp), intent(in) :: bs_volatility, maturity
      type(market_spec), intent(in) :: market
      integer, intent(out), optional :: status
      type(market_spec) :: no_default
      type(option_value) :: target_value
      real(dp) :: jump_volatility

      no_default = market
      no_default%default_intensity = 0.0_dp
      no_default%use_hazard_link = .false.
      no_default%volatility = bs_volatility
      no_default%use_vol_curve = .false.
      target_value = black_scholes_on_term_structures(put_option, 1.0_dp, &
                                                       1.0_dp, maturity, &
                                                       no_default)
      jump_volatility = implied_volatility_with_term_struct( &
         target_value%price, put_option, 1.0_dp, 1.0_dp, maturity, market, &
         status=status)
   end function equivalent_jump_vola_to_bs

   function equivalent_bs_vola_to_jump(jump_volatility, maturity, market, &
                                       status) result(bs_volatility)
      real(dp), intent(in) :: jump_volatility, maturity
      type(market_spec), intent(in) :: market
      integer, intent(out), optional :: status
      type(market_spec) :: jump_market, no_default
      type(option_value) :: target_value
      real(dp) :: bs_volatility

      jump_market = market
      jump_market%volatility = jump_volatility
      jump_market%use_vol_curve = .false.
      target_value = black_scholes_on_term_structures(put_option, 1.0_dp, &
                                                       1.0_dp, maturity, &
                                                       jump_market)
      no_default = market
      no_default%default_intensity = 0.0_dp
      no_default%use_hazard_link = .false.
      bs_volatility = implied_volatility_with_term_struct( &
         target_value%price, put_option, 1.0_dp, 1.0_dp, maturity, &
         no_default, status=status)
   end function equivalent_bs_vola_to_jump

   function implied_jump_process_volatility(instrument_price, spot, inst, &
                                            market, min_steps, status) &
                                            result(vol)
      real(dp), intent(in) :: instrument_price, spot
      type(instrument_spec), intent(in) :: inst
      type(market_spec), intent(in) :: market
      integer, intent(in), optional :: min_steps
      integer, intent(out), optional :: status
      type(market_spec) :: trial_market
      real(dp) :: vol, lo, hi, mid, px
      integer :: i, ns, st

      ns = 75
      if (present(min_steps)) ns = min_steps
      lo = 1.0e-6_dp
      hi = 4.0_dp
      st = ragtop_ok
      do i = 1, 100
         mid = 0.5_dp*(lo+hi)
         trial_market = market
         trial_market%volatility = mid
         trial_market%use_vol_curve = .false.
         px = find_present_value(spot, ns, inst, trial_market)
         if (abs(px-instrument_price) <= &
             1.0e-5_dp*max(1.0_dp, abs(instrument_price))) exit
         if (px > instrument_price) then
            hi = mid
         else
            lo = mid
         end if
      end do
      vol = 0.5_dp*(lo+hi)
      if (i > 100) st = ragtop_max_iterations
      if (present(status)) status = st
   end function implied_jump_process_volatility

   subroutine fit_variance_cumulation(spot, options, prices, market, curve, &
                                      status)
      real(dp), intent(in) :: spot
      type(instrument_spec), intent(in) :: options(:)
      real(dp), intent(in) :: prices(:)
      type(market_spec), intent(in) :: market
      type(volatility_curve), intent(out) :: curve
      integer, intent(out), optional :: status
      real(dp), allocatable :: time(:), volatility(:)
      integer :: i, st, overall

      overall = ragtop_ok
      allocate(time(size(options)), volatility(size(options)))
      do i = 1, size(options)
         time(i) = options(i)%maturity
         volatility(i) = implied_volatility_with_term_struct( &
            prices(i), options(i)%callput, spot, options(i)%strike, time(i), &
            market, status=st)
         if (st /= ragtop_ok) overall = st
      end do
      call initialize_volatility_curve(curve, time, volatility, st)
      if (st /= ragtop_ok) overall = st
      if (present(status)) status = overall
   end subroutine fit_variance_cumulation

   function price_with_intensity_link(power, constant_fraction, base_hazard, &
                                      spot, inst, market, min_steps) &
                                      result(price)
      real(dp), intent(in) :: power, constant_fraction, base_hazard, spot
      type(instrument_spec), intent(in) :: inst
      type(market_spec), intent(in) :: market
      integer, intent(in), optional :: min_steps
      real(dp) :: price
      type(market_spec) :: trial_market
      integer :: ns

      ns = 75
      if (present(min_steps)) ns = min_steps
      trial_market = market
      trial_market%use_hazard_link = .true.
      trial_market%hazard%base_intensity = max(0.0_dp, base_hazard)
      trial_market%hazard%constant_fraction = &
         min(1.0_dp, max(0.0_dp, constant_fraction))
      trial_market%hazard%power = max(0.0_dp, power)
      trial_market%hazard%reference_spot = spot
      price = find_present_value(spot, ns, inst, trial_market)
   end function price_with_intensity_link

   function penalty_with_intensity_link(power, constant_fraction, &
                                        base_hazard, spot, instruments, &
                                        prices, spreads, market, min_steps) &
                                        result(penalty)
      real(dp), intent(in) :: power, constant_fraction, base_hazard, spot
      type(instrument_spec), intent(in) :: instruments(:)
      real(dp), intent(in) :: prices(:)
      real(dp), intent(in), optional :: spreads(:)
      type(market_spec), intent(in) :: market
      integer, intent(in), optional :: min_steps
      real(dp) :: penalty, error_value, scale, model_price
      integer :: i

      penalty = 0.0_dp
      do i = 1, size(instruments)
         scale = 1.0_dp
         if (present(spreads)) then
            scale = max(spreads(i), sqrt(epsilon(1.0_dp)))
         end if
         model_price = price_with_intensity_link(power, constant_fraction, &
                                                  base_hazard, spot, &
                                                  instruments(i), market, &
                                                  min_steps)
         error_value = (model_price-prices(i))/scale
         penalty = penalty+error_value*error_value
      end do
   end function penalty_with_intensity_link

   subroutine fit_to_option_market(spot, instruments, prices, market, power, &
                                   constant_fraction, base_hazard, spreads, &
                                   min_steps, status)
      real(dp), intent(in) :: spot
      type(instrument_spec), intent(in) :: instruments(:)
      real(dp), intent(in) :: prices(:)
      type(market_spec), intent(in) :: market
      real(dp), intent(out) :: power, constant_fraction, base_hazard
      real(dp), intent(in), optional :: spreads(:)
      integer, intent(in), optional :: min_steps
      integer, intent(out), optional :: status
      real(dp) :: parameters(3), step(3), candidate(3)
      real(dp) :: best, trial
      integer :: iter, j

      parameters = [1.0_dp, 0.9_dp, &
                    max(0.01_dp, market%default_intensity)]
      step = [0.5_dp, 0.1_dp, 0.01_dp]
      best = penalty_with_intensity_link( &
         parameters(1), parameters(2), parameters(3), spot, instruments, &
         prices, spreads, market, min_steps)
      do iter = 1, 80
         do j = 1, 3
            candidate = parameters
            candidate(j) = candidate(j)+step(j)
            call constrain_intensity_parameters(candidate)
            trial = penalty_with_intensity_link( &
               candidate(1), candidate(2), candidate(3), spot, instruments, &
               prices, spreads, market, min_steps)
            if (trial < best) then
               parameters = candidate
               best = trial
            else
               candidate = parameters
               candidate(j) = candidate(j)-step(j)
               call constrain_intensity_parameters(candidate)
               trial = penalty_with_intensity_link( &
                  candidate(1), candidate(2), candidate(3), spot, &
                  instruments, prices, spreads, market, min_steps)
               if (trial < best) then
                  parameters = candidate
                  best = trial
               end if
            end if
         end do
         step = step*0.85_dp
      end do
      power = parameters(1)
      constant_fraction = parameters(2)
      base_hazard = parameters(3)
      if (present(status)) status = ragtop_ok
   end subroutine fit_to_option_market

   pure subroutine constrain_intensity_parameters(parameters)
      real(dp), intent(inout) :: parameters(3)
      parameters(1) = max(0.0_dp, parameters(1))
      parameters(2) = min(1.0_dp, max(0.0_dp, parameters(2)))
      parameters(3) = max(0.0_dp, parameters(3))
   end subroutine constrain_intensity_parameters

end module ragtop_calibration
