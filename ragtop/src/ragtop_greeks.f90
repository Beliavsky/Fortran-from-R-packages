! SPDX-License-Identifier: GPL-2.0-or-later
module ragtop_greeks
   use ragtop_kinds, only : dp
   use ragtop_constants
   use ragtop_types, only : market_spec, instrument_spec, greek_result, price_grid
   use ragtop_pde, only : find_present_value, form_present_value_grid, american
   use ragtop_math, only : cubic_spline_interp
   implicit none
   private
   public :: construct_descending_bumps, resolve_bump
   public :: greek_by_fd, robust_greek, grid_delta_gamma, find_greeks
   public :: bump_rate, bump_volatility, bump_hazard

contains

   function construct_descending_bumps(level, divisor, n) result(bumps)
      real(dp), intent(in) :: level
      real(dp), intent(in), optional :: divisor
      integer, intent(in), optional :: n
      real(dp), allocatable :: bumps(:)
      real(dp) :: div
      integer :: count, i
      div = 40.0_dp
      if (present(divisor)) div = divisor
      count = 3
      if (present(n)) count = n
      allocate(bumps(count))
      do i = 1, count
         bumps(i) = max(abs(level),1.0_dp)/div**i
      end do
   end function construct_descending_bumps

   pure real(dp) function resolve_bump(user_bump, level, default_divisor) &
                                       result(bump)
      real(dp), intent(in), optional :: user_bump
      real(dp), intent(in) :: level
      real(dp), intent(in), optional :: default_divisor
      real(dp) :: divisor
      divisor = 1000.0_dp
      if (present(default_divisor)) divisor = default_divisor
      bump = max(abs(level),1.0_dp)/divisor
      if (present(user_bump)) then
         if (user_bump > 0.0_dp) bump = user_bump
      end if
   end function resolve_bump

   pure real(dp) function greek_by_fd(v0, v_up, v_down, bump) result(greek)
      real(dp), intent(in) :: v0, v_up, v_down, bump
      greek = (v_up-v_down)/(2.0_dp*bump)+0.0_dp*v0
   end function greek_by_fd

   function robust_greek(v0, values_up, values_down, bumps) result(greek)
      real(dp), intent(in) :: v0, values_up(:), values_down(:), bumps(:)
      real(dp) :: greek
      real(dp), allocatable :: estimates(:)
      integer :: n
      n = min(size(values_up),min(size(values_down),size(bumps)))
      allocate(estimates(n))
      estimates = (values_up(1:n)-values_down(1:n))/(2.0_dp*bumps(1:n))
      if (n == 1) then
         greek = estimates(1)
      else
         greek = estimates(n)+(estimates(n)-estimates(n-1))/3.0_dp
      end if
      greek = greek+0.0_dp*v0
   end function robust_greek

   pure function bump_rate(market, bump) result(bumped)
      type(market_spec), intent(in) :: market
      real(dp), intent(in) :: bump
      type(market_spec) :: bumped
      bumped = market
      if (bumped%use_rate_curve .and. allocated(bumped%rates%rate)) then
         bumped%rates%rate = bumped%rates%rate+bump
         bumped%rates%df = exp(-bumped%rates%time*bumped%rates%rate)
         call refresh_forward_rates(bumped)
      else
         bumped%short_rate = bumped%short_rate+bump
      end if
   end function bump_rate

   pure subroutine refresh_forward_rates(market)
      type(market_spec), intent(inout) :: market
      integer :: i, n
      if (.not. allocated(market%rates%time)) return
      n = size(market%rates%time)
      if (n == 1) then
         market%rates%forward_rate(1) = market%rates%rate(1)
      else
         do i = 1, n-1
            market%rates%forward_rate(i) = &
               -log(market%rates%df(i+1)/market%rates%df(i))/ &
               (market%rates%time(i+1)-market%rates%time(i))
         end do
         market%rates%forward_rate(n) = market%rates%forward_rate(n-1)
      end if
   end subroutine refresh_forward_rates

   pure function bump_volatility(market, bump) result(bumped)
      type(market_spec), intent(in) :: market
      real(dp), intent(in) :: bump
      type(market_spec) :: bumped
      real(dp), allocatable :: cumulative(:)
      bumped = market
      if (bumped%use_vol_curve .and. allocated(bumped%vols%volatility)) then
         bumped%vols%volatility = max(0.0_dp,bumped%vols%volatility+bump)
         cumulative = bumped%vols%time*bumped%vols%volatility**2
         bumped%vols%cumulative_variance = cumulative
         call refresh_forward_volatilities(bumped)
      else
         bumped%volatility = max(0.0_dp,bumped%volatility+bump)
      end if
   end function bump_volatility

   pure subroutine refresh_forward_volatilities(market)
      type(market_spec), intent(inout) :: market
      real(dp) :: previous_time, previous_variance, dv, dt
      integer :: i
      previous_time = 0.0_dp
      previous_variance = 0.0_dp
      do i = 1, size(market%vols%time)
         dt = market%vols%time(i)-previous_time
         dv = market%vols%cumulative_variance(i)-previous_variance
         market%vols%forward_volatility(i) = sqrt(max(0.0_dp,dv/dt))
         previous_time = market%vols%time(i)
         previous_variance = market%vols%cumulative_variance(i)
      end do
   end subroutine refresh_forward_volatilities

   pure function bump_hazard(market, bump) result(bumped)
      type(market_spec), intent(in) :: market
      real(dp), intent(in) :: bump
      type(market_spec) :: bumped
      bumped = market
      if (bumped%use_hazard_link) then
         bumped%hazard%base_intensity = &
            max(0.0_dp,bumped%hazard%base_intensity+bump)
      else
         bumped%default_intensity = max(0.0_dp,bumped%default_intensity+bump)
      end if
   end function bump_hazard

   subroutine grid_delta_gamma(grid, spot, delta, gamma, status)
      type(price_grid), intent(in) :: grid
      real(dp), intent(in) :: spot
      real(dp), intent(out) :: delta, gamma
      integer, intent(out), optional :: status
      real(dp) :: h, v0, vu, vd
      if (.not. allocated(grid%stock) .or. size(grid%stock) < 3) then
         delta = 0.0_dp
         gamma = 0.0_dp
         if (present(status)) status = ragtop_invalid_argument
         return
      end if
      h = max(1.0e-4_dp*spot,1.0e-6_dp)
      v0 = cubic_spline_interp(grid%stock,grid%value,spot)
      vu = cubic_spline_interp(grid%stock,grid%value,spot+h)
      vd = cubic_spline_interp(grid%stock,grid%value,spot-h)
      delta = (vu-vd)/(2.0_dp*h)
      gamma = (vu-2.0_dp*v0+vd)/(h*h)
      if (present(status)) status = ragtop_ok
   end subroutine grid_delta_gamma

   function price_for_greeks(spot, min_steps, inst, market) result(price)
      real(dp), intent(in) :: spot
      integer, intent(in) :: min_steps
      type(instrument_spec), intent(in) :: inst
      type(market_spec), intent(in) :: market
      real(dp) :: price
      if (inst%kind == instrument_american_option) then
         price = american(inst%callput,spot,inst%strike,inst%maturity, &
                          market,min_steps)
      else
         price = find_present_value(spot,min_steps,inst,market)
      end if
   end function price_for_greeks

   function find_greeks(spot, min_steps, inst, market, spot_bump, vol_bump, &
                        rate_bump, hazard_bump, time_bump) result(greeks)
      real(dp), intent(in) :: spot
      integer, intent(in) :: min_steps
      type(instrument_spec), intent(in) :: inst
      type(market_spec), intent(in) :: market
      real(dp), intent(in), optional :: spot_bump, vol_bump, rate_bump
      real(dp), intent(in), optional :: hazard_bump, time_bump
      type(greek_result) :: greeks
      type(market_spec) :: up_market, down_market
      type(instrument_spec) :: shorter
      real(dp) :: ds, dv, dr, dh, dt
      real(dp) :: up, down
      type(price_grid) :: grid

      ds = resolve_bump(spot_bump,spot,1000.0_dp)
      dv = resolve_bump(vol_bump,market%volatility,100.0_dp)
      dr = resolve_bump(rate_bump,market%short_rate,10000.0_dp)
      dh = resolve_bump(hazard_bump,market%default_intensity,10000.0_dp)
      dt = resolve_bump(time_bump,inst%maturity,365.0_dp)

      greeks%price = price_for_greeks(spot,min_steps,inst,market)
      call form_present_value_grid(spot,min_steps,inst,market,grid)
      call grid_delta_gamma(grid,spot,greeks%delta,greeks%gamma,greeks%status)

      up_market = bump_volatility(market,dv)
      down_market = bump_volatility(market,-dv)
      up = price_for_greeks(spot,min_steps,inst,up_market)
      down = price_for_greeks(spot,min_steps,inst,down_market)
      greeks%vega = (up-down)/(2.0_dp*dv)

      up_market = bump_rate(market,dr)
      down_market = bump_rate(market,-dr)
      up = price_for_greeks(spot,min_steps,inst,up_market)
      down = price_for_greeks(spot,min_steps,inst,down_market)
      greeks%rho = (up-down)/(2.0_dp*dr)

      up_market = bump_hazard(market,dh)
      down_market = bump_hazard(market,-dh)
      up = price_for_greeks(spot,min_steps,inst,up_market)
      down = price_for_greeks(spot,min_steps,inst,down_market)
      greeks%hazard_sensitivity = (up-down)/(2.0_dp*dh)

      shorter = inst
      shorter%maturity = max(dt,inst%maturity-dt)
      up = price_for_greeks(spot,min_steps,shorter,market)
      greeks%theta = (up-greeks%price)/dt
      greeks%status = ragtop_ok
   end function find_greeks

end module ragtop_greeks
