! SPDX-License-Identifier: GPL-2.0-or-later
module ragtop_pde
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
   use ragtop_kinds, only : dp
   use ragtop_constants
   use ragtop_types, only : instrument_spec, market_spec, dividend_schedule, grid_spec, price_grid, option_value
   use ragtop_math, only : solve_tridiagonal, cubic_spline_interp, sort_unique_real
   use ragtop_term_structures, only : discount_factor, cumulative_variance, default_intensity, &
                                      effective_rate, effective_volatility
   use ragtop_cashflows, only : adjust_for_dividends
   use ragtop_instruments, only : terminal_values, recovery_values, apply_optionality, &
                                  instrument_cashflow_between, critical_times, EuropeanOption, AmericanOption
   use ragtop_black_scholes, only : black_scholes_on_term_structures
   implicit none
   private
   public :: construct_implicit_grid_structure, construct_tridiagonals
   public :: infer_conforming_time_grid, pde_matrix_solve
   public :: form_present_value_grid, find_present_value, american

contains

   subroutine construct_implicit_grid_structure(tenor, min_steps, spot, center, drift_rate, volatility, &
                                                 structure_constant, std_devs_width, grid, status)
      real(dp), intent(in) :: tenor, spot, center, drift_rate, volatility
      integer, intent(in) :: min_steps
      real(dp), intent(in), optional :: structure_constant, std_devs_width
      type(grid_spec), intent(out) :: grid
      integer, intent(out), optional :: status
      real(dp) :: sc, width_mult, dt, dz, width
      integer :: half_n, i, st
      st = ragtop_ok
      sc = 2.0_dp
      width_mult = 3.0_dp
      if (present(structure_constant)) sc = structure_constant
      if (present(std_devs_width)) width_mult = std_devs_width
      if (tenor <= 0.0_dp .or. min_steps < 1 .or. spot <= 0.0_dp .or. center <= 0.0_dp .or. &
          volatility < 0.0_dp .or. sc <= 0.0_dp .or. width_mult <= 0.0_dp) then
         st = ragtop_invalid_argument
         if (present(status)) status = st
         return
      end if
      dt = tenor/real(min_steps,dp)
      dz = sqrt(dt/sc)
      width = max(width_mult*volatility*sqrt(tenor),4.0_dp*dz)
      half_n = max(2,ceiling(width/dz))
      grid%n_time_steps = min_steps
      grid%n_space = 2*half_n+1
      grid%t_max = tenor
      grid%dt = dt
      grid%dz = dz
      grid%z0 = log(spot/center)+(drift_rate-0.5_dp*volatility**2)*tenor
      grid%z_width = width
      allocate(grid%z(grid%n_space))
      do i = 1, grid%n_space
         grid%z(i) = grid%z0+dz*real(i-half_n-1,dp)
      end do
      if (present(status)) status = st
   end subroutine construct_implicit_grid_structure

   subroutine construct_tridiagonals(volatility, structure_constant, drift, sub, diag, super, status)
      real(dp), intent(in) :: volatility, structure_constant, drift(:)
      real(dp), intent(out) :: sub(:), diag(:), super(:)
      integer, intent(out), optional :: status
      real(dp) :: sig2sc
      integer :: n, i, st
      n = size(drift)
      st = ragtop_ok
      if (n < 2 .or. size(diag) /= n .or. size(sub) /= n-1 .or. size(super) /= n-1) then
         st = ragtop_invalid_argument
         if (present(status)) status = st
         return
      end if
      sig2sc = volatility*volatility*structure_constant
      diag = 1.0_dp+sig2sc
      do i = 1, n-1
         sub(i) = -0.5_dp*(sig2sc-drift(i+1))
         super(i) = -0.5_dp*(sig2sc+drift(i))
      end do
      diag(1) = 1.0_dp+drift(1)
      super(1) = -drift(1)
      diag(n) = 1.0_dp+drift(n)
      sub(n-1) = -drift(n)
      if (any(.not.ieee_is_finite(diag)) .or. any(.not.ieee_is_finite(sub)) .or. &
          any(.not.ieee_is_finite(super))) then
         st = ragtop_invalid_argument
      end if
      if (present(status)) status = st
   end subroutine construct_tridiagonals

   subroutine pde_matrix_solve(sub, diag, super, rhs, x, status)
      real(dp), intent(in) :: sub(:), diag(:), super(:), rhs(:)
      real(dp), intent(out) :: x(:)
      integer, intent(out), optional :: status
      call solve_tridiagonal(sub,diag,super,rhs,x,status)
   end subroutine pde_matrix_solve

   subroutine infer_conforming_time_grid(min_steps, tmax, inst, times)
      integer, intent(in) :: min_steps
      real(dp), intent(in) :: tmax
      type(instrument_spec), intent(in), optional :: inst
      real(dp), allocatable, intent(out) :: times(:)
      real(dp), allocatable :: all_times(:), extra(:), unique_times(:)
      integer :: i, nbase, n
      nbase = min_steps+1
      if (present(inst)) then
         call critical_times(inst,extra)
         allocate(all_times(nbase+size(extra)))
      else
         allocate(all_times(nbase))
      end if
      do i = 1, nbase
         all_times(i) = tmax*real(i-1,dp)/real(min_steps,dp)
      end do
      if (present(inst)) all_times(nbase+1:) = extra
      all_times = max(0.0_dp,min(tmax,all_times))
      call sort_unique_real(all_times,unique_times,1.0e-12_dp)
      n = count(unique_times >= 0.0_dp .and. unique_times <= tmax)
      allocate(times(n))
      times = pack(unique_times,unique_times >= 0.0_dp .and. unique_times <= tmax)
   end subroutine infer_conforming_time_grid

   pure subroutine stock_levels(z, t, tmax, center, drift_rate, volatility, stock)
      real(dp), intent(in) :: z(:), t, tmax, center, drift_rate, volatility
      real(dp), intent(out) :: stock(:)
      stock = center*exp(z-(drift_rate-0.5_dp*volatility**2)*(tmax-t))
   end subroutine stock_levels

   subroutine form_present_value_grid(spot, min_steps, inst, market, result, dividends, &
                                      structure_constant, std_devs_width, grid_center)
      real(dp), intent(in) :: spot
      integer, intent(in) :: min_steps
      type(instrument_spec), intent(in) :: inst
      type(market_spec), intent(in) :: market
      type(price_grid), intent(out) :: result
      type(dividend_schedule), intent(in), optional :: dividends
      real(dp), intent(in), optional :: structure_constant, std_devs_width, grid_center
      type(grid_spec) :: grid
      type(dividend_schedule) :: no_dividends
      real(dp), allocatable :: times(:), stock(:), stock_next(:), values(:), adjusted(:,:), solved(:)
      real(dp), allocatable :: recovery(:), optional_values(:), h(:), surv(:), drift(:)
      real(dp), allocatable :: sub(:), diag(:), super(:)
      real(dp) :: t0, t1, dt, full_df, prev_df, local_df, local_rate, local_vol
      real(dp) :: r_eff, vol_eff, c, center, sc, width, cf
      integer :: m, n, st
      sc = 2.0_dp
      width = 3.0_dp
      if (present(structure_constant)) sc = structure_constant
      if (present(std_devs_width)) width = std_devs_width
      center = spot
      if (inst%strike > 0.0_dp) center = inst%strike
      if (present(grid_center)) center = grid_center
      result%status = ragtop_ok
      if (spot <= 0.0_dp .or. min_steps < 1 .or. inst%maturity <= 0.0_dp) then
         result%status = ragtop_invalid_argument
         return
      end if
      r_eff = effective_rate(market,inst%maturity)
      vol_eff = effective_volatility(market,inst%maturity)
      c = r_eff-market%dividend_rate-market%borrow_cost
      call construct_implicit_grid_structure(inst%maturity,min_steps,spot,center,c,vol_eff,sc,width,grid,st)
      if (st /= ragtop_ok) then
         result%status = st
         return
      end if
      call infer_conforming_time_grid(min_steps,inst%maturity,inst,times)
      n = grid%n_space
      allocate(stock(n),stock_next(n),values(n),solved(n),recovery(n),optional_values(n),h(n),surv(n),drift(n))
      allocate(sub(n-1),diag(n),super(n-1),adjusted(n,1))
      call stock_levels(grid%z,inst%maturity,inst%maturity,center,c,vol_eff,stock_next)
      call terminal_values(inst,stock_next,market,values)
      values = discount_factor(market,inst%maturity,0.0_dp)*values
      do m = size(times)-1, 1, -1
         t0 = times(m)
         t1 = times(m+1)
         dt = t1-t0
         if (dt <= 0.0_dp) cycle
         full_df = discount_factor(market,t0,0.0_dp)
         prev_df = discount_factor(market,t1,0.0_dp)
         local_df = discount_factor(market,t1,t0)
         local_rate = -log(max(local_df,tiny(1.0_dp)))/dt
         local_vol = sqrt(max(0.0_dp,cumulative_variance(market,t1,t0)/dt))
         call stock_levels(grid%z,t0,inst%maturity,center,c,vol_eff,stock)
         h = default_intensity(market,t0,stock)
         adjusted(:,1) = values
         if (present(dividends)) then
            call adjust_for_dividends(reshape(values,[n,1]),t0,dt,local_rate,h,stock,spot,dividends,adjusted)
         else
            call adjust_for_dividends(reshape(values,[n,1]),t0,dt,local_rate,h,stock,spot,no_dividends,adjusted)
         end if
         values = adjusted(:,1)
         cf = instrument_cashflow_between(inst,t0,t1,market)
         if (abs(cf) > 0.0_dp) values = values+cf*prev_df
         surv = exp(-h*dt)
         drift = h*dt/grid%dz
         call construct_tridiagonals(local_vol,dt/(grid%dz*grid%dz),drift,sub,diag,super,st)
         if (st /= ragtop_ok) then
            result%status = st
            return
         end if
         call solve_tridiagonal(sub,diag,super,values,solved,st)
         if (st /= ragtop_ok) then
            result%status = st
            return
         end if
         solved = max(0.0_dp,solved)
         call recovery_values(inst,solved/max(full_df,sqrt(tiny(1.0_dp))),stock,t0,market,recovery)
         values = surv*solved+(1.0_dp-surv)*local_df*full_df*recovery
         call apply_optionality(inst,values/max(full_df,sqrt(tiny(1.0_dp))),stock,t0,market,optional_values)
         values = full_df*optional_values
         stock_next = stock
      end do
      allocate(result%stock(n),result%value(n))
      result%stock = stock_next
      result%value = values
   end subroutine form_present_value_grid

   function find_present_value(spot, min_steps, inst, market, dividends, structure_constant, &
                               std_devs_width, grid_center, status) result(value)
      real(dp), intent(in) :: spot
      integer, intent(in) :: min_steps
      type(instrument_spec), intent(in) :: inst
      type(market_spec), intent(in) :: market
      type(dividend_schedule), intent(in), optional :: dividends
      real(dp), intent(in), optional :: structure_constant, std_devs_width, grid_center
      integer, intent(out), optional :: status
      real(dp) :: value
      type(price_grid) :: grid
      if (present(dividends)) then
         call form_present_value_grid(spot,min_steps,inst,market,grid,dividends,structure_constant,std_devs_width,grid_center)
      else
         call form_present_value_grid(spot,min_steps,inst,market,grid,structure_constant=structure_constant, &
                                      std_devs_width=std_devs_width,grid_center=grid_center)
      end if
      if (grid%status == ragtop_ok) then
         value = cubic_spline_interp(grid%stock,grid%value,spot)
      else
         value = 0.0_dp
      end if
      if (present(status)) status = grid%status
   end function find_present_value

   function american(callput, spot, strike, maturity, market, min_steps, dividends, &
                     structure_constant, std_devs_width, status) result(value)
      integer, intent(in) :: callput
      real(dp), intent(in) :: spot, strike, maturity
      type(market_spec), intent(in) :: market
      integer, intent(in), optional :: min_steps
      type(dividend_schedule), intent(in), optional :: dividends
      real(dp), intent(in), optional :: structure_constant, std_devs_width
      integer, intent(out), optional :: status
      real(dp) :: value, pa, pe, exact, early
      type(option_value) :: bs
      integer :: nsteps, sta, ste
      type(instrument_spec) :: ai, ei
      nsteps = 100
      if (present(min_steps)) nsteps = min_steps
      ai = AmericanOption(maturity,strike,callput,'American')
      ei = EuropeanOption(maturity,strike,callput,'European')
      if (present(dividends)) then
         pa = find_present_value(spot,nsteps,ai,market,dividends,structure_constant,std_devs_width,status=sta)
         pe = find_present_value(spot,nsteps,ei,market,dividends,structure_constant,std_devs_width,status=ste)
         bs = black_scholes_on_term_structures(callput,spot,strike,maturity,market,dividends)
         exact = bs%price
      else
         pa = find_present_value(spot,nsteps,ai,market,structure_constant=structure_constant, &
                                 std_devs_width=std_devs_width,status=sta)
         pe = find_present_value(spot,nsteps,ei,market,structure_constant=structure_constant, &
                                 std_devs_width=std_devs_width,status=ste)
         bs = black_scholes_on_term_structures(callput,spot,strike,maturity,market)
         exact = bs%price
      end if
      value = pa+(exact-pe)
      early = max(0.0_dp,real(callput,dp)*(spot-strike))
      value = max(value,early)
      if (present(status)) then
         if (sta == ragtop_ok .and. ste == ragtop_ok) then
            status = ragtop_ok
         else
            status = max(sta,ste)
         end if
      end if
   end function american

end module ragtop_pde
