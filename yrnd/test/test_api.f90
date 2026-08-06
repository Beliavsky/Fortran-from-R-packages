program test_api
   use yrnd
   implicit none

   integer, parameter :: nopt = 7
   real(dp) :: strikes(nopt), calls(nopt), puts(nopt)
   type(lognormal_mixture_t) :: known
   type(density_result_t) :: future_result
   type(transformed_density_t) :: rate_result, yield_result
   type(ctd_probability_result_t) :: p_now, p_future
   type(spread_result_t) :: spread
   type(bond_t) :: bonds(2)
   type(date_t) :: start_date, option_date, future_date
   integer :: i, failures

   failures = 0
   strikes = [(92.0_dp + 3.0_dp * real(i - 1, dp), i = 1, nopt)]
   known%n_components = 2
   known%meanlog(1) = log(98.0_dp) - 0.5_dp * 0.08_dp ** 2
   known%meanlog(2) = log(102.0_dp) - 0.5_dp * 0.15_dp ** 2
   known%sdlog(1:2) = [0.08_dp, 0.15_dp]
   known%weight(1:2) = [0.50_dp, 0.50_dp]
   call mixture_option_prices(known, strikes, strikes, 0.025_dp, 0.5_dp, &
      option_futures_margin, calls, puts)

   start_date = date_t(2026, 8, 5)
   option_date = date_t(2027, 2, 5)
   future_date = date_t(2027, 3, 15)

   bonds(1)%id = "A"
   bonds(1)%coupon = 0.04_dp
   bonds(1)%coupon_frequency = 2
   bonds(1)%maturity = date_t(2031, 6, 15)
   bonds(1)%conversion_factor = 0.92_dp
   bonds(1)%current_yield = 0.041_dp
   bonds(1)%nominal = 100.0_dp
   bonds(2) = bonds(1)
   bonds(2)%id = "B"
   bonds(2)%coupon = 0.0525_dp
   bonds(2)%maturity = date_t(2034, 12, 15)
   bonds(2)%conversion_factor = 0.97_dp
   bonds(2)%current_yield = 0.046_dp

   call stir_future_price(calls, strikes, puts, strikes, 2, 0.025_dp, dc_act_365, &
      option_futures_margin, known%mean(), future_date, option_date, start_date, &
      future_result, grid_step=0.10_dp, max_iter=300)
   call check(future_result%objective < 1.0e-2_dp, "stir_future_price")

   call stir_rate(calls, strikes, puts, strikes, 2, 0.025_dp, dc_act_365, &
      option_futures_margin, known%mean(), future_date, option_date, start_date, &
      rate_result, grid_step=0.10_dp, max_iter=300)
   call check(allocated(rate_result%density), "stir_rate")

   call ctd_bond_yield(calls, strikes, puts, strikes, 2, 0.025_dp, 0.027_dp, &
      dc_act_365, option_futures_margin, known%mean(), future_date, option_date, &
      start_date, bonds(1), 1, yield_result, grid_step=0.10_dp, max_iter=300)
   call check(allocated(yield_result%density), "ctd_bond_yield")

   call proba_ctd(calls, strikes, puts, strikes, 2, 0.025_dp, dc_act_365, &
      option_futures_margin, known%mean(), future_date, option_date, start_date, &
      bonds, 1, p_now, grid_step=0.20_dp, max_iter=300)
   call check(abs(sum(p_now%probability) - 1.0_dp) < 1.0e-10_dp, "proba_ctd")

   call proba_ctd_opt(calls, strikes, puts, strikes, 2, 0.025_dp, 0.027_dp, &
      dc_act_365, option_futures_margin, known%mean(), future_date, option_date, &
      start_date, bonds, 1, p_future, grid_step=0.20_dp, max_iter=300)
   call check(abs(sum(p_future%probability) - 1.0_dp) < 1.0e-10_dp, "proba_ctd_opt")

   call bond_yield_spread(calls, strikes, puts, strikes, calls, strikes, puts, strikes, &
      0.025_dp, 0.027_dp, dc_act_365, dc_act_365, option_futures_margin, &
      option_futures_margin, known%mean(), known%mean(), future_date, option_date, &
      start_date, bonds, 1, bonds, 1, 0.4_dp, 100, spread, seed=2026, max_iter=250)
   call check(size(spread%samples) == 100, "bond_yield_spread")

   if (failures > 0) error stop 1
   write(*, '(a)') "All yrnd API tests passed."

contains

   subroutine check(condition, name)
      logical, intent(in) :: condition
      character(len=*), intent(in) :: name
      if (condition) then
         write(*, '(a)') "PASS: " // trim(name)
      else
         write(*, '(a)') "FAIL: " // trim(name)
         failures = failures + 1
      end if
   end subroutine check

end program test_api
