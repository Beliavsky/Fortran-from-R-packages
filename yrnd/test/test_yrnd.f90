program test_yrnd
   use yrnd
   implicit none

   integer, parameter :: nopt = 9
   real(dp) :: strikes(nopt), calls(nopt), puts(nopt)
   real(dp) :: area, y, dirty, objective3
   type(lognormal_mixture_t) :: known, known3, fit3
   type(density_result_t) :: fit_result
   type(transformed_density_t) :: rate_result, yield_result
   type(bond_t) :: bond, bonds(2)
   type(bond_context_t) :: context
   type(ctd_probability_result_t) :: probabilities
   type(spread_result_t) :: spread
   type(date_t) :: start_date, option_date, future_date
   integer :: i, status, failures, convergence3

   failures = 0
   strikes = [(90.0_dp + 2.5_dp * real(i - 1, dp), i = 1, nopt)]
   known%n_components = 2
   known%meanlog(1) = log(97.0_dp) - 0.5_dp * 0.07_dp ** 2
   known%meanlog(2) = log(103.0_dp) - 0.5_dp * 0.16_dp ** 2
   known%sdlog(1:2) = [0.07_dp, 0.16_dp]
   known%weight(1:2) = [0.45_dp, 0.55_dp]
   call mixture_option_prices(known, strikes, strikes, 0.03_dp, 0.5_dp, &
      option_european, calls, puts)

   start_date = date_t(2026, 8, 5)
   option_date = date_t(2027, 2, 5)
   future_date = date_t(2027, 3, 15)

   call bond_future_price(calls, strikes, puts, strikes, 2, 0.03_dp, dc_act_365, &
      option_european, known%mean(), future_date, option_date, start_date, &
      fit_result, grid_step=0.05_dp, max_iter=450)
   area = trapezoid_integral(fit_result%domain, fit_result%density)
   call check(abs(area - 1.0_dp) < 2.0e-5_dp, "future density integrates to one")
   call check(abs(fit_result%model%mean() - known%mean()) < 0.2_dp, "fitted forward mean")
   call check(fit_result%objective < 5.0e-3_dp, "mixture option fit objective")
   call check(all(fit_result%quantiles(2:) >= fit_result%quantiles(:12)), "ordered quantiles")

   known3%n_components = 3
   known3%meanlog(1) = log(94.0_dp) - 0.5_dp * 0.05_dp ** 2
   known3%meanlog(2) = log(100.0_dp) - 0.5_dp * 0.11_dp ** 2
   known3%meanlog(3) = log(108.0_dp) - 0.5_dp * 0.20_dp ** 2
   known3%sdlog(1:3) = [0.05_dp, 0.11_dp, 0.20_dp]
   known3%weight(1:3) = [0.20_dp, 0.50_dp, 0.30_dp]
   call mixture_option_prices(known3, strikes, strikes, 0.03_dp, 0.5_dp, &
      option_european, calls, puts)
   call fit_lognormal_mixture(calls, strikes, puts, strikes, 3, 0.03_dp, 0.5_dp, &
      option_european, known3%mean(), fit3, objective3, convergence3, max_iter=700)
   call check(objective3 < 1.0e-3_dp, "three-component mixture fit")
   call check(abs(fit3%mean() - known3%mean()) < 0.1_dp, "three-component forward mean")

   call mixture_option_prices(known, strikes, strikes, 0.03_dp, 0.5_dp, &
      option_european, calls, puts)

   call future_price_to_stir_rate(fit_result, rate_result)
   call check(abs(rate_result%moments(1) - (100.0_dp - fit_result%moments(1))) < 1.0e-4_dp, &
      "STIR rate transformation mean")
   call check(rate_result%domain(1) < rate_result%domain(size(rate_result%domain)), &
      "STIR domain ascending")

   bond%id = "BOND_A"
   bond%coupon = 0.04_dp
   bond%coupon_frequency = 2
   bond%maturity = date_t(2031, 6, 15)
   bond%conversion_factor = 0.91_dp
   bond%current_yield = 0.042_dp
   bond%nominal = 100.0_dp
   context = bond_context_t(start_date, option_date, future_date, dc_act_365, 1, 0.03_dp, 0.032_dp)

   dirty = dirty_price_from_yield(bond, bond%current_yield, option_date, 1, dc_act_365)
   call check(dirty > 50.0_dp .and. dirty < 150.0_dp, "bond dirty price")
   y = yield_from_future_price(bond, 110.0_dp, context, .true., status)
   call check(status == 0, "bond yield root convergence")
   call check(y > -0.5_dp .and. y < 1.0_dp, "bond yield plausible")

   call future_price_to_bond_yield(fit_result, bond, context, yield_result, .true.)
   area = trapezoid_integral(yield_result%domain, yield_result%density)
   call check(abs(area - 1.0_dp) < 2.0e-5_dp, "bond yield density integrates to one")
   call check(all(yield_result%domain(2:) >= yield_result%domain(:size(yield_result%domain)-1)), &
      "bond yield domain ascending")

   bonds(1) = bond
   bonds(2) = bond
   bonds(2)%id = "BOND_B"
   bonds(2)%coupon = 0.055_dp
   bonds(2)%maturity = date_t(2033, 12, 15)
   bonds(2)%conversion_factor = 0.96_dp
   bonds(2)%current_yield = 0.046_dp
   call ctd_probabilities(bonds, fit_result%domain, fit_result%density, context, probabilities, .true.)
   call check(abs(sum(probabilities%probability) - 1.0_dp) < 1.0e-10_dp, "CTD probabilities sum to one")
   call check(all(probabilities%probability >= 0.0_dp), "CTD probabilities nonnegative")

   call simulate_bond_yield_spread(fit_result%model, known%mean(), bonds, context, &
      fit_result%model, known%mean(), bonds, context, 0.5_dp, 300, spread, seed=12345)
   call check(size(spread%samples) == 300, "spread simulation size")
   call check(spread%stddev >= 0.0_dp, "spread standard deviation")
   area = trapezoid_integral(spread%domain, spread%density)
   call check(abs(area - 1.0_dp) < 1.0e-4_dp, "spread KDE integrates to one")

   if (failures /= 0) then
      write(*, '(a,i0)') "FAILED TESTS: ", failures
      error stop 1
   end if
   write(*, '(a)') "All yrnd tests passed."

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

end program test_yrnd
