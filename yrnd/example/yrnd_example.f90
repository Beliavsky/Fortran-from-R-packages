program yrnd_example
   use yrnd
   implicit none

   real(dp) :: strikes(7), calls(7), puts(7)
   type(lognormal_mixture_t) :: source
   type(density_result_t) :: fitted
   type(transformed_density_t) :: rates
   type(date_t) :: today, option_date, futures_date
   integer :: i

   strikes = [(92.0_dp + 3.0_dp * real(i - 1, dp), i = 1, size(strikes))]
   source%n_components = 2
   source%meanlog(1) = log(98.0_dp) - 0.5_dp * 0.08_dp ** 2
   source%meanlog(2) = log(102.0_dp) - 0.5_dp * 0.15_dp ** 2
   source%sdlog(1:2) = [0.08_dp, 0.15_dp]
   source%weight(1:2) = [0.50_dp, 0.50_dp]

   call mixture_option_prices(source, strikes, strikes, 0.025_dp, 0.5_dp, &
      option_futures_margin, calls, puts)

   today = date_t(2026, 8, 5)
   option_date = date_t(2027, 2, 5)
   futures_date = date_t(2027, 3, 15)

   call bond_future_price(calls, strikes, puts, strikes, 2, 0.025_dp, dc_act_365, &
      option_futures_margin, source%mean(), futures_date, option_date, today, &
      fitted, grid_step=0.05_dp)
   call future_price_to_stir_rate(fitted, rates)

   write(*, '(a,f10.4)') "Future-price mean: ", fitted%moments(1)
   write(*, '(a,f10.4)') "Future-price standard deviation: ", fitted%moments(2)
   write(*, '(a,f10.4)') "Implied rate mean: ", rates%moments(1)
   write(*, '(a,f10.4)') "Median future price: ", fitted%quantiles(7)
end program yrnd_example
