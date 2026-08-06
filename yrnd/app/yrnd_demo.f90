program yrnd_demo
   use yrnd
   implicit none

   real(dp) :: strikes(5), calls(5), puts(5)
   type(lognormal_mixture_t) :: source
   type(density_result_t) :: result
   type(date_t) :: start_date, option_date, futures_date

   strikes = [94.0_dp, 97.0_dp, 100.0_dp, 103.0_dp, 106.0_dp]
   source%n_components = 2
   source%meanlog(1:2) = [log(98.0_dp), log(102.0_dp)]
   source%sdlog(1:2) = [0.08_dp, 0.15_dp]
   source%weight(1:2) = [0.5_dp, 0.5_dp]
   call mixture_option_prices(source, strikes, strikes, 0.025_dp, 0.5_dp, &
      option_european, calls, puts)

   start_date = date_t(2026, 8, 5)
   option_date = date_t(2027, 2, 5)
   futures_date = date_t(2027, 3, 15)
   call bond_future_price(calls, strikes, puts, strikes, 2, 0.025_dp, dc_act_365, &
      option_european, source%mean(), futures_date, option_date, start_date, &
      result, grid_step=0.05_dp)

   write(*, '(a,f10.5)') "Objective: ", result%objective
   write(*, '(a,f10.5)') "Mean: ", result%moments(1)
   write(*, '(a,f10.5)') "Standard deviation: ", result%moments(2)
end program yrnd_demo
