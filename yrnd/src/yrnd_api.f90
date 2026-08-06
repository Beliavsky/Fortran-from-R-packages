! SPDX-License-Identifier: GPL-3.0-only
module yrnd_api
   use yrnd_kinds, only : dp
   use yrnd_dates, only : date_t, year_fraction, operator(<), operator(<=), operator(>)
   use yrnd_mixture, only : density_result_t, lognormal_mixture_t, fit_lognormal_mixture, build_density_result
   use yrnd_bonds, only : bond_t, bond_context_t, ctd_probability_result_t, ctd_probabilities
   use yrnd_transforms, only : transformed_density_t, spread_result_t, &
      future_price_to_stir_rate, future_price_to_bond_yield, simulate_bond_yield_spread
   implicit none
   private

   public :: bond_future_price, stir_future_price, stir_rate
   public :: ctd_bond_yield, proba_ctd, proba_ctd_opt, bond_yield_spread

contains

   subroutine bond_future_price(call_prices, call_strikes, put_prices, put_strikes, &
                                n_components, risk_free_rate, day_count, option_style, &
                                futures_price, futures_maturity, option_maturity, start_date, &
                                result, grid_step, max_iter)
      real(dp), intent(in) :: call_prices(:), call_strikes(:), put_prices(:), put_strikes(:)
      integer, intent(in) :: n_components, day_count, option_style
      real(dp), intent(in) :: risk_free_rate, futures_price
      type(date_t), intent(in) :: futures_maturity, option_maturity, start_date
      type(density_result_t), intent(out) :: result
      real(dp), intent(in), optional :: grid_step
      integer, intent(in), optional :: max_iter
      real(dp) :: term, objective
      integer :: convergence
      type(lognormal_mixture_t) :: model

      if (.not. (start_date < option_maturity) .or. .not. (option_maturity <= futures_maturity)) then
         error stop "bond_future_price: inconsistent dates"
      end if
      term = year_fraction(start_date, option_maturity, day_count)
      call fit_lognormal_mixture(call_prices, call_strikes, put_prices, put_strikes, &
         n_components, risk_free_rate, term, option_style, futures_price, &
         model, objective, convergence, max_iter)
      call build_density_result(model, call_strikes, put_strikes, risk_free_rate, term, &
         option_style, result, grid_step, objective, convergence)
   end subroutine bond_future_price

   subroutine stir_future_price(call_prices, call_strikes, put_prices, put_strikes, &
                                n_components, risk_free_rate, day_count, option_style, &
                                futures_price, futures_maturity, option_maturity, start_date, &
                                result, grid_step, max_iter)
      real(dp), intent(in) :: call_prices(:), call_strikes(:), put_prices(:), put_strikes(:)
      integer, intent(in) :: n_components, day_count, option_style
      real(dp), intent(in) :: risk_free_rate, futures_price
      type(date_t), intent(in) :: futures_maturity, option_maturity, start_date
      type(density_result_t), intent(out) :: result
      real(dp), intent(in), optional :: grid_step
      integer, intent(in), optional :: max_iter
      call bond_future_price(call_prices, call_strikes, put_prices, put_strikes, n_components, &
         risk_free_rate, day_count, option_style, futures_price, futures_maturity, &
         option_maturity, start_date, result, grid_step, max_iter)
   end subroutine stir_future_price

   subroutine stir_rate(call_prices, call_strikes, put_prices, put_strikes, &
                        n_components, risk_free_rate, day_count, option_style, &
                        futures_price, futures_maturity, option_maturity, start_date, &
                        result, grid_step, max_iter)
      real(dp), intent(in) :: call_prices(:), call_strikes(:), put_prices(:), put_strikes(:)
      integer, intent(in) :: n_components, day_count, option_style
      real(dp), intent(in) :: risk_free_rate, futures_price
      type(date_t), intent(in) :: futures_maturity, option_maturity, start_date
      type(transformed_density_t), intent(out) :: result
      real(dp), intent(in), optional :: grid_step
      integer, intent(in), optional :: max_iter
      type(density_result_t) :: future_result
      call stir_future_price(call_prices, call_strikes, put_prices, put_strikes, n_components, &
         risk_free_rate, day_count, option_style, futures_price, futures_maturity, &
         option_maturity, start_date, future_result, grid_step, max_iter)
      call future_price_to_stir_rate(future_result, result)
   end subroutine stir_rate

   subroutine ctd_bond_yield(call_prices, call_strikes, put_prices, put_strikes, &
                             n_components, option_rate, futures_rate, day_count, option_style, &
                             futures_price, futures_maturity, option_maturity, start_date, &
                             ctd_bond, settlement_days, result, grid_step, max_iter)
      real(dp), intent(in) :: call_prices(:), call_strikes(:), put_prices(:), put_strikes(:)
      integer, intent(in) :: n_components, day_count, option_style, settlement_days
      real(dp), intent(in) :: option_rate, futures_rate, futures_price
      type(date_t), intent(in) :: futures_maturity, option_maturity, start_date
      type(bond_t), intent(in) :: ctd_bond
      type(transformed_density_t), intent(out) :: result
      real(dp), intent(in), optional :: grid_step
      integer, intent(in), optional :: max_iter
      type(density_result_t) :: future_result
      type(bond_context_t) :: context
      call bond_future_price(call_prices, call_strikes, put_prices, put_strikes, n_components, &
         option_rate, day_count, option_style, futures_price, futures_maturity, &
         option_maturity, start_date, future_result, grid_step, max_iter)
      context%calibration_date = start_date
      context%option_date = option_maturity
      context%futures_date = futures_maturity
      context%day_count = day_count
      context%settlement_days = settlement_days
      context%option_zero_rate = option_rate
      context%futures_zero_rate = futures_rate
      call future_price_to_bond_yield(future_result, ctd_bond, context, result, .true.)
   end subroutine ctd_bond_yield

   subroutine proba_ctd(call_prices, call_strikes, put_prices, put_strikes, &
                        n_components, risk_free_rate, day_count, option_style, &
                        futures_price, futures_maturity, option_maturity, start_date, &
                        bonds, settlement_days, result, grid_step, max_iter)
      real(dp), intent(in) :: call_prices(:), call_strikes(:), put_prices(:), put_strikes(:)
      integer, intent(in) :: n_components, day_count, option_style, settlement_days
      real(dp), intent(in) :: risk_free_rate, futures_price
      type(date_t), intent(in) :: futures_maturity, option_maturity, start_date
      type(bond_t), intent(in) :: bonds(:)
      type(ctd_probability_result_t), intent(out) :: result
      real(dp), intent(in), optional :: grid_step
      integer, intent(in), optional :: max_iter
      type(density_result_t) :: future_result
      type(bond_context_t) :: context
      call bond_future_price(call_prices, call_strikes, put_prices, put_strikes, n_components, &
         risk_free_rate, day_count, option_style, futures_price, futures_maturity, &
         option_maturity, start_date, future_result, grid_step, max_iter)
      context%calibration_date = start_date
      context%option_date = option_maturity
      context%futures_date = option_maturity
      context%day_count = day_count
      context%settlement_days = settlement_days
      context%option_zero_rate = risk_free_rate
      context%futures_zero_rate = risk_free_rate
      call ctd_probabilities(bonds, future_result%domain, future_result%density, context, result, .false.)
   end subroutine proba_ctd

   subroutine proba_ctd_opt(call_prices, call_strikes, put_prices, put_strikes, &
                            n_components, option_rate, futures_rate, day_count, option_style, &
                            futures_price, futures_maturity, option_maturity, start_date, &
                            bonds, settlement_days, result, grid_step, max_iter)
      real(dp), intent(in) :: call_prices(:), call_strikes(:), put_prices(:), put_strikes(:)
      integer, intent(in) :: n_components, day_count, option_style, settlement_days
      real(dp), intent(in) :: option_rate, futures_rate, futures_price
      type(date_t), intent(in) :: futures_maturity, option_maturity, start_date
      type(bond_t), intent(in) :: bonds(:)
      type(ctd_probability_result_t), intent(out) :: result
      real(dp), intent(in), optional :: grid_step
      integer, intent(in), optional :: max_iter
      type(density_result_t) :: future_result
      type(bond_context_t) :: context
      call bond_future_price(call_prices, call_strikes, put_prices, put_strikes, n_components, &
         option_rate, day_count, option_style, futures_price, futures_maturity, &
         option_maturity, start_date, future_result, grid_step, max_iter)
      context%calibration_date = start_date
      context%option_date = option_maturity
      context%futures_date = futures_maturity
      context%day_count = day_count
      context%settlement_days = settlement_days
      context%option_zero_rate = option_rate
      context%futures_zero_rate = futures_rate
      call ctd_probabilities(bonds, future_result%domain, future_result%density, context, result, .true.)
   end subroutine proba_ctd_opt

   subroutine bond_yield_spread(call_prices1, call_strikes1, put_prices1, put_strikes1, &
                                call_prices2, call_strikes2, put_prices2, put_strikes2, &
                                option_rate, futures_rate, day_count1, day_count2, &
                                option_style1, option_style2, futures_price1, futures_price2, &
                                futures_maturity, option_maturity, start_date, &
                                bonds1, settlement_days1, bonds2, settlement_days2, &
                                correlation, n_simulations, result, seed, max_iter)
      real(dp), intent(in) :: call_prices1(:), call_strikes1(:), put_prices1(:), put_strikes1(:)
      real(dp), intent(in) :: call_prices2(:), call_strikes2(:), put_prices2(:), put_strikes2(:)
      real(dp), intent(in) :: option_rate, futures_rate, futures_price1, futures_price2, correlation
      integer, intent(in) :: day_count1, day_count2, option_style1, option_style2
      integer, intent(in) :: settlement_days1, settlement_days2, n_simulations
      type(date_t), intent(in) :: futures_maturity, option_maturity, start_date
      type(bond_t), intent(in) :: bonds1(:), bonds2(:)
      type(spread_result_t), intent(out) :: result
      integer, intent(in), optional :: seed, max_iter
      type(density_result_t) :: density1, density2
      type(bond_context_t) :: context1, context2

      call bond_future_price(call_prices1, call_strikes1, put_prices1, put_strikes1, 2, &
         option_rate, day_count1, option_style1, futures_price1, futures_maturity, &
         option_maturity, start_date, density1, grid_step=0.01_dp, max_iter=max_iter)
      call bond_future_price(call_prices2, call_strikes2, put_prices2, put_strikes2, 2, &
         option_rate, day_count2, option_style2, futures_price2, futures_maturity, &
         option_maturity, start_date, density2, grid_step=0.01_dp, max_iter=max_iter)
      context1 = bond_context_t(start_date, option_maturity, futures_maturity, day_count1, &
         settlement_days1, option_rate, futures_rate)
      context2 = bond_context_t(start_date, option_maturity, futures_maturity, day_count2, &
         settlement_days2, option_rate, futures_rate)
      call simulate_bond_yield_spread(density1%model, futures_price1, bonds1, context1, &
         density2%model, futures_price2, bonds2, context2, correlation, n_simulations, result, seed)
   end subroutine bond_yield_spread

end module yrnd_api
