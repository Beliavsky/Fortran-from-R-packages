! SPDX-License-Identifier: GPL-3.0-or-later
module quant_bond_curves
   use qbc_kinds, only : dp
   use qbc_status, only : qbc_success, qbc_invalid_argument, qbc_size_mismatch, &
      qbc_no_convergence, qbc_singular, qbc_infeasible, qbc_status_message
   use qbc_dates, only : qbc_date, make_date, parse_date, date_string, valid_date, &
      days_between, add_days, add_months, add_years, is_leap_year, days_in_month, &
      day_of_week, is_business_day, adjust_business_day, discount_time, year_fraction, &
      operator(==), operator(/=), operator(<), operator(<=), operator(>), operator(>=)
   use qbc_types
   use qbc_cashflows, only : coupon_dates, coupon_dates_from_years, coupon_cashflows, coupons => coupon_cashflows, accrued_interest
   use qbc_curves, only : curve_rate, interpolate_curve, spot_to_forward, spot2forward => spot_to_forward, &
      forward_to_spot, fwd2spot => forward_to_spot, &
      discount_factor, discount_factors, curve_from_market_yields
   use qbc_bonds, only : valuation_bonds, dirty_to_clean, price_dirty2clean, bond_price_to_rate, bond_price2rate, &
      bond_sensitivity, sens_bonds, average_life, accrued_interests
   use qbc_swaps, only : valuation_swaps, swap_leg_value
   use qbc_optimization, only : qbc_optimizer_result, bounded_nelder_mead
   use qbc_calibration, only : curve_calibration, curve_calculation, bootstrap_curve, basis_curve
   implicit none
   public
end module quant_bond_curves
