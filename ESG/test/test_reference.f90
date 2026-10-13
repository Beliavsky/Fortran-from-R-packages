! SPDX-License-Identifier: GPL-2.0-or-later
program test_reference
   use esg
   use esg_reference_values
   implicit none
   type(esg_parameters) :: p
   type(esg_scenario) :: scenario
   real(dp) :: zero(85), innovations(3,5,5), flat(75), deviation
   real(dp), allocatable :: value(:), forward(:), curve(:), paths(:,:), means(:), joined(:)
   integer :: i, status
   p = esg_parameters(rate_vol = 0.017_dp, k = 0.42_dp, stock_vol = 0.23_dp, estate_vol = 0.12_dp, &
      stock0 = 1.1_dp, estate0 = 1.4_dp, rho = -0.35_dp, eta = 0.005_dp, liquidity0 = 0.011_dp, &
      default0 = 0.035_dp, default_vol = 0.018_dp, alpha = 0.6_dp, beta = 0.05_dp)
   do i = 1, 85
      zero(i) = 0.018_dp + 0.00003_dp * real(i, dp) + sin(real(i, dp) * 0.2_dp) * 0.0002_dp
   end do
   do i = 1, 75
      flat(i) = sin(real(i, dp) * 0.7_dp) * 0.7_dp
   end do
   innovations = reshape(flat, shape(innovations))
   call esg_forward_rates(zero, 5, forward, status)
   call check(status == 0, 'setForwardRates status')
   call near(forward, ref_forward, 'setForwardRates')
   call esg_zero_rates(zero, 5, curve, status)
   call near(curve, ref_zero, 'setZCRates helper')
   call esg_generate(p, zero, innovations, scenario, status)
   call check(status == 0, 'customPathsGeneration status')
   call near(reshape(scenario%rate,[18]), ref_rate, 'rate paths')
   call near(reshape(scenario%stock,[18]), ref_stock, 'stock paths')
   call near(reshape(scenario%estate,[18]), ref_estate, 'estate paths')
   call near(reshape(scenario%liquidity,[18]), ref_liquidity, 'liquidity paths')
   call near(reshape(scenario%default,[18]), ref_default, 'default paths')
   joined = [reshape(scenario%rate,[18]), reshape(scenario%stock,[18]), reshape(scenario%estate,[18]), &
      reshape(scenario%liquidity,[18]), reshape(scenario%default,[18])]
   call near(joined, ref_rAllRisksFactors, 'rAllRisksFactors')
   call esg_short_rate(p, forward, innovations(:,:,1), paths, status)
   call near(reshape(paths,[18]), ref_rShortRate, 'rShortRate')
   call esg_stock(p, scenario%rate, innovations(:,:,1), innovations(:,:,2), paths, status)
   call near(reshape(paths,[18]), ref_rStock, 'rStock')
   call esg_real_estate(p, scenario%rate, innovations(:,:,3), paths, status)
   call near(reshape(paths,[18]), ref_rRealEstate, 'rRealEstate')
   call esg_liquidity(p, innovations(:,:,4), paths, status)
   call near(reshape(paths,[18]), ref_rLiquiditySpread, 'rLiquiditySpread')
   call esg_default(p, innovations(:,:,5), paths, status)
   call near(reshape(paths,[18]), ref_rDefaultSpread, 'rDefaultSpread')
   call esg_martingale(scenario, means, deviation, status)
   call near([deviation], ref_martingale, 'MartingaleTest')
   call esg_zc_bond(scenario, 0, 3.25_dp, value, status)
   call check(status == 0, 'Zero-Coupon at 0 status')
   call near(value, ref_price_0_1, 'Zero-Coupon at 0')
   call esg_asset_price(scenario, 'Zero-Coupon', 0, 3.25_dp, value, status)
   call near(value, ref_dispatch_0_1, 'dispatch Zero-Coupon at 0')
   call esg_bond(scenario, 0, 3.25_dp, 2, 0.025_dp, value, status)
   call check(status == 0, 'Bond at 0 status')
   call near(value, ref_price_0_2, 'Bond at 0')
   call esg_asset_price(scenario, 'Bond', 0, 3.25_dp, value, status, &
      coupons=2, coupon_rate=0.025_dp)
   call near(value, ref_dispatch_0_2, 'dispatch Bond at 0')
   call esg_corporate_bond(scenario, 0, 3.25_dp, 2, 0.025_dp, 0.4_dp, value, status)
   call check(status == 0, 'CBond at 0 status')
   call near(value, ref_price_0_3, 'CBond at 0')
   call esg_asset_price(scenario, 'CBond', 0, 3.25_dp, value, status, &
      coupons=2, coupon_rate=0.025_dp, recovery=0.4_dp)
   call near(value, ref_dispatch_0_3, 'dispatch CBond at 0')
   call esg_convertible_bond(scenario, 0, 3.25_dp, 2, 0.025_dp, value, status)
   call check(status == 0, 'ConvBond at 0 status')
   call near(value, ref_price_0_4, 'ConvBond at 0')
   call esg_asset_price(scenario, 'ConvBond', 0, 3.25_dp, value, status, &
      coupons=2, coupon_rate=0.025_dp)
   call near(value, ref_dispatch_0_4, 'dispatch ConvBond at 0')
   call esg_stock_call(scenario, 0, 3.25_dp, 1.05_dp, value, status)
   call check(status == 0, 'EuroCall_Stock at 0 status')
   call near(value, ref_price_0_5, 'EuroCall_Stock at 0')
   call esg_asset_price(scenario, 'EuroCall_Stock', 0, 3.25_dp, value, status, &
      strike=1.05_dp)
   call near(value, ref_dispatch_0_5, 'dispatch EuroCall_Stock at 0')
   call esg_stock_put(scenario, 0, 3.25_dp, 1.05_dp, value, status)
   call check(status == 0, 'EuroPut_Stock at 0 status')
   call near(value, ref_price_0_6, 'EuroPut_Stock at 0')
   call esg_asset_price(scenario, 'EuroPut_Stock', 0, 3.25_dp, value, status, &
      strike=1.05_dp)
   call near(value, ref_dispatch_0_6, 'dispatch EuroPut_Stock at 0')
   call esg_zc_call(scenario, 0, 3.25_dp, 4.5_dp, 0.93_dp, value, status)
   call check(status == 0, 'EuroCall_ZC at 0 status')
   call near(value, ref_price_0_7, 'EuroCall_ZC at 0')
   call esg_asset_price(scenario, 'EuroCall_ZC', 0, 3.25_dp, value, status, &
      underlying=4.5_dp, strike=0.93_dp)
   call near(value, ref_dispatch_0_7, 'dispatch EuroCall_ZC at 0')
   call esg_zc_put(scenario, 0, 3.25_dp, 4.5_dp, 0.93_dp, value, status)
   call check(status == 0, 'EuroPut_ZC at 0 status')
   call near(value, ref_price_0_8, 'EuroPut_ZC at 0')
   call esg_asset_price(scenario, 'EuroPut_ZC', 0, 3.25_dp, value, status, &
      underlying=4.5_dp, strike=0.93_dp)
   call near(value, ref_dispatch_0_8, 'dispatch EuroPut_ZC at 0')
   call esg_cds(scenario, 0, 3.25_dp, 0.4_dp, value, status)
   call check(status == 0, 'CDS at 0 status')
   call near(value, ref_price_0_9, 'CDS at 0')
   call esg_asset_price(scenario, 'CDS', 0, 3.25_dp, value, status, &
      recovery=0.4_dp)
   call near(value, ref_dispatch_0_9, 'dispatch CDS at 0')
   call esg_zc_bond(scenario, 1, 3.25_dp, value, status)
   call check(status == 0, 'Zero-Coupon at 1 status')
   call near(value, ref_price_1_1, 'Zero-Coupon at 1')
   call esg_asset_price(scenario, 'Zero-Coupon', 1, 3.25_dp, value, status)
   call near(value, ref_dispatch_1_1, 'dispatch Zero-Coupon at 1')
   call esg_bond(scenario, 1, 3.25_dp, 2, 0.025_dp, value, status)
   call check(status == 0, 'Bond at 1 status')
   call near(value, ref_price_1_2, 'Bond at 1')
   call esg_asset_price(scenario, 'Bond', 1, 3.25_dp, value, status, &
      coupons=2, coupon_rate=0.025_dp)
   call near(value, ref_dispatch_1_2, 'dispatch Bond at 1')
   call esg_corporate_bond(scenario, 1, 3.25_dp, 2, 0.025_dp, 0.4_dp, value, status)
   call check(status == 0, 'CBond at 1 status')
   call near(value, ref_price_1_3, 'CBond at 1')
   call esg_asset_price(scenario, 'CBond', 1, 3.25_dp, value, status, &
      coupons=2, coupon_rate=0.025_dp, recovery=0.4_dp)
   call near(value, ref_dispatch_1_3, 'dispatch CBond at 1')
   call esg_convertible_bond(scenario, 1, 3.25_dp, 2, 0.025_dp, value, status)
   call check(status == 0, 'ConvBond at 1 status')
   call near(value, ref_price_1_4, 'ConvBond at 1')
   call esg_asset_price(scenario, 'ConvBond', 1, 3.25_dp, value, status, &
      coupons=2, coupon_rate=0.025_dp)
   call near(value, ref_dispatch_1_4, 'dispatch ConvBond at 1')
   call esg_stock_call(scenario, 1, 3.25_dp, 1.05_dp, value, status)
   call check(status == 0, 'EuroCall_Stock at 1 status')
   call near(value, ref_price_1_5, 'EuroCall_Stock at 1')
   call esg_asset_price(scenario, 'EuroCall_Stock', 1, 3.25_dp, value, status, &
      strike=1.05_dp)
   call near(value, ref_dispatch_1_5, 'dispatch EuroCall_Stock at 1')
   call esg_stock_put(scenario, 1, 3.25_dp, 1.05_dp, value, status)
   call check(status == 0, 'EuroPut_Stock at 1 status')
   call near(value, ref_price_1_6, 'EuroPut_Stock at 1')
   call esg_asset_price(scenario, 'EuroPut_Stock', 1, 3.25_dp, value, status, &
      strike=1.05_dp)
   call near(value, ref_dispatch_1_6, 'dispatch EuroPut_Stock at 1')
   call esg_zc_call(scenario, 1, 3.25_dp, 4.5_dp, 0.93_dp, value, status)
   call check(status == 0, 'EuroCall_ZC at 1 status')
   call near(value, ref_price_1_7, 'EuroCall_ZC at 1')
   call esg_asset_price(scenario, 'EuroCall_ZC', 1, 3.25_dp, value, status, &
      underlying=4.5_dp, strike=0.93_dp)
   call near(value, ref_dispatch_1_7, 'dispatch EuroCall_ZC at 1')
   call esg_zc_put(scenario, 1, 3.25_dp, 4.5_dp, 0.93_dp, value, status)
   call check(status == 0, 'EuroPut_ZC at 1 status')
   call near(value, ref_price_1_8, 'EuroPut_ZC at 1')
   call esg_asset_price(scenario, 'EuroPut_ZC', 1, 3.25_dp, value, status, &
      underlying=4.5_dp, strike=0.93_dp)
   call near(value, ref_dispatch_1_8, 'dispatch EuroPut_ZC at 1')
   call esg_cds(scenario, 1, 3.25_dp, 0.4_dp, value, status)
   call check(status == 0, 'CDS at 1 status')
   call near(value, ref_price_1_9, 'CDS at 1')
   call esg_asset_price(scenario, 'CDS', 1, 3.25_dp, value, status, &
      recovery=0.4_dp)
   call near(value, ref_dispatch_1_9, 'dispatch CDS at 1')
   call esg_zc_bond(scenario, 2, 3.25_dp, value, status)
   call check(status == 0, 'Zero-Coupon at 2 status')
   call near(value, ref_price_2_1, 'Zero-Coupon at 2')
   call esg_asset_price(scenario, 'Zero-Coupon', 2, 3.25_dp, value, status)
   call near(value, ref_dispatch_2_1, 'dispatch Zero-Coupon at 2')
   call esg_bond(scenario, 2, 3.25_dp, 2, 0.025_dp, value, status)
   call check(status == 0, 'Bond at 2 status')
   call near(value, ref_price_2_2, 'Bond at 2')
   call esg_asset_price(scenario, 'Bond', 2, 3.25_dp, value, status, &
      coupons=2, coupon_rate=0.025_dp)
   call near(value, ref_dispatch_2_2, 'dispatch Bond at 2')
   call esg_corporate_bond(scenario, 2, 3.25_dp, 2, 0.025_dp, 0.4_dp, value, status)
   call check(status == 0, 'CBond at 2 status')
   call near(value, ref_price_2_3, 'CBond at 2')
   call esg_asset_price(scenario, 'CBond', 2, 3.25_dp, value, status, &
      coupons=2, coupon_rate=0.025_dp, recovery=0.4_dp)
   call near(value, ref_dispatch_2_3, 'dispatch CBond at 2')
   call esg_convertible_bond(scenario, 2, 3.25_dp, 2, 0.025_dp, value, status)
   call check(status == 0, 'ConvBond at 2 status')
   call near(value, ref_price_2_4, 'ConvBond at 2')
   call esg_asset_price(scenario, 'ConvBond', 2, 3.25_dp, value, status, &
      coupons=2, coupon_rate=0.025_dp)
   call near(value, ref_dispatch_2_4, 'dispatch ConvBond at 2')
   call esg_stock_call(scenario, 2, 3.25_dp, 1.05_dp, value, status)
   call check(status == 0, 'EuroCall_Stock at 2 status')
   call near(value, ref_price_2_5, 'EuroCall_Stock at 2')
   call esg_asset_price(scenario, 'EuroCall_Stock', 2, 3.25_dp, value, status, &
      strike=1.05_dp)
   call near(value, ref_dispatch_2_5, 'dispatch EuroCall_Stock at 2')
   call esg_stock_put(scenario, 2, 3.25_dp, 1.05_dp, value, status)
   call check(status == 0, 'EuroPut_Stock at 2 status')
   call near(value, ref_price_2_6, 'EuroPut_Stock at 2')
   call esg_asset_price(scenario, 'EuroPut_Stock', 2, 3.25_dp, value, status, &
      strike=1.05_dp)
   call near(value, ref_dispatch_2_6, 'dispatch EuroPut_Stock at 2')
   call esg_zc_call(scenario, 2, 3.25_dp, 4.5_dp, 0.93_dp, value, status)
   call check(status == 0, 'EuroCall_ZC at 2 status')
   call near(value, ref_price_2_7, 'EuroCall_ZC at 2')
   call esg_asset_price(scenario, 'EuroCall_ZC', 2, 3.25_dp, value, status, &
      underlying=4.5_dp, strike=0.93_dp)
   call near(value, ref_dispatch_2_7, 'dispatch EuroCall_ZC at 2')
   call esg_zc_put(scenario, 2, 3.25_dp, 4.5_dp, 0.93_dp, value, status)
   call check(status == 0, 'EuroPut_ZC at 2 status')
   call near(value, ref_price_2_8, 'EuroPut_ZC at 2')
   call esg_asset_price(scenario, 'EuroPut_ZC', 2, 3.25_dp, value, status, &
      underlying=4.5_dp, strike=0.93_dp)
   call near(value, ref_dispatch_2_8, 'dispatch EuroPut_ZC at 2')
   call esg_cds(scenario, 2, 3.25_dp, 0.4_dp, value, status)
   call check(status == 0, 'CDS at 2 status')
   call near(value, ref_price_2_9, 'CDS at 2')
   call esg_asset_price(scenario, 'CDS', 2, 3.25_dp, value, status, &
      recovery=0.4_dp)
   call near(value, ref_dispatch_2_9, 'dispatch CDS at 2')
   call esg_simulate_asset(p, zero, innovations(:,:4,:), 'Zero-Coupon', 1, 3.25_dp, value, status)
   call near(value, ref_rAssetDistribution, 'rAssetDistribution')
   call esg_stock_call(scenario, 2, 2.0_dp, 1.05_dp, value, status)
   call near(value, spread(ref_maturity_stock(1), 1, 3), 'stock maturity convention')
   call esg_cds(scenario, 2, 2.0_dp, 0.4_dp, value, status)
   call near(value, spread(ref_maturity_cds(1), 1, 3), 'CDS maturity convention')
   print *, 'PASS: genuine ESG 1.3 reference cases, 20 exported operations'
contains
   subroutine check(condition, label)
      logical, intent(in) :: condition !! Required true assertion for this reference case.
      character(*), intent(in) :: label !! Human-readable case identifying a failed assertion.
      if (.not. condition) then
         print *, 'FAIL ', label
         error stop 1
      end if
   end subroutine
   subroutine near(actual, expected, label)
      real(dp), intent(in) :: actual(:) !! Fortran result flattened in column-major order.
      real(dp), intent(in) :: expected(:) !! Genuine upstream R result in matching order.
      character(*), intent(in) :: label !! Name of the tested R operation or fixture.
      call check(size(actual) == size(expected), label // ' shape')
      if (.not. all(abs(actual - expected) <= 2.0e-10_dp * max(1.0_dp, abs(expected)))) then
         print *, label, ' actual ', actual, ' expected ', expected
         error stop 2
      end if
   end subroutine
end program
