! SPDX-License-Identifier: GPL-3.0-or-later
program bond_pricing_example
   use quant_bond_curves
   implicit none
   type(qbc_bond) :: bond
   type(qbc_bond_sensitivity) :: risk
   real(dp) :: price, yield

   bond%analysis_date = make_date(2022,6,1)
   bond%maturity = make_date(2026,6,1)
   bond%coupon_rate = 0.06_dp
   bond%principal = 100.0_dp
   bond%frequency = 1
   bond%asset_type = qbc_asset_tes
   bond%business_convention = 'N'

   price = valuation_bonds(bond, [bond%coupon_rate], [0.08_dp])
   yield = bond_price2rate(bond, [bond%coupon_rate], price)
   risk = bond_sensitivity(bond, [bond%coupon_rate], yield)

   write(*,'(a,f12.6)') 'dirty price:       ', price
   write(*,'(a,f12.6)') 'yield:             ', yield
   write(*,'(a,f12.6)') 'modified duration: ', risk%modified_duration
   write(*,'(a,f12.6)') 'convexity:         ', risk%convexity
   write(*,'(a,f12.6)') 'DV01:              ', risk%dv01
end program bond_pricing_example
