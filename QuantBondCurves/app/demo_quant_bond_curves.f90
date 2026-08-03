! SPDX-License-Identifier: GPL-3.0-or-later
program demo_quant_bond_curves
   use quant_bond_curves
   implicit none
   type(qbc_calibration_result) :: market_curve
   type(qbc_bond) :: bond
   type(qbc_bond_sensitivity) :: risk
   real(dp) :: nodes(6), terms(5), yields(5), price, ytm
   integer :: i

   terms=[0.25_dp,0.5_dp,1.0_dp,2.0_dp,5.0_dp]
   yields=[0.030_dp,0.032_dp,0.035_dp,0.039_dp,0.044_dp]
   nodes=[0.0_dp,0.5_dp,1.0_dp,2.0_dp,3.0_dp,5.0_dp]
   call curve_calibration(terms,yields,nodes,market_curve,approximation=2,rate_type=qbc_rate_continuous)

   write(*,'(a)') 'Calibrated/interpolated continuous spot curve'
   do i=1,size(nodes)
      write(*,'(f7.2,f12.7)') market_curve%curve%terms(i),market_curve%curve%rates(i)
   end do

   bond%analysis_date=make_date(2024,1,1)
   bond%maturity=make_date(2029,1,1)
   bond%coupon_rate=0.04_dp
   bond%principal=100.0_dp
   bond%frequency=2
   bond%asset_type=qbc_asset_fixed
   bond%rate_type=qbc_rate_continuous
   price=valuation_bonds(bond,[bond%coupon_rate],market_curve%curve)
   ytm=bond_price2rate(bond,[bond%coupon_rate],price)
   risk=bond_sensitivity(bond,[bond%coupon_rate],ytm)
   write(*,'(/,a,f12.6)') 'bond price: ',price
   write(*,'(a,f12.6)') 'flat-equivalent yield: ',ytm
   write(*,'(a,f12.6)') 'modified duration: ',risk%modified_duration
end program demo_quant_bond_curves
