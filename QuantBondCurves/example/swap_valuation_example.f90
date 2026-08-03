! SPDX-License-Identifier: GPL-3.0-or-later
program swap_valuation_example
   use quant_bond_curves
   implicit none
   type(qbc_curve) :: curve
   type(qbc_swap) :: swap
   real(dp) :: value

   allocate(curve%terms(4),curve%rates(4))
   curve%terms=[0.0_dp,1.0_dp,3.0_dp,5.0_dp]
   curve%rates=[0.035_dp,0.037_dp,0.040_dp,0.042_dp]
   curve%approximation=2
   curve%rate_type=qbc_rate_continuous

   swap%analysis_date=make_date(2024,1,1)
   swap%maturity=make_date(2029,1,1)
   swap%frequency=2
   swap%rate_type=qbc_rate_continuous
   swap%legs=qbc_leg_fixed_float
   swap%coupon_rate1=0.045_dp
   swap%spread2=0.001_dp
   swap%principal1=1.0e6_dp
   swap%principal2=1.0e6_dp
   value=valuation_swaps(swap,curve,curve)
   write(*,'(a,f16.2)') 'fixed-minus-floating value: ',value
end program swap_valuation_example
