! SPDX-License-Identifier: GPL-3.0-or-later
program basis_curve_example
   use quant_bond_curves
   implicit none
   type(qbc_swap) :: swaps(3)
   type(qbc_curve) :: local_curve, foreign_curve, true_basis
   type(qbc_calibration_result) :: fit
   real(dp) :: terms(3), market_values(3)
   integer :: i

   terms=[1.0_dp,2.0_dp,3.0_dp]
   allocate(local_curve%terms(3),local_curve%rates(3),foreign_curve%terms(3),foreign_curve%rates(3), &
            true_basis%terms(3),true_basis%rates(3))
   local_curve%terms=terms; local_curve%rates=0.04_dp; local_curve%approximation=2; local_curve%rate_type=qbc_rate_continuous
   foreign_curve=local_curve
   true_basis%terms=terms; true_basis%rates=[0.02_dp,0.025_dp,0.03_dp]
   true_basis%approximation=2; true_basis%rate_type=qbc_rate_continuous
   do i=1,3
      swaps(i)%analysis_date=make_date(2024,1,1)
      swaps(i)%maturity=make_date(2024+i,1,1)
      swaps(i)%frequency=0
      swaps(i)%rate_type=qbc_rate_continuous
      swaps(i)%principal1=1.0_dp; swaps(i)%principal2=1.0_dp
      market_values(i)=valuation_swaps(swaps(i),local_curve,foreign_curve,true_basis)
   end do
   call basis_curve(swaps,local_curve,foreign_curve,market_values,terms,[0.01_dp,0.01_dp,0.01_dp],fit)
   write(*,'(a,es14.6)') 'objective: ',fit%objective
   do i=1,3
      write(*,'(f6.2,f12.7)') terms(i),fit%curve%rates(i)
   end do
end program basis_curve_example
