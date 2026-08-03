! SPDX-License-Identifier: GPL-3.0-or-later
program curve_calibration_example
   use quant_bond_curves
   implicit none
   type(qbc_bond) :: bonds(4)
   type(qbc_calibration_result) :: fit
   real(dp) :: terms(4), true_rates(4), prices(4), initial(4)
   integer :: i

   terms = [1.0_dp,2.0_dp,3.0_dp,5.0_dp]
   true_rates = [0.030_dp,0.035_dp,0.040_dp,0.045_dp]
   initial = 0.025_dp
   do i=1,4
      bonds(i)%analysis_date=make_date(2024,1,1)
      bonds(i)%maturity=make_date(2024+nint(terms(i)),1,1)
      bonds(i)%frequency=0
      bonds(i)%asset_type=qbc_asset_fixed
      bonds(i)%rate_type=qbc_rate_continuous
      prices(i)=exp(-true_rates(i)*terms(i))
   end do
   call bootstrap_curve(bonds,prices,terms,initial,fit)
   write(*,'(a,es14.6)') 'objective: ',fit%objective
   do i=1,size(terms)
      write(*,'(f6.2,2f12.7)') terms(i),true_rates(i),fit%curve%rates(i)
   end do
end program curve_calibration_example
