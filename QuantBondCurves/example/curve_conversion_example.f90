! SPDX-License-Identifier: GPL-3.0-or-later
program curve_conversion_example
   use quant_bond_curves
   implicit none
   real(dp), allocatable :: forward(:), recovered(:), t1(:), t2(:)
   real(dp) :: terms(6), spot(6)
   integer :: i

   terms = [0.0_dp,0.5_dp,1.0_dp,2.0_dp,5.0_dp,10.0_dp]
   spot = [0.025_dp,0.027_dp,0.030_dp,0.034_dp,0.040_dp,0.043_dp]
   call spot2forward(terms,spot,2,forward,t1)
   call fwd2spot(t1,forward,1,recovered,t2)
   write(*,'(a)') ' term       spot    forward  recovered'
   do i=1,size(t1)
      write(*,'(f6.2,3f11.6)') t1(i),spot(i),forward(i),recovered(i)
   end do
end program curve_conversion_example
