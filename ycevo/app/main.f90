program ycevo_example
   use ycevo, only : dp, ycevo_success, bond_panel_t, yield_curve_t
   use ycevo, only : simulate_bond_panel, estimate_yield, write_yield_curve_csv
   implicit none

   type(bond_panel_t) :: bonds
   type(yield_curve_t) :: curve
   real(dp) :: tau(7), ht(7)
   integer :: i, status
   character(len=256) :: message

   call simulate_bond_panel(bonds, nday=40, n_bonds=60, seed=12345, &
                            max_maturity_years=8.0_dp)
   tau = [0.5_dp,1.0_dp,2.0_dp,3.0_dp,4.0_dp,5.0_dp,7.0_dp]
   ht = [0.35_dp,0.45_dp,0.60_dp,0.70_dp,0.80_dp,0.90_dp,1.00_dp]
   call estimate_yield(bonds, 0.5_dp, 0.35_dp, tau, ht, curve, status, message)
   if (status /= ycevo_success) then
      write(*, '(a)') 'Estimation failed: '//trim(message)
      error stop 1
   end if

   write(*, '(a)') ' tau          discount       yield'
   do i = 1, size(tau)
      write(*, '(f7.3,2(2x,f12.8))') curve%tau(i), curve%discount(i), curve%yield(i)
   end do
   call write_yield_curve_csv('yield_curve.csv', curve, status)
end program ycevo_example
