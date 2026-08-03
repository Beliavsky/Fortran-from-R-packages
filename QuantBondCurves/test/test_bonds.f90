! SPDX-License-Identifier: GPL-3.0-or-later
program test_bonds
   use quant_bond_curves
   implicit none
   type(qbc_bond) :: bond
   type(qbc_bond_sensitivity) :: sens
   real(dp) :: price, expected, recovered, wal
   integer :: st

   bond%analysis_date = make_date(2022,6,1)
   bond%maturity = make_date(2026,6,1)
   bond%coupon_rate = 0.06_dp
   bond%principal = 1.0_dp
   bond%frequency = 1
   bond%asset_type = qbc_asset_tes
   bond%rate_type = qbc_rate_discrete
   bond%business_convention = 'N'
   expected = sum([0.06_dp/1.08_dp, 0.06_dp/1.08_dp**2, 0.06_dp/1.08_dp**3, 1.06_dp/1.08_dp**4])
   price = valuation_bonds(bond, [0.06_dp], [0.08_dp], status=st)
   call assert_close(price, expected, 1.0e-12_dp, 'bond price')
   recovered = bond_price_to_rate(bond, [0.06_dp], price, status=st)
   call assert_close(recovered, 0.08_dp, 1.0e-9_dp, 'yield inversion')

   sens = bond_sensitivity(bond, [0.06_dp], 0.08_dp, bump=1.0e-5_dp, status=st)
   call assert_true(sens%modified_duration > 3.0_dp .and. sens%modified_duration < 4.0_dp, 'duration range')
   call assert_true(sens%convexity > 10.0_dp, 'positive convexity')
   call assert_close(sens%dv01, sens%modified_duration*price*1.0e-4_dp, 1.0e-14_dp, 'DV01 identity')

   wal = average_life(bond, [0.06_dp], [0.08_dp], discounted=.false., status=st)
   call assert_true(wal > 3.5_dp .and. wal < 4.0_dp, 'weighted average life')

   bond%analysis_date = make_date(2022,9,1)
   call assert_true(dirty_to_clean(bond, 1.0_dp, 0.06_dp, st) < 1.0_dp, 'clean below dirty between coupons')

   print '(a)', 'test_bonds: PASS'
contains
   subroutine assert_true(condition, label)
      logical, intent(in) :: condition
      character(len=*), intent(in) :: label
      if (.not. condition) error stop label
   end subroutine
   subroutine assert_close(x, y, tol, label)
      real(dp), intent(in) :: x, y, tol
      character(len=*), intent(in) :: label
      if (abs(x-y) > tol) then
         write(*,*) label, x, y
         error stop 'assert_close'
      end if
   end subroutine
end program test_bonds
