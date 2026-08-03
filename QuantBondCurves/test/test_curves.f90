! SPDX-License-Identifier: GPL-3.0-or-later
program test_curves
   use quant_bond_curves
   implicit none
   real(dp), allocatable :: f(:), s(:), t1(:), t2(:), vals(:), dfs(:)
   real(dp) :: terms(5), spot(5)
   type(qbc_curve) :: curve
   integer :: st

   terms = [0.0_dp, 0.5_dp, 1.0_dp, 2.0_dp, 3.0_dp]
   spot = [0.02_dp, 0.022_dp, 0.025_dp, 0.03_dp, 0.035_dp]
   call spot_to_forward(terms, spot, 2, f, t1, st)
   call assert_true(st == qbc_success, 'spot_to_forward status')
   call forward_to_spot(t1, f, 1, s, t2, st)
   call assert_close(maxval(abs(s-spot)), 0.0_dp, 1.0e-12_dp, 'linear spot/constant forward round trip')

   allocate(curve%terms(3), curve%rates(3))
   curve%terms = [0.0_dp, 1.0_dp, 2.0_dp]
   curve%rates = [0.02_dp, 0.04_dp, 0.06_dp]
   curve%approximation = 2
   call interpolate_curve(curve, [0.5_dp, 1.5_dp, 3.0_dp], vals, st)
   call assert_close(maxval(abs(vals-[0.03_dp,0.05_dp,0.06_dp])), 0.0_dp, 1.0e-14_dp, 'linear interpolation and flat extrapolation')

   call discount_factors([make_date(2025,1,1), make_date(2026,1,1)], [0.05_dp], make_date(2024,1,1), &
                         qbc_rate_continuous, 1, dfs, st)
   call assert_close(dfs(1), exp(-0.05_dp), 1.0e-14_dp, 'continuous discount')
   call assert_close(dfs(2), exp(-0.10_dp), 1.0e-14_dp, 'continuous discount 2')

   print '(a)', 'test_curves: PASS'
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
end program test_curves
