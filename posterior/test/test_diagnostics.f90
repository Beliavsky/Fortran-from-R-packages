! SPDX-License-Identifier: BSD-3-Clause
program test_diagnostics
   use, intrinsic :: ieee_arithmetic, only: ieee_is_nan, ieee_is_finite, ieee_value, &
      ieee_quiet_nan, ieee_positive_inf, ieee_negative_inf
   use posterior_diagnostics
   implicit none
   integer :: unit, ncase, case_id, n, m, j, checks
   character(80) :: label
   real(dp), allocatable :: x(:, :)
   real(dp) :: actual(20), expected(20), ranks(5), inf, nan, q(3), ac(4)
   checks = 0
   open(newunit = unit, file = 'test/reference_diagnostics.dat', status = 'old', action = 'read')
   read(unit, *) ncase
   do case_id = 1, ncase
      read(unit, '(a)') label
      read(unit, *) n, m
      allocate(x(n, m))
      read(unit, *) x
      read(unit, *) expected
      actual(1) = rhat_basic(x)
      actual(2) = rhat_basic(x, split = .false.)
      actual(3) = rhat(x)
      actual(4) = ess_basic(x)
      actual(5) = ess_basic(x, split = .false.)
      actual(6) = ess_bulk(x)
      actual(7) = ess_tail(x)
      actual(8) = ess_mean(x)
      actual(9) = ess_median(x)
      actual(10) = ess_sd(x)
      actual(11) = mcse_mean(x)
      actual(12) = mcse_median(x)
      actual(13) = mcse_sd(x)
      actual(14:16) = ess_quantile(x, [0.0_dp, 0.3_dp, 1.0_dp])
      actual(17:19) = mcse_quantile(x, [0.0_dp, 0.3_dp, 1.0_dp])
      actual(20) = rhat_nested(x, [(j, j = 1, m)])
      do j = 1, 20
         call check_number(actual(j), expected(j), 2.0e-10_dp, trim(label), j)
      end do
      deallocate(x)
   end do
   close(unit)
   nan = ieee_value(0.0_dp, ieee_quiet_nan)
   inf = ieee_value(0.0_dp, ieee_positive_inf)
   ranks = r_scale([3.0_dp, 1.0_dp, 1.0_dp, inf, nan])
   call check_number(ranks(1), 3.0_dp, 0.0_dp, 'tie ranks', 1)
   call check_number(ranks(2), 1.5_dp, 0.0_dp, 'tie ranks', 2)
   call check_number(ranks(4), 4.0_dp, 0.0_dp, 'infinite rank', 4)
   call check_number(ranks(5), nan, 0.0_dp, 'missing rank', 5)
   q = quantile2([-inf, 1.0_dp, inf], [0.0_dp, 0.5_dp, 1.0_dp])
   call check_number(q(1), -inf, 0.0_dp, 'quantile endpoint', 1)
   call check_number(q(2), 1.0_dp, 0.0_dp, 'quantile endpoint', 2)
   call check_number(q(3), inf, 0.0_dp, 'quantile endpoint', 3)
   ac = autocovariance([1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp])
   call check_number(ac(1), 1.25_dp, 1.0e-14_dp, 'analytic covariance', 1)
   call check_number(ac(2), 0.3125_dp, 1.0e-14_dp, 'analytic covariance', 2)
   call check_number(ac(3), -0.375_dp, 1.0e-14_dp, 'analytic covariance', 3)
   call check_number(ac(4), -0.5625_dp, 1.0e-14_dp, 'analytic covariance', 4)
   if (.not. is_constant([2.0_dp, 2.0_dp])) error stop 'Constant test failed'
   if (is_constant([inf, inf])) error stop 'Infinite repeated values must not be constant'
   if (is_constant([nan, nan])) error stop 'Missing repeated values must not be constant'
   if (is_constant([huge(1.0_dp), -huge(1.0_dp)])) error stop 'Overflowing range must not be constant'
   if (is_constant([2.0_dp, 2.0_dp], tol = -1.0_dp)) error stop 'Negative tolerance must be rejected'
   if (is_constant([2.0_dp, 2.0_dp], tol = inf)) error stop 'Infinite tolerance must be rejected'
   if (is_constant([2.0_dp, 2.0_dp], tol = nan)) error stop 'Missing tolerance must be rejected'
   if (is_constant([2.0_dp, 2.0_dp], tol = 0.0_dp)) error stop 'Constant comparison must be strict'

   print '(a,i0,a)', 'posterior diagnostics: ', checks, ' R-oracle and analytic checks passed'
contains
   subroutine check_number(actual, expected, tolerance, context, index)
      real(dp), intent(in) :: actual !! Computed diagnostic or numerical property.
      real(dp), intent(in) :: expected !! R oracle or exact analytic expected value.
      real(dp), intent(in) :: tolerance !! Relative error tolerance with unit absolute floor.
      character(*), intent(in) :: context !! Test-case description for failure reporting.
      integer, intent(in) :: index !! Position of the checked diagnostic within the case.
      logical :: ok
      if (ieee_is_nan(expected)) then
         ok = ieee_is_nan(actual)
      else if (.not. ieee_is_finite(expected)) then
         ok = actual == expected
      else
         ok = ieee_is_finite(actual)
         if (ok) ok = abs(actual - expected) <= tolerance * max(1.0_dp, abs(expected))
      end if
      if (.not. ok) then
         print *, context, index, ' actual=', actual, ' expected=', expected
         error stop 'Numerical comparison failed'
      end if
      checks = checks + 1
   end subroutine check_number
end program test_diagnostics
