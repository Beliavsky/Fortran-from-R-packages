! SPDX-License-Identifier: GPL-2.0-or-later
! Derived from the R package roll by Jason Foster; see PROVENANCE.md and NOTICE.md.
program test_roll
   use, intrinsic :: ieee_arithmetic, only : ieee_is_nan, ieee_value, ieee_quiet_nan
   use roll, only : dp, roll_na_logical, roll_lm_result
   use roll, only : roll_any, roll_all, roll_sum, roll_prod, roll_mean, roll_min, roll_max
   use roll, only : roll_idxmin, roll_idxmax, roll_median, roll_quantile, roll_var, roll_sd
   use roll, only : roll_scale, roll_cov, roll_cor, roll_crossprod, roll_lm
   implicit none

   real(dp) :: x(5), y(5), w(3), nanv
   real(dp) :: xm(5, 2), ym(5, 1)
   real(dp), allocatable :: r(:), m(:, :), cube(:, :, :)
   integer :: logical_x(5)
   integer, allocatable :: ir(:)
   type(roll_lm_result) :: fit

   nanv = ieee_value(0.0_dp, ieee_quiet_nan)
   x = [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 5.0_dp]
   y = 1.0_dp + 2.0_dp * x
   w = [0.25_dp, 0.5_dp, 1.0_dp]
   logical_x = [1, 0, roll_na_logical, 0, 1]

   call roll_sum(x, 3, r, min_obs=1)
   call check_vec(r, [1.0_dp, 3.0_dp, 6.0_dp, 9.0_dp, 12.0_dp], 1.0e-13_dp, 'sum')

   call roll_sum(x, 3, r, weights=w, min_obs=1)
   call assert_close(r(3), 4.25_dp, 1.0e-13_dp, 'weighted sum')

   call roll_prod(x, 3, r, min_obs=1)
   call check_vec(r, [1.0_dp, 2.0_dp, 6.0_dp, 24.0_dp, 60.0_dp], 1.0e-13_dp, 'product')

   call roll_mean(x, 3, r, min_obs=1)
   call check_vec(r, [1.0_dp, 1.5_dp, 2.0_dp, 3.0_dp, 4.0_dp], 1.0e-13_dp, 'mean')

   call roll_min(x, 3, r, min_obs=1)
   call check_vec(r, [1.0_dp, 1.0_dp, 1.0_dp, 2.0_dp, 3.0_dp], 1.0e-13_dp, 'min')
   call roll_max(x, 3, r, min_obs=1)
   call check_vec(r, [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 5.0_dp], 1.0e-13_dp, 'max')

   call roll_idxmin(x, 3, ir, min_obs=1)
   call check_int(ir, [1, 1, 1, 1, 1], 'idxmin')
   call roll_idxmax(x, 3, ir, min_obs=1)
   call check_int(ir, [1, 2, 3, 3, 3], 'idxmax')

   call roll_median(x, 3, r, min_obs=1)
   call check_vec(r, [1.0_dp, 1.5_dp, 2.0_dp, 3.0_dp, 4.0_dp], 1.0e-13_dp, 'median')
   call roll_quantile(x, 3, 0.25_dp, r, min_obs=1)
   call check_vec(r, [1.0_dp, 1.0_dp, 1.0_dp, 2.0_dp, 3.0_dp], 1.0e-13_dp, 'quantile')

   call roll_var(x, 3, r, min_obs=1)
   call assert_true(ieee_is_nan(r(1)), 'variance one observation is NA')
   call check_vec(r(2:5), [0.5_dp, 1.0_dp, 1.0_dp, 1.0_dp], 1.0e-13_dp, 'variance')
   call roll_sd(x, 3, r, min_obs=1)
   call assert_close(r(2), sqrt(0.5_dp), 1.0e-13_dp, 'sd')
   call roll_scale(x, 3, r, min_obs=1)
   call assert_true(ieee_is_nan(r(1)), 'scale one observation is NA')
   call assert_close(r(3), 1.0_dp, 1.0e-13_dp, 'scale')

   call roll_any(logical_x, 3, ir, min_obs=1)
   call check_int(ir, [1, 1, 1, roll_na_logical, 1], 'any')
   call roll_all(logical_x, 3, ir, min_obs=1)
   call check_int(ir, [1, 0, 0, 0, 0], 'all')

   call roll_cov(x, 3, r, y=2.0_dp*x, min_obs=1)
   call assert_true(ieee_is_nan(r(1)), 'cov one observation is NA')
   call check_vec(r(2:5), [1.0_dp, 2.0_dp, 2.0_dp, 2.0_dp], 1.0e-13_dp, 'covariance')
   call roll_cor(x, 3, r, y=2.0_dp*x, min_obs=1)
   call check_vec(r(2:5), [1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp], 1.0e-13_dp, 'correlation')
   call roll_crossprod(x, 3, r, y=2.0_dp*x, min_obs=1)
   call check_vec(r, [2.0_dp, 10.0_dp, 28.0_dp, 58.0_dp, 100.0_dp], 1.0e-13_dp, 'crossproduct')

   xm(:, 1) = x
   xm(:, 2) = 10.0_dp + x
   xm(3, 2) = nanv
   call roll_mean(xm, 3, m, min_obs=1, complete_obs=.true.)
   call assert_close(m(3, 1), 1.5_dp, 1.0e-13_dp, 'complete-row matrix mean')
   call roll_mean(xm, 3, m, min_obs=1, complete_obs=.false.)
   call assert_close(m(3, 1), 2.0_dp, 1.0e-13_dp, 'pairwise matrix mean')

   ym(:, 1) = y
   call roll_cov(reshape(x, [5, 1]), 3, cube, y=ym, min_obs=1)
   call assert_close(cube(1, 1, 5), 2.0_dp, 1.0e-13_dp, 'matrix covariance cube')

   call roll_lm(x, y, 4, fit, min_obs=3)
   call assert_close(fit%coefficients(3, 1, 1), 1.0_dp, 1.0e-11_dp, 'lm intercept')
   call assert_close(fit%coefficients(3, 2, 1), 2.0_dp, 1.0e-11_dp, 'lm slope')
   call assert_close(fit%r_squared(3, 1), 1.0_dp, 1.0e-11_dp, 'lm r squared')
   call assert_close(fit%std_error(3, 2, 1), 0.0_dp, 1.0e-10_dp, 'lm standard error')

   write (*, '(a)') 'All roll tests passed.'

contains

   subroutine assert_true(condition, label)
      logical, intent(in) :: condition !! Condition that must evaluate true.
      character(len=*), intent(in) :: label !! Human-readable assertion label.
      if (.not. condition) then
         write (*, '(a)') 'FAIL: ' // trim(label)
         error stop 1
      end if
   end subroutine assert_true

   subroutine assert_close(actual, expected, tolerance, label)
      real(dp), intent(in) :: actual !! Computed scalar value.
      real(dp), intent(in) :: expected !! Expected scalar value.
      real(dp), intent(in) :: tolerance !! Maximum absolute error allowed.
      character(len=*), intent(in) :: label !! Human-readable assertion label.
      if (ieee_is_nan(actual) .or. abs(actual - expected) > tolerance) then
         write (*, '(a,2(1x,es24.16))') 'FAIL: ' // trim(label), actual, expected
         error stop 1
      end if
   end subroutine assert_close

   subroutine check_vec(actual, expected, tolerance, label)
      real(dp), intent(in) :: actual(:) !! Computed vector.
      real(dp), intent(in) :: expected(:) !! Expected vector of the same size.
      real(dp), intent(in) :: tolerance !! Maximum elementwise absolute error allowed.
      character(len=*), intent(in) :: label !! Human-readable assertion label.
      integer :: i
      call assert_true(size(actual) == size(expected), trim(label) // ' size')
      do i = 1, size(actual)
         call assert_close(actual(i), expected(i), tolerance, trim(label) // ' value')
      end do
   end subroutine check_vec

   subroutine check_int(actual, expected, label)
      integer, intent(in) :: actual(:) !! Computed integer vector.
      integer, intent(in) :: expected(:) !! Expected integer vector of the same size.
      character(len=*), intent(in) :: label !! Human-readable assertion label.
      call assert_true(size(actual) == size(expected), trim(label) // ' size')
      call assert_true(all(actual == expected), trim(label) // ' values')
   end subroutine check_int

end program test_roll
