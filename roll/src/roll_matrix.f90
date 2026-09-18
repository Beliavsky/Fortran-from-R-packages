! SPDX-License-Identifier: GPL-2.0-or-later
! Derived from the R package roll by Jason Foster; see PROVENANCE.md and NOTICE.md.
module roll_matrix
   use, intrinsic :: ieee_arithmetic, only : ieee_is_nan, ieee_value, ieee_quiet_nan
   use r_kinds, only : dp
   use roll_univariate, only : roll_na_logical
   use roll_univariate, only : roll_any_vec, roll_all_vec, roll_sum_vec, roll_prod_vec, roll_mean_vec
   use roll_univariate, only : roll_min_vec, roll_max_vec, roll_idxmin_vec, roll_idxmax_vec
   use roll_univariate, only : roll_median_vec, roll_quantile_vec, roll_var_vec, roll_sd_vec, roll_scale_vec
   implicit none
   private
   public :: roll_any_mat, roll_all_mat, roll_sum_mat, roll_prod_mat, roll_mean_mat
   public :: roll_min_mat, roll_max_mat, roll_idxmin_mat, roll_idxmax_mat
   public :: roll_median_mat, roll_quantile_mat, roll_var_mat, roll_sd_mat, roll_scale_mat

contains

   pure function na_real() result(value)
      real(dp) :: value
      value = ieee_value(0.0_dp, ieee_quiet_nan)
   end function na_real

   pure subroutine masked_real_column(x, column, complete_obs, values)
      real(dp), intent(in) :: x(:, :) !! Source matrix with observations in rows and variables in columns.
      integer, intent(in) :: column !! One-based source column to copy.
      logical, intent(in) :: complete_obs !! Exclude an entire row when any column is missing when true.
      real(dp), allocatable, intent(out) :: values(:) !! Selected column with complete-row exclusions encoded as NaNs.
      integer :: i

      allocate(values(size(x, 1)))
      values = x(:, column)
      if (.not. complete_obs) return
      do i = 1, size(x, 1)
         if (any(ieee_is_nan(x(i, :)))) values(i) = na_real()
      end do
   end subroutine masked_real_column

   pure subroutine masked_logical_column(x, column, complete_obs, values)
      integer, intent(in) :: x(:, :) !! Logical matrix encoded as zero/one/`roll_na_logical`.
      integer, intent(in) :: column !! One-based source column to copy.
      logical, intent(in) :: complete_obs !! Exclude an entire row when any column is missing when true.
      integer, allocatable, intent(out) :: values(:) !! Selected column with complete-row exclusions encoded as missing.
      integer :: i

      allocate(values(size(x, 1)))
      values = x(:, column)
      if (.not. complete_obs) return
      do i = 1, size(x, 1)
         if (any(x(i, :) == roll_na_logical)) values(i) = roll_na_logical
      end do
   end subroutine masked_logical_column

   pure subroutine roll_any_mat(x, width, result, min_obs, complete_obs, na_restore, online)
      integer, intent(in) :: x(:, :) !! Logical matrix encoded as zero/one/`roll_na_logical`.
      integer, intent(in) :: width !! Rolling window width; must be at least one.
      integer, allocatable, intent(out) :: result(:, :) !! Rolling logical any with the same shape as `x`.
      integer, intent(in), optional :: min_obs !! Minimum nonmissing observations; defaults to `width`.
      logical, intent(in), optional :: complete_obs !! Exclude rows missing in any column when true.
      logical, intent(in), optional :: na_restore !! Restore current-cell missing values when true.
      logical, intent(in), optional :: online !! Accepted for API parity; deterministic batch evaluation is used.
      integer :: j
      integer, allocatable :: column(:), temp(:)
      logical :: complete, restore

      complete = .false.
      if (present(complete_obs)) complete = complete_obs
      restore = .false.
      if (present(na_restore)) restore = na_restore
      allocate(result(size(x, 1), size(x, 2)))
      do j = 1, size(x, 2)
         call masked_logical_column(x, j, complete, column)
         call roll_any_vec(column, width, temp, min_obs=min_obs, na_restore=.false., online=online)
         result(:, j) = temp
         if (restore) where (x(:, j) == roll_na_logical) result(:, j) = roll_na_logical
      end do
   end subroutine roll_any_mat

   pure subroutine roll_all_mat(x, width, result, min_obs, complete_obs, na_restore, online)
      integer, intent(in) :: x(:, :) !! Logical matrix encoded as zero/one/`roll_na_logical`.
      integer, intent(in) :: width !! Rolling window width; must be at least one.
      integer, allocatable, intent(out) :: result(:, :) !! Rolling logical all with the same shape as `x`.
      integer, intent(in), optional :: min_obs !! Minimum nonmissing observations; defaults to `width`.
      logical, intent(in), optional :: complete_obs !! Exclude rows missing in any column when true.
      logical, intent(in), optional :: na_restore !! Restore current-cell missing values when true.
      logical, intent(in), optional :: online !! Accepted for API parity; deterministic batch evaluation is used.
      integer :: j
      integer, allocatable :: column(:), temp(:)
      logical :: complete, restore

      complete = .false.
      if (present(complete_obs)) complete = complete_obs
      restore = .false.
      if (present(na_restore)) restore = na_restore
      allocate(result(size(x, 1), size(x, 2)))
      do j = 1, size(x, 2)
         call masked_logical_column(x, j, complete, column)
         call roll_all_vec(column, width, temp, min_obs=min_obs, na_restore=.false., online=online)
         result(:, j) = temp
         if (restore) where (x(:, j) == roll_na_logical) result(:, j) = roll_na_logical
      end do
   end subroutine roll_all_mat

   pure subroutine roll_sum_mat(x, width, result, weights, min_obs, complete_obs, na_restore, online)
      real(dp), intent(in) :: x(:, :) !! Numeric matrix with observations in rows and variables in columns.
      integer, intent(in) :: width !! Rolling window width; must be at least one.
      real(dp), allocatable, intent(out) :: result(:, :) !! Rolling weighted sums with the same shape as `x`.
      real(dp), intent(in), optional :: weights(:) !! Weights whose final element aligns with the current row.
      integer, intent(in), optional :: min_obs !! Minimum usable observations; defaults to `width`.
      logical, intent(in), optional :: complete_obs !! Exclude rows missing in any column when true.
      logical, intent(in), optional :: na_restore !! Restore current-cell NaNs when true.
      logical, intent(in), optional :: online !! Accepted for API parity; deterministic batch evaluation is used.
      call unary_real_matrix(1, x, width, result, weights, min_obs, complete_obs, na_restore, online)
   end subroutine roll_sum_mat

   pure subroutine roll_prod_mat(x, width, result, weights, min_obs, complete_obs, na_restore, online)
      real(dp), intent(in) :: x(:, :) !! Numeric matrix with observations in rows and variables in columns.
      integer, intent(in) :: width !! Rolling window width; must be at least one.
      real(dp), allocatable, intent(out) :: result(:, :) !! Rolling products with the same shape as `x`.
      real(dp), intent(in), optional :: weights(:) !! Weights whose final element aligns with the current row.
      integer, intent(in), optional :: min_obs !! Minimum usable observations; defaults to `width`.
      logical, intent(in), optional :: complete_obs !! Exclude rows missing in any column when true.
      logical, intent(in), optional :: na_restore !! Restore current-cell NaNs when true.
      logical, intent(in), optional :: online !! Accepted for API parity; deterministic batch evaluation is used.
      call unary_real_matrix(2, x, width, result, weights, min_obs, complete_obs, na_restore, online)
   end subroutine roll_prod_mat

   pure subroutine roll_mean_mat(x, width, result, weights, min_obs, complete_obs, na_restore, online)
      real(dp), intent(in) :: x(:, :) !! Numeric matrix with observations in rows and variables in columns.
      integer, intent(in) :: width !! Rolling window width; must be at least one.
      real(dp), allocatable, intent(out) :: result(:, :) !! Rolling weighted means with the same shape as `x`.
      real(dp), intent(in), optional :: weights(:) !! Weights whose final element aligns with the current row.
      integer, intent(in), optional :: min_obs !! Minimum usable observations; defaults to `width`.
      logical, intent(in), optional :: complete_obs !! Exclude rows missing in any column when true.
      logical, intent(in), optional :: na_restore !! Restore current-cell NaNs when true.
      logical, intent(in), optional :: online !! Accepted for API parity; deterministic batch evaluation is used.
      call unary_real_matrix(3, x, width, result, weights, min_obs, complete_obs, na_restore, online)
   end subroutine roll_mean_mat

   pure subroutine roll_min_mat(x, width, result, weights, min_obs, complete_obs, na_restore, online)
      real(dp), intent(in) :: x(:, :) !! Numeric matrix with observations in rows and variables in columns.
      integer, intent(in) :: width !! Rolling window width; must be at least one.
      real(dp), allocatable, intent(out) :: result(:, :) !! Rolling minima with the same shape as `x`.
      real(dp), intent(in), optional :: weights(:) !! Accepted for API parity; positive weights do not alter minima.
      integer, intent(in), optional :: min_obs !! Minimum usable observations; defaults to `width`.
      logical, intent(in), optional :: complete_obs !! Exclude rows missing in any column when true.
      logical, intent(in), optional :: na_restore !! Restore current-cell NaNs when true.
      logical, intent(in), optional :: online !! Accepted for API parity; deterministic batch evaluation is used.
      call unary_real_matrix(4, x, width, result, weights, min_obs, complete_obs, na_restore, online)
   end subroutine roll_min_mat

   pure subroutine roll_max_mat(x, width, result, weights, min_obs, complete_obs, na_restore, online)
      real(dp), intent(in) :: x(:, :) !! Numeric matrix with observations in rows and variables in columns.
      integer, intent(in) :: width !! Rolling window width; must be at least one.
      real(dp), allocatable, intent(out) :: result(:, :) !! Rolling maxima with the same shape as `x`.
      real(dp), intent(in), optional :: weights(:) !! Accepted for API parity; positive weights do not alter maxima.
      integer, intent(in), optional :: min_obs !! Minimum usable observations; defaults to `width`.
      logical, intent(in), optional :: complete_obs !! Exclude rows missing in any column when true.
      logical, intent(in), optional :: na_restore !! Restore current-cell NaNs when true.
      logical, intent(in), optional :: online !! Accepted for API parity; deterministic batch evaluation is used.
      call unary_real_matrix(5, x, width, result, weights, min_obs, complete_obs, na_restore, online)
   end subroutine roll_max_mat

   pure subroutine roll_median_mat(x, width, result, weights, min_obs, complete_obs, na_restore, online)
      real(dp), intent(in) :: x(:, :) !! Numeric matrix with observations in rows and variables in columns.
      integer, intent(in) :: width !! Rolling window width; must be at least one.
      real(dp), allocatable, intent(out) :: result(:, :) !! Rolling weighted medians with the same shape as `x`.
      real(dp), intent(in), optional :: weights(:) !! Positive observation weights aligned from the end.
      integer, intent(in), optional :: min_obs !! Minimum usable observations; defaults to `width`.
      logical, intent(in), optional :: complete_obs !! Exclude rows missing in any column when true.
      logical, intent(in), optional :: na_restore !! Restore current-cell NaNs when true.
      logical, intent(in), optional :: online !! Accepted for API parity; deterministic batch evaluation is used.
      call unary_real_matrix(6, x, width, result, weights, min_obs, complete_obs, na_restore, online)
   end subroutine roll_median_mat

   pure subroutine unary_real_matrix(code, x, width, result, weights, min_obs, complete_obs, na_restore, online)
      integer, intent(in) :: code !! Internal statistic selector used by public matrix wrappers.
      real(dp), intent(in) :: x(:, :) !! Numeric matrix with observations in rows.
      integer, intent(in) :: width !! Rolling window width.
      real(dp), allocatable, intent(out) :: result(:, :) !! Matrix result with the same shape as `x`.
      real(dp), intent(in), optional :: weights(:) !! Observation weights aligned from the end.
      integer, intent(in), optional :: min_obs !! Minimum usable observations.
      logical, intent(in), optional :: complete_obs !! Exclude rows missing in any column when true.
      logical, intent(in), optional :: na_restore !! Restore current-cell NaNs when true.
      logical, intent(in), optional :: online !! Accepted compatibility selector.
      integer :: j
      real(dp), allocatable :: column(:), temp(:)
      logical :: complete, restore

      complete = .false.
      if (present(complete_obs)) complete = complete_obs
      restore = .false.
      if (present(na_restore)) restore = na_restore
      allocate(result(size(x, 1), size(x, 2)))
      do j = 1, size(x, 2)
         call masked_real_column(x, j, complete, column)
         select case (code)
         case (1); call roll_sum_vec(column, width, temp, weights, min_obs, na_restore=.false., online=online)
         case (2); call roll_prod_vec(column, width, temp, weights, min_obs, na_restore=.false., online=online)
         case (3); call roll_mean_vec(column, width, temp, weights, min_obs, na_restore=.false., online=online)
         case (4); call roll_min_vec(column, width, temp, weights, min_obs, na_restore=.false., online=online)
         case (5); call roll_max_vec(column, width, temp, weights, min_obs, na_restore=.false., online=online)
         case (6); call roll_median_vec(column, width, temp, weights, min_obs, na_restore=.false., online=online)
         end select
         result(:, j) = temp
         if (restore) where (ieee_is_nan(x(:, j))) result(:, j) = na_real()
      end do
   end subroutine unary_real_matrix

   pure subroutine roll_idxmin_mat(x, width, result, weights, min_obs, complete_obs, na_restore, online)
      real(dp), intent(in) :: x(:, :) !! Numeric matrix with observations in rows and variables in columns.
      integer, intent(in) :: width !! Rolling window width; must be at least one.
      integer, allocatable, intent(out) :: result(:, :) !! Window-relative one-based indices of minima.
      real(dp), intent(in), optional :: weights(:) !! Accepted for R API parity; positive weights do not alter the index.
      integer, intent(in), optional :: min_obs !! Minimum usable observations; defaults to `width`.
      logical, intent(in), optional :: complete_obs !! Exclude rows missing in any column when true.
      logical, intent(in), optional :: na_restore !! Restore missing current cells when true.
      logical, intent(in), optional :: online !! Accepted for API parity; deterministic batch evaluation is used.
      call index_matrix(.true., x, width, result, weights, min_obs, complete_obs, na_restore, online)
   end subroutine roll_idxmin_mat

   pure subroutine roll_idxmax_mat(x, width, result, weights, min_obs, complete_obs, na_restore, online)
      real(dp), intent(in) :: x(:, :) !! Numeric matrix with observations in rows and variables in columns.
      integer, intent(in) :: width !! Rolling window width; must be at least one.
      integer, allocatable, intent(out) :: result(:, :) !! Window-relative one-based indices of maxima.
      real(dp), intent(in), optional :: weights(:) !! Accepted for R API parity; positive weights do not alter the index.
      integer, intent(in), optional :: min_obs !! Minimum usable observations; defaults to `width`.
      logical, intent(in), optional :: complete_obs !! Exclude rows missing in any column when true.
      logical, intent(in), optional :: na_restore !! Restore missing current cells when true.
      logical, intent(in), optional :: online !! Accepted for API parity; deterministic batch evaluation is used.
      call index_matrix(.false., x, width, result, weights, min_obs, complete_obs, na_restore, online)
   end subroutine roll_idxmax_mat

   pure subroutine index_matrix(want_min, x, width, result, weights, min_obs, complete_obs, na_restore, online)
      logical, intent(in) :: want_min !! Select minima when true and maxima when false.
      real(dp), intent(in) :: x(:, :) !! Numeric matrix with observations in rows.
      integer, intent(in) :: width !! Rolling window width.
      integer, allocatable, intent(out) :: result(:, :) !! Window-relative extrema indices.
      real(dp), intent(in), optional :: weights(:) !! Accepted weights argument for API compatibility.
      integer, intent(in), optional :: min_obs !! Minimum usable observations.
      logical, intent(in), optional :: complete_obs !! Exclude complete rows when true.
      logical, intent(in), optional :: na_restore !! Restore current-cell missing values when true.
      logical, intent(in), optional :: online !! Accepted compatibility selector.
      integer :: j
      real(dp), allocatable :: column(:)
      integer, allocatable :: temp(:)
      logical :: complete, restore

      complete = .false.
      if (present(complete_obs)) complete = complete_obs
      restore = .false.
      if (present(na_restore)) restore = na_restore
      allocate(result(size(x, 1), size(x, 2)))
      do j = 1, size(x, 2)
         call masked_real_column(x, j, complete, column)
         if (want_min) then
            call roll_idxmin_vec(column, width, temp, weights, min_obs, na_restore=.false., online=online)
         else
            call roll_idxmax_vec(column, width, temp, weights, min_obs, na_restore=.false., online=online)
         end if
         result(:, j) = temp
         if (restore) where (ieee_is_nan(x(:, j))) result(:, j) = roll_na_logical
      end do
   end subroutine index_matrix

   pure subroutine roll_quantile_mat(x, width, p, result, weights, min_obs, complete_obs, na_restore, online)
      real(dp), intent(in) :: x(:, :) !! Numeric matrix with observations in rows and variables in columns.
      integer, intent(in) :: width !! Rolling window width; must be at least one.
      real(dp), intent(in) :: p !! Quantile probability in the closed interval `[0,1]`.
      real(dp), allocatable, intent(out) :: result(:, :) !! Rolling weighted type-2 quantiles.
      real(dp), intent(in), optional :: weights(:) !! Positive observation weights aligned from the end.
      integer, intent(in), optional :: min_obs !! Minimum usable observations; defaults to `width`.
      logical, intent(in), optional :: complete_obs !! Exclude rows missing in any column when true.
      logical, intent(in), optional :: na_restore !! Restore current-cell NaNs when true.
      logical, intent(in), optional :: online !! Accepted for API parity; deterministic batch evaluation is used.
      integer :: j
      real(dp), allocatable :: column(:), temp(:)
      logical :: complete, restore

      complete = .false.
      if (present(complete_obs)) complete = complete_obs
      restore = .false.
      if (present(na_restore)) restore = na_restore
      allocate(result(size(x, 1), size(x, 2)))
      do j = 1, size(x, 2)
         call masked_real_column(x, j, complete, column)
         call roll_quantile_vec(column, width, p, temp, weights, min_obs, na_restore=.false., online=online)
         result(:, j) = temp
         if (restore) where (ieee_is_nan(x(:, j))) result(:, j) = na_real()
      end do
   end subroutine roll_quantile_mat

   pure subroutine roll_var_mat(x, width, result, weights, center, min_obs, complete_obs, na_restore, online)
      real(dp), intent(in) :: x(:, :) !! Numeric matrix with observations in rows and variables in columns.
      integer, intent(in) :: width !! Rolling window width; must be at least one.
      real(dp), allocatable, intent(out) :: result(:, :) !! Unbiased weighted rolling variances.
      real(dp), intent(in), optional :: weights(:) !! Observation weights aligned from the end.
      logical, intent(in), optional :: center !! Use weighted rolling means when true.
      integer, intent(in), optional :: min_obs !! Minimum usable observations; defaults to `width`.
      logical, intent(in), optional :: complete_obs !! Exclude rows missing in any column when true.
      logical, intent(in), optional :: na_restore !! Restore current-cell NaNs when true.
      logical, intent(in), optional :: online !! Accepted for API parity; deterministic batch evaluation is used.
      call moment_matrix(1, x, width, result, weights, center, .true., min_obs, complete_obs, na_restore, online)
   end subroutine roll_var_mat

   pure subroutine roll_sd_mat(x, width, result, weights, center, min_obs, complete_obs, na_restore, online)
      real(dp), intent(in) :: x(:, :) !! Numeric matrix with observations in rows and variables in columns.
      integer, intent(in) :: width !! Rolling window width; must be at least one.
      real(dp), allocatable, intent(out) :: result(:, :) !! Rolling standard deviations.
      real(dp), intent(in), optional :: weights(:) !! Observation weights aligned from the end.
      logical, intent(in), optional :: center !! Use weighted rolling means when true.
      integer, intent(in), optional :: min_obs !! Minimum usable observations; defaults to `width`.
      logical, intent(in), optional :: complete_obs !! Exclude rows missing in any column when true.
      logical, intent(in), optional :: na_restore !! Restore current-cell NaNs when true.
      logical, intent(in), optional :: online !! Accepted for API parity; deterministic batch evaluation is used.
      call moment_matrix(2, x, width, result, weights, center, .true., min_obs, complete_obs, na_restore, online)
   end subroutine roll_sd_mat

   pure subroutine roll_scale_mat(x, width, result, weights, center, scale, min_obs, complete_obs, na_restore, online)
      real(dp), intent(in) :: x(:, :) !! Numeric matrix with observations in rows and variables in columns.
      integer, intent(in) :: width !! Rolling window width; must be at least one.
      real(dp), allocatable, intent(out) :: result(:, :) !! Latest finite observations after rolling centering/scaling.
      real(dp), intent(in), optional :: weights(:) !! Observation weights aligned from the end.
      logical, intent(in), optional :: center !! Subtract weighted rolling means when true.
      logical, intent(in), optional :: scale !! Divide by weighted rolling standard deviations when true.
      integer, intent(in), optional :: min_obs !! Minimum usable observations; defaults to `width`.
      logical, intent(in), optional :: complete_obs !! Exclude rows missing in any column when true.
      logical, intent(in), optional :: na_restore !! Restore current-cell NaNs when true.
      logical, intent(in), optional :: online !! Accepted for API parity; deterministic batch evaluation is used.
      call moment_matrix(3, x, width, result, weights, center, scale, min_obs, complete_obs, na_restore, online)
   end subroutine roll_scale_mat

   pure subroutine moment_matrix(code, x, width, result, weights, center, scale, min_obs, complete_obs, na_restore, online)
      integer, intent(in) :: code !! Internal selector: one variance, two standard deviation, three scaling.
      real(dp), intent(in) :: x(:, :) !! Numeric matrix with observations in rows.
      integer, intent(in) :: width !! Rolling window width.
      real(dp), allocatable, intent(out) :: result(:, :) !! Matrix result with the same shape as `x`.
      real(dp), intent(in), optional :: weights(:) !! Observation weights aligned from the end.
      logical, intent(in), optional :: center !! Use weighted rolling means when true.
      logical, intent(in), optional :: scale !! Scale by rolling standard deviation when true.
      integer, intent(in), optional :: min_obs !! Minimum usable observations.
      logical, intent(in), optional :: complete_obs !! Exclude rows missing in any column when true.
      logical, intent(in), optional :: na_restore !! Restore current-cell NaNs when true.
      logical, intent(in), optional :: online !! Accepted compatibility selector.
      integer :: j
      real(dp), allocatable :: column(:), temp(:)
      logical :: complete, restore

      complete = .false.
      if (present(complete_obs)) complete = complete_obs
      restore = .false.
      if (present(na_restore)) restore = na_restore
      allocate(result(size(x, 1), size(x, 2)))
      do j = 1, size(x, 2)
         call masked_real_column(x, j, complete, column)
         select case (code)
         case (1); call roll_var_vec(column, width, temp, weights, center, min_obs, na_restore=.false., online=online)
         case (2); call roll_sd_vec(column, width, temp, weights, center, min_obs, na_restore=.false., online=online)
         case (3); call roll_scale_vec(column, width, temp, weights, center, scale, min_obs, na_restore=.false., online=online)
         end select
         result(:, j) = temp
         if (restore) where (ieee_is_nan(x(:, j))) result(:, j) = na_real()
      end do
   end subroutine moment_matrix

end module roll_matrix
