! SPDX-License-Identifier: Artistic-2.0
! Translation of matrixStats 1.5.0; attribution and material differences in NOTICE.md.
module matrixstats
   use matrixstats_common
   use matrixstats_order
   use r_ordering, only: r_sort_values_in_place
   use, intrinsic :: ieee_arithmetic, only: ieee_is_nan, ieee_is_finite
   implicit none
   private
   public :: dp, ms_na_index
   public :: sum2
   public :: mean2
   public :: all_value
   public :: any_value
   public :: count_value
   public :: log_sum_exp
   public :: iqr
   public :: product
   public :: any_missing
   public :: diff2
   public :: weighted_mean
   public :: weighted_median
   public :: weighted_var
   public :: weighted_sd
   public :: weighted_mad
   public :: row_sums2
   public :: row_means2
   public :: row_vars
   public :: row_sds
   public :: row_medians
   public :: row_mads
   public :: row_iqrs
   public :: row_mins
   public :: row_maxs
   public :: row_log_sum_exps
   public :: row_counts
   public :: row_alls
   public :: row_anys
   public :: row_prods
   public :: row_any_missings
   public :: row_any_nas
   public :: row_weighted_means
   public :: row_weighted_medians
   public :: row_weighted_vars
   public :: row_weighted_sds
   public :: row_weighted_mads
   public :: row_ranges
   public :: row_quantiles
   public :: row_order_stats
   public :: row_ranks
   public :: row_cumsums
   public :: row_cumprods
   public :: row_cummins
   public :: row_cummaxs
   public :: row_diffs
   public :: col_sums2
   public :: col_means2
   public :: col_vars
   public :: col_sds
   public :: col_medians
   public :: col_mads
   public :: col_iqrs
   public :: col_mins
   public :: col_maxs
   public :: col_log_sum_exps
   public :: col_counts
   public :: col_alls
   public :: col_anys
   public :: col_prods
   public :: col_any_missings
   public :: col_any_nas
   public :: col_weighted_means
   public :: col_weighted_medians
   public :: col_weighted_vars
   public :: col_weighted_sds
   public :: col_weighted_mads
   public :: col_ranges
   public :: col_quantiles
   public :: col_order_stats
   public :: col_ranks
   public :: col_cumsums
   public :: col_cumprods
   public :: col_cummins
   public :: col_cummaxs
   public :: col_diffs
   public :: var_diff
   public :: row_var_diffs
   public :: col_var_diffs
   public :: sd_diff
   public :: row_sd_diffs
   public :: col_sd_diffs
   public :: mad_diff
   public :: row_mad_diffs
   public :: col_mad_diffs
   public :: iqr_diff
   public :: row_iqr_diffs
   public :: col_iqr_diffs
   public :: bin_counts
   public :: bin_means
   public :: sign_tabulate
   public :: index_by_row
   public :: row_collapse
   public :: col_collapse
   public :: row_tabulates
   public :: col_tabulates
   public :: row_avgs_per_col_set
   public :: col_avgs_per_row_set
   public :: x_op_y
   public :: t_tx_op_y
contains

   pure function sum2(x, idxs, na_rm) result(y)
      real(dp), intent(in) :: x(:) !! Numeric vector; NaN represents both R NA and NaN.
      integer, intent(in), optional :: idxs(:) !! Optional R-style vector indices; negative excludes, zero drops.
      logical, intent(in), optional :: na_rm !! Remove NaN observations when true; default false.
      real(dp) :: y
      y = ms_sum(ms_vector(x, idxs), na_rm = na_rm)
   end function sum2

   pure function mean2(x, idxs, na_rm, refine) result(y)
      real(dp), intent(in) :: x(:) !! Numeric vector; NaN represents both R NA and NaN.
      integer, intent(in), optional :: idxs(:) !! Optional R-style vector indices; negative excludes, zero drops.
      logical, intent(in), optional :: na_rm !! Remove NaN observations when true; default false.
      logical, intent(in), optional :: refine !! Refine the mean using a residual correction; default true.
      real(dp) :: y
      y = ms_mean(ms_vector(x, idxs), na_rm = na_rm, refine = refine)
   end function mean2

   pure function all_value(x, idxs, na_rm, value) result(y)
      real(dp), intent(in) :: x(:) !! Numeric vector; NaN represents both R NA and NaN.
      integer, intent(in), optional :: idxs(:) !! Optional R-style vector indices; negative excludes, zero drops.
      logical, intent(in), optional :: na_rm !! Remove NaN observations when true; default false.
      real(dp), intent(in), optional :: value !! Explicit numeric target; omitted tests nonzero, and NaN matches missingness.
      real(dp) :: y
      y = ms_all(ms_vector(x, idxs), na_rm = na_rm, value = value)
   end function all_value

   pure function any_value(x, idxs, na_rm, value) result(y)
      real(dp), intent(in) :: x(:) !! Numeric vector; NaN represents both R NA and NaN.
      integer, intent(in), optional :: idxs(:) !! Optional R-style vector indices; negative excludes, zero drops.
      logical, intent(in), optional :: na_rm !! Remove NaN observations when true; default false.
      real(dp), intent(in), optional :: value !! Explicit numeric target; omitted tests nonzero, and NaN matches missingness.
      real(dp) :: y
      y = ms_any(ms_vector(x, idxs), na_rm = na_rm, value = value)
   end function any_value

   pure function count_value(x, idxs, na_rm, value) result(y)
      real(dp), intent(in) :: x(:) !! Numeric vector; NaN represents both R NA and NaN.
      integer, intent(in), optional :: idxs(:) !! Optional R-style vector indices; negative excludes, zero drops.
      logical, intent(in), optional :: na_rm !! Remove NaN observations when true; default false.
      real(dp), intent(in), optional :: value !! Exact target value; default one, and NaN matches missingness.
      real(dp) :: y
      y = ms_count(ms_vector(x, idxs), na_rm = na_rm, value = value)
   end function count_value

   pure function log_sum_exp(x, idxs, na_rm) result(y)
      real(dp), intent(in) :: x(:) !! Numeric vector; NaN represents both R NA and NaN.
      integer, intent(in), optional :: idxs(:) !! Optional R-style vector indices; negative excludes, zero drops.
      logical, intent(in), optional :: na_rm !! Remove NaN observations when true; default false.
      real(dp) :: y
      y = ms_log_sum_exp(ms_vector(x, idxs), na_rm = na_rm)
   end function log_sum_exp

   pure function iqr(x, idxs, na_rm, type) result(y)
      real(dp), intent(in) :: x(:) !! Numeric vector; NaN represents both R NA and NaN.
      integer, intent(in), optional :: idxs(:) !! Optional R-style vector indices; negative excludes, zero drops.
      logical, intent(in), optional :: na_rm !! Remove NaN observations when true; default false.
      integer, intent(in), optional :: type !! Hyndman-Fan quantile type in 1 through 9; default seven.
      real(dp) :: y
      y = ms_iqr(ms_vector(x, idxs), na_rm = na_rm, type = type)
   end function iqr

   pure function product(x, idxs, na_rm) result(y)
      real(dp), intent(in) :: x(:) !! Numeric vector; NaN represents both R NA and NaN.
      integer, intent(in), optional :: idxs(:) !! Optional R-style vector indices; negative excludes, zero drops.
      logical, intent(in), optional :: na_rm !! Remove NaN observations when true; default false.
      real(dp) :: y
      y = ms_product(ms_vector(x, idxs), na_rm, 'expSumLog')
   end function product

   pure function any_missing(x, idxs) result(y)
      real(dp), intent(in) :: x(:) !! Numeric vector; NaN represents both R NA and NaN.
      integer, intent(in), optional :: idxs(:) !! Optional R-style vector indices; negative excludes, zero drops.
      logical :: y
      y = any(ieee_is_nan(ms_vector(x, idxs)))
   end function any_missing

   pure function diff2(x, idxs, lag, differences) result(y)
      real(dp), intent(in) :: x(:) !! Numeric vector; NaN represents both R NA and NaN.
      integer, intent(in), optional :: idxs(:) !! Optional R-style vector indices; negative excludes, zero drops.
      integer, intent(in), optional :: lag !! Positive difference lag in observations; default one.
      integer, intent(in), optional :: differences !! Nonnegative difference order; default one.
      real(dp), allocatable :: y(:)
      y = ms_difference(ms_vector(x, idxs), lag, differences)
   end function diff2

   pure function weighted_mean(x, w, idxs, na_rm, refine) result(y)
      real(dp), intent(in) :: x(:) !! Numeric vector; NaN represents both R NA and NaN.
      real(dp), intent(in), optional :: w(:) !! Weights, same length as unselected observations; default equal weights.
      integer, intent(in), optional :: idxs(:) !! Optional R-style vector indices; negative excludes, zero drops.
      logical, intent(in), optional :: na_rm !! Remove NaN observations when true; default false.
      logical, intent(in), optional :: refine !! Refine the mean using a residual correction; default false.
      real(dp) :: y
      real(dp), allocatable :: a(:), weights(:)
      a = ms_vector(x, idxs)
      if (present(w)) then
         weights = ms_vector(w, idxs)
         y = ms_weighted_mean(a, weights, na_rm, refine = refine)
      else
         y = ms_weighted_mean(a, na_rm = na_rm, refine = refine)
      end if
   end function weighted_mean

   pure function weighted_median(x, w, idxs, na_rm, interpolate, ties) result(y)
      real(dp), intent(in) :: x(:) !! Numeric vector; NaN represents both R NA and NaN.
      real(dp), intent(in), optional :: w(:) !! Weights, same length as unselected observations; default equal weights.
      integer, intent(in), optional :: idxs(:) !! Optional R-style vector indices; negative excludes, zero drops.
      logical, intent(in), optional :: na_rm !! Remove NaN observations when true; default false.
      logical, intent(in), optional :: interpolate !! Interpolate between observations; default true unless ties is supplied.
      character(*), intent(in), optional :: ties !! Half-weight tie rule: weighted (default), min, max, or mean.
      real(dp) :: y
      real(dp), allocatable :: a(:), weights(:)
      a = ms_vector(x, idxs)
      if (present(w)) then
         weights = ms_vector(w, idxs)
         y = ms_weighted_median(a, weights, na_rm, interpolate = interpolate, ties = ties)
      else
         y = ms_weighted_median(a, na_rm = na_rm, interpolate = interpolate, ties = ties)
      end if
   end function weighted_median

   pure function weighted_var(x, w, idxs, na_rm, center) result(y)
      real(dp), intent(in) :: x(:) !! Numeric vector; NaN represents both R NA and NaN.
      real(dp), intent(in), optional :: w(:) !! Weights, same length as unselected observations; default equal weights.
      integer, intent(in), optional :: idxs(:) !! Optional R-style vector indices; negative excludes, zero drops.
      logical, intent(in), optional :: na_rm !! Remove NaN observations when true; default false.
      real(dp), intent(in), optional :: center !! Prespecified weighted mean; omitted means estimate internally.
      real(dp) :: y
      real(dp), allocatable :: a(:), weights(:)
      a = ms_vector(x, idxs)
      if (present(w)) then
         weights = ms_vector(w, idxs)
         y = ms_weighted_var(a, weights, na_rm, center = center)
      else
         y = ms_weighted_var(a, na_rm = na_rm, center = center)
      end if
   end function weighted_var

   pure function weighted_sd(x, w, idxs, na_rm, center) result(y)
      real(dp), intent(in) :: x(:) !! Numeric vector; NaN represents both R NA and NaN.
      real(dp), intent(in), optional :: w(:) !! Weights, same length as unselected observations; default equal weights.
      integer, intent(in), optional :: idxs(:) !! Optional R-style vector indices; negative excludes, zero drops.
      logical, intent(in), optional :: na_rm !! Remove NaN observations when true; default false.
      real(dp), intent(in), optional :: center !! Prespecified weighted mean; omitted means estimate internally.
      real(dp) :: y
      real(dp), allocatable :: a(:), weights(:)
      a = ms_vector(x, idxs)
      if (present(w)) then
         weights = ms_vector(w, idxs)
         y = ms_weighted_sd(a, weights, na_rm, center = center)
      else
         y = ms_weighted_sd(a, na_rm = na_rm, center = center)
      end if
   end function weighted_sd

   pure function weighted_mad(x, w, idxs, na_rm, constant, center) result(y)
      real(dp), intent(in) :: x(:) !! Numeric vector; NaN represents both R NA and NaN.
      real(dp), intent(in), optional :: w(:) !! Weights, same length as unselected observations; default equal weights.
      integer, intent(in), optional :: idxs(:) !! Optional R-style vector indices; negative excludes, zero drops.
      logical, intent(in), optional :: na_rm !! Remove NaN observations when true; default false.
      real(dp), intent(in), optional :: constant !! MAD scale multiplier; default 1.4826.
      real(dp), intent(in), optional :: center !! Prespecified weighted median; omitted means estimate internally.
      real(dp) :: y
      real(dp), allocatable :: a(:), weights(:)
      a = ms_vector(x, idxs)
      if (present(w)) then
         weights = ms_vector(w, idxs)
         y = ms_weighted_mad(a, weights, na_rm, constant = constant, center = center)
      else
         y = ms_weighted_mad(a, na_rm = na_rm, constant = constant, center = center)
      end if
   end function weighted_mad

   pure function row_sums2(x, rows, cols, na_rm) result(y)
      real(dp), intent(in) :: x(:, :) !! Numeric matrix; NaN represents both R NA and NaN.
      integer, intent(in), optional :: rows(:) !! Optional R-style row indices; negative excludes, zero drops.
      integer, intent(in), optional :: cols(:) !! Optional R-style column indices; negative excludes, zero drops.
      logical, intent(in), optional :: na_rm !! Remove NaN observations when true; default false.
      real(dp), allocatable :: y(:)
      real(dp), allocatable :: a(:, :), centers(:)
      integer :: j
      a = ms_select(x, rows, cols)
      allocate(y(size(a, 1)))
      do j = 1, size(a, 1)
         y(j) = ms_sum(a(j, :), na_rm = na_rm)
      end do
   end function row_sums2

   pure function row_means2(x, rows, cols, na_rm, refine) result(y)
      real(dp), intent(in) :: x(:, :) !! Numeric matrix; NaN represents both R NA and NaN.
      integer, intent(in), optional :: rows(:) !! Optional R-style row indices; negative excludes, zero drops.
      integer, intent(in), optional :: cols(:) !! Optional R-style column indices; negative excludes, zero drops.
      logical, intent(in), optional :: na_rm !! Remove NaN observations when true; default false.
      logical, intent(in), optional :: refine !! Refine the mean using a residual correction; default true.
      real(dp), allocatable :: y(:)
      real(dp), allocatable :: a(:, :), centers(:)
      integer :: j
      a = ms_select(x, rows, cols)
      allocate(y(size(a, 1)))
      do j = 1, size(a, 1)
         y(j) = ms_mean(a(j, :), na_rm = na_rm, refine = refine)
      end do
   end function row_means2

   pure function row_vars(x, rows, cols, na_rm, refine, center) result(y)
      real(dp), intent(in) :: x(:, :) !! Numeric matrix; NaN represents both R NA and NaN.
      integer, intent(in), optional :: rows(:) !! Optional R-style row indices; negative excludes, zero drops.
      integer, intent(in), optional :: cols(:) !! Optional R-style column indices; negative excludes, zero drops.
      logical, intent(in), optional :: na_rm !! Remove NaN observations when true; default false.
      logical, intent(in), optional :: refine !! Refine the mean using a residual correction; default true.
      real(dp), intent(in), optional :: center(:) !! Means for each original row or column, before selection.
      real(dp), allocatable :: y(:)
      real(dp), allocatable :: a(:, :), centers(:)
      integer :: j
      a = ms_select(x, rows, cols)
      allocate(y(size(a, 1)))
      if (present(center)) then
         centers = ms_vector(center, rows)
         if (size(centers) /= size(a, 1)) then
            y = ms_nan()
            return
         end if
      end if
      do j = 1, size(a, 1)
         if (present(center)) then
            y(j) = ms_var(a(j, :), na_rm = na_rm, refine = refine, center = centers(j))
         else
            y(j) = ms_var(a(j, :), na_rm = na_rm, refine = refine)
         end if
      end do
   end function row_vars

   pure function row_sds(x, rows, cols, na_rm, refine, center) result(y)
      real(dp), intent(in) :: x(:, :) !! Numeric matrix; NaN represents both R NA and NaN.
      integer, intent(in), optional :: rows(:) !! Optional R-style row indices; negative excludes, zero drops.
      integer, intent(in), optional :: cols(:) !! Optional R-style column indices; negative excludes, zero drops.
      logical, intent(in), optional :: na_rm !! Remove NaN observations when true; default false.
      logical, intent(in), optional :: refine !! Refine the mean using a residual correction; default true.
      real(dp), intent(in), optional :: center(:) !! Means for each original row or column, before selection.
      real(dp), allocatable :: y(:)
      real(dp), allocatable :: a(:, :), centers(:)
      integer :: j
      a = ms_select(x, rows, cols)
      allocate(y(size(a, 1)))
      if (present(center)) then
         centers = ms_vector(center, rows)
         if (size(centers) /= size(a, 1)) then
            y = ms_nan()
            return
         end if
      end if
      do j = 1, size(a, 1)
         if (present(center)) then
            y(j) = ms_var(a(j, :), na_rm = na_rm, refine = refine, center = centers(j))
         else
            y(j) = ms_var(a(j, :), na_rm = na_rm, refine = refine)
         end if
         y(j) = sqrt(y(j))
      end do
   end function row_sds

   pure function row_medians(x, rows, cols, na_rm) result(y)
      real(dp), intent(in) :: x(:, :) !! Numeric matrix; NaN represents both R NA and NaN.
      integer, intent(in), optional :: rows(:) !! Optional R-style row indices; negative excludes, zero drops.
      integer, intent(in), optional :: cols(:) !! Optional R-style column indices; negative excludes, zero drops.
      logical, intent(in), optional :: na_rm !! Remove NaN observations when true; default false.
      real(dp), allocatable :: y(:)
      real(dp), allocatable :: a(:, :), centers(:)
      integer :: j
      a = ms_select(x, rows, cols)
      allocate(y(size(a, 1)))
      do j = 1, size(a, 1)
         y(j) = ms_median(a(j, :), na_rm = na_rm)
      end do
   end function row_medians

   pure function row_mads(x, rows, cols, na_rm, constant, center) result(y)
      real(dp), intent(in) :: x(:, :) !! Numeric matrix; NaN represents both R NA and NaN.
      integer, intent(in), optional :: rows(:) !! Optional R-style row indices; negative excludes, zero drops.
      integer, intent(in), optional :: cols(:) !! Optional R-style column indices; negative excludes, zero drops.
      logical, intent(in), optional :: na_rm !! Remove NaN observations when true; default false.
      real(dp), intent(in), optional :: constant !! MAD scale multiplier; default 1.4826.
      real(dp), intent(in), optional :: center(:) !! Medians for each original row or column, before selection.
      real(dp), allocatable :: y(:)
      real(dp), allocatable :: a(:, :), centers(:)
      integer :: j
      a = ms_select(x, rows, cols)
      allocate(y(size(a, 1)))
      if (present(center)) then
         centers = ms_vector(center, rows)
         if (size(centers) /= size(a, 1)) then
            y = ms_nan()
            return
         end if
      end if
      do j = 1, size(a, 1)
         if (present(center)) then
            y(j) = ms_mad(a(j, :), na_rm = na_rm, constant = constant, center = centers(j))
         else
            y(j) = ms_mad(a(j, :), na_rm = na_rm, constant = constant)
         end if
      end do
   end function row_mads

   pure function row_iqrs(x, rows, cols, na_rm, type) result(y)
      real(dp), intent(in) :: x(:, :) !! Numeric matrix; NaN represents both R NA and NaN.
      integer, intent(in), optional :: rows(:) !! Optional R-style row indices; negative excludes, zero drops.
      integer, intent(in), optional :: cols(:) !! Optional R-style column indices; negative excludes, zero drops.
      logical, intent(in), optional :: na_rm !! Remove NaN observations when true; default false.
      integer, intent(in), optional :: type !! Hyndman-Fan quantile type in 1 through 9; default seven.
      real(dp), allocatable :: y(:)
      real(dp), allocatable :: a(:, :), centers(:)
      integer :: j
      a = ms_select(x, rows, cols)
      allocate(y(size(a, 1)))
      do j = 1, size(a, 1)
         y(j) = ms_quantile(a(j, :), 0.75_dp, na_rm, type) - ms_quantile(a(j, :), 0.25_dp, na_rm, type)
      end do
   end function row_iqrs

   pure function row_mins(x, rows, cols, na_rm) result(y)
      real(dp), intent(in) :: x(:, :) !! Numeric matrix; NaN represents both R NA and NaN.
      integer, intent(in), optional :: rows(:) !! Optional R-style row indices; negative excludes, zero drops.
      integer, intent(in), optional :: cols(:) !! Optional R-style column indices; negative excludes, zero drops.
      logical, intent(in), optional :: na_rm !! Remove NaN observations when true; default false.
      real(dp), allocatable :: y(:)
      real(dp), allocatable :: a(:, :), centers(:)
      integer :: j
      a = ms_select(x, rows, cols)
      allocate(y(size(a, 1)))
      do j = 1, size(a, 1)
         y(j) = ms_extreme(a(j, :), na_rm = na_rm, maximum = .false.)
      end do
   end function row_mins

   pure function row_maxs(x, rows, cols, na_rm) result(y)
      real(dp), intent(in) :: x(:, :) !! Numeric matrix; NaN represents both R NA and NaN.
      integer, intent(in), optional :: rows(:) !! Optional R-style row indices; negative excludes, zero drops.
      integer, intent(in), optional :: cols(:) !! Optional R-style column indices; negative excludes, zero drops.
      logical, intent(in), optional :: na_rm !! Remove NaN observations when true; default false.
      real(dp), allocatable :: y(:)
      real(dp), allocatable :: a(:, :), centers(:)
      integer :: j
      a = ms_select(x, rows, cols)
      allocate(y(size(a, 1)))
      do j = 1, size(a, 1)
         y(j) = ms_extreme(a(j, :), na_rm = na_rm, maximum = .true.)
      end do
   end function row_maxs

   pure function row_log_sum_exps(x, rows, cols, na_rm) result(y)
      real(dp), intent(in) :: x(:, :) !! Numeric matrix; NaN represents both R NA and NaN.
      integer, intent(in), optional :: rows(:) !! Optional R-style row indices; negative excludes, zero drops.
      integer, intent(in), optional :: cols(:) !! Optional R-style column indices; negative excludes, zero drops.
      logical, intent(in), optional :: na_rm !! Remove NaN observations when true; default false.
      real(dp), allocatable :: y(:)
      real(dp), allocatable :: a(:, :), centers(:)
      integer :: j
      a = ms_select(x, rows, cols)
      allocate(y(size(a, 1)))
      do j = 1, size(a, 1)
         y(j) = ms_log_sum_exp(a(j, :), na_rm = na_rm)
      end do
   end function row_log_sum_exps

   pure function row_counts(x, rows, cols, na_rm, value) result(y)
      real(dp), intent(in) :: x(:, :) !! Numeric matrix; NaN represents both R NA and NaN.
      integer, intent(in), optional :: rows(:) !! Optional R-style row indices; negative excludes, zero drops.
      integer, intent(in), optional :: cols(:) !! Optional R-style column indices; negative excludes, zero drops.
      logical, intent(in), optional :: na_rm !! Remove NaN observations when true; default false.
      real(dp), intent(in), optional :: value !! Exact target value; default one, and NaN matches missingness.
      real(dp), allocatable :: y(:)
      real(dp), allocatable :: a(:, :), centers(:)
      integer :: j
      a = ms_select(x, rows, cols)
      allocate(y(size(a, 1)))
      do j = 1, size(a, 1)
         y(j) = ms_count(a(j, :), na_rm = na_rm, value = value)
      end do
   end function row_counts

   pure function row_alls(x, rows, cols, na_rm, value) result(y)
      real(dp), intent(in) :: x(:, :) !! Numeric matrix; NaN represents both R NA and NaN.
      integer, intent(in), optional :: rows(:) !! Optional R-style row indices; negative excludes, zero drops.
      integer, intent(in), optional :: cols(:) !! Optional R-style column indices; negative excludes, zero drops.
      logical, intent(in), optional :: na_rm !! Remove NaN observations when true; default false.
      real(dp), intent(in), optional :: value !! Explicit numeric target; omitted tests nonzero, and NaN matches missingness.
      real(dp), allocatable :: y(:)
      real(dp), allocatable :: a(:, :), centers(:)
      integer :: j
      a = ms_select(x, rows, cols)
      allocate(y(size(a, 1)))
      do j = 1, size(a, 1)
         y(j) = ms_all(a(j, :), na_rm = na_rm, value = value)
      end do
   end function row_alls

   pure function row_anys(x, rows, cols, na_rm, value) result(y)
      real(dp), intent(in) :: x(:, :) !! Numeric matrix; NaN represents both R NA and NaN.
      integer, intent(in), optional :: rows(:) !! Optional R-style row indices; negative excludes, zero drops.
      integer, intent(in), optional :: cols(:) !! Optional R-style column indices; negative excludes, zero drops.
      logical, intent(in), optional :: na_rm !! Remove NaN observations when true; default false.
      real(dp), intent(in), optional :: value !! Explicit numeric target; omitted tests nonzero, and NaN matches missingness.
      real(dp), allocatable :: y(:)
      real(dp), allocatable :: a(:, :), centers(:)
      integer :: j
      a = ms_select(x, rows, cols)
      allocate(y(size(a, 1)))
      do j = 1, size(a, 1)
         y(j) = ms_any(a(j, :), na_rm = na_rm, value = value)
      end do
   end function row_anys

   pure function row_prods(x, rows, cols, na_rm, method) result(y)
      real(dp), intent(in) :: x(:, :) !! Numeric matrix; NaN represents both R NA and NaN.
      integer, intent(in), optional :: rows(:) !! Optional R-style row indices; negative excludes, zero drops.
      integer, intent(in), optional :: cols(:) !! Optional R-style column indices; negative excludes, zero drops.
      logical, intent(in), optional :: na_rm !! Remove NaN observations when true; default false.
      character(*), intent(in), optional :: method !! 'direct' product by default, or 'expSumLog'.
      real(dp), allocatable :: y(:)
      real(dp), allocatable :: a(:, :), centers(:)
      integer :: j
      a = ms_select(x, rows, cols)
      allocate(y(size(a, 1)))
      do j = 1, size(a, 1)
         y(j) = ms_product(a(j, :), na_rm = na_rm, method = method)
      end do
   end function row_prods

   pure function row_any_missings(x, rows, cols) result(y)
      real(dp), intent(in) :: x(:, :) !! Numeric matrix; NaN represents both R NA and NaN.
      integer, intent(in), optional :: rows(:) !! Optional R-style row indices; negative excludes, zero drops.
      integer, intent(in), optional :: cols(:) !! Optional R-style column indices; negative excludes, zero drops.
      logical, allocatable :: y(:)
      real(dp), allocatable :: a(:, :)
      integer :: j
      a = ms_select(x, rows, cols)
      allocate(y(size(a, 1)))
      do j = 1, size(a, 1)
         y(j) = any(ieee_is_nan(a(j, :)))
      end do
   end function row_any_missings

   pure function row_any_nas(x, rows, cols) result(y)
      real(dp), intent(in) :: x(:, :) !! Numeric matrix; NaN represents both R NA and NaN.
      integer, intent(in), optional :: rows(:) !! Optional R-style row indices; negative excludes, zero drops.
      integer, intent(in), optional :: cols(:) !! Optional R-style column indices; negative excludes, zero drops.
      logical, allocatable :: y(:)
      real(dp), allocatable :: a(:, :)
      integer :: j
      a = ms_select(x, rows, cols)
      allocate(y(size(a, 1)))
      do j = 1, size(a, 1)
         y(j) = any(ieee_is_nan(a(j, :)))
      end do
   end function row_any_nas

   pure function row_weighted_means(x, w, rows, cols, na_rm) result(y)
      real(dp), intent(in) :: x(:, :) !! Numeric matrix; NaN represents both R NA and NaN.
      real(dp), intent(in), optional :: w(:) !! Weights, same length as unselected observations; default equal weights.
      integer, intent(in), optional :: rows(:) !! Optional R-style row indices; negative excludes, zero drops.
      integer, intent(in), optional :: cols(:) !! Optional R-style column indices; negative excludes, zero drops.
      logical, intent(in), optional :: na_rm !! Remove NaN observations when true; default false.
      real(dp), allocatable :: y(:)
      real(dp), allocatable :: a(:, :), weights(:), centers(:)
      integer :: j
      a = ms_select(x, rows, cols)
      allocate(y(size(a, 1)))
      if (present(w)) then
         weights = ms_vector(w, cols)
      else
         allocate(weights(size(a, 2)))
         weights = 1.0_dp
      end if
      if (size(weights) /= size(a, 2)) then
         y = ms_nan()
         return
      end if
      if (any(weights < 0.0_dp)) then
         y = ms_nan()
         return
      end if
      do j = 1, size(a, 1)
         if (.not. present(w)) then
            y(j) = ms_mean(a(j, :), na_rm)
            cycle
         end if
         y(j) = matrix_weighted_mean(a(j, :), weights, na_rm)
      end do
   end function row_weighted_means

   pure function row_weighted_medians(x, w, rows, cols, na_rm, interpolate, ties) result(y)
      real(dp), intent(in) :: x(:, :) !! Numeric matrix; NaN represents both R NA and NaN.
      real(dp), intent(in), optional :: w(:) !! Weights, same length as unselected observations; default equal weights.
      integer, intent(in), optional :: rows(:) !! Optional R-style row indices; negative excludes, zero drops.
      integer, intent(in), optional :: cols(:) !! Optional R-style column indices; negative excludes, zero drops.
      logical, intent(in), optional :: na_rm !! Remove NaN observations when true; default false.
      logical, intent(in), optional :: interpolate !! Interpolate between observations; default true unless ties is supplied.
      character(*), intent(in), optional :: ties !! Half-weight median tie rule: weighted, min, max, or mean.
      real(dp), allocatable :: y(:)
      real(dp), allocatable :: a(:, :), weights(:), centers(:)
      integer :: j
      a = ms_select(x, rows, cols)
      allocate(y(size(a, 1)))
      if (present(w)) then
         weights = ms_vector(w, cols)
      else
         allocate(weights(size(a, 2)))
         weights = 1.0_dp
      end if
      if (size(weights) /= size(a, 2)) then
         y = ms_nan()
         return
      end if
      if (any(weights < 0.0_dp)) then
         y = ms_nan()
         return
      end if
      do j = 1, size(a, 1)
         if (.not. present(w)) then
            y(j) = ms_median(a(j, :), na_rm)
            cycle
         end if
         y(j) = ms_weighted_median(a(j, :), weights, na_rm, interpolate = interpolate, ties = ties)
      end do
   end function row_weighted_medians

   pure function row_weighted_vars(x, w, rows, cols, na_rm, center) result(y)
      real(dp), intent(in) :: x(:, :) !! Numeric matrix; NaN represents both R NA and NaN.
      real(dp), intent(in), optional :: w(:) !! Weights, same length as unselected observations; default equal weights.
      integer, intent(in), optional :: rows(:) !! Optional R-style row indices; negative excludes, zero drops.
      integer, intent(in), optional :: cols(:) !! Optional R-style column indices; negative excludes, zero drops.
      logical, intent(in), optional :: na_rm !! Remove NaN observations when true; default false.
      real(dp), intent(in), optional :: center(:) !! Weighted means for each original row or column.
      real(dp), allocatable :: y(:)
      real(dp), allocatable :: a(:, :), weights(:), centers(:)
      integer :: j
      a = ms_select(x, rows, cols)
      allocate(y(size(a, 1)))
      if (present(w)) then
         weights = ms_vector(w, cols)
      else
         allocate(weights(size(a, 2)))
         weights = 1.0_dp
      end if
      if (size(weights) /= size(a, 2)) then
         y = ms_nan()
         return
      end if
      if (present(center)) then
         centers = ms_vector(center, rows)
         if (size(centers) /= size(a, 1)) then
            y = ms_nan()
            return
         end if
      end if
      do j = 1, size(a, 1)
         if (present(center)) then
            y(j) = ms_weighted_var(a(j, :), weights, na_rm, center = centers(j))
         else
            y(j) = ms_weighted_var(a(j, :), weights, na_rm)
         end if
      end do
   end function row_weighted_vars

   pure function row_weighted_sds(x, w, rows, cols, na_rm, center) result(y)
      real(dp), intent(in) :: x(:, :) !! Numeric matrix; NaN represents both R NA and NaN.
      real(dp), intent(in), optional :: w(:) !! Weights, same length as unselected observations; default equal weights.
      integer, intent(in), optional :: rows(:) !! Optional R-style row indices; negative excludes, zero drops.
      integer, intent(in), optional :: cols(:) !! Optional R-style column indices; negative excludes, zero drops.
      logical, intent(in), optional :: na_rm !! Remove NaN observations when true; default false.
      real(dp), intent(in), optional :: center(:) !! Weighted means for each original row or column.
      real(dp), allocatable :: y(:)
      real(dp), allocatable :: a(:, :), weights(:), centers(:)
      integer :: j
      a = ms_select(x, rows, cols)
      allocate(y(size(a, 1)))
      if (present(w)) then
         weights = ms_vector(w, cols)
      else
         allocate(weights(size(a, 2)))
         weights = 1.0_dp
      end if
      if (size(weights) /= size(a, 2)) then
         y = ms_nan()
         return
      end if
      if (present(center)) then
         centers = ms_vector(center, rows)
         if (size(centers) /= size(a, 1)) then
            y = ms_nan()
            return
         end if
      end if
      do j = 1, size(a, 1)
         if (present(center)) then
            y(j) = ms_weighted_sd(a(j, :), weights, na_rm, center = centers(j))
         else
            y(j) = ms_weighted_sd(a(j, :), weights, na_rm)
         end if
      end do
   end function row_weighted_sds

   pure function row_weighted_mads(x, w, rows, cols, na_rm, constant, center) result(y)
      real(dp), intent(in) :: x(:, :) !! Numeric matrix; NaN represents both R NA and NaN.
      real(dp), intent(in), optional :: w(:) !! Weights, same length as unselected observations; default equal weights.
      integer, intent(in), optional :: rows(:) !! Optional R-style row indices; negative excludes, zero drops.
      integer, intent(in), optional :: cols(:) !! Optional R-style column indices; negative excludes, zero drops.
      logical, intent(in), optional :: na_rm !! Remove NaN observations when true; default false.
      real(dp), intent(in), optional :: constant !! MAD scale multiplier; default 1.4826.
      real(dp), intent(in), optional :: center(:) !! Weighted medians for each original row or column.
      real(dp), allocatable :: y(:)
      real(dp), allocatable :: a(:, :), weights(:), centers(:)
      integer :: j
      a = ms_select(x, rows, cols)
      allocate(y(size(a, 1)))
      if (present(w)) then
         weights = ms_vector(w, cols)
      else
         allocate(weights(size(a, 2)))
         weights = 1.0_dp
      end if
      if (size(weights) /= size(a, 2)) then
         y = ms_nan()
         return
      end if
      if (present(center)) then
         centers = ms_vector(center, rows)
         if (size(centers) /= size(a, 1)) then
            y = ms_nan()
            return
         end if
      end if
      do j = 1, size(a, 1)
         if (.not. present(w)) then
            if (present(center)) then
               y(j) = ms_weighted_mad(a(j, :), na_rm = na_rm, constant = constant, center = centers(j))
            else
               y(j) = ms_weighted_mad(a(j, :), na_rm = na_rm, constant = constant)
            end if
            cycle
         end if
         if (present(center)) then
            y(j) = ms_weighted_mad(a(j, :), weights, na_rm, constant = constant, center = centers(j))
         else
            y(j) = ms_weighted_mad(a(j, :), weights, na_rm, constant = constant)
         end if
      end do
   end function row_weighted_mads

   pure function row_ranges(x, rows, cols, na_rm) result(y)
      real(dp), intent(in) :: x(:, :) !! Numeric matrix; NaN represents both R NA and NaN.
      integer, intent(in), optional :: rows(:) !! Optional R-style row indices; negative excludes, zero drops.
      integer, intent(in), optional :: cols(:) !! Optional R-style column indices; negative excludes, zero drops.
      logical, intent(in), optional :: na_rm !! Remove NaN observations when true; default false.
      real(dp), allocatable :: y(:, :)
      real(dp), allocatable :: a(:, :)
      integer :: j
      a = ms_select(x, rows, cols)
      allocate(y(size(a, 1), 2))
      do j = 1, size(a, 1)
         y(j, 1) = ms_extreme(a(j, :), .false., na_rm)
         y(j, 2) = ms_extreme(a(j, :), .true., na_rm)
      end do
   end function row_ranges

   pure function row_quantiles(x, probs, rows, cols, na_rm, type) result(y)
      real(dp), intent(in) :: x(:, :) !! Numeric matrix; NaN represents both R NA and NaN.
      real(dp), intent(in) :: probs(:) !! Quantile probabilities, each in the closed interval zero to one.
      integer, intent(in), optional :: rows(:) !! Optional R-style row indices; negative excludes, zero drops.
      integer, intent(in), optional :: cols(:) !! Optional R-style column indices; negative excludes, zero drops.
      logical, intent(in), optional :: na_rm !! Remove NaN observations when true; default false.
      integer, intent(in), optional :: type !! Hyndman-Fan quantile type in 1 through 9; default seven.
      real(dp), allocatable :: y(:, :)
      real(dp), allocatable :: a(:, :)
      integer :: j, k
      a = ms_select(x, rows, cols)
      allocate(y(size(a, 1), size(probs)))
      do k = 1, size(probs)
         do j = 1, size(a, 1)
            y(j, k) = ms_quantile(a(j, :), probs(k), na_rm, type)
         end do
      end do
   end function row_quantiles

   pure function row_order_stats(x, which, rows, cols) result(y)
      real(dp), intent(in) :: x(:, :) !! Numeric matrix; NaN represents both R NA and NaN.
      integer, intent(in) :: which !! One-based order statistic rank within the selected observations.
      integer, intent(in), optional :: rows(:) !! Optional R-style row indices; negative excludes, zero drops.
      integer, intent(in), optional :: cols(:) !! Optional R-style column indices; negative excludes, zero drops.
      real(dp), allocatable :: y(:)
      real(dp), allocatable :: a(:, :)
      integer :: j
      a = ms_select(x, rows, cols)
      allocate(y(size(a, 1)))
      do j = 1, size(a, 1)
         y(j) = ms_order_stat(a(j, :), which)
      end do
   end function row_order_stats

   pure function row_ranks(x, rows, cols, ties_method) result(y)
      real(dp), intent(in) :: x(:, :) !! Numeric matrix; NaN represents both R NA and NaN.
      integer, intent(in), optional :: rows(:) !! Optional R-style row indices; negative excludes, zero drops.
      integer, intent(in), optional :: cols(:) !! Optional R-style column indices; negative excludes, zero drops.
      character(*), intent(in), optional :: ties_method !! Rank ties: max (default), min, average, first, last, or dense.
      real(dp), allocatable :: y(:, :)
      real(dp), allocatable :: a(:, :)
      integer :: j
      a = ms_select(x, rows, cols)
      allocate(y(size(a, 1), size(a, 2)))
      do j = 1, size(a, 1)
         y(j, :) = ms_ranks(a(j, :), ties_method)
      end do
   end function row_ranks

   pure function row_cumsums(x, rows, cols) result(y)
      real(dp), intent(in) :: x(:, :) !! Numeric matrix; NaN represents both R NA and NaN.
      integer, intent(in), optional :: rows(:) !! Optional R-style row indices; negative excludes, zero drops.
      integer, intent(in), optional :: cols(:) !! Optional R-style column indices; negative excludes, zero drops.
      real(dp), allocatable :: y(:, :)
      real(dp), allocatable :: a(:, :)
      integer :: j
      a = ms_select(x, rows, cols)
      allocate(y(size(a, 1), size(a, 2)))
      do j = 1, size(a, 1)
         y(j, :) = ms_cumulative(a(j, :), 'sum')
      end do
   end function row_cumsums

   pure function row_cumprods(x, rows, cols) result(y)
      real(dp), intent(in) :: x(:, :) !! Numeric matrix; NaN represents both R NA and NaN.
      integer, intent(in), optional :: rows(:) !! Optional R-style row indices; negative excludes, zero drops.
      integer, intent(in), optional :: cols(:) !! Optional R-style column indices; negative excludes, zero drops.
      real(dp), allocatable :: y(:, :)
      real(dp), allocatable :: a(:, :)
      integer :: j
      a = ms_select(x, rows, cols)
      allocate(y(size(a, 1), size(a, 2)))
      do j = 1, size(a, 1)
         y(j, :) = ms_cumulative(a(j, :), 'prod')
      end do
   end function row_cumprods

   pure function row_cummins(x, rows, cols) result(y)
      real(dp), intent(in) :: x(:, :) !! Numeric matrix; NaN represents both R NA and NaN.
      integer, intent(in), optional :: rows(:) !! Optional R-style row indices; negative excludes, zero drops.
      integer, intent(in), optional :: cols(:) !! Optional R-style column indices; negative excludes, zero drops.
      real(dp), allocatable :: y(:, :)
      real(dp), allocatable :: a(:, :)
      integer :: j
      a = ms_select(x, rows, cols)
      allocate(y(size(a, 1), size(a, 2)))
      do j = 1, size(a, 1)
         y(j, :) = ms_cumulative(a(j, :), 'min')
      end do
   end function row_cummins

   pure function row_cummaxs(x, rows, cols) result(y)
      real(dp), intent(in) :: x(:, :) !! Numeric matrix; NaN represents both R NA and NaN.
      integer, intent(in), optional :: rows(:) !! Optional R-style row indices; negative excludes, zero drops.
      integer, intent(in), optional :: cols(:) !! Optional R-style column indices; negative excludes, zero drops.
      real(dp), allocatable :: y(:, :)
      real(dp), allocatable :: a(:, :)
      integer :: j
      a = ms_select(x, rows, cols)
      allocate(y(size(a, 1), size(a, 2)))
      do j = 1, size(a, 1)
         y(j, :) = ms_cumulative(a(j, :), 'max')
      end do
   end function row_cummaxs

   pure function row_diffs(x, rows, cols, lag, differences) result(y)
      real(dp), intent(in) :: x(:, :) !! Numeric matrix; NaN represents both R NA and NaN.
      integer, intent(in), optional :: rows(:) !! Optional R-style row indices; negative excludes, zero drops.
      integer, intent(in), optional :: cols(:) !! Optional R-style column indices; negative excludes, zero drops.
      integer, intent(in), optional :: lag !! Positive difference lag in observations; default one.
      integer, intent(in), optional :: differences !! Nonnegative difference order; default one.
      real(dp), allocatable :: y(:, :)
      real(dp), allocatable :: a(:, :), d(:)
      integer :: j, n, l, order
      a = ms_select(x, rows, cols)
      l = 1
      order = 1
      if (present(lag)) l = lag
      if (present(differences)) order = differences
      n = max(0, size(a, 2) - max(0, l) * max(0, order))
      if (l < 1 .or. order < 0) n = 0
      allocate(y(size(a, 1), n))
      do j = 1, size(a, 1)
         d = ms_difference(a(j, :), lag, differences)
         y(j, :) = d
      end do
   end function row_diffs

   pure function col_sums2(x, rows, cols, na_rm) result(y)
      real(dp), intent(in) :: x(:, :) !! Numeric matrix; NaN represents both R NA and NaN.
      integer, intent(in), optional :: rows(:) !! Optional R-style row indices; negative excludes, zero drops.
      integer, intent(in), optional :: cols(:) !! Optional R-style column indices; negative excludes, zero drops.
      logical, intent(in), optional :: na_rm !! Remove NaN observations when true; default false.
      real(dp), allocatable :: y(:)
      real(dp), allocatable :: a(:, :), centers(:)
      integer :: j
      a = transpose(ms_select(x, rows, cols))
      allocate(y(size(a, 1)))
      do j = 1, size(a, 1)
         y(j) = ms_sum(a(j, :), na_rm = na_rm)
      end do
   end function col_sums2

   pure function col_means2(x, rows, cols, na_rm, refine) result(y)
      real(dp), intent(in) :: x(:, :) !! Numeric matrix; NaN represents both R NA and NaN.
      integer, intent(in), optional :: rows(:) !! Optional R-style row indices; negative excludes, zero drops.
      integer, intent(in), optional :: cols(:) !! Optional R-style column indices; negative excludes, zero drops.
      logical, intent(in), optional :: na_rm !! Remove NaN observations when true; default false.
      logical, intent(in), optional :: refine !! Refine the mean using a residual correction; default true.
      real(dp), allocatable :: y(:)
      real(dp), allocatable :: a(:, :), centers(:)
      integer :: j
      a = transpose(ms_select(x, rows, cols))
      allocate(y(size(a, 1)))
      do j = 1, size(a, 1)
         y(j) = ms_mean(a(j, :), na_rm = na_rm, refine = refine)
      end do
   end function col_means2

   pure function col_vars(x, rows, cols, na_rm, refine, center) result(y)
      real(dp), intent(in) :: x(:, :) !! Numeric matrix; NaN represents both R NA and NaN.
      integer, intent(in), optional :: rows(:) !! Optional R-style row indices; negative excludes, zero drops.
      integer, intent(in), optional :: cols(:) !! Optional R-style column indices; negative excludes, zero drops.
      logical, intent(in), optional :: na_rm !! Remove NaN observations when true; default false.
      logical, intent(in), optional :: refine !! Refine the mean using a residual correction; default true.
      real(dp), intent(in), optional :: center(:) !! Means for each original row or column, before selection.
      real(dp), allocatable :: y(:)
      real(dp), allocatable :: a(:, :), centers(:)
      integer :: j
      a = transpose(ms_select(x, rows, cols))
      allocate(y(size(a, 1)))
      if (present(center)) then
         centers = ms_vector(center, cols)
         if (size(centers) /= size(a, 1)) then
            y = ms_nan()
            return
         end if
      end if
      do j = 1, size(a, 1)
         if (present(center)) then
            y(j) = ms_var(a(j, :), na_rm = na_rm, refine = refine, center = centers(j))
         else
            y(j) = ms_var(a(j, :), na_rm = na_rm, refine = refine)
         end if
      end do
   end function col_vars

   pure function col_sds(x, rows, cols, na_rm, refine, center) result(y)
      real(dp), intent(in) :: x(:, :) !! Numeric matrix; NaN represents both R NA and NaN.
      integer, intent(in), optional :: rows(:) !! Optional R-style row indices; negative excludes, zero drops.
      integer, intent(in), optional :: cols(:) !! Optional R-style column indices; negative excludes, zero drops.
      logical, intent(in), optional :: na_rm !! Remove NaN observations when true; default false.
      logical, intent(in), optional :: refine !! Refine the mean using a residual correction; default true.
      real(dp), intent(in), optional :: center(:) !! Means for each original row or column, before selection.
      real(dp), allocatable :: y(:)
      real(dp), allocatable :: a(:, :), centers(:)
      integer :: j
      a = transpose(ms_select(x, rows, cols))
      allocate(y(size(a, 1)))
      if (present(center)) then
         centers = ms_vector(center, cols)
         if (size(centers) /= size(a, 1)) then
            y = ms_nan()
            return
         end if
      end if
      do j = 1, size(a, 1)
         if (present(center)) then
            y(j) = ms_var(a(j, :), na_rm = na_rm, refine = refine, center = centers(j))
         else
            y(j) = ms_var(a(j, :), na_rm = na_rm, refine = refine)
         end if
         y(j) = sqrt(y(j))
      end do
   end function col_sds

   pure function col_medians(x, rows, cols, na_rm) result(y)
      real(dp), intent(in) :: x(:, :) !! Numeric matrix; NaN represents both R NA and NaN.
      integer, intent(in), optional :: rows(:) !! Optional R-style row indices; negative excludes, zero drops.
      integer, intent(in), optional :: cols(:) !! Optional R-style column indices; negative excludes, zero drops.
      logical, intent(in), optional :: na_rm !! Remove NaN observations when true; default false.
      real(dp), allocatable :: y(:)
      real(dp), allocatable :: a(:, :), centers(:)
      integer :: j
      a = transpose(ms_select(x, rows, cols))
      allocate(y(size(a, 1)))
      do j = 1, size(a, 1)
         y(j) = ms_median(a(j, :), na_rm = na_rm)
      end do
   end function col_medians

   pure function col_mads(x, rows, cols, na_rm, constant, center) result(y)
      real(dp), intent(in) :: x(:, :) !! Numeric matrix; NaN represents both R NA and NaN.
      integer, intent(in), optional :: rows(:) !! Optional R-style row indices; negative excludes, zero drops.
      integer, intent(in), optional :: cols(:) !! Optional R-style column indices; negative excludes, zero drops.
      logical, intent(in), optional :: na_rm !! Remove NaN observations when true; default false.
      real(dp), intent(in), optional :: constant !! MAD scale multiplier; default 1.4826.
      real(dp), intent(in), optional :: center(:) !! Medians for each original row or column, before selection.
      real(dp), allocatable :: y(:)
      real(dp), allocatable :: a(:, :), centers(:)
      integer :: j
      a = transpose(ms_select(x, rows, cols))
      allocate(y(size(a, 1)))
      if (present(center)) then
         centers = ms_vector(center, cols)
         if (size(centers) /= size(a, 1)) then
            y = ms_nan()
            return
         end if
      end if
      do j = 1, size(a, 1)
         if (present(center)) then
            y(j) = ms_mad(a(j, :), na_rm = na_rm, constant = constant, center = centers(j))
         else
            y(j) = ms_mad(a(j, :), na_rm = na_rm, constant = constant)
         end if
      end do
   end function col_mads

   pure function col_iqrs(x, rows, cols, na_rm, type) result(y)
      real(dp), intent(in) :: x(:, :) !! Numeric matrix; NaN represents both R NA and NaN.
      integer, intent(in), optional :: rows(:) !! Optional R-style row indices; negative excludes, zero drops.
      integer, intent(in), optional :: cols(:) !! Optional R-style column indices; negative excludes, zero drops.
      logical, intent(in), optional :: na_rm !! Remove NaN observations when true; default false.
      integer, intent(in), optional :: type !! Hyndman-Fan quantile type in 1 through 9; default seven.
      real(dp), allocatable :: y(:)
      real(dp), allocatable :: a(:, :), centers(:)
      integer :: j
      a = transpose(ms_select(x, rows, cols))
      allocate(y(size(a, 1)))
      do j = 1, size(a, 1)
         y(j) = ms_quantile(a(j, :), 0.75_dp, na_rm, type) - ms_quantile(a(j, :), 0.25_dp, na_rm, type)
      end do
   end function col_iqrs

   pure function col_mins(x, rows, cols, na_rm) result(y)
      real(dp), intent(in) :: x(:, :) !! Numeric matrix; NaN represents both R NA and NaN.
      integer, intent(in), optional :: rows(:) !! Optional R-style row indices; negative excludes, zero drops.
      integer, intent(in), optional :: cols(:) !! Optional R-style column indices; negative excludes, zero drops.
      logical, intent(in), optional :: na_rm !! Remove NaN observations when true; default false.
      real(dp), allocatable :: y(:)
      real(dp), allocatable :: a(:, :), centers(:)
      integer :: j
      a = transpose(ms_select(x, rows, cols))
      allocate(y(size(a, 1)))
      do j = 1, size(a, 1)
         y(j) = ms_extreme(a(j, :), na_rm = na_rm, maximum = .false.)
      end do
   end function col_mins

   pure function col_maxs(x, rows, cols, na_rm) result(y)
      real(dp), intent(in) :: x(:, :) !! Numeric matrix; NaN represents both R NA and NaN.
      integer, intent(in), optional :: rows(:) !! Optional R-style row indices; negative excludes, zero drops.
      integer, intent(in), optional :: cols(:) !! Optional R-style column indices; negative excludes, zero drops.
      logical, intent(in), optional :: na_rm !! Remove NaN observations when true; default false.
      real(dp), allocatable :: y(:)
      real(dp), allocatable :: a(:, :), centers(:)
      integer :: j
      a = transpose(ms_select(x, rows, cols))
      allocate(y(size(a, 1)))
      do j = 1, size(a, 1)
         y(j) = ms_extreme(a(j, :), na_rm = na_rm, maximum = .true.)
      end do
   end function col_maxs

   pure function col_log_sum_exps(x, rows, cols, na_rm) result(y)
      real(dp), intent(in) :: x(:, :) !! Numeric matrix; NaN represents both R NA and NaN.
      integer, intent(in), optional :: rows(:) !! Optional R-style row indices; negative excludes, zero drops.
      integer, intent(in), optional :: cols(:) !! Optional R-style column indices; negative excludes, zero drops.
      logical, intent(in), optional :: na_rm !! Remove NaN observations when true; default false.
      real(dp), allocatable :: y(:)
      real(dp), allocatable :: a(:, :), centers(:)
      integer :: j
      a = transpose(ms_select(x, rows, cols))
      allocate(y(size(a, 1)))
      do j = 1, size(a, 1)
         y(j) = ms_log_sum_exp(a(j, :), na_rm = na_rm)
      end do
   end function col_log_sum_exps

   pure function col_counts(x, rows, cols, na_rm, value) result(y)
      real(dp), intent(in) :: x(:, :) !! Numeric matrix; NaN represents both R NA and NaN.
      integer, intent(in), optional :: rows(:) !! Optional R-style row indices; negative excludes, zero drops.
      integer, intent(in), optional :: cols(:) !! Optional R-style column indices; negative excludes, zero drops.
      logical, intent(in), optional :: na_rm !! Remove NaN observations when true; default false.
      real(dp), intent(in), optional :: value !! Exact target value; default one, and NaN matches missingness.
      real(dp), allocatable :: y(:)
      real(dp), allocatable :: a(:, :), centers(:)
      integer :: j
      a = transpose(ms_select(x, rows, cols))
      allocate(y(size(a, 1)))
      do j = 1, size(a, 1)
         y(j) = ms_count(a(j, :), na_rm = na_rm, value = value)
      end do
   end function col_counts

   pure function col_alls(x, rows, cols, na_rm, value) result(y)
      real(dp), intent(in) :: x(:, :) !! Numeric matrix; NaN represents both R NA and NaN.
      integer, intent(in), optional :: rows(:) !! Optional R-style row indices; negative excludes, zero drops.
      integer, intent(in), optional :: cols(:) !! Optional R-style column indices; negative excludes, zero drops.
      logical, intent(in), optional :: na_rm !! Remove NaN observations when true; default false.
      real(dp), intent(in), optional :: value !! Explicit numeric target; omitted tests nonzero, and NaN matches missingness.
      real(dp), allocatable :: y(:)
      real(dp), allocatable :: a(:, :), centers(:)
      integer :: j
      a = transpose(ms_select(x, rows, cols))
      allocate(y(size(a, 1)))
      do j = 1, size(a, 1)
         y(j) = ms_all(a(j, :), na_rm = na_rm, value = value)
      end do
   end function col_alls

   pure function col_anys(x, rows, cols, na_rm, value) result(y)
      real(dp), intent(in) :: x(:, :) !! Numeric matrix; NaN represents both R NA and NaN.
      integer, intent(in), optional :: rows(:) !! Optional R-style row indices; negative excludes, zero drops.
      integer, intent(in), optional :: cols(:) !! Optional R-style column indices; negative excludes, zero drops.
      logical, intent(in), optional :: na_rm !! Remove NaN observations when true; default false.
      real(dp), intent(in), optional :: value !! Explicit numeric target; omitted tests nonzero, and NaN matches missingness.
      real(dp), allocatable :: y(:)
      real(dp), allocatable :: a(:, :), centers(:)
      integer :: j
      a = transpose(ms_select(x, rows, cols))
      allocate(y(size(a, 1)))
      do j = 1, size(a, 1)
         y(j) = ms_any(a(j, :), na_rm = na_rm, value = value)
      end do
   end function col_anys

   pure function col_prods(x, rows, cols, na_rm, method) result(y)
      real(dp), intent(in) :: x(:, :) !! Numeric matrix; NaN represents both R NA and NaN.
      integer, intent(in), optional :: rows(:) !! Optional R-style row indices; negative excludes, zero drops.
      integer, intent(in), optional :: cols(:) !! Optional R-style column indices; negative excludes, zero drops.
      logical, intent(in), optional :: na_rm !! Remove NaN observations when true; default false.
      character(*), intent(in), optional :: method !! 'direct' product by default, or 'expSumLog'.
      real(dp), allocatable :: y(:)
      real(dp), allocatable :: a(:, :), centers(:)
      integer :: j
      a = transpose(ms_select(x, rows, cols))
      allocate(y(size(a, 1)))
      do j = 1, size(a, 1)
         y(j) = ms_product(a(j, :), na_rm = na_rm, method = method)
      end do
   end function col_prods

   pure function col_any_missings(x, rows, cols) result(y)
      real(dp), intent(in) :: x(:, :) !! Numeric matrix; NaN represents both R NA and NaN.
      integer, intent(in), optional :: rows(:) !! Optional R-style row indices; negative excludes, zero drops.
      integer, intent(in), optional :: cols(:) !! Optional R-style column indices; negative excludes, zero drops.
      logical, allocatable :: y(:)
      real(dp), allocatable :: a(:, :)
      integer :: j
      a = transpose(ms_select(x, rows, cols))
      allocate(y(size(a, 1)))
      do j = 1, size(a, 1)
         y(j) = any(ieee_is_nan(a(j, :)))
      end do
   end function col_any_missings

   pure function col_any_nas(x, rows, cols) result(y)
      real(dp), intent(in) :: x(:, :) !! Numeric matrix; NaN represents both R NA and NaN.
      integer, intent(in), optional :: rows(:) !! Optional R-style row indices; negative excludes, zero drops.
      integer, intent(in), optional :: cols(:) !! Optional R-style column indices; negative excludes, zero drops.
      logical, allocatable :: y(:)
      real(dp), allocatable :: a(:, :)
      integer :: j
      a = transpose(ms_select(x, rows, cols))
      allocate(y(size(a, 1)))
      do j = 1, size(a, 1)
         y(j) = any(ieee_is_nan(a(j, :)))
      end do
   end function col_any_nas

   pure function col_weighted_means(x, w, rows, cols, na_rm) result(y)
      real(dp), intent(in) :: x(:, :) !! Numeric matrix; NaN represents both R NA and NaN.
      real(dp), intent(in), optional :: w(:) !! Weights, same length as unselected observations; default equal weights.
      integer, intent(in), optional :: rows(:) !! Optional R-style row indices; negative excludes, zero drops.
      integer, intent(in), optional :: cols(:) !! Optional R-style column indices; negative excludes, zero drops.
      logical, intent(in), optional :: na_rm !! Remove NaN observations when true; default false.
      real(dp), allocatable :: y(:)
      real(dp), allocatable :: a(:, :), weights(:), centers(:)
      integer :: j
      a = transpose(ms_select(x, rows, cols))
      allocate(y(size(a, 1)))
      if (present(w)) then
         weights = ms_vector(w, rows)
      else
         allocate(weights(size(a, 2)))
         weights = 1.0_dp
      end if
      if (size(weights) /= size(a, 2)) then
         y = ms_nan()
         return
      end if
      if (any(weights < 0.0_dp)) then
         y = ms_nan()
         return
      end if
      do j = 1, size(a, 1)
         if (.not. present(w)) then
            y(j) = ms_mean(a(j, :), na_rm)
            cycle
         end if
         y(j) = matrix_weighted_mean(a(j, :), weights, na_rm)
      end do
   end function col_weighted_means

   pure function col_weighted_medians(x, w, rows, cols, na_rm, interpolate, ties) result(y)
      real(dp), intent(in) :: x(:, :) !! Numeric matrix; NaN represents both R NA and NaN.
      real(dp), intent(in), optional :: w(:) !! Weights, same length as unselected observations; default equal weights.
      integer, intent(in), optional :: rows(:) !! Optional R-style row indices; negative excludes, zero drops.
      integer, intent(in), optional :: cols(:) !! Optional R-style column indices; negative excludes, zero drops.
      logical, intent(in), optional :: na_rm !! Remove NaN observations when true; default false.
      logical, intent(in), optional :: interpolate !! Interpolate between observations; default true unless ties is supplied.
      character(*), intent(in), optional :: ties !! Half-weight median tie rule: weighted, min, max, or mean.
      real(dp), allocatable :: y(:)
      real(dp), allocatable :: a(:, :), weights(:), centers(:)
      integer :: j
      a = transpose(ms_select(x, rows, cols))
      allocate(y(size(a, 1)))
      if (present(w)) then
         weights = ms_vector(w, rows)
      else
         allocate(weights(size(a, 2)))
         weights = 1.0_dp
      end if
      if (size(weights) /= size(a, 2)) then
         y = ms_nan()
         return
      end if
      if (any(weights < 0.0_dp)) then
         y = ms_nan()
         return
      end if
      do j = 1, size(a, 1)
         if (.not. present(w)) then
            y(j) = ms_median(a(j, :), na_rm)
            cycle
         end if
         y(j) = ms_weighted_median(a(j, :), weights, na_rm, interpolate = interpolate, ties = ties)
      end do
   end function col_weighted_medians

   pure function col_weighted_vars(x, w, rows, cols, na_rm, center) result(y)
      real(dp), intent(in) :: x(:, :) !! Numeric matrix; NaN represents both R NA and NaN.
      real(dp), intent(in), optional :: w(:) !! Weights, same length as unselected observations; default equal weights.
      integer, intent(in), optional :: rows(:) !! Optional R-style row indices; negative excludes, zero drops.
      integer, intent(in), optional :: cols(:) !! Optional R-style column indices; negative excludes, zero drops.
      logical, intent(in), optional :: na_rm !! Remove NaN observations when true; default false.
      real(dp), intent(in), optional :: center(:) !! Weighted means for each original row or column.
      real(dp), allocatable :: y(:)
      real(dp), allocatable :: a(:, :), weights(:), centers(:)
      integer :: j
      a = transpose(ms_select(x, rows, cols))
      allocate(y(size(a, 1)))
      if (present(w)) then
         weights = ms_vector(w, rows)
      else
         allocate(weights(size(a, 2)))
         weights = 1.0_dp
      end if
      if (size(weights) /= size(a, 2)) then
         y = ms_nan()
         return
      end if
      if (present(center)) then
         centers = ms_vector(center, cols)
         if (size(centers) /= size(a, 1)) then
            y = ms_nan()
            return
         end if
      end if
      do j = 1, size(a, 1)
         if (present(center)) then
            y(j) = ms_weighted_var(a(j, :), weights, na_rm, center = centers(j))
         else
            y(j) = ms_weighted_var(a(j, :), weights, na_rm)
         end if
      end do
   end function col_weighted_vars

   pure function col_weighted_sds(x, w, rows, cols, na_rm, center) result(y)
      real(dp), intent(in) :: x(:, :) !! Numeric matrix; NaN represents both R NA and NaN.
      real(dp), intent(in), optional :: w(:) !! Weights, same length as unselected observations; default equal weights.
      integer, intent(in), optional :: rows(:) !! Optional R-style row indices; negative excludes, zero drops.
      integer, intent(in), optional :: cols(:) !! Optional R-style column indices; negative excludes, zero drops.
      logical, intent(in), optional :: na_rm !! Remove NaN observations when true; default false.
      real(dp), intent(in), optional :: center(:) !! Weighted means for each original row or column.
      real(dp), allocatable :: y(:)
      real(dp), allocatable :: a(:, :), weights(:), centers(:)
      integer :: j
      a = transpose(ms_select(x, rows, cols))
      allocate(y(size(a, 1)))
      if (present(w)) then
         weights = ms_vector(w, rows)
      else
         allocate(weights(size(a, 2)))
         weights = 1.0_dp
      end if
      if (size(weights) /= size(a, 2)) then
         y = ms_nan()
         return
      end if
      if (present(center)) then
         centers = ms_vector(center, cols)
         if (size(centers) /= size(a, 1)) then
            y = ms_nan()
            return
         end if
      end if
      do j = 1, size(a, 1)
         if (present(center)) then
            y(j) = ms_weighted_sd(a(j, :), weights, na_rm, center = centers(j))
         else
            y(j) = ms_weighted_sd(a(j, :), weights, na_rm)
         end if
      end do
   end function col_weighted_sds

   pure function col_weighted_mads(x, w, rows, cols, na_rm, constant, center) result(y)
      real(dp), intent(in) :: x(:, :) !! Numeric matrix; NaN represents both R NA and NaN.
      real(dp), intent(in), optional :: w(:) !! Weights, same length as unselected observations; default equal weights.
      integer, intent(in), optional :: rows(:) !! Optional R-style row indices; negative excludes, zero drops.
      integer, intent(in), optional :: cols(:) !! Optional R-style column indices; negative excludes, zero drops.
      logical, intent(in), optional :: na_rm !! Remove NaN observations when true; default false.
      real(dp), intent(in), optional :: constant !! MAD scale multiplier; default 1.4826.
      real(dp), intent(in), optional :: center(:) !! Weighted medians for each original row or column.
      real(dp), allocatable :: y(:)
      real(dp), allocatable :: a(:, :), weights(:), centers(:)
      integer :: j
      a = transpose(ms_select(x, rows, cols))
      allocate(y(size(a, 1)))
      if (present(w)) then
         weights = ms_vector(w, rows)
      else
         allocate(weights(size(a, 2)))
         weights = 1.0_dp
      end if
      if (size(weights) /= size(a, 2)) then
         y = ms_nan()
         return
      end if
      if (present(center)) then
         centers = ms_vector(center, cols)
         if (size(centers) /= size(a, 1)) then
            y = ms_nan()
            return
         end if
      end if
      do j = 1, size(a, 1)
         if (.not. present(w)) then
            if (present(center)) then
               y(j) = ms_weighted_mad(a(j, :), na_rm = na_rm, constant = constant, center = centers(j))
            else
               y(j) = ms_weighted_mad(a(j, :), na_rm = na_rm, constant = constant)
            end if
            cycle
         end if
         if (present(center)) then
            y(j) = ms_weighted_mad(a(j, :), weights, na_rm, constant = constant, center = centers(j))
         else
            y(j) = ms_weighted_mad(a(j, :), weights, na_rm, constant = constant)
         end if
      end do
   end function col_weighted_mads

   pure function col_ranges(x, rows, cols, na_rm) result(y)
      real(dp), intent(in) :: x(:, :) !! Numeric matrix; NaN represents both R NA and NaN.
      integer, intent(in), optional :: rows(:) !! Optional R-style row indices; negative excludes, zero drops.
      integer, intent(in), optional :: cols(:) !! Optional R-style column indices; negative excludes, zero drops.
      logical, intent(in), optional :: na_rm !! Remove NaN observations when true; default false.
      real(dp), allocatable :: y(:, :)
      real(dp), allocatable :: a(:, :)
      integer :: j
      a = transpose(ms_select(x, rows, cols))
      allocate(y(size(a, 1), 2))
      do j = 1, size(a, 1)
         y(j, 1) = ms_extreme(a(j, :), .false., na_rm)
         y(j, 2) = ms_extreme(a(j, :), .true., na_rm)
      end do
   end function col_ranges

   pure function col_quantiles(x, probs, rows, cols, na_rm, type) result(y)
      real(dp), intent(in) :: x(:, :) !! Numeric matrix; NaN represents both R NA and NaN.
      real(dp), intent(in) :: probs(:) !! Quantile probabilities, each in the closed interval zero to one.
      integer, intent(in), optional :: rows(:) !! Optional R-style row indices; negative excludes, zero drops.
      integer, intent(in), optional :: cols(:) !! Optional R-style column indices; negative excludes, zero drops.
      logical, intent(in), optional :: na_rm !! Remove NaN observations when true; default false.
      integer, intent(in), optional :: type !! Hyndman-Fan quantile type in 1 through 9; default seven.
      real(dp), allocatable :: y(:, :)
      real(dp), allocatable :: a(:, :)
      integer :: j, k
      a = transpose(ms_select(x, rows, cols))
      allocate(y(size(a, 1), size(probs)))
      do k = 1, size(probs)
         do j = 1, size(a, 1)
            y(j, k) = ms_quantile(a(j, :), probs(k), na_rm, type)
         end do
      end do
   end function col_quantiles

   pure function col_order_stats(x, which, rows, cols) result(y)
      real(dp), intent(in) :: x(:, :) !! Numeric matrix; NaN represents both R NA and NaN.
      integer, intent(in) :: which !! One-based order statistic rank within the selected observations.
      integer, intent(in), optional :: rows(:) !! Optional R-style row indices; negative excludes, zero drops.
      integer, intent(in), optional :: cols(:) !! Optional R-style column indices; negative excludes, zero drops.
      real(dp), allocatable :: y(:)
      real(dp), allocatable :: a(:, :)
      integer :: j
      a = transpose(ms_select(x, rows, cols))
      allocate(y(size(a, 1)))
      do j = 1, size(a, 1)
         y(j) = ms_order_stat(a(j, :), which)
      end do
   end function col_order_stats

   pure function col_ranks(x, rows, cols, ties_method, preserve_shape) result(y)
      real(dp), intent(in) :: x(:, :) !! Numeric matrix; NaN represents both R NA and NaN.
      integer, intent(in), optional :: rows(:) !! Optional R-style row indices; negative excludes, zero drops.
      integer, intent(in), optional :: cols(:) !! Optional R-style column indices; negative excludes, zero drops.
      character(*), intent(in), optional :: ties_method !! Rank ties: max (default), min, average, first, last, or dense.
      logical, intent(in), optional :: preserve_shape !! True preserves input shape; default returns ranked columns as rows.
      real(dp), allocatable :: y(:, :)
      real(dp), allocatable :: a(:, :)
      integer :: j
      a = transpose(ms_select(x, rows, cols))
      allocate(y(size(a, 1), size(a, 2)))
      do j = 1, size(a, 1)
         y(j, :) = ms_ranks(a(j, :), ties_method)
      end do
      if (present(preserve_shape)) then
         if (preserve_shape) y = transpose(y)
      end if
   end function col_ranks

   pure function col_cumsums(x, rows, cols) result(y)
      real(dp), intent(in) :: x(:, :) !! Numeric matrix; NaN represents both R NA and NaN.
      integer, intent(in), optional :: rows(:) !! Optional R-style row indices; negative excludes, zero drops.
      integer, intent(in), optional :: cols(:) !! Optional R-style column indices; negative excludes, zero drops.
      real(dp), allocatable :: y(:, :)
      real(dp), allocatable :: a(:, :)
      integer :: j
      a = transpose(ms_select(x, rows, cols))
      allocate(y(size(a, 1), size(a, 2)))
      do j = 1, size(a, 1)
         y(j, :) = ms_cumulative(a(j, :), 'sum')
      end do
      y = transpose(y)
   end function col_cumsums

   pure function col_cumprods(x, rows, cols) result(y)
      real(dp), intent(in) :: x(:, :) !! Numeric matrix; NaN represents both R NA and NaN.
      integer, intent(in), optional :: rows(:) !! Optional R-style row indices; negative excludes, zero drops.
      integer, intent(in), optional :: cols(:) !! Optional R-style column indices; negative excludes, zero drops.
      real(dp), allocatable :: y(:, :)
      real(dp), allocatable :: a(:, :)
      integer :: j
      a = transpose(ms_select(x, rows, cols))
      allocate(y(size(a, 1), size(a, 2)))
      do j = 1, size(a, 1)
         y(j, :) = ms_cumulative(a(j, :), 'prod')
      end do
      y = transpose(y)
   end function col_cumprods

   pure function col_cummins(x, rows, cols) result(y)
      real(dp), intent(in) :: x(:, :) !! Numeric matrix; NaN represents both R NA and NaN.
      integer, intent(in), optional :: rows(:) !! Optional R-style row indices; negative excludes, zero drops.
      integer, intent(in), optional :: cols(:) !! Optional R-style column indices; negative excludes, zero drops.
      real(dp), allocatable :: y(:, :)
      real(dp), allocatable :: a(:, :)
      integer :: j
      a = transpose(ms_select(x, rows, cols))
      allocate(y(size(a, 1), size(a, 2)))
      do j = 1, size(a, 1)
         y(j, :) = ms_cumulative(a(j, :), 'min')
      end do
      y = transpose(y)
   end function col_cummins

   pure function col_cummaxs(x, rows, cols) result(y)
      real(dp), intent(in) :: x(:, :) !! Numeric matrix; NaN represents both R NA and NaN.
      integer, intent(in), optional :: rows(:) !! Optional R-style row indices; negative excludes, zero drops.
      integer, intent(in), optional :: cols(:) !! Optional R-style column indices; negative excludes, zero drops.
      real(dp), allocatable :: y(:, :)
      real(dp), allocatable :: a(:, :)
      integer :: j
      a = transpose(ms_select(x, rows, cols))
      allocate(y(size(a, 1), size(a, 2)))
      do j = 1, size(a, 1)
         y(j, :) = ms_cumulative(a(j, :), 'max')
      end do
      y = transpose(y)
   end function col_cummaxs

   pure function col_diffs(x, rows, cols, lag, differences) result(y)
      real(dp), intent(in) :: x(:, :) !! Numeric matrix; NaN represents both R NA and NaN.
      integer, intent(in), optional :: rows(:) !! Optional R-style row indices; negative excludes, zero drops.
      integer, intent(in), optional :: cols(:) !! Optional R-style column indices; negative excludes, zero drops.
      integer, intent(in), optional :: lag !! Positive difference lag in observations; default one.
      integer, intent(in), optional :: differences !! Nonnegative difference order; default one.
      real(dp), allocatable :: y(:, :)
      real(dp), allocatable :: a(:, :), d(:)
      integer :: j, n, l, order
      a = transpose(ms_select(x, rows, cols))
      l = 1
      order = 1
      if (present(lag)) l = lag
      if (present(differences)) order = differences
      n = max(0, size(a, 2) - max(0, l) * max(0, order))
      if (l < 1 .or. order < 0) n = 0
      allocate(y(size(a, 1), n))
      do j = 1, size(a, 1)
         d = ms_difference(a(j, :), lag, differences)
         y(j, :) = d
      end do
      y = transpose(y)
   end function col_diffs

   pure function var_diff(x, idxs, na_rm, diff, trim) result(y)
      real(dp), intent(in) :: x(:) !! Numeric vector; NaN represents both R NA and NaN.
      integer, intent(in), optional :: idxs(:) !! Optional R-style vector indices; negative excludes, zero drops.
      logical, intent(in), optional :: na_rm !! Remove NaN observations when true; default false.
      integer, intent(in), optional :: diff !! Nonnegative repeated difference order; default one.
      real(dp), intent(in), optional :: trim !! Fraction removed from each sorted tail, between zero and one half.
      real(dp) :: y
      y = scale_difference(ms_vector(x, idxs), 'var', na_rm, diff, trim)
   end function var_diff

   pure function row_var_diffs(x, rows, cols, na_rm, diff, trim) result(y)
      real(dp), intent(in) :: x(:, :) !! Numeric matrix; NaN represents both R NA and NaN.
      integer, intent(in), optional :: rows(:) !! Optional R-style row indices; negative excludes, zero drops.
      integer, intent(in), optional :: cols(:) !! Optional R-style column indices; negative excludes, zero drops.
      logical, intent(in), optional :: na_rm !! Remove NaN observations when true; default false.
      integer, intent(in), optional :: diff !! Nonnegative repeated difference order; default one.
      real(dp), intent(in), optional :: trim !! Fraction removed from each sorted tail, between zero and one half.
      real(dp), allocatable :: y(:)
      real(dp), allocatable :: a(:, :)
      integer :: j
      a = ms_select(x, rows, cols)
      allocate(y(size(a, 1)))
      do j = 1, size(a, 1)
         y(j) = scale_difference(a(j, :), 'var', na_rm, diff, trim)
      end do
   end function row_var_diffs

   pure function col_var_diffs(x, rows, cols, na_rm, diff, trim) result(y)
      real(dp), intent(in) :: x(:, :) !! Numeric matrix; NaN represents both R NA and NaN.
      integer, intent(in), optional :: rows(:) !! Optional R-style row indices; negative excludes, zero drops.
      integer, intent(in), optional :: cols(:) !! Optional R-style column indices; negative excludes, zero drops.
      logical, intent(in), optional :: na_rm !! Remove NaN observations when true; default false.
      integer, intent(in), optional :: diff !! Nonnegative repeated difference order; default one.
      real(dp), intent(in), optional :: trim !! Fraction removed from each sorted tail, between zero and one half.
      real(dp), allocatable :: y(:)
      real(dp), allocatable :: a(:, :)
      integer :: j
      a = transpose(ms_select(x, rows, cols))
      allocate(y(size(a, 1)))
      do j = 1, size(a, 1)
         y(j) = scale_difference(a(j, :), 'var', na_rm, diff, trim)
      end do
   end function col_var_diffs

   pure function sd_diff(x, idxs, na_rm, diff, trim) result(y)
      real(dp), intent(in) :: x(:) !! Numeric vector; NaN represents both R NA and NaN.
      integer, intent(in), optional :: idxs(:) !! Optional R-style vector indices; negative excludes, zero drops.
      logical, intent(in), optional :: na_rm !! Remove NaN observations when true; default false.
      integer, intent(in), optional :: diff !! Nonnegative repeated difference order; default one.
      real(dp), intent(in), optional :: trim !! Fraction removed from each sorted tail, between zero and one half.
      real(dp) :: y
      y = scale_difference(ms_vector(x, idxs), 'sd', na_rm, diff, trim)
   end function sd_diff

   pure function row_sd_diffs(x, rows, cols, na_rm, diff, trim) result(y)
      real(dp), intent(in) :: x(:, :) !! Numeric matrix; NaN represents both R NA and NaN.
      integer, intent(in), optional :: rows(:) !! Optional R-style row indices; negative excludes, zero drops.
      integer, intent(in), optional :: cols(:) !! Optional R-style column indices; negative excludes, zero drops.
      logical, intent(in), optional :: na_rm !! Remove NaN observations when true; default false.
      integer, intent(in), optional :: diff !! Nonnegative repeated difference order; default one.
      real(dp), intent(in), optional :: trim !! Fraction removed from each sorted tail, between zero and one half.
      real(dp), allocatable :: y(:)
      real(dp), allocatable :: a(:, :)
      integer :: j
      a = ms_select(x, rows, cols)
      allocate(y(size(a, 1)))
      do j = 1, size(a, 1)
         y(j) = scale_difference(a(j, :), 'sd', na_rm, diff, trim)
      end do
   end function row_sd_diffs

   pure function col_sd_diffs(x, rows, cols, na_rm, diff, trim) result(y)
      real(dp), intent(in) :: x(:, :) !! Numeric matrix; NaN represents both R NA and NaN.
      integer, intent(in), optional :: rows(:) !! Optional R-style row indices; negative excludes, zero drops.
      integer, intent(in), optional :: cols(:) !! Optional R-style column indices; negative excludes, zero drops.
      logical, intent(in), optional :: na_rm !! Remove NaN observations when true; default false.
      integer, intent(in), optional :: diff !! Nonnegative repeated difference order; default one.
      real(dp), intent(in), optional :: trim !! Fraction removed from each sorted tail, between zero and one half.
      real(dp), allocatable :: y(:)
      real(dp), allocatable :: a(:, :)
      integer :: j
      a = transpose(ms_select(x, rows, cols))
      allocate(y(size(a, 1)))
      do j = 1, size(a, 1)
         y(j) = scale_difference(a(j, :), 'sd', na_rm, diff, trim)
      end do
   end function col_sd_diffs

   pure function mad_diff(x, idxs, na_rm, diff, trim, constant) result(y)
      real(dp), intent(in) :: x(:) !! Numeric vector; NaN represents both R NA and NaN.
      integer, intent(in), optional :: idxs(:) !! Optional R-style vector indices; negative excludes, zero drops.
      logical, intent(in), optional :: na_rm !! Remove NaN observations when true; default false.
      integer, intent(in), optional :: diff !! Nonnegative repeated difference order; default one.
      real(dp), intent(in), optional :: trim !! Fraction removed from each sorted tail, between zero and one half.
      real(dp), intent(in), optional :: constant !! MAD scale multiplier; default 1.4826.
      real(dp) :: y
      y = scale_difference(ms_vector(x, idxs), 'mad', na_rm, diff, trim, constant = constant)
   end function mad_diff

   pure function row_mad_diffs(x, rows, cols, na_rm, diff, trim, constant) result(y)
      real(dp), intent(in) :: x(:, :) !! Numeric matrix; NaN represents both R NA and NaN.
      integer, intent(in), optional :: rows(:) !! Optional R-style row indices; negative excludes, zero drops.
      integer, intent(in), optional :: cols(:) !! Optional R-style column indices; negative excludes, zero drops.
      logical, intent(in), optional :: na_rm !! Remove NaN observations when true; default false.
      integer, intent(in), optional :: diff !! Nonnegative repeated difference order; default one.
      real(dp), intent(in), optional :: trim !! Fraction removed from each sorted tail, between zero and one half.
      real(dp), intent(in), optional :: constant !! MAD scale multiplier; default 1.4826.
      real(dp), allocatable :: y(:)
      real(dp), allocatable :: a(:, :)
      integer :: j
      a = ms_select(x, rows, cols)
      allocate(y(size(a, 1)))
      do j = 1, size(a, 1)
         y(j) = scale_difference(a(j, :), 'mad', na_rm, diff, trim, constant = constant)
      end do
   end function row_mad_diffs

   pure function col_mad_diffs(x, rows, cols, na_rm, diff, trim, constant) result(y)
      real(dp), intent(in) :: x(:, :) !! Numeric matrix; NaN represents both R NA and NaN.
      integer, intent(in), optional :: rows(:) !! Optional R-style row indices; negative excludes, zero drops.
      integer, intent(in), optional :: cols(:) !! Optional R-style column indices; negative excludes, zero drops.
      logical, intent(in), optional :: na_rm !! Remove NaN observations when true; default false.
      integer, intent(in), optional :: diff !! Nonnegative repeated difference order; default one.
      real(dp), intent(in), optional :: trim !! Fraction removed from each sorted tail, between zero and one half.
      real(dp), intent(in), optional :: constant !! MAD scale multiplier; default 1.4826.
      real(dp), allocatable :: y(:)
      real(dp), allocatable :: a(:, :)
      integer :: j
      a = transpose(ms_select(x, rows, cols))
      allocate(y(size(a, 1)))
      do j = 1, size(a, 1)
         y(j) = scale_difference(a(j, :), 'mad', na_rm, diff, trim, constant = constant)
      end do
   end function col_mad_diffs

   pure function iqr_diff(x, idxs, na_rm, diff, trim, type) result(y)
      real(dp), intent(in) :: x(:) !! Numeric vector; NaN represents both R NA and NaN.
      integer, intent(in), optional :: idxs(:) !! Optional R-style vector indices; negative excludes, zero drops.
      logical, intent(in), optional :: na_rm !! Remove NaN observations when true; default false.
      integer, intent(in), optional :: diff !! Nonnegative repeated difference order; default one.
      real(dp), intent(in), optional :: trim !! Fraction removed from each sorted tail, between zero and one half.
      integer, intent(in), optional :: type !! Hyndman-Fan quantile type in 1 through 9; default seven.
      real(dp) :: y
      y = scale_difference(ms_vector(x, idxs), 'iqr', na_rm, diff, trim, type = type)
   end function iqr_diff

   pure function row_iqr_diffs(x, rows, cols, na_rm, diff, trim, type) result(y)
      real(dp), intent(in) :: x(:, :) !! Numeric matrix; NaN represents both R NA and NaN.
      integer, intent(in), optional :: rows(:) !! Optional R-style row indices; negative excludes, zero drops.
      integer, intent(in), optional :: cols(:) !! Optional R-style column indices; negative excludes, zero drops.
      logical, intent(in), optional :: na_rm !! Remove NaN observations when true; default false.
      integer, intent(in), optional :: diff !! Nonnegative repeated difference order; default one.
      real(dp), intent(in), optional :: trim !! Fraction removed from each sorted tail, between zero and one half.
      integer, intent(in), optional :: type !! Hyndman-Fan quantile type in 1 through 9; default seven.
      real(dp), allocatable :: y(:)
      real(dp), allocatable :: a(:, :)
      integer :: j
      a = ms_select(x, rows, cols)
      allocate(y(size(a, 1)))
      do j = 1, size(a, 1)
         y(j) = scale_difference(a(j, :), 'iqr', na_rm, diff, trim, type = type)
      end do
   end function row_iqr_diffs

   pure function col_iqr_diffs(x, rows, cols, na_rm, diff, trim, type) result(y)
      real(dp), intent(in) :: x(:, :) !! Numeric matrix; NaN represents both R NA and NaN.
      integer, intent(in), optional :: rows(:) !! Optional R-style row indices; negative excludes, zero drops.
      integer, intent(in), optional :: cols(:) !! Optional R-style column indices; negative excludes, zero drops.
      logical, intent(in), optional :: na_rm !! Remove NaN observations when true; default false.
      integer, intent(in), optional :: diff !! Nonnegative repeated difference order; default one.
      real(dp), intent(in), optional :: trim !! Fraction removed from each sorted tail, between zero and one half.
      integer, intent(in), optional :: type !! Hyndman-Fan quantile type in 1 through 9; default seven.
      real(dp), allocatable :: y(:)
      real(dp), allocatable :: a(:, :)
      integer :: j
      a = transpose(ms_select(x, rows, cols))
      allocate(y(size(a, 1)))
      do j = 1, size(a, 1)
         y(j) = scale_difference(a(j, :), 'iqr', na_rm, diff, trim, type = type)
      end do
   end function col_iqr_diffs
   pure real(dp) function matrix_weighted_mean(x, w, na_rm) result(y)
      real(dp), intent(in) :: x(:) !! One selected matrix row or column.
      real(dp), intent(in) :: w(:) !! Nonnegative matching weights, possibly NaN or infinity.
      logical, intent(in), optional :: na_rm !! Remove missing x observations before normalizing weights.
      real(dp), allocatable :: a(:), weights(:)
      logical, allocatable :: keep(:)
      real(dp) :: total
      allocate(keep(size(x)))
      keep = ieee_is_nan(w) .or. w /= 0.0_dp
      if (ms_remove(na_rm)) keep = keep .and. .not. ieee_is_nan(x)
      a = pack(x, keep)
      weights = pack(w, keep)
      total = ms_sum(weights)
      y = ms_nan()
      if (total == 0.0_dp .or. ieee_is_nan(total)) return
      y = ms_sum((weights / total) * a)
   end function matrix_weighted_mean

   pure real(dp) function scale_difference(x, operation, na_rm, diff, trim, constant, type) result(y)
      real(dp), intent(in) :: x(:) !! Observations, in original sequential order.
      character(*), intent(in) :: operation !! One of var, sd, mad, or iqr.
      logical, intent(in), optional :: na_rm !! Remove NaNs before differencing when true.
      integer, intent(in), optional :: diff !! Nonnegative repeated difference order; default one.
      real(dp), intent(in), optional :: trim !! Fraction trimmed from each sorted tail; default zero.
      real(dp), intent(in), optional :: constant !! Normal-consistency multiplier for MAD; default 1.4826.
      integer, intent(in), optional :: type !! IQR Hyndman-Fan quantile type; default seven.
      real(dp), allocatable :: a(:)
      real(dp) :: fraction
      integer :: d, n, k
      d = 1
      if (present(diff)) d = diff
      fraction = 0.0_dp
      if (present(trim)) fraction = trim
      y = ms_nan()
      if (d < 0 .or. fraction < 0.0_dp .or. fraction > 0.5_dp .or. ieee_is_nan(fraction)) return
      a = ms_vector(x, na_rm = na_rm)
      n = size(a)
      if (any(ieee_is_nan(a)) .or. n == 0) return
      if (operation == 'iqr' .and. n == 1) then
         y = 0.0_dp
         return
      end if
      if ((operation == 'var' .or. operation == 'sd') .and. n < 2) return
      if (d > 0) a = ms_difference(a, differences = d)
      n = size(a)
      if (d > 0 .and. n == 1) then
         if (operation == 'iqr') y = 0.0_dp
         return
      end if
      if (fraction > 0.0_dp .and. n > 0) then
         call r_sort_values_in_place(a)
         k = floor(real(n, dp) * fraction)
         ! Upstream trim = 0.5 uses the middle pair for even-length samples.
         if (2 * k == n) k = max(0, k - 1)
         a = a(k + 1:n - k)
      end if
      select case (operation)
      case ('var')
         y = ms_var(a) / (2.0_dp**d)
      case ('sd')
         y = sqrt(ms_var(a)) / (sqrt(2.0_dp)**d)
      case ('mad')
         y = ms_mad(a, constant = constant) / (sqrt(2.0_dp)**d)
      case ('iqr')
         y = (ms_quantile(a, 0.75_dp, type = type) - ms_quantile(a, 0.25_dp, type = type)) / &
            (sqrt(2.0_dp)**d)
      end select
   end function scale_difference

   pure function bin_counts(x, bx, idxs, right) result(y)
      real(dp), intent(in) :: x(:) !! Numeric positions; missing and out-of-range positions are ignored.
      real(dp), intent(in) :: bx(:) !! Finite nondecreasing bin boundaries, at least two elements.
      integer, intent(in), optional :: idxs(:) !! Optional R-style position selection.
      logical, intent(in), optional :: right !! Use (left,right] intervals if true; default [left,right).
      integer, allocatable :: y(:)
      real(dp), allocatable :: a(:)
      integer :: i, k
      allocate(y(max(0, size(bx) - 1)))
      y = 0
      a = ms_vector(x, idxs)
      do i = 1, size(a)
         k = bin_index(a(i), bx, right)
         if (k > 0) y(k) = y(k) + 1
      end do
   end function bin_counts

   pure integer function bin_index(x, bx, right) result(k)
      real(dp), intent(in) :: x !! Coordinate to classify; NaN maps to no bin.
      real(dp), intent(in) :: bx(:) !! Finite nondecreasing boundaries.
      logical, intent(in), optional :: right !! Select right-closed rather than left-closed bins.
      integer :: j
      k = 0
      if (ieee_is_nan(x)) return
      if (any(.not. ieee_is_finite(bx))) return
      do j = 1, size(bx) - 1
         if (ms_remove(right)) then
            if (x > bx(j) .and. x <= bx(j + 1)) k = j
         else
            if (x >= bx(j) .and. x < bx(j + 1)) k = j
         end if
         if (k > 0) return
      end do
   end function bin_index

   pure subroutine bin_means(y, x, bx, means, counts, idxs, na_rm, right)
      real(dp), intent(in) :: y(:) !! Values averaged in bins; length must equal x.
      real(dp), intent(in) :: x(:) !! Bin coordinates; NaN coordinates are always removed.
      real(dp), intent(in) :: bx(:) !! Finite nondecreasing bin boundaries.
      real(dp), allocatable, intent(out) :: means(:) !! Means of each bin; empty bins contain NaN.
      integer, allocatable, intent(out), optional :: counts(:) !! Number of retained values in each bin.
      integer, intent(in), optional :: idxs(:) !! Optional R-style joint selection of x and y.
      logical, intent(in), optional :: na_rm !! Remove missing y values; default true, unlike most other routines.
      logical, intent(in), optional :: right !! Use right-closed intervals when true.
      real(dp), allocatable :: a(:), positions(:), values(:)
      integer, allocatable :: n(:)
      integer :: i, k, j
      logical :: remove
      allocate(means(max(0, size(bx) - 1)), n(max(0, size(bx) - 1)))
      means = ms_nan()
      n = 0
      if (size(x) /= size(y)) then
         if (present(counts)) counts = n
         return
      end if
      a = ms_vector(y, idxs)
      positions = ms_vector(x, idxs)
      remove = .true.
      if (present(na_rm)) remove = na_rm
      do j = 1, size(means)
         allocate(values(size(a)))
         do i = 1, size(a)
            if (remove .and. ieee_is_nan(a(i))) cycle
            k = bin_index(positions(i), bx, right)
            if (k /= j) cycle
            n(j) = n(j) + 1
            values(n(j)) = a(i)
         end do
         means(j) = ms_mean(values(:n(j)))
         deallocate(values)
      end do
      if (present(counts)) counts = n
   end subroutine bin_means

   pure function sign_tabulate(x, idxs) result(y)
      real(dp), intent(in) :: x(:) !! Numeric observations, including missing or infinite values.
      integer, intent(in), optional :: idxs(:) !! Optional R-style vector selection.
      integer :: y(6)
      real(dp), allocatable :: a(:)
      a = ms_vector(x, idxs)
      y(1) = count(a < 0.0_dp)
      y(2) = count(a == 0.0_dp)
      y(3) = count(a > 0.0_dp)
      y(4) = count(ieee_is_nan(a))
      y(5) = count(a == ms_inf(.true.))
      y(6) = count(a == ms_inf())
   end function sign_tabulate

   pure function index_by_row(dim, idxs) result(y)
      integer, intent(in) :: dim(2) !! Nonnegative matrix dimensions: row count and column count.
      integer, intent(in), optional :: idxs(:) !! Optional positive row-major positions to convert.
      integer, allocatable :: y(:)
      integer, allocatable :: indices(:)
      integer :: n, i, p
      if (any(dim < 0)) then
         allocate(y(0))
         return
      end if
      if (dim(1) > huge(n) / max(1, dim(2))) then
         allocate(y(0))
         return
      end if
      n = dim(1) * dim(2)
      if (present(idxs)) then
         indices = idxs
      else
         indices = [(i, i = 1, n)]
      end if
      allocate(y(size(indices)))
      y = ms_na_index
      do i = 1, size(indices)
         p = indices(i)
         if (p < 1 .or. p > n) cycle
         y(i) = (p - 1) / dim(2) + dim(1) * mod(p - 1, dim(2)) + 1
      end do
   end function index_by_row

   pure function row_collapse(x, idxs, rows) result(y)
      real(dp), intent(in) :: x(:, :) !! Matrix from which one column entry is chosen in each row.
      integer, intent(in) :: idxs(:) !! One-based column choices, recycled across original rows.
      integer, intent(in), optional :: rows(:) !! Optional R-style selection of original rows.
      real(dp), allocatable :: y(:)
      integer, allocatable :: ir(:)
      integer :: j, k, r
      ir = ms_indices(size(x, 1), rows)
      allocate(y(size(ir)))
      y = ms_nan()
      if (size(idxs) == 0) return
      do j = 1, size(ir)
         r = ir(j)
         if (r == ms_na_index) cycle
         k = idxs(mod(r - 1, size(idxs)) + 1)
         if (k >= 1 .and. k <= size(x, 2)) y(j) = x(r, k)
      end do
   end function row_collapse

   pure function col_collapse(x, idxs, cols) result(y)
      real(dp), intent(in) :: x(:, :) !! Matrix from which one row entry is chosen in each column.
      integer, intent(in) :: idxs(:) !! One-based row choices, recycled across original columns.
      integer, intent(in), optional :: cols(:) !! Optional R-style selection of original columns.
      real(dp), allocatable :: y(:)
      y = row_collapse(transpose(x), idxs, cols)
   end function col_collapse

   pure function tabulate_values(x, values) result(y)
      real(dp), intent(in) :: x(:, :) !! Selected matrix used to infer tabulation categories.
      real(dp), intent(in), optional :: values(:) !! Explicit categories; absent means sorted unique values with NaN last.
      real(dp), allocatable :: y(:)
      real(dp), allocatable :: a(:), distinct(:)
      integer :: i, k
      if (present(values)) then
         y = values
         return
      end if
      a = pack(x, .not. ieee_is_nan(x))
      call r_sort_values_in_place(a)
      allocate(distinct(size(a) + 1))
      k = 0
      do i = 1, size(a)
         if (k > 0) then
            if (a(i) == distinct(k)) cycle
         end if
         k = k + 1
         distinct(k) = a(i)
      end do
      if (any(ieee_is_nan(x))) then
         k = k + 1
         distinct(k) = ms_nan()
      end if
      y = distinct(:k)
   end function tabulate_values

   pure function row_tabulates(x, values, rows, cols) result(y)
      real(dp), intent(in) :: x(:, :) !! Numeric representation of integer or logical categories; NaN represents NA.
      real(dp), intent(in), optional :: values(:) !! Categories to count; default sorted unique selected values, NaN last.
      integer, intent(in), optional :: rows(:) !! Optional R-style row selection.
      integer, intent(in), optional :: cols(:) !! Optional R-style column selection.
      integer, allocatable :: y(:, :)
      real(dp), allocatable :: a(:, :), v(:)
      integer :: i, j
      a = ms_select(x, rows, cols)
      v = tabulate_values(a, values)
      allocate(y(size(a, 1), size(v)))
      do j = 1, size(v)
         do i = 1, size(a, 1)
            y(i, j) = int(ms_count(a(i, :), v(j), .true.))
         end do
      end do
   end function row_tabulates

   pure function col_tabulates(x, values, rows, cols) result(y)
      real(dp), intent(in) :: x(:, :) !! Numeric representation of integer or logical categories.
      real(dp), intent(in), optional :: values(:) !! Categories to count; default sorted unique values.
      integer, intent(in), optional :: rows(:) !! Optional R-style row selection.
      integer, intent(in), optional :: cols(:) !! Optional R-style column selection.
      integer, allocatable :: y(:, :)
      y = row_tabulates(transpose(x), values, cols, rows)
   end function col_tabulates

   pure function row_avgs_per_col_set(x, sets, rows, na_rm, w) result(y)
      real(dp), intent(in) :: x(:, :) !! Matrix whose rows are averaged over column sets.
      integer, intent(in) :: sets(:, :) !! Each column lists one set of positive column indices; nonpositive entries are padding.
      integer, intent(in), optional :: rows(:) !! Optional R-style row selection.
      logical, intent(in), optional :: na_rm !! Remove NaNs; default true when x contains any missing value.
      real(dp), intent(in), optional :: w(:, :) !! Optional nonnegative weights with the same shape as x.
      real(dp), allocatable :: y(:, :)
      real(dp), allocatable :: a(:, :), weights(:, :)
      integer, allocatable :: members(:)
      integer :: i, j
      logical :: remove
      a = ms_select(x, rows)
      if (present(w)) then
         if (any(shape(w) /= shape(x))) then
            allocate(y(size(a, 1), size(sets, 2)))
            y = ms_nan()
            return
         end if
         weights = ms_select(w, rows)
      end if
      remove = any(ieee_is_nan(a))
      if (present(na_rm)) remove = na_rm
      allocate(y(size(a, 1), size(sets, 2)))
      y = ms_nan()
      do j = 1, size(sets, 2)
         members = pack(sets(:, j), sets(:, j) > 0 .and. sets(:, j) <= size(x, 2))
         do i = 1, size(a, 1)
            if (present(w)) then
               y(i, j) = matrix_weighted_mean(a(i, members), weights(i, members), remove)
            else
               y(i, j) = ms_mean(a(i, members), remove)
            end if
         end do
      end do
   end function row_avgs_per_col_set

   pure function col_avgs_per_row_set(x, sets, cols, na_rm, w) result(y)
      real(dp), intent(in) :: x(:, :) !! Matrix whose columns are averaged over row sets.
      integer, intent(in) :: sets(:, :) !! Each column lists one set of positive row indices; nonpositive entries are padding.
      integer, intent(in), optional :: cols(:) !! Optional R-style column selection.
      logical, intent(in), optional :: na_rm !! Remove NaNs; default true if x contains missing values.
      real(dp), intent(in), optional :: w(:, :) !! Optional nonnegative weights with the same shape as x.
      real(dp), allocatable :: y(:, :)
      if (present(w)) then
         y = transpose(row_avgs_per_col_set(transpose(x), sets, cols, na_rm, transpose(w)))
      else
         y = transpose(row_avgs_per_col_set(transpose(x), sets, cols, na_rm))
      end if
   end function col_avgs_per_row_set

   pure function x_op_y(x, y, op, xrows, xcols, yidxs, commute, na_rm) result(z)
      real(dp), intent(in) :: x(:, :) !! Matrix operand; selections are applied before vector recycling.
      real(dp), intent(in) :: y(:) !! Vector operand, recycled in column-major order.
      character(*), intent(in) :: op !! Arithmetic operation: +, -, *, or /.
      integer, intent(in), optional :: xrows(:) !! Optional R-style selection of matrix rows.
      integer, intent(in), optional :: xcols(:) !! Optional R-style selection of matrix columns.
      integer, intent(in), optional :: yidxs(:) !! Optional R-style selection of vector entries.
      logical, intent(in), optional :: commute !! Compute y OP x rather than x OP y when true.
      logical, intent(in), optional :: na_rm !! For + and *, one missing operand is replaced by the other; default false.
      real(dp), allocatable :: z(:, :)
      real(dp), allocatable :: a(:, :), b(:)
      integer :: i, j, k
      a = ms_select(x, xrows, xcols)
      b = ms_vector(y, yidxs)
      allocate(z(size(a, 1), size(a, 2)))
      z = ms_nan()
      if (size(b) == 0) return
      k = 0
      do j = 1, size(a, 2)
         do i = 1, size(a, 1)
            k = mod(k, size(b)) + 1
            z(i, j) = arithmetic(a(i, j), b(k), op, commute, na_rm)
         end do
      end do
   end function x_op_y

   pure function t_tx_op_y(x, y, op, xrows, xcols, yidxs, commute, na_rm) result(z)
      real(dp), intent(in) :: x(:, :) !! Matrix operand, returned in its selected original orientation.
      real(dp), intent(in) :: y(:) !! Vector operand recycled in row-major order across the selected matrix.
      character(*), intent(in) :: op !! Arithmetic operation: +, -, *, or /.
      integer, intent(in), optional :: xrows(:) !! Optional R-style matrix row selection.
      integer, intent(in), optional :: xcols(:) !! Optional R-style matrix column selection.
      integer, intent(in), optional :: yidxs(:) !! Optional R-style vector selection.
      logical, intent(in), optional :: commute !! Reverse operand order when true.
      logical, intent(in), optional :: na_rm !! Missing-operand removal for addition or multiplication only.
      real(dp), allocatable :: z(:, :)
      z = transpose(x_op_y(transpose(x), y, op, xcols, xrows, yidxs, commute, na_rm))
   end function t_tx_op_y

   pure real(dp) function arithmetic(x, y, op, commute, na_rm) result(z)
      real(dp), intent(in) :: x !! First scalar operand.
      real(dp), intent(in) :: y !! Second scalar operand.
      character(*), intent(in) :: op !! One arithmetic operator: +, -, *, or /.
      logical, intent(in), optional :: commute !! Reverse operands when true.
      logical, intent(in), optional :: na_rm !! Remove a missing operand for + or * when true.
      real(dp) :: a, b
      a = x
      b = y
      if (ms_remove(commute)) then
         a = y
         b = x
      end if
      z = ms_nan()
      if ((op == '+' .or. op == '*') .and. ms_remove(na_rm)) then
         if (ieee_is_nan(a)) then
            z = b
            return
         else if (ieee_is_nan(b)) then
            z = a
            return
         end if
      end if
      select case (op)
      case ('+')
         z = a + b
      case ('-')
         z = a - b
      case ('*')
         z = a * b
      case ('/')
         z = a / b
      end select
   end function arithmetic

end module matrixstats
