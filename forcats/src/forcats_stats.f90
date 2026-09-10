! SPDX-License-Identifier: MIT
! SPDX-FileComment: Factor statistics translated from the R forcats package.
module forcats_stats
   !! Implements factor counting, lumping, and data-driven reordering.
   use forcats_levels, only: lvls_reorder
   use forcats_transform, only: fct_other, fct_relevel, fct_unique
   use forcats_types, only: dp, factor_count_type, factor_type, new_factor
   implicit none
   private

   public :: fct_count, fct_lump, fct_lump_lowfreq, fct_lump_min
   public :: fct_lump_n, fct_lump_prop, fct_reorder, fct_reorder2
   public :: first2, last2

contains

   pure elemental function fct_count(factor, sort, prop) result(table)
      !! Counts observations for every level and for missing values when present.
      type(factor_type), intent(in) :: factor !! Factor to count.
      logical, intent(in), optional :: sort   !! Sort rows by decreasing count.
      logical, intent(in), optional :: prop   !! Compute marginal proportions.
      type(factor_count_type) :: table
      integer, allocatable :: counts(:), permutation(:)
      logical :: include_proportion, sort_rows
      integer :: i, n, total

      n = factor%nlevels() + merge(1, 0, any(factor%missing))
      allocate (counts(n), source=0)
      do i = 1, factor%size()
         if (factor%missing(i)) then
            counts(n) = counts(n) + 1
         else
            counts(factor%codes(i)) = counts(factor%codes(i)) + 1
         end if
      end do
      table%factor = fct_unique(factor)
      sort_rows = .false.
      if (present(sort)) sort_rows = sort
      if (sort_rows) then
         permutation = stable_count_order(counts)
         table%factor = reorder_observations(table%factor, permutation)
         table%count = counts(permutation)
      else
         table%count = counts
      end if
      include_proportion = .false.
      if (present(prop)) include_proportion = prop
      if (include_proportion) then
         total = sum(table%count)
         allocate (table%proportion(n), source=0.0_dp)
         if (total > 0) table%proportion = real(table%count, dp) / real(total, dp)
      end if
   end function fct_count

   pure function fct_lump(factor, n, prop, weights, other_level) result(output)
      !! Applies the historical lump dispatcher using n, prop, or low frequency.
      type(factor_type), intent(in) :: factor !! Factor to lump.
      integer, intent(in), optional :: n      !! Number of frequent levels to retain.
      real(dp), intent(in), optional :: prop  !! Proportion threshold.
      real(dp), intent(in), optional :: weights(:) !! Observation weights.
      character(len=*), intent(in), optional :: other_level !! Collapsed label.
      type(factor_type) :: output

      if (present(n) .and. present(prop)) error stop "fct_lump: supply n or prop, not both"
      if (present(n)) then
         output = fct_lump_n(factor, n, weights, other_level)
      else if (present(prop)) then
         output = fct_lump_prop(factor, prop, weights, other_level)
      else
         output = fct_lump_lowfreq(factor, weights, other_level)
      end if
   end function fct_lump

   pure function fct_lump_min(factor, minimum, weights, other_level) result(output)
      !! Preserves levels whose weighted frequency is at least a threshold.
      type(factor_type), intent(in) :: factor !! Factor to lump.
      real(dp), intent(in) :: minimum         !! Minimum retained frequency.
      real(dp), intent(in), optional :: weights(:) !! Observation weights.
      character(len=*), intent(in), optional :: other_level !! Collapsed label.
      type(factor_type) :: output
      character(len=:), allocatable :: keep(:), replacement
      real(dp), allocatable :: frequencies(:)

      if (minimum < 0.0_dp) error stop "fct_lump_min: minimum must be nonnegative"
      frequencies = level_weights(factor, weights)
      keep = pack(factor%levels, frequencies >= minimum)
      replacement = "Other"
      if (present(other_level)) replacement = other_level
      output = fct_other(factor, keep, replacement)
   end function fct_lump_min

   pure function fct_lump_prop(factor, proportion, weights, other_level) result(output)
      !! Preserves levels on the selected side of a marginal proportion threshold.
      type(factor_type), intent(in) :: factor !! Factor to lump.
      real(dp), intent(in) :: proportion      !! Signed proportion threshold.
      real(dp), intent(in), optional :: weights(:) !! Observation weights.
      character(len=*), intent(in), optional :: other_level !! Collapsed label.
      type(factor_type) :: output
      character(len=:), allocatable :: keep(:), replacement
      real(dp), allocatable :: frequencies(:), proportions(:)
      real(dp) :: total

      frequencies = level_weights(factor, weights)
      if (present(weights)) then
         total = sum(weights)
      else
         total = real(factor%size(), dp)
      end if
      allocate (proportions(size(frequencies)), source=0.0_dp)
      if (total > 0.0_dp) proportions = frequencies / total
      if (proportion < 0.0_dp) then
         keep = pack(factor%levels, proportions <= -proportion)
      else
         keep = pack(factor%levels, proportions > proportion)
      end if
      replacement = "Other"
      if (present(other_level)) replacement = other_level
      output = fct_other(factor, keep, replacement)
   end function fct_lump_prop

   pure function fct_lump_n(factor, n, weights, other_level) result(output)
      !! Preserves the n most frequent levels, or least frequent when n is negative.
      type(factor_type), intent(in) :: factor !! Factor to lump.
      integer, intent(in) :: n               !! Signed number of levels to retain.
      real(dp), intent(in), optional :: weights(:) !! Observation weights.
      character(len=*), intent(in), optional :: other_level !! Collapsed label.
      type(factor_type) :: output
      character(len=:), allocatable :: keep(:), replacement
      real(dp), allocatable :: frequencies(:)
      logical, allocatable :: retained(:)
      real(dp) :: threshold
      integer :: requested

      frequencies = level_weights(factor, weights)
      requested = min(abs(n), factor%nlevels())
      allocate (retained(factor%nlevels()), source=.false.)
      if (requested > 0) then
         threshold = nth_frequency(frequencies, requested, descending=n >= 0)
         if (n >= 0) then
            retained = frequencies >= threshold
         else
            retained = frequencies <= threshold
         end if
      end if
      keep = pack(factor%levels, retained)
      replacement = "Other"
      if (present(other_level)) replacement = other_level
      output = fct_other(factor, keep, replacement)
   end function fct_lump_n

   pure function fct_lump_lowfreq(factor, weights, other_level) result(output)
      !! Lumps the smallest groups while keeping their collective group smallest.
      type(factor_type), intent(in) :: factor !! Factor to lump.
      real(dp), intent(in), optional :: weights(:) !! Observation weights.
      character(len=*), intent(in), optional :: other_level !! Collapsed label.
      type(factor_type) :: output
      character(len=:), allocatable :: keep(:), replacement
      integer, allocatable :: order(:)
      logical, allocatable :: retain(:)
      real(dp), allocatable :: frequencies(:)
      real(dp) :: remaining
      integer :: i, cutoff

      frequencies = level_weights(factor, weights)
      order = stable_real_order(frequencies, descending=.true.)
      remaining = sum(frequencies)
      cutoff = factor%nlevels() + 1
      do i = 1, factor%nlevels()
         remaining = remaining - frequencies(order(i))
         if (frequencies(order(i)) > remaining) then
            cutoff = i + 1
            exit
         end if
      end do
      allocate (retain(factor%nlevels()), source=.true.)
      do i = cutoff, factor%nlevels()
         retain(order(i)) = .false.
      end do
      keep = pack(factor%levels, retain)
      replacement = "Other"
      if (present(other_level)) replacement = other_level
      output = fct_other(factor, keep, replacement)
   end function fct_lump_lowfreq

   pure function fct_reorder(factor, x, statistic, descending, missing_x, default) result(output)
      !! Reorders levels using a built-in scalar summary of an aligned real vector.
      type(factor_type), intent(in) :: factor !! Factor defining groups.
      real(dp), intent(in) :: x(:)            !! Values summarized within levels.
      character(len=*), intent(in), optional :: statistic !! mean, median, min, or max.
      logical, intent(in), optional :: descending !! Sort summaries in descending order.
      logical, intent(in), optional :: missing_x(:) !! Missing-value mask for x.
      real(dp), intent(in), optional :: default !! Summary for empty levels.
      type(factor_type) :: output
      character(len=:), allocatable :: operation
      integer, allocatable :: permutation(:)
      real(dp), allocatable :: summaries(:), values(:)
      logical, allocatable :: mask(:)
      logical :: reverse
      integer :: level

      if (size(x) /= factor%size()) error stop "fct_reorder: x has the wrong size"
      allocate (mask(size(x)), source=.false.)
      if (present(missing_x)) then
         if (size(missing_x) /= size(x)) error stop "fct_reorder: missing mask has the wrong size"
         mask = missing_x
      end if
      operation = "median"
      if (present(statistic)) operation = trim(statistic)
      allocate (summaries(factor%nlevels()))
      summaries = huge(0.0_dp)
      if (present(default)) summaries = default
      do level = 1, factor%nlevels()
         values = pack(x, .not. factor%missing .and. .not. mask .and. factor%codes == level)
         if (size(values) > 0) summaries(level) = summarize(values, operation)
      end do
      reverse = .false.
      if (present(descending)) reverse = descending
      permutation = stable_real_order(summaries, reverse)
      output = lvls_reorder(factor, permutation)
   end function fct_reorder

   pure function fct_reorder2(factor, x, y, last, descending, missing_x, missing_y, default) result(output)
      !! Reorders levels by the first or last y after ordering aligned x values.
      type(factor_type), intent(in) :: factor !! Factor defining groups.
      real(dp), intent(in) :: x(:)            !! Within-group ordering variable.
      real(dp), intent(in) :: y(:)            !! Within-group terminal value.
      logical, intent(in), optional :: last   !! Use last terminal value; defaults true.
      logical, intent(in), optional :: descending !! Sort summaries descending; defaults true.
      logical, intent(in), optional :: missing_x(:), missing_y(:) !! Missing masks.
      real(dp), intent(in), optional :: default !! Summary for empty levels.
      type(factor_type) :: output
      integer, allocatable :: permutation(:), locations(:)
      real(dp), allocatable :: summaries(:)
      logical, allocatable :: mask(:)
      logical :: reverse, use_last
      integer :: i, level, selected

      if (size(x) /= factor%size() .or. size(y) /= factor%size()) then
         error stop "fct_reorder2: aligned arrays have the wrong size"
      end if
      allocate (mask(factor%size()), source=.false.)
      if (present(missing_x)) mask = mask .or. missing_x
      if (present(missing_y)) mask = mask .or. missing_y
      allocate (summaries(factor%nlevels()))
      summaries = -huge(0.0_dp)
      if (present(default)) summaries = default
      use_last = .true.
      if (present(last)) use_last = last
      do level = 1, factor%nlevels()
         locations = pack([(i, i=1, factor%size())], &
            .not. factor%missing .and. .not. mask .and. factor%codes == level)
         if (size(locations) == 0) cycle
         selected = locations(1)
         do i = 2, size(locations)
            if (use_last .and. x(locations(i)) > x(selected)) selected = locations(i)
            if (.not. use_last .and. x(locations(i)) < x(selected)) selected = locations(i)
         end do
         summaries(level) = y(selected)
      end do
      reverse = .true.
      if (present(descending)) reverse = descending
      permutation = stable_real_order(summaries, reverse)
      output = lvls_reorder(factor, permutation)
   end function fct_reorder2

   pure function last2(x, y) result(value)
      !! Returns y at the largest x, retaining first occurrence on ties.
      real(dp), intent(in) :: x(:) !! Ordering values.
      real(dp), intent(in) :: y(:) !! Values from which to select.
      real(dp) :: value

      if (size(x) /= size(y) .or. size(x) == 0) error stop "last2: arrays must be equal and nonempty"
      value = y(maxloc(x, dim=1))
   end function last2

   pure function first2(x, y) result(value)
      !! Returns y at the smallest x, retaining first occurrence on ties.
      real(dp), intent(in) :: x(:) !! Ordering values.
      real(dp), intent(in) :: y(:) !! Values from which to select.
      real(dp) :: value

      if (size(x) /= size(y) .or. size(x) == 0) error stop "first2: arrays must be equal and nonempty"
      value = y(minloc(x, dim=1))
   end function first2

   pure function level_weights(factor, weights) result(frequencies)
      !! Computes nonnegative observation weights for each level.
      type(factor_type), intent(in) :: factor    !! Factor defining groups.
      real(dp), intent(in), optional :: weights(:) !! Observation weights.
      real(dp), allocatable :: frequencies(:)
      integer :: i

      allocate (frequencies(factor%nlevels()), source=0.0_dp)
      if (present(weights)) then
         if (size(weights) /= factor%size()) error stop "factor weights have the wrong size"
         if (any(weights < 0.0_dp)) error stop "factor weights must be nonnegative"
         do i = 1, factor%size()
            if (.not. factor%missing(i)) frequencies(factor%codes(i)) = frequencies(factor%codes(i)) + weights(i)
         end do
      else
         do i = 1, factor%size()
            if (.not. factor%missing(i)) frequencies(factor%codes(i)) = frequencies(factor%codes(i)) + 1.0_dp
         end do
      end if
   end function level_weights

   pure real(dp) function nth_frequency(values, n, descending) result(value)
      !! Returns the nth sorted frequency.
      real(dp), intent(in) :: values(:) !! Frequencies to order.
      integer, intent(in) :: n          !! One-based rank.
      logical, intent(in) :: descending !! Sort direction.
      integer, allocatable :: order(:)

      order = stable_real_order(values, descending)
      value = values(order(n))
   end function nth_frequency

   pure function stable_real_order(values, descending) result(indices)
      !! Returns stable ascending or descending order indices.
      real(dp), intent(in) :: values(:) !! Values to order.
      logical, intent(in), optional :: descending !! Reverse ordering when true.
      integer, allocatable :: indices(:)
      logical :: reverse
      integer :: i, j, index

      reverse = .false.
      if (present(descending)) reverse = descending
      allocate (indices(size(values)))
      do i = 1, size(values)
         indices(i) = i
      end do
      do i = 2, size(indices)
         index = indices(i)
         j = i - 1
         do while (j >= 1)
            if (reverse) then
               if (values(indices(j)) >= values(index)) exit
            else
               if (values(indices(j)) <= values(index)) exit
            end if
            indices(j + 1) = indices(j)
            j = j - 1
         end do
         indices(j + 1) = index
      end do
   end function stable_real_order

   pure function stable_count_order(values) result(indices)
      !! Returns stable decreasing order indices for integer counts.
      integer, intent(in) :: values(:) !! Counts to order.
      integer, allocatable :: indices(:)
      integer :: i, j, index

      allocate (indices(size(values)))
      do i = 1, size(values)
         indices(i) = i
      end do
      do i = 2, size(indices)
         index = indices(i)
         j = i - 1
         do while (j >= 1)
            if (values(indices(j)) >= values(index)) exit
            indices(j + 1) = indices(j)
            j = j - 1
         end do
         indices(j + 1) = index
      end do
   end function stable_count_order

   pure function reorder_observations(factor, permutation) result(output)
      !! Reorders observations without changing the factor level order.
      type(factor_type), intent(in) :: factor !! Factor to slice.
      integer, intent(in) :: permutation(:)   !! Observation permutation.
      type(factor_type) :: output

      output = new_factor(factor%codes(permutation), factor%levels, factor%missing(permutation), factor%ordered)
   end function reorder_observations

   pure real(dp) function summarize(values, statistic) result(value)
      !! Applies one supported scalar summary to a nonempty array.
      real(dp), intent(in) :: values(:)       !! Values to summarize.
      character(len=*), intent(in) :: statistic !! Summary name.
      real(dp), allocatable :: sorted(:)
      integer :: middle

      select case (statistic)
      case ("mean")
         value = sum(values) / real(size(values), dp)
      case ("min")
         value = minval(values)
      case ("max")
         value = maxval(values)
      case ("median")
         sorted = values(stable_real_order(values))
         middle = size(sorted) / 2
         if (mod(size(sorted), 2) == 0) then
            value = 0.5_dp * (sorted(middle) + sorted(middle + 1))
         else
            value = sorted(middle + 1)
         end if
      case default
         error stop "fct_reorder: unsupported statistic"
      end select
   end function summarize

end module forcats_stats
