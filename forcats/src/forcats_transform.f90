! SPDX-License-Identifier: MIT
! SPDX-FileComment: Factor transformations translated from the R forcats package.
module forcats_transform
   !! Implements deterministic factor reordering, recoding, and combining.
   use forcats_levels, only: lvls_expand, lvls_reorder, lvls_revalue, lvls_union
   use forcats_types, only: dp, factor_type, fct, new_factor
   implicit none
   private

   public :: fct_c, fct_collapse, fct_cross, fct_drop, fct_expand, fct_explicit_na
   public :: fct_infreq, fct_inorder, fct_inseq, fct_match, fct_other
   public :: fct_na_level_to_value, fct_na_value_to_level, fct_recode
   public :: fct_relevel, fct_rev, fct_shift, fct_unify, fct_unique

contains

   pure elemental function fct_explicit_na(factor, na_level) result(output)
      !! Provides the deprecated upstream alias for converting missing values to a level.
      type(factor_type), intent(in) :: factor       !! Factor to transform.
      character(len=*), intent(in), optional :: na_level !! Missing label; defaults to (Missing).
      type(factor_type) :: output

      if (present(na_level)) then
         output = fct_na_value_to_level(factor, na_level)
      else
         output = fct_na_value_to_level(factor, "(Missing)")
      end if
   end function fct_explicit_na

   pure function fct_relevel(factor, first_levels, after) result(output)
      !! Moves named levels to a requested position while retaining other order.
      type(factor_type), intent(in) :: factor          !! Factor to reorder.
      character(len=*), intent(in) :: first_levels(:) !! Levels to move.
      integer, intent(in), optional :: after           !! Number of retained levels preceding moved levels.
      type(factor_type) :: output
      integer, allocatable :: selected(:), retained(:), permutation(:)
      integer :: i, insertion, number_retained, number_selected

      allocate (selected(factor%nlevels()), retained(factor%nlevels()))
      number_selected = 0
      do i = 1, size(first_levels)
         if (find_character(factor%levels, first_levels(i)) == 0) cycle
         if (find_integer(selected(:number_selected), find_character(factor%levels, first_levels(i))) > 0) cycle
         number_selected = number_selected + 1
         selected(number_selected) = find_character(factor%levels, first_levels(i))
      end do
      number_retained = 0
      do i = 1, factor%nlevels()
         if (find_integer(selected(:number_selected), i) == 0) then
            number_retained = number_retained + 1
            retained(number_retained) = i
         end if
      end do
      insertion = 0
      if (present(after)) insertion = max(0, min(after, number_retained))
      permutation = [retained(:insertion), selected(:number_selected), retained(insertion + 1:number_retained)]
      output = lvls_reorder(factor, permutation)
   end function fct_relevel

   pure elemental function fct_inorder(factor, ordered) result(output)
      !! Reorders levels by first nonmissing appearance, then unused levels.
      type(factor_type), intent(in) :: factor   !! Factor to reorder.
      logical, intent(in), optional :: ordered !! Replacement ordered status.
      type(factor_type) :: output
      integer, allocatable :: permutation(:)
      integer :: i, n

      allocate (permutation(factor%nlevels()))
      n = 0
      do i = 1, factor%size()
         if (factor%missing(i)) cycle
         if (find_integer(permutation(:n), factor%codes(i)) == 0) then
            n = n + 1
            permutation(n) = factor%codes(i)
         end if
      end do
      do i = 1, factor%nlevels()
         if (find_integer(permutation(:n), i) == 0) then
            n = n + 1
            permutation(n) = i
         end if
      end do
      output = lvls_reorder(factor, permutation, ordered)
   end function fct_inorder

   pure function fct_infreq(factor, weights, ordered) result(output)
      !! Reorders levels by decreasing weighted frequency using stable ties.
      type(factor_type), intent(in) :: factor    !! Factor to reorder.
      real(dp), intent(in), optional :: weights(:) !! Nonnegative observation weights.
      logical, intent(in), optional :: ordered   !! Replacement ordered status.
      type(factor_type) :: output
      integer, allocatable :: permutation(:)
      real(dp), allocatable :: frequencies(:)

      frequencies = level_weights(factor, weights)
      permutation = stable_order(frequencies, descending=.true.)
      output = lvls_reorder(factor, permutation, ordered)
   end function fct_infreq

   pure elemental function fct_inseq(factor, ordered) result(output)
      !! Reorders levels by their numeric interpretation.
      type(factor_type), intent(in) :: factor   !! Factor to reorder.
      logical, intent(in), optional :: ordered !! Replacement ordered status.
      type(factor_type) :: output
      integer, allocatable :: permutation(:)
      real(dp), allocatable :: numeric_levels(:)
      integer :: i, status

      allocate (numeric_levels(factor%nlevels()))
      do i = 1, factor%nlevels()
         read (factor%levels(i), *, iostat=status) numeric_levels(i)
         if (status /= 0) error stop "fct_inseq: every level must be numeric"
      end do
      permutation = stable_order(numeric_levels)
      output = lvls_reorder(factor, permutation, ordered)
   end function fct_inseq

   pure elemental function fct_rev(factor) result(output)
      !! Reverses factor level order.
      type(factor_type), intent(in) :: factor !! Factor to reorder.
      type(factor_type) :: output
      integer, allocatable :: permutation(:)
      integer :: i

      allocate (permutation(factor%nlevels()))
      do i = 1, factor%nlevels()
         permutation(i) = factor%nlevels() - i + 1
      end do
      output = lvls_reorder(factor, permutation)
   end function fct_rev

   pure elemental function fct_shift(factor, n) result(output)
      !! Cyclically shifts levels left for positive n and right for negative n.
      type(factor_type), intent(in) :: factor !! Factor to reorder.
      integer, intent(in), optional :: n      !! Signed shift; defaults to one.
      type(factor_type) :: output
      integer, allocatable :: permutation(:)
      integer :: i, shift

      if (factor%nlevels() == 0) then
         output = factor
         return
      end if
      shift = 1
      if (present(n)) shift = n
      allocate (permutation(factor%nlevels()))
      do i = 1, factor%nlevels()
         permutation(i) = modulo(i - 1 + shift, factor%nlevels()) + 1
      end do
      output = lvls_reorder(factor, permutation)
   end function fct_shift

   pure function fct_expand(factor, additional_levels, after) result(output)
      !! Adds previously absent levels at a requested position.
      type(factor_type), intent(in) :: factor              !! Factor to expand.
      character(len=*), intent(in) :: additional_levels(:) !! Candidate levels to add.
      integer, intent(in), optional :: after               !! Existing levels preceding additions.
      type(factor_type) :: output
      character(len=:), allocatable :: combined(:), additions(:)
      integer :: i, insertion, n

      allocate (character(len=max(len(factor%levels), len(additional_levels))) :: additions(size(additional_levels)))
      n = 0
      do i = 1, size(additional_levels)
         if (find_character(factor%levels, additional_levels(i)) > 0) cycle
         if (find_character(additions(:n), additional_levels(i)) > 0) cycle
         n = n + 1
         additions(n) = additional_levels(i)
      end do
      insertion = factor%nlevels()
      if (present(after)) insertion = max(0, min(after, factor%nlevels()))
      allocate (character(len=max(len(factor%levels), len(additions))) :: combined(factor%nlevels() + n))
      combined(:insertion) = factor%levels(:insertion)
      combined(insertion + 1:insertion + n) = additions(:n)
      combined(insertion + n + 1:) = factor%levels(insertion + 1:)
      output = lvls_expand(factor, combined)
   end function fct_expand

   pure function fct_drop(factor, only) result(output)
      !! Drops unused levels, optionally restricting which labels may be removed.
      type(factor_type), intent(in) :: factor    !! Factor to compact.
      character(len=*), intent(in), optional :: only(:) !! Eligible unused levels.
      type(factor_type) :: output
      character(len=len(factor%levels)), allocatable :: retained(:)
      integer, allocatable :: mapping(:), codes(:)
      logical :: eligible
      integer :: i, n

      allocate (retained(factor%nlevels()), mapping(factor%nlevels()))
      n = 0
      do i = 1, factor%nlevels()
         eligible = .true.
         if (present(only)) eligible = find_character(only, factor%levels(i)) > 0
         if (eligible .and. .not. level_is_used(factor, i)) then
            mapping(i) = 0
         else
            n = n + 1
            retained(n) = factor%levels(i)
            mapping(i) = n
         end if
      end do
      codes = factor%codes
      do i = 1, factor%size()
         if (.not. factor%missing(i)) codes(i) = mapping(factor%codes(i))
      end do
      output = new_factor(codes, retained(:n), factor%missing, factor%ordered)
   end function fct_drop

   pure function fct_recode(factor, old_levels, new_levels) result(output)
      !! Renames selected old levels and collapses equal replacements.
      type(factor_type), intent(in) :: factor       !! Factor to recode.
      character(len=*), intent(in) :: old_levels(:) !! Existing labels to replace.
      character(len=*), intent(in) :: new_levels(:) !! Replacement labels.
      type(factor_type) :: output
      character(len=:), allocatable :: replacements(:)
      integer :: i, index

      if (size(old_levels) /= size(new_levels)) error stop "fct_recode: mapping arrays differ in size"
      allocate (character(len=max(len(factor%levels), len(new_levels))) :: replacements(factor%nlevels()))
      replacements(:) = factor%levels
      do i = 1, size(old_levels)
         index = find_character(factor%levels, old_levels(i))
         if (index > 0) replacements(index) = new_levels(i)
      end do
      output = lvls_revalue(factor, replacements)
   end function fct_recode

   pure function fct_collapse(factor, old_levels, group_levels, other_level) result(output)
      !! Collapses selected levels into named groups and optionally groups all others.
      type(factor_type), intent(in) :: factor         !! Factor to collapse.
      character(len=*), intent(in) :: old_levels(:)   !! Existing labels to map.
      character(len=*), intent(in) :: group_levels(:) !! Group label for each old label.
      character(len=*), intent(in), optional :: other_level !! Optional label for unmapped levels.
      type(factor_type) :: output
      logical, allocatable :: mapped(:)
      integer :: i, index

      if (size(old_levels) /= size(group_levels)) error stop "fct_collapse: mapping arrays differ in size"
      output = fct_recode(factor, old_levels, group_levels)
      if (.not. present(other_level)) return
      allocate (mapped(output%nlevels()), source=.false.)
      do i = 1, size(group_levels)
         index = find_character(output%levels, group_levels(i))
         if (index > 0) mapped(index) = .true.
      end do
      output = other_by_mask(output, mapped, other_level)
   end function fct_collapse

   pure function fct_other(factor, keep, other_level) result(output)
      !! Keeps selected labels and replaces every other level with one label.
      type(factor_type), intent(in) :: factor   !! Factor to transform.
      character(len=*), intent(in) :: keep(:)  !! Labels to preserve.
      character(len=*), intent(in), optional :: other_level !! Replacement label.
      type(factor_type) :: output
      logical, allocatable :: keep_mask(:)
      character(len=:), allocatable :: replacement
      integer :: i

      replacement = "Other"
      if (present(other_level)) replacement = other_level
      allocate (keep_mask(factor%nlevels()))
      do i = 1, factor%nlevels()
         keep_mask(i) = find_character(keep, factor%levels(i)) > 0
      end do
      output = other_by_mask(factor, keep_mask, replacement)
   end function fct_other

   pure elemental function fct_na_value_to_level(factor, level) result(output)
      !! Converts missing observations to an explicit factor level.
      type(factor_type), intent(in) :: factor !! Factor to transform.
      character(len=*), intent(in), optional :: level !! Missing-level label.
      type(factor_type) :: output
      character(len=:), allocatable :: replacement, expanded(:)
      integer, allocatable :: codes(:)
      logical, allocatable :: missing(:)
      integer :: index

      replacement = "NA"
      if (present(level)) replacement = level
      index = find_character(factor%levels, replacement)
      if (index == 0) then
         allocate (character(len=max(len(factor%levels), len(replacement))) :: expanded(factor%nlevels() + 1))
         expanded(:factor%nlevels()) = factor%levels
         expanded(factor%nlevels() + 1) = replacement
         index = size(expanded)
      else
         expanded = factor%levels
      end if
      codes = factor%codes
      where (factor%missing) codes = index
      allocate (missing(factor%size()), source=.false.)
      output = new_factor(codes, expanded, missing, factor%ordered)
   end function fct_na_value_to_level

   pure function fct_na_level_to_value(factor, levels) result(output)
      !! Converts selected explicit levels to missing observations and removes them.
      type(factor_type), intent(in) :: factor       !! Factor to transform.
      character(len=*), intent(in), optional :: levels(:) !! Labels treated as missing; defaults to NA.
      type(factor_type) :: output
      character(len=:), allocatable :: targets(:)
      character(len=len(factor%levels)), allocatable :: retained(:)
      integer, allocatable :: mapping(:), codes(:)
      logical, allocatable :: missing(:)
      integer :: i, n

      if (present(levels)) then
         targets = levels
      else
         targets = [character(len=2) :: "NA"]
      end if
      allocate (retained(factor%nlevels()), mapping(factor%nlevels()))
      n = 0
      do i = 1, factor%nlevels()
         if (find_character(targets, factor%levels(i)) > 0) then
            mapping(i) = 0
         else
            n = n + 1
            retained(n) = factor%levels(i)
            mapping(i) = n
         end if
      end do
      codes = factor%codes
      missing = factor%missing
      do i = 1, factor%size()
         if (missing(i)) cycle
         if (mapping(factor%codes(i)) == 0) then
            missing(i) = .true.
            codes(i) = 0
         else
            codes(i) = mapping(factor%codes(i))
         end if
      end do
      output = new_factor(codes, retained(:n), missing, factor%ordered)
   end function fct_na_level_to_value

   pure function fct_match(factor, levels) result(matches)
      !! Tests observations against validated factor levels.
      type(factor_type), intent(in) :: factor    !! Factor to test.
      character(len=*), intent(in) :: levels(:) !! Existing labels to match.
      logical, allocatable :: matches(:)
      integer :: i

      do i = 1, size(levels)
         if (find_character(factor%levels, levels(i)) == 0) error stop "fct_match: unknown level"
      end do
      allocate (matches(factor%size()), source=.false.)
      do i = 1, factor%size()
         if (.not. factor%missing(i)) matches(i) = find_character(levels, factor%levels(factor%codes(i))) > 0
      end do
   end function fct_match

   pure elemental function fct_unique(factor) result(output)
      !! Returns one value for every possible level and one missing value if present.
      type(factor_type), intent(in) :: factor !! Factor whose possibilities are returned.
      type(factor_type) :: output
      integer, allocatable :: codes(:)
      logical, allocatable :: missing(:)
      integer :: i, n

      n = factor%nlevels() + merge(1, 0, any(factor%missing))
      allocate (codes(n), source=0)
      allocate (missing(n), source=.false.)
      do i = 1, factor%nlevels()
         codes(i) = i
      end do
      if (n > factor%nlevels()) missing(n) = .true.
      output = new_factor(codes, factor%levels, missing, factor%ordered)
   end function fct_unique

   pure function fct_c(factors) result(output)
      !! Concatenates factors after taking the stable union of their levels.
      type(factor_type), intent(in) :: factors(:) !! Factors to concatenate.
      type(factor_type) :: output
      character(len=:), allocatable :: levels(:)
      integer, allocatable :: codes(:)
      logical, allocatable :: missing(:)
      integer :: i, j, offset, total

      levels = lvls_union(factors)
      total = 0
      do i = 1, size(factors)
         total = total + factors(i)%size()
      end do
      allocate (codes(total), missing(total))
      offset = 0
      do i = 1, size(factors)
         do j = 1, factors(i)%size()
            missing(offset + j) = factors(i)%missing(j)
            if (missing(offset + j)) then
               codes(offset + j) = 0
            else
               codes(offset + j) = find_character(levels, factors(i)%levels(factors(i)%codes(j)))
            end if
         end do
         offset = offset + factors(i)%size()
      end do
      output = new_factor(codes, levels, missing)
   end function fct_c

   pure function fct_unify(factors, levels) result(output)
      !! Gives every factor the same supplied or inferred stable level union.
      type(factor_type), intent(in) :: factors(:) !! Factors to harmonize.
      character(len=*), intent(in), optional :: levels(:) !! Common level set.
      type(factor_type), allocatable :: output(:)
      character(len=:), allocatable :: common_levels(:)
      integer :: i

      if (present(levels)) then
         common_levels = levels
      else
         common_levels = lvls_union(factors)
      end if
      allocate (output(size(factors)))
      do i = 1, size(factors)
         output(i) = lvls_expand(factors(i), common_levels)
      end do
   end function fct_unify

   pure elemental function fct_cross(left, right, separator, keep_empty) result(output)
      !! Combines two factors using Cartesian-order levels and rowwise values.
      type(factor_type), intent(in) :: left, right !! Same-sized factors to combine.
      character(len=*), intent(in), optional :: separator !! Label separator; defaults to colon.
      logical, intent(in), optional :: keep_empty !! Retain unobserved combinations.
      type(factor_type) :: output
      character(len=:), allocatable :: sep, all_levels(:), used_levels(:), row_values(:)
      logical, allocatable :: missing(:), used(:)
      integer :: i, j, k

      if (left%size() /= right%size()) error stop "fct_cross: factors differ in size"
      sep = ":"
      if (present(separator)) sep = separator
      allocate (character(len=len(left%levels) + len(sep) + len(right%levels)) :: &
         all_levels(left%nlevels() * right%nlevels()), row_values(left%size()))
      k = 0
      do i = 1, left%nlevels()
         do j = 1, right%nlevels()
            k = k + 1
            all_levels(k) = trim(left%levels(i)) // sep // trim(right%levels(j))
         end do
      end do
      missing = left%missing .or. right%missing
      row_values(:) = ""
      do i = 1, left%size()
         if (.not. missing(i)) then
            row_values(i) = trim(left%levels(left%codes(i))) // sep // trim(right%levels(right%codes(i)))
         end if
      end do
      if (present(keep_empty)) then
         if (keep_empty) then
            output = fct(row_values, all_levels, missing)
            return
         end if
      end if
      allocate (used(size(all_levels)), source=.false.)
      do i = 1, size(all_levels)
         used(i) = any(.not. missing .and. row_values == all_levels(i))
      end do
      used_levels = pack(all_levels, used)
      output = fct(row_values, used_levels, missing)
   end function fct_cross

   pure function other_by_mask(factor, keep, other_level) result(output)
      !! Collapses levels selected by a logical complement into a final level.
      type(factor_type), intent(in) :: factor     !! Factor to transform.
      logical, intent(in) :: keep(:)              !! One preservation flag per level.
      character(len=*), intent(in) :: other_level !! Replacement label.
      type(factor_type) :: output
      character(len=:), allocatable :: replacements(:)

      if (size(keep) /= factor%nlevels()) error stop "fct_other: keep mask has the wrong size"
      if (all(keep)) then
         output = factor
         return
      end if
      allocate (character(len=max(len(factor%levels), len(other_level))) :: replacements(factor%nlevels()))
      replacements(:) = factor%levels
      where (.not. keep) replacements = other_level
      output = lvls_revalue(factor, replacements)
      output = fct_relevel(output, [other_level], after=output%nlevels())
   end function other_by_mask

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

   pure function stable_order(values, descending) result(indices)
      !! Returns stable ascending or descending order indices.
      real(dp), intent(in) :: values(:)          !! Values to rank.
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
   end function stable_order

   pure logical function level_is_used(factor, level) result(used)
      !! Reports whether a level has at least one nonmissing observation.
      type(factor_type), intent(in) :: factor !! Factor to inspect.
      integer, intent(in) :: level           !! One-based level index.

      used = any(.not. factor%missing .and. factor%codes == level)
   end function level_is_used

   pure integer function find_character(values, target) result(index)
      !! Finds a string, returning zero when absent.
      character(len=*), intent(in) :: values(:) !! Values to search.
      character(len=*), intent(in) :: target    !! Value to find.
      integer :: i

      index = 0
      do i = 1, size(values)
         if (values(i) == target) then
            index = i
            return
         end if
      end do
   end function find_character

   pure integer function find_integer(values, target) result(index)
      !! Finds an integer, returning zero when absent.
      integer, intent(in) :: values(:) !! Values to search.
      integer, intent(in) :: target    !! Value to find.
      integer :: i

      index = 0
      do i = 1, size(values)
         if (values(i) == target) then
            index = i
            return
         end if
      end do
   end function find_integer

end module forcats_transform
