! SPDX-License-Identifier: MIT
! SPDX-FileComment: Level algorithms translated from the R forcats package.
module forcats_levels
   !! Implements low-level operations on factor levels and codes.
   use forcats_types, only: factor_type, new_factor
   implicit none
   private

   public :: lvls_expand, lvls_reorder, lvls_revalue, lvls_union

contains

   pure function lvls_reorder(factor, indices, ordered) result(output)
      !! Reorders levels while preserving each observation's decoded value.
      type(factor_type), intent(in) :: factor   !! Factor to transform.
      integer, intent(in) :: indices(:)        !! Permutation of existing level indices.
      logical, intent(in), optional :: ordered !! Replacement ordered status.
      type(factor_type) :: output
      character(len=:), allocatable :: reordered_levels(:)
      integer, allocatable :: codes(:)
      logical :: output_ordered
      integer :: i, old_code

      if (size(indices) /= factor%nlevels()) error stop "lvls_reorder: wrong number of indices"
      do i = 1, size(indices)
         if (count(indices == i) /= 1) error stop "lvls_reorder: indices must be a permutation"
      end do

      codes = factor%codes
      do i = 1, factor%size()
         if (factor%missing(i)) cycle
         old_code = factor%codes(i)
         codes(i) = find_integer(indices, old_code)
      end do
      output_ordered = factor%ordered
      if (present(ordered)) output_ordered = ordered
      reordered_levels = factor%levels(indices)
      output = new_factor(codes, reordered_levels, factor%missing, output_ordered)
   end function lvls_reorder

   pure function lvls_revalue(factor, new_levels) result(output)
      !! Relabels levels and collapses duplicated replacement labels.
      type(factor_type), intent(in) :: factor    !! Factor to transform.
      character(len=*), intent(in) :: new_levels(:) !! One replacement per level.
      type(factor_type) :: output
      character(len=len(new_levels)), allocatable :: unique_levels(:)
      integer, allocatable :: mapping(:), codes(:)
      integer :: i, j, number_unique

      if (size(new_levels) /= factor%nlevels()) error stop "lvls_revalue: wrong number of levels"
      allocate (unique_levels(size(new_levels)), mapping(size(new_levels)))
      number_unique = 0
      do i = 1, size(new_levels)
         j = find_character(unique_levels(:number_unique), new_levels(i))
         if (j == 0) then
            number_unique = number_unique + 1
            unique_levels(number_unique) = new_levels(i)
            j = number_unique
         end if
         mapping(i) = j
      end do

      codes = factor%codes
      do i = 1, factor%size()
         if (.not. factor%missing(i)) codes(i) = mapping(factor%codes(i))
      end do
      output = new_factor(codes, unique_levels(:number_unique), factor%missing, factor%ordered)
   end function lvls_revalue

   pure function lvls_expand(factor, new_levels) result(output)
      !! Replaces the level set with a superset and preserves decoded values.
      type(factor_type), intent(in) :: factor    !! Factor to transform.
      character(len=*), intent(in) :: new_levels(:) !! Superset of existing levels.
      type(factor_type) :: output
      integer, allocatable :: codes(:), mapping(:)
      integer :: i

      allocate (mapping(factor%nlevels()))
      do i = 1, factor%nlevels()
         mapping(i) = find_character(new_levels, factor%levels(i))
         if (mapping(i) == 0) error stop "lvls_expand: new levels omit an existing level"
      end do
      if (has_duplicates(new_levels)) error stop "lvls_expand: new levels must be unique"

      codes = factor%codes
      do i = 1, factor%size()
         if (.not. factor%missing(i)) codes(i) = mapping(factor%codes(i))
      end do
      output = new_factor(codes, new_levels, factor%missing, factor%ordered)
   end function lvls_expand

   pure function lvls_union(factors) result(levels)
      !! Returns the stable union of levels from a factor array.
      type(factor_type), intent(in) :: factors(:) !! Factors whose levels are combined.
      character(len=:), allocatable :: levels(:)
      character(len=:), allocatable :: workspace(:)
      integer :: i, j, maximum_width, number_of_levels

      maximum_width = 1
      number_of_levels = 0
      do i = 1, size(factors)
         if (allocated(factors(i)%levels)) then
            maximum_width = max(maximum_width, len(factors(i)%levels))
            number_of_levels = number_of_levels + factors(i)%nlevels()
         end if
      end do
      allocate (character(len=maximum_width) :: workspace(number_of_levels))
      number_of_levels = 0
      do i = 1, size(factors)
         do j = 1, factors(i)%nlevels()
            if (find_character(workspace(:number_of_levels), factors(i)%levels(j)) == 0) then
               number_of_levels = number_of_levels + 1
               workspace(number_of_levels) = factors(i)%levels(j)
            end if
         end do
      end do
      levels = workspace(:number_of_levels)
   end function lvls_union

   pure integer function find_integer(values, target) result(index)
      !! Finds an integer in an array, returning zero when absent.
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

   pure integer function find_character(values, target) result(index)
      !! Finds a character value in an array, returning zero when absent.
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

   pure logical function has_duplicates(values) result(found)
      !! Reports whether a character array contains duplicate values.
      character(len=*), intent(in) :: values(:) !! Values to inspect.
      integer :: i

      found = .false.
      do i = 2, size(values)
         if (any(values(:i - 1) == values(i))) then
            found = .true.
            return
         end if
      end do
   end function has_duplicates

end module forcats_levels
