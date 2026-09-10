! SPDX-License-Identifier: MIT
! SPDX-FileComment: Categorical container inspired by the R forcats package.
module forcats_types
   !! Defines the categorical value and count-table types used by forcats.
   use, intrinsic :: ieee_arithmetic, only: ieee_is_nan
   use, intrinsic :: iso_fortran_env, only: real64
   implicit none
   private

   integer, parameter, public :: dp = real64

   type, public :: factor_type
      !! Stores factor codes, ordered levels, missingness, and ordered status.
      integer, allocatable :: codes(:)
      character(len=:), allocatable :: levels(:)
      logical, allocatable :: missing(:)
      logical :: ordered = .false.
   contains
      procedure :: size => factor_size
      procedure :: nlevels => factor_nlevels
      procedure :: values => factor_values
      procedure :: valid => factor_valid
   end type factor_type

   type, public :: factor_count_type
      !! Stores the result of counting every level and any missing values.
      type(factor_type) :: factor
      integer, allocatable :: count(:)
      real(dp), allocatable :: proportion(:)
   end type factor_count_type

   public :: as_factor, fct, new_factor

   interface as_factor
      module procedure as_factor_character, as_factor_integer
      module procedure as_factor_real, as_factor_logical
   end interface as_factor

   interface fct
      module procedure fct_character
   end interface fct

contains

   pure function new_factor(codes, levels, missing, ordered) result(factor)
      !! Constructs a factor from validated one-based codes and explicit levels.
      integer, intent(in) :: codes(:)             !! One-based level codes.
      character(len=*), intent(in) :: levels(:)   !! Ordered level labels.
      logical, intent(in), optional :: missing(:) !! Missing-value mask.
      logical, intent(in), optional :: ordered    !! Whether levels are ordinal.
      type(factor_type) :: factor

      if (present(missing)) then
         if (size(missing) /= size(codes)) error stop "new_factor: missing mask has the wrong size"
      end if
      if (any(codes < 1 .or. codes > size(levels))) then
         if (.not. present(missing)) error stop "new_factor: code is outside the level range"
         if (any((codes < 1 .or. codes > size(levels)) .and. .not. missing)) then
            error stop "new_factor: nonmissing code is outside the level range"
         end if
      end if
      if (has_duplicates(levels)) error stop "new_factor: levels must be unique"

      factor%codes = codes
      factor%levels = levels
      allocate (factor%missing(size(codes)), source=.false.)
      if (present(missing)) factor%missing = missing
      if (present(ordered)) factor%ordered = ordered
   end function new_factor

   pure function fct_character(x, levels, missing, ordered) result(factor)
      !! Creates a factor, preserving first appearance when levels are omitted.
      character(len=*), intent(in) :: x(:)        !! Categorical observations.
      character(len=*), intent(in), optional :: levels(:) !! Explicit level order.
      logical, intent(in), optional :: missing(:) !! Missing-value mask.
      logical, intent(in), optional :: ordered    !! Whether levels are ordinal.
      type(factor_type) :: factor
      character(len=len(x)), allocatable :: inferred(:)
      logical, allocatable :: mask(:)
      integer, allocatable :: codes(:)
      integer :: i, j, number_of_levels

      allocate (mask(size(x)), source=.false.)
      if (present(missing)) then
         if (size(missing) /= size(x)) error stop "fct: missing mask has the wrong size"
         mask = missing
      end if

      if (present(levels)) then
         factor%levels = levels
      else
         allocate (inferred(size(x)))
         number_of_levels = 0
         do i = 1, size(x)
            if (mask(i)) cycle
            if (find_level(inferred(:number_of_levels), x(i)) == 0) then
               number_of_levels = number_of_levels + 1
               inferred(number_of_levels) = x(i)
            end if
         end do
         factor%levels = inferred(:number_of_levels)
      end if

      if (has_duplicates(factor%levels)) error stop "fct: levels must be unique"
      allocate (codes(size(x)), source=0)
      do i = 1, size(x)
         if (mask(i)) cycle
         j = find_level(factor%levels, x(i))
         if (j == 0) then
            mask(i) = .true.
         else
            codes(i) = j
         end if
      end do
      factor%codes = codes
      factor%missing = mask
      if (present(ordered)) factor%ordered = ordered
   end function fct_character

   pure function as_factor_character(x, missing, ordered) result(factor)
      !! Converts character values using deterministic first-appearance levels.
      character(len=*), intent(in) :: x(:)        !! Values to encode.
      logical, intent(in), optional :: missing(:) !! Missing-value mask.
      logical, intent(in), optional :: ordered    !! Whether levels are ordinal.
      type(factor_type) :: factor

      factor = fct_character(x, missing=missing, ordered=ordered)
   end function as_factor_character

   pure function as_factor_integer(x, missing, ordered) result(factor)
      !! Converts integers using ascending numeric levels.
      integer, intent(in) :: x(:)                 !! Values to encode.
      logical, intent(in), optional :: missing(:) !! Missing-value mask.
      logical, intent(in), optional :: ordered    !! Whether levels are ordinal.
      type(factor_type) :: factor
      character(len=32), allocatable :: labels(:)
      integer, allocatable :: values(:)
      logical, allocatable :: mask(:)
      integer :: i, n

      call make_mask(size(x), missing, mask)
      allocate (values(size(x)))
      n = 0
      do i = 1, size(x)
         if (mask(i) .or. any(values(:n) == x(i))) cycle
         n = n + 1
         values(n) = x(i)
      end do
      call sort_integer(values(:n))
      allocate (labels(n))
      do i = 1, n
         write (labels(i), '(i0)') values(i)
      end do
      factor = fct_character(integer_labels(x), labels, mask, ordered)
   end function as_factor_integer

   pure function as_factor_real(x, missing, ordered) result(factor)
      !! Converts real values using ascending numeric levels.
      real(dp), intent(in) :: x(:)                !! Values to encode.
      logical, intent(in), optional :: missing(:) !! Missing-value mask.
      logical, intent(in), optional :: ordered    !! Whether levels are ordinal.
      type(factor_type) :: factor
      character(len=48), allocatable :: labels(:), all_labels(:)
      real(dp), allocatable :: values(:)
      logical, allocatable :: mask(:)
      integer :: i, n

      call make_mask(size(x), missing, mask)
      mask = mask .or. ieee_is_nan(x)
      allocate (values(size(x)))
      n = 0
      do i = 1, size(x)
         if (mask(i) .or. any(values(:n) <= x(i) .and. values(:n) >= x(i))) cycle
         n = n + 1
         values(n) = x(i)
      end do
      call sort_real(values(:n))
      allocate (labels(n))
      do i = 1, n
         write (labels(i), '(g0)') values(i)
      end do
      allocate (all_labels(size(x)))
      do i = 1, size(x)
         write (all_labels(i), '(g0)') x(i)
      end do
      factor = fct_character(all_labels, labels, mask, ordered)
   end function as_factor_real

   pure function as_factor_logical(x, missing, ordered) result(factor)
      !! Converts logical values with FALSE before TRUE.
      logical, intent(in) :: x(:)                 !! Values to encode.
      logical, intent(in), optional :: missing(:) !! Missing-value mask.
      logical, intent(in), optional :: ordered    !! Whether levels are ordinal.
      type(factor_type) :: factor
      character(len=5), allocatable :: labels(:)
      integer :: i

      allocate (labels(size(x)))
      do i = 1, size(x)
         labels(i) = merge("TRUE ", "FALSE", x(i))
      end do
      factor = fct_character(labels, [character(len=5) :: "FALSE", "TRUE"], missing, ordered)
   end function as_factor_logical

   pure elemental integer function factor_size(self) result(n)
      !! Returns the number of observations.
      class(factor_type), intent(in) :: self !! Factor to inspect.

      n = 0
      if (allocated(self%codes)) n = size(self%codes)
   end function factor_size

   pure elemental integer function factor_nlevels(self) result(n)
      !! Returns the number of possible levels.
      class(factor_type), intent(in) :: self !! Factor to inspect.

      n = 0
      if (allocated(self%levels)) n = size(self%levels)
   end function factor_nlevels

   pure function factor_values(self, missing_value) result(values)
      !! Decodes factor observations to character values.
      class(factor_type), intent(in) :: self       !! Factor to decode.
      character(len=*), intent(in), optional :: missing_value !! Missing placeholder.
      character(len=:), allocatable :: values(:)
      character(len=:), allocatable :: replacement
      integer :: i, width

      replacement = "NA"
      if (present(missing_value)) replacement = missing_value
      width = len(replacement)
      if (allocated(self%levels)) width = max(width, len(self%levels))
      allocate (character(len=width) :: values(self%size()))
      do i = 1, self%size()
         if (self%missing(i)) then
            values(i) = replacement
         else
            values(i) = self%levels(self%codes(i))
         end if
      end do
   end function factor_values

   pure logical function factor_valid(self) result(valid)
      !! Reports whether all components satisfy the factor representation invariants.
      class(factor_type), intent(in) :: self !! Factor to validate.

      valid = allocated(self%codes) .and. allocated(self%levels) .and. allocated(self%missing)
      if (.not. valid) return
      valid = size(self%codes) == size(self%missing) .and. .not. has_duplicates(self%levels)
      if (.not. valid) return
      valid = all(self%missing .or. (self%codes >= 1 .and. self%codes <= size(self%levels)))
   end function factor_valid

   pure integer function find_level(levels, value) result(index)
      !! Finds an exact level label, returning zero when absent.
      character(len=*), intent(in) :: levels(:) !! Labels to search.
      character(len=*), intent(in) :: value     !! Label to find.
      integer :: i

      index = 0
      do i = 1, size(levels)
         if (levels(i) == value) then
            index = i
            return
         end if
      end do
   end function find_level

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

   pure subroutine make_mask(n, missing, mask)
      !! Creates and validates an optional missing-value mask.
      integer, intent(in) :: n                    !! Required mask size.
      logical, intent(in), optional :: missing(:) !! Optional source mask.
      logical, allocatable, intent(out) :: mask(:) !! Resulting mask.

      allocate (mask(n), source=.false.)
      if (present(missing)) then
         if (size(missing) /= n) error stop "as_factor: missing mask has the wrong size"
         mask = missing
      end if
   end subroutine make_mask

   pure function integer_labels(x) result(labels)
      !! Formats integer values for factor matching.
      integer, intent(in) :: x(:) !! Integers to format.
      character(len=32), allocatable :: labels(:)
      integer :: i

      allocate (labels(size(x)))
      do i = 1, size(x)
         write (labels(i), '(i0)') x(i)
      end do
   end function integer_labels

   pure subroutine sort_integer(x)
      !! Sorts a short integer array in ascending order.
      integer, intent(inout) :: x(:) !! Values to sort.
      integer :: i, j, value

      do i = 2, size(x)
         value = x(i)
         j = i - 1
         do while (j >= 1)
            if (x(j) <= value) exit
            x(j + 1) = x(j)
            j = j - 1
         end do
         x(j + 1) = value
      end do
   end subroutine sort_integer

   pure subroutine sort_real(x)
      !! Sorts a short real array in ascending order.
      real(dp), intent(inout) :: x(:) !! Values to sort.
      real(dp) :: value
      integer :: i, j

      do i = 2, size(x)
         value = x(i)
         j = i - 1
         do while (j >= 1)
            if (x(j) <= value) exit
            x(j + 1) = x(j)
            j = j - 1
         end do
         x(j + 1) = value
      end do
   end subroutine sort_real

end module forcats_types
