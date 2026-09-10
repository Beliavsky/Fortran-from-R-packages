! SPDX-License-Identifier: MIT
! SPDX-FileComment: Vector operations for the Fortran stringr translation.
module stringr_vector
   !! Implements joining, splitting, ordering, ranking, subsets, and uniqueness.
   use stringr_scalar, only: str_detect
   use stringr_types, only: split_result_type
   implicit none
   private

   public :: str_c, str_flatten, str_flatten_comma, str_order, str_rank
   public :: str_sort, str_split, str_split_1, str_split_fixed, str_split_i
   public :: str_subset, str_unique, str_which

contains

   pure function str_c(strings, separator) result(output)
      !! Concatenates a vector using an optional separator.
      character(len=*), intent(in) :: strings(:) !! Input strings.
      character(len=*), intent(in), optional :: separator !! Separator; empty by default.
      character(len=:), allocatable :: output
      character(len=:), allocatable :: sep
      integer :: i
      sep = ''
      if (present(separator)) sep = separator
      output = ''
      do i = 1, size(strings)
         if (i > 1) output = output // sep
         output = output // trim(strings(i))
      end do
   end function str_c

   pure function str_flatten(strings, collapse, last) result(output)
      !! Collapses strings, optionally using a distinct final separator.
      character(len=*), intent(in) :: strings(:) !! Input strings.
      character(len=*), intent(in), optional :: collapse !! Separator; empty by default.
      character(len=*), intent(in), optional :: last !! Final separator.
      character(len=:), allocatable :: output
      character(len=:), allocatable :: separator
      integer :: i
      separator = ''
      if (present(collapse)) separator = collapse
      output = ''
      do i = 1, size(strings)
         if (i > 1) then
            if (present(last) .and. i == size(strings)) then
               output = output // last
            else
               output = output // separator
            end if
         end if
         output = output // trim(strings(i))
      end do
   end function str_flatten

   pure function str_flatten_comma(strings, last) result(output)
      !! Collapses strings with comma-space separators.
      character(len=*), intent(in) :: strings(:) !! Input strings.
      character(len=*), intent(in), optional :: last !! Final separator.
      character(len=:), allocatable :: output
      if (present(last)) then
         output = str_flatten(strings, ', ', last)
      else
         output = str_flatten(strings, ', ')
      end if
   end function str_flatten_comma

   pure function str_order(strings, descending) result(order)
      !! Returns a stable bytewise lexical ordering.
      character(len=*), intent(in) :: strings(:) !! Input strings.
      logical, intent(in), optional :: descending !! Reverse order when true.
      integer, allocatable :: order(:)
      logical :: reverse
      integer :: i, j, key
      reverse = .false.
      if (present(descending)) reverse = descending
      allocate (order(size(strings)))
      do i = 1, size(order)
         order(i) = i
      end do
      do i = 2, size(order)
         key = order(i)
         j = i - 1
         do while (j >= 1)
            if (reverse) then
               if (strings(key) <= strings(order(j))) exit
            else
               if (strings(key) >= strings(order(j))) exit
            end if
            order(j + 1) = order(j)
            j = j - 1
         end do
         order(j + 1) = key
      end do
   end function str_order

   pure function str_sort(strings, descending) result(sorted)
      !! Sorts strings using stable bytewise lexical order.
      character(len=*), intent(in) :: strings(:) !! Input strings.
      logical, intent(in), optional :: descending !! Reverse order when true.
      character(len=len(strings)), allocatable :: sorted(:)
      sorted = strings(str_order(strings, descending))
   end function str_sort

   pure function str_rank(strings, descending) result(rank)
      !! Returns ordinal ranks in original input order.
      character(len=*), intent(in) :: strings(:) !! Input strings.
      logical, intent(in), optional :: descending !! Reverse order when true.
      integer, allocatable :: rank(:)
      integer, allocatable :: order(:)
      integer :: i
      order = str_order(strings, descending)
      allocate (rank(size(strings)))
      do i = 1, size(order)
         rank(order(i)) = i
      end do
   end function str_rank

   pure function str_unique(strings) result(unique)
      !! Retains the first occurrence of each distinct string.
      character(len=*), intent(in) :: strings(:) !! Input strings.
      character(len=len(strings)), allocatable :: unique(:)
      logical, allocatable :: keep(:)
      integer :: i, j
      allocate (keep(size(strings)), source=.true.)
      do i = 1, size(strings)
         do j = 1, i - 1
            if (strings(i) == strings(j)) then
               keep(i) = .false.
               exit
            end if
         end do
      end do
      unique = pack(strings, keep)
   end function str_unique

   pure function str_which(strings, pattern) result(indices)
      !! Returns indices of strings containing a literal pattern.
      character(len=*), intent(in) :: strings(:) !! Input strings.
      character(len=*), intent(in) :: pattern    !! Literal pattern.
      integer, allocatable :: indices(:)
      logical, allocatable :: keep(:)
      integer :: i
      allocate (keep(size(strings)))
      do i = 1, size(strings)
         keep(i) = str_detect(strings(i), pattern)
      end do
      indices = pack([(i, i = 1, size(strings))], keep)
   end function str_which

   pure function str_subset(strings, pattern, negate) result(subset)
      !! Selects strings by literal containment.
      character(len=*), intent(in) :: strings(:) !! Input strings.
      character(len=*), intent(in) :: pattern    !! Literal pattern.
      logical, intent(in), optional :: negate    !! Select nonmatches when true.
      character(len=len(strings)), allocatable :: subset(:)
      logical, allocatable :: keep(:)
      logical :: invert
      integer :: i
      invert = .false.
      if (present(negate)) invert = negate
      allocate (keep(size(strings)))
      do i = 1, size(strings)
         keep(i) = str_detect(strings(i), pattern) .neqv. invert
      end do
      subset = pack(strings, keep)
   end function str_subset

   pure function str_split(strings, pattern, n) result(parts)
      !! Splits each string by a literal delimiter into a rectangular result.
      character(len=*), intent(in) :: strings(:) !! Input strings.
      character(len=*), intent(in) :: pattern    !! Nonempty literal delimiter.
      integer, intent(in), optional :: n         !! Maximum pieces per string.
      type(split_result_type) :: parts
      integer :: columns, i, j, limit
      if (len(pattern) == 0) error stop 'str_split: pattern cannot be empty'
      limit = huge(0)
      if (present(n)) limit = max(1, n)
      allocate (parts%counts(size(strings)))
      columns = 0
      do i = 1, size(strings)
         parts%counts(i) = min(piece_count(strings(i), pattern), limit)
         columns = max(columns, parts%counts(i))
      end do
      allocate (character(len=max(1, len(strings))) :: parts%values(size(strings), columns))
      parts%values = ''
      do i = 1, size(strings)
         do j = 1, parts%counts(i)
            parts%values(i, j) = split_piece(strings(i), pattern, j, parts%counts(i))
         end do
      end do
   end function str_split

   pure function str_split_1(string, pattern, n) result(parts)
      !! Splits one string and returns its fields as a vector.
      character(len=*), intent(in) :: string  !! Input string.
      character(len=*), intent(in) :: pattern !! Literal delimiter.
      integer, intent(in), optional :: n      !! Maximum fields.
      character(len=len(string)), allocatable :: parts(:)
      type(split_result_type) :: matrix
      matrix = str_split([string], pattern, n)
      allocate (parts(matrix%counts(1)))
      parts = matrix%values(1, :matrix%counts(1))
   end function str_split_1

   pure function str_split_fixed(strings, pattern, n) result(parts)
      !! Splits strings into exactly `n` rectangular columns, blank-filling absent fields.
      character(len=*), intent(in) :: strings(:) !! Input strings.
      character(len=*), intent(in) :: pattern    !! Literal delimiter.
      integer, intent(in) :: n                   !! Output column count.
      character(len=len(strings)), allocatable :: parts(:,:)
      type(split_result_type) :: result
      result = str_split(strings, pattern, n)
      allocate (parts(size(strings), n))
      parts = ''
      if (size(result%values, 2) > 0) then
         parts(:, :min(n, size(result%values, 2))) = result%values(:, :min(n, size(result%values, 2)))
      end if
   end function str_split_fixed

   pure function str_split_i(strings, pattern, i) result(parts)
      !! Extracts one field from each literal split, blank-filling absent fields.
      character(len=*), intent(in) :: strings(:) !! Input strings.
      character(len=*), intent(in) :: pattern    !! Literal delimiter.
      integer, intent(in) :: i                   !! One-based field index.
      character(len=len(strings)), allocatable :: parts(:)
      type(split_result_type) :: result
      integer :: row
      result = str_split(strings, pattern, i)
      allocate (parts(size(strings)))
      parts = ''
      do row = 1, size(strings)
         if (result%counts(row) >= i) parts(row) = result%values(row, i)
      end do
   end function str_split_i

   pure elemental integer function piece_count(text, delimiter) result(count)
      !! Counts literal-delimiter-separated fields.
      character(len=*), intent(in) :: text      !! Input text.
      character(len=*), intent(in) :: delimiter !! Literal delimiter.
      integer :: offset, position
      count = 1
      offset = 1
      do while (offset <= len_trim(text))
         position = index(text(offset:len_trim(text)), delimiter)
         if (position == 0) exit
         count = count + 1
         offset = offset + position - 1 + len(delimiter)
      end do
   end function piece_count

   pure subroutine split_one(text, delimiter, fields, count)
      !! Splits one string into caller-sized fields.
      character(len=*), intent(in) :: text       !! Input text.
      character(len=*), intent(in) :: delimiter  !! Literal delimiter.
      character(len=*), intent(out) :: fields(:) !! Output fields.
      integer, intent(in) :: count               !! Maximum output fields.
      integer :: field, offset, position
      fields = ''
      offset = 1
      do field = 1, count
         if (field == count) then
            fields(field) = text(offset:len_trim(text))
            exit
         end if
         position = index(text(offset:len_trim(text)), delimiter)
         if (position == 0) then
            fields(field) = text(offset:len_trim(text))
            exit
         end if
         fields(field) = text(offset:offset + position - 2)
         offset = offset + position - 1 + len(delimiter)
      end do
   end subroutine split_one

   pure elemental function split_piece(text, delimiter, wanted, limit) result(piece)
      !! Extracts one field from a literal split without an array-section temporary.
      character(len=*), intent(in) :: text      !! Input text.
      character(len=*), intent(in) :: delimiter !! Literal delimiter.
      integer, intent(in) :: wanted             !! One-based requested field.
      integer, intent(in) :: limit              !! Maximum split-field count.
      character(len=len(text)) :: piece
      integer :: field, offset, position
      piece = ''
      offset = 1
      do field = 1, wanted
         if (field == limit) then
            if (field == wanted) piece = text(offset:len_trim(text))
            return
         end if
         position = index(text(offset:len_trim(text)), delimiter)
         if (position == 0) then
            if (field == wanted) piece = text(offset:len_trim(text))
            return
         end if
         if (field == wanted) then
            piece = text(offset:offset + position - 2)
            return
         end if
         offset = offset + position - 1 + len(delimiter)
      end do
   end function split_piece

end module stringr_vector
