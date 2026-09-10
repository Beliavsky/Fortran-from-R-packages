! SPDX-License-Identifier: MIT
! SPDX-FileComment: Character-column operations for the Fortran tidyr translation.
module tidyr_strings
   !! Separates, lengthens, and unites rectangular character data.
   use tibble, only: add_column, drop_columns, make_column, new_tibble, replace_column
   use tibble, only: tibble_character, tibble_column, tibble_type
   use stringr, only: str_split_1
   use vctrs, only: vec_slice
   implicit none
   private

   public :: extract_numeric, separate, separate_longer_delim, separate_rows
   public :: separate_wider_delim, unite

contains

   function separate_wider_delim(data, column, names, delim, too_few) result(output)
      !! Splits one character column into a fixed set of character columns.
      type(tibble_type), intent(in) :: data !! Input table.
      character(len=*), intent(in) :: column !! Character column to split.
      character(len=*), intent(in) :: names(:) !! Output column names.
      character(len=*), intent(in) :: delim !! Literal nonempty delimiter.
      character(len=*), intent(in), optional :: too_few !! `error` or `align_start`.
      type(tibble_type) :: output
      character(len=256), allocatable :: values(:,:)
      logical, allocatable :: missing(:,:)
      character(len=:), allocatable :: policy
      integer :: i, index, count
      character(len=:), allocatable :: pieces(:)

      if (len(delim) == 0) error stop 'separate_wider_delim: delimiter cannot be empty'
      index = data%column_index(column)
      if (index == 0) error stop 'separate_wider_delim: column not found'
      if (data%columns(index)%type_code /= tibble_character) then
         error stop 'separate_wider_delim: column must be character'
      end if
      policy = 'error'
      if (present(too_few)) policy = trim(too_few)
      allocate (values(data%nrow(), size(names)))
      values = ''
      allocate (missing(data%nrow(), size(names)), source=.false.)
      do i = 1, data%nrow()
         if (data%columns(index)%missing(i)) then
            missing(i, :) = .true.
         else
            pieces = str_split_1(data%columns(index)%character_values(i), delim)
            count = size(pieces)
            values(i, :min(count, size(names))) = pieces(:min(count, size(names)))
            if (count < size(names)) then
               if (policy == 'error') error stop 'separate_wider_delim: too few pieces'
               missing(i, count + 1:) = .true.
            end if
            if (count > size(names)) error stop 'separate_wider_delim: too many pieces'
         end if
      end do
      output = drop_columns(data, [character(len=len(column)) :: column])
      do i = 1, size(names)
         output = add_column(output, make_column(names(i), values(:, i), missing(:, i)))
      end do
   end function separate_wider_delim

   function separate(data, column, into, sep) result(output)
      !! Legacy alias for delimiter-based separation.
      type(tibble_type), intent(in) :: data !! Input table.
      character(len=*), intent(in) :: column !! Character column to split.
      character(len=*), intent(in) :: into(:) !! Output names.
      character(len=*), intent(in) :: sep !! Literal delimiter.
      type(tibble_type) :: output
      output = separate_wider_delim(data, column, into, sep)
   end function separate

   function unite(data, new_name, columns, sep, na_rm) result(output)
      !! Joins selected columns as one character column.
      type(tibble_type), intent(in) :: data !! Input table.
      character(len=*), intent(in) :: new_name !! Output column name.
      character(len=*), intent(in) :: columns(:) !! Columns to join.
      character(len=*), intent(in), optional :: sep !! Separator; underscore by default.
      logical, intent(in), optional :: na_rm !! Omit missing components when true.
      type(tibble_type) :: output
      character(len=512), allocatable :: values(:)
      logical, allocatable :: missing(:)
      character(len=:), allocatable :: separator
      character(len=128) :: piece
      logical :: omit_missing
      integer :: i, j, index

      separator = '_'
      if (present(sep)) separator = sep
      omit_missing = .false.
      if (present(na_rm)) omit_missing = na_rm
      allocate (values(data%nrow()))
      values = ''
      allocate (missing(data%nrow()), source=.false.)
      do i = 1, data%nrow()
         do j = 1, size(columns)
            index = data%column_index(columns(j))
            if (index == 0) error stop 'unite: column not found'
            if (data%columns(index)%missing(i)) then
               if (.not. omit_missing) missing(i) = .true.
               cycle
            end if
            piece = element_text(data%columns(index), i)
            if (len_trim(values(i)) > 0) values(i) = trim(values(i)) // separator
            values(i) = trim(values(i)) // trim(piece)
         end do
      end do
      output = drop_columns(data, columns)
      output = add_column(output, make_column(new_name, values, missing))
   end function unite

   function separate_longer_delim(data, column, delim) result(output)
      !! Splits each character cell into rows while repeating other columns.
      type(tibble_type), intent(in) :: data !! Input table.
      character(len=*), intent(in) :: column !! Character column to split.
      character(len=*), intent(in) :: delim !! Literal nonempty delimiter.
      type(tibble_type) :: output
      type(tibble_column), allocatable :: columns(:)
      character(len=256), allocatable :: tokens(:)
      logical, allocatable :: token_missing(:)
      integer, allocatable :: rows(:)
      integer :: i, index, j, k, number, total

      if (len(delim) == 0) error stop 'separate_longer_delim: delimiter cannot be empty'
      index = data%column_index(column)
      if (index == 0) error stop 'separate_longer_delim: column not found'
      if (data%columns(index)%type_code /= tibble_character) then
         error stop 'separate_longer_delim: column must be character'
      end if
      total = 0
      do i = 1, data%nrow()
         if (data%columns(index)%missing(i)) then
            total = total + 1
         else
            total = total + token_count(data%columns(index)%character_values(i), delim)
         end if
      end do
      allocate (rows(total))
      allocate (tokens(total))
      tokens = ''
      allocate (token_missing(total), source=.false.)
      k = 0
      do i = 1, data%nrow()
         if (data%columns(index)%missing(i)) then
            k = k + 1
            rows(k) = i
            tokens(k) = ''
            token_missing(k) = .true.
         else
            number = token_count(data%columns(index)%character_values(i), delim)
            call split_fixed(data%columns(index)%character_values(i), delim, tokens(k + 1:k + number), j)
            rows(k + 1:k + number) = i
            k = k + number
         end if
      end do
      allocate (columns(data%ncol()))
      do j = 1, data%ncol()
         columns(j) = vec_slice(data%columns(j), rows)
      end do
      columns(index) = make_column(column, tokens, token_missing)
      output = new_tibble(columns, nrow=total)
   end function separate_longer_delim

   function separate_rows(data, column, sep) result(output)
      !! Legacy alias for `separate_longer_delim`.
      type(tibble_type), intent(in) :: data !! Input table.
      character(len=*), intent(in) :: column !! Character column to split.
      character(len=*), intent(in) :: sep !! Literal delimiter.
      type(tibble_type) :: output
      output = separate_longer_delim(data, column, sep)
   end function separate_rows

   pure function extract_numeric(text) result(value)
      !! Extracts characters that can form one signed decimal number.
      character(len=*), intent(in) :: text !! Text containing a number.
      character(len=:), allocatable :: value
      character(len=len(text)) :: buffer
      integer :: i, n
      n = 0
      buffer = ''
      do i = 1, len_trim(text)
         if (index('0123456789+-.eE', text(i:i)) > 0) then
            n = n + 1
            buffer(n:n) = text(i:i)
         end if
      end do
      value = buffer(:n)
   end function extract_numeric

   pure subroutine split_fixed(text, delim, pieces, count)
      !! Splits text into caller-sized fixed-length fields.
      character(len=*), intent(in) :: text !! Text to split.
      character(len=*), intent(in) :: delim !! Literal delimiter.
      character(len=*), intent(out) :: pieces(:) !! Extracted pieces.
      integer, intent(out) :: count !! Actual piece count, possibly exceeding capacity.
      integer :: position, relative, start
      pieces = ''
      count = 0
      start = 1
      do
         relative = index(text(start:), delim)
         count = count + 1
         if (count <= size(pieces)) then
            if (relative == 0) then
               pieces(count) = text(start:)
            else
               position = start + relative - 1
               pieces(count) = text(start:position - 1)
            end if
         end if
         if (relative == 0) exit
         start = start + relative - 1 + len(delim)
         if (start > len(text) + 1) exit
      end do
   end subroutine split_fixed

   pure integer function token_count(text, delim) result(count)
      !! Counts literal-delimiter-separated fields, including empty fields.
      character(len=*), intent(in) :: text !! Text to inspect.
      character(len=*), intent(in) :: delim !! Literal delimiter.
      integer :: relative, start
      count = 1
      start = 1
      do
         relative = index(text(start:), delim)
         if (relative == 0) exit
         count = count + 1
         start = start + relative - 1 + len(delim)
         if (start > len(text)) exit
      end do
   end function token_count

   function element_text(vector, row) result(text)
      !! Formats one nonmissing vector element compactly.
      type(tibble_column), intent(in) :: vector !! Source column.
      integer, intent(in) :: row                !! Source row.
      character(len=128) :: text
      text = ''
      select case (vector%type_code)
      case (1)
         write (text, '(i0)') vector%integer_values(row)
      case (2)
         write (text, '(g0)') vector%real_values(row)
      case (3)
         text = merge('TRUE ', 'FALSE', vector%logical_values(row))
      case (4)
         text = vector%character_values(row)
      case default
         error stop 'unite: unsupported column type'
      end select
   end function element_text

end module tidyr_strings
