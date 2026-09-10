! SPDX-License-Identifier: MIT
! SPDX-FileComment: Rectangular writers translated from the R readr package behavior.
module readr_writer
   !! Formats and writes typed tibbles as delimited text.
   use readr_types, only: formatted_lines_type
   use tibble, only: tibble_character, tibble_integer, tibble_logical, tibble_real, tibble_type
   implicit none
   private

   public :: format_csv, format_csv2, format_delim, format_tsv
   public :: write_csv, write_csv2, write_delim, write_tsv

contains

   pure function format_delim(table, delimiter, na, col_names) result(output)
      !! Formats a tibble into delimited records with necessary quoting.
      type(tibble_type), intent(in) :: table     !! Table to format.
      character(len=1), intent(in) :: delimiter !! Field delimiter.
      character(len=*), intent(in), optional :: na !! Missing token; defaults NA.
      logical, intent(in), optional :: col_names !! Include header; defaults true.
      type(formatted_lines_type) :: output
      character(len=:), allocatable :: line, missing_token, token
      character(len=:), allocatable :: workspace(:)
      logical :: header
      integer :: i, j, line_number, number_of_lines, width

      header = .true.
      if (present(col_names)) header = col_names
      missing_token = "NA"
      if (present(na)) missing_token = na
      number_of_lines = table%nrow() + merge(1, 0, header)
      width = max(1, table%ncol() * (maximum_field_width(table, missing_token) + 3))
      allocate (character(len=width) :: workspace(number_of_lines))
      workspace(:) = ""
      line_number = 0
      if (header) then
         line = ""
         do j = 1, table%ncol()
            token = quote_field(table%columns(j)%name, delimiter, missing_token, character_field=.true.)
            call append_field(line, token, delimiter, j)
         end do
         line_number = 1
         workspace(line_number) = line
      end if
      do i = 1, table%nrow()
         line = ""
         do j = 1, table%ncol()
            token = encoded_value(table, i, j, delimiter, missing_token)
            call append_field(line, token, delimiter, j)
         end do
         line_number = line_number + 1
         workspace(line_number) = line
      end do
      allocate (character(len=width) :: output%line(number_of_lines))
      output%line = workspace
   end function format_delim

   pure elemental function format_csv(table, na, col_names) result(output)
      !! Formats a tibble as comma-delimited records.
      type(tibble_type), intent(in) :: table !! Table to format.
      character(len=*), intent(in), optional :: na !! Missing token.
      logical, intent(in), optional :: col_names !! Include header.
      type(formatted_lines_type) :: output

      output = format_delim(table, ",", na, col_names)
   end function format_csv

   pure elemental function format_csv2(table, na, col_names) result(output)
      !! Formats a tibble as semicolon-delimited records.
      type(tibble_type), intent(in) :: table !! Table to format.
      character(len=*), intent(in), optional :: na !! Missing token.
      logical, intent(in), optional :: col_names !! Include header.
      type(formatted_lines_type) :: output

      output = format_delim(table, ";", na, col_names)
   end function format_csv2

   pure elemental function format_tsv(table, na, col_names) result(output)
      !! Formats a tibble as tab-delimited records.
      type(tibble_type), intent(in) :: table !! Table to format.
      character(len=*), intent(in), optional :: na !! Missing token.
      logical, intent(in), optional :: col_names !! Include header.
      type(formatted_lines_type) :: output

      output = format_delim(table, achar(9), na, col_names)
   end function format_tsv

   subroutine write_delim(table, file, delimiter, na, col_names, append)
      !! Writes a tibble to a delimited text file.
      type(tibble_type), intent(in) :: table     !! Table to write.
      character(len=*), intent(in) :: file      !! Output path.
      character(len=1), intent(in) :: delimiter !! Field delimiter.
      character(len=*), intent(in), optional :: na !! Missing token.
      logical, intent(in), optional :: col_names !! Include header.
      logical, intent(in), optional :: append   !! Append to an existing file.
      type(formatted_lines_type) :: formatted
      logical :: append_file
      integer :: i, status, unit

      formatted = format_delim(table, delimiter, na, col_names)
      append_file = .false.
      if (present(append)) append_file = append
      if (append_file) then
         open (newunit=unit, file=file, status="unknown", position="append", action="write", iostat=status)
      else
         open (newunit=unit, file=file, status="replace", action="write", iostat=status)
      end if
      if (status /= 0) error stop "write_delim: cannot open output file"
      do i = 1, size(formatted%line)
         write (unit, '(a)', iostat=status) trim(formatted%line(i))
         if (status /= 0) exit
      end do
      close (unit)
      if (status /= 0) error stop "write_delim: failed while writing output file"
   end subroutine write_delim

   subroutine write_csv(table, file, na, col_names, append)
      !! Writes a tibble as comma-delimited text.
      type(tibble_type), intent(in) :: table !! Table to write.
      character(len=*), intent(in) :: file  !! Output path.
      character(len=*), intent(in), optional :: na !! Missing token.
      logical, intent(in), optional :: col_names !! Include header.
      logical, intent(in), optional :: append !! Append mode.

      call write_delim(table, file, ",", na, col_names, append)
   end subroutine write_csv

   subroutine write_csv2(table, file, na, col_names, append)
      !! Writes a tibble as semicolon-delimited text.
      type(tibble_type), intent(in) :: table !! Table to write.
      character(len=*), intent(in) :: file  !! Output path.
      character(len=*), intent(in), optional :: na !! Missing token.
      logical, intent(in), optional :: col_names !! Include header.
      logical, intent(in), optional :: append !! Append mode.

      call write_delim(table, file, ";", na, col_names, append)
   end subroutine write_csv2

   subroutine write_tsv(table, file, na, col_names, append)
      !! Writes a tibble as tab-delimited text.
      type(tibble_type), intent(in) :: table !! Table to write.
      character(len=*), intent(in) :: file  !! Output path.
      character(len=*), intent(in), optional :: na !! Missing token.
      logical, intent(in), optional :: col_names !! Include header.
      logical, intent(in), optional :: append !! Append mode.

      call write_delim(table, file, achar(9), na, col_names, append)
   end subroutine write_tsv

   pure function encoded_value(table, row, column, delimiter, na) result(token)
      !! Formats one typed table value for delimited output.
      type(tibble_type), intent(in) :: table !! Source table.
      integer, intent(in) :: row, column     !! One-based value coordinates.
      character(len=1), intent(in) :: delimiter !! Field delimiter.
      character(len=*), intent(in) :: na     !! Missing token.
      character(len=:), allocatable :: token
      character(len=64) :: buffer

      if (table%columns(column)%missing(row)) then
         token = na
         return
      end if
      select case (table%columns(column)%type_code)
      case (tibble_integer)
         write (buffer, '(i0)') table%columns(column)%integer_values(row)
         token = trim(buffer)
      case (tibble_real)
         write (buffer, '(g0)') table%columns(column)%real_values(row)
         token = trim(buffer)
      case (tibble_logical)
         token = merge("TRUE ", "FALSE", table%columns(column)%logical_values(row))
         token = trim(token)
      case (tibble_character)
         token = quote_field(table%columns(column)%character_values(row), delimiter, na, character_field=.true.)
      case default
         error stop "format_delim: unsupported column type"
      end select
   end function encoded_value

   pure function quote_field(value, delimiter, na, character_field) result(output)
      !! Applies necessary CSV quoting and doubles embedded quote characters.
      character(len=*), intent(in) :: value     !! Field value.
      character(len=1), intent(in) :: delimiter !! Field delimiter.
      character(len=*), intent(in) :: na        !! Missing token.
      logical, intent(in) :: character_field    !! Whether NA ambiguity requires quoting.
      character(len=:), allocatable :: output
      character(len=2 * len_trim(value) + 2) :: escaped
      logical :: needs_quotes
      integer :: i, n

      needs_quotes = index(value, delimiter) > 0 .or. index(value, '"') > 0 .or. &
         index(value, achar(10)) > 0 .or. index(value, achar(13)) > 0
      if (character_field) needs_quotes = needs_quotes .or. trim(value) == trim(na) .or. len_trim(value) == 0
      if (.not. needs_quotes) then
         output = trim(value)
         return
      end if
      n = 1
      escaped(n:n) = '"'
      do i = 1, len_trim(value)
         n = n + 1
         escaped(n:n) = value(i:i)
         if (value(i:i) == '"') then
            n = n + 1
            escaped(n:n) = '"'
         end if
      end do
      n = n + 1
      escaped(n:n) = '"'
      output = escaped(:n)
   end function quote_field

   pure subroutine append_field(line, field, delimiter, column)
      !! Appends one encoded field to a growing record.
      character(len=:), allocatable, intent(inout) :: line !! Record under construction.
      character(len=*), intent(in) :: field       !! Encoded field.
      character(len=1), intent(in) :: delimiter   !! Field separator.
      integer, intent(in) :: column               !! One-based field number.

      if (column == 1) then
         line = field
      else
         line = line // delimiter // field
      end if
   end subroutine append_field

   pure integer function maximum_field_width(table, na) result(width)
      !! Returns a conservative maximum encoded field width.
      type(tibble_type), intent(in) :: table !! Table to inspect.
      character(len=*), intent(in) :: na     !! Missing token.
      integer :: j

      width = max(64, len(na))
      do j = 1, table%ncol()
         width = max(width, len(table%columns(j)%name))
         if (table%columns(j)%type_code == tibble_character) then
            if (allocated(table%columns(j)%character_values)) then
               width = max(width, 2 * len(table%columns(j)%character_values) + 2)
            end if
         end if
      end do
   end function maximum_field_width

end module readr_writer
