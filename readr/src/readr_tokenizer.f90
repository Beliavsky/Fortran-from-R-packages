! SPDX-License-Identifier: MIT
! SPDX-FileComment: Quote-aware tokenization translated from the R readr package behavior.
module readr_tokenizer
   !! Tokenizes delimited text, including doubled quotes and quoted newlines.
   use readr_types, only: token_table_type
   implicit none
   private

   public :: count_fields, read_file, read_lines, tokenize, write_file, write_lines

contains

   function read_file(file) result(text)
      !! Reads a complete file as an uninterpreted character sequence.
      character(len=*), intent(in) :: file !! Input path.
      character(len=:), allocatable :: text
      integer :: file_size, status, unit

      inquire (file=file, size=file_size, iostat=status)
      if (status /= 0) error stop "read_file: cannot inspect input file"
      allocate (character(len=file_size) :: text)
      open (newunit=unit, file=file, access="stream", form="unformatted", status="old", action="read", iostat=status)
      if (status /= 0) error stop "read_file: cannot open input file"
      if (file_size > 0) read (unit, iostat=status) text
      close (unit)
      if (status /= 0 .and. file_size > 0) error stop "read_file: failed while reading input file"
   end function read_file

   function read_lines(file, skip_empty_rows) result(lines)
      !! Reads logical text lines, optionally omitting empty physical lines.
      character(len=*), intent(in) :: file     !! Input path.
      logical, intent(in), optional :: skip_empty_rows !! Omit empty lines; defaults true.
      character(len=:), allocatable :: lines(:)
      character(len=:), allocatable :: text
      logical :: skip_empty
      integer :: i, line_length, maximum_length, number_of_lines, start

      text = read_file(file)
      skip_empty = .true.
      if (present(skip_empty_rows)) skip_empty = skip_empty_rows
      number_of_lines = 0
      maximum_length = 1
      start = 1
      do i = 1, len(text) + 1
         if (i <= len(text)) then
            if (text(i:i) /= achar(10) .and. text(i:i) /= achar(13)) cycle
         end if
         line_length = i - start
         if (.not. skip_empty .or. line_length > 0) then
            number_of_lines = number_of_lines + 1
            maximum_length = max(maximum_length, line_length)
         end if
         if (i <= len(text) .and. text(i:i) == achar(13)) then
            if (i < len(text)) then
               if (text(i + 1:i + 1) == achar(10)) cycle
            end if
         end if
         start = i + 1
      end do
      allocate (character(len=maximum_length) :: lines(number_of_lines))
      number_of_lines = 0
      start = 1
      i = 1
      do while (i <= len(text) + 1)
         if (i <= len(text)) then
            if (text(i:i) /= achar(10) .and. text(i:i) /= achar(13)) then
               i = i + 1
               cycle
            end if
         end if
         line_length = i - start
         if (.not. skip_empty .or. line_length > 0) then
            number_of_lines = number_of_lines + 1
            if (line_length > 0) lines(number_of_lines) = text(start:i - 1)
         end if
         if (i <= len(text) .and. text(i:i) == achar(13) .and. i < len(text)) then
            if (text(i + 1:i + 1) == achar(10)) i = i + 1
         end if
         start = i + 1
         i = i + 1
      end do
   end function read_lines

   subroutine write_file(text, file, append)
      !! Writes a complete character sequence to a file.
      character(len=*), intent(in) :: text !! Text to write.
      character(len=*), intent(in) :: file !! Output path.
      logical, intent(in), optional :: append !! Append instead of replace.
      integer :: status, unit

      if (present(append)) then
         if (append) then
            open (newunit=unit, file=file, access="stream", form="unformatted", status="unknown", &
               position="append", action="write", iostat=status)
         else
            open (newunit=unit, file=file, access="stream", form="unformatted", status="replace", &
               action="write", iostat=status)
         end if
      else
         open (newunit=unit, file=file, access="stream", form="unformatted", status="replace", &
            action="write", iostat=status)
      end if
      if (status /= 0) error stop "write_file: cannot open output file"
      if (len(text) > 0) write (unit, iostat=status) text
      close (unit)
      if (status /= 0) error stop "write_file: failed while writing output file"
   end subroutine write_file

   subroutine write_lines(lines, file)
      !! Writes character values as newline-terminated records.
      character(len=*), intent(in) :: lines(:) !! Lines to write.
      character(len=*), intent(in) :: file     !! Output path.
      integer :: i, status, unit

      open (newunit=unit, file=file, status="replace", action="write", iostat=status)
      if (status /= 0) error stop "write_lines: cannot open output file"
      do i = 1, size(lines)
         write (unit, '(a)', iostat=status) trim(lines(i))
         if (status /= 0) exit
      end do
      close (unit)
      if (status /= 0) error stop "write_lines: failed while writing output file"
   end subroutine write_lines

   pure function tokenize(text, delimiter, quote, comment) result(tokens)
      !! Splits delimited text with RFC-4180-style doubled quote escaping.
      character(len=*), intent(in) :: text      !! Text to tokenize.
      character(len=1), intent(in), optional :: delimiter !! Field delimiter; defaults comma.
      character(len=1), intent(in), optional :: quote     !! Quote character; defaults double quote.
      character(len=1), intent(in), optional :: comment   !! Optional full-line comment marker.
      type(token_table_type) :: tokens
      character(len=1) :: delim, quote_mark
      integer :: maximum_columns, maximum_width, number_of_rows

      delim = ","
      if (present(delimiter)) delim = delimiter
      quote_mark = '"'
      if (present(quote)) quote_mark = quote
      call token_dimensions(text, delim, quote_mark, comment, number_of_rows, maximum_columns, maximum_width, &
         tokens%valid, tokens%message)
      if (.not. tokens%valid) then
         allocate (character(len=1) :: tokens%field(0, 0))
         allocate (tokens%quoted(0, 0), tokens%fields_per_row(0))
         return
      end if
      allocate (character(len=max(1, maximum_width)) :: tokens%field(number_of_rows, maximum_columns))
      allocate (tokens%quoted(number_of_rows, maximum_columns), source=.false.)
      allocate (tokens%fields_per_row(number_of_rows), source=0)
      tokens%field(:,:) = ""
      call fill_tokens(text, delim, quote_mark, comment, tokens)
   end function tokenize

   function count_fields(file, delimiter) result(counts)
      !! Counts fields in every nonempty logical record of a delimited file.
      character(len=*), intent(in) :: file      !! Input path.
      character(len=1), intent(in), optional :: delimiter !! Field delimiter.
      integer, allocatable :: counts(:)
      type(token_table_type) :: tokens

      tokens = tokenize(read_file(file), delimiter)
      if (.not. tokens%valid) error stop "count_fields: " // tokens%message
      counts = tokens%fields_per_row
   end function count_fields

   pure subroutine token_dimensions(text, delimiter, quote, comment, rows, columns, width, valid, message)
      !! Determines token storage dimensions and validates quote balance.
      character(len=*), intent(in) :: text      !! Text to scan.
      character(len=1), intent(in) :: delimiter !! Field delimiter.
      character(len=1), intent(in) :: quote     !! Quote character.
      character(len=1), intent(in), optional :: comment !! Comment marker.
      integer, intent(out) :: rows, columns, width !! Required dimensions.
      logical, intent(out) :: valid             !! Whether quotes are balanced.
      character(len=:), allocatable, intent(out) :: message !! Validation message.
      logical :: in_quotes, row_has_data
      integer :: field_length, fields, i

      rows = 0
      columns = 0
      width = 0
      fields = 1
      field_length = 0
      in_quotes = .false.
      row_has_data = .false.
      i = 1
      do while (i <= len(text))
         if (.not. in_quotes .and. present(comment)) then
            if (.not. row_has_data .and. field_length == 0 .and. fields == 1 .and. text(i:i) == comment) then
               do while (i <= len(text))
                  if (is_newline(text(i:i))) exit
                  i = i + 1
               end do
               cycle
            end if
         end if
         if (text(i:i) == quote) then
            row_has_data = .true.
            if (in_quotes .and. i < len(text)) then
               if (text(i + 1:i + 1) == quote) then
                  field_length = field_length + 1
                  i = i + 2
                  cycle
               end if
            end if
            in_quotes = .not. in_quotes
         else if (.not. in_quotes .and. text(i:i) == delimiter) then
            width = max(width, field_length)
            field_length = 0
            fields = fields + 1
            row_has_data = .true.
         else if (.not. in_quotes .and. is_newline(text(i:i))) then
            if (row_has_data .or. field_length > 0) then
               rows = rows + 1
               columns = max(columns, fields)
               width = max(width, field_length)
            end if
            fields = 1
            field_length = 0
            row_has_data = .false.
            if (text(i:i) == achar(13) .and. i < len(text)) then
               if (text(i + 1:i + 1) == achar(10)) i = i + 1
            end if
         else
            field_length = field_length + 1
            row_has_data = .true.
         end if
         i = i + 1
      end do
      valid = .not. in_quotes
      if (.not. valid) then
         message = "unterminated quoted field"
         return
      end if
      if (row_has_data .or. field_length > 0) then
         rows = rows + 1
         columns = max(columns, fields)
         width = max(width, field_length)
      end if
      message = ""
   end subroutine token_dimensions

   pure elemental subroutine fill_tokens(text, delimiter, quote, comment, tokens)
      !! Performs tokenization into preallocated rectangular storage.
      character(len=*), intent(in) :: text      !! Text to scan.
      character(len=1), intent(in) :: delimiter !! Field delimiter.
      character(len=1), intent(in) :: quote     !! Quote character.
      character(len=1), intent(in), optional :: comment !! Comment marker.
      type(token_table_type), intent(inout) :: tokens !! Destination token table.
      logical :: in_quotes, row_has_data
      integer :: column, i, position, row

      row = 1
      column = 1
      position = 0
      in_quotes = .false.
      row_has_data = .false.
      i = 1
      do while (i <= len(text))
         if (.not. in_quotes .and. present(comment)) then
            if (.not. row_has_data .and. position == 0 .and. column == 1 .and. text(i:i) == comment) then
               do while (i <= len(text))
                  if (is_newline(text(i:i))) exit
                  i = i + 1
               end do
               cycle
            end if
         end if
         if (text(i:i) == quote) then
            row_has_data = .true.
            tokens%quoted(row, column) = .true.
            if (in_quotes .and. i < len(text)) then
               if (text(i + 1:i + 1) == quote) then
                  position = position + 1
                  tokens%field(row, column)(position:position) = quote
                  i = i + 2
                  cycle
               end if
            end if
            in_quotes = .not. in_quotes
         else if (.not. in_quotes .and. text(i:i) == delimiter) then
            column = column + 1
            position = 0
            row_has_data = .true.
         else if (.not. in_quotes .and. is_newline(text(i:i))) then
            if (row_has_data .or. position > 0) then
               tokens%fields_per_row(row) = column
               row = row + 1
            end if
            column = 1
            position = 0
            row_has_data = .false.
            if (text(i:i) == achar(13) .and. i < len(text)) then
               if (text(i + 1:i + 1) == achar(10)) i = i + 1
            end if
         else
            position = position + 1
            tokens%field(row, column)(position:position) = text(i:i)
            row_has_data = .true.
         end if
         i = i + 1
      end do
      if (row <= tokens%nrow() .and. (row_has_data .or. position > 0)) tokens%fields_per_row(row) = column
   end subroutine fill_tokens

   pure elemental logical function is_newline(character) result(is_line_end)
      !! Reports whether a character is CR or LF.
      character(len=1), intent(in) :: character !! Character to inspect.

      is_line_end = character == achar(10) .or. character == achar(13)
   end function is_newline

end module readr_tokenizer
