! SPDX-License-Identifier: MIT
! SPDX-FileComment: Rectangular readers translated from the R readr package behavior.
module readr_reader
   !! Converts tokenized delimited files into typed tibble columns.
   use readr_parsers, only: parse_vector
   use readr_tokenizer, only: read_file, tokenize
   use readr_types, only: col_character_code, col_double_code, col_factor_code, col_guess, &
      col_guess_code, col_integer_code, col_logical_code, col_skip_code, col_spec_type, &
      collector_type, locale, locale_type, parse_problem_type, parse_result_type, &
      read_result_type, token_table_type
   use tibble, only: make_column, new_tibble, tibble_character, tibble_column, tibble_integer, &
      tibble_logical, tibble_real
   implicit none
   private

   public :: problems, read_csv, read_csv2, read_delim, read_table, read_tsv
   public :: spec, spec_csv, spec_csv2, spec_delim, spec_table, spec_tsv, stop_for_problems

contains

   function read_delim(file, delimiter, col_types, col_names, na, comment, skip, settings) result(result)
      !! Reads an eager delimited text file into a typed tibble and diagnostics.
      character(len=*), intent(in) :: file      !! Input path.
      character(len=1), intent(in) :: delimiter !! Field delimiter.
      type(col_spec_type), intent(in), optional :: col_types !! Column specification.
      logical, intent(in), optional :: col_names !! Use first unskipped row as names; defaults true.
      character(len=*), intent(in), optional :: na(:) !! Missing tokens.
      character(len=1), intent(in), optional :: comment !! Full-line comment marker.
      integer, intent(in), optional :: skip      !! Number of token rows to skip.
      type(locale_type), intent(in), optional :: settings !! Numeric locale.
      type(read_result_type) :: result
      type(token_table_type) :: tokens

      tokens = tokenize(read_file(file), delimiter, comment=comment)
      result = tokens_to_result(tokens, col_types, col_names, na, skip, settings)
   end function read_delim

   function read_csv(file, col_types, col_names, na, comment, skip, settings) result(result)
      !! Reads a comma-delimited file.
      character(len=*), intent(in) :: file      !! Input path.
      type(col_spec_type), intent(in), optional :: col_types !! Column specification.
      logical, intent(in), optional :: col_names !! Whether the first row is a header.
      character(len=*), intent(in), optional :: na(:) !! Missing tokens.
      character(len=1), intent(in), optional :: comment !! Comment marker.
      integer, intent(in), optional :: skip      !! Token rows to skip.
      type(locale_type), intent(in), optional :: settings !! Numeric locale.
      type(read_result_type) :: result

      result = read_delim(file, ",", col_types, col_names, na, comment, skip, settings)
   end function read_csv

   function read_csv2(file, col_types, col_names, na, comment, skip, settings) result(result)
      !! Reads semicolon-delimited data using comma as the default decimal mark.
      character(len=*), intent(in) :: file      !! Input path.
      type(col_spec_type), intent(in), optional :: col_types !! Column specification.
      logical, intent(in), optional :: col_names !! Whether the first row is a header.
      character(len=*), intent(in), optional :: na(:) !! Missing tokens.
      character(len=1), intent(in), optional :: comment !! Comment marker.
      integer, intent(in), optional :: skip      !! Token rows to skip.
      type(locale_type), intent(in), optional :: settings !! Numeric locale override.
      type(read_result_type) :: result
      type(locale_type) :: active

      active = locale(decimal_mark=",", grouping_mark=".")
      if (present(settings)) active = settings
      result = read_delim(file, ";", col_types, col_names, na, comment, skip, active)
   end function read_csv2

   function read_tsv(file, col_types, col_names, na, comment, skip, settings) result(result)
      !! Reads a tab-delimited file.
      character(len=*), intent(in) :: file      !! Input path.
      type(col_spec_type), intent(in), optional :: col_types !! Column specification.
      logical, intent(in), optional :: col_names !! Whether the first row is a header.
      character(len=*), intent(in), optional :: na(:) !! Missing tokens.
      character(len=1), intent(in), optional :: comment !! Comment marker.
      integer, intent(in), optional :: skip      !! Token rows to skip.
      type(locale_type), intent(in), optional :: settings !! Numeric locale.
      type(read_result_type) :: result

      result = read_delim(file, achar(9), col_types, col_names, na, comment, skip, settings)
   end function read_tsv

   function read_table(file, col_types, col_names, na, comment, skip, settings) result(result)
      !! Reads fields separated by runs of spaces or tabs outside quoted fields.
      character(len=*), intent(in) :: file      !! Input path.
      type(col_spec_type), intent(in), optional :: col_types !! Column specification.
      logical, intent(in), optional :: col_names !! Whether the first row is a header.
      character(len=*), intent(in), optional :: na(:) !! Missing tokens.
      character(len=1), intent(in), optional :: comment !! Comment marker.
      integer, intent(in), optional :: skip      !! Token rows to skip.
      type(locale_type), intent(in), optional :: settings !! Numeric locale.
      type(read_result_type) :: result
      type(token_table_type) :: tokens

      tokens = tokenize(whitespace_delimiters(read_file(file)), achar(9), comment=comment)
      result = tokens_to_result(tokens, col_types, col_names, na, skip, settings)
   end function read_table

   function tokens_to_result(tokens, col_types, col_names, na, skip, settings) result(result)
      !! Converts rectangular tokens to typed columns and aggregates diagnostics.
      type(token_table_type), intent(in) :: tokens !! Tokenized text.
      type(col_spec_type), intent(in), optional :: col_types !! Column specification.
      logical, intent(in), optional :: col_names !! Whether first row is a header.
      character(len=*), intent(in), optional :: na(:) !! Missing tokens.
      integer, intent(in), optional :: skip      !! Token rows to skip.
      type(locale_type), intent(in), optional :: settings !! Numeric locale.
      type(read_result_type) :: result
      character(len=:), allocatable :: names(:)
      character(len=:), allocatable :: column_tokens(:)
      type(collector_type), allocatable :: collectors(:)
      type(tibble_column), allocatable :: columns(:)
      type(parse_result_type) :: parsed
      logical :: header
      integer :: data_start, i, j, number_output, skipped

      if (.not. tokens%valid) then
         call add_problem(result%problem, 0, 0, "balanced quotes", "end of input", tokens%message)
         allocate (columns(0))
         result%data = new_tibble(columns, nrow=0)
         return
      end if
      skipped = 0
      if (present(skip)) skipped = max(0, skip)
      header = .true.
      if (present(col_names)) header = col_names
      if (skipped >= tokens%nrow()) then
         allocate (columns(0))
         result%data = new_tibble(columns, nrow=0)
         return
      end if

      allocate (character(len=max(1, len(tokens%field))) :: names(tokens%ncol()))
      if (header) then
         names = tokens%field(skipped + 1, :)
         data_start = skipped + 2
      else
         do j = 1, tokens%ncol()
            names(j) = "X" // integer_text(j)
         end do
         data_start = skipped + 1
      end if
      collectors = effective_collectors(names, col_types)
      result%specification%names = names
      result%specification%collectors = collectors
      result%specification%default = col_guess()
      if (present(col_types)) result%specification%default = col_types%default

      do i = data_start, tokens%nrow()
         if (tokens%fields_per_row(i) /= tokens%ncol()) then
            call add_problem(result%problem, i, tokens%fields_per_row(i), integer_text(tokens%ncol()) // " columns", &
               integer_text(tokens%fields_per_row(i)) // " columns", "ragged row")
         end if
      end do

      number_output = count([(collectors(j)%type_code /= col_skip_code, j=1, size(collectors))])
      allocate (columns(number_output))
      number_output = 0
      do j = 1, tokens%ncol()
         if (collectors(j)%type_code == col_skip_code) cycle
         column_tokens = tokens%field(data_start:, j)
         parsed = parse_vector(column_tokens, collectors(j), na, settings)
         parsed%values%name = trim(names(j))
         number_output = number_output + 1
         columns(number_output) = parsed%values
         call append_problems(result%problem, parsed%problem, data_start - 1, j)
         if (collectors(j)%type_code == col_guess_code) collectors(j)%type_code = collector_code(parsed%values%type_code)
      end do
      result%specification%collectors = collectors
      result%data = new_tibble(columns, nrow=max(0, tokens%nrow() - data_start + 1), name_repair="unique")
   end function tokens_to_result

   pure function effective_collectors(names, specification) result(collectors)
      !! Resolves named collectors against input column names.
      character(len=*), intent(in) :: names(:) !! Input column names.
      type(col_spec_type), intent(in), optional :: specification !! Requested specification.
      type(collector_type), allocatable :: collectors(:)
      integer :: i, j

      allocate (collectors(size(names)))
      collectors = col_guess()
      if (.not. present(specification)) return
      collectors = specification%default
      if (.not. allocated(specification%names)) return
      do i = 1, size(specification%names)
         do j = 1, size(names)
            if (trim(specification%names(i)) == trim(names(j))) collectors(j) = specification%collectors(i)
         end do
      end do
   end function effective_collectors

   pure elemental integer function collector_code(tibble_code) result(code)
      !! Converts a vctrs/tibble type code to a readr collector code.
      integer, intent(in) :: tibble_code !! Source type code.

      select case (tibble_code)
      case (tibble_integer)
         code = col_integer_code
      case (tibble_real)
         code = col_double_code
      case (tibble_logical)
         code = col_logical_code
      case (tibble_character)
         code = col_character_code
      case default
         code = col_character_code
      end select
   end function collector_code

   pure function problems(result) result(problem)
      !! Returns all structured diagnostics associated with a read result.
      type(read_result_type), intent(in) :: result !! Result to inspect.
      type(parse_problem_type), allocatable :: problem(:)

      if (allocated(result%problem)) then
         problem = result%problem
      else
         allocate (problem(0))
      end if
   end function problems

   pure elemental subroutine stop_for_problems(result)
      !! Stops execution if a read result contains any parsing diagnostic.
      type(read_result_type), intent(in) :: result !! Result to validate.

      if (.not. result%ok()) error stop "stop_for_problems: input has parsing problems"
   end subroutine stop_for_problems

   pure elemental function spec(result) result(specification)
      !! Returns the effective column specification used for a read.
      type(read_result_type), intent(in) :: result !! Result to inspect.
      type(col_spec_type) :: specification

      specification = result%specification
   end function spec

   function spec_delim(file, delimiter, col_names, na, comment, skip, settings) result(specification)
      !! Guesses a column specification for a delimited file.
      character(len=*), intent(in) :: file      !! Input path.
      character(len=1), intent(in) :: delimiter !! Field delimiter.
      logical, intent(in), optional :: col_names !! Whether first row is a header.
      character(len=*), intent(in), optional :: na(:) !! Missing tokens.
      character(len=1), intent(in), optional :: comment !! Comment marker.
      integer, intent(in), optional :: skip      !! Token rows to skip.
      type(locale_type), intent(in), optional :: settings !! Numeric locale.
      type(col_spec_type) :: specification
      type(read_result_type) :: result

      result = read_delim(file, delimiter, col_names=col_names, na=na, comment=comment, skip=skip, settings=settings)
      specification = result%specification
   end function spec_delim

   function spec_csv(file) result(specification)
      !! Guesses a specification for a comma-delimited file.
      character(len=*), intent(in) :: file !! Input path.
      type(col_spec_type) :: specification

      specification = spec_delim(file, ",")
   end function spec_csv

   function spec_csv2(file) result(specification)
      !! Guesses a specification for a semicolon-delimited decimal-comma file.
      character(len=*), intent(in) :: file !! Input path.
      type(col_spec_type) :: specification

      specification = spec_delim(file, ";", settings=locale(decimal_mark=",", grouping_mark="."))
   end function spec_csv2

   function spec_tsv(file) result(specification)
      !! Guesses a specification for a tab-delimited file.
      character(len=*), intent(in) :: file !! Input path.
      type(col_spec_type) :: specification

      specification = spec_delim(file, achar(9))
   end function spec_tsv

   function spec_table(file) result(specification)
      !! Guesses a specification for a whitespace-delimited file.
      character(len=*), intent(in) :: file !! Input path.
      type(col_spec_type) :: specification
      type(read_result_type) :: result

      result = read_table(file)
      specification = result%specification
   end function spec_table

   pure function whitespace_delimiters(text) result(output)
      !! Replaces runs of unquoted horizontal whitespace with one tab delimiter.
      character(len=*), intent(in) :: text !! Whitespace-delimited text.
      character(len=:), allocatable :: output
      character(len=len(text)) :: workspace
      logical :: in_quotes, previous_space
      integer :: i, n

      in_quotes = .false.
      previous_space = .false.
      n = 0
      do i = 1, len(text)
         if (text(i:i) == '"') in_quotes = .not. in_quotes
         if (.not. in_quotes .and. (text(i:i) == " " .or. text(i:i) == achar(9))) then
            if (previous_space .or. n == 0) cycle
            n = n + 1
            workspace(n:n) = achar(9)
            previous_space = .true.
         else
            if (previous_space .and. is_line_end(text(i:i))) n = n - 1
            n = n + 1
            workspace(n:n) = text(i:i)
            previous_space = .false.
         end if
      end do
      output = workspace(:n)
   end function whitespace_delimiters

   pure elemental logical function is_line_end(character) result(line_end)
      !! Reports whether a character is a line terminator.
      character(len=1), intent(in) :: character !! Character to inspect.

      line_end = character == achar(10) .or. character == achar(13)
   end function is_line_end

   pure function integer_text(value) result(text)
      !! Formats a nonnegative integer without external I/O.
      integer, intent(in) :: value !! Value to format.
      character(len=:), allocatable :: text
      character(len=32) :: buffer

      write (buffer, '(i0)') value
      text = trim(buffer)
   end function integer_text

   pure subroutine append_problems(destination, source, row_offset, column)
      !! Appends parser diagnostics while converting their table coordinates.
      type(parse_problem_type), allocatable, intent(inout) :: destination(:) !! Aggregate diagnostics.
      type(parse_problem_type), allocatable, intent(in) :: source(:) !! Column diagnostics.
      integer, intent(in) :: row_offset, column !! Coordinate adjustment.
      integer :: i

      if (.not. allocated(source)) return
      do i = 1, size(source)
         call add_problem(destination, source(i)%row + row_offset, column, source(i)%expected, &
            source(i)%actual, source(i)%message)
      end do
   end subroutine append_problems

   pure subroutine add_problem(problem, row, column, expected, actual, message)
      !! Appends one structured read diagnostic.
      type(parse_problem_type), allocatable, intent(inout) :: problem(:) !! Diagnostic list.
      integer, intent(in) :: row, column !! One-based location.
      character(len=*), intent(in) :: expected, actual, message !! Diagnostic details.
      type(parse_problem_type), allocatable :: temporary(:)
      integer :: n

      n = 0
      if (allocated(problem)) n = size(problem)
      allocate (temporary(n + 1))
      if (n > 0) temporary(:n) = problem
      temporary(n + 1)%row = row
      temporary(n + 1)%column = column
      temporary(n + 1)%expected = expected
      temporary(n + 1)%actual = actual
      temporary(n + 1)%message = message
      call move_alloc(temporary, problem)
   end subroutine add_problem

end module readr_reader
