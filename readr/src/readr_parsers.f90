! SPDX-License-Identifier: MIT
! SPDX-FileComment: Column parsers translated from the R readr package behavior.
module readr_parsers
   !! Parses character vectors into typed columns with structured diagnostics.
   use forcats, only: as_factor, factor_type, fct
   use readr_types, only: col_character_code, col_double_code, col_factor_code, col_guess_code, &
      col_integer_code, col_logical_code, col_number_code, collector_type, dp, factor_parse_result_type, &
      locale_type, parse_problem_type, parse_result_type
   use tibble, only: make_column
   implicit none
   private

   public :: guess_parser, parse_character, parse_double, parse_factor
   public :: parse_guess, parse_integer, parse_logical, parse_number, parse_vector

contains

   pure function parse_character(x, na) result(result)
      !! Parses character fields, marking configured tokens as missing.
      character(len=*), intent(in) :: x(:)   !! Input fields.
      character(len=*), intent(in), optional :: na(:) !! Missing tokens.
      type(parse_result_type) :: result
      character(len=len(x)), allocatable :: values(:)
      logical, allocatable :: missing(:)
      integer :: i

      allocate (values(size(x)))
      allocate (missing(size(x)), source=.false.)
      do i = 1, size(x)
         values(i) = trim(adjustl(x(i)))
         missing(i) = is_missing_token(values(i), na)
      end do
      result%values = make_column("", values, missing)
   end function parse_character

   pure function parse_logical(x, na) result(result)
      !! Parses TRUE/FALSE and T/F fields without locale dependence.
      character(len=*), intent(in) :: x(:)   !! Input fields.
      character(len=*), intent(in), optional :: na(:) !! Missing tokens.
      type(parse_result_type) :: result
      logical, allocatable :: missing(:), values(:)
      character(len=:), allocatable :: token
      integer :: i

      allocate (values(size(x)), source=.false.)
      allocate (missing(size(x)), source=.false.)
      do i = 1, size(x)
         token = uppercase(trim(adjustl(x(i))))
         if (is_missing_token(trim(adjustl(x(i))), na)) then
            missing(i) = .true.
         else if (token == "TRUE" .or. token == "T") then
            values(i) = .true.
         else if (token == "FALSE" .or. token == "F") then
            values(i) = .false.
         else
            missing(i) = .true.
            call add_problem(result%problem, i, 1, "a logical", trim(adjustl(x(i))), "invalid logical")
         end if
      end do
      result%values = make_column("", values, missing)
   end function parse_logical

   pure function parse_integer(x, na) result(result)
      !! Parses base-ten integers and records conversion failures.
      character(len=*), intent(in) :: x(:)   !! Input fields.
      character(len=*), intent(in), optional :: na(:) !! Missing tokens.
      type(parse_result_type) :: result
      integer, allocatable :: values(:)
      logical, allocatable :: missing(:)
      character(len=:), allocatable :: token
      integer :: i, status

      allocate (values(size(x)), source=0)
      allocate (missing(size(x)), source=.false.)
      do i = 1, size(x)
         token = trim(adjustl(x(i)))
         if (is_missing_token(token, na)) then
            missing(i) = .true.
         else if (.not. integer_syntax(token)) then
            missing(i) = .true.
            call add_problem(result%problem, i, 1, "an integer", token, "invalid integer")
         else
            read (token, *, iostat=status) values(i)
            if (status /= 0) then
               missing(i) = .true.
               call add_problem(result%problem, i, 1, "an integer", token, "integer outside supported range")
            end if
         end if
      end do
      result%values = make_column("", values, missing)
   end function parse_integer

   pure function parse_double(x, na, settings) result(result)
      !! Parses complete floating-point fields using configured decimal marks.
      character(len=*), intent(in) :: x(:)   !! Input fields.
      character(len=*), intent(in), optional :: na(:) !! Missing tokens.
      type(locale_type), intent(in), optional :: settings !! Numeric locale.
      type(parse_result_type) :: result
      real(dp), allocatable :: values(:)
      logical, allocatable :: missing(:)
      character(len=:), allocatable :: token
      integer :: i, status

      allocate (values(size(x)), source=0.0_dp)
      allocate (missing(size(x)), source=.false.)
      do i = 1, size(x)
         token = localized_number(trim(adjustl(x(i))), settings, remove_grouping=.false.)
         if (is_missing_token(trim(adjustl(x(i))), na)) then
            missing(i) = .true.
         else if (.not. double_syntax(token)) then
            missing(i) = .true.
            call add_problem(result%problem, i, 1, "a double", trim(adjustl(x(i))), "invalid double")
         else
            read (token, *, iostat=status) values(i)
            if (status /= 0) then
               missing(i) = .true.
               call add_problem(result%problem, i, 1, "a double", trim(adjustl(x(i))), "double outside supported range")
            end if
         end if
      end do
      result%values = make_column("", values, missing)
   end function parse_double

   pure function parse_number(x, na, settings) result(result)
      !! Extracts the first flexible numeric token, ignoring grouping and adornment.
      character(len=*), intent(in) :: x(:)   !! Input fields.
      character(len=*), intent(in), optional :: na(:) !! Missing tokens.
      type(locale_type), intent(in), optional :: settings !! Numeric locale.
      type(parse_result_type) :: result
      character(len=len(x)), allocatable :: tokens(:)
      integer :: i

      allocate (tokens(size(x)))
      do i = 1, size(x)
         tokens(i) = extract_number(trim(adjustl(x(i))), settings)
      end do
      result = parse_double(tokens, na, settings)
      do i = 1, size(x)
         if (.not. is_missing_token(trim(adjustl(x(i))), na) .and. len_trim(tokens(i)) == 0) then
            if (.not. result%values%missing(i)) result%values%missing(i) = .true.
         end if
      end do
   end function parse_number

   pure function parse_factor(x, levels, ordered, na) result(result)
      !! Parses categorical fields with first-appearance or explicit levels.
      character(len=*), intent(in) :: x(:)      !! Input fields.
      character(len=*), intent(in), optional :: levels(:) !! Permitted levels.
      logical, intent(in), optional :: ordered  !! Whether the factor is ordinal.
      character(len=*), intent(in), optional :: na(:) !! Missing tokens.
      type(factor_parse_result_type) :: result
      character(len=len(x)), allocatable :: values(:)
      logical, allocatable :: missing(:)
      integer :: i

      allocate (values(size(x)), missing(size(x)))
      do i = 1, size(x)
         values(i) = trim(adjustl(x(i)))
         missing(i) = is_missing_token(values(i), na)
         if (present(levels)) then
            if (.not. missing(i) .and. find_character(levels, values(i)) == 0) then
               missing(i) = .true.
               call add_problem(result%problem, i, 1, "a value in level set", values(i), "unknown factor level")
            end if
         end if
      end do
      if (present(levels)) then
         result%values = fct(values, levels, missing, ordered)
      else
         result%values = as_factor(values, missing, ordered)
      end if
   end function parse_factor

   pure function parse_guess(x, na, settings) result(result)
      !! Guesses the narrowest common logical, integer, double, or character type.
      character(len=*), intent(in) :: x(:)   !! Input fields.
      character(len=*), intent(in), optional :: na(:) !! Missing tokens.
      type(locale_type), intent(in), optional :: settings !! Numeric locale.
      type(parse_result_type) :: result
      character(len=:), allocatable :: parser

      parser = guess_parser(x, na, settings)
      select case (parser)
      case ("logical")
         result = parse_logical(x, na)
      case ("integer")
         result = parse_integer(x, na)
      case ("double")
         result = parse_double(x, na, settings)
      case default
         result = parse_character(x, na)
      end select
   end function parse_guess

   pure function parse_vector(x, collector, na, settings) result(result)
      !! Parses a character vector according to one collector.
      character(len=*), intent(in) :: x(:)       !! Input fields.
      type(collector_type), intent(in) :: collector !! Requested type.
      character(len=*), intent(in), optional :: na(:) !! Missing tokens.
      type(locale_type), intent(in), optional :: settings !! Numeric locale.
      type(parse_result_type) :: result
      type(factor_parse_result_type) :: factor_result

      select case (collector%type_code)
      case (col_guess_code)
         result = parse_guess(x, na, settings)
      case (col_logical_code)
         result = parse_logical(x, na)
      case (col_integer_code)
         result = parse_integer(x, na)
      case (col_double_code)
         result = parse_double(x, na, settings)
      case (col_number_code)
         result = parse_number(x, na, settings)
      case (col_character_code)
         result = parse_character(x, na)
      case (col_factor_code)
         if (allocated(collector%levels)) then
            factor_result = parse_factor(x, collector%levels, collector%ordered, na)
         else
            factor_result = parse_factor(x, ordered=collector%ordered, na=na)
         end if
         result = parse_character(factor_result%values%values(), na)
         result%problem = factor_result%problem
      case default
         error stop "parse_vector: unsupported collector"
      end select
   end function parse_vector

   pure function guess_parser(x, na, settings) result(parser)
      !! Returns the narrowest parser accepting every nonmissing field.
      character(len=*), intent(in) :: x(:)   !! Input fields.
      character(len=*), intent(in), optional :: na(:) !! Missing tokens.
      type(locale_type), intent(in), optional :: settings !! Numeric locale.
      character(len=:), allocatable :: parser
      type(parse_result_type) :: attempt

      attempt = parse_logical(x, na)
      if (attempt%ok()) then
         parser = "logical"
         return
      end if
      attempt = parse_integer(x, na)
      if (attempt%ok()) then
         parser = "integer"
         return
      end if
      attempt = parse_double(x, na, settings)
      if (attempt%ok()) then
         parser = "double"
         return
      end if
      parser = "character"
   end function guess_parser

   pure logical function is_missing_token(token, na) result(missing)
      !! Reports whether a token is empty, NA, or in the caller's missing set.
      character(len=*), intent(in) :: token !! Token to inspect.
      character(len=*), intent(in), optional :: na(:) !! Additional missing tokens.

      missing = len_trim(token) == 0 .or. trim(token) == "NA"
      if (present(na)) missing = missing .or. any(trim(token) == na)
   end function is_missing_token

   pure elemental logical function integer_syntax(token) result(valid)
      !! Validates a complete signed decimal integer token.
      character(len=*), intent(in) :: token !! Token to validate.
      integer :: i, start

      valid = len(token) > 0
      if (.not. valid) return
      start = 1
      if (token(1:1) == "+" .or. token(1:1) == "-") start = 2
      if (start > len(token)) then
         valid = .false.
         return
      end if
      do i = start, len(token)
         if (token(i:i) < "0" .or. token(i:i) > "9") then
            valid = .false.
            return
         end if
      end do
   end function integer_syntax

   pure elemental logical function double_syntax(token) result(valid)
      !! Validates a complete decimal token with optional exponent.
      character(len=*), intent(in) :: token !! Token to validate.
      logical :: exponent_seen, digit_seen, decimal_seen
      integer :: i

      valid = len(token) > 0
      exponent_seen = .false.
      decimal_seen = .false.
      digit_seen = .false.
      do i = 1, len(token)
         select case (token(i:i))
         case ("0":"9")
            digit_seen = .true.
         case ("+", "-")
            if (i /= 1) then
               if (token(i - 1:i - 1) /= "e" .and. token(i - 1:i - 1) /= "E" .and. &
                   token(i - 1:i - 1) /= "d" .and. token(i - 1:i - 1) /= "D") valid = .false.
            end if
         case (".")
            if (decimal_seen .or. exponent_seen) valid = .false.
            decimal_seen = .true.
         case ("e", "E", "d", "D")
            if (exponent_seen .or. .not. digit_seen .or. i == len(token)) valid = .false.
            exponent_seen = .true.
            digit_seen = .false.
         case default
            valid = .false.
         end select
         if (.not. valid) return
      end do
      valid = valid .and. digit_seen
   end function double_syntax

   pure function localized_number(token, settings, remove_grouping) result(output)
      !! Normalizes configured decimal and optional grouping marks.
      character(len=*), intent(in) :: token !! Numeric token.
      type(locale_type), intent(in), optional :: settings !! Numeric locale.
      logical, intent(in) :: remove_grouping !! Remove grouping marks.
      character(len=:), allocatable :: output
      type(locale_type) :: active
      character(len=len(token)) :: workspace
      integer :: i, n

      active = locale_type()
      if (present(settings)) active = settings
      n = 0
      do i = 1, len(token)
         if (remove_grouping .and. token(i:i) == active%grouping_mark) cycle
         n = n + 1
         if (token(i:i) == active%decimal_mark) then
            workspace(n:n) = "."
         else
            workspace(n:n) = token(i:i)
         end if
      end do
      output = workspace(:n)
   end function localized_number

   pure function extract_number(token, settings) result(number)
      !! Extracts a permissive numeric substring and normalizes its locale.
      character(len=*), intent(in) :: token !! Decorated input token.
      type(locale_type), intent(in), optional :: settings !! Numeric locale.
      character(len=:), allocatable :: number
      type(locale_type) :: active
      character(len=len(token)) :: workspace
      logical :: started
      integer :: i, n

      active = locale_type()
      if (present(settings)) active = settings
      started = .false.
      n = 0
      do i = 1, len(token)
         if (.not. started) then
            if ((token(i:i) >= "0" .and. token(i:i) <= "9") .or. token(i:i) == "+" .or. &
                token(i:i) == "-" .or. token(i:i) == active%decimal_mark) started = .true.
         end if
         if (.not. started) cycle
         if ((token(i:i) >= "0" .and. token(i:i) <= "9") .or. token(i:i) == "+" .or. &
             token(i:i) == "-" .or. token(i:i) == active%decimal_mark .or. &
             token(i:i) == active%grouping_mark .or. index("eEdD", token(i:i)) > 0) then
            if (token(i:i) /= active%grouping_mark) then
               n = n + 1
               if (token(i:i) == active%decimal_mark) then
                  workspace(n:n) = "."
               else
                  workspace(n:n) = token(i:i)
               end if
            end if
         else
            exit
         end if
      end do
      number = workspace(:n)
   end function extract_number

   pure elemental function uppercase(text) result(output)
      !! Converts ASCII letters to uppercase.
      character(len=*), intent(in) :: text !! Text to convert.
      character(len=len(text)) :: output
      integer :: code, i

      output = text
      do i = 1, len(text)
         code = iachar(output(i:i))
         if (code >= iachar("a") .and. code <= iachar("z")) output(i:i) = achar(code - 32)
      end do
   end function uppercase

   pure integer function find_character(values, target) result(index)
      !! Finds a character value, returning zero when absent.
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

   pure subroutine add_problem(problem, row, column, expected, actual, message)
      !! Appends one structured parse diagnostic.
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

end module readr_parsers
