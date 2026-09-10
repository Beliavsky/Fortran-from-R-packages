! SPDX-License-Identifier: MIT
module dplyr_expression
   use dplyr_vector_ops
   use tibble, only: tibble_type
   use vctrs, only: dp, new_vctr, vctr_type, vctrs_logical
   implicit none
   private

   type :: parser_state
      !! Tracks a position in one filter expression.
      character(len=:), allocatable :: source
      integer :: position = 1
   end type parser_state

   public :: evaluate_expression

contains

   function evaluate_expression(data, expression) result(value)
      !! Parses and evaluates a scalar-broadcasting vector expression against a table.
      type(tibble_type), intent(in) :: data !! Table supplying named columns.
      character(len=*), intent(in) :: expression !! Expression text.
      type(vctr_type) :: value
      type(parser_state) :: state

      state%source = trim(expression)
      if (len(state%source) == 0) error stop "filter expression is empty"
      value = parse_or(state, data)
      call skip_spaces(state)
      if (state%position <= len(state%source)) then
         error stop "unexpected text in filter expression near: " // state%source(state%position:)
      end if
   end function evaluate_expression

   recursive function parse_or(state, data) result(value)
      !! Parses left-associative logical OR expressions.
      type(parser_state), intent(inout) :: state !! Parser cursor.
      type(tibble_type), intent(in) :: data      !! Column environment.
      type(vctr_type) :: value, right

      value = parse_and(state, data)
      do while (match_or_operator(state))
         right = parse_and(state, data)
         value = value .or. right
      end do
   end function parse_or

   recursive function parse_and(state, data) result(value)
      !! Parses left-associative logical AND expressions.
      type(parser_state), intent(inout) :: state !! Parser cursor.
      type(tibble_type), intent(in) :: data      !! Column environment.
      type(vctr_type) :: value, right

      value = parse_comparison(state, data)
      do while (match_and_operator(state))
         right = parse_comparison(state, data)
         value = value .and. right
      end do
   end function parse_and

   recursive function parse_comparison(state, data) result(value)
      !! Parses at most one relational comparison.
      type(parser_state), intent(inout) :: state !! Parser cursor.
      type(tibble_type), intent(in) :: data      !! Column environment.
      type(vctr_type) :: value, right

      value = parse_additive(state, data)
      if (match_text(state, "==")) then
         right = parse_additive(state, data)
         value = value == right
      else if (match_not_equal_operator(state)) then
         right = parse_additive(state, data)
         value = value /= right
      else if (match_text(state, "<=")) then
         right = parse_additive(state, data)
         value = value <= right
      else if (match_text(state, ">=")) then
         right = parse_additive(state, data)
         value = value >= right
      else if (match_text(state, "<")) then
         right = parse_additive(state, data)
         value = value < right
      else if (match_text(state, ">")) then
         right = parse_additive(state, data)
         value = value > right
      end if
   end function parse_comparison

   recursive function parse_additive(state, data) result(value)
      !! Parses addition and subtraction.
      type(parser_state), intent(inout) :: state !! Parser cursor.
      type(tibble_type), intent(in) :: data      !! Column environment.
      type(vctr_type) :: value, right

      value = parse_multiplicative(state, data)
      do
         if (match_text(state, "+")) then
            right = parse_multiplicative(state, data)
            value = value + right
         else if (match_text(state, "-")) then
            right = parse_multiplicative(state, data)
            value = value - right
         else
            exit
         end if
      end do
   end function parse_additive

   recursive function parse_multiplicative(state, data) result(value)
      !! Parses multiplication and division.
      type(parser_state), intent(inout) :: state !! Parser cursor.
      type(tibble_type), intent(in) :: data      !! Column environment.
      type(vctr_type) :: value, right

      value = parse_unary(state, data)
      do
         if (match_text(state, "*")) then
            right = parse_unary(state, data)
            value = value * right
         else if (match_text(state, "/")) then
            right = parse_unary(state, data)
            value = value / right
         else
            exit
         end if
      end do
   end function parse_multiplicative

   recursive function parse_unary(state, data) result(value)
      !! Parses logical negation and numeric unary signs.
      type(parser_state), intent(inout) :: state !! Parser cursor.
      type(tibble_type), intent(in) :: data      !! Column environment.
      type(vctr_type) :: value

      if (match_not_operator(state)) then
         value = .not. parse_unary(state, data)
      else if (match_text(state, "-")) then
         value = -parse_unary(state, data)
      else if (match_text(state, "+")) then
         value = parse_unary(state, data)
      else
         value = parse_primary(state, data)
      end if
   end function parse_unary

   recursive function parse_primary(state, data) result(value)
      !! Parses parenthesized expressions, literals, columns, and supported functions.
      type(parser_state), intent(inout) :: state !! Parser cursor.
      type(tibble_type), intent(in) :: data      !! Column environment.
      type(vctr_type) :: value, argument
      character(len=:), allocatable :: identifier, lowered, text
      character(len=1) :: current

      call skip_spaces(state)
      if (state%position > len(state%source)) error stop "expected an expression operand"
      current = state%source(state%position:state%position)
      if (current == "(") then
         state%position = state%position + 1
         value = parse_or(state, data)
         call require_text(state, ")")
      else if (current == "'" .or. current == '"') then
         text = parse_string(state)
         value = new_vctr("", text, 1)
      else if (is_digit(current) .or. current == ".") then
         value = parse_number(state)
      else if (is_identifier_start(current)) then
         identifier = parse_identifier(state)
         lowered = lower_ascii(identifier)
         select case (lowered)
         case ("true")
            value = new_vctr("", .true., 1)
         case ("false")
            value = new_vctr("", .false., 1)
         case ("missing", "abs")
            call require_text(state, "(")
            argument = parse_or(state, data)
            call require_text(state, ")")
            if (lowered == "missing") then
               value = is_missing(argument)
            else
               value = vector_abs(argument)
            end if
         case default
            value = col(data, identifier)
         end select
      else
         error stop "invalid expression operand near: " // state%source(state%position:)
      end if
   end function parse_primary

   function parse_number(state) result(value)
      !! Parses a default integer or double-precision literal.
      type(parser_state), intent(inout) :: state !! Parser cursor.
      type(vctr_type) :: value
      character(len=:), allocatable :: token
      integer :: finish, integer_value, io_status
      real(dp) :: real_value
      logical :: is_real

      call skip_spaces(state)
      finish = state%position
      is_real = .false.
      do while (finish <= len(state%source))
         if (index("0123456789.eEdD+-", state%source(finish:finish)) == 0) exit
         if (index(".eEdD", state%source(finish:finish)) > 0) is_real = .true.
         if ((state%source(finish:finish) == "+" .or. state%source(finish:finish) == "-") .and. &
             finish > state%position) then
            if (index("eEdD", state%source(finish - 1:finish - 1)) == 0) exit
         end if
         finish = finish + 1
      end do
      token = state%source(state%position:finish - 1)
      state%position = finish
      if (is_real) then
         read (token, *, iostat=io_status) real_value
         if (io_status /= 0) error stop "invalid real literal: " // token
         value = new_vctr("", real_value, 1)
      else
         read (token, *, iostat=io_status) integer_value
         if (io_status /= 0) error stop "invalid integer literal: " // token
         value = new_vctr("", integer_value, 1)
      end if
   end function parse_number

   function parse_string(state) result(value)
      !! Parses a single- or double-quoted character literal with doubled-quote escapes.
      type(parser_state), intent(inout) :: state !! Parser cursor.
      character(len=:), allocatable :: value
      character(len=1) :: quote

      quote = state%source(state%position:state%position)
      state%position = state%position + 1
      value = ""
      do while (state%position <= len(state%source))
         if (state%source(state%position:state%position) == quote) then
            if (state%position < len(state%source)) then
               if (state%source(state%position + 1:state%position + 1) == quote) then
                  value = value // quote
                  state%position = state%position + 2
                  cycle
               end if
            end if
            state%position = state%position + 1
            return
         end if
         value = value // state%source(state%position:state%position)
         state%position = state%position + 1
      end do
      error stop "unterminated character literal"
   end function parse_string

   function parse_identifier(state) result(identifier)
      !! Parses an unquoted Fortran-like column or function name.
      type(parser_state), intent(inout) :: state !! Parser cursor.
      character(len=:), allocatable :: identifier
      integer :: finish

      call skip_spaces(state)
      finish = state%position + 1
      do while (finish <= len(state%source))
         if (.not. is_identifier_part(state%source(finish:finish))) exit
         finish = finish + 1
      end do
      identifier = state%source(state%position:finish - 1)
      state%position = finish
   end function parse_identifier

   logical function match_text(state, token) result(matched)
      !! Consumes exact punctuation when it occurs at the current nonblank position.
      type(parser_state), intent(inout) :: state !! Parser cursor.
      character(len=*), intent(in) :: token      !! Token to match.

      call skip_spaces(state)
      matched = .false.
      if (state%position + len(token) - 1 > len(state%source)) return
      if (token == "/" .and. state%position < len(state%source)) then
         if (state%source(state%position + 1:state%position + 1) == "=") return
      end if
      if (state%source(state%position:state%position + len(token) - 1) == token) then
         state%position = state%position + len(token)
         matched = .true.
      end if
   end function match_text

   logical function match_keyword(state, keyword) result(matched)
      !! Consumes a case-insensitive keyword bounded by nonidentifier characters.
      type(parser_state), intent(inout) :: state !! Parser cursor.
      character(len=*), intent(in) :: keyword    !! Lowercase keyword.
      integer :: finish

      call skip_spaces(state)
      matched = .false.
      finish = state%position + len(keyword) - 1
      if (finish > len(state%source)) return
      if (lower_ascii(state%source(state%position:finish)) /= keyword) return
      if (finish < len(state%source)) then
         if (is_identifier_part(state%source(finish + 1:finish + 1))) return
      end if
      state%position = finish + 1
      matched = .true.
   end function match_keyword

   logical function match_or_operator(state) result(matched)
      !! Consumes either the `or` keyword or its vertical-bar alias.
      type(parser_state), intent(inout) :: state !! Parser cursor.

      matched = match_keyword(state, "or")
      if (.not. matched) matched = match_text(state, "|")
   end function match_or_operator

   logical function match_and_operator(state) result(matched)
      !! Consumes either the `and` keyword or its ampersand alias.
      type(parser_state), intent(inout) :: state !! Parser cursor.

      matched = match_keyword(state, "and")
      if (.not. matched) matched = match_text(state, "&")
   end function match_and_operator

   logical function match_not_equal_operator(state) result(matched)
      !! Consumes either `/=` or the R-style `!=` alias.
      type(parser_state), intent(inout) :: state !! Parser cursor.

      matched = match_text(state, "/=")
      if (.not. matched) matched = match_text(state, "!=")
   end function match_not_equal_operator

   logical function match_not_operator(state) result(matched)
      !! Consumes either the `not` keyword or its exclamation-mark alias.
      type(parser_state), intent(inout) :: state !! Parser cursor.

      matched = match_keyword(state, "not")
      if (.not. matched) matched = match_text(state, "!")
   end function match_not_operator

   subroutine require_text(state, token)
      !! Consumes required punctuation or raises a parse error.
      type(parser_state), intent(inout) :: state !! Parser cursor.
      character(len=*), intent(in) :: token      !! Required token.

      if (.not. match_text(state, token)) error stop "expected token: " // token
   end subroutine require_text

   pure elemental subroutine skip_spaces(state)
      !! Advances over ASCII spaces and horizontal tabs.
      type(parser_state), intent(inout) :: state !! Parser cursor.

      do while (state%position <= len(state%source))
         if (state%source(state%position:state%position) /= " " .and. &
             state%source(state%position:state%position) /= achar(9)) exit
         state%position = state%position + 1
      end do
   end subroutine skip_spaces

   pure elemental logical function is_digit(character) result(answer)
      !! Reports whether a character is an ASCII decimal digit.
      character(len=1), intent(in) :: character !! Character to inspect.
      answer = character >= "0" .and. character <= "9"
   end function is_digit

   pure elemental logical function is_identifier_start(character) result(answer)
      !! Reports whether a character can begin an unquoted identifier.
      character(len=1), intent(in) :: character !! Character to inspect.
      answer = (character >= "a" .and. character <= "z") .or. &
         (character >= "A" .and. character <= "Z") .or. character == "_"
   end function is_identifier_start

   pure elemental logical function is_identifier_part(character) result(answer)
      !! Reports whether a character can occur after the first identifier character.
      character(len=1), intent(in) :: character !! Character to inspect.
      answer = is_identifier_start(character) .or. is_digit(character) .or. character == "."
   end function is_identifier_part

   pure elemental function lower_ascii(text) result(lowered)
      !! Converts ASCII uppercase letters to lowercase without locale dependence.
      character(len=*), intent(in) :: text !! Text to convert.
      character(len=len(text)) :: lowered
      integer :: code, i

      lowered = text
      do i = 1, len(text)
         code = iachar(lowered(i:i))
         if (code >= iachar("A") .and. code <= iachar("Z")) lowered(i:i) = achar(code + 32)
      end do
   end function lower_ascii

end module dplyr_expression
