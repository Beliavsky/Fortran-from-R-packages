module car_hypothesis
   use car_kinds, only: dp
   implicit none
   private

   public :: make_hypothesis

contains

   subroutine make_hypothesis(coefficient_names, hypotheses, hypothesis_matrix, rhs, status, rhs_values)
      character(len=*), intent(in) :: coefficient_names(:) !! Coefficient names allowed in the linear hypotheses.
      character(len=*), intent(in) :: hypotheses(:) !! Linear expressions using +, -, optional *, and at most one equals sign per row.
      real(dp), allocatable, intent(out) :: hypothesis_matrix(:, :) !! Numeric L matrix for hypotheses L beta = rhs.
      real(dp), allocatable, intent(out) :: rhs(:) !! Numeric right-hand side for each hypothesis.
      integer, intent(out) :: status !! Zero on success; nonzero for malformed expressions, unsupported names, or rhs shape errors.
      real(dp), intent(in), optional :: rhs_values(:) !! Optional explicit right-hand sides; length one is recycled, otherwise length must match hypotheses.
      real(dp), allocatable :: left_coefficients(:)
      real(dp), allocatable :: right_coefficients(:)
      real(dp) :: left_constant
      real(dp) :: right_constant
      character(len=:), allocatable :: expression
      character(len=:), allocatable :: left_text
      character(len=:), allocatable :: right_text
      integer :: equals_position
      integer :: i
      integer :: parse_status
      integer :: p

      p = size(coefficient_names)
      allocate(hypothesis_matrix(size(hypotheses), p), rhs(size(hypotheses)))
      hypothesis_matrix = 0.0_dp
      rhs = 0.0_dp
      status = 0
      if (p <= 0 .or. size(hypotheses) <= 0) then
         status = 1
         return
      end if
      if (any_unsupported_name(coefficient_names)) then
         status = 2
         return
      end if
      if (present(rhs_values)) then
         if (size(rhs_values) /= 1 .and. size(rhs_values) /= size(hypotheses)) then
            status = 3
            return
         end if
      end if

      allocate(left_coefficients(p), right_coefficients(p))
      do i = 1, size(hypotheses)
         expression = compact_text(hypotheses(i))
         if (len(expression) == 0) then
            status = 10 + i
            return
         end if
         equals_position = index(expression, '=')
         if (present(rhs_values)) then
            if (equals_position > 0) then
               status = 20 + i
               return
            end if
            left_text = expression
            right_text = ''
            call parse_linear_expression(left_text, coefficient_names, left_coefficients, left_constant, parse_status)
            if (parse_status /= 0) then
               status = 30 + i
               return
            end if
            hypothesis_matrix(i, :) = left_coefficients
            if (size(rhs_values) == 1) then
               rhs(i) = rhs_values(1) - left_constant
            else
               rhs(i) = rhs_values(i) - left_constant
            end if
         else
            if (equals_position > 0) then
               if (index(expression(equals_position + 1:), '=') > 0) then
                  status = 40 + i
                  return
               end if
               left_text = expression(:equals_position - 1)
               right_text = expression(equals_position + 1:)
            else
               left_text = expression
               right_text = ''
            end if
            call parse_linear_expression(left_text, coefficient_names, left_coefficients, left_constant, parse_status)
            if (parse_status /= 0) then
               status = 50 + i
               return
            end if
            call parse_linear_expression(right_text, coefficient_names, right_coefficients, right_constant, parse_status)
            if (parse_status /= 0) then
               status = 60 + i
               return
            end if
            hypothesis_matrix(i, :) = left_coefficients - right_coefficients
            rhs(i) = right_constant - left_constant
         end if
      end do
   end subroutine make_hypothesis

   subroutine parse_linear_expression(expression, names, coefficients, constant, status)
      character(len=*), intent(in) :: expression !! Compact linear expression on one side of a hypothesis equation.
      character(len=*), intent(in) :: names(:) !! Supported coefficient names used to identify symbolic terms.
      real(dp), intent(out) :: coefficients(:) !! Parsed coefficient multiplier for each name.
      real(dp), intent(out) :: constant !! Parsed numeric constant on this side of the equation.
      integer, intent(out) :: status !! Zero on success; nonzero for malformed or ambiguous terms.
      character(len=:), allocatable :: term
      real(dp) :: numeric_value
      integer :: end_position
      integer :: ios
      integer :: name_index
      integer :: position
      integer :: sign_value

      coefficients = 0.0_dp
      constant = 0.0_dp
      status = 0
      if (len_trim(expression) == 0) return

      position = 1
      do while (position <= len_trim(expression))
         sign_value = 1
         if (expression(position:position) == '+') then
            position = position + 1
         else if (expression(position:position) == '-') then
            sign_value = -1
            position = position + 1
         end if
         if (position > len_trim(expression)) then
            status = 1
            return
         end if
         end_position = position
         do while (end_position <= len_trim(expression))
            if (expression(end_position:end_position) == '+' .or. expression(end_position:end_position) == '-') exit
            end_position = end_position + 1
         end do
         term = expression(position:end_position - 1)
         if (len_trim(term) == 0) then
            status = 1
            return
         end if

         call identify_name_term(term, names, name_index, numeric_value, ios)
         if (ios == 0 .and. name_index > 0) then
            coefficients(name_index) = coefficients(name_index) + real(sign_value, dp) * numeric_value
         else if (ios == 0) then
            read(term, *, iostat=ios) numeric_value
            if (ios /= 0) then
               status = 2
               return
            end if
            constant = constant + real(sign_value, dp) * numeric_value
         else
            status = 3
            return
         end if
         position = end_position
      end do
   end subroutine parse_linear_expression

   subroutine identify_name_term(term, names, name_index, multiplier, status)
      character(len=*), intent(in) :: term !! Single unsigned additive term to classify as a coefficient term or constant.
      character(len=*), intent(in) :: names(:) !! Supported coefficient names, none containing +, -, *, or =.
      integer, intent(out) :: name_index !! One-based matched name index, or zero when the term is numeric-only.
      real(dp), intent(out) :: multiplier !! Numeric multiplier preceding a matched coefficient name.
      integer, intent(out) :: status !! Zero for a valid name or possible constant; nonzero for an invalid numeric prefix.
      character(len=:), allocatable :: prefix
      character(len=:), allocatable :: term_trimmed
      integer :: best_length
      integer :: i
      integer :: ios
      integer :: name_length
      integer :: term_length

      term_trimmed = trim(term)
      term_length = len_trim(term_trimmed)
      name_index = 0
      best_length = -1
      multiplier = 0.0_dp
      status = 0
      do i = 1, size(names)
         name_length = len_trim(names(i))
         if (name_length <= 0 .or. name_length > term_length) cycle
         if (term_trimmed(term_length - name_length + 1:term_length) == trim(names(i))) then
            if (name_length > best_length) then
               best_length = name_length
               name_index = i
            end if
         end if
      end do
      if (name_index == 0) return

      name_length = len_trim(names(name_index))
      if (term_length == name_length) then
         multiplier = 1.0_dp
         return
      end if
      prefix = term_trimmed(:term_length - name_length)
      if (len(prefix) > 0) then
         if (prefix(len(prefix):len(prefix)) == '*') prefix = prefix(:len(prefix) - 1)
      end if
      if (len_trim(prefix) == 0) then
         multiplier = 1.0_dp
         return
      end if
      read(prefix, *, iostat=ios) multiplier
      if (ios /= 0) status = 1
   end subroutine identify_name_term

   pure function compact_text(text) result(compact)
      character(len=*), intent(in) :: text !! Source hypothesis text from which spaces and tabs are removed.
      character(len=:), allocatable :: compact
      character(len=len(text)) :: work
      integer :: i
      integer :: j

      j = 0
      do i = 1, len_trim(text)
         if (text(i:i) /= ' ' .and. iachar(text(i:i)) /= 9) then
            j = j + 1
            work(j:j) = text(i:i)
         end if
      end do
      if (j == 0) then
         compact = ''
      else
         compact = work(:j)
      end if
   end function compact_text

   pure function any_unsupported_name(names) result(unsupported)
      character(len=*), intent(in) :: names(:) !! Coefficient names checked for operators unsupported by this partial parser.
      logical :: unsupported
      integer :: i

      unsupported = .false.
      do i = 1, size(names)
         if (len_trim(names(i)) == 0 .or. index(names(i), '+') > 0 .or. index(names(i), '-') > 0 .or. &
             index(names(i), '*') > 0 .or. index(names(i), '=') > 0) then
            unsupported = .true.
            return
         end if
      end do
   end function any_unsupported_name

end module car_hypothesis
