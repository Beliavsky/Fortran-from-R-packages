! SPDX-License-Identifier: MIT
module dplyr_vector_ops
   use, intrinsic :: ieee_arithmetic, only: ieee_is_nan
   use tibble, only: tibble_type
   use vctrs, only: dp, new_vctr, vctr_type, vec_cast, vec_ptype_common, &
      vec_recycle_common, vctrs_character, vctrs_integer, vctrs_logical, vctrs_real
   implicit none
   private

   public :: col, column, is_missing, vector_abs
   public :: operator(+), operator(-), operator(*), operator(/)
   public :: operator(==), operator(/=), operator(<), operator(<=), operator(>), operator(>=)
   public :: operator(.and.), operator(.or.), operator(.not.)

   interface column
      module procedure named_vector
   end interface column

   interface operator(+)
      module procedure add_vectors, add_vector_integer, add_vector_real
   end interface

   interface operator(-)
      module procedure subtract_vectors, subtract_vector_integer, subtract_vector_real
      module procedure negate_vector
   end interface

   interface operator(*)
      module procedure multiply_vectors, multiply_vector_integer, multiply_vector_real
   end interface

   interface operator(/)
      module procedure divide_vectors, divide_vector_integer, divide_vector_real
   end interface

   interface operator(==)
      module procedure equal_vectors, equal_vector_integer, equal_vector_real
      module procedure equal_vector_character, equal_vector_logical
   end interface

   interface operator(/=)
      module procedure not_equal_vectors, not_equal_vector_integer, not_equal_vector_real
      module procedure not_equal_vector_character, not_equal_vector_logical
   end interface

   interface operator(<)
      module procedure less_vectors, less_vector_integer, less_vector_real, less_vector_character
   end interface

   interface operator(<=)
      module procedure less_equal_vectors, less_equal_vector_integer, less_equal_vector_real
      module procedure less_equal_vector_character
   end interface

   interface operator(>)
      module procedure greater_vectors, greater_vector_integer, greater_vector_real
      module procedure greater_vector_character
   end interface

   interface operator(>=)
      module procedure greater_equal_vectors, greater_equal_vector_integer
      module procedure greater_equal_vector_real
      module procedure greater_equal_vector_character
   end interface

   interface operator(.and.)
      module procedure and_vectors
   end interface

   interface operator(.or.)
      module procedure or_vectors
   end interface

   interface operator(.not.)
      module procedure not_vector
   end interface

contains

   pure elemental function col(data, name) result(vector)
      !! Returns an exact named column as a typed vector.
      type(tibble_type), intent(in) :: data !! Table containing the column.
      character(len=*), intent(in) :: name  !! Exact column name.
      type(vctr_type) :: vector
      integer :: index

      index = data%column_index(name)
      if (index == 0) error stop "col: column not found: " // name
      vector = data%columns(index)
   end function col

   pure elemental function named_vector(name, vector) result(out)
      !! Assigns a semantic column name to an expression result.
      character(len=*), intent(in) :: name !! Nonempty output name.
      type(vctr_type), intent(in) :: vector !! Vector expression result.
      type(vctr_type) :: out

      if (len_trim(name) == 0) error stop "column: name cannot be empty"
      out = vector
      out%name = trim(name)
   end function named_vector

   pure elemental function add_vectors(left, right) result(out)
      !! Adds two numeric vectors with size recycling and type promotion.
      type(vctr_type), intent(in) :: left  !! Left operand.
      type(vctr_type), intent(in) :: right !! Right operand.
      type(vctr_type) :: out

      out = arithmetic_vectors(left, right, "+")
   end function add_vectors

   pure elemental function subtract_vectors(left, right) result(out)
      !! Subtracts two numeric vectors with size recycling and type promotion.
      type(vctr_type), intent(in) :: left  !! Left operand.
      type(vctr_type), intent(in) :: right !! Right operand.
      type(vctr_type) :: out

      out = arithmetic_vectors(left, right, "-")
   end function subtract_vectors

   pure elemental function multiply_vectors(left, right) result(out)
      !! Multiplies two numeric vectors with size recycling and type promotion.
      type(vctr_type), intent(in) :: left  !! Left operand.
      type(vctr_type), intent(in) :: right !! Right operand.
      type(vctr_type) :: out

      out = arithmetic_vectors(left, right, "*")
   end function multiply_vectors

   pure elemental function divide_vectors(left, right) result(out)
      !! Divides two numeric vectors and returns double precision.
      type(vctr_type), intent(in) :: left  !! Left operand.
      type(vctr_type), intent(in) :: right !! Right operand.
      type(vctr_type) :: out

      out = arithmetic_vectors(left, right, "/")
   end function divide_vectors

   pure elemental function arithmetic_vectors(left, right, operation) result(out)
      !! Evaluates one binary arithmetic operation on recycled numeric vectors.
      type(vctr_type), intent(in) :: left  !! Left operand.
      type(vctr_type), intent(in) :: right !! Right operand.
      character(len=1), intent(in) :: operation !! Arithmetic operator.
      type(vctr_type) :: out
      type(vctr_type) :: inputs(2), values(2)
      integer :: result_type

      inputs = [left, right]
      result_type = vec_ptype_common(inputs)
      if (result_type /= vctrs_integer .and. result_type /= vctrs_real) then
         error stop "arithmetic: operands must be integer or real"
      end if
      if (operation == "/") result_type = vctrs_real
      values = vec_recycle_common(inputs)
      values(1) = vec_cast(values(1), result_type)
      values(2) = vec_cast(values(2), result_type)
      out = allocate_result(result_type, values(1)%size())
      out%missing = values(1)%missing .or. values(2)%missing
      if (result_type == vctrs_integer) then
         select case (operation)
         case ("+")
            out%integer_values = values(1)%integer_values + values(2)%integer_values
         case ("-")
            out%integer_values = values(1)%integer_values - values(2)%integer_values
         case ("*")
            out%integer_values = values(1)%integer_values * values(2)%integer_values
         end select
      else
         select case (operation)
         case ("+")
            out%real_values = values(1)%real_values + values(2)%real_values
         case ("-")
            out%real_values = values(1)%real_values - values(2)%real_values
         case ("*")
            out%real_values = values(1)%real_values * values(2)%real_values
         case ("/")
            out%real_values = values(1)%real_values / values(2)%real_values
         end select
      end if
   end function arithmetic_vectors

   pure elemental function negate_vector(vector) result(out)
      !! Negates an integer or real vector while retaining its missing mask.
      type(vctr_type), intent(in) :: vector !! Numeric operand.
      type(vctr_type) :: out

      out = vector
      out%name = ""
      select case (vector%type_code)
      case (vctrs_integer)
         out%integer_values = -out%integer_values
      case (vctrs_real)
         out%real_values = -out%real_values
      case default
         error stop "unary minus: operand must be integer or real"
      end select
   end function negate_vector

   pure elemental function compare_vectors(left, right, operation) result(out)
      !! Compares recycled compatible vectors and propagates missingness.
      type(vctr_type), intent(in) :: left  !! Left operand.
      type(vctr_type), intent(in) :: right !! Right operand.
      character(len=*), intent(in) :: operation !! Comparison operator.
      type(vctr_type) :: out
      type(vctr_type) :: inputs(2), values(2)
      integer :: i, result_type

      inputs = [left, right]
      result_type = vec_ptype_common(inputs)
      if (result_type == 0) error stop "comparison: operands have incompatible types"
      values = vec_recycle_common(inputs)
      values(1) = vec_cast(values(1), result_type)
      values(2) = vec_cast(values(2), result_type)
      out = allocate_result(vctrs_logical, values(1)%size())
      out%missing = values(1)%missing .or. values(2)%missing
      do i = 1, out%size()
         if (.not. out%missing(i)) out%logical_values(i) = compare_element(values, i, operation)
      end do
   end function compare_vectors

   pure logical function compare_element(values, index, operation) result(answer)
      !! Compares one pair of same-type nonmissing values.
      type(vctr_type), intent(in) :: values(2) !! Recycled and cast operands.
      integer, intent(in) :: index             !! Element position.
      character(len=*), intent(in) :: operation !! Comparison operator.

      select case (values(1)%type_code)
      case (vctrs_integer)
         answer = compare_integer(values(1)%integer_values(index), &
            values(2)%integer_values(index), operation)
      case (vctrs_real)
         answer = compare_real( &
            values(1)%real_values(index), values(2)%real_values(index), operation)
      case (vctrs_logical)
         answer = compare_logical(values(1)%logical_values(index), &
            values(2)%logical_values(index), operation)
      case (vctrs_character)
         answer = compare_character(values(1)%character_values(index), &
            values(2)%character_values(index), operation)
      case default
         error stop "comparison: unsupported operand type"
      end select
   end function compare_element

   pure elemental logical function compare_integer(left, right, operation) result(answer)
      !! Compares two integer scalars for a selected operation.
      integer, intent(in) :: left  !! Left scalar operand.
      integer, intent(in) :: right !! Right scalar operand.
      character(len=*), intent(in) :: operation !! Comparison operator.

      select case (operation)
      case ("==")
      answer = left == right
      case ("/=")
      answer = left /= right
      case ("<")
      answer = left < right
      case ("<=")
      answer = left <= right
      case (">")
      answer = left > right
      case (">=")
      answer = left >= right
      case default
      error stop "comparison: unknown operator"
      end select
   end function compare_integer

   pure elemental logical function compare_real(left, right, operation) result(answer)
      !! Compares two real scalars, with NaN comparisons returning false except inequality.
      real(dp), intent(in) :: left  !! Left scalar operand.
      real(dp), intent(in) :: right !! Right scalar operand.
      character(len=*), intent(in) :: operation !! Comparison operator.

      if (ieee_is_nan(left) .or. ieee_is_nan(right)) then
         answer = operation == "/="
         return
      end if
      select case (operation)
      case ("==")
      answer = left <= right .and. right <= left
      case ("/=")
      answer = left < right .or. right < left
      case ("<")
      answer = left < right
      case ("<=")
      answer = left <= right
      case (">")
      answer = left > right
      case (">=")
      answer = left >= right
      case default
      error stop "comparison: unknown operator"
      end select
   end function compare_real

   pure elemental logical function compare_logical(left, right, operation) result(answer)
      !! Compares logical scalars for equality or inequality.
      logical, intent(in) :: left  !! Left scalar operand.
      logical, intent(in) :: right !! Right scalar operand.
      character(len=*), intent(in) :: operation !! Comparison operator.

      select case (operation)
      case ("==")
      answer = left .eqv. right
      case ("/=")
      answer = left .neqv. right
      case default
      error stop "comparison: logical values support only == and /="
      end select
   end function compare_logical

   pure elemental logical function compare_character(left, right, operation) result(answer)
      !! Compares character scalars lexically using Fortran collation.
      character(len=*), intent(in) :: left  !! Left scalar operand.
      character(len=*), intent(in) :: right !! Right scalar operand.
      character(len=*), intent(in) :: operation !! Comparison operator.

      select case (operation)
      case ("==")
      answer = left == right
      case ("/=")
      answer = left /= right
      case ("<")
      answer = left < right
      case ("<=")
      answer = left <= right
      case (">")
      answer = left > right
      case (">=")
      answer = left >= right
      case default
      error stop "comparison: unknown operator"
      end select
   end function compare_character

   pure elemental function logical_vectors(left, right, operation) result(out)
      !! Applies three-valued AND or OR to recycled logical vectors.
      type(vctr_type), intent(in) :: left  !! Left logical operand.
      type(vctr_type), intent(in) :: right !! Right logical operand.
      character(len=*), intent(in) :: operation !! `and` or `or`.
      type(vctr_type) :: out
      type(vctr_type) :: inputs(2), values(2)
      integer :: i

      if (left%type_code /= vctrs_logical .or. right%type_code /= vctrs_logical) then
         error stop "logical operation: operands must be logical"
      end if
      inputs = [left, right]
      values = vec_recycle_common(inputs)
      out = allocate_result(vctrs_logical, values(1)%size())
      do i = 1, out%size()
         call logical_element(values(1), values(2), i, operation, &
            out%logical_values(i), out%missing(i))
      end do
   end function logical_vectors

   pure elemental subroutine logical_element(left, right, index, operation, value, missing)
      !! Evaluates one three-valued logical operation.
      type(vctr_type), intent(in) :: left  !! Left logical operand.
      type(vctr_type), intent(in) :: right !! Right logical operand.
      integer, intent(in) :: index              !! Element position.
      character(len=*), intent(in) :: operation !! `and` or `or`.
      logical, intent(out) :: value             !! Logical result storage.
      logical, intent(out) :: missing           !! Result missingness.

      value = .false.
      if (operation == "and") then
         if ((.not. left%missing(index) .and. .not. left%logical_values(index)) .or. &
             (.not. right%missing(index) .and. .not. right%logical_values(index))) then
            missing = .false.
         else if (left%missing(index) .or. right%missing(index)) then
            missing = .true.
         else
            value = .true.
            missing = .false.
         end if
      else
         if ((.not. left%missing(index) .and. left%logical_values(index)) .or. &
             (.not. right%missing(index) .and. right%logical_values(index))) then
            value = .true.
            missing = .false.
         else if (left%missing(index) .or. right%missing(index)) then
            missing = .true.
         else
            missing = .false.
         end if
      end if
   end subroutine logical_element

   pure elemental function and_vectors(left, right) result(out)
      !! Computes three-valued elementwise conjunction.
      type(vctr_type), intent(in) :: left  !! Left logical operand.
      type(vctr_type), intent(in) :: right !! Right logical operand.
      type(vctr_type) :: out
      out = logical_vectors(left, right, "and")
   end function and_vectors

   pure elemental function or_vectors(left, right) result(out)
      !! Computes three-valued elementwise disjunction.
      type(vctr_type), intent(in) :: left  !! Left logical operand.
      type(vctr_type), intent(in) :: right !! Right logical operand.
      type(vctr_type) :: out
      out = logical_vectors(left, right, "or")
   end function or_vectors

   pure elemental function not_vector(vector) result(out)
      !! Negates a logical vector while preserving missingness.
      type(vctr_type), intent(in) :: vector !! Logical operand.
      type(vctr_type) :: out

      if (vector%type_code /= vctrs_logical) error stop "not: operand must be logical"
      out = vector
      out%name = ""
      out%logical_values = .not. out%logical_values
   end function not_vector

   pure elemental function is_missing(vector) result(out)
      !! Returns a nonmissing logical vector identifying missing elements.
      type(vctr_type), intent(in) :: vector !! Vector to inspect.
      type(vctr_type) :: out

      out = new_vctr("", vector%missing)
   end function is_missing

   pure elemental function vector_abs(vector) result(out)
      !! Returns absolute values for an integer or real vector.
      type(vctr_type), intent(in) :: vector !! Numeric operand.
      type(vctr_type) :: out

      out = vector
      out%name = ""
      select case (vector%type_code)
      case (vctrs_integer)
         out%integer_values = abs(out%integer_values)
      case (vctrs_real)
         out%real_values = abs(out%real_values)
      case default
         error stop "abs: operand must be integer or real"
      end select
   end function vector_abs

   pure elemental function allocate_result(type_code, n) result(out)
      !! Allocates an unnamed vector result for a supported primitive type.
      integer, intent(in) :: type_code !! Requested vctrs type code.
      integer, intent(in) :: n         !! Nonnegative vector size.
      type(vctr_type) :: out

      out%name = ""
      out%type_code = type_code
      allocate (out%missing(n), source=.false.)
      select case (type_code)
      case (vctrs_integer)
         allocate (out%integer_values(n), source=0)
      case (vctrs_real)
         allocate (out%real_values(n), source=0.0_dp)
      case (vctrs_logical)
         allocate (out%logical_values(n), source=.false.)
      case default
         error stop "allocate_result: unsupported type"
      end select
   end function allocate_result

   pure elemental function scalar_integer(value) result(vector)
      !! Wraps an integer scalar as a size-one vector.
      integer, intent(in) :: value !! Scalar to wrap.
      type(vctr_type) :: vector
      vector = new_vctr("", value, 1)
   end function scalar_integer

   pure elemental function scalar_real(value) result(vector)
      !! Wraps a real scalar as a size-one vector.
      real(dp), intent(in) :: value !! Scalar to wrap.
      type(vctr_type) :: vector
      vector = new_vctr("", value, 1)
   end function scalar_real

   pure elemental function scalar_character(value) result(vector)
      !! Wraps a character scalar as a size-one vector.
      character(len=*), intent(in) :: value !! Scalar to wrap.
      type(vctr_type) :: vector
      vector = new_vctr("", value, 1)
   end function scalar_character

   pure elemental function scalar_logical(value) result(vector)
      !! Wraps a logical scalar as a size-one vector.
      logical, intent(in) :: value !! Scalar to wrap.
      type(vctr_type) :: vector
      vector = new_vctr("", value, 1)
   end function scalar_logical

   pure elemental function add_vector_integer(left, right) result(out)
      !! Adds an integer scalar to a numeric vector.
      type(vctr_type), intent(in) :: left !! Vector operand.
      integer, intent(in) :: right        !! Scalar operand.
      type(vctr_type) :: out
      out = add_vectors(left, scalar_integer(right))
   end function add_vector_integer

   pure elemental function add_vector_real(left, right) result(out)
      !! Adds a real scalar to a numeric vector.
      type(vctr_type), intent(in) :: left !! Vector operand.
      real(dp), intent(in) :: right       !! Scalar operand.
      type(vctr_type) :: out
      out = add_vectors(left, scalar_real(right))
   end function add_vector_real

   pure elemental function subtract_vector_integer(left, right) result(out)
      !! Subtracts an integer scalar from a numeric vector.
      type(vctr_type), intent(in) :: left !! Vector operand.
      integer, intent(in) :: right        !! Scalar operand.
      type(vctr_type) :: out
      out = subtract_vectors(left, scalar_integer(right))
   end function subtract_vector_integer

   pure elemental function subtract_vector_real(left, right) result(out)
      !! Subtracts a real scalar from a numeric vector.
      type(vctr_type), intent(in) :: left !! Vector operand.
      real(dp), intent(in) :: right       !! Scalar operand.
      type(vctr_type) :: out
      out = subtract_vectors(left, scalar_real(right))
   end function subtract_vector_real

   pure elemental function multiply_vector_integer(left, right) result(out)
      !! Multiplies a numeric vector by an integer scalar.
      type(vctr_type), intent(in) :: left !! Vector operand.
      integer, intent(in) :: right        !! Scalar operand.
      type(vctr_type) :: out
      out = multiply_vectors(left, scalar_integer(right))
   end function multiply_vector_integer

   pure elemental function multiply_vector_real(left, right) result(out)
      !! Multiplies a numeric vector by a real scalar.
      type(vctr_type), intent(in) :: left !! Vector operand.
      real(dp), intent(in) :: right       !! Scalar operand.
      type(vctr_type) :: out
      out = multiply_vectors(left, scalar_real(right))
   end function multiply_vector_real

   pure elemental function divide_vector_integer(left, right) result(out)
      !! Divides a numeric vector by an integer scalar.
      type(vctr_type), intent(in) :: left !! Vector operand.
      integer, intent(in) :: right        !! Scalar operand.
      type(vctr_type) :: out
      out = divide_vectors(left, scalar_integer(right))
   end function divide_vector_integer

   pure elemental function divide_vector_real(left, right) result(out)
      !! Divides a numeric vector by a real scalar.
      type(vctr_type), intent(in) :: left !! Vector operand.
      real(dp), intent(in) :: right       !! Scalar operand.
      type(vctr_type) :: out
      out = divide_vectors(left, scalar_real(right))
   end function divide_vector_real

   pure elemental function equal_vectors(left, right) result(out)
      !! Tests recycled compatible vectors for equality.
      type(vctr_type), intent(in) :: left  !! Left operand.
      type(vctr_type), intent(in) :: right !! Right operand.
      type(vctr_type) :: out
      out = compare_vectors(left, right, "==")
   end function equal_vectors

   pure elemental function not_equal_vectors(left, right) result(out)
      !! Tests recycled compatible vectors for inequality.
      type(vctr_type), intent(in) :: left  !! Left operand.
      type(vctr_type), intent(in) :: right !! Right operand.
      type(vctr_type) :: out
      out = compare_vectors(left, right, "/=")
   end function not_equal_vectors

   pure elemental function less_vectors(left, right) result(out)
      !! Tests whether left vector values precede right values.
      type(vctr_type), intent(in) :: left  !! Left operand.
      type(vctr_type), intent(in) :: right !! Right operand.
      type(vctr_type) :: out
      out = compare_vectors(left, right, "<")
   end function less_vectors

   pure elemental function less_equal_vectors(left, right) result(out)
      !! Tests whether left vector values do not exceed right values.
      type(vctr_type), intent(in) :: left  !! Left operand.
      type(vctr_type), intent(in) :: right !! Right operand.
      type(vctr_type) :: out
      out = compare_vectors(left, right, "<=")
   end function less_equal_vectors

   pure elemental function greater_vectors(left, right) result(out)
      !! Tests whether left vector values exceed right values.
      type(vctr_type), intent(in) :: left  !! Left operand.
      type(vctr_type), intent(in) :: right !! Right operand.
      type(vctr_type) :: out
      out = compare_vectors(left, right, ">")
   end function greater_vectors

   pure elemental function greater_equal_vectors(left, right) result(out)
      !! Tests whether left vector values do not precede right values.
      type(vctr_type), intent(in) :: left  !! Left operand.
      type(vctr_type), intent(in) :: right !! Right operand.
      type(vctr_type) :: out
      out = compare_vectors(left, right, ">=")
   end function greater_equal_vectors

   pure elemental function equal_vector_integer(left, right) result(out)
      !! Tests a vector for equality with an integer scalar.
      type(vctr_type), intent(in) :: left !! Vector operand.
      integer, intent(in) :: right        !! Integer scalar operand.
      type(vctr_type) :: out
      out = equal_vectors(left, scalar_integer(right))
   end function equal_vector_integer

   pure elemental function equal_vector_real(left, right) result(out)
      !! Tests a vector for equality with a real scalar.
      type(vctr_type), intent(in) :: left !! Vector operand.
      real(dp), intent(in) :: right       !! Real scalar operand.
      type(vctr_type) :: out
      out = equal_vectors(left, scalar_real(right))
   end function equal_vector_real

   pure elemental function equal_vector_character(left, right) result(out)
      !! Tests a vector for equality with a character scalar.
      type(vctr_type), intent(in) :: left       !! Vector operand.
      character(len=*), intent(in) :: right     !! Character scalar operand.
      type(vctr_type) :: out
      out = equal_vectors(left, scalar_character(right))
   end function equal_vector_character

   pure elemental function equal_vector_logical(left, right) result(out)
      !! Tests a vector for equality with a logical scalar.
      type(vctr_type), intent(in) :: left !! Vector operand.
      logical, intent(in) :: right        !! Logical scalar operand.
      type(vctr_type) :: out
      out = equal_vectors(left, scalar_logical(right))
   end function equal_vector_logical

   pure elemental function not_equal_vector_integer(left, right) result(out)
      !! Tests a vector for inequality with an integer scalar.
      type(vctr_type), intent(in) :: left !! Vector operand.
      integer, intent(in) :: right        !! Integer scalar operand.
      type(vctr_type) :: out
      out = not_equal_vectors(left, scalar_integer(right))
   end function not_equal_vector_integer

   pure elemental function not_equal_vector_real(left, right) result(out)
      !! Tests a vector for inequality with a real scalar.
      type(vctr_type), intent(in) :: left !! Vector operand.
      real(dp), intent(in) :: right       !! Real scalar operand.
      type(vctr_type) :: out
      out = not_equal_vectors(left, scalar_real(right))
   end function not_equal_vector_real

   pure elemental function not_equal_vector_character(left, right) result(out)
      !! Tests a vector for inequality with a character scalar.
      type(vctr_type), intent(in) :: left   !! Vector operand.
      character(len=*), intent(in) :: right !! Character scalar operand.
      type(vctr_type) :: out
      out = not_equal_vectors(left, scalar_character(right))
   end function not_equal_vector_character

   pure elemental function not_equal_vector_logical(left, right) result(out)
      !! Tests a vector for inequality with a logical scalar.
      type(vctr_type), intent(in) :: left !! Vector operand.
      logical, intent(in) :: right        !! Logical scalar operand.
      type(vctr_type) :: out
      out = not_equal_vectors(left, scalar_logical(right))
   end function not_equal_vector_logical

   pure elemental function less_vector_integer(left, right) result(out)
      !! Tests whether vector values are less than an integer scalar.
      type(vctr_type), intent(in) :: left !! Vector operand.
      integer, intent(in) :: right        !! Integer scalar operand.
      type(vctr_type) :: out
      out = less_vectors(left, scalar_integer(right))
   end function less_vector_integer

   pure elemental function less_vector_real(left, right) result(out)
      !! Tests whether vector values are less than a real scalar.
      type(vctr_type), intent(in) :: left !! Vector operand.
      real(dp), intent(in) :: right       !! Real scalar operand.
      type(vctr_type) :: out
      out = less_vectors(left, scalar_real(right))
   end function less_vector_real

   pure elemental function less_vector_character(left, right) result(out)
      !! Tests whether vector values are less than a character scalar.
      type(vctr_type), intent(in) :: left   !! Vector operand.
      character(len=*), intent(in) :: right !! Character scalar operand.
      type(vctr_type) :: out
      out = less_vectors(left, scalar_character(right))
   end function less_vector_character

   pure elemental function less_equal_vector_integer(left, right) result(out)
      !! Tests whether vector values are at most an integer scalar.
      type(vctr_type), intent(in) :: left !! Vector operand.
      integer, intent(in) :: right        !! Integer scalar operand.
      type(vctr_type) :: out
      out = less_equal_vectors(left, scalar_integer(right))
   end function less_equal_vector_integer

   pure elemental function less_equal_vector_real(left, right) result(out)
      !! Tests whether vector values are at most a real scalar.
      type(vctr_type), intent(in) :: left !! Vector operand.
      real(dp), intent(in) :: right       !! Real scalar operand.
      type(vctr_type) :: out
      out = less_equal_vectors(left, scalar_real(right))
   end function less_equal_vector_real

   pure elemental function less_equal_vector_character(left, right) result(out)
      !! Tests whether vector values are at most a character scalar.
      type(vctr_type), intent(in) :: left   !! Vector operand.
      character(len=*), intent(in) :: right !! Character scalar operand.
      type(vctr_type) :: out
      out = less_equal_vectors(left, scalar_character(right))
   end function less_equal_vector_character

   pure elemental function greater_vector_integer(left, right) result(out)
      !! Tests whether vector values exceed an integer scalar.
      type(vctr_type), intent(in) :: left !! Vector operand.
      integer, intent(in) :: right        !! Integer scalar operand.
      type(vctr_type) :: out
      out = greater_vectors(left, scalar_integer(right))
   end function greater_vector_integer

   pure elemental function greater_vector_real(left, right) result(out)
      !! Tests whether vector values exceed a real scalar.
      type(vctr_type), intent(in) :: left !! Vector operand.
      real(dp), intent(in) :: right       !! Real scalar operand.
      type(vctr_type) :: out
      out = greater_vectors(left, scalar_real(right))
   end function greater_vector_real

   pure elemental function greater_vector_character(left, right) result(out)
      !! Tests whether vector values exceed a character scalar.
      type(vctr_type), intent(in) :: left   !! Vector operand.
      character(len=*), intent(in) :: right !! Character scalar operand.
      type(vctr_type) :: out
      out = greater_vectors(left, scalar_character(right))
   end function greater_vector_character

   pure elemental function greater_equal_vector_integer(left, right) result(out)
      !! Tests whether vector values are at least an integer scalar.
      type(vctr_type), intent(in) :: left !! Vector operand.
      integer, intent(in) :: right        !! Integer scalar operand.
      type(vctr_type) :: out
      out = greater_equal_vectors(left, scalar_integer(right))
   end function greater_equal_vector_integer

   pure elemental function greater_equal_vector_real(left, right) result(out)
      !! Tests whether vector values are at least a real scalar.
      type(vctr_type), intent(in) :: left !! Vector operand.
      real(dp), intent(in) :: right       !! Real scalar operand.
      type(vctr_type) :: out
      out = greater_equal_vectors(left, scalar_real(right))
   end function greater_equal_vector_real

   pure elemental function greater_equal_vector_character(left, right) result(out)
      !! Tests whether vector values are at least a character scalar.
      type(vctr_type), intent(in) :: left   !! Vector operand.
      character(len=*), intent(in) :: right !! Character scalar operand.
      type(vctr_type) :: out
      out = greater_equal_vectors(left, scalar_character(right))
   end function greater_equal_vector_character

end module dplyr_vector_ops
