! SPDX-License-Identifier: MIT
! SPDX-FileComment: Shared assertion helpers for tests of the Fortran stats translation.
module stats_test_assertions
   use iso_fortran_env, only: real64
   implicit none
   private

   integer, parameter :: dp = real64
   real(dp), parameter :: default_tolerance = 1.0e-10_dp

   public :: assert_approx, assert_close, assert_close_tolerance, assert_equal
   public :: assert_integer, assert_integer_matrix, assert_integer_vector
   public :: assert_matrix_close, assert_true, assert_vector_approx
   public :: assert_vector_close, assert_vector_close_tolerance

   interface assert_close
      module procedure assert_close_label, assert_close_tolerance_first
   end interface assert_close

   interface assert_vector_close
      module procedure assert_vector_close_label, assert_vector_close_tolerance_first
   end interface assert_vector_close

   interface assert_matrix_close
      module procedure assert_matrix_close_label, assert_matrix_close_tolerance_first
   end interface assert_matrix_close

   interface assert_approx
      module procedure assert_close_relative
   end interface assert_approx

   interface assert_close_tolerance
      module procedure assert_close_tolerance_first
   end interface assert_close_tolerance

   interface assert_vector_approx
      module procedure assert_vector_close_relative
   end interface assert_vector_approx

   interface assert_vector_close_tolerance
      module procedure assert_vector_close_label
   end interface assert_vector_close_tolerance

   interface assert_equal
      module procedure assert_integer, assert_integer_matrix, assert_integer_vector
   end interface assert_equal

contains

   subroutine assert_true(condition, label)
      !! Stops the test unless a condition is true.
      logical, intent(in) :: condition !! Condition expected to be true.
      character(len=*), intent(in) :: label !! Description of the assertion.

      if (.not. condition) error stop label
   end subroutine assert_true

   subroutine assert_close_label(actual, expected, label, allowed_error)
      !! Stops the test when two real scalars differ beyond a tolerance.
      real(dp), intent(in) :: actual !! Computed value.
      real(dp), intent(in) :: expected !! Expected value.
      character(len=*), intent(in) :: label !! Description of the assertion.
      real(dp), intent(in), optional :: allowed_error !! Absolute tolerance.
      real(dp) :: tolerance

      tolerance = default_tolerance
      if (present(allowed_error)) tolerance = allowed_error
      if (abs(actual - expected) > tolerance) then
         print *, label, actual, expected, abs(actual - expected)
         error stop label
      end if
   end subroutine assert_close_label

   subroutine assert_close_tolerance_first(actual_value, expected_value, tolerance_value, assertion_label)
      !! Supports the legacy assertion argument order with tolerance before label.
      real(dp), intent(in) :: actual_value !! Computed value.
      real(dp), intent(in) :: expected_value !! Expected value.
      real(dp), intent(in) :: tolerance_value !! Absolute tolerance.
      character(len=*), intent(in) :: assertion_label !! Description of the assertion.

      call assert_close_label(actual_value, expected_value, assertion_label, tolerance_value)
   end subroutine assert_close_tolerance_first

   subroutine assert_close_relative(actual, expected, allowed_error, label)
      !! Stops the test when scalar error exceeds a magnitude-scaled tolerance.
      real(dp), intent(in) :: actual !! Computed value.
      real(dp), intent(in) :: expected !! Expected value.
      real(dp), intent(in) :: allowed_error !! Relative-plus-absolute tolerance.
      character(len=*), intent(in) :: label !! Description of the assertion.

      if (abs(actual - expected) > allowed_error*(1.0_dp + abs(expected))) then
         print *, label, actual, expected, abs(actual - expected)
         error stop label
      end if
   end subroutine assert_close_relative

   subroutine assert_vector_close_label(actual, expected, label, allowed_error)
      !! Stops the test when real vectors have different shapes or values.
      real(dp), intent(in) :: actual(:) !! Computed vector.
      real(dp), intent(in) :: expected(:) !! Expected vector.
      character(len=*), intent(in) :: label !! Description of the assertion.
      real(dp), intent(in), optional :: allowed_error !! Maximum absolute error.
      real(dp) :: tolerance

      if (size(actual) /= size(expected)) error stop label//" size"
      tolerance = default_tolerance
      if (present(allowed_error)) tolerance = allowed_error
      if (any(abs(actual - expected) > tolerance)) then
         print *, label, maxval(abs(actual - expected))
         error stop label
      end if
   end subroutine assert_vector_close_label

   subroutine assert_vector_close_tolerance_first(actual_value, expected_value, tolerance_value, &
                                                  assertion_label)
      !! Supports the legacy vector assertion order with tolerance before label.
      real(dp), intent(in) :: actual_value(:) !! Computed vector.
      real(dp), intent(in) :: expected_value(:) !! Expected vector.
      real(dp), intent(in) :: tolerance_value !! Maximum absolute error.
      character(len=*), intent(in) :: assertion_label !! Description of the assertion.

      call assert_vector_close_label(actual_value, expected_value, assertion_label, tolerance_value)
   end subroutine assert_vector_close_tolerance_first

   subroutine assert_vector_close_relative(actual, expected, allowed_error, label)
      !! Stops the test when vector error exceeds a magnitude-scaled tolerance.
      real(dp), intent(in) :: actual(:) !! Computed vector.
      real(dp), intent(in) :: expected(:) !! Expected vector.
      real(dp), intent(in) :: allowed_error !! Relative-plus-absolute tolerance.
      character(len=*), intent(in) :: label !! Description of the assertion.

      if (size(actual) /= size(expected)) error stop label//" size"
      if (any(abs(actual - expected) > allowed_error*(1.0_dp + abs(expected)))) then
         print *, label, maxval(abs(actual - expected))
         error stop label
      end if
   end subroutine assert_vector_close_relative

   subroutine assert_matrix_close_label(actual, expected, label, allowed_error)
      !! Stops the test when real matrices have different shapes or values.
      real(dp), intent(in) :: actual(:, :) !! Computed matrix.
      real(dp), intent(in) :: expected(:, :) !! Expected matrix.
      character(len=*), intent(in) :: label !! Description of the assertion.
      real(dp), intent(in), optional :: allowed_error !! Maximum absolute error.
      real(dp) :: tolerance

      if (any(shape(actual) /= shape(expected))) error stop label//" shape"
      tolerance = default_tolerance
      if (present(allowed_error)) tolerance = allowed_error
      if (any(abs(actual - expected) > tolerance)) then
         print *, label, maxval(abs(actual - expected))
         error stop label
      end if
   end subroutine assert_matrix_close_label

   subroutine assert_matrix_close_tolerance_first(actual_value, expected_value, tolerance_value, &
                                                  assertion_label)
      !! Supports the legacy matrix assertion order with tolerance before label.
      real(dp), intent(in) :: actual_value(:, :) !! Computed matrix.
      real(dp), intent(in) :: expected_value(:, :) !! Expected matrix.
      real(dp), intent(in) :: tolerance_value !! Maximum absolute error.
      character(len=*), intent(in) :: assertion_label !! Description of the assertion.

      call assert_matrix_close_label(actual_value, expected_value, assertion_label, tolerance_value)
   end subroutine assert_matrix_close_tolerance_first

   subroutine assert_integer(actual, expected, label)
      !! Stops the test when two integer scalars differ.
      integer, intent(in) :: actual !! Computed integer value.
      integer, intent(in) :: expected !! Expected integer value.
      character(len=*), intent(in) :: label !! Description of the assertion.

      if (actual /= expected) error stop label
   end subroutine assert_integer

   subroutine assert_integer_vector(actual, expected, label)
      !! Stops the test when integer vectors have different shapes or values.
      integer, intent(in) :: actual(:) !! Computed integer vector.
      integer, intent(in) :: expected(:) !! Expected integer vector.
      character(len=*), intent(in) :: label !! Description of the assertion.

      if (size(actual) /= size(expected)) error stop label//" size"
      if (any(actual /= expected)) error stop label
   end subroutine assert_integer_vector

   subroutine assert_integer_matrix(actual, expected, label)
      !! Stops the test when integer matrices have different shapes or values.
      integer, intent(in) :: actual(:, :) !! Computed integer matrix.
      integer, intent(in) :: expected(:, :) !! Expected integer matrix.
      character(len=*), intent(in) :: label !! Description of the assertion.

      if (any(shape(actual) /= shape(expected))) error stop label//" shape"
      if (any(actual /= expected)) error stop label
   end subroutine assert_integer_matrix

end module stats_test_assertions
