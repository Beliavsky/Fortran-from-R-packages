! SPDX-License-Identifier: MIT
! SPDX-FileComment: Callback contracts for the Fortran purrr translation.
module purrr_callbacks
   !! Defines explicit procedure interfaces used by typed functional operations.
   use, intrinsic :: iso_fortran_env, only: real64
   implicit none
   private
   integer, parameter, public :: dp = real64

   public :: integer_binary, integer_predicate, integer_unary, integer_walker
   public :: real_binary, real_predicate, real_unary, real_walker

   abstract interface
      pure integer function integer_unary(value)
         !! Transforms one integer.
         integer, intent(in) :: value !! Input value.
      end function integer_unary

      pure integer function integer_binary(left, right)
         !! Combines two integers.
         integer, intent(in) :: left  !! Left input.
         integer, intent(in) :: right !! Right input.
      end function integer_binary

      pure logical function integer_predicate(value)
         !! Tests one integer.
         integer, intent(in) :: value !! Value to test.
      end function integer_predicate

      subroutine integer_walker(value)
         !! Visits one integer for side effects.
         integer, intent(in) :: value !! Value supplied for its side effect.
      end subroutine integer_walker

      pure real(dp) function real_unary(value)
         !! Transforms one real value.
         import dp
         real(dp), intent(in) :: value !! Input value.
      end function real_unary

      pure real(dp) function real_binary(left, right)
         !! Combines two real values.
         import dp
         real(dp), intent(in) :: left  !! Left input.
         real(dp), intent(in) :: right !! Right input.
      end function real_binary

      pure logical function real_predicate(value)
         !! Tests one real value.
         import dp
         real(dp), intent(in) :: value !! Value to test.
      end function real_predicate

      subroutine real_walker(value)
         !! Visits one real value for side effects.
         import dp
         real(dp), intent(in) :: value !! Value supplied for its side effect.
      end subroutine real_walker
   end interface
end module purrr_callbacks
