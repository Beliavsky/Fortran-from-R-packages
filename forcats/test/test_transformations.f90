! SPDX-License-Identifier: MIT
! SPDX-FileComment: Transformation tests for the Fortran forcats translation.
program test_transformations
   !! Tests high-level deterministic factor transformations.
   use forcats
   implicit none

   type(factor_type) :: factor, output, right, crossed, combined
   type(factor_type) :: factors(2)
   logical, allocatable :: matched(:)

   factor = fct([character(len=1) :: "b", "b", "a", "c", "c", "c"], &
      [character(len=1) :: "a", "b", "c", "d"])

   output = fct_inorder(factor)
   call assert_levels(output, [character(len=1) :: "b", "a", "c", "d"], "inorder")
   output = fct_infreq(factor)
   call assert_levels(output, [character(len=1) :: "c", "b", "a", "d"], "infreq")
   output = fct_rev(factor)
   call assert_levels(output, [character(len=1) :: "d", "c", "b", "a"], "reverse")
   output = fct_shift(factor, -1)
   call assert_levels(output, [character(len=1) :: "d", "a", "b", "c"], "shift")
   output = fct_relevel(factor, [character(len=1) :: "c", "a"], after=1)
   call assert_levels(output, [character(len=1) :: "b", "c", "a", "d"], "relevel")

   output = fct_expand(factor, [character(len=1) :: "e", "b", "f"], after=1)
   call assert_levels(output, [character(len=1) :: "a", "e", "f", "b", "c", "d"], "expand")
   output = fct_drop(factor)
   call assert_levels(output, [character(len=1) :: "a", "b", "c"], "drop")

   output = fct_recode(factor, [character(len=1) :: "a", "c"], [character(len=5) :: "alpha", "alpha"])
   call assert_levels(output, [character(len=5) :: "alpha", "b", "d"], "recode collapse")
   output = fct_other(factor, [character(len=1) :: "b", "c"])
   call assert_levels(output, [character(len=5) :: "b", "c", "Other"], "other")

   matched = fct_match(factor, [character(len=1) :: "a", "c"])
   call assert_true(all(matched .eqv. [.false., .false., .true., .true., .true., .true.]), "match")

   output = fct_na_value_to_level(fct([character(len=1) :: "a", "x", "b"], &
      missing=[.false., .true., .false.]), "Missing")
   call assert_levels(output, [character(len=7) :: "a", "b", "Missing"], "missing to level")
   call assert_true(.not. any(output%missing), "missing resolved")
   output = fct_na_level_to_value(output, [character(len=7) :: "Missing"])
   call assert_levels(output, [character(len=1) :: "a", "b"], "level to missing")
   call assert_true(output%missing(2), "missing restored")

   factors(1) = fct([character(len=1) :: "a", "b"])
   factors(2) = fct([character(len=1) :: "b", "c"])
   combined = fct_c(factors)
   call assert_levels(combined, [character(len=1) :: "a", "b", "c"], "concatenate levels")
   call assert_integer(combined%codes, [1, 2, 2, 3], "concatenate codes")

   right = fct([character(len=5) :: "green", "green", "red", "green"])
   crossed = fct_cross(fct([character(len=5) :: "apple", "kiwi", "apple", "apple"]), right)
   call assert_levels(crossed, [character(len=11) :: "apple:green", "apple:red", "kiwi:green"], "cross")
   call assert_integer(crossed%codes, [1, 3, 2, 1], "cross codes")

contains

   pure elemental subroutine assert_true(condition, message)
      !! Stops the test when a logical condition is false.
      logical, intent(in) :: condition        !! Condition that must hold.
      character(len=*), intent(in) :: message !! Failure description.

      if (.not. condition) error stop "FAIL: " // message
   end subroutine assert_true

   pure subroutine assert_levels(actual, expected, message)
      !! Stops the test when factor levels differ.
      type(factor_type), intent(in) :: actual  !! Computed factor.
      character(len=*), intent(in) :: expected(:) !! Required levels.
      character(len=*), intent(in) :: message  !! Failure description.
      integer :: i

      call assert_true(actual%nlevels() == size(expected), message // " size")
      do i = 1, size(expected)
         call assert_true(trim(actual%levels(i)) == trim(expected(i)), message)
      end do
   end subroutine assert_levels

   pure subroutine assert_integer(actual, expected, message)
      !! Stops the test when integer arrays differ.
      integer, intent(in) :: actual(:)         !! Computed values.
      integer, intent(in) :: expected(:)       !! Required values.
      character(len=*), intent(in) :: message  !! Failure description.

      call assert_true(size(actual) == size(expected), message // " size")
      call assert_true(all(actual == expected), message)
   end subroutine assert_integer

end program test_transformations
