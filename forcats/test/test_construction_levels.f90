! SPDX-License-Identifier: MIT
! SPDX-FileComment: Construction tests for the Fortran forcats translation.
program test_construction_levels
   !! Tests construction and low-level level transformations.
   use forcats
   implicit none

   type(factor_type) :: factor, output
   character(len=:), allocatable :: values(:)
   logical :: missing(6)

   missing = [.false., .false., .false., .true., .false., .false.]
   factor = fct([character(len=1) :: "b", "a", "b", "x", "c", "a"], missing=missing)
   call assert_character(factor%levels, [character(len=1) :: "b", "a", "c"], "first appearance")
   call assert_integer(factor%codes, [1, 2, 1, 0, 3, 2], "initial codes")
   call assert_true(all(factor%missing .eqv. missing), "initial missing mask")
   call assert_true(factor%valid(), "factor invariants")

   output = lvls_reorder(factor, [3, 1, 2], ordered=.true.)
   call assert_character(output%levels, [character(len=1) :: "c", "b", "a"], "reordered levels")
   call assert_integer(output%codes, [2, 3, 2, 0, 1, 3], "reordered codes")
   call assert_true(output%ordered, "ordered override")

   output = lvls_revalue(factor, [character(len=5) :: "other", "fruit", "other"])
   call assert_character(output%levels, [character(len=5) :: "other", "fruit"], "collapsed levels")
   call assert_integer(output%codes, [1, 2, 1, 0, 1, 2], "collapsed codes")

   output = lvls_expand(factor, [character(len=1) :: "z", "b", "a", "c", "q"])
   call assert_character(output%levels, [character(len=1) :: "z", "b", "a", "c", "q"], "expanded levels")
   call assert_integer(output%codes, [2, 3, 2, 0, 4, 3], "expanded codes")

   output = as_factor([3, 1, 2, 1])
   call assert_character(output%levels, [character(len=1) :: "1", "2", "3"], "integer levels")
   call assert_integer(output%codes, [3, 1, 2, 1], "integer codes")

   output = as_factor([.true., .false., .true.])
   call assert_character(output%levels, [character(len=5) :: "FALSE", "TRUE"], "logical levels")
   call assert_integer(output%codes, [2, 1, 2], "logical codes")

   values = factor%values("missing")
   call assert_character(values, [character(len=7) :: "b", "a", "b", "missing", "c", "a"], "decoded values")

contains

   pure elemental subroutine assert_true(condition, message)
      !! Stops the test when a logical condition is false.
      logical, intent(in) :: condition   !! Condition that must hold.
      character(len=*), intent(in) :: message !! Failure description.

      if (.not. condition) error stop "FAIL: " // message
   end subroutine assert_true

   pure subroutine assert_integer(actual, expected, message)
      !! Stops the test when integer arrays differ.
      integer, intent(in) :: actual(:)   !! Computed values.
      integer, intent(in) :: expected(:) !! Required values.
      character(len=*), intent(in) :: message !! Failure description.

      call assert_true(size(actual) == size(expected), message // " size")
      call assert_true(all(actual == expected), message)
   end subroutine assert_integer

   pure subroutine assert_character(actual, expected, message)
      !! Stops the test when character arrays differ after trimming.
      character(len=*), intent(in) :: actual(:)   !! Computed values.
      character(len=*), intent(in) :: expected(:) !! Required values.
      character(len=*), intent(in) :: message     !! Failure description.
      integer :: i

      call assert_true(size(actual) == size(expected), message // " size")
      do i = 1, size(actual)
         call assert_true(trim(actual(i)) == trim(expected(i)), message)
      end do
   end subroutine assert_character

end program test_construction_levels
