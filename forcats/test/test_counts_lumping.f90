! SPDX-License-Identifier: MIT
! SPDX-FileComment: Counting and lumping tests for the Fortran forcats translation.
program test_counts_lumping
   !! Tests counts, lumping, and data-driven reorder operations.
   use forcats
   implicit none

   type(factor_type) :: factor, output
   type(factor_count_type) :: table
   real(dp) :: x(9), y(9)

   factor = fct([character(len=1) :: "a", "a", "a", "a", "b", "b", "c", "d", "e"], &
      [character(len=1) :: "a", "b", "c", "d", "e", "z"])
   table = fct_count(factor, prop=.true.)
   call assert_integer(table%count, [4, 2, 1, 1, 1, 0], "counts")
   call assert_true(abs(sum(table%proportion) - 1.0_dp) < 1.0e-14_dp, "proportions")

   output = fct_lump_min(factor, 2.0_dp)
   call assert_levels(output, [character(len=5) :: "a", "b", "Other"], "lump minimum")
   call assert_integer(output%codes, [1, 1, 1, 1, 2, 2, 3, 3, 3], "lump minimum codes")
   output = fct_lump_n(factor, 2)
   call assert_levels(output, [character(len=5) :: "a", "b", "Other"], "lump n")
   output = fct_lump_prop(factor, 0.15_dp)
   call assert_levels(output, [character(len=5) :: "a", "b", "Other"], "lump proportion")

   x = [4.0_dp, 6.0_dp, 8.0_dp, 10.0_dp, 1.0_dp, 3.0_dp, 20.0_dp, 12.0_dp, 16.0_dp]
   output = fct_reorder(factor, x, statistic="mean")
   call assert_levels(output, [character(len=1) :: "b", "a", "d", "e", "c", "z"], "reorder mean")

   x = [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 1.0_dp, 2.0_dp, 1.0_dp, 1.0_dp, 1.0_dp]
   y = [2.0_dp, 4.0_dp, 8.0_dp, 9.0_dp, 7.0_dp, 6.0_dp, 5.0_dp, 3.0_dp, 1.0_dp]
   output = fct_reorder2(factor, x, y)
   call assert_levels(output, [character(len=1) :: "a", "b", "c", "d", "e", "z"], "reorder2")
   call assert_true(abs(last2([1.0_dp, 3.0_dp, 2.0_dp], [10.0_dp, 30.0_dp, 20.0_dp]) - 30.0_dp) < &
      1.0e-14_dp, "last2")
   call assert_true(abs(first2([1.0_dp, 3.0_dp, 2.0_dp], [10.0_dp, 30.0_dp, 20.0_dp]) - 10.0_dp) < &
      1.0e-14_dp, "first2")

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

end program test_counts_lumping
