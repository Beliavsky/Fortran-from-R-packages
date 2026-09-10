! SPDX-License-Identifier: MIT
! SPDX-FileComment: Deterministic tests for the Fortran purrr translation.
program test_purrr
   !! Verifies typed maps, folds, filtering, predicates, and detection.
   use purrr
   implicit none
   integer, allocatable :: integers(:)
   real(dp), allocatable :: reals(:)

   integers = map_int([1, 2, 3], square_integer)
   if (any(integers /= [1, 4, 9])) error stop 'map_int failed'
   reals = map_dbl([1.0_dp, 2.0_dp], square_real)
   if (any(abs(reals - [1.0_dp, 4.0_dp]) > 1.0e-12_dp)) error stop 'map_dbl failed'
   integers = map2_int([1, 2], [3, 4], add_integer)
   if (any(integers /= [4, 6])) error stop 'map2_int failed'
   if (reduce([1, 2, 3, 4], add_integer) /= 10) error stop 'reduce failed'
   integers = accumulate([1, 2, 3, 4], add_integer)
   if (any(integers /= [1, 3, 6, 10])) error stop 'accumulate failed'
   integers = accumulate([1, 2, 3], add_integer, initial=10)
   if (any(integers /= [10, 11, 13, 16])) error stop 'initial accumulate failed'
   integers = keep([1, 2, 3, 4], is_even)
   if (any(integers /= [2, 4])) error stop 'keep failed'
   integers = discard([1, 2, 3, 4], is_even)
   if (any(integers /= [1, 3])) error stop 'discard failed'
   if (.not. some([1, 2, 3], is_even)) error stop 'some failed'
   if (every([1, 2, 3], is_even)) error stop 'every failed'
   if (.not. none([1, 3, 5], is_even)) error stop 'none failed'
   if (detect_index([1, 3, 4, 6], is_even) /= 3) error stop 'detect_index failed'
   if (detect([1, 3, 4, 6], is_even) /= 4) error stop 'detect failed'
   integers = head_while([2, 4, 5, 6], is_even)
   if (any(integers /= [2, 4])) error stop 'head_while failed'
   integers = tail_while([1, 2, 4, 6], is_even)
   if (any(integers /= [2, 4, 6])) error stop 'tail_while failed'
   print '(a)', 'All purrr tests passed.'

contains

   pure integer function square_integer(value)
      !! Squares one integer.
      integer, intent(in) :: value !! Value to square.
      square_integer = value * value
   end function square_integer

   pure real(dp) function square_real(value)
      !! Squares one real value.
      real(dp), intent(in) :: value !! Value to square.
      square_real = value * value
   end function square_real

   pure integer function add_integer(left, right)
      !! Adds two integers.
      integer, intent(in) :: left  !! Left addend.
      integer, intent(in) :: right !! Right addend.
      add_integer = left + right
   end function add_integer

   pure logical function is_even(value)
      !! Tests whether an integer is even.
      integer, intent(in) :: value !! Value to test.
      is_even = modulo(value, 2) == 0
   end function is_even
end program test_purrr
