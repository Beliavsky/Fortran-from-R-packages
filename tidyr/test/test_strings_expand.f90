! SPDX-License-Identifier: MIT
! SPDX-FileComment: String and expansion tests for the Fortran tidyr translation.
program test_strings_expand
   use tidyr, only: crossing, expand_grid, make_column, new_tibble, separate_wider_delim, &
      tibble_column, tibble_type, unite
   implicit none

   type(tibble_column) :: columns(2), dimensions(2)
   type(tibble_type) :: data, grid, result

   columns(1) = make_column("id", [1, 2])
   columns(2) = make_column("code", ["a-b", "c-d"])
   data = new_tibble(columns)
   result = separate_wider_delim(data, "code", [character(len=5) :: "left", "right"], "-")
   call check(result%ncol() == 3, "separate column count")
   call check(all(result%columns(2)%character_values == ["a", "c"]), "separate left")
   call check(all(result%columns(3)%character_values == ["b", "d"]), "separate right")

   result = unite(result, "code", [character(len=5) :: "left", "right"], sep="-")
   call check(result%ncol() == 2, "unite column count")
   call check(all(result%columns(2)%character_values == ["a-b", "c-d"]), "unite values")

   dimensions(1) = make_column("x", [2, 1])
   dimensions(2) = make_column("y", ["b", "a"])
   grid = expand_grid(dimensions)
   call check(grid%nrow() == 4, "expand_grid row count")
   call check(all(grid%columns(1)%integer_values == [2, 2, 1, 1]), "expand_grid first")
   call check(all(grid%columns(2)%character_values == ["b", "a", "b", "a"]), &
      "expand_grid second")

   dimensions(1) = make_column("x", [2, 1, 2])
   dimensions(2) = make_column("y", ["b", "a"])
   grid = crossing(dimensions)
   call check(all(grid%columns(1)%integer_values == [1, 1, 2, 2]), "crossing first")
   call check(all(grid%columns(2)%character_values == ["a", "b", "a", "b"]), &
      "crossing second")

contains

   subroutine check(condition, message)
      !! Stops the test with a focused diagnostic when an invariant fails.
      logical, intent(in) :: condition    !! Assertion result.
      character(len=*), intent(in) :: message !! Assertion label.

      if (.not. condition) error stop "test_strings_expand: " // message
   end subroutine check

end program test_strings_expand
