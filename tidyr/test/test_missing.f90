! SPDX-License-Identifier: MIT
! SPDX-FileComment: Missing-value tests for the Fortran tidyr translation.
program test_missing
   use tidyr, only: drop_na, fill, make_column, new_tibble, replace_na, tibble_column, &
      tibble_type, uncount
   implicit none

   type(tibble_column) :: columns(3), replacements(1)
   type(tibble_type) :: data, result

   columns(1) = make_column("a", [1, 2, 3], [.false., .true., .false.])
   columns(2) = make_column("b", ["x", "y", "z"], [.false., .false., .true.])
   columns(3) = make_column("w", [1, 0, 2])
   data = new_tibble(columns)

   result = drop_na(data)
   call check(result%nrow() == 1, "drop_na all")
   result = drop_na(data, [character(len=1) :: "a"])
   call check(result%nrow() == 2, "drop_na selected")
   call check(all(result%columns(1)%integer_values == [1, 3]), "drop_na row values")

   replacements(1) = make_column("a", 9, 1)
   result = replace_na(data, replacements)
   call check(all(result%columns(1)%integer_values == [1, 9, 3]), "replace_na values")
   call check(.not. any(result%columns(1)%missing), "replace_na mask")

   result = fill(data, [character(len=1) :: "a", "b"], direction="downup")
   call check(all(result%columns(1)%integer_values == [1, 1, 3]), "fill down integer")
   call check(all(result%columns(2)%character_values == ["x", "y", "y"]), &
      "fill down character")
   call check(.not. any(result%columns(1)%missing) .and. .not. any(result%columns(2)%missing), &
      "fill masks")

   result = uncount(data, "w")
   call check(result%nrow() == 3 .and. result%ncol() == 2, "uncount shape")
   call check(all(result%columns(1)%integer_values == [1, 3, 3]), "uncount rows")

contains

   subroutine check(condition, message)
      !! Stops the test with a focused diagnostic when an invariant fails.
      logical, intent(in) :: condition    !! Assertion result.
      character(len=*), intent(in) :: message !! Assertion label.

      if (.not. condition) error stop "test_missing: " // message
   end subroutine check

end program test_missing
