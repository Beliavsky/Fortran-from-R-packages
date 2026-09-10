! SPDX-License-Identifier: MIT
! SPDX-FileComment: Pivot tests for the Fortran tidyr translation.
program test_pivot
   use tidyr, only: dp, make_column, new_tibble, pivot_longer, pivot_wider, tibble_column, &
      tibble_type
   implicit none

   type(tibble_column) :: columns(3)
   type(tibble_type) :: data, long, slow, sparse, wide

   columns(1) = make_column("id", [1, 2])
   columns(2) = make_column("x", [10.0_dp, 20.0_dp])
   columns(3) = make_column("y", [11.0_dp, 21.0_dp])
   data = new_tibble(columns)

   long = pivot_longer(data, [character(len=1) :: "x", "y"])
   call check(long%nrow() == 4, "pivot_longer row count")
   call check(long%ncol() == 3, "pivot_longer column count")
   call check(all(long%columns(1)%integer_values == [1, 1, 2, 2]), "fast identifier order")
   call check(all(long%columns(2)%character_values == ["x", "y", "x", "y"]), &
      "fast name order")
   call check(all(abs(long%columns(3)%real_values - [10.0_dp, 11.0_dp, 20.0_dp, 21.0_dp]) &
      < 1.0e-12_dp), "fast value order")

   slow = pivot_longer(data, [character(len=1) :: "x", "y"], cols_vary="slowest")
   call check(all(slow%columns(1)%integer_values == [1, 2, 1, 2]), "slow identifier order")
   call check(all(slow%columns(2)%character_values == ["x", "x", "y", "y"]), &
      "slow name order")
   call check(all(abs(slow%columns(3)%real_values - [10.0_dp, 20.0_dp, 11.0_dp, 21.0_dp]) &
      < 1.0e-12_dp), "slow value order")

   wide = pivot_wider(long, "name", "value", [character(len=2) :: "id"])
   call check(wide%nrow() == 2 .and. wide%ncol() == 3, "pivot_wider shape")
   call check(all(wide%columns(1)%integer_values == [1, 2]), "wide identifiers")
   call check(all(abs(wide%columns(2)%real_values - [10.0_dp, 20.0_dp]) < 1.0e-12_dp), &
      "wide x values")
   call check(all(abs(wide%columns(3)%real_values - [11.0_dp, 21.0_dp]) < 1.0e-12_dp), &
      "wide y values")
   call check(.not. any(wide%columns(2)%missing) .and. .not. any(wide%columns(3)%missing), &
      "wide complete mask")

   columns(1) = make_column("id", [1, 1, 2])
   columns(2) = make_column("name", ["x", "y", "x"])
   columns(3) = make_column("value", [10.0_dp, 11.0_dp, 20.0_dp])
   sparse = pivot_wider(new_tibble(columns), "name", "value", [character(len=2) :: "id"])
   call check(sparse%nrow() == 2 .and. sparse%ncol() == 3, "sparse wide shape")
   call check(.not. sparse%columns(2)%missing(2), "sparse x is present")
   call check(sparse%columns(3)%missing(2), "sparse y is missing")

contains

   subroutine check(condition, message)
      !! Stops the test with a focused diagnostic when an invariant fails.
      logical, intent(in) :: condition    !! Assertion result.
      character(len=*), intent(in) :: message !! Assertion label.

      if (.not. condition) error stop "test_pivot: " // message
   end subroutine check

end program test_pivot
