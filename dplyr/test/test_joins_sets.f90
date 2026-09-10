program test_joins_sets
   use dplyr
   implicit none

   type(tibble_column) :: x_columns(2), y_columns(2), one_column(1)
   type(tibble_type) :: joined, set_x, set_y, x, y

   x_columns(1) = make_column("id", [1, 2, 2])
   x_columns(2) = make_column("label", ["a", "b", "c"])
   y_columns(1) = make_column("id", [2, 2, 3])
   y_columns(2) = make_column("score", [10.0_dp, 20.0_dp, 30.0_dp])
   x = new_tibble(x_columns)
   y = new_tibble(y_columns)

   joined = inner_join(x, y, [character(len=2) :: "id"])
   call check(joined%nrow() == 4 .and. joined%ncol() == 3, "inner join shape")
   joined = left_join(x, y, [character(len=2) :: "id"])
   call check(joined%nrow() == 5 .and. joined%columns(3)%missing(1), "left join unmatched")
   joined = right_join(x, y, [character(len=2) :: "id"])
   call check(joined%nrow() == 5 .and. joined%columns(2)%missing(5), "right join unmatched")
   call check(joined%columns(1)%integer_values(5) == 3, "right join coalesced key")
   joined = full_join(x, y, [character(len=2) :: "id"])
   call check(joined%nrow() == 6, "full join shape")
   joined = semi_join(x, y, [character(len=2) :: "id"])
   call check(joined%nrow() == 2, "semi join")
   joined = anti_join(x, y, [character(len=2) :: "id"])
   call check(joined%nrow() == 1, "anti join")
   joined = cross_join(slice_head(x, 2), slice_head(y, 2))
   call check(joined%nrow() == 4, "cross join")

   one_column(1) = make_column("value", [1, 2, 2])
   set_x = new_tibble(one_column)
   one_column(1) = make_column("value", [2, 3])
   set_y = new_tibble(one_column)
   joined = union(set_x, set_y)
   call check(all(joined%columns(1)%integer_values == [1, 2, 3]), "union")
   joined = intersect(set_x, set_y)
   call check(all(joined%columns(1)%integer_values == [2]), "intersect")
   joined = setdiff(set_x, set_y)
   call check(all(joined%columns(1)%integer_values == [1]), "setdiff")
   joined = symdiff(set_x, set_y)
   call check(all(joined%columns(1)%integer_values == [1, 3]), "symdiff")
   call check(setequal(union(set_x, set_y), new_tibble([make_column("value", [3, 2, 1])])), &
      "setequal")

contains

   subroutine check(condition, message)
      !! Stops the test with a focused diagnostic when an invariant fails.
      logical, intent(in) :: condition    !! Assertion result.
      character(len=*), intent(in) :: message !! Assertion label.
      if (.not. condition) error stop "test_joins_sets: " // message
   end subroutine check

end program test_joins_sets
