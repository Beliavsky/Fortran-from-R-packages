program test_groups_helpers
   use dplyr
   implicit none

   type(tibble_column) :: columns(2)
   type(tibble_type) :: counts, data, result
   type(grouped_df) :: grouped
   type(summary_spec) :: specifications(3)
   type(vctr_type) :: values(2), vector
   integer, allocatable :: ids(:)

   columns(1) = make_column("g", ["a", "a", "b"])
   columns(2) = make_column("x", [1, 3, 4])
   data = new_tibble(columns)
   grouped = group_by(data, [character(len=1) :: "g"])
   call check(n_groups(grouped) == 2, "number of groups")
   call check(all(group_indices(grouped) == [1, 1, 2]), "group indices")
   call check(all(group_size(grouped) == [2, 1]), "group sizes")

   specifications(1) = summary("n", "", "n")
   specifications(2) = summary("total", "x", "sum")
   specifications(3) = summary("average", "x", "mean")
   result = summarise(grouped, specifications)
   call check(all(result%columns(2)%integer_values == [2, 1]), "summary counts")
   call check(all(result%columns(3)%integer_values == [4, 4]), "summary sums")
   call check(all(abs(result%columns(4)%real_values - [2.0_dp, 4.0_dp]) < 1.0e-12_dp), &
      "summary means")
   counts = count_rows(data, [character(len=1) :: "g"])
   call check(all(counts%columns(2)%integer_values == [2, 1]), "count rows")

   vector = lag(col(data, "x"))
   call check(vector%missing(1) .and. all(vector%integer_values(2:) == [1, 3]), "lag")
   vector = lead(col(data, "x"))
   call check(vector%missing(3) .and. all(vector%integer_values(:2) == [3, 4]), "lead")
   vector = between(col(data, "x"), 2.0_dp, 4.0_dp)
   call check(all(vector%logical_values .eqv. [.false., .true., .true.]), "between")
   vector = near(col(data, "x"), col(data, "x") + 1.0e-10_dp)
   call check(all(vector%logical_values), "near")
   call check(n_distinct(col(data, "g")) == 2, "n_distinct")

   values(1) = make_column("x", [1, 2, 3], [.false., .true., .false.])
   values(2) = make_column("x", [9])
   vector = coalesce(values)
   call check(all(vector%integer_values == [1, 9, 3]), "coalesce")
   vector = na_if(col(data, "x"), make_column("", 3, 1))
   call check(vector%missing(2), "na_if")
   ids = consecutive_id([make_column("", [1, 1, 2, 2]), make_column("", ["a", "b", "b", "b"])])
   call check(all(ids == [1, 2, 3, 3]), "consecutive_id")

contains

   subroutine check(condition, message)
      !! Stops the test with a focused diagnostic when an invariant fails.
      logical, intent(in) :: condition    !! Assertion result.
      character(len=*), intent(in) :: message !! Assertion label.
      if (.not. condition) error stop "test_groups_helpers: " // message
   end subroutine check

end program test_groups_helpers
