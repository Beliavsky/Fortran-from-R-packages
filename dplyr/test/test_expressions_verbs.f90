program test_expressions_verbs
   use dplyr
   implicit none

   type(tibble_column) :: columns(3)
   type(tibble_type) :: data, result, typed
   type(vctr_type) :: computed

   columns(1) = make_column("x", [1.0_dp, 4.0_dp, 6.0_dp, 8.0_dp], &
      [.false., .false., .false., .true.])
   columns(2) = make_column("y", [5, 12, 8, 1])
   columns(3) = make_column("sector", ["energy", "tech  ", "energy", "energy"])
   data = new_tibble(columns)

   result = filter(data, "x > 3 and y <= 10")
   typed = filter(data, (col(data, "x") > 3.0_dp) .and. (col(data, "y") <= 10))
   call check(result%nrow() == 1, "string filter row count")
   call check(result%columns(2)%integer_values(1) == 8, "string filter value")
   call check(all(result%columns(2)%integer_values == typed%columns(2)%integer_values), &
      "typed and string filters agree")
   result = filter(data, "sector == 'energy'")
   call check(result%nrow() == 3, "character expression")
   result = filter(data, "missing(x)")
   call check(result%nrow() == 1, "missing expression")
   result = filter(data, "(x + y) >= 14 and not missing(x)")
   call check(result%nrow() == 2, "parenthesized arithmetic expression")
   result = filter(data, "abs(x - 5) <= 1")
   call check(result%nrow() == 2, "absolute-value expression")
   result = filter(data, "sector /= 'tech'")
   call check(result%nrow() == 3, "not-equal expression")
   result = filter(data, "sector != 'tech' & !missing(x)")
   call check(result%nrow() == 2, "R-style logical expression")
   result = filter(data, "true")
   call check(result%nrow() == 4, "scalar condition recycling")

   computed = column("z", col(data, "x") + col(data, "y") * 2)
   result = mutate(data, computed)
   call check(result%ncol() == 4, "mutate adds column")
   call check(abs(result%columns(4)%real_values(3) - 22.0_dp) < 1.0e-12_dp, &
      "mutate arithmetic")
   call check(result%columns(4)%missing(4), "mutate propagates missingness")

   result = arrange(data, [character(len=1) :: "y"], [.true.])
   call check(all(result%columns(2)%integer_values == [12, 8, 5, 1]), "arrange descending")
   result = distinct(data, [character(len=6) :: "sector"], keep_all=.true.)
   call check(result%nrow() == 2, "distinct selected key")
   result = relocate(data, [character(len=6) :: "sector"])
   call check(result%columns(1)%name == "sector", "relocate to front")
   result = rename(data, [character(len=6) :: "sector"], [character(len=5) :: "group"])
   call check(result%columns(3)%name == "group", "rename")
   result = select(data, [character(len=1) :: "y", "x"])
   call check(result%ncol() == 2 .and. result%columns(1)%name == "y", "select order")

contains

   subroutine check(condition, message)
      !! Stops the test with a focused diagnostic when an invariant fails.
      logical, intent(in) :: condition    !! Assertion result.
      character(len=*), intent(in) :: message !! Assertion label.
      if (.not. condition) error stop "test_expressions_verbs: " // message
   end subroutine check

end program test_expressions_verbs
