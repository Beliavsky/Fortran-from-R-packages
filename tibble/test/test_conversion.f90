program test_conversion
   use tibble, only: as_tibble, deframe_real, dp, enframe, get_character, get_integer, &
      make_column, new_tibble, tibble_column, tibble_type, validate_tibble
   implicit none

   type(tibble_type) :: table, framed, repaired, empty
   type(tibble_column) :: duplicate_columns(4)
   character(len=:), allocatable :: names(:), message
   character(len=3) :: character_matrix(2, 2)
   real(dp), allocatable :: values(:)
   integer :: matrix(2, 2)

   matrix = reshape([1, 2, 10, 20], [2, 2])
   table = as_tibble([character(len=5) :: "id", "value"], matrix)
   if (any(get_integer(table, "id") /= [1, 2])) error stop "matrix conversion failed"
   character_matrix = reshape([character(len=3) :: "a", "b", "c", "d"], [2, 2])
   table = as_tibble([character(len=5) :: "left", "right"], character_matrix)
   if (any(get_character(table, "right") /= [character(len=3) :: "c", "d"])) &
      error stop "character matrix conversion failed"

   framed = enframe([3.0_dp, 2.0_dp, 1.0_dp], &
      [character(len=5) :: "alpha", "beta", "gamma"], name="term", value="estimate")
   if (any(get_character(framed, "term") /= [character(len=5) :: "alpha", "beta", "gamma"])) &
      error stop "enframe names failed"
   values = deframe_real(framed, names)
   if (any(abs(values - [3.0_dp, 2.0_dp, 1.0_dp]) > 1.0e-12_dp)) error stop "deframe values failed"
   if (any(names /= [character(len=5) :: "alpha", "beta", "gamma"])) then
      error stop "deframe names failed"
   end if

   framed = enframe([8, 9])
   if (any(get_integer(framed, "name") /= [1, 2])) error stop "unnamed enframe sequence failed"

   duplicate_columns(1) = make_column("x...3", [1, 2])
   duplicate_columns(2) = make_column("x", [3, 4])
   duplicate_columns(3) = make_column("x", [5, 6])
   duplicate_columns(4) = make_column("", [7, 8])
   repaired = new_tibble(duplicate_columns, name_repair="unique")
   if (repaired%columns(1)%name /= "x...3") error stop "first unique name changed"
   if (repaired%columns(2)%name /= "x") error stop "unique name changed"
   if (repaired%columns(3)%name /= "x...4") error stop "conflicting suffix repair failed"
   if (repaired%columns(4)%name /= "...") error stop "empty name repair failed"
   if (.not. validate_tibble(repaired, message)) error stop "repaired table invalid: " // message

   empty = new_tibble([tibble_column ::], nrow=4)
   if (empty%nrow() /= 4 .or. empty%ncol() /= 0) error stop "zero-column table dimensions failed"
end program test_conversion
