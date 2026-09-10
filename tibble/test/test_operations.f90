program test_operations
   use tibble, only: add_column, add_rows, dp, drop_columns, filter_rows, get_character, &
      get_integer, get_real, make_column, new_tibble, replace_column, select_columns, &
      slice_rows, tibble_column, tibble_real, tibble_type
   implicit none

   type(tibble_column) :: columns(3), inserted_columns(3)
   type(tibble_type) :: table, selected, filtered, sliced, changed, inserted
   character(len=:), allocatable :: labels(:)
   integer, allocatable :: identifiers(:)
   real(dp), allocatable :: weights(:)

   columns(1) = make_column("id", [1, 2, 3])
   columns(2) = make_column("label", [character(len=5) :: "one", "two", "three"])
   columns(3) = make_column("weight", [1.0_dp, 2.0_dp, 3.0_dp])
   table = new_tibble(columns)

   filtered = filter_rows(table, [.true., .false., .true.])
   if (filtered%nrow() /= 2) error stop "filter row count failed"
   if (any(get_integer(filtered, "id") /= [1, 3])) error stop "filter values failed"

   sliced = slice_rows(table, [3, 1, 1])
   if (any(get_integer(sliced, "id") /= [3, 1, 1])) error stop "row slicing failed"

   selected = select_columns(table, [character(len=6) :: "weight", "label"])
   if (selected%ncol() /= 2) error stop "selection count failed"
   if (selected%columns(1)%name /= "weight") error stop "selection order failed"
   selected = drop_columns(table, [character(len=6) :: "weight"])
   if (selected%ncol() /= 2 .or. selected%has_name("weight")) error stop "column drop failed"

   changed = add_column(table, make_column("constant", 7, table%nrow()), before=2)
   if (changed%columns(2)%name /= "constant") error stop "column insertion position failed"
   if (any(get_integer(changed, "constant") /= 7)) error stop "scalar recycling failed"
   changed = replace_column(changed, make_column("constant", [0.5_dp, 1.5_dp, 2.5_dp]))
   if (any(abs(get_real(changed, "constant") - [0.5_dp, 1.5_dp, 2.5_dp]) > 1.0e-12_dp)) &
      error stop "typed column replacement failed"

   inserted_columns(1) = make_column("id", [9])
   inserted_columns(2) = make_column("label", [character(len=11) :: "nine-longer"])
   inserted_columns(3) = make_column("weight", [9.0_dp])
   inserted = add_rows(table, new_tibble(inserted_columns), before=2)
   identifiers = get_integer(inserted, "id")
   labels = get_character(inserted, "label")
   weights = get_real(inserted, "weight")
   if (any(identifiers /= [1, 9, 2, 3])) error stop "row insertion order failed"
   if (labels(2) /= "nine-longer") error stop "character width preservation failed"
   if (any(abs(weights - [1.0_dp, 9.0_dp, 2.0_dp, 3.0_dp]) > 1.0e-12_dp)) &
      error stop "row insertion values failed"

   inserted = add_rows(table, new_tibble(inserted_columns), before=1)
   if (any(get_integer(inserted, "id") /= [9, 1, 2, 3])) error stop "front row insertion failed"
   inserted = add_rows(table, new_tibble(inserted_columns))
   if (any(get_integer(inserted, "id") /= [1, 2, 3, 9])) error stop "appended row insertion failed"

   inserted_columns(1) = make_column("id", [9.5_dp])
   inserted = add_rows(table, new_tibble(inserted_columns))
   if (inserted%columns(1)%type_code /= tibble_real) error stop "row insertion promotion failed"
   weights = get_real(inserted, "id")
   if (abs(weights(4) - 9.5_dp) > 0.0_dp) error stop "promoted row insertion value failed"
end program test_operations
