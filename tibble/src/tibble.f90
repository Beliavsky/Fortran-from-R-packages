! SPDX-License-Identifier: MIT
module tibble
   !! Public umbrella module for the restricted modern-Fortran tibble API.
   use tibble_types, only: dp, make_column, tibble_character, tibble_column, &
      tibble_integer, tibble_logical, tibble_real, tibble_type
   use tibble_operations, only: add_column, add_rows, as_tibble, deframe_character, &
      deframe_integer, deframe_logical, deframe_real, drop_columns, enframe, &
      filter_rows, get_character, get_integer, get_logical, get_missing, get_real, glimpse, &
      new_tibble, print_tibble, repair_names, replace_column, select_columns, &
      slice_columns, slice_rows, validate_tibble
   implicit none
   public
end module tibble
