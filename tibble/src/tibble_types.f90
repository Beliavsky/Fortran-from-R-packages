! SPDX-License-Identifier: MIT
module tibble_types
   !! Re-exports vctrs storage using the established tibble-facing names.
   use vctrs, only: dp, make_column => new_vctr, &
      tibble_character => vctrs_character, tibble_column => vctr_type, &
      tibble_integer => vctrs_integer, tibble_logical => vctrs_logical, &
      tibble_real => vctrs_real, tibble_type => vctrs_data_frame
   implicit none
   private

   public :: dp, make_column, tibble_character, tibble_column, tibble_integer
   public :: tibble_logical, tibble_real, tibble_type
end module tibble_types
