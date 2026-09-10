! SPDX-License-Identifier: MIT
! SPDX-FileComment: Public facade for the Fortran tidyr translation.
module tidyr
   !! Re-exports supported tidy-data reshaping operations.
   use tibble, only: as_tibble, dp, make_column, new_tibble, print_tibble
   use tibble, only: tibble_character, tibble_column, tibble_integer
   use tibble, only: tibble_logical, tibble_real, tibble_type
   use tidyr_expand
   use tidyr_missing
   use tidyr_pivot
   use tidyr_strings
   implicit none
   public
end module tidyr
