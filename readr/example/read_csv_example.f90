! SPDX-License-Identifier: MIT
! SPDX-FileComment: Example for the Fortran readr translation.
program read_csv_example
   !! Reads the bundled mixed-type CSV file and prints its inferred structure.
   use readr
   use tibble, only: print_tibble
   implicit none

   type(read_result_type) :: result

   result = read_csv("data/sample.csv")
   call stop_for_problems(result)
   call print_tibble(result%data)
end program read_csv_example
