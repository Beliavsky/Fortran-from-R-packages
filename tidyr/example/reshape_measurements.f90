! SPDX-License-Identifier: MIT
! SPDX-FileComment: Measurement reshape example for the Fortran tidyr translation.
program reshape_measurements
   use tidyr, only: dp, make_column, new_tibble, pivot_longer, tibble_column, tibble_type
   implicit none

   type(tibble_column) :: columns(3)
   type(tibble_type) :: long, wide

   columns(1) = make_column("station", [1, 2])
   columns(2) = make_column("morning", [12.5_dp, 14.0_dp])
   columns(3) = make_column("evening", [10.0_dp, 11.5_dp])
   wide = new_tibble(columns)
   long = pivot_longer(wide, [character(len=7) :: "morning", "evening"], &
      names_to="period", values_to="temperature")

   print "(a,i0,a,i0)", "Reshaped ", wide%nrow(), " rows to ", long%nrow()
   print "(*(g0,1x))", long%columns(3)%real_values
end program reshape_measurements
