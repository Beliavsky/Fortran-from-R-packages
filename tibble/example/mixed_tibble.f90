program mixed_tibble
   use tibble, only: dp, glimpse, make_column, new_tibble, print_tibble, tibble_column, tibble_type
   implicit none

   type(tibble_column) :: columns(4)
   type(tibble_type) :: observations

   columns(1) = make_column("id", [1, 2, 3, 4])
   columns(2) = make_column("asset", [character(len=4) :: "SPY", "EFA", "IEF", "GLD"])
   columns(3) = make_column("return", [0.012_dp, -0.004_dp, 0.003_dp, 0.009_dp])
   columns(4) = make_column("selected", [.true., .false., .true., .true.])
   observations = new_tibble(columns)

   call print_tibble(observations)
   write (*, *)
   call glimpse(observations)
end program mixed_tibble
