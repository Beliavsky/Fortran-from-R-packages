program filter_and_summarise
   use dplyr
   implicit none

   type(tibble_column) :: columns(3)
   type(tibble_type) :: data, result

   columns(1) = make_column("sector", ["energy", "energy", "tech  ", "tech  "])
   columns(2) = make_column("return", [0.03_dp, -0.01_dp, 0.02_dp, 0.04_dp])
   columns(3) = make_column("volume", [12, 20, 8, 15])
   data = new_tibble(columns)

   result = filter(data, "return > 0 and volume >= 10")
   print "(a,i0)", "Rows passing the expression: ", result%nrow()
   print "(*(g0,1x))", result%columns(2)%real_values
end program filter_and_summarise
