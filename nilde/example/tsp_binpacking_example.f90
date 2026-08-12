! SPDX-License-Identifier: GPL-2.0-or-later
program tsp_binpacking_example
   use nilde
   implicit none
   type(bin_packing_result_t) :: bp
   type(tsp_result_t) :: tsp
   integer(i8) :: c(4,4)

   bp = bin_packing([70_i8,60_i8,50_i8,40_i8,30_i8,20_i8,10_i8],100_i8)
   print '(a,i0,a,i0)', 'minimum bins = ',bp%min_bins,', optimal packings = ',bp%nsol

   c = reshape([0_i8,1_i8,10_i8,10_i8, &
                1_i8,0_i8,10_i8,10_i8, &
                10_i8,10_i8,0_i8,1_i8, &
                10_i8,10_i8,1_i8,0_i8],[4,4])
   tsp = tsp_solver(c)
   print '(a,i0)', 'TSP optimum = ',tsp%tour_length
   print '(a,i0)', 'optimal tours = ',tsp%ntours
end program tsp_binpacking_example
