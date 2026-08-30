! SPDX-License-Identifier: GPL-2.0-or-later
! Computational translation derived from R package spdep 1.4-2.
! Translated/modified 2026-08-30 for modern free-form Fortran.
program moran_example
   use spdep, only : dp, neighbor_list, spatial_weights, spatial_test_result, &
      cell2nb, nb2listw, moran_test
   implicit none

   type(neighbor_list) :: nb
   type(spatial_weights) :: listw
   type(spatial_test_result) :: res
   real(dp) :: x(6)

   nb = cell2nb(1, 6)
   listw = nb2listw(nb, "W")
   x = [1.0_dp, 2.0_dp, 3.0_dp, 5.0_dp, 8.0_dp, 13.0_dp]
   res = moran_test(x, listw)

   print '(a,f12.6)', "Moran I:     ", res%statistic
   print '(a,f12.6)', "Expectation: ", res%expectation
   print '(a,f12.6)', "Z score:     ", res%z_score
   print '(a,f12.6)', "Two-sided p: ", res%p_value
end program moran_example
