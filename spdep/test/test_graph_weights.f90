! SPDX-License-Identifier: GPL-2.0-or-later
! Computational translation derived from R package spdep 1.4-2.
! Translated/modified 2026-08-30 for modern free-form Fortran.
program test_graph_weights
   use spdep, only : dp, neighbor_list, spatial_weights, knn_result, weights_constants, &
      cell2nb, card, nb2listw, lag_listw, spweights_constants, knearneigh, knn2nb, &
      dnearneigh, connected_components, is_symmetric_nb, include_self, remove_self, &
      graph_distance_matrix, great_circle_distance
   implicit none

   type(neighbor_list) :: nb
   type(neighbor_list) :: nb2
   type(spatial_weights) :: listw
   type(knn_result) :: knn
   type(weights_constants) :: c
   integer, allocatable :: counts(:)
   integer, allocatable :: components(:)
   real(dp), allocatable :: lagged(:)
   real(dp), allocatable :: distances(:, :)
   real(dp) :: coords(4, 2)
   real(dp) :: x(4)
   integer :: failures

   failures = 0
   nb = cell2nb(1, 4)
   counts = card(nb)
   call assert_true(all(counts == [1, 2, 2, 1]), "cell2nb/card", failures)
   call assert_true(is_symmetric_nb(nb), "rook chain symmetry", failures)

   listw = nb2listw(nb, "W")
   x = [1.0_dp, 2.0_dp, 4.0_dp, 8.0_dp]
   lagged = lag_listw(listw, x)
   call assert_close(lagged(1), 2.0_dp, 1.0e-12_dp, "lag first", failures)
   call assert_close(lagged(2), 2.5_dp, 1.0e-12_dp, "lag second", failures)
   call assert_close(lagged(3), 5.0_dp, 1.0e-12_dp, "lag third", failures)
   call assert_close(lagged(4), 4.0_dp, 1.0e-12_dp, "lag fourth", failures)

   c = spweights_constants(listw)
   call assert_close(c%s0, 4.0_dp, 1.0e-12_dp, "S0", failures)
   call assert_close(c%s1, 5.5_dp, 1.0e-12_dp, "S1", failures)
   call assert_close(c%s2, 17.0_dp, 1.0e-12_dp, "S2", failures)

   nb2 = include_self(nb)
   call assert_true(nb2%self_included, "include_self flag", failures)
   nb2 = remove_self(nb2)
   call assert_true(.not. nb2%self_included, "remove_self flag", failures)

   coords(:, 1) = [0.0_dp, 1.0_dp, 2.0_dp, 4.0_dp]
   coords(:, 2) = 0.0_dp
   knn = knearneigh(coords, 1)
   call assert_true(all(knn%index(:, 1) == [2, 1, 2, 3]), "knearest indices", failures)
   nb2 = knn2nb(knn, .true.)
   call assert_true(is_symmetric_nb(nb2), "symmetric knn graph", failures)

   nb2 = dnearneigh(coords, 0.0_dp, 1.01_dp)
   components = connected_components(nb2)
   call assert_true(maxval(components) == 2, "distance-band components", failures)
   distances = graph_distance_matrix(nb)
   call assert_close(distances(1, 4), 3.0_dp, 1.0e-12_dp, "graph distance", failures)

   call assert_close(great_circle_distance(0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp), &
      0.0_dp, 1.0e-12_dp, "great-circle identity", failures)

   if (failures /= 0) error stop "test_graph_weights failed"
   print '(a)', "test_graph_weights: PASS"

contains

   subroutine assert_true(condition, label, failures)
      logical, intent(in) :: condition !! Boolean condition that must be true for the test to pass.
      character(len=*), intent(in) :: label !! Human-readable test label printed when the assertion fails.
      integer, intent(inout) :: failures !! Running number of failed assertions, incremented on failure.

      if (.not. condition) then
         failures = failures + 1
         print '(a,1x,a)', "FAIL:", trim(label)
      end if
   end subroutine assert_true

   subroutine assert_close(actual, expected, tolerance, label, failures)
      real(dp), intent(in) :: actual !! Computed scalar value being checked.
      real(dp), intent(in) :: expected !! Reference scalar value expected from the deterministic case.
      real(dp), intent(in) :: tolerance !! Maximum permitted absolute difference between actual and expected.
      character(len=*), intent(in) :: label !! Human-readable test label printed when the assertion fails.
      integer, intent(inout) :: failures !! Running number of failed assertions, incremented on failure.

      if (abs(actual - expected) > tolerance) then
         failures = failures + 1
         print '(a,1x,a,2(1x,es24.16))', "FAIL:", trim(label), actual, expected
      end if
   end subroutine assert_close

end program test_graph_weights
