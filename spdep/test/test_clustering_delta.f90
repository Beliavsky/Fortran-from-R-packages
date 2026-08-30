! SPDX-License-Identifier: GPL-2.0-or-later
! Computational translation derived from R package spdep 1.4-2.
! Translated/modified 2026-08-30 for modern free-form Fortran.
program test_clustering_delta
   use spdep, only : dp, neighbor_list, real_vector, mst_result, spatial_delta_result, &
      cell2nb, nb_adjacency_matrix, nbcosts, mstree, prunecost, skater_groups, &
      metropolis_hastings_weights, linearised_diffusive_weights, &
      iterative_proportional_fitting_weights, graph_distance_weights, &
      spatialdelta, cornish_fisher
   use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
   implicit none

   type(neighbor_list) :: nb
   type(real_vector), allocatable :: costs(:)
   type(mst_result) :: tree
   type(spatial_delta_result) :: delta
   type(spatial_delta_result) :: corrected
   integer, allocatable :: groups(:)
   logical, allocatable :: ladj(:, :)
   real(dp), allocatable :: adjacency(:, :)
   real(dp), allocatable :: w(:, :)
   real(dp), allocatable :: d(:, :)
   real(dp), allocatable :: gains(:)
   real(dp) :: data(5, 1)
   real(dp) :: rw(5)
   integer :: i
   integer :: j
   integer :: failures

   failures = 0
   nb = cell2nb(1, 5)
   data(:, 1) = [0.0_dp, 0.1_dp, 0.2_dp, 10.0_dp, 10.1_dp]
   costs = nbcosts(nb, data)
   tree = mstree(nb, costs)
   call assert_close(tree%total_cost, 10.1_dp, 1.0e-12_dp, "MST total cost", failures)
   gains = prunecost(tree, data)
   call assert_true(maxloc(gains, dim = 1) == 3, "SKATER pruning gain", failures)
   groups = skater_groups(tree, data, 1)
   call assert_true(groups(1) == groups(2) .and. groups(2) == groups(3), &
      "SKATER first cluster", failures)
   call assert_true(groups(4) == groups(5) .and. groups(3) /= groups(4), &
      "SKATER second cluster", failures)

   ladj = nb_adjacency_matrix(nb)
   allocate(adjacency(5, 5))
   adjacency = merge(1.0_dp, 0.0_dp, ladj)
   rw = 0.2_dp

   w = metropolis_hastings_weights(adjacency, rw)
   call assert_rowsum_one(w, "Metropolis-Hastings row sums", failures)
   call assert_stationary(w, rw, "Metropolis-Hastings stationary masses", failures)

   w = linearised_diffusive_weights(adjacency, rw)
   call assert_rowsum_one(w, "diffusive row sums", failures)
   call assert_stationary(w, rw, "diffusive stationary masses", failures)

   w = iterative_proportional_fitting_weights(adjacency, rw)
   call assert_rowsum_one(w, "IPFP row sums", failures)
   call assert_stationary(w, rw, "IPFP stationary masses", failures)

   w = graph_distance_weights(adjacency, rw)
   call assert_rowsum_one(w, "graph-distance row sums", failures)
   call assert_stationary(w, rw, "graph-distance stationary masses", failures)

   allocate(d(5, 5))
   do i = 1, 5
      do j = 1, 5
         d(i, j) = (data(i, 1) - data(j, 1)) ** 2
      end do
   end do
   delta = spatialdelta(d, w, rw, "two.sided")
   call assert_true(.not. ieee_is_nan(delta%delta), "spatial delta finite", failures)
   call assert_true(delta%p_value >= 0.0_dp .and. delta%p_value <= 1.0_dp, &
      "spatial delta p-value", failures)
   corrected = cornish_fisher(delta, "two.sided")
   call assert_true(.not. ieee_is_nan(corrected%z_score), "Cornish-Fisher finite", failures)

   if (failures /= 0) error stop "test_clustering_delta failed"
   print '(a)', "test_clustering_delta: PASS"

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

   subroutine assert_rowsum_one(w, label, failures)
      real(dp), intent(in) :: w(:, :) !! Square transition matrix whose row sums should equal one.
      character(len=*), intent(in) :: label !! Human-readable test label printed when the assertion fails.
      integer, intent(inout) :: failures !! Running number of failed assertions, incremented on failure.

      call assert_true(maxval(abs(sum(w, dim = 2) - 1.0_dp)) < 1.0e-8_dp, label, failures)
   end subroutine assert_rowsum_one

   subroutine assert_stationary(w, rw, label, failures)
      real(dp), intent(in) :: w(:, :) !! Square transition matrix tested against the supplied stationary masses.
      real(dp), intent(in) :: rw(:) !! Regional stationary masses expected to satisfy transpose(W)*rw = rw.
      character(len=*), intent(in) :: label !! Human-readable test label printed when the assertion fails.
      integer, intent(inout) :: failures !! Running number of failed assertions, incremented on failure.
      real(dp), allocatable :: fitted(:)

      fitted = matmul(transpose(w), rw)
      call assert_true(maxval(abs(fitted - rw)) < 1.0e-8_dp, label, failures)
   end subroutine assert_stationary

end program test_clustering_delta
