! SPDX-License-Identifier: GPL-3.0-or-later
module loo_folds
   use loo_types, only: dp
   use r_ordering, only: r_order
   implicit none
   private
   public :: kfold_split_random, kfold_split_stratified, kfold_split_grouped
contains
   pure function valid_permutation(permutation, n) result(valid)
      integer, intent(in) :: permutation(:) !! Candidate one-based permutation of n distinct indexes.
      integer, intent(in) :: n !! Required permutation length and upper index limit.
      logical :: valid
      integer :: i
      valid = .false.
      if (size(permutation) /= n .or. n < 1) return
      if (any(permutation < 1) .or. any(permutation > n)) return
      do i = 1, n
         if (count(permutation == i) /= 1) return
      end do
      valid = .true.
   end function valid_permutation

   pure function kfold_split_random(k, permutation) result(folds)
      integer, intent(in) :: k !! Number of folds, at least two and no greater than the observation count.
      integer, intent(in) :: permutation(:) !! Caller-supplied random permutation of indexes 1 through N; replaces R's RNG draw.
      integer :: folds(size(permutation)), n, i, boundary, j
      n = size(permutation)
      folds = 0
      if (k < 2 .or. k > n) return
      if (.not. valid_permutation(permutation, n)) return
      do i = 1, n
         folds(i) = 1
         do j = 1, k - 1
            boundary = ceiling(1.0_dp + real(n - 1, dp) * real(j, dp) / real(k, dp))
            if (permutation(i) >= boundary) folds(i) = j + 1
         end do
      end do
   end function kfold_split_random

   pure function kfold_split_stratified(k, categories, permutation) result(folds)
      integer, intent(in) :: k !! Number of balanced folds, between two and the number of observations.
      integer, intent(in) :: categories(:) !! Integer categorical labels, ordered numerically as R factor levels.
      integer, intent(in) :: permutation(:) !! Random order of observation indexes; within-category relative order is preserved.
      integer :: folds(size(categories)), n, i, j, at, category
      integer, allocatable :: order(:)
      n = size(categories)
      folds = 0
      if (k < 2 .or. k > n) return
      if (.not. valid_permutation(permutation, n)) return
      call r_order(real(categories, dp), order)
      at = 0
      do i = 1, n
         if (i > 1) then
            if (categories(order(i)) == categories(order(i - 1))) cycle
         end if
         category = categories(order(i))
         do j = 1, n
            if (categories(permutation(j)) /= category) cycle
            folds(permutation(j)) = mod(at, k) + 1
            at = at + 1
         end do
      end do
   end function kfold_split_stratified

   pure function kfold_split_grouped(k, groups, permutation) result(folds)
      integer, intent(in) :: k !! Number of folds, at least two and at most the number of distinct groups.
      integer, intent(in) :: groups(:) !! Integer group labels; every observation in a group is assigned together.
      integer, intent(in) :: permutation(:) !! Permutation of sorted unique group indexes; ignored only when groups equal k.
      integer :: folds(size(groups)), levels(size(groups)), labels(size(groups)), count_levels
      integer :: i, j, big, n_big, boundary, bin
      integer, allocatable :: order(:)
      folds = 0
      if (size(groups) < k .or. k < 2) return
      call r_order(real(groups, dp), order)
      count_levels = 0
      do i = 1, size(groups)
         if (i == 1) then
            count_levels = 1
            levels(count_levels) = groups(order(i))
         else if (groups(order(i)) /= groups(order(i - 1))) then
            count_levels = count_levels + 1
            levels(count_levels) = groups(order(i))
         end if
         labels(order(i)) = count_levels
      end do
      if (count_levels < k) return
      if (count_levels == k) then
         folds = labels
         return
      end if
      if (.not. valid_permutation(permutation, count_levels)) return
      big = ceiling(real(count_levels, dp) / real(k, dp))
      n_big = k - (big * k - count_levels)
      do i = 1, size(groups)
         j = permutation(labels(i))
         boundary = 0
         bin = 1
         do while (bin < k)
            if (bin <= n_big) then
               boundary = boundary + big
            else
               boundary = boundary + big - 1
            end if
            if (j <= boundary) exit
            bin = bin + 1
         end do
         folds(i) = bin
      end do
   end function kfold_split_grouped
end module loo_folds
