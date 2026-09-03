! SPDX-License-Identifier: GPL-2.0-only OR GPL-3.0-only
!
! Neighbor-joining follows ape src/nj.c (Copyright 2006-2023 Emmanuel Paradis).
! BIONJ follows ape src/BIONJ.c (Copyright 2007-2008 Olivier Gascuel,
! Hoa Sien Cuong; R port by Vincent Lefort and Emmanuel Paradis).
! MVR follows ape src/mvr.c (Copyright 2011-2012 Andrei-Alin Popescu).
! NJ* follows ape src/njs.c (Copyright 2011-2013 Andrei-Alin Popescu).
! BIONJ* follows ape src/bionjs.c (Copyright 2011-2014 Andrei-Alin Popescu).
! MVR* follows ape src/mvrs.c (Copyright 2011-2012 Andrei-Alin Popescu).
! Missing-distance completion follows ape src/additive.c and src/ultrametric.c
! (Copyright 2011 Andrei-Alin Popescu).
module ape_reconstruction
   use ieee_arithmetic, only : ieee_is_nan
   use r_kinds, only : dp
   use ape_types, only : phylo_tree, make_phylo_tree
   implicit none
   private

   public :: nj
   public :: bionj
   public :: mvr
   public :: njs
   public :: bionjs
   public :: mvrs
   public :: additive_completion
   public :: ultrametric_completion

contains

   pure subroutine nj(distance, tree, info)
      !! Reconstructs an unrooted tree by the neighbor-joining algorithm.
      real(dp), intent(in) :: distance(:, :) !! Symmetric taxon-distance matrix; diagonal entries are ignored.
      type(phylo_tree), intent(out) :: tree !! Reconstructed ape-style tree with a trifurcating root representation.
      integer, intent(out) :: info !! Status code: zero on success, nonzero for invalid input.
      real(dp), allocatable :: d(:, :)
      real(dp), allocatable :: sums(:)
      real(dp), allocatable :: lengths(:)
      integer, allocatable :: labels(:)
      integer, allocatable :: edges(:, :)
      logical, allocatable :: active(:)
      integer :: a
      integer :: b
      integer :: c
      integer :: cur_node
      integer :: edge_pos
      integer :: i
      integer :: j
      integer :: n
      integer :: r
      real(dp) :: best_q
      real(dp) :: dab
      real(dp) :: q
      real(dp) :: delta

      info = 0
      n = size(distance, 1)
      if (n < 3 .or. size(distance, 2) /= n) then
         info = 1
         return
      end if
      if (.not. valid_distance_matrix(distance)) then
         info = 2
         return
      end if

      allocate(d(n, n), sums(n), labels(n), active(n))
      allocate(edges(2 * n - 3, 2), lengths(2 * n - 3))
      d = distance
      labels = [(i, i = 1, n)]
      active = .true.
      cur_node = 2 * n - 2
      edge_pos = 0
      r = n

      do while (r > 3)
         sums = 0.0_dp
         do i = 1, n
            if (.not. active(i)) cycle
            do j = 1, n
               if (.not. active(j) .or. j == i) cycle
               sums(i) = sums(i) + d(i, j)
            end do
         end do

         best_q = huge(1.0_dp)
         a = 0
         b = 0
         do i = 1, n - 1
            if (.not. active(i)) cycle
            do j = i + 1, n
               if (.not. active(j)) cycle
               q = real(r - 2, dp) * d(i, j) - sums(i) - sums(j)
               if (q < best_q) then
                  best_q = q
                  a = i
                  b = j
               end if
            end do
         end do
         if (a == 0 .or. b == 0) then
            info = 3
            return
         end if

         dab = d(a, b)
         delta = (sums(a) - sums(b)) / real(r - 2, dp)
         edge_pos = edge_pos + 1
         edges(edge_pos, :) = [cur_node, labels(a)]
         lengths(edge_pos) = 0.5_dp * (dab + delta)
         edge_pos = edge_pos + 1
         edges(edge_pos, :) = [cur_node, labels(b)]
         lengths(edge_pos) = 0.5_dp * (dab - delta)

         do i = 1, n
            if (.not. active(i) .or. i == a .or. i == b) cycle
            d(a, i) = 0.5_dp * (d(a, i) + d(b, i) - dab)
            d(i, a) = d(a, i)
         end do
         labels(a) = cur_node
         active(b) = .false.
         cur_node = cur_node - 1
         r = r - 1
      end do

      call final_three(active, a, b, c, info)
      if (info /= 0) return
      edge_pos = edge_pos + 1
      edges(edge_pos, :) = [cur_node, labels(a)]
      lengths(edge_pos) = 0.5_dp * (d(a, b) + d(a, c) - d(b, c))
      edge_pos = edge_pos + 1
      edges(edge_pos, :) = [cur_node, labels(b)]
      lengths(edge_pos) = 0.5_dp * (d(a, b) + d(b, c) - d(a, c))
      edge_pos = edge_pos + 1
      edges(edge_pos, :) = [cur_node, labels(c)]
      lengths(edge_pos) = 0.5_dp * (d(a, c) + d(b, c) - d(a, b))

      tree = make_phylo_tree(n, edges, lengths)
   end subroutine nj

   pure subroutine bionj(distance, tree, info)
      !! Reconstructs an unrooted tree by Gascuel's BIONJ variance-weighted algorithm.
      real(dp), intent(in) :: distance(:, :) !! Symmetric taxon-distance matrix; diagonal entries are ignored.
      type(phylo_tree), intent(out) :: tree !! Reconstructed ape-style tree with a trifurcating root representation.
      integer, intent(out) :: info !! Status code: zero on success, nonzero for invalid input.
      real(dp), allocatable :: d(:, :)
      real(dp), allocatable :: variance(:, :)
      real(dp), allocatable :: sums(:)
      real(dp), allocatable :: lengths(:)
      integer, allocatable :: labels(:)
      integer, allocatable :: edges(:, :)
      logical, allocatable :: active(:)
      integer :: a
      integer :: b
      integer :: c
      integer :: cur_node
      integer :: edge_pos
      integer :: i
      integer :: j
      integer :: n
      integer :: r
      real(dp) :: best_q
      real(dp) :: la
      real(dp) :: lb
      real(dp) :: lambda
      real(dp) :: q
      real(dp) :: vab
      real(dp) :: lambda_sum

      info = 0
      n = size(distance, 1)
      if (n < 3 .or. size(distance, 2) /= n) then
         info = 1
         return
      end if
      if (.not. valid_distance_matrix(distance)) then
         info = 2
         return
      end if

      allocate(d(n, n), variance(n, n), sums(n), labels(n), active(n))
      allocate(edges(2 * n - 3, 2), lengths(2 * n - 3))
      d = distance
      variance = distance
      do i = 1, n
         variance(i, i) = 0.0_dp
      end do
      labels = [(i, i = 1, n)]
      active = .true.
      cur_node = 2 * n - 2
      edge_pos = 0
      r = n

      do while (r > 3)
         sums = 0.0_dp
         do i = 1, n
            if (.not. active(i)) cycle
            do j = 1, n
               if (.not. active(j) .or. j == i) cycle
               sums(i) = sums(i) + d(i, j)
            end do
         end do

         best_q = huge(1.0_dp)
         a = 0
         b = 0
         do i = 1, n
            if (.not. active(i)) cycle
            do j = 1, i - 1
               if (.not. active(j)) cycle
               q = real(r - 2, dp) * d(i, j) - sums(i) - sums(j)
               if (q < best_q - 1.0e-6_dp) then
                  best_q = q
                  a = i
                  b = j
               end if
            end do
         end do
         if (a == 0 .or. b == 0) then
            info = 3
            return
         end if

         vab = variance(a, b)
         la = 0.5_dp * (d(a, b) + (sums(a) - sums(b)) / real(r - 2, dp))
         lb = 0.5_dp * (d(a, b) + (sums(b) - sums(a)) / real(r - 2, dp))
         if (abs(vab) <= tiny(1.0_dp)) then
            lambda = 0.5_dp
         else
            lambda_sum = 0.0_dp
            do i = 1, n
               if (.not. active(i) .or. i == a .or. i == b) cycle
               lambda_sum = lambda_sum + variance(b, i) - variance(a, i)
            end do
            lambda = 0.5_dp + lambda_sum / (2.0_dp * real(r - 2, dp) * vab)
            lambda = max(0.0_dp, min(1.0_dp, lambda))
         end if

         edge_pos = edge_pos + 1
         edges(edge_pos, :) = [cur_node, labels(a)]
         lengths(edge_pos) = la
         edge_pos = edge_pos + 1
         edges(edge_pos, :) = [cur_node, labels(b)]
         lengths(edge_pos) = lb

         do i = 1, n
            if (.not. active(i) .or. i == a .or. i == b) cycle
            d(a, i) = lambda * (d(a, i) - la) + (1.0_dp - lambda) * (d(b, i) - lb)
            d(i, a) = d(a, i)
            variance(a, i) = lambda * variance(a, i) + (1.0_dp - lambda) * variance(b, i) &
               - lambda * (1.0_dp - lambda) * vab
            variance(i, a) = variance(a, i)
         end do
         labels(a) = cur_node
         active(b) = .false.
         cur_node = cur_node - 1
         r = r - 1
      end do

      call final_three(active, a, b, c, info)
      if (info /= 0) return
      edge_pos = edge_pos + 1
      edges(edge_pos, :) = [cur_node, labels(a)]
      lengths(edge_pos) = 0.5_dp * (d(a, b) + d(a, c) - d(b, c))
      edge_pos = edge_pos + 1
      edges(edge_pos, :) = [cur_node, labels(b)]
      lengths(edge_pos) = 0.5_dp * (d(a, b) + d(b, c) - d(a, c))
      edge_pos = edge_pos + 1
      edges(edge_pos, :) = [cur_node, labels(c)]
      lengths(edge_pos) = 0.5_dp * (d(a, c) + d(b, c) - d(a, b))

      tree = make_phylo_tree(n, edges, lengths)
   end subroutine bionj

   pure subroutine mvr(distance, variance_input, tree, info)
      !! Reconstructs a tree with ape's minimum-variance reduction algorithm.
      real(dp), intent(in) :: distance(:, :) !! Symmetric taxon-distance matrix; diagonal entries are ignored.
      real(dp), intent(in) :: variance_input(:, :) !! Symmetric positive pairwise variance matrix matching `distance`.
      type(phylo_tree), intent(out) :: tree !! Reconstructed ape-style tree with a trifurcating root representation.
      integer, intent(out) :: info !! Status code: zero on success, nonzero for invalid dimensions or variances.
      real(dp), allocatable :: d(:, :)
      real(dp), allocatable :: variance(:, :)
      real(dp), allocatable :: sums(:)
      real(dp), allocatable :: lengths(:)
      integer, allocatable :: labels(:)
      integer, allocatable :: edges(:, :)
      logical, allocatable :: active(:)
      integer :: a
      integer :: b
      integer :: c
      integer :: cur_node
      integer :: edge_pos
      integer :: i
      integer :: j
      integer :: n
      integer :: r
      real(dp) :: best_q
      real(dp) :: branch_a
      real(dp) :: branch_b
      real(dp) :: denominator
      real(dp) :: inv_sum
      real(dp) :: mu
      real(dp) :: q
      real(dp) :: weight

      info = 0
      n = size(distance, 1)
      if (n < 3 .or. size(distance, 2) /= n) then
         info = 1
         return
      end if
      if (size(variance_input, 1) /= n .or. size(variance_input, 2) /= n) then
         info = 1
         return
      end if
      if (.not. valid_distance_matrix(distance)) then
         info = 2
         return
      end if
      do i = 1, n - 1
         do j = i + 1, n
            if (variance_input(i, j) <= 0.0_dp .or. variance_input(j, i) <= 0.0_dp) then
               info = 3
               return
            end if
            if (abs(variance_input(i, j) - variance_input(j, i)) > 100.0_dp * epsilon(1.0_dp) * &
               max(1.0_dp, abs(variance_input(i, j)), abs(variance_input(j, i)))) then
               info = 3
               return
            end if
         end do
      end do

      allocate(d(n, n), variance(n, n), sums(n), labels(n), active(n))
      allocate(edges(2 * n - 3, 2), lengths(2 * n - 3))
      d = distance
      variance = variance_input
      labels = [(i, i = 1, n)]
      active = .true.
      cur_node = 2 * n - 2
      edge_pos = 0
      r = n

      do while (r > 3)
         sums = 0.0_dp
         do i = 1, n
            if (.not. active(i)) cycle
            do j = 1, n
               if (.not. active(j) .or. i == j) cycle
               sums(i) = sums(i) + d(i, j)
            end do
         end do

         best_q = huge(1.0_dp)
         a = 0
         b = 0
         do i = 1, n - 1
            if (.not. active(i)) cycle
            do j = i + 1, n
               if (.not. active(j)) cycle
               q = real(r - 2, dp) * d(i, j) - sums(i) - sums(j)
               if (q < best_q) then
                  best_q = q
                  a = i
                  b = j
               end if
            end do
         end do
         if (a == 0 .or. b == 0) then
            info = 4
            return
         end if

         inv_sum = 0.0_dp
         do i = 1, n
            if (.not. active(i) .or. i == a .or. i == b) cycle
            denominator = variance(i, a) + variance(i, b)
            if (denominator <= 0.0_dp) then
               info = 5
               return
            end if
            inv_sum = inv_sum + 1.0_dp / denominator
         end do
         if (inv_sum <= 0.0_dp) then
            info = 5
            return
         end if
         mu = 0.5_dp / inv_sum
         branch_a = 0.5_dp * d(a, b)
         do i = 1, n
            if (.not. active(i) .or. i == a .or. i == b) cycle
            denominator = variance(i, a) + variance(i, b)
            weight = mu / denominator
            branch_a = branch_a + weight * (d(i, a) - d(i, b))
         end do
         branch_b = d(a, b) - branch_a

         edge_pos = edge_pos + 1
         edges(edge_pos, :) = [cur_node, labels(a)]
         lengths(edge_pos) = branch_a
         edge_pos = edge_pos + 1
         edges(edge_pos, :) = [cur_node, labels(b)]
         lengths(edge_pos) = branch_b

         do i = 1, n
            if (.not. active(i) .or. i == a .or. i == b) cycle
            denominator = variance(i, a) + variance(i, b)
            if (denominator <= 0.0_dp) then
               info = 5
               return
            end if
            weight = variance(i, b) / denominator
            d(a, i) = weight * (d(i, a) - branch_a) + (1.0_dp - weight) * (d(i, b) - branch_b)
            d(i, a) = d(a, i)
            variance(a, i) = variance(i, a) * variance(i, b) / denominator
            variance(i, a) = variance(a, i)
         end do
         labels(a) = cur_node
         active(b) = .false.
         cur_node = cur_node - 1
         r = r - 1
      end do

      call final_three(active, a, b, c, info)
      if (info /= 0) return
      edge_pos = edge_pos + 1
      edges(edge_pos, :) = [cur_node, labels(a)]
      lengths(edge_pos) = 0.5_dp * (d(a, b) + d(a, c) - d(b, c))
      edge_pos = edge_pos + 1
      edges(edge_pos, :) = [cur_node, labels(b)]
      lengths(edge_pos) = 0.5_dp * (d(a, b) + d(b, c) - d(a, c))
      edge_pos = edge_pos + 1
      edges(edge_pos, :) = [cur_node, labels(c)]
      lengths(edge_pos) = 0.5_dp * (d(a, c) + d(b, c) - d(a, b))

      tree = make_phylo_tree(n, edges, lengths)
   end subroutine mvr

   pure subroutine njs(distance, tree, info, fs)
      !! Reconstructs an unrooted tree from incomplete distances with ape's NJ* algorithm.
      real(dp), intent(in) :: distance(:, :) !! Symmetric distance matrix; NaN or negative off-diagonal entries are missing.
      type(phylo_tree), intent(out) :: tree !! Reconstructed ape-style tree with a trifurcating root representation.
      integer, intent(out) :: info !! Status code: zero on success; nonzero indicates invalid or insufficient distance data.
      integer, intent(in), optional :: fs !! Number of top agglomeration candidates retained before NJ* tie breaking; default is 15.
      real(dp), allocatable :: d(:, :)
      real(dp), allocatable :: dnew(:, :)
      real(dp), allocatable :: lengths(:)
      real(dp), allocatable :: pair_r(:, :)
      integer, allocatable :: pair_s(:, :)
      integer, allocatable :: labels(:)
      integer, allocatable :: labels_new(:)
      integer, allocatable :: edges(:, :)
      integer, allocatable :: survivors(:)
      integer :: a
      integer :: b
      integer :: cur_node
      integer :: edge_pos
      integer :: fs_value
      integer :: i
      integer :: j
      integer :: n
      integer :: r
      logical :: incomplete
      real(dp) :: branch_a
      real(dp) :: branch_b
      real(dp) :: d12
      real(dp) :: d13
      real(dp) :: d23
      real(dp) :: dxy
      real(dp) :: difference_sum
      integer :: common_count

      info = 0
      n = size(distance, 1)
      if (n < 3 .or. size(distance, 2) /= n) then
         info = 1
         return
      end if
      fs_value = 15
      if (present(fs)) fs_value = fs
      if (fs_value < 1) then
         info = 2
         return
      end if
      if (.not. valid_incomplete_distance_matrix(distance)) then
         info = 3
         return
      end if

      allocate(d(n, n), labels(n))
      allocate(edges(2 * n - 3, 2), lengths(2 * n - 3))
      d = distance
      do i = 1, n
         d(i, i) = 0.0_dp
         do j = i + 1, n
            if (is_missing(d(i, j))) then
               d(i, j) = -1.0_dp
               d(j, i) = -1.0_dp
            end if
         end do
      end do
      labels = [(i, i = 1, n)]
      cur_node = 2 * n - 2
      edge_pos = 0
      r = n
      incomplete = .true.

      do while (r > 3)
         if (incomplete) then
            allocate(pair_r(r, r), pair_s(r, r))
            call incomplete_pair_statistics(d, pair_r, pair_s)
            call choose_njs_pair(d, pair_r, pair_s, fs_value, a, b, incomplete, info)
         else
            call choose_complete_nj_pair(d, a, b, info)
         end if
         if (info /= 0) return

         if (incomplete) then
            common_count = pair_s(a, b) - 2
         else
            common_count = r - 2
         end if
         if (common_count <= 0) then
            info = 4
            return
         end if

         difference_sum = 0.0_dp
         do i = 1, r
            if (i == a .or. i == b) cycle
            if (is_missing(d(a, i)) .or. is_missing(d(b, i))) cycle
            difference_sum = difference_sum + d(a, i) - d(b, i)
         end do
         difference_sum = difference_sum / (2.0_dp * real(common_count, dp))
         dxy = 0.5_dp * d(a, b)
         branch_a = dxy + difference_sum
         branch_b = dxy - difference_sum

         edge_pos = edge_pos + 1
         edges(edge_pos, :) = [cur_node, labels(a)]
         lengths(edge_pos) = branch_a
         edge_pos = edge_pos + 1
         edges(edge_pos, :) = [cur_node, labels(b)]
         lengths(edge_pos) = branch_b

         allocate(survivors(r - 2), labels_new(r - 1), dnew(r - 1, r - 1))
         call survivor_indices(r, a, b, survivors)
         dnew = 0.0_dp
         labels_new(1) = cur_node
         do i = 1, r - 2
            labels_new(i + 1) = labels(survivors(i))
            if (.not. is_missing(d(a, survivors(i))) .and. .not. is_missing(d(b, survivors(i)))) then
               dnew(1, i + 1) = 0.5_dp * (d(a, survivors(i)) - branch_a + &
                  d(b, survivors(i)) - branch_b)
            else if (.not. is_missing(d(a, survivors(i)))) then
               dnew(1, i + 1) = d(a, survivors(i)) - branch_a
            else if (.not. is_missing(d(b, survivors(i)))) then
               dnew(1, i + 1) = d(b, survivors(i)) - branch_b
            else
               dnew(1, i + 1) = -1.0_dp
            end if
            dnew(i + 1, 1) = dnew(1, i + 1)
         end do
         do i = 1, r - 3
            do j = i + 1, r - 2
               dnew(i + 1, j + 1) = d(survivors(i), survivors(j))
               dnew(j + 1, i + 1) = dnew(i + 1, j + 1)
            end do
         end do

         call move_alloc(dnew, d)
         call move_alloc(labels_new, labels)
         if (allocated(pair_r)) deallocate(pair_r)
         if (allocated(pair_s)) deallocate(pair_s)
         deallocate(survivors)
         r = r - 1
         cur_node = cur_node - 1
      end do

      d12 = d(1, 2)
      d13 = d(1, 3)
      d23 = d(2, 3)
      call complete_final_three(d12, d13, d23)
      edge_pos = edge_pos + 1
      edges(edge_pos, :) = [cur_node, labels(1)]
      lengths(edge_pos) = 0.5_dp * (d12 + d13 - d23)
      edge_pos = edge_pos + 1
      edges(edge_pos, :) = [cur_node, labels(2)]
      lengths(edge_pos) = 0.5_dp * (d12 + d23 - d13)
      edge_pos = edge_pos + 1
      edges(edge_pos, :) = [cur_node, labels(3)]
      lengths(edge_pos) = 0.5_dp * (d13 + d23 - d12)

      tree = make_phylo_tree(n, edges, lengths)
   end subroutine njs

   pure subroutine bionjs(distance, tree, info, fs)
      !! Reconstructs an unrooted tree from incomplete distances with ape's BIONJ* algorithm.
      real(dp), intent(in) :: distance(:, :) !! Symmetric distance matrix; NaN or negative off-diagonal entries are missing.
      type(phylo_tree), intent(out) :: tree !! Reconstructed ape-style tree with a trifurcating root representation.
      integer, intent(out) :: info !! Status code: zero on success; nonzero indicates invalid or insufficient distance data.
      integer, intent(in), optional :: fs !! Number of top agglomeration candidates retained before NJ* tie breaking; default is 15.
      real(dp), allocatable :: d(:, :)
      real(dp), allocatable :: dnew(:, :)
      real(dp), allocatable :: lengths(:)
      real(dp), allocatable :: pair_r(:, :)
      real(dp), allocatable :: variance(:, :)
      real(dp), allocatable :: variance_new(:, :)
      integer, allocatable :: pair_s(:, :)
      integer, allocatable :: labels(:)
      integer, allocatable :: labels_new(:)
      integer, allocatable :: edges(:, :)
      integer, allocatable :: survivors(:)
      integer :: a
      integer :: b
      integer :: common_count
      integer :: cur_node
      integer :: edge_pos
      integer :: fs_value
      integer :: i
      integer :: j
      integer :: n
      integer :: r
      logical :: incomplete
      real(dp) :: branch_a
      real(dp) :: branch_b
      real(dp) :: d12
      real(dp) :: d13
      real(dp) :: d23
      real(dp) :: dxy
      real(dp) :: difference_sum
      real(dp) :: lambda
      real(dp) :: lambda_sum
      real(dp) :: vab

      info = 0
      n = size(distance, 1)
      if (n < 3 .or. size(distance, 2) /= n) then
         info = 1
         return
      end if
      fs_value = 15
      if (present(fs)) fs_value = fs
      if (fs_value < 1) then
         info = 2
         return
      end if
      if (.not. valid_incomplete_distance_matrix(distance)) then
         info = 3
         return
      end if

      allocate(d(n, n), variance(n, n), labels(n))
      allocate(edges(2 * n - 3, 2), lengths(2 * n - 3))
      d = distance
      do i = 1, n
         d(i, i) = 0.0_dp
         do j = i + 1, n
            if (is_missing(d(i, j))) then
               d(i, j) = -1.0_dp
               d(j, i) = -1.0_dp
            end if
         end do
      end do
      variance = d
      labels = [(i, i = 1, n)]
      cur_node = 2 * n - 2
      edge_pos = 0
      r = n
      incomplete = .true.

      do while (r > 3)
         if (incomplete) then
            allocate(pair_r(r, r), pair_s(r, r))
            call incomplete_pair_statistics(d, pair_r, pair_s)
            call choose_njs_pair(d, pair_r, pair_s, fs_value, a, b, incomplete, info)
         else
            call choose_complete_nj_pair(d, a, b, info)
         end if
         if (info /= 0) return

         if (incomplete) then
            common_count = pair_s(a, b) - 2
         else
            common_count = r - 2
         end if
         if (common_count <= 0) then
            info = 4
            return
         end if

         difference_sum = 0.0_dp
         lambda_sum = 0.0_dp
         do i = 1, r
            if (i == a .or. i == b) cycle
            if (is_missing(d(a, i)) .or. is_missing(d(b, i))) cycle
            difference_sum = difference_sum + d(a, i) - d(b, i)
            lambda_sum = lambda_sum + variance(b, i) - variance(a, i)
         end do
         difference_sum = difference_sum / (2.0_dp * real(common_count, dp))
         dxy = 0.5_dp * d(a, b)
         branch_a = dxy + difference_sum
         branch_b = dxy - difference_sum

         vab = variance(a, b)
         if (abs(vab) > 0.0_dp) then
            lambda = 0.5_dp + lambda_sum / (2.0_dp * real(common_count, dp) * vab)
         else
            lambda = 0.5_dp
         end if
         lambda = max(0.0_dp, min(1.0_dp, lambda))

         edge_pos = edge_pos + 1
         edges(edge_pos, :) = [cur_node, labels(a)]
         lengths(edge_pos) = branch_a
         edge_pos = edge_pos + 1
         edges(edge_pos, :) = [cur_node, labels(b)]
         lengths(edge_pos) = branch_b

         allocate(survivors(r - 2), labels_new(r - 1), dnew(r - 1, r - 1), variance_new(r - 1, r - 1))
         call survivor_indices(r, a, b, survivors)
         dnew = 0.0_dp
         variance_new = 0.0_dp
         labels_new(1) = cur_node
         do i = 1, r - 2
            labels_new(i + 1) = labels(survivors(i))
            if (.not. is_missing(d(a, survivors(i))) .and. .not. is_missing(d(b, survivors(i)))) then
               dnew(1, i + 1) = lambda * (d(a, survivors(i)) - branch_a) + &
                  (1.0_dp - lambda) * (d(b, survivors(i)) - branch_b)
               variance_new(1, i + 1) = lambda * variance(a, survivors(i)) + &
                  (1.0_dp - lambda) * variance(b, survivors(i)) - &
                  lambda * (1.0_dp - lambda) * variance(a, b)
            else if (.not. is_missing(d(a, survivors(i)))) then
               dnew(1, i + 1) = d(a, survivors(i)) - branch_a
               variance_new(1, i + 1) = variance(a, survivors(i))
            else if (.not. is_missing(d(b, survivors(i)))) then
               dnew(1, i + 1) = d(b, survivors(i)) - branch_b
               variance_new(1, i + 1) = variance(b, survivors(i))
            else
               dnew(1, i + 1) = -1.0_dp
               variance_new(1, i + 1) = -1.0_dp
            end if
            dnew(i + 1, 1) = dnew(1, i + 1)
            variance_new(i + 1, 1) = variance_new(1, i + 1)
         end do
         do i = 1, r - 3
            do j = i + 1, r - 2
               dnew(i + 1, j + 1) = d(survivors(i), survivors(j))
               dnew(j + 1, i + 1) = dnew(i + 1, j + 1)
               variance_new(i + 1, j + 1) = variance(survivors(i), survivors(j))
               variance_new(j + 1, i + 1) = variance_new(i + 1, j + 1)
            end do
         end do

         call move_alloc(dnew, d)
         call move_alloc(variance_new, variance)
         call move_alloc(labels_new, labels)
         if (allocated(pair_r)) deallocate(pair_r)
         if (allocated(pair_s)) deallocate(pair_s)
         deallocate(survivors)
         r = r - 1
         cur_node = cur_node - 1
      end do

      d12 = d(1, 2)
      d13 = d(1, 3)
      d23 = d(2, 3)
      call complete_final_three(d12, d13, d23)
      edge_pos = edge_pos + 1
      edges(edge_pos, :) = [cur_node, labels(1)]
      lengths(edge_pos) = 0.5_dp * (d12 + d13 - d23)
      edge_pos = edge_pos + 1
      edges(edge_pos, :) = [cur_node, labels(2)]
      lengths(edge_pos) = 0.5_dp * (d12 + d23 - d13)
      edge_pos = edge_pos + 1
      edges(edge_pos, :) = [cur_node, labels(3)]
      lengths(edge_pos) = 0.5_dp * (d13 + d23 - d12)

      tree = make_phylo_tree(n, edges, lengths)
   end subroutine bionjs

   pure subroutine mvrs(distance, variance_input, tree, info, fs)
      !! Reconstructs an unrooted tree from incomplete distances with ape's MVR* algorithm.
      real(dp), intent(in) :: distance(:, :) !! Symmetric distance matrix; NaN or negative off-diagonal entries are missing.
      real(dp), intent(in) :: variance_input(:, :) !! Symmetric positive pairwise variance matrix matching `distance`.
      type(phylo_tree), intent(out) :: tree !! Reconstructed ape-style tree with a trifurcating root representation.
      integer, intent(out) :: info !! Status code: zero on success; nonzero indicates invalid or insufficient input data.
      integer, intent(in), optional :: fs !! Number of top agglomeration candidates retained before NJ* tie breaking; default is 15.
      real(dp), allocatable :: d(:, :)
      real(dp), allocatable :: dnew(:, :)
      real(dp), allocatable :: lengths(:)
      real(dp), allocatable :: pair_r(:, :)
      real(dp), allocatable :: variance(:, :)
      real(dp), allocatable :: variance_new(:, :)
      integer, allocatable :: pair_s(:, :)
      integer, allocatable :: labels(:)
      integer, allocatable :: labels_new(:)
      integer, allocatable :: edges(:, :)
      integer, allocatable :: survivors(:)
      integer :: a
      integer :: b
      integer :: cur_node
      integer :: edge_pos
      integer :: fs_value
      integer :: i
      integer :: j
      integer :: n
      integer :: r
      logical :: incomplete
      real(dp) :: branch_a
      real(dp) :: branch_b
      real(dp) :: d12
      real(dp) :: d13
      real(dp) :: d23
      real(dp) :: denominator
      real(dp) :: inverse_sum
      real(dp) :: mu
      real(dp) :: weight

      info = 0
      n = size(distance, 1)
      if (n < 3 .or. size(distance, 2) /= n) then
         info = 1
         return
      end if
      if (size(variance_input, 1) /= n .or. size(variance_input, 2) /= n) then
         info = 1
         return
      end if
      fs_value = 15
      if (present(fs)) fs_value = fs
      if (fs_value < 1) then
         info = 2
         return
      end if
      if (.not. valid_incomplete_distance_matrix(distance)) then
         info = 3
         return
      end if
      if (.not. valid_variance_matrix(variance_input)) then
         info = 4
         return
      end if

      allocate(d(n, n), variance(n, n), labels(n))
      allocate(edges(2 * n - 3, 2), lengths(2 * n - 3))
      d = distance
      do i = 1, n
         d(i, i) = 0.0_dp
         do j = i + 1, n
            if (is_missing(d(i, j))) then
               d(i, j) = -1.0_dp
               d(j, i) = -1.0_dp
            end if
         end do
      end do
      variance = variance_input
      labels = [(i, i = 1, n)]
      cur_node = 2 * n - 2
      edge_pos = 0
      r = n
      incomplete = .true.

      do while (r > 3)
         if (incomplete) then
            allocate(pair_r(r, r), pair_s(r, r))
            call incomplete_pair_statistics(d, pair_r, pair_s)
            call choose_njs_pair(d, pair_r, pair_s, fs_value, a, b, incomplete, info)
         else
            call choose_complete_nj_pair(d, a, b, info)
         end if
         if (info /= 0) return

         inverse_sum = 0.0_dp
         do i = 1, r
            if (i == a .or. i == b) cycle
            if (is_missing(d(a, i)) .or. is_missing(d(b, i))) cycle
            denominator = variance(i, a) + variance(i, b)
            if (denominator <= 0.0_dp) then
               info = 5
               return
            end if
            inverse_sum = inverse_sum + 1.0_dp / denominator
         end do
         if (inverse_sum <= 0.0_dp) then
            info = 5
            return
         end if
         mu = 0.5_dp / inverse_sum
         branch_a = 0.5_dp * d(a, b)
         branch_b = 0.5_dp * d(a, b)
         do i = 1, r
            if (i == a .or. i == b) cycle
            if (is_missing(d(a, i)) .or. is_missing(d(b, i))) cycle
            denominator = variance(i, a) + variance(i, b)
            weight = mu / denominator
            branch_a = branch_a + weight * (d(i, a) - d(i, b))
            branch_b = branch_b + weight * (d(i, b) - d(i, a))
         end do

         edge_pos = edge_pos + 1
         edges(edge_pos, :) = [cur_node, labels(a)]
         lengths(edge_pos) = branch_a
         edge_pos = edge_pos + 1
         edges(edge_pos, :) = [cur_node, labels(b)]
         lengths(edge_pos) = branch_b

         allocate(survivors(r - 2), labels_new(r - 1), dnew(r - 1, r - 1), variance_new(r - 1, r - 1))
         call survivor_indices(r, a, b, survivors)
         dnew = 0.0_dp
         variance_new = 0.0_dp
         labels_new(1) = cur_node
         do i = 1, r - 2
            labels_new(i + 1) = labels(survivors(i))
            if (.not. is_missing(d(a, survivors(i))) .and. .not. is_missing(d(b, survivors(i)))) then
               denominator = variance(survivors(i), b) + variance(survivors(i), a)
               if (denominator <= 0.0_dp) then
                  info = 5
                  return
               end if
               weight = variance(survivors(i), b) / denominator
               dnew(1, i + 1) = weight * (d(a, survivors(i)) - branch_a) + &
                  (1.0_dp - weight) * (d(b, survivors(i)) - branch_b)
               variance_new(1, i + 1) = variance(survivors(i), b) * &
                  variance(survivors(i), a) / denominator
            else if (.not. is_missing(d(a, survivors(i)))) then
               dnew(1, i + 1) = d(a, survivors(i)) - branch_a
               variance_new(1, i + 1) = variance(a, survivors(i))
            else if (.not. is_missing(d(b, survivors(i)))) then
               dnew(1, i + 1) = d(b, survivors(i)) - branch_b
               variance_new(1, i + 1) = variance(b, survivors(i))
            else
               dnew(1, i + 1) = -1.0_dp
               variance_new(1, i + 1) = -1.0_dp
            end if
            dnew(i + 1, 1) = dnew(1, i + 1)
            variance_new(i + 1, 1) = variance_new(1, i + 1)
         end do
         do i = 1, r - 3
            do j = i + 1, r - 2
               dnew(i + 1, j + 1) = d(survivors(i), survivors(j))
               dnew(j + 1, i + 1) = dnew(i + 1, j + 1)
               variance_new(i + 1, j + 1) = variance(survivors(i), survivors(j))
               variance_new(j + 1, i + 1) = variance_new(i + 1, j + 1)
            end do
         end do

         call move_alloc(dnew, d)
         call move_alloc(variance_new, variance)
         call move_alloc(labels_new, labels)
         if (allocated(pair_r)) deallocate(pair_r)
         if (allocated(pair_s)) deallocate(pair_s)
         deallocate(survivors)
         r = r - 1
         cur_node = cur_node - 1
      end do

      d12 = d(1, 2)
      d13 = d(1, 3)
      d23 = d(2, 3)
      call complete_final_three(d12, d13, d23)
      edge_pos = edge_pos + 1
      edges(edge_pos, :) = [cur_node, labels(1)]
      lengths(edge_pos) = 0.5_dp * (d12 + d13 - d23)
      edge_pos = edge_pos + 1
      edges(edge_pos, :) = [cur_node, labels(2)]
      lengths(edge_pos) = 0.5_dp * (d12 + d23 - d13)
      edge_pos = edge_pos + 1
      edges(edge_pos, :) = [cur_node, labels(3)]
      lengths(edge_pos) = 0.5_dp * (d13 + d23 - d12)

      tree = make_phylo_tree(n, edges, lengths)
   end subroutine mvrs

   pure logical function valid_variance_matrix(variance) result(ok)
      !! Validates a symmetric positive off-diagonal variance matrix for MVR and MVR* reconstruction.
      real(dp), intent(in) :: variance(:, :) !! Pairwise variances; off-diagonal entries must be positive and finite.
      integer :: i
      integer :: j
      integer :: n
      real(dp) :: tolerance

      n = size(variance, 1)
      ok = size(variance, 2) == n
      if (.not. ok) return
      tolerance = 64.0_dp * epsilon(1.0_dp)
      do i = 1, n - 1
         do j = i + 1, n
            if (ieee_is_nan(variance(i, j)) .or. ieee_is_nan(variance(j, i))) then
               ok = .false.
               return
            end if
            if (variance(i, j) <= 0.0_dp .or. variance(j, i) <= 0.0_dp) then
               ok = .false.
               return
            end if
            if (abs(variance(i, j) - variance(j, i)) > tolerance * &
               max(1.0_dp, abs(variance(i, j)), abs(variance(j, i)))) then
               ok = .false.
               return
            end if
         end do
      end do
   end function valid_variance_matrix

   pure subroutine choose_complete_nj_pair(distance, a, b, info)
      !! Selects the first pair maximizing ape's complete-matrix NJ agglomeration criterion.
      real(dp), intent(in) :: distance(:, :) !! Current complete symmetric inter-cluster distance matrix.
      integer, intent(out) :: a !! Index of the first selected cluster.
      integer, intent(out) :: b !! Index of the second selected cluster.
      integer, intent(out) :: info !! Zero on success; nonzero if no pair can be selected.
      real(dp), allocatable :: sums(:)
      real(dp) :: best
      real(dp) :: criterion
      integer :: i
      integer :: j
      integer :: n

      n = size(distance, 1)
      allocate(sums(n))
      sums = 0.0_dp
      do i = 1, n
         do j = 1, n
            if (i == j) cycle
            sums(i) = sums(i) + distance(i, j)
         end do
      end do
      best = -huge(1.0_dp)
      a = 0
      b = 0
      do i = 1, n - 1
         do j = i + 1, n
            criterion = sums(i) + sums(j) - real(n - 2, dp) * distance(i, j)
            if (criterion > best) then
               best = criterion
               a = i
               b = j
            end if
         end do
      end do
      if (a == 0 .or. b == 0) then
         info = 1
      else
         info = 0
      end if
   end subroutine choose_complete_nj_pair

   pure subroutine incomplete_pair_statistics(distance, pair_r, pair_s)
      !! Computes NJ* R(x,y) sums and S(x,y) cardinalities for known pairs.
      real(dp), intent(in) :: distance(:, :) !! Current symmetric distance matrix with negative entries marking missing values.
      real(dp), intent(out) :: pair_r(:, :) !! Pairwise NJ* R sums for the current agglomeration state.
      integer, intent(out) :: pair_s(:, :) !! Pairwise NJ* S cardinalities for the current agglomeration state.
      integer :: i
      integer :: j
      integer :: k
      integer :: n

      n = size(distance, 1)
      pair_r = 0.0_dp
      pair_s = 0
      do i = 1, n - 1
         do j = i + 1, n
            if (is_missing(distance(i, j))) cycle
            do k = 1, n
               if (k == i .or. k == j) then
                  pair_s(i, j) = pair_s(i, j) + 1
                  if (i /= k) pair_r(i, j) = pair_r(i, j) + distance(i, k)
                  if (j /= k) pair_r(i, j) = pair_r(i, j) + distance(j, k)
               else if (.not. is_missing(distance(i, k)) .and. .not. is_missing(distance(j, k))) then
                  pair_s(i, j) = pair_s(i, j) + 1
                  pair_r(i, j) = pair_r(i, j) + distance(i, k) + distance(j, k)
               end if
            end do
            pair_r(j, i) = pair_r(i, j)
            pair_s(j, i) = pair_s(i, j)
         end do
      end do
   end subroutine incomplete_pair_statistics

   pure subroutine choose_njs_pair(distance, pair_r, pair_s, fs, a, b, incomplete, info)
      !! Selects an NJ* agglomeration pair using ape's Q shortlist and missing-data tie breakers.
      real(dp), intent(in) :: distance(:, :) !! Current symmetric distance matrix with negative entries denoting missing values.
      real(dp), intent(in) :: pair_r(:, :) !! NJ* R sums corresponding to `distance`.
      integer, intent(in) :: pair_s(:, :) !! NJ* S cardinalities corresponding to `distance`.
      integer, intent(in) :: fs !! Positive maximum number of Q-ranked candidate pairs retained for tie breaking.
      integer, intent(out) :: a !! Index of the first selected cluster in the current matrix.
      integer, intent(out) :: b !! Index of the second selected cluster in the current matrix.
      logical, intent(out) :: incomplete !! True when at least one current inter-cluster distance is missing.
      integer, intent(out) :: info !! Zero on success; nonzero when insufficient information prevents pair selection.
      real(dp), allocatable :: scores(:)
      integer, allocatable :: first(:)
      integer, allocatable :: second(:)
      real(dp) :: best
      real(dp) :: criterion
      real(dp) :: value
      integer :: candidate_count
      integer :: i
      integer :: j
      integer :: k
      integer :: n
      integer :: position

      n = size(distance, 1)
      allocate(scores(fs), first(fs), second(fs))
      scores = -huge(1.0_dp)
      first = 0
      second = 0
      incomplete = .false.

      do i = 1, n - 1
         do j = i + 1, n
            if (is_missing(distance(i, j))) then
               incomplete = .true.
               cycle
            end if
            if (pair_s(i, j) <= 2) cycle
            criterion = pair_r(i, j) / real(pair_s(i, j) - 2, dp) - distance(i, j)
            position = 1
            do while (position <= fs)
               if (.not. (scores(position) > criterion)) exit
               position = position + 1
            end do
            if (position <= fs) then
               do k = fs, position + 1, -1
                  scores(k) = scores(k - 1)
                  first(k) = first(k - 1)
                  second(k) = second(k - 1)
               end do
               scores(position) = criterion
               first(position) = i
               second(position) = j
            end if
         end do
      end do

      if (first(1) == 0 .or. second(1) == 0) then
         a = 0
         b = 0
         info = 1
         return
      end if
      if (.not. incomplete) then
         a = first(1)
         b = second(1)
         info = 0
         return
      end if

      candidate_count = count(first /= 0)
      call retain_best_njs_candidates(distance, first, second, candidate_count, 1, info)
      if (info /= 0) then
         a = 0
         b = 0
         return
      end if
      if (candidate_count > 1) call retain_best_njs_candidates(distance, first, second, candidate_count, 2, info)
      if (candidate_count > 1) call retain_best_njs_candidates(distance, first, second, candidate_count, 3, info)
      if (candidate_count > 1) then
         best = -huge(1.0_dp)
         position = 0
         do i = 1, candidate_count
            value = njs_cnxy(first(i), second(i), distance)
            if (value > best) then
               best = value
               position = i
            end if
         end do
         if (position == 0) then
            a = 0
            b = 0
            info = 1
            return
         end if
         a = first(position)
         b = second(position)
      else
         a = first(1)
         b = second(1)
      end if
      info = 0
   end subroutine choose_njs_pair

   pure subroutine retain_best_njs_candidates(distance, first, second, candidate_count, stage, info)
      !! Retains candidates tied at the maximum value of one NJ* missing-data tie-break statistic.
      real(dp), intent(in) :: distance(:, :) !! Current NJ* distance matrix used to evaluate candidate diagnostics.
      integer, intent(inout) :: first(:) !! First cluster index for each retained candidate; compacted in place.
      integer, intent(inout) :: second(:) !! Second cluster index for each retained candidate; compacted in place.
      integer, intent(inout) :: candidate_count !! Number of active entries in `first` and `second`, updated after filtering.
      integer, intent(in) :: stage !! Tie-break stage: 1 for nxy, 2 for cxy, or 3 for mxy.
      integer, intent(out) :: info !! Zero on success; nonzero if no candidate can be evaluated.
      real(dp), allocatable :: values(:)
      real(dp) :: best
      integer :: i
      integer :: kept

      allocate(values(candidate_count))
      best = -huge(1.0_dp)
      do i = 1, candidate_count
         select case (stage)
         case (1)
            values(i) = njs_nxy(first(i), second(i), distance)
         case (2)
            values(i) = real(njs_cxy(first(i), second(i), distance), dp)
         case (3)
            values(i) = real(njs_mxy(first(i), second(i), distance), dp)
         case default
            info = 1
            return
         end select
         best = max(best, values(i))
      end do
      kept = 0
      do i = 1, candidate_count
         if (abs(values(i) - best) <= 0.0_dp) then
            kept = kept + 1
            first(kept) = first(i)
            second(kept) = second(i)
         end if
      end do
      candidate_count = kept
      if (candidate_count == 0) then
         info = 1
      else
         info = 0
      end if
   end subroutine retain_best_njs_candidates

   pure real(dp) function njs_nxy(x, y, distance) result(value)
      !! Computes ape NJ*'s normalized four-point consistency score N*(x,y).
      integer, intent(in) :: x !! First candidate cluster index.
      integer, intent(in) :: y !! Second candidate cluster index.
      real(dp), intent(in) :: distance(:, :) !! Current symmetric NJ* distance matrix.
      real(dp) :: n1
      real(dp) :: n2
      real(dp) :: total
      integer :: count_valid
      integer :: i
      integer :: j
      integer :: n

      n = size(distance, 1)
      count_valid = 0
      total = 0.0_dp
      do i = 1, n
         do j = 1, n
            if (i == j) cycle
            if ((i == x .and. j == y) .or. (j == x .and. i == y)) cycle
            n1 = 0.0_dp
            n2 = 0.0_dp
            if (i /= x) n1 = distance(i, x)
            if (j /= y) n2 = distance(j, y)
            if (is_missing(n1) .or. is_missing(n2) .or. is_missing(distance(i, j))) cycle
            count_valid = count_valid + 1
            if (n1 + n2 - distance(x, y) - distance(i, j) >= -1.0e-10_dp) total = total + 1.0_dp
         end do
      end do
      if (count_valid == 0) then
         value = 0.0_dp
      else
         value = total / real(count_valid, dp)
      end if
   end function njs_nxy

   pure integer function njs_cxy(x, y, distance) result(value)
      !! Counts evaluable four-point comparisons for an NJ* candidate pair.
      integer, intent(in) :: x !! First candidate cluster index.
      integer, intent(in) :: y !! Second candidate cluster index.
      real(dp), intent(in) :: distance(:, :) !! Current symmetric NJ* distance matrix.
      real(dp) :: n1
      real(dp) :: n2
      integer :: i
      integer :: j
      integer :: n

      n = size(distance, 1)
      value = 0
      do i = 1, n
         do j = 1, n
            if (i == j) cycle
            if ((i == x .and. j == y) .or. (j == x .and. i == y)) cycle
            n1 = 0.0_dp
            n2 = 0.0_dp
            if (i /= x) n1 = distance(i, x)
            if (j /= y) n2 = distance(j, y)
            if (is_missing(n1) .or. is_missing(n2) .or. is_missing(distance(i, j))) cycle
            value = value + 1
         end do
      end do
   end function njs_cxy

   pure integer function njs_mxy(x, y, distance) result(value)
      !! Counts complementary missing-distance incidences for an NJ* candidate pair.
      integer, intent(in) :: x !! First candidate cluster index.
      integer, intent(in) :: y !! Second candidate cluster index.
      real(dp), intent(in) :: distance(:, :) !! Current symmetric NJ* distance matrix.
      logical, allocatable :: missing_x(:)
      logical, allocatable :: missing_y(:)
      integer :: i
      integer :: n

      n = size(distance, 1)
      allocate(missing_x(n), missing_y(n))
      missing_x = .false.
      missing_y = .false.
      do i = 1, n
         if (i /= x .and. is_missing(distance(x, i))) missing_x(i) = .true.
         if (i /= y .and. is_missing(distance(y, i))) missing_y(i) = .true.
      end do
      value = 0
      do i = 1, n
         if (i /= x .and. missing_x(i) .and. .not. missing_y(i)) value = value + 1
         if (i /= y .and. missing_y(i) .and. .not. missing_x(i)) value = value + 1
      end do
   end function njs_mxy

   pure real(dp) function njs_cnxy(x, y, distance) result(value)
      !! Computes ape NJ*'s final unnormalized four-point tie-break statistic.
      integer, intent(in) :: x !! First candidate cluster index.
      integer, intent(in) :: y !! Second candidate cluster index.
      real(dp), intent(in) :: distance(:, :) !! Current symmetric NJ* distance matrix.
      real(dp) :: n1
      real(dp) :: n2
      integer :: i
      integer :: j
      integer :: n

      n = size(distance, 1)
      value = 0.0_dp
      do i = 1, n
         do j = 1, n
            if (i == j) cycle
            if ((i == x .and. j == y) .or. (j == x .and. i == y)) cycle
            n1 = 0.0_dp
            n2 = 0.0_dp
            if (i /= x) n1 = distance(i, x)
            if (j /= y) n2 = distance(j, y)
            if (is_missing(n1) .or. is_missing(n2) .or. is_missing(distance(i, j))) cycle
            value = value + n1 + n2 - distance(x, y) - distance(i, j)
         end do
      end do
   end function njs_cnxy

   pure subroutine survivor_indices(n, a, b, survivors)
      !! Lists current cluster indices not selected for agglomeration, in ascending order.
      integer, intent(in) :: n !! Current number of clusters.
      integer, intent(in) :: a !! First excluded cluster index.
      integer, intent(in) :: b !! Second excluded cluster index.
      integer, intent(out) :: survivors(:) !! Ascending indices of all clusters other than `a` and `b`.
      integer :: i
      integer :: k

      k = 0
      do i = 1, n
         if (i == a .or. i == b) cycle
         k = k + 1
         survivors(k) = i
      end do
   end subroutine survivor_indices

   pure subroutine complete_final_three(d12, d13, d23)
      !! Applies ape NJ*'s deterministic completion rule to the last three inter-cluster distances.
      real(dp), intent(inout) :: d12 !! Distance between final clusters 1 and 2; negative means missing on entry.
      real(dp), intent(inout) :: d13 !! Distance between final clusters 1 and 3; negative means missing on entry.
      real(dp), intent(inout) :: d23 !! Distance between final clusters 2 and 3; negative means missing on entry.
      real(dp) :: values(3)
      real(dp) :: maximum
      integer :: i
      integer :: known
      integer :: known_index
      integer :: missing_index

      values = [d12, d13, d23]
      known = 0
      known_index = 0
      missing_index = 0
      do i = 1, 3
         if (is_missing(values(i))) then
            missing_index = i
         else
            known = known + 1
            known_index = i
         end if
      end do
      select case (known)
      case (2)
         maximum = -huge(1.0_dp)
         do i = 1, 3
            if (.not. is_missing(values(i))) maximum = max(maximum, values(i))
         end do
         values(missing_index) = maximum
      case (1)
         do i = 1, 3
            if (is_missing(values(i))) values(i) = values(known_index)
         end do
      case (0)
         values = 1.0_dp
      end select
      d12 = values(1)
      d13 = values(2)
      d23 = values(3)
   end subroutine complete_final_three

   pure logical function valid_incomplete_distance_matrix(distance) result(ok)
      !! Validates symmetry while permitting NaN or negative off-diagonal entries as missing distances.
      real(dp), intent(in) :: distance(:, :) !! Candidate square symmetric matrix with optional missing off-diagonal values.
      integer :: i
      integer :: j
      integer :: n
      real(dp) :: tolerance

      n = size(distance, 1)
      ok = size(distance, 2) == n
      if (.not. ok) return
      tolerance = 64.0_dp * epsilon(1.0_dp)
      do i = 1, n - 1
         do j = i + 1, n
            if (is_missing(distance(i, j)) .or. is_missing(distance(j, i))) then
               if (.not. (is_missing(distance(i, j)) .and. is_missing(distance(j, i)))) then
                  ok = .false.
                  return
               end if
            else if (abs(distance(i, j) - distance(j, i)) > tolerance * &
               max(1.0_dp, abs(distance(i, j)), abs(distance(j, i)))) then
               ok = .false.
               return
            end if
         end do
      end do
   end function valid_incomplete_distance_matrix

   pure subroutine additive_completion(distance, completed, unresolved)
      !! Completes missing dissimilarities using ape's four-point additive minimax rule.
      real(dp), intent(in) :: distance(:, :) !! Square distance matrix; NaN or negative off-diagonal values denote missing entries.
      real(dp), allocatable, intent(out) :: completed(:, :) !! Completed symmetric matrix; unresolved entries remain NaN.
      integer, intent(out) :: unresolved !! Number of unresolved pairs in the upper triangle after iteration.
      real(dp) :: best
      real(dp) :: candidate
      real(dp) :: maximum
      integer :: i
      integer :: j
      integer :: k
      integer :: l
      integer :: n
      logical :: changed
      logical :: found

      n = size(distance, 1)
      allocate(completed(n, size(distance, 2)))
      completed = distance
      unresolved = 0
      if (size(distance, 2) /= n .or. n == 0) then
         unresolved = -1
         return
      end if
      do i = 1, n
         completed(i, i) = 0.0_dp
      end do
      maximum = 0.0_dp
      do i = 1, n - 1
         do j = i + 1, n
            if (.not. is_missing(completed(i, j))) maximum = max(maximum, completed(i, j))
         end do
      end do

      do
         changed = .false.
         do i = 1, n - 1
            do j = i + 1, n
               if (.not. is_missing(completed(i, j))) cycle
               best = maximum
               found = .false.
               do k = 1, n
                  if (is_missing(completed(i, k)) .or. is_missing(completed(j, k))) cycle
                  do l = 1, n
                     if (k == l) cycle
                     if (is_missing(completed(k, l))) cycle
                     if (is_missing(completed(i, l)) .or. is_missing(completed(j, l))) cycle
                     candidate = max(completed(i, k) + completed(j, l), &
                        completed(i, l) + completed(j, k)) - completed(k, l)
                     best = min(best, candidate)
                     found = .true.
                  end do
               end do
               if (found) then
                  completed(i, j) = best
                  completed(j, i) = best
                  changed = .true.
               end if
            end do
         end do
         if (.not. changed) exit
      end do
      call normalize_unresolved(completed, unresolved)
   end subroutine additive_completion

   pure subroutine ultrametric_completion(distance, completed, unresolved)
      !! Completes missing dissimilarities using ape's ultrametric minimax rule.
      real(dp), intent(in) :: distance(:, :) !! Square distance matrix; NaN or negative off-diagonal values denote missing entries.
      real(dp), allocatable, intent(out) :: completed(:, :) !! Completed symmetric matrix; unresolved entries remain NaN.
      integer, intent(out) :: unresolved !! Number of unresolved pairs in the upper triangle after iteration.
      real(dp) :: best
      real(dp) :: candidate
      real(dp) :: maximum
      integer :: i
      integer :: j
      integer :: k
      integer :: n
      logical :: changed
      logical :: found

      n = size(distance, 1)
      allocate(completed(n, size(distance, 2)))
      completed = distance
      unresolved = 0
      if (size(distance, 2) /= n .or. n == 0) then
         unresolved = -1
         return
      end if
      do i = 1, n
         completed(i, i) = 0.0_dp
      end do
      maximum = 0.0_dp
      do i = 1, n - 1
         do j = i + 1, n
            if (.not. is_missing(completed(i, j))) maximum = max(maximum, completed(i, j))
         end do
      end do

      do
         changed = .false.
         do i = 1, n - 1
            do j = i + 1, n
               if (.not. is_missing(completed(i, j))) cycle
               best = maximum
               found = .false.
               do k = 1, n
                  if (is_missing(completed(i, k)) .or. is_missing(completed(j, k))) cycle
                  candidate = max(completed(i, k), completed(j, k))
                  best = min(best, candidate)
                  found = .true.
               end do
               if (found) then
                  completed(i, j) = best
                  completed(j, i) = best
                  changed = .true.
               end if
            end do
         end do
         if (.not. changed) exit
      end do
      call normalize_unresolved(completed, unresolved)
   end subroutine ultrametric_completion

   pure logical function valid_distance_matrix(distance) result(ok)
      real(dp), intent(in) :: distance(:, :) !! Candidate symmetric distance matrix to validate.
      integer :: i
      integer :: j
      integer :: n
      real(dp) :: tolerance

      n = size(distance, 1)
      ok = size(distance, 2) == n
      if (.not. ok) return
      tolerance = 64.0_dp * epsilon(1.0_dp)
      do i = 1, n
         if (ieee_is_nan(distance(i, i))) then
            ok = .false.
            return
         end if
         do j = i + 1, n
            if (ieee_is_nan(distance(i, j)) .or. ieee_is_nan(distance(j, i))) then
               ok = .false.
               return
            end if
            if (distance(i, j) < 0.0_dp .or. distance(j, i) < 0.0_dp) then
               ok = .false.
               return
            end if
            if (abs(distance(i, j) - distance(j, i)) > tolerance * &
               max(1.0_dp, abs(distance(i, j)), abs(distance(j, i)))) then
               ok = .false.
               return
            end if
         end do
      end do
   end function valid_distance_matrix

   pure subroutine final_three(active, a, b, c, info)
      logical, intent(in) :: active(:) !! Mask identifying the three remaining BIONJ/NJ clusters.
      integer, intent(out) :: a !! Index of the first remaining cluster.
      integer, intent(out) :: b !! Index of the second remaining cluster.
      integer, intent(out) :: c !! Index of the third remaining cluster.
      integer, intent(out) :: info !! Zero if exactly three active clusters are found, otherwise nonzero.
      integer :: i
      integer :: count

      a = 0
      b = 0
      c = 0
      count = 0
      do i = 1, size(active)
         if (.not. active(i)) cycle
         count = count + 1
         if (count == 1) a = i
         if (count == 2) b = i
         if (count == 3) c = i
      end do
      if (count == 3) then
         info = 0
      else
         info = 4
      end if
   end subroutine final_three

   pure elemental logical function is_missing(value) result(missing)
      real(dp), intent(in) :: value !! Distance value; NaN or a negative value is treated as missing.

      missing = ieee_is_nan(value) .or. value < 0.0_dp
   end function is_missing

   pure subroutine normalize_unresolved(matrix, unresolved)
      real(dp), intent(inout) :: matrix(:, :) !! Symmetric completion matrix whose missing-pair count is requested.
      integer, intent(out) :: unresolved !! Number of still-missing unordered off-diagonal pairs.
      integer :: i
      integer :: j

      unresolved = 0
      do i = 1, size(matrix, 1) - 1
         do j = i + 1, size(matrix, 1)
            if (is_missing(matrix(i, j))) unresolved = unresolved + 1
         end do
      end do
   end subroutine normalize_unresolved

end module ape_reconstruction
