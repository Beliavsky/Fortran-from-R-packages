! SPDX-License-Identifier: GPL-2.0-or-later
! Computational translation derived from R package spdep 1.4-2.
! Translated/modified 2026-08-30 for modern free-form Fortran.
module spdep_clustering
   use spdep_kinds, only : dp
   use spdep_types, only : real_vector, neighbor_list, mst_result
   use spdep_math, only : distance_metric, safe_nan
   implicit none
   private

   public :: nbcosts
   public :: ssw
   public :: mstree
   public :: prunecost
   public :: skater_groups

contains

   pure function nbcosts(nb, data, method, p) result(costs)
      type(neighbor_list), intent(in) :: nb !! Neighbor graph whose links receive feature-space dissimilarities.
      real(dp), intent(in) :: data(:, :) !! Feature matrix with one region per row and one variable per column.
      character(len=*), intent(in), optional :: method !! Distance metric accepted by distance_metric; default is euclidean.
      real(dp), intent(in), optional :: p !! Positive Minkowski power used only for Minkowski distance.
      type(real_vector), allocatable :: costs(:)
      integer :: i
      integer :: j
      integer :: k
      integer :: n

      n = nb%size()
      allocate(costs(n))
      if (size(data, 1) /= n) then
         do i = 1, n
            allocate(costs(i)%values(0))
         end do
         return
      end if
      do i = 1, n
         allocate(costs(i)%values(size(nb%neighbors(i)%values)))
         do k = 1, size(nb%neighbors(i)%values)
            j = nb%neighbors(i)%values(k)
            costs(i)%values(k) = distance_metric(data(i, :), data(j, :), method, p)
         end do
      end do
   end function nbcosts

   pure real(dp) function ssw(data, ids, method, p) result(value)
      real(dp), intent(in) :: data(:, :) !! Feature matrix with one region per row.
      integer, intent(in) :: ids(:) !! One-based row indices belonging to the group whose within-group dispersion is required.
      character(len=*), intent(in), optional :: method !! Distance metric from each observation to the group centroid.
      real(dp), intent(in), optional :: p !! Positive Minkowski power used only for Minkowski distance.
      real(dp), allocatable :: centroid(:)
      integer :: i

      if (size(ids) == 0) then
         value = 0.0_dp
         return
      end if
      if (any(ids < 1) .or. any(ids > size(data, 1))) then
         value = safe_nan()
         return
      end if
      allocate(centroid(size(data, 2)))
      centroid = 0.0_dp
      do i = 1, size(ids)
         centroid = centroid + data(ids(i), :)
      end do
      centroid = centroid / real(size(ids), dp)
      value = 0.0_dp
      do i = 1, size(ids)
         value = value + distance_metric(data(ids(i), :), centroid, method, p)
      end do
   end function ssw

   pure function mstree(nb, costs, ini) result(tree)
      type(neighbor_list), intent(in) :: nb !! Neighbor graph constraining the edges available to Prim's minimum spanning tree.
      type(real_vector), intent(in) :: costs(:) !! Edge costs conformable with nb neighbor vectors.
      integer, intent(in), optional :: ini !! One-based starting region; default is 1 for deterministic behavior.
      type(mst_result) :: tree
      logical, allocatable :: chosen(:)
      real(dp), allocatable :: best_cost(:)
      integer, allocatable :: parent(:)
      integer :: n
      integer :: start
      integer :: i
      integer :: j
      integer :: k
      integer :: next
      integer :: edge
      real(dp) :: candidate

      n = nb%size()
      if (n <= 1 .or. size(costs) /= n) then
         allocate(tree%from(0), tree%to(0), tree%cost(0))
         tree%total_cost = 0.0_dp
         return
      end if
      start = 1
      if (present(ini)) start = ini
      if (start < 1 .or. start > n) start = 1
      allocate(tree%from(n - 1), tree%to(n - 1), tree%cost(n - 1))
      allocate(chosen(n), best_cost(n), parent(n))
      chosen = .false.
      best_cost = huge(1.0_dp)
      parent = 0
      chosen(start) = .true.
      do k = 1, size(nb%neighbors(start)%values)
         j = nb%neighbors(start)%values(k)
         if (k <= size(costs(start)%values)) then
            best_cost(j) = costs(start)%values(k)
            parent(j) = start
         end if
      end do
      tree%total_cost = 0.0_dp
      do edge = 1, n - 1
         next = 0
         candidate = huge(1.0_dp)
         do i = 1, n
            if (.not. chosen(i) .and. best_cost(i) < candidate) then
               candidate = best_cost(i)
               next = i
            end if
         end do
         if (next == 0 .or. parent(next) == 0) then
            tree%from(edge:) = 0
            tree%to(edge:) = 0
            tree%cost(edge:) = huge(1.0_dp)
            tree%total_cost = huge(1.0_dp)
            return
         end if
         tree%from(edge) = parent(next)
         tree%to(edge) = next
         tree%cost(edge) = candidate
         tree%total_cost = tree%total_cost + candidate
         chosen(next) = .true.
         do k = 1, size(nb%neighbors(next)%values)
            j = nb%neighbors(next)%values(k)
            if (chosen(j)) cycle
            if (k > size(costs(next)%values)) cycle
            if (costs(next)%values(k) < best_cost(j)) then
               best_cost(j) = costs(next)%values(k)
               parent(j) = next
            end if
         end do
      end do
   end function mstree

   pure function prunecost(tree, data, method, p, min_group_size) result(gain)
      type(mst_result), intent(in) :: tree !! Minimum spanning tree whose individual edges are considered for removal.
      real(dp), intent(in) :: data(:, :) !! Feature matrix with one region per row.
      character(len=*), intent(in), optional :: method !! Distance metric used by within-group dispersion.
      real(dp), intent(in), optional :: p !! Positive Minkowski power used only for Minkowski distance.
      integer, intent(in), optional :: min_group_size !! Minimum size required for both groups created by a cut; default is 1.
      real(dp), allocatable :: gain(:)
      logical, allocatable :: active(:)
      integer, allocatable :: labels(:)
      integer, allocatable :: ids1(:)
      integer, allocatable :: ids2(:)
      integer, allocatable :: all_ids(:)
      integer :: n
      integer :: e
      integer :: min_size
      integer :: n1
      integer :: n2
      integer :: idx
      real(dp) :: total

      n = size(data, 1)
      min_size = 1
      if (present(min_group_size)) min_size = max(1, min_group_size)
      allocate(gain(size(tree%from)), active(size(tree%from)), all_ids(n))
      all_ids = [(e, e = 1, n)]
      total = ssw(data, all_ids, method, p)
      active = .true.
      do e = 1, size(tree%from)
         active(e) = .false.
         labels = tree_components(tree, active, n)
         n1 = count(labels == labels(tree%from(e)))
         n2 = count(labels == labels(tree%to(e)))
         if (n1 < min_size .or. n2 < min_size) then
            gain(e) = -huge(1.0_dp)
         else
            ids1 = pack([(idx, idx = 1, n)], labels == labels(tree%from(e)))
            ids2 = pack([(idx, idx = 1, n)], labels == labels(tree%to(e)))
            gain(e) = total - ssw(data, ids1, method, p) - ssw(data, ids2, method, p)
         end if
         active(e) = .true.
      end do
   end function prunecost

   pure function skater_groups(tree, data, ncuts, method, p, min_group_size) result(groups)
      type(mst_result), intent(in) :: tree !! Minimum spanning tree to prune into contiguous feature-homogeneous groups.
      real(dp), intent(in) :: data(:, :) !! Feature matrix with one region per row.
      integer, intent(in) :: ncuts !! Number of tree edges to remove; resulting group count is at most ncuts+1.
      character(len=*), intent(in), optional :: method !! Distance metric used by within-group dispersion.
      real(dp), intent(in), optional :: p !! Positive Minkowski power used only for Minkowski distance.
      integer, intent(in), optional :: min_group_size !! Minimum component size allowed after each cut; default is 1.
      integer, allocatable :: groups(:)
      logical, allocatable :: active(:)
      integer, allocatable :: candidate_groups(:)
      integer :: n
      integer :: e
      integer :: cut
      integer :: best_edge
      integer :: min_size
      real(dp) :: best_gain
      real(dp) :: current_obj
      real(dp) :: new_obj

      n = size(data, 1)
      min_size = 1
      if (present(min_group_size)) min_size = max(1, min_group_size)
      allocate(active(size(tree%from)))
      active = .true.
      groups = tree_components(tree, active, n)
      current_obj = grouped_ssw(data, groups, method, p)
      do cut = 1, max(0, ncuts)
         best_edge = 0
         best_gain = -huge(1.0_dp)
         do e = 1, size(tree%from)
            if (.not. active(e)) cycle
            active(e) = .false.
            candidate_groups = tree_components(tree, active, n)
            if (minimum_group_count(candidate_groups) >= min_size) then
               new_obj = grouped_ssw(data, candidate_groups, method, p)
               if (current_obj - new_obj > best_gain) then
                  best_gain = current_obj - new_obj
                  best_edge = e
               end if
            end if
            active(e) = .true.
         end do
         if (best_edge == 0) exit
         active(best_edge) = .false.
         groups = tree_components(tree, active, n)
         current_obj = grouped_ssw(data, groups, method, p)
      end do
   end function skater_groups

   pure function tree_components(tree, active, n) result(labels)
      type(mst_result), intent(in) :: tree !! Tree edge list whose active subset defines the component structure.
      logical, intent(in) :: active(:) !! Per-edge mask; true edges remain connected.
      integer, intent(in) :: n !! Number of tree vertices to label.
      integer, allocatable :: labels(:)
      logical, allocatable :: adj(:, :)
      integer, allocatable :: stack(:)
      logical, allocatable :: seen(:)
      integer :: e
      integer :: i
      integer :: j
      integer :: v
      integer :: top
      integer :: label

      allocate(labels(n), adj(n, n), stack(max(1, n)), seen(n))
      labels = 0
      adj = .false.
      seen = .false.
      do e = 1, min(size(active), size(tree%from))
         if (.not. active(e)) cycle
         i = tree%from(e)
         j = tree%to(e)
         if (i < 1 .or. i > n .or. j < 1 .or. j > n) cycle
         adj(i, j) = .true.
         adj(j, i) = .true.
      end do
      label = 0
      do i = 1, n
         if (seen(i)) cycle
         label = label + 1
         top = 1
         stack(1) = i
         seen(i) = .true.
         do while (top > 0)
            v = stack(top)
            top = top - 1
            labels(v) = label
            do j = 1, n
               if (adj(v, j) .and. .not. seen(j)) then
                  top = top + 1
                  stack(top) = j
                  seen(j) = .true.
               end if
            end do
         end do
      end do
   end function tree_components

   pure real(dp) function grouped_ssw(data, groups, method, p) result(value)
      real(dp), intent(in) :: data(:, :) !! Feature matrix with one observation per row.
      integer, intent(in) :: groups(:) !! Positive component label for each observation.
      character(len=*), intent(in), optional :: method !! Distance metric used within each component.
      real(dp), intent(in), optional :: p !! Positive Minkowski power used only for Minkowski distance.
      integer, allocatable :: ids(:)
      integer :: g
      integer :: i
      integer :: max_group

      if (size(groups) /= size(data, 1) .or. size(groups) == 0) then
         value = 0.0_dp
         return
      end if
      max_group = maxval(groups)
      value = 0.0_dp
      do g = 1, max_group
         ids = pack([(i, i = 1, size(groups))], groups == g)
         value = value + ssw(data, ids, method, p)
      end do
   end function grouped_ssw

   pure integer function minimum_group_count(groups) result(min_count)
      integer, intent(in) :: groups(:) !! Positive component label for each observation.
      integer :: g
      integer :: max_group

      if (size(groups) == 0) then
         min_count = 0
         return
      end if
      max_group = maxval(groups)
      min_count = size(groups)
      do g = 1, max_group
         min_count = min(min_count, count(groups == g))
      end do
   end function minimum_group_count

end module spdep_clustering
