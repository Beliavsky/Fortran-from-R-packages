! SPDX-License-Identifier: MIT
! SPDX-FileComment: Clustering adapters corresponding to selected R stats functions.
module r_stats_clustering
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use fastcluster_core, only: fast_hclust_condensed => hclust_condensed, &
                               fast_hclust_matrix => hclust_matrix
   use fastcluster_distances, only: matrix_to_condensed, pairwise_distances
   use fastcluster_types, only: hclust_result_t => hclust_result
   use r_kinds, only: dp
   use r_optional, only: optval
   use r_stats_types, only: kmeans_result_t
   implicit none
   private

   public :: cutree, dist, dist_matrix, hclust, hclust_condensed, hclust_result_t
   public :: kmeans, kmeans_matrix, kmeans_vector

   interface kmeans
      module procedure kmeans_matrix, kmeans_vector
   end interface kmeans

contains

   pure function dist(x, method, p) result(condensed)
      !! Returns pairwise distances in R's column-wise lower-triangle layout.
      real(dp), intent(in) :: x(:, :) !! Observations by rows and variables by columns.
      character(len=*), intent(in), optional :: method !! Distance method; defaults to Euclidean.
      real(dp), intent(in), optional :: p !! Positive Minkowski exponent.
      real(dp), allocatable :: condensed(:)
      real(dp), allocatable :: distances(:, :)
      character(len=:), allocatable :: metric
      integer :: status

      metric = "euclidean"
      if (present(method)) metric = trim(adjustl(method))
      call pairwise_distances(x, metric, distances, p=p, status=status)
      if (status /= 0) then
         allocate (condensed(0))
         return
      end if
      call matrix_to_condensed(distances, condensed, status=status)
      if (status /= 0) then
         if (allocated(condensed)) deallocate (condensed)
         allocate (condensed(0))
      end if
   end function dist

   pure function dist_matrix(x, method, p) result(distances)
      !! Returns a symmetric matrix of pairwise distances.
      real(dp), intent(in) :: x(:, :) !! Observations by rows and variables by columns.
      character(len=*), intent(in), optional :: method !! Distance method; defaults to Euclidean.
      real(dp), intent(in), optional :: p !! Positive Minkowski exponent.
      real(dp), allocatable :: distances(:, :)
      character(len=:), allocatable :: metric

      metric = "euclidean"
      if (present(method)) metric = trim(adjustl(method))
      call pairwise_distances(x, metric, distances, p=p)
   end function dist_matrix

   pure function kmeans_matrix(x, centers, iter_max) result(fit)
      !! Fits Lloyd k-means deterministically from explicitly supplied initial centers.
      real(dp), intent(in) :: x(:, :) !! Observations with shape `(n, p)`.
      real(dp), intent(in) :: centers(:, :) !! Initial centers with shape `(k, p)`.
      integer, intent(in), optional :: iter_max !! Positive iteration limit; defaults to ten.
      type(kmeans_result_t) :: fit
      integer, allocatable :: assignments(:), counts(:), previous_assignments(:)
      real(dp), allocatable :: current_centers(:, :), new_centers(:, :), overall_center(:)
      integer :: cluster_count, i, iteration, j, limit, observation_count, variable_count
      logical :: converged

      observation_count = size(x, 1)
      variable_count = size(x, 2)
      cluster_count = size(centers, 1)
      if (observation_count == 0 .or. variable_count == 0 .or. cluster_count == 0) then
         fit%status = 1
         return
      end if
      if (size(centers, 2) /= variable_count .or. cluster_count > observation_count) then
         fit%status = 1
         return
      end if
      if (any(.not. ieee_is_finite(x)) .or. any(.not. ieee_is_finite(centers))) then
         fit%status = 1
         return
      end if

      limit = max(1, optval(iter_max, 10))
      allocate (assignments(observation_count), previous_assignments(observation_count), source=0)
      allocate (counts(cluster_count), new_centers(cluster_count, variable_count))
      current_centers = centers
      converged = .false.
      do iteration = 1, limit
         counts = 0
         new_centers = 0.0_dp
         do i = 1, observation_count
            assignments(i) = nearest_center(x(i, :), current_centers)
            counts(assignments(i)) = counts(assignments(i)) + 1
            new_centers(assignments(i), :) = new_centers(assignments(i), :) + x(i, :)
         end do
         if (any(counts == 0)) then
            fit%status = 3
            return
         end if
         do j = 1, cluster_count
            new_centers(j, :) = new_centers(j, :) / real(counts(j), dp)
         end do
         converged = all(assignments == previous_assignments)
         previous_assignments = assignments
         current_centers = new_centers
         fit%iter = iteration
         if (converged) exit
      end do
      if (.not. converged) fit%status = 2

      fit%centers = current_centers
      fit%cluster = assignments
      fit%size = counts
      allocate (fit%withinss(cluster_count), source=0.0_dp)
      do i = 1, observation_count
         j = fit%cluster(i)
         fit%withinss(j) = fit%withinss(j) + sum((x(i, :) - fit%centers(j, :))**2)
      end do
      overall_center = sum(x, dim=1) / real(observation_count, dp)
      fit%totss = 0.0_dp
      do i = 1, observation_count
         fit%totss = fit%totss + sum((x(i, :) - overall_center)**2)
      end do
      fit%tot_withinss = sum(fit%withinss)
      fit%betweenss = fit%totss - fit%tot_withinss
   end function kmeans_matrix

   pure function kmeans_vector(x, centers, iter_max) result(fit)
      !! Fits deterministic Lloyd k-means to one-dimensional observations.
      real(dp), intent(in) :: x(:) !! One-dimensional observations.
      real(dp), intent(in) :: centers(:) !! Explicit initial cluster centers.
      integer, intent(in), optional :: iter_max !! Positive iteration limit; defaults to ten.
      type(kmeans_result_t) :: fit
      real(dp), allocatable :: center_matrix(:, :), data_matrix(:, :)

      allocate (data_matrix(size(x), 1), center_matrix(size(centers), 1))
      data_matrix(:, 1) = x
      center_matrix(:, 1) = centers
      fit = kmeans_matrix(data_matrix, center_matrix, iter_max)
   end function kmeans_vector

   pure function nearest_center(observation, centers) result(index)
      !! Returns the first center having minimum squared Euclidean distance.
      real(dp), intent(in) :: observation(:) !! Observation with shape `(p)`.
      real(dp), intent(in) :: centers(:, :) !! Candidate centers with shape `(k, p)`.
      integer :: index
      real(dp) :: candidate, squared_distance
      integer :: j

      index = 1
      squared_distance = sum((observation - centers(1, :))**2)
      do j = 2, size(centers, 1)
         candidate = sum((observation - centers(j, :))**2)
         if (candidate < squared_distance) then
            squared_distance = candidate
            index = j
         end if
      end do
   end function nearest_center

   pure function hclust(distances, method) result(tree)
      !! Performs hierarchical clustering on a symmetric distance matrix.
      real(dp), intent(in) :: distances(:, :) !! Symmetric pairwise distance matrix.
      character(len=*), intent(in), optional :: method !! Linkage method; defaults to complete.
      type(hclust_result_t) :: tree
      character(len=:), allocatable :: linkage

      linkage = "complete"
      if (present(method)) linkage = trim(adjustl(method))
      call fast_hclust_matrix(distances, linkage, tree)
   end function hclust

   pure function hclust_condensed(distances, n, method) result(tree)
      !! Performs hierarchical clustering directly from R-layout condensed distances.
      real(dp), intent(in) :: distances(:) !! Condensed lower-triangle distances.
      integer, intent(in) :: n !! Number of observations represented by `distances`.
      character(len=*), intent(in), optional :: method !! Linkage method; defaults to complete.
      type(hclust_result_t) :: tree
      character(len=:), allocatable :: linkage

      linkage = "complete"
      if (present(method)) linkage = trim(adjustl(method))
      call fast_hclust_condensed(distances, n, linkage, tree)
   end function hclust_condensed

   pure function cutree(tree, k, h) result(groups)
      !! Cuts a hierarchical tree at a requested cluster count or merge height.
      type(hclust_result_t), intent(in) :: tree !! Hierarchical clustering result.
      integer, intent(in), optional :: k !! Requested cluster count.
      real(dp), intent(in), optional :: h !! Requested cut height.
      integer, allocatable :: groups(:)
      integer, allocatable :: merge_representative(:), parent(:), representative_label(:)
      integer :: first, first_root, i, merges_to_apply, n, next_label, second, second_root, root

      n = tree%n
      allocate (groups(max(0, n)))
      if (n <= 0) return
      allocate (parent(n), representative_label(n), merge_representative(max(0, n - 1)))
      do i = 1, n
         parent(i) = i
      end do
      merge_representative = 0
      if (present(h)) then
         merges_to_apply = count(tree%height <= h)
      else
         merges_to_apply = n - max(1, min(optval(k, 1), n))
      end if
      merges_to_apply = min(merges_to_apply, size(tree%merge, 1))

      do i = 1, merges_to_apply
         first = tree%merge(i, 1)
         second = tree%merge(i, 2)
         if (first < 0) then
            first_root = -first
         else
            first_root = merge_representative(first)
         end if
         if (second < 0) then
            second_root = -second
         else
            second_root = merge_representative(second)
         end if
         first_root = find_root(parent, first_root)
         second_root = find_root(parent, second_root)
         if (first_root /= second_root) parent(second_root) = first_root
         merge_representative(i) = first_root
      end do

      representative_label = 0
      next_label = 0
      do i = 1, n
         root = find_root(parent, i)
         if (representative_label(root) == 0) then
            next_label = next_label + 1
            representative_label(root) = next_label
         end if
         groups(i) = representative_label(root)
      end do
   end function cutree

   pure function find_root(parent, node) result(root)
      !! Returns the representative of a union-find node without modifying the parent array.
      integer, intent(in) :: parent(:) !! Union-find parent links.
      integer, intent(in) :: node !! One-based node whose representative is requested.
      integer :: root

      root = node
      do while (parent(root) /= root)
         root = parent(root)
      end do
   end function find_root

end module r_stats_clustering
