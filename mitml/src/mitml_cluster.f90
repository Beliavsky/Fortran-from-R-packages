! SPDX-License-Identifier: GPL-2.0-or-later
! Upstream mitml 0.4-5 (2023-03-08), authored by Simon Grund,
! Alexander Robitzsch, and Oliver Luedtke; upstream license GPL (>= 2).
! Modern free-form Fortran translation for Fortran-from-R-packages.
! Numeric clusterMeans translation from mitml.
module mitml_cluster
   use, intrinsic :: ieee_arithmetic, only : ieee_is_nan, ieee_quiet_nan, ieee_value
   use r_kinds, only : dp
   use mitml_types, only : MITML_ERR_DIMENSION, MITML_OK
   implicit none
   private

   public :: cluster_means
   public :: cluster_means_matrix

contains

   subroutine cluster_means(x, cluster, means, status, adjusted, group)
      real(dp), intent(in) :: x(:) !! Values whose within-cluster means are requested; NaNs are treated as missing.
      integer, intent(in) :: cluster(:) !! Cluster label for each observation; labels need not be consecutive.
      real(dp), intent(out) :: means(:) !! Observation-aligned cluster means or leave-one-out means.
      integer, intent(out) :: status !! MITML_OK on success or a dimension error code.
      logical, intent(in), optional :: adjusted !! If true, return leave-one-out cluster means.
      integer, intent(in), optional :: group(:) !! Optional group label; clusters are nested within group for averaging.
      integer, allocatable :: count_obs(:)
      integer, allocatable :: membership(:)
      real(dp), allocatable :: sums(:)
      real(dp) :: nan
      integer :: i
      integer :: j
      integer :: n
      integer :: n_cluster
      logical :: leave_one_out

      status = MITML_OK
      n = size(x)
      if (size(cluster) /= n .or. size(means) /= n) then
         status = MITML_ERR_DIMENSION
         return
      end if
      if (present(group)) then
         if (size(group) /= n) then
            status = MITML_ERR_DIMENSION
            return
         end if
      end if
      leave_one_out = .false.
      if (present(adjusted)) leave_one_out = adjusted
      nan = ieee_value(0.0_dp, ieee_quiet_nan)

      allocate(membership(n))
      membership = 0
      n_cluster = 0
      do i = 1, n
         do j = 1, i - 1
            if (cluster(j) /= cluster(i)) cycle
            if (present(group)) then
               if (group(j) /= group(i)) cycle
            end if
            membership(i) = membership(j)
            exit
         end do
         if (membership(i) == 0) then
            n_cluster = n_cluster + 1
            membership(i) = n_cluster
         end if
      end do

      allocate(count_obs(n_cluster), sums(n_cluster))
      count_obs = 0
      sums = 0.0_dp
      do i = 1, n
         if (.not. ieee_is_nan(x(i))) then
            count_obs(membership(i)) = count_obs(membership(i)) + 1
            sums(membership(i)) = sums(membership(i)) + x(i)
         end if
      end do

      do i = 1, n
         j = membership(i)
         if (.not. leave_one_out) then
            if (count_obs(j) > 0) then
               means(i) = sums(j) / real(count_obs(j), dp)
            else
               means(i) = nan
            end if
         else
            if (ieee_is_nan(x(i)) .or. count_obs(j) <= 1) then
               means(i) = nan
            else
               means(i) = (sums(j) - x(i)) / real(count_obs(j) - 1, dp)
            end if
         end if
      end do
   end subroutine cluster_means

   subroutine cluster_means_matrix(x, cluster, means, status, adjusted, group)
      real(dp), intent(in) :: x(:, :) !! Data matrix with observations in rows and variables in columns.
      integer, intent(in) :: cluster(:) !! Cluster label for each row of x.
      real(dp), intent(out) :: means(:, :) !! Cluster means with the same shape as x.
      integer, intent(out) :: status !! MITML_OK on success or a dimension error code.
      logical, intent(in), optional :: adjusted !! If true, compute leave-one-out means for every column.
      integer, intent(in), optional :: group(:) !! Optional group label nested with cluster labels.
      integer :: j
      integer :: local_status

      status = MITML_OK
      if (any(shape(means) /= shape(x)) .or. size(cluster) /= size(x, 1)) then
         status = MITML_ERR_DIMENSION
         return
      end if
      if (present(group)) then
         if (size(group) /= size(x, 1)) then
            status = MITML_ERR_DIMENSION
            return
         end if
      end if
      do j = 1, size(x, 2)
         if (present(group) .and. present(adjusted)) then
            call cluster_means(x(:, j), cluster, means(:, j), local_status, adjusted, group)
         else if (present(group)) then
            call cluster_means(x(:, j), cluster, means(:, j), local_status, group=group)
         else if (present(adjusted)) then
            call cluster_means(x(:, j), cluster, means(:, j), local_status, adjusted)
         else
            call cluster_means(x(:, j), cluster, means(:, j), local_status)
         end if
         if (local_status /= MITML_OK) then
            status = local_status
            return
         end if
      end do
   end subroutine cluster_means_matrix

end module mitml_cluster
