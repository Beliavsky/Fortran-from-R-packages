! SPDX-License-Identifier: GPL-2.0-only OR GPL-3.0-only
module sandwich_panel
   use sandwich_kinds, only : dp
   use sandwich_status, only : SANDWICH_SUCCESS, SANDWICH_INVALID_ARGUMENT, &
      SANDWICH_DIMENSION_MISMATCH, SANDWICH_INSUFFICIENT_DATA
   use sandwich_utils, only : sort_index_integer, unique_values_integer
   use sandwich_kernels, only : kernel_weights
   use sandwich_linalg, only : outer_product, project_psd
   use sandwich_core, only : sandwich_covariance
   implicit none
   private

   integer, parameter, public :: PANEL_LAG_MAX = -1
   integer, parameter, public :: PANEL_LAG_NW1987 = -2
   integer, parameter, public :: PANEL_LAG_NW1994 = -3

   public :: meat_panel_longitudinal, vcov_panel_longitudinal
   public :: meat_panel_corrected, vcov_panel_corrected

contains

   subroutine meat_panel_longitudinal(scores, cluster, time, meat_matrix, status, kernel, &
      lag, bandwidth, adjust, aggregate)
      real(dp), intent(in) :: scores(:, :)
      integer, intent(in) :: cluster(:), time(:)
      real(dp), allocatable, intent(out) :: meat_matrix(:, :)
      integer, intent(out), optional :: status
      character(len=*), intent(in), optional :: kernel
      integer, intent(in), optional :: lag
      real(dp), intent(in), optional :: bandwidth
      logical, intent(in), optional :: adjust, aggregate
      real(dp), allocatable :: ordered_scores(:, :), ef(:, :), grid(:, :, :)
      real(dp), allocatable :: weight_arguments(:), weights(:), term(:, :)
      integer, allocatable :: index(:), ordered_time(:), ordered_cluster(:)
      integer, allocatable :: unique_time(:), time_label(:), unique_cluster(:), cluster_label(:)
      logical, allocatable :: present_pair(:, :)
      integer :: n, k, nt, nc, selected_lag, nw, i, j, t, c, info
      real(dp) :: bw
      logical :: use_adjust, use_aggregate
      character(len=32) :: selected_kernel

      n = size(scores, 1)
      k = size(scores, 2)
      if (n <= 0 .or. k <= 0 .or. size(cluster) /= n .or. size(time) /= n) then
         allocate(meat_matrix(0, 0))
         if (present(status)) status = SANDWICH_DIMENSION_MISMATCH
         return
      end if
      use_adjust = .true.
      if (present(adjust)) use_adjust = adjust
      use_aggregate = .true.
      if (present(aggregate)) use_aggregate = aggregate
      if (use_adjust .and. n <= k) then
         allocate(meat_matrix(0, 0))
         if (present(status)) status = SANDWICH_INVALID_ARGUMENT
         return
      end if
      selected_kernel = 'Bartlett'
      if (present(kernel)) selected_kernel = kernel

      call sort_index_integer(time, index)
      allocate(ordered_scores(n, k), ordered_time(n), ordered_cluster(n))
      do i = 1, n
         ordered_scores(i, :) = scores(index(i), :)
         ordered_time(i) = time(index(i))
         ordered_cluster(i) = cluster(index(i))
      end do
      call unique_values_integer(ordered_time, unique_time, time_label)
      call unique_values_integer(ordered_cluster, unique_cluster, cluster_label)
      nt = size(unique_time)
      nc = size(unique_cluster)

      if (present(lag)) then
         selected_lag = lag
      else
         selected_lag = PANEL_LAG_NW1987
      end if
      select case (selected_lag)
      case (PANEL_LAG_MAX)
         selected_lag = nt - 1
      case (PANEL_LAG_NW1987)
         selected_lag = floor(real(nt, dp)**0.25_dp)
      case (PANEL_LAG_NW1994)
         selected_lag = floor(4.0_dp * (real(nt, dp) / 100.0_dp)**(2.0_dp / 9.0_dp))
      case default
         continue
      end select
      if (selected_lag < 0 .or. selected_lag >= nt) then
         allocate(meat_matrix(0, 0))
         if (present(status)) status = SANDWICH_INVALID_ARGUMENT
         return
      end if
      if (present(bandwidth)) then
         bw = bandwidth
      else
         bw = real(selected_lag + 1, dp)
      end if
      if (bw <= 0.0_dp) then
         allocate(meat_matrix(0, 0))
         if (present(status)) status = SANDWICH_INVALID_ARGUMENT
         return
      end if

      allocate(weight_arguments(nt))
      do i = 1, nt
         weight_arguments(i) = real(i - 1, dp) / bw
      end do
      call kernel_weights(weight_arguments, selected_kernel, weights, info)
      if (info /= SANDWICH_SUCCESS) then
         allocate(meat_matrix(0, 0))
         if (present(status)) status = info
         return
      end if
      nw = 1
      do i = 1, nt
         if (abs(weights(i)) > tiny(1.0_dp)) nw = i
      end do
      nw = min(nw, selected_lag + 1)

      if (use_aggregate) then
         allocate(ef(nt, k))
         ef = 0.0_dp
         do i = 1, n
            ef(time_label(i), :) = ef(time_label(i), :) + ordered_scores(i, :)
         end do
         meat_matrix = 0.5_dp * weights(1) * matmul(transpose(ef), ef)
         do j = 1, nw - 1
            term = matmul(transpose(ef(1:nt - j, :)), ef(1 + j:nt, :))
            meat_matrix = meat_matrix + weights(j + 1) * term
         end do
      else
         allocate(grid(nt, nc, k), present_pair(nt, nc))
         grid = 0.0_dp
         present_pair = .false.
         do i = 1, n
            t = time_label(i)
            c = cluster_label(i)
            grid(t, c, :) = grid(t, c, :) + ordered_scores(i, :)
            present_pair(t, c) = .true.
         end do
         allocate(meat_matrix(k, k))
         meat_matrix = 0.0_dp
         do t = 1, nt
            do c = 1, nc
               if (present_pair(t, c)) then
                  meat_matrix = meat_matrix + 0.5_dp * weights(1) * &
                     outer_product(grid(t, c, :), grid(t, c, :))
               end if
            end do
         end do
         do j = 1, nw - 1
            do t = 1, nt - j
               do c = 1, nc
                  if (present_pair(t, c) .and. present_pair(t + j, c)) then
                     meat_matrix = meat_matrix + weights(j + 1) * &
                        outer_product(grid(t, c, :), grid(t + j, c, :))
                  end if
               end do
            end do
         end do
      end if

      meat_matrix = meat_matrix + transpose(meat_matrix)
      if (use_adjust) meat_matrix = real(n, dp) / real(n - k, dp) * meat_matrix
      meat_matrix = meat_matrix / real(n, dp)
      meat_matrix = 0.5_dp * (meat_matrix + transpose(meat_matrix))
      if (present(status)) status = SANDWICH_SUCCESS
   end subroutine meat_panel_longitudinal

   subroutine vcov_panel_longitudinal(scores, cluster, time, bread, covariance, status, &
      kernel, lag, bandwidth, adjust, aggregate, fix)
      real(dp), intent(in) :: scores(:, :), bread(:, :)
      integer, intent(in) :: cluster(:), time(:)
      real(dp), allocatable, intent(out) :: covariance(:, :)
      integer, intent(out), optional :: status
      character(len=*), intent(in), optional :: kernel
      integer, intent(in), optional :: lag
      real(dp), intent(in), optional :: bandwidth
      logical, intent(in), optional :: adjust, aggregate, fix
      real(dp), allocatable :: meat_matrix(:, :), fixed_covariance(:, :)
      character(len=32) :: selected_kernel
      integer :: selected_lag, info
      real(dp) :: bw
      logical :: use_adjust, use_aggregate, use_fix

      selected_kernel = 'Bartlett'
      if (present(kernel)) selected_kernel = kernel
      selected_lag = PANEL_LAG_NW1987
      if (present(lag)) selected_lag = lag
      use_adjust = .true.
      if (present(adjust)) use_adjust = adjust
      use_aggregate = .true.
      if (present(aggregate)) use_aggregate = aggregate
      use_fix = .false.
      if (present(fix)) use_fix = fix
      if (present(bandwidth)) then
         bw = bandwidth
         call meat_panel_longitudinal(scores, cluster, time, meat_matrix, info, &
            selected_kernel, selected_lag, bw, use_adjust, use_aggregate)
      else
         call meat_panel_longitudinal(scores, cluster, time, meat_matrix, info, &
            selected_kernel, selected_lag, adjust = use_adjust, aggregate = use_aggregate)
      end if
      if (info /= SANDWICH_SUCCESS) then
         allocate(covariance(0, 0))
         if (present(status)) status = info
         return
      end if
      call sandwich_covariance(bread, meat_matrix, size(scores, 1), covariance, info)
      if (info == SANDWICH_SUCCESS .and. use_fix) then
         call project_psd(covariance, fixed_covariance, info)
         if (info == SANDWICH_SUCCESS) call move_alloc(fixed_covariance, covariance)
      end if
      if (present(status)) status = info
   end subroutine vcov_panel_longitudinal

   subroutine meat_panel_corrected(x, residuals, cluster, time, meat_matrix, status, pairwise)
      real(dp), intent(in) :: x(:, :), residuals(:)
      integer, intent(in) :: cluster(:), time(:)
      real(dp), allocatable, intent(out) :: meat_matrix(:, :)
      integer, intent(out), optional :: status
      logical, intent(in), optional :: pairwise
      real(dp), allocatable :: egrid(:, :), xgrid(:, :, :), sigma(:, :)
      integer, allocatable :: unique_cluster(:), cluster_label(:), unique_time(:), time_label(:)
      logical, allocatable :: present_pair(:, :), complete_time(:)
      real(dp) :: numerator
      integer :: n, k, nt, nc, i, t, c, d, denominator
      logical :: use_pairwise

      n = size(x, 1)
      k = size(x, 2)
      if (n <= 0 .or. k <= 0 .or. size(residuals) /= n .or. &
          size(cluster) /= n .or. size(time) /= n) then
         allocate(meat_matrix(0, 0))
         if (present(status)) status = SANDWICH_DIMENSION_MISMATCH
         return
      end if
      use_pairwise = .false.
      if (present(pairwise)) use_pairwise = pairwise
      call unique_values_integer(cluster, unique_cluster, cluster_label)
      call unique_values_integer(time, unique_time, time_label)
      nc = size(unique_cluster)
      nt = size(unique_time)
      if (nc <= 0 .or. nt <= 0) then
         allocate(meat_matrix(0, 0))
         if (present(status)) status = SANDWICH_INSUFFICIENT_DATA
         return
      end if

      allocate(egrid(nt, nc), xgrid(nt, nc, k), present_pair(nt, nc), complete_time(nt))
      egrid = 0.0_dp
      xgrid = 0.0_dp
      present_pair = .false.
      do i = 1, n
         t = time_label(i)
         c = cluster_label(i)
         egrid(t, c) = residuals(i)
         xgrid(t, c, :) = x(i, :)
         present_pair(t, c) = .true.
      end do
      do t = 1, nt
         complete_time(t) = all(present_pair(t, :))
      end do

      allocate(sigma(nc, nc))
      sigma = 0.0_dp
      do c = 1, nc
         do d = 1, nc
            numerator = 0.0_dp
            denominator = 0
            do t = 1, nt
               if (use_pairwise) then
                  if (present_pair(t, c) .and. present_pair(t, d)) then
                     numerator = numerator + egrid(t, c) * egrid(t, d)
                     denominator = denominator + 1
                  end if
               else if (complete_time(t)) then
                  numerator = numerator + egrid(t, c) * egrid(t, d)
                  denominator = denominator + 1
               end if
            end do
            if (denominator <= 0) then
               allocate(meat_matrix(0, 0))
               if (present(status)) status = SANDWICH_INSUFFICIENT_DATA
               return
            end if
            sigma(c, d) = numerator / real(denominator, dp)
         end do
      end do

      allocate(meat_matrix(k, k))
      meat_matrix = 0.0_dp
      do t = 1, nt
         do c = 1, nc
            if (.not. present_pair(t, c)) cycle
            do d = 1, nc
               if (.not. present_pair(t, d)) cycle
               meat_matrix = meat_matrix + sigma(c, d) * &
                  outer_product(xgrid(t, c, :), xgrid(t, d, :))
            end do
         end do
      end do
      meat_matrix = meat_matrix / real(n, dp)
      meat_matrix = 0.5_dp * (meat_matrix + transpose(meat_matrix))
      if (present(status)) status = SANDWICH_SUCCESS
   end subroutine meat_panel_corrected

   subroutine vcov_panel_corrected(x, residuals, cluster, time, bread, covariance, status, &
      pairwise, fix)
      real(dp), intent(in) :: x(:, :), residuals(:), bread(:, :)
      integer, intent(in) :: cluster(:), time(:)
      real(dp), allocatable, intent(out) :: covariance(:, :)
      integer, intent(out), optional :: status
      logical, intent(in), optional :: pairwise, fix
      real(dp), allocatable :: meat_matrix(:, :), fixed_covariance(:, :)
      integer :: info
      logical :: use_pairwise, use_fix

      use_pairwise = .false.
      if (present(pairwise)) use_pairwise = pairwise
      use_fix = .false.
      if (present(fix)) use_fix = fix
      call meat_panel_corrected(x, residuals, cluster, time, meat_matrix, info, use_pairwise)
      if (info /= SANDWICH_SUCCESS) then
         allocate(covariance(0, 0))
         if (present(status)) status = info
         return
      end if
      call sandwich_covariance(bread, meat_matrix, size(x, 1), covariance, info)
      if (info == SANDWICH_SUCCESS .and. use_fix) then
         call project_psd(covariance, fixed_covariance, info)
         if (info == SANDWICH_SUCCESS) call move_alloc(fixed_covariance, covariance)
      end if
      if (present(status)) status = info
   end subroutine vcov_panel_corrected

end module sandwich_panel
