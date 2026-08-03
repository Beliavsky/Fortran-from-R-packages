! SPDX-License-Identifier: GPL-2.0-only OR GPL-3.0-only
module sandwich_cluster
   use sandwich_kinds, only : dp
   use sandwich_status, only : SANDWICH_SUCCESS, SANDWICH_INVALID_ARGUMENT, &
      SANDWICH_DIMENSION_MISMATCH
   use sandwich_utils, only : lowercase, count_bits, group_subset, aggregate_rows
   use sandwich_linalg, only : inverse_matrix, identity_matrix, symmetric_matrix_power, project_psd
   use sandwich_regression, only : ols_hatvalues
   use sandwich_core, only : sandwich_covariance
   implicit none
   private

   public :: meat_cluster, vcov_cluster

contains

   subroutine meat_cluster(scores, cluster, meat_matrix, status, type, cadjust, multi0, &
      x, residuals, hat)
      real(dp), intent(in) :: scores(:, :)
      integer, intent(in) :: cluster(:, :)
      real(dp), allocatable, intent(out) :: meat_matrix(:, :)
      integer, intent(out), optional :: status
      character(len=*), intent(in), optional :: type
      logical, intent(in), optional :: cadjust, multi0
      real(dp), intent(in), optional :: x(:, :), residuals(:), hat(:)
      real(dp), allocatable :: efi(:, :), aggregated(:, :), invinfo(:, :), leverage(:)
      real(dp), allocatable :: xg(:, :), hg(:, :), adjustment_matrix(:, :), resg(:), adjusted_res(:)
      real(dp) :: adj, scale
      integer, allocatable :: labels(:), indices(:)
      integer :: n, k, p, full_mask, mask, nset, sign_value, g, group, i, j, info
      integer :: group_size, position
      logical :: use_cadjust, use_multi0
      character(len=8) :: kind

      n = size(scores, 1)
      k = size(scores, 2)
      p = size(cluster, 2)
      if (n <= 0 .or. k <= 0 .or. p <= 0 .or. size(cluster, 1) /= n) then
         allocate(meat_matrix(0, 0))
         if (present(status)) status = SANDWICH_DIMENSION_MISMATCH
         return
      end if
      kind = 'hc0'
      if (present(type)) kind = trim(lowercase(type))
      if (trim(kind) == 'hc') kind = 'hc0'
      if (trim(kind) /= 'hc0' .and. trim(kind) /= 'hc1' .and. &
          trim(kind) /= 'hc2' .and. trim(kind) /= 'hc3') then
         allocate(meat_matrix(0, 0))
         if (present(status)) status = SANDWICH_INVALID_ARGUMENT
         return
      end if
      use_cadjust = .true.
      if (present(cadjust)) use_cadjust = cadjust
      use_multi0 = .false.
      if (present(multi0)) use_multi0 = multi0

      if (trim(kind) == 'hc1' .and. n <= k) then
         allocate(meat_matrix(0, 0))
         if (present(status)) status = SANDWICH_INVALID_ARGUMENT
         return
      end if

      if (trim(kind) == 'hc2' .or. trim(kind) == 'hc3') then
         if (.not. present(x) .or. .not. present(residuals)) then
            allocate(meat_matrix(0, 0))
            if (present(status)) status = SANDWICH_INVALID_ARGUMENT
            return
         end if
         if (size(x, 1) /= n .or. size(x, 2) /= k .or. size(residuals) /= n) then
            allocate(meat_matrix(0, 0))
            if (present(status)) status = SANDWICH_DIMENSION_MISMATCH
            return
         end if
         call inverse_matrix(matmul(transpose(x), x), invinfo, info)
         if (info /= SANDWICH_SUCCESS) then
            allocate(meat_matrix(0, 0))
            if (present(status)) status = info
            return
         end if
         if (present(hat)) then
            if (size(hat) /= n) then
               allocate(meat_matrix(0, 0))
               if (present(status)) status = SANDWICH_DIMENSION_MISMATCH
               return
            end if
            allocate(leverage(n))
            leverage = hat
         else
            call ols_hatvalues(x, leverage, info)
            if (info /= SANDWICH_SUCCESS) then
               allocate(meat_matrix(0, 0))
               if (present(status)) status = info
               return
            end if
         end if
      end if

      allocate(meat_matrix(k, k))
      meat_matrix = 0.0_dp
      full_mask = 2**p - 1
      do mask = 1, full_mask
         nset = count_bits(mask, p)
         if (mod(nset, 2) == 1) then
            sign_value = 1
         else
            sign_value = -1
         end if

         if (use_multi0 .and. mask == full_mask) then
            allocate(labels(n))
            do i = 1, n
               labels(i) = i
            end do
            g = n
         else
            call group_subset(cluster, mask, labels, g)
         end if
         if (g <= 1 .and. use_cadjust .and. .not. (use_multi0 .and. mask == full_mask)) then
            deallocate(meat_matrix)
            allocate(meat_matrix(0, 0))
            if (present(status)) status = SANDWICH_INVALID_ARGUMENT
            return
         end if

         allocate(efi(n, k))
         efi = scores
         if (trim(kind) == 'hc2' .or. trim(kind) == 'hc3') then
            if (g == n) then
               do i = 1, n
                  if (leverage(i) >= 1.0_dp) then
                     deallocate(meat_matrix)
                     allocate(meat_matrix(0, 0))
                     if (present(status)) status = SANDWICH_INVALID_ARGUMENT
                     return
                  end if
                  if (trim(kind) == 'hc2') then
                     efi(i, :) = efi(i, :) / sqrt(1.0_dp - leverage(i))
                  else
                     efi(i, :) = efi(i, :) / (1.0_dp - leverage(i))
                  end if
               end do
            else
               do group = 1, g
                  group_size = count(labels == group)
                  allocate(indices(group_size), xg(group_size, k), resg(group_size))
                  position = 0
                  do i = 1, n
                     if (labels(i) == group) then
                        position = position + 1
                        indices(position) = i
                        xg(position, :) = x(i, :)
                        resg(position) = residuals(i)
                     end if
                  end do
                  hg = matmul(xg, matmul(invinfo, transpose(xg)))
                  hg = identity_matrix(group_size) - hg
                  if (trim(kind) == 'hc2') then
                     call symmetric_matrix_power(hg, -0.5_dp, adjustment_matrix, info)
                  else
                     call inverse_matrix(hg, adjustment_matrix, info)
                  end if
                  if (info /= SANDWICH_SUCCESS) then
                     deallocate(meat_matrix)
                     allocate(meat_matrix(0, 0))
                     if (present(status)) status = info
                     return
                  end if
                  adjusted_res = matmul(adjustment_matrix, resg)
                  do j = 1, group_size
                     efi(indices(j), :) = adjusted_res(j) * xg(j, :)
                  end do
                  deallocate(indices, xg, resg, hg, adjustment_matrix, adjusted_res)
               end do
            end if
            scale = sqrt(real(g - 1, dp) / real(g, dp))
            efi = scale * efi
         end if

         call aggregate_rows(efi, labels, g, aggregated)
         if (use_multi0 .and. mask == full_mask) then
            if (trim(kind) == 'hc1') then
               adj = real(n - k, dp) / real(n - 1, dp)
            else
               adj = 1.0_dp
            end if
         else if (use_cadjust) then
            adj = real(g, dp) / real(g - 1, dp)
         else
            adj = 1.0_dp
         end if
         meat_matrix = meat_matrix + real(sign_value, dp) * adj * &
            matmul(transpose(aggregated), aggregated) / real(n, dp)

         deallocate(labels, efi, aggregated)
      end do

      if (trim(kind) == 'hc1') then
         meat_matrix = real(n - 1, dp) / real(n - k, dp) * meat_matrix
      end if
      meat_matrix = 0.5_dp * (meat_matrix + transpose(meat_matrix))
      if (present(status)) status = SANDWICH_SUCCESS
   end subroutine meat_cluster

   subroutine vcov_cluster(scores, cluster, bread, covariance, status, type, cadjust, multi0, &
      fix, x, residuals, hat)
      real(dp), intent(in) :: scores(:, :), bread(:, :)
      integer, intent(in) :: cluster(:, :)
      real(dp), allocatable, intent(out) :: covariance(:, :)
      integer, intent(out), optional :: status
      character(len=*), intent(in), optional :: type
      logical, intent(in), optional :: cadjust, multi0, fix
      real(dp), intent(in), optional :: x(:, :), residuals(:), hat(:)
      real(dp), allocatable :: meat_matrix(:, :), fixed_covariance(:, :)
      character(len=8) :: selected_type
      logical :: use_cadjust, use_multi0, use_fix
      integer :: info

      selected_type = 'hc0'
      if (present(type)) selected_type = type
      use_cadjust = .true.
      if (present(cadjust)) use_cadjust = cadjust
      use_multi0 = .false.
      if (present(multi0)) use_multi0 = multi0
      use_fix = .false.
      if (present(fix)) use_fix = fix

      if (present(x) .and. present(residuals) .and. present(hat)) then
         call meat_cluster(scores, cluster, meat_matrix, info, selected_type, use_cadjust, &
            use_multi0, x, residuals, hat)
      else if (present(x) .and. present(residuals)) then
         call meat_cluster(scores, cluster, meat_matrix, info, selected_type, use_cadjust, &
            use_multi0, x, residuals)
      else
         call meat_cluster(scores, cluster, meat_matrix, info, selected_type, use_cadjust, use_multi0)
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
   end subroutine vcov_cluster

end module sandwich_cluster
