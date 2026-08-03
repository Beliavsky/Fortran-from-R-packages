! SPDX-License-Identifier: GPL-2.0-only OR GPL-3.0-only
module sandwich_bootstrap
   use, intrinsic :: iso_fortran_env, only : int64
   use sandwich_kinds, only : dp
   use sandwich_status, only : SANDWICH_SUCCESS, SANDWICH_INVALID_ARGUMENT, &
      SANDWICH_DIMENSION_MISMATCH, SANDWICH_NUMERICAL_FAILURE
   use sandwich_utils, only : lowercase, count_bits, group_subset
   use sandwich_linalg, only : covariance_matrix, outer_product, project_psd
   use sandwich_regression, only : ols_model, fit_ols
   implicit none
   private

   real(dp), parameter :: pi = acos(-1.0_dp)

   public :: bootstrap_covariance, jackknife_covariance
   public :: vcov_bootstrap_ols, set_bootstrap_seed

contains

   subroutine set_bootstrap_seed(seed)
      integer, intent(in) :: seed
      integer, allocatable :: put(:)
      integer :: n, i
      integer(int64) :: value, modulus

      call random_seed(size = n)
      allocate(put(n))
      modulus = int(huge(1) - 1, int64)
      value = abs(int(seed, int64)) + 104729_int64
      do i = 1, n
         value = modulo(1664525_int64 * value + 1013904223_int64, modulus)
         put(i) = int(max(value, 1_int64))
      end do
      call random_seed(put = put)
   end subroutine set_bootstrap_seed

   subroutine bootstrap_covariance(replicates, covariance, status)
      real(dp), intent(in) :: replicates(:, :)
      real(dp), allocatable, intent(out) :: covariance(:, :)
      integer, intent(out), optional :: status
      integer :: info

      call covariance_matrix(replicates, covariance, info)
      if (present(status)) status = info
   end subroutine bootstrap_covariance

   subroutine jackknife_covariance(estimates, covariance, status, center_estimate)
      real(dp), intent(in) :: estimates(:, :)
      real(dp), allocatable, intent(out) :: covariance(:, :)
      integer, intent(out), optional :: status
      real(dp), intent(in), optional :: center_estimate(:)
      real(dp), allocatable :: center(:), deviations(:, :)
      integer :: r, k, i

      r = size(estimates, 1)
      k = size(estimates, 2)
      if (r < 2 .or. k <= 0) then
         allocate(covariance(0, 0))
         if (present(status)) status = SANDWICH_INVALID_ARGUMENT
         return
      end if
      allocate(center(k), deviations(r, k))
      if (present(center_estimate)) then
         if (size(center_estimate) /= k) then
            allocate(covariance(0, 0))
            if (present(status)) status = SANDWICH_DIMENSION_MISMATCH
            return
         end if
         center = center_estimate
      else
         center = sum(estimates, dim = 1) / real(r, dp)
      end if
      do i = 1, r
         deviations(i, :) = estimates(i, :) - center
      end do
      covariance = real(r - 1, dp) / real(r, dp) * &
         matmul(transpose(deviations), deviations)
      covariance = 0.5_dp * (covariance + transpose(covariance))
      if (present(status)) status = SANDWICH_SUCCESS
   end subroutine jackknife_covariance

   subroutine vcov_bootstrap_ols(x, y, cluster, covariance, status, replications, type, &
      weights, offset, center, fix, seed)
      real(dp), intent(in) :: x(:, :), y(:)
      integer, intent(in) :: cluster(:, :)
      real(dp), allocatable, intent(out) :: covariance(:, :)
      integer, intent(out), optional :: status
      integer, intent(in), optional :: replications, seed
      character(len=*), intent(in), optional :: type, center
      real(dp), intent(in), optional :: weights(:), offset(:)
      logical, intent(in), optional :: fix
      type(ols_model) :: base_model, boot_model
      real(dp), allocatable :: w(:), off(:), replicates(:, :), component(:, :)
      real(dp), allocatable :: xb(:, :), yb(:), wb(:), ob(:), yboot(:), fw(:), multipliers(:)
      real(dp), allocatable :: fixed_covariance(:, :)
      integer, allocatable :: labels(:), group_sizes(:), source_groups(:), row_indices(:)
      integer :: n, k, p, full_mask, mask, nset, sign_value, g, r_requested, r_actual
      integer :: i, j, group, source_group, target_group, total_rows, position, info
      logical :: use_fix, center_on_estimate
      character(len=16) :: selected_type, selected_center

      n = size(x, 1)
      k = size(x, 2)
      p = size(cluster, 2)
      if (n <= k .or. k <= 0 .or. size(y) /= n .or. size(cluster, 1) /= n .or. p <= 0) then
         allocate(covariance(0, 0))
         if (present(status)) status = SANDWICH_DIMENSION_MISMATCH
         return
      end if
      allocate(w(n), off(n))
      w = 1.0_dp
      off = 0.0_dp
      if (present(weights)) then
         if (size(weights) /= n .or. any(weights < 0.0_dp)) then
            allocate(covariance(0, 0))
            if (present(status)) status = SANDWICH_INVALID_ARGUMENT
            return
         end if
         w = weights
      end if
      if (present(offset)) then
         if (size(offset) /= n) then
            allocate(covariance(0, 0))
            if (present(status)) status = SANDWICH_DIMENSION_MISMATCH
            return
         end if
         off = offset
      end if

      selected_type = 'xy'
      if (present(type)) selected_type = trim(lowercase(type))
      if (trim(selected_type) == 'wild') selected_type = 'rademacher'
      if (trim(selected_type) == 'wild-rademacher') selected_type = 'rademacher'
      if (trim(selected_type) == 'wild-mammen') selected_type = 'mammen'
      if (trim(selected_type) == 'wild-norm') selected_type = 'norm'
      if (trim(selected_type) == 'wild-webb') selected_type = 'webb'
      if (trim(selected_type) /= 'xy' .and. trim(selected_type) /= 'jackknife' .and. &
          trim(selected_type) /= 'fractional' .and. trim(selected_type) /= 'residual' .and. &
          trim(selected_type) /= 'rademacher' .and. trim(selected_type) /= 'mammen' .and. &
          trim(selected_type) /= 'norm' .and. trim(selected_type) /= 'webb') then
         allocate(covariance(0, 0))
         if (present(status)) status = SANDWICH_INVALID_ARGUMENT
         return
      end if
      selected_center = 'mean'
      if (present(center)) selected_center = trim(lowercase(center))
      center_on_estimate = trim(selected_center) == 'estimate'
      use_fix = .false.
      if (present(fix)) use_fix = fix
      r_requested = 250
      if (present(replications)) r_requested = replications
      if (r_requested < 2 .and. trim(selected_type) /= 'jackknife') then
         allocate(covariance(0, 0))
         if (present(status)) status = SANDWICH_INVALID_ARGUMENT
         return
      end if
      if (present(seed)) call set_bootstrap_seed(seed)

      call fit_ols(x, y, base_model, info, w, off)
      if (info /= SANDWICH_SUCCESS) then
         allocate(covariance(0, 0))
         if (present(status)) status = info
         return
      end if

      allocate(covariance(k, k))
      covariance = 0.0_dp
      full_mask = 2**p - 1
      do mask = 1, full_mask
         nset = count_bits(mask, p)
         if (mod(nset, 2) == 1) then
            sign_value = 1
         else
            sign_value = -1
         end if
         call group_subset(cluster, mask, labels, g)
         if (g < 2) then
            deallocate(covariance)
            allocate(covariance(0, 0))
            if (present(status)) status = SANDWICH_INVALID_ARGUMENT
            return
         end if
         allocate(group_sizes(g))
         do group = 1, g
            group_sizes(group) = count(labels == group)
         end do

         if (trim(selected_type) == 'jackknife') then
            r_actual = g
         else
            r_actual = r_requested
         end if
         allocate(replicates(r_actual, k))

         do i = 1, r_actual
            select case (trim(selected_type))
            case ('xy')
               allocate(source_groups(g))
               do j = 1, g
                  source_groups(j) = random_integer(g)
               end do
               total_rows = 0
               do j = 1, g
                  total_rows = total_rows + group_sizes(source_groups(j))
               end do
               allocate(xb(total_rows, k), yb(total_rows), wb(total_rows), ob(total_rows))
               position = 0
               do j = 1, g
                  source_group = source_groups(j)
                  do group = 1, n
                     if (labels(group) == source_group) then
                        position = position + 1
                        xb(position, :) = x(group, :)
                        yb(position) = y(group)
                        wb(position) = w(group)
                        ob(position) = off(group)
                     end if
                  end do
               end do
               call fit_ols(xb, yb, boot_model, info, wb, ob)
               deallocate(source_groups, xb, yb, wb, ob)
            case ('jackknife')
               total_rows = n - group_sizes(i)
               allocate(xb(total_rows, k), yb(total_rows), wb(total_rows), ob(total_rows))
               position = 0
               do group = 1, n
                  if (labels(group) /= i) then
                     position = position + 1
                     xb(position, :) = x(group, :)
                     yb(position) = y(group)
                     wb(position) = w(group)
                     ob(position) = off(group)
                  end if
               end do
               call fit_ols(xb, yb, boot_model, info, wb, ob)
               deallocate(xb, yb, wb, ob)
            case ('fractional')
               allocate(fw(g), wb(n))
               do group = 1, g
                  fw(group) = random_exponential()
               end do
               fw = fw / (sum(fw) / real(g, dp))
               do group = 1, n
                  wb(group) = w(group) * fw(labels(group))
               end do
               call fit_ols(x, y, boot_model, info, wb, off)
               deallocate(fw, wb)
            case ('residual')
               if (any(group_sizes /= group_sizes(1))) then
                  info = SANDWICH_INVALID_ARGUMENT
               else
                  allocate(source_groups(g), yboot(n))
                  do target_group = 1, g
                     source_groups(target_group) = random_integer(g)
                  end do
                  yboot = base_model%fitted
                  do target_group = 1, g
                     source_group = source_groups(target_group)
                     call group_rows(labels, target_group, row_indices)
                     call assign_residual_block(labels, source_group, base_model%residuals, &
                        yboot, row_indices)
                     deallocate(row_indices)
                  end do
                  call fit_ols(x, yboot, boot_model, info, w, off)
                  deallocate(source_groups, yboot)
               end if
            case default
               allocate(multipliers(g), yboot(n))
               do group = 1, g
                  multipliers(group) = wild_multiplier(trim(selected_type))
               end do
               do group = 1, n
                  yboot(group) = base_model%fitted(group) + &
                     base_model%residuals(group) * multipliers(labels(group))
               end do
               call fit_ols(x, yboot, boot_model, info, w, off)
               deallocate(multipliers, yboot)
            end select
            if (info /= SANDWICH_SUCCESS) then
               deallocate(covariance)
               allocate(covariance(0, 0))
               if (present(status)) status = SANDWICH_NUMERICAL_FAILURE
               return
            end if
            replicates(i, :) = boot_model%coefficients
         end do

         if (trim(selected_type) == 'jackknife') then
            if (center_on_estimate) then
               call jackknife_covariance(replicates, component, info, base_model%coefficients)
            else
               call jackknife_covariance(replicates, component, info)
            end if
         else
            call bootstrap_covariance(replicates, component, info)
         end if
         if (info /= SANDWICH_SUCCESS) then
            deallocate(covariance)
            allocate(covariance(0, 0))
            if (present(status)) status = info
            return
         end if
         covariance = covariance + real(sign_value, dp) * component
         deallocate(labels, group_sizes, replicates, component)
      end do

      covariance = 0.5_dp * (covariance + transpose(covariance))
      if (use_fix) then
         call project_psd(covariance, fixed_covariance, info)
         if (info == SANDWICH_SUCCESS) call move_alloc(fixed_covariance, covariance)
      end if
      if (present(status)) status = SANDWICH_SUCCESS
   contains

      integer function random_integer(maximum) result(value)
         integer, intent(in) :: maximum
         real(dp) :: u
         call random_number(u)
         value = min(maximum, int(u * real(maximum, dp)) + 1)
      end function random_integer

      real(dp) function random_exponential() result(value)
         real(dp) :: u
         call random_number(u)
         u = max(u, tiny(1.0_dp))
         value = -log(u)
      end function random_exponential

      real(dp) function random_normal() result(value)
         real(dp) :: u1, u2
         call random_number(u1)
         call random_number(u2)
         u1 = max(u1, tiny(1.0_dp))
         value = sqrt(-2.0_dp * log(u1)) * cos(2.0_dp * pi * u2)
      end function random_normal

      real(dp) function wild_multiplier(kind) result(value)
         character(len=*), intent(in) :: kind
         real(dp) :: u, sqrt5
         integer :: choice

         call random_number(u)
         select case (kind)
         case ('rademacher')
            if (u < 0.5_dp) then
               value = -1.0_dp
            else
               value = 1.0_dp
            end if
         case ('mammen')
            sqrt5 = sqrt(5.0_dp)
            if (u < (sqrt5 + 1.0_dp) / (2.0_dp * sqrt5)) then
               value = -(sqrt5 - 1.0_dp) / 2.0_dp
            else
               value = (sqrt5 + 1.0_dp) / 2.0_dp
            end if
         case ('norm')
            value = random_normal()
         case ('webb')
            choice = min(6, int(6.0_dp * u) + 1)
            select case (choice)
            case (1)
               value = -sqrt(1.5_dp)
            case (2)
               value = -1.0_dp
            case (3)
               value = -sqrt(0.5_dp)
            case (4)
               value = sqrt(0.5_dp)
            case (5)
               value = 1.0_dp
            case default
               value = sqrt(1.5_dp)
            end select
         case default
            value = 0.0_dp
         end select
      end function wild_multiplier

      subroutine group_rows(group_labels, selected_group, rows)
         integer, intent(in) :: group_labels(:), selected_group
         integer, allocatable, intent(out) :: rows(:)
         integer :: row, pos

         allocate(rows(count(group_labels == selected_group)))
         pos = 0
         do row = 1, size(group_labels)
            if (group_labels(row) == selected_group) then
               pos = pos + 1
               rows(pos) = row
            end if
         end do
      end subroutine group_rows

      subroutine assign_residual_block(group_labels, source, residuals, response, target_rows)
         integer, intent(in) :: group_labels(:), source, target_rows(:)
         real(dp), intent(in) :: residuals(:)
         real(dp), intent(inout) :: response(:)
         integer, allocatable :: source_rows(:)
         integer :: row

         call group_rows(group_labels, source, source_rows)
         do row = 1, size(target_rows)
            response(target_rows(row)) = response(target_rows(row)) + residuals(source_rows(row))
         end do
      end subroutine assign_residual_block

   end subroutine vcov_bootstrap_ols

end module sandwich_bootstrap
