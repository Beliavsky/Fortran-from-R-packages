! SPDX-License-Identifier: GPL-2.0-only
! Modern Fortran translation of computational code from the Spillover R package.
module spillover_fevd
   use spillover_kinds, only : dp
   use spillover_status, only : spillover_success, spillover_invalid_argument, &
      spillover_not_positive_definite, set_status
   use spillover_linalg, only : cholesky_lower
   use spillover_var, only : var_model, ma_coefficients
   implicit none
   private

   public :: generalized_fevd
   public :: orthogonalized_fevd

contains

   subroutine generalized_fevd(model, horizon, normalized, fevd, info, message)
      type(var_model), intent(in) :: model
      integer, intent(in) :: horizon
      logical, intent(in) :: normalized
      real(dp), allocatable, intent(out) :: fevd(:, :, :)
      integer, intent(out), optional :: info
      character(len=*), intent(out), optional :: message

      real(dp), allocatable :: phi(:, :, :), numerator(:, :), denominator(:)
      real(dp), allocatable :: response(:, :), variance_term(:, :)
      real(dp) :: row_sum
      integer :: k, h, i, j, phi_info
      character(len=160) :: phi_message

      call set_status(info, message, spillover_success, 'success')
      k = model%k
      if (horizon < 1 .or. k < 1 .or. .not. allocated(model%sigma)) then
         allocate(fevd(0, 0, 0))
         call set_status(info, message, spillover_invalid_argument, &
            'invalid model or forecast horizon')
         return
      end if
      if (any([(model%sigma(i, i) <= 0.0_dp, i = 1, k)])) then
         allocate(fevd(0, 0, 0))
         call set_status(info, message, spillover_not_positive_definite, &
            'innovation covariance has a nonpositive diagonal')
         return
      end if

      call ma_coefficients(model, horizon, phi, phi_info, phi_message)
      if (phi_info /= spillover_success) then
         allocate(fevd(0, 0, 0))
         call set_status(info, message, phi_info, trim(phi_message))
         return
      end if

      allocate(fevd(k, k, horizon), numerator(k, k), denominator(k))
      allocate(response(k, k), variance_term(k, k))
      numerator = 0.0_dp
      denominator = 0.0_dp
      fevd = 0.0_dp

      do h = 1, horizon
         response = matmul(phi(:, :, h), model%sigma)
         variance_term = matmul(response, transpose(phi(:, :, h)))
         do i = 1, k
            denominator(i) = denominator(i) + variance_term(i, i)
            do j = 1, k
               numerator(i, j) = numerator(i, j) + &
                  response(i, j) * response(i, j) / model%sigma(j, j)
            end do
         end do

         do i = 1, k
            if (denominator(i) <= tiny(1.0_dp)) then
               call set_status(info, message, spillover_invalid_argument, &
                  'zero forecast-error variance encountered')
               return
            end if
            fevd(i, :, h) = numerator(i, :) / denominator(i)
            if (normalized) then
               row_sum = sum(fevd(i, :, h))
               if (row_sum <= tiny(1.0_dp)) then
                  call set_status(info, message, spillover_invalid_argument, &
                     'zero generalized FEVD row sum encountered')
                  return
               end if
               fevd(i, :, h) = fevd(i, :, h) / row_sum
            end if
         end do
      end do
   end subroutine generalized_fevd

   subroutine orthogonalized_fevd(model, horizon, fevd, perm, info, message)
      type(var_model), intent(in) :: model
      integer, intent(in) :: horizon
      real(dp), allocatable, intent(out) :: fevd(:, :, :)
      integer, intent(in), optional :: perm(:)
      integer, intent(out), optional :: info
      character(len=*), intent(out), optional :: message

      real(dp), allocatable :: phi(:, :, :), phi_p(:, :, :), sigma_p(:, :)
      real(dp), allocatable :: l(:, :), response(:, :), numerator(:, :), denominator(:)
      real(dp), allocatable :: fevd_p(:, :, :)
      integer, allocatable :: order(:)
      integer :: k, h, i, j, phi_info, chol_info
      character(len=160) :: local_message

      call set_status(info, message, spillover_success, 'success')
      k = model%k
      if (horizon < 1 .or. k < 1 .or. .not. allocated(model%sigma)) then
         allocate(fevd(0, 0, 0))
         call set_status(info, message, spillover_invalid_argument, &
            'invalid model or forecast horizon')
         return
      end if

      allocate(order(k))
      order = [(i, i = 1, k)]
      if (present(perm)) then
         if (.not. valid_permutation(perm, k)) then
            allocate(fevd(0, 0, 0))
            call set_status(info, message, spillover_invalid_argument, &
               'perm must contain each integer from 1 to k exactly once')
            return
         end if
         order = perm
      end if

      call ma_coefficients(model, horizon, phi, phi_info, local_message)
      if (phi_info /= spillover_success) then
         allocate(fevd(0, 0, 0))
         call set_status(info, message, phi_info, trim(local_message))
         return
      end if

      allocate(phi_p(k, k, horizon), sigma_p(k, k))
      do i = 1, k
         do j = 1, k
            sigma_p(i, j) = model%sigma(order(i), order(j))
            phi_p(i, j, :) = phi(order(i), order(j), :)
         end do
      end do

      call cholesky_lower(sigma_p, l, chol_info, local_message)
      if (chol_info /= spillover_success) then
         allocate(fevd(0, 0, 0))
         call set_status(info, message, chol_info, trim(local_message))
         return
      end if

      allocate(fevd_p(k, k, horizon), numerator(k, k), denominator(k), response(k, k))
      numerator = 0.0_dp
      denominator = 0.0_dp
      fevd_p = 0.0_dp
      do h = 1, horizon
         response = matmul(phi_p(:, :, h), l)
         do i = 1, k
            denominator(i) = denominator(i) + sum(response(i, :) ** 2)
            numerator(i, :) = numerator(i, :) + response(i, :) ** 2
            if (denominator(i) <= tiny(1.0_dp)) then
               call set_status(info, message, spillover_invalid_argument, &
                  'zero orthogonalized forecast-error variance encountered')
               return
            end if
            fevd_p(i, :, h) = numerator(i, :) / denominator(i)
         end do
      end do

      allocate(fevd(k, k, horizon))
      fevd = 0.0_dp
      do h = 1, horizon
         do i = 1, k
            do j = 1, k
               fevd(order(i), order(j), h) = fevd_p(i, j, h)
            end do
         end do
      end do
   end subroutine orthogonalized_fevd

   pure logical function valid_permutation(perm, k)
      integer, intent(in) :: perm(:)
      integer, intent(in) :: k
      logical :: seen(k)
      integer :: i

      valid_permutation = .false.
      if (size(perm) /= k) return
      seen = .false.
      do i = 1, k
         if (perm(i) < 1 .or. perm(i) > k) return
         if (seen(perm(i))) return
         seen(perm(i)) = .true.
      end do
      valid_permutation = .true.
   end function valid_permutation

end module spillover_fevd
