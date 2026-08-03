! SPDX-License-Identifier: GPL-2.0-only OR GPL-3.0-only
module sandwich_hc
   use sandwich_kinds, only : dp
   use sandwich_status, only : SANDWICH_SUCCESS, SANDWICH_INVALID_ARGUMENT, &
      SANDWICH_DIMENSION_MISMATCH, SANDWICH_NUMERICAL_FAILURE
   use sandwich_utils, only : lowercase
   use sandwich_regression, only : ols_hatvalues
   use sandwich_core, only : sandwich_covariance
   implicit none
   private

   public :: hc_weights, meat_hc, vcov_hc

contains

   subroutine hc_weights(residuals, hat, df_residual, type, omega, status)
      real(dp), intent(in) :: residuals(:), hat(:)
      integer, intent(in) :: df_residual
      character(len=*), intent(in) :: type
      real(dp), allocatable, intent(out) :: omega(:)
      integer, intent(out), optional :: status
      character(len=:), allocatable :: kind
      real(dp), allocatable :: delta(:)
      real(dp) :: p, nreal, leverage_ratio, cap, gamma1, gamma2
      integer :: n, i

      n = size(residuals)
      if (n <= 0 .or. size(hat) /= n .or. df_residual <= 0) then
         allocate(omega(0))
         if (present(status)) status = SANDWICH_DIMENSION_MISMATCH
         return
      end if

      kind = trim(lowercase(type))
      if (kind == 'hc') kind = 'hc0'
      allocate(omega(n))
      nreal = real(n, dp)

      select case (kind)
      case ('const')
         omega = sum(residuals**2) / real(df_residual, dp)
      case ('hc0')
         omega = residuals**2
      case ('hc1')
         omega = residuals**2 * nreal / real(df_residual, dp)
      case ('hc2')
         if (any(hat >= 1.0_dp)) then
            omega = huge(1.0_dp)
            if (present(status)) status = SANDWICH_NUMERICAL_FAILURE
            return
         end if
         omega = residuals**2 / (1.0_dp - hat)
      case ('hc3')
         if (any(hat >= 1.0_dp)) then
            omega = huge(1.0_dp)
            if (present(status)) status = SANDWICH_NUMERICAL_FAILURE
            return
         end if
         omega = residuals**2 / (1.0_dp - hat)**2
      case ('hc4', 'hc4m', 'hc5')
         if (any(hat >= 1.0_dp)) then
            omega = huge(1.0_dp)
            if (present(status)) status = SANDWICH_NUMERICAL_FAILURE
            return
         end if
         p = real(nint(sum(hat)), dp)
         if (p <= 0.0_dp) then
            omega = 0.0_dp
            if (present(status)) status = SANDWICH_INVALID_ARGUMENT
            return
         end if
         allocate(delta(n))
         select case (kind)
         case ('hc4')
            do i = 1, n
               leverage_ratio = nreal * hat(i) / p
               delta(i) = min(4.0_dp, leverage_ratio)
            end do
            omega = residuals**2 / (1.0_dp - hat)**delta
         case ('hc4m')
            gamma1 = 1.0_dp
            gamma2 = 1.5_dp
            do i = 1, n
               leverage_ratio = nreal * hat(i) / p
               delta(i) = min(gamma1, leverage_ratio) + min(gamma2, leverage_ratio)
            end do
            omega = residuals**2 / (1.0_dp - hat)**delta
         case ('hc5')
            cap = max(4.0_dp, nreal * 0.7_dp * maxval(hat) / p)
            do i = 1, n
               leverage_ratio = nreal * hat(i) / p
               delta(i) = min(leverage_ratio, cap)
            end do
            omega = residuals**2 / sqrt((1.0_dp - hat)**delta)
         end select
      case default
         omega = 0.0_dp
         if (present(status)) status = SANDWICH_INVALID_ARGUMENT
         return
      end select

      if (present(status)) status = SANDWICH_SUCCESS
   end subroutine hc_weights

   subroutine meat_hc(x, residuals, type, meat_matrix, status, hat, omega)
      real(dp), intent(in) :: x(:, :), residuals(:)
      character(len=*), intent(in) :: type
      real(dp), allocatable, intent(out) :: meat_matrix(:, :)
      integer, intent(out), optional :: status
      real(dp), intent(in), optional :: hat(:), omega(:)
      real(dp), allocatable :: leverage(:), variance_weights(:), weighted_x(:, :)
      integer :: n, k, i, info
      character(len=:), allocatable :: kind

      n = size(x, 1)
      k = size(x, 2)
      if (n <= k .or. k <= 0 .or. size(residuals) /= n) then
         allocate(meat_matrix(0, 0))
         if (present(status)) status = SANDWICH_DIMENSION_MISMATCH
         return
      end if

      if (present(omega)) then
         if (size(omega) /= n .or. any(omega < 0.0_dp)) then
            allocate(meat_matrix(0, 0))
            if (present(status)) status = SANDWICH_INVALID_ARGUMENT
            return
         end if
         allocate(variance_weights(n))
         variance_weights = omega
      else
         kind = trim(lowercase(type))
         if (kind == 'const' .or. kind == 'hc0' .or. kind == 'hc' .or. kind == 'hc1') then
            allocate(leverage(n))
            leverage = 0.0_dp
         else if (present(hat)) then
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
         call hc_weights(residuals, leverage, n - k, type, variance_weights, info)
         if (info /= SANDWICH_SUCCESS) then
            allocate(meat_matrix(0, 0))
            if (present(status)) status = info
            return
         end if
      end if

      allocate(weighted_x(n, k))
      do i = 1, n
         weighted_x(i, :) = sqrt(variance_weights(i)) * x(i, :)
      end do
      meat_matrix = matmul(transpose(weighted_x), weighted_x) / real(n, dp)
      meat_matrix = 0.5_dp * (meat_matrix + transpose(meat_matrix))
      if (present(status)) status = SANDWICH_SUCCESS
   end subroutine meat_hc

   subroutine vcov_hc(x, residuals, bread, type, covariance, status, hat, omega)
      real(dp), intent(in) :: x(:, :), residuals(:), bread(:, :)
      character(len=*), intent(in) :: type
      real(dp), allocatable, intent(out) :: covariance(:, :)
      integer, intent(out), optional :: status
      real(dp), intent(in), optional :: hat(:), omega(:)
      real(dp), allocatable :: meat_matrix(:, :)
      integer :: info

      if (present(hat) .and. present(omega)) then
         call meat_hc(x, residuals, type, meat_matrix, info, hat, omega)
      else if (present(hat)) then
         call meat_hc(x, residuals, type, meat_matrix, info, hat = hat)
      else if (present(omega)) then
         call meat_hc(x, residuals, type, meat_matrix, info, omega = omega)
      else
         call meat_hc(x, residuals, type, meat_matrix, info)
      end if
      if (info /= SANDWICH_SUCCESS) then
         allocate(covariance(0, 0))
         if (present(status)) status = info
         return
      end if
      call sandwich_covariance(bread, meat_matrix, size(x, 1), covariance, info)
      if (present(status)) status = info
   end subroutine vcov_hc

end module sandwich_hc
