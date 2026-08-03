! SPDX-License-Identifier: GPL-2.0-only OR GPL-3.0-only
module sandwich_kernels
   use sandwich_kinds, only : dp
   use sandwich_status, only : SANDWICH_SUCCESS, SANDWICH_INVALID_ARGUMENT
   use sandwich_utils, only : lowercase
   implicit none
   private

   real(dp), parameter :: pi = acos(-1.0_dp)

   public :: kernel_weight, kernel_weights

contains

   real(dp) function kernel_weight(x, kernel, normalize, status) result(weight)
      real(dp), intent(in) :: x
      character(len=*), intent(in) :: kernel
      logical, intent(in), optional :: normalize
      integer, intent(out), optional :: status
      character(len=:), allocatable :: kind
      real(dp) :: ca, z, y, ay
      logical :: use_normalize

      use_normalize = .false.
      if (present(normalize)) use_normalize = normalize
      kind = trim(lowercase(kernel))
      select case (kind)
      case ('truncated', 'uniform')
         ca = 1.0_dp
         if (use_normalize) ca = 2.0_dp
         z = ca * x
         if (abs(z) > 1.0_dp) then
            weight = 0.0_dp
         else
            weight = 1.0_dp
         end if
      case ('bartlett', 'triangular')
         ca = 1.0_dp
         if (use_normalize) ca = 2.0_dp / 3.0_dp
         z = ca * x
         if (abs(z) > 1.0_dp) then
            weight = 0.0_dp
         else
            weight = 1.0_dp - abs(z)
         end if
      case ('parzen')
         ca = 1.0_dp
         if (use_normalize) ca = 0.539285_dp
         z = ca * x
         ay = abs(z)
         if (ay > 1.0_dp) then
            weight = 0.0_dp
         else if (ay < 0.5_dp) then
            weight = 1.0_dp - 6.0_dp * z * z + 6.0_dp * ay**3
         else
            weight = 2.0_dp * (1.0_dp - ay)**3
         end if
      case ('tukey-hanning', 'tukey hanning', 'tukey_hanning')
         ca = 1.0_dp
         if (use_normalize) ca = 0.75_dp
         z = ca * x
         if (abs(z) > 1.0_dp) then
            weight = 0.0_dp
         else
            weight = 0.5_dp * (1.0_dp + cos(pi * z))
         end if
      case ('quadratic spectral', 'quadratic-spectral', 'quadratic_spectral', 'qs')
         z = x
         y = 6.0_dp * pi * z / 5.0_dp
         if (abs(y) < 1.0e-4_dp) then
            weight = 1.0_dp - y * y / 10.0_dp + y**4 / 280.0_dp
         else
            weight = 3.0_dp / (y * y) * (sin(y) / y - cos(y))
         end if
      case default
         weight = 0.0_dp
         if (present(status)) status = SANDWICH_INVALID_ARGUMENT
         return
      end select
      if (present(status)) status = SANDWICH_SUCCESS
   end function kernel_weight

   subroutine kernel_weights(x, kernel, weights, status, normalize)
      real(dp), intent(in) :: x(:)
      character(len=*), intent(in) :: kernel
      real(dp), allocatable, intent(out) :: weights(:)
      integer, intent(out), optional :: status
      logical, intent(in), optional :: normalize
      integer :: i, info
      logical :: use_normalize

      use_normalize = .false.
      if (present(normalize)) use_normalize = normalize
      allocate(weights(size(x)))
      do i = 1, size(x)
         weights(i) = kernel_weight(x(i), kernel, use_normalize, info)
         if (info /= SANDWICH_SUCCESS) then
            weights = 0.0_dp
            if (present(status)) status = info
            return
         end if
      end do
      if (present(status)) status = SANDWICH_SUCCESS
   end subroutine kernel_weights

end module sandwich_kernels
