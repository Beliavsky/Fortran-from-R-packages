! SPDX-License-Identifier: MIT
module r_time_series
   use r_descriptive, only : r_mean
   use r_kinds, only : dp
   use r_missing, only : r_is_finite
   use r_status, only : r_invalid_input, r_ok, r_singular
   implicit none
   private

   public :: r_autocovariance, r_autocorrelation
   public :: r_cross_covariance, r_cross_correlation

contains

   pure subroutine r_autocovariance(x, values, lag_max, demean, status)
      real(dp), intent(in) :: x(:)
      real(dp), allocatable, intent(out) :: values(:)
      integer, intent(in), optional :: lag_max
      logical, intent(in), optional :: demean
      integer, intent(out), optional :: status
      real(dp) :: center
      integer :: k, m, n
      logical :: center_data

      if (present(status)) status = r_ok
      n = size(x)
      if (n < 1 .or. any(.not. r_is_finite(x))) then
         if (present(status)) status = r_invalid_input
         return
      end if

      if (n == 1) then
         m = 0
      else
         m = min(n - 1, max(1, int(10.0_dp * log10(real(n, dp)))))
      end if
      if (present(lag_max)) m = lag_max
      if (m < 0 .or. m >= n) then
         if (present(status)) status = r_invalid_input
         return
      end if

      center_data = .true.
      if (present(demean)) center_data = demean
      center = 0.0_dp
      if (center_data) center = r_mean(x)

      allocate(values(0:m))
      do k = 0, m
         values(k) = dot_product(x(1:n-k) - center, x(1+k:n) - center) / real(n, dp)
      end do
   end subroutine r_autocovariance

   pure subroutine r_autocorrelation(x, values, lag_max, demean, status)
      real(dp), intent(in) :: x(:)
      real(dp), allocatable, intent(out) :: values(:)
      integer, intent(in), optional :: lag_max
      logical, intent(in), optional :: demean
      integer, intent(out), optional :: status
      real(dp), allocatable :: covariance(:)
      integer :: local_status

      call r_autocovariance(x, covariance, lag_max, demean, local_status)
      if (local_status /= r_ok) then
         if (present(status)) status = local_status
         return
      end if
      if (covariance(0) <= 0.0_dp) then
         if (present(status)) status = r_singular
         return
      end if
      allocate(values(lbound(covariance, 1):ubound(covariance, 1)))
      values = covariance / covariance(0)
      if (present(status)) status = r_ok
   end subroutine r_autocorrelation

   pure subroutine r_cross_covariance(x, values, lag_max, demean, status)
      real(dp), intent(in) :: x(:,:)
      real(dp), allocatable, intent(out) :: values(:,:,:)
      integer, intent(in), optional :: lag_max
      logical, intent(in), optional :: demean
      integer, intent(out), optional :: status
      real(dp), allocatable :: centered(:,:), means(:)
      integer :: i, j, k, m, n, p
      logical :: center_data

      if (present(status)) status = r_ok
      n = size(x, 1)
      p = size(x, 2)
      if (n < 1 .or. p < 1 .or. any(.not. r_is_finite(x))) then
         if (present(status)) status = r_invalid_input
         return
      end if

      if (n == 1) then
         m = 0
      else
         m = min(n - 1, max(1, int(10.0_dp * log10(real(n, dp)))))
      end if
      if (present(lag_max)) m = lag_max
      if (m < 0 .or. m >= n) then
         if (present(status)) status = r_invalid_input
         return
      end if

      center_data = .true.
      if (present(demean)) center_data = demean
      allocate(centered(n, p), means(p))
      means = 0.0_dp
      if (center_data) then
         do j = 1, p
            means(j) = r_mean(x(:, j))
         end do
      end if
      do j = 1, p
         centered(:, j) = x(:, j) - means(j)
      end do

      allocate(values(0:m, p, p))
      do k = 0, m
         do i = 1, p
            do j = 1, p
               values(k, i, j) = dot_product(centered(1:n-k, i), centered(1+k:n, j)) / real(n, dp)
            end do
         end do
      end do
   end subroutine r_cross_covariance

   pure subroutine r_cross_correlation(x, values, lag_max, demean, status)
      real(dp), intent(in) :: x(:,:)
      real(dp), allocatable, intent(out) :: values(:,:,:)
      integer, intent(in), optional :: lag_max
      logical, intent(in), optional :: demean
      integer, intent(out), optional :: status
      real(dp), allocatable :: covariance(:,:,:), scale(:)
      integer :: i, j, local_status, p

      call r_cross_covariance(x, covariance, lag_max, demean, local_status)
      if (local_status /= r_ok) then
         if (present(status)) status = local_status
         return
      end if

      p = size(x, 2)
      allocate(scale(p))
      do i = 1, p
         if (covariance(0, i, i) <= 0.0_dp) then
            if (present(status)) status = r_singular
            return
         end if
         scale(i) = sqrt(covariance(0, i, i))
      end do
      allocate(values(lbound(covariance, 1):ubound(covariance, 1), p, p))
      do i = 1, p
         do j = 1, p
            values(:, i, j) = covariance(:, i, j) / (scale(i) * scale(j))
         end do
      end do
      if (present(status)) status = r_ok
   end subroutine r_cross_correlation

end module r_time_series
