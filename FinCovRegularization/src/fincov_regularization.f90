! SPDX-License-Identifier: GPL-2.0-only
module fincov_regularization
   use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
   use fincov_kinds, only : dp
   use fincov_status, only : fincov_ok, fincov_invalid_input, fincov_size_mismatch, fincov_no_convergence
   use fincov_utils, only : lowercase
   use fincov_linalg, only : sample_covariance, minimum_eigenvalue
   implicit none
   private

   public :: banding, tapering, hard_thresholding, soft_thresholding
   public :: ind_cov, threshold_min
contains
   function banding(sigma, k, status) result(regularized)
      real(dp), intent(in) :: sigma(:,:)
      integer, intent(in), optional :: k
      integer, intent(out), optional :: status
      real(dp) :: regularized(size(sigma,1),size(sigma,2))
      integer :: bandwidth, i, j, p

      regularized = sigma
      p = size(sigma,1)
      bandwidth = 0
      if (present(k)) bandwidth = k
      if (p < 1 .or. size(sigma,2) /= p .or. bandwidth < 0) then
         regularized = ieee_value(0.0_dp, ieee_quiet_nan)
         if (present(status)) status = fincov_invalid_input
         return
      end if
      do i = 1, p
         do j = 1, p
            if (abs(i-j) > bandwidth) regularized(i,j) = 0.0_dp
         end do
      end do
      regularized = 0.5_dp * (regularized + transpose(regularized))
      if (present(status)) status = fincov_ok
   end function banding

   function tapering(sigma, l, h, status) result(regularized)
      real(dp), intent(in) :: sigma(:,:)
      real(dp), intent(in) :: l
      real(dp), intent(in), optional :: h
      integer, intent(out), optional :: status
      real(dp) :: regularized(size(sigma,1),size(sigma,2))
      real(dp) :: ratio, distance, weight
      integer :: i, j, p

      p = size(sigma,1)
      ratio = 0.5_dp
      if (present(h)) ratio = h
      if (p < 1 .or. size(sigma,2) /= p .or. l < 0.0_dp .or. ratio <= 0.0_dp) then
         regularized = ieee_value(0.0_dp, ieee_quiet_nan)
         if (present(status)) status = fincov_invalid_input
         return
      end if

      regularized = 0.0_dp
      do i = 1, p
         do j = 1, p
            distance = real(abs(i-j),dp)
            if (distance <= l*ratio) then
               weight = 1.0_dp
            else if (distance < l) then
               weight = 2.0_dp - distance/(l*ratio)
            else
               weight = 0.0_dp
            end if
            regularized(i,j) = sigma(i,j)*weight
         end do
      end do
      regularized = 0.5_dp * (regularized + transpose(regularized))
      if (present(status)) status = fincov_ok
   end function tapering

   function hard_thresholding(sigma, threshold, status) result(regularized)
      real(dp), intent(in) :: sigma(:,:)
      real(dp), intent(in), optional :: threshold
      integer, intent(out), optional :: status
      real(dp) :: regularized(size(sigma,1),size(sigma,2))
      real(dp) :: cutoff
      integer :: i, j, p

      p = size(sigma,1)
      cutoff = 0.5_dp
      if (present(threshold)) cutoff = threshold
      if (p < 1 .or. size(sigma,2) /= p .or. cutoff < 0.0_dp) then
         regularized = ieee_value(0.0_dp, ieee_quiet_nan)
         if (present(status)) status = fincov_invalid_input
         return
      end if
      regularized = sigma
      do i = 1, p
         do j = 1, p
            if (i /= j .and. abs(regularized(i,j)) < cutoff) regularized(i,j) = 0.0_dp
         end do
      end do
      regularized = 0.5_dp * (regularized + transpose(regularized))
      if (present(status)) status = fincov_ok
   end function hard_thresholding

   function soft_thresholding(sigma, threshold, status) result(regularized)
      real(dp), intent(in) :: sigma(:,:)
      real(dp), intent(in), optional :: threshold
      integer, intent(out), optional :: status
      real(dp) :: regularized(size(sigma,1),size(sigma,2))
      real(dp) :: cutoff
      integer :: i, j, p

      p = size(sigma,1)
      cutoff = 0.5_dp
      if (present(threshold)) cutoff = threshold
      if (p < 1 .or. size(sigma,2) /= p .or. cutoff < 0.0_dp) then
         regularized = ieee_value(0.0_dp, ieee_quiet_nan)
         if (present(status)) status = fincov_invalid_input
         return
      end if
      regularized = 0.0_dp
      do i = 1, p
         regularized(i,i) = sigma(i,i)
         do j = 1, p
            if (i /= j) regularized(i,j) = sign(max(abs(sigma(i,j)) - cutoff, 0.0_dp), sigma(i,j))
         end do
      end do
      regularized = 0.5_dp * (regularized + transpose(regularized))
      if (present(status)) status = fincov_ok
   end function soft_thresholding

   function ind_cov(sigma, status) result(regularized)
      real(dp), intent(in) :: sigma(:,:)
      integer, intent(out), optional :: status
      real(dp) :: regularized(size(sigma,1),size(sigma,2))
      real(dp), allocatable :: covariance(:,:)
      integer :: i, p, local_status

      p = size(sigma,1)
      regularized = 0.0_dp
      if (p < 2 .or. size(sigma,2) /= p) then
         regularized = ieee_value(0.0_dp, ieee_quiet_nan)
         if (present(status)) status = fincov_size_mismatch
         return
      end if
      call sample_covariance(sigma, covariance, local_status)
      if (local_status == fincov_ok) then
         do i = 1, p
            regularized(i,i) = covariance(i,i)
         end do
      else
         regularized = ieee_value(0.0_dp, ieee_quiet_nan)
      end if
      if (present(status)) status = local_status
   end function ind_cov

   function threshold_min(sigma, method, tolerance, status) result(value)
      real(dp), intent(in) :: sigma(:,:)
      character(len=*), intent(in), optional :: method
      real(dp), intent(in), optional :: tolerance
      integer, intent(out), optional :: status
      real(dp) :: value
      character(len=:), allocatable :: method_name
      real(dp) :: lower, upper, middle, f_lower, f_upper, f_middle, tol
      integer :: local_status, iteration

      method_name = 'hard'
      if (present(method)) method_name = trim(lowercase(method))
      if (method_name /= 'hard' .and. method_name /= 'soft') then
         value = ieee_value(0.0_dp, ieee_quiet_nan)
         if (present(status)) status = fincov_invalid_input
         return
      end if
      if (size(sigma,1) < 1 .or. size(sigma,2) /= size(sigma,1)) then
         value = ieee_value(0.0_dp, ieee_quiet_nan)
         if (present(status)) status = fincov_size_mismatch
         return
      end if

      lower = 0.0_dp
      upper = maxval(sigma)
      tol = sqrt(epsilon(1.0_dp))
      if (present(tolerance)) tol = max(tolerance, epsilon(1.0_dp))
      f_lower = threshold_eigenvalue(sigma, lower, method_name, local_status)
      if (local_status /= fincov_ok .and. local_status /= fincov_no_convergence) then
         value = ieee_value(0.0_dp, ieee_quiet_nan)
         if (present(status)) status = local_status
         return
      end if
      if (upper <= 0.0_dp) then
         value = 0.0_dp
         if (present(status)) status = fincov_ok
         return
      end if
      f_upper = threshold_eigenvalue(sigma, upper, method_name, local_status)
      if (local_status /= fincov_ok .and. local_status /= fincov_no_convergence) then
         value = ieee_value(0.0_dp, ieee_quiet_nan)
         if (present(status)) status = local_status
         return
      end if

      if (f_lower*f_upper >= 0.0_dp) then
         value = 0.0_dp
         if (present(status)) status = fincov_ok
         return
      end if

      do iteration = 1, 200
         middle = 0.5_dp*(lower + upper)
         f_middle = threshold_eigenvalue(sigma, middle, method_name, local_status)
         if (abs(f_middle) <= tol .or. abs(upper-lower) <= tol*max(1.0_dp,abs(middle))) exit
         if (f_lower*f_middle <= 0.0_dp) then
            upper = middle
            f_upper = f_middle
         else
            lower = middle
            f_lower = f_middle
         end if
      end do
      value = max(0.0_dp, 0.5_dp*(lower + upper))
      if (present(status)) then
         if (iteration > 200) then
            status = fincov_no_convergence
         else
            status = fincov_ok
         end if
      end if
   end function threshold_min

   function threshold_eigenvalue(sigma, threshold, method, status) result(value)
      real(dp), intent(in) :: sigma(:,:), threshold
      character(len=*), intent(in) :: method
      integer, intent(out) :: status
      real(dp) :: value
      real(dp) :: regularized(size(sigma,1),size(sigma,2))

      value = ieee_value(0.0_dp, ieee_quiet_nan)
      if (method == 'hard') then
         regularized = hard_thresholding(sigma, threshold, status)
      else
         regularized = soft_thresholding(sigma, threshold, status)
      end if
      if (status == fincov_ok) value = minimum_eigenvalue(regularized, status)
   end function threshold_eigenvalue
end module fincov_regularization
