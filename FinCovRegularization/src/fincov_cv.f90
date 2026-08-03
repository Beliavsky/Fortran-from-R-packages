! SPDX-License-Identifier: GPL-2.0-only
module fincov_cv
   use fincov_kinds, only : dp
   use fincov_status, only : fincov_ok, fincov_invalid_input, fincov_size_mismatch, fincov_no_convergence
   use fincov_types, only : cv_result
   use fincov_utils, only : lowercase
   use fincov_rng, only : sample_without_replacement
   use fincov_linalg, only : sample_covariance, frobenius_norm_squared, spectral_norm_squared
   use fincov_regularization, only : banding, tapering, hard_thresholding, soft_thresholding, threshold_min
   implicit none
   private

   public :: banding_cv, tapering_cv, threshold_cv
contains
   function banding_cv(data, n_cv, norm, seed) result(result)
      real(dp), intent(in) :: data(:,:)
      integer, intent(in), optional :: n_cv, seed
      character(len=*), intent(in), optional :: norm
      type(cv_result) :: result
      real(dp), allocatable :: part1(:,:), part2(:,:), cov1(:,:), cov2(:,:), regularized(:,:), difference(:,:)
      character(len=1) :: norm_code
      integer :: n, p, n1, ncv_use, seed_use, repetition, k, local_status

      result%regularization = 'Banding'
      result%method = ''
      call cv_settings(data, n_cv, norm, seed, n, p, n1, ncv_use, seed_use, norm_code, local_status)
      call set_result_metadata(result, ncv_use, seed_use, norm_code)
      if (local_status /= fincov_ok) then
         result%status = local_status
         return
      end if
      allocate(result%parameter_grid(p), result%cv_error(p))
      result%parameter_grid = [(real(k,dp), k=0,p-1)]
      result%cv_error = 0.0_dp

      do repetition = 1, ncv_use
         call split_sample(data, n1, seed_use + repetition, part1, part2, local_status)
         if (local_status /= fincov_ok) exit
         call sample_covariance(part1, cov1, local_status)
         if (local_status /= fincov_ok) exit
         call sample_covariance(part2, cov2, local_status)
         if (local_status /= fincov_ok) exit
         do k = 0, p - 1
            regularized = banding(cov1, k, local_status)
            if (local_status /= fincov_ok) exit
            difference = regularized - cov2
            result%cv_error(k+1) = result%cv_error(k+1) + matrix_loss(difference, norm_code, local_status)
            if (local_status /= fincov_ok .and. local_status /= fincov_no_convergence) exit
         end do
         if (local_status /= fincov_ok .and. local_status /= fincov_no_convergence) exit
      end do
      if (local_status /= fincov_ok .and. local_status /= fincov_no_convergence) then
         result%status = local_status
         return
      end if
      result%cv_error = result%cv_error/real(ncv_use,dp)
      result%parameter_index = minloc(result%cv_error,dim=1)
      result%parameter_opt = result%parameter_grid(result%parameter_index)
      result%status = local_status
   end function banding_cv

   function tapering_cv(data, h, n_cv, norm, seed) result(result)
      real(dp), intent(in) :: data(:,:)
      real(dp), intent(in), optional :: h
      integer, intent(in), optional :: n_cv, seed
      character(len=*), intent(in), optional :: norm
      type(cv_result) :: result
      real(dp), allocatable :: part1(:,:), part2(:,:), cov1(:,:), cov2(:,:), regularized(:,:), difference(:,:)
      real(dp) :: ratio
      character(len=1) :: norm_code
      integer :: n, p, n1, ncv_use, seed_use, repetition, l_index, local_status

      ratio = 0.5_dp
      if (present(h)) ratio = h
      result%regularization = 'Tapering'
      result%method = ''
      result%h = ratio
      call cv_settings(data, n_cv, norm, seed, n, p, n1, ncv_use, seed_use, norm_code, local_status)
      call set_result_metadata(result, ncv_use, seed_use, norm_code)
      if (local_status /= fincov_ok .or. ratio <= 0.0_dp) then
         if (local_status == fincov_ok) local_status = fincov_invalid_input
         result%status = local_status
         return
      end if
      allocate(result%parameter_grid(p), result%cv_error(p))
      result%parameter_grid = [(real(l_index,dp), l_index=0,p-1)]
      result%cv_error = 0.0_dp

      do repetition = 1, ncv_use
         call split_sample(data, n1, seed_use + repetition, part1, part2, local_status)
         if (local_status /= fincov_ok) exit
         call sample_covariance(part1, cov1, local_status)
         if (local_status /= fincov_ok) exit
         call sample_covariance(part2, cov2, local_status)
         if (local_status /= fincov_ok) exit
         do l_index = 0, p - 1
            regularized = tapering(cov1, real(l_index,dp), ratio, local_status)
            if (local_status /= fincov_ok) exit
            difference = regularized - cov2
            result%cv_error(l_index+1) = result%cv_error(l_index+1) + matrix_loss(difference, norm_code, local_status)
            if (local_status /= fincov_ok .and. local_status /= fincov_no_convergence) exit
         end do
         if (local_status /= fincov_ok .and. local_status /= fincov_no_convergence) exit
      end do
      if (local_status /= fincov_ok .and. local_status /= fincov_no_convergence) then
         result%status = local_status
         return
      end if
      result%cv_error = result%cv_error/real(ncv_use,dp)
      result%parameter_index = minloc(result%cv_error,dim=1)
      result%parameter_opt = result%parameter_grid(result%parameter_index)
      result%status = local_status
   end function tapering_cv

   function threshold_cv(data, method, thresh_len, n_cv, norm, seed) result(result)
      real(dp), intent(in) :: data(:,:)
      character(len=*), intent(in), optional :: method
      integer, intent(in), optional :: thresh_len, n_cv, seed
      character(len=*), intent(in), optional :: norm
      type(cv_result) :: result
      real(dp), allocatable :: part1(:,:), part2(:,:), cov1(:,:), cov2(:,:), regularized(:,:), difference(:,:)
      real(dp), allocatable :: full_covariance(:,:)
      character(len=:), allocatable :: method_name
      character(len=1) :: norm_code
      real(dp) :: minimum_threshold, maximum_threshold
      integer :: n, p, n1, ncv_use, seed_use, grid_length, repetition, grid_index, i, j, local_status

      method_name = 'hard'
      if (present(method)) method_name = trim(lowercase(method))
      grid_length = 20
      if (present(thresh_len)) grid_length = thresh_len
      if (method_name == 'hard') then
         result%regularization = 'Hard Thresholding'
      else if (method_name == 'soft') then
         result%regularization = 'Soft Thresholding'
      else
         result%regularization = 'Thresholding'
      end if
      result%method = method_name
      call cv_settings(data, n_cv, norm, seed, n, p, n1, ncv_use, seed_use, norm_code, local_status)
      call set_result_metadata(result, ncv_use, seed_use, norm_code)
      if (local_status /= fincov_ok .or. grid_length < 1 .or. &
          (method_name /= 'hard' .and. method_name /= 'soft')) then
         if (local_status == fincov_ok) local_status = fincov_invalid_input
         result%status = local_status
         return
      end if

      call sample_covariance(data, full_covariance, local_status)
      if (local_status /= fincov_ok) then
         result%status = local_status
         return
      end if
      minimum_threshold = threshold_min(full_covariance, method_name, status=local_status)
      if (local_status /= fincov_ok .and. local_status /= fincov_no_convergence) then
         result%status = local_status
         return
      end if
      maximum_threshold = 0.0_dp
      do i = 1, p
         do j = 1, p
            if (i /= j) maximum_threshold = max(maximum_threshold, abs(full_covariance(i,j)))
         end do
      end do
      if (maximum_threshold <= tiny(1.0_dp)) maximum_threshold = minimum_threshold

      allocate(result%parameter_grid(grid_length), result%cv_error(grid_length))
      if (grid_length == 1) then
         result%parameter_grid(1) = minimum_threshold
      else
         do grid_index = 1, grid_length
            result%parameter_grid(grid_index) = minimum_threshold + &
               real(grid_index-1,dp)*(maximum_threshold-minimum_threshold)/real(grid_length-1,dp)
         end do
      end if
      result%cv_error = 0.0_dp

      do repetition = 1, ncv_use
         call split_sample(data, n1, seed_use + repetition, part1, part2, local_status)
         if (local_status /= fincov_ok) exit
         call sample_covariance(part1, cov1, local_status)
         if (local_status /= fincov_ok) exit
         call sample_covariance(part2, cov2, local_status)
         if (local_status /= fincov_ok) exit
         do grid_index = 1, grid_length
            if (method_name == 'hard') then
               regularized = hard_thresholding(cov1, result%parameter_grid(grid_index), local_status)
            else
               regularized = soft_thresholding(cov1, result%parameter_grid(grid_index), local_status)
            end if
            if (local_status /= fincov_ok) exit
            difference = regularized - cov2
            result%cv_error(grid_index) = result%cv_error(grid_index) + matrix_loss(difference, norm_code, local_status)
            if (local_status /= fincov_ok .and. local_status /= fincov_no_convergence) exit
         end do
         if (local_status /= fincov_ok .and. local_status /= fincov_no_convergence) exit
      end do
      if (local_status /= fincov_ok .and. local_status /= fincov_no_convergence) then
         result%status = local_status
         return
      end if
      result%cv_error = result%cv_error/real(ncv_use,dp)
      result%parameter_index = minloc(result%cv_error,dim=1)
      result%parameter_opt = result%parameter_grid(result%parameter_index)
      result%status = local_status
   end function threshold_cv

   subroutine cv_settings(data, n_cv, norm, seed, n, p, n1, ncv_use, seed_use, norm_code, status)
      real(dp), intent(in) :: data(:,:)
      integer, intent(in), optional :: n_cv, seed
      character(len=*), intent(in), optional :: norm
      integer, intent(out) :: n, p, n1, ncv_use, seed_use, status
      character(len=1), intent(out) :: norm_code
      character(len=:), allocatable :: norm_name
      integer :: n2

      n = size(data,1)
      p = size(data,2)
      ncv_use = 10
      seed_use = 142857
      if (present(n_cv)) ncv_use = n_cv
      if (present(seed)) seed_use = seed
      norm_name = 'f'
      if (present(norm)) norm_name = trim(lowercase(norm))
      if (len_trim(norm_name) > 0) then
         norm_code = norm_name(1:1)
      else
         norm_code = '?'
      end if

      if (n < 4 .or. p < 1 .or. ncv_use < 1 .or. (norm_code /= 'f' .and. norm_code /= 'o')) then
         n1 = 0
         status = fincov_invalid_input
         return
      end if
      n1 = ceiling(real(n,dp)*(1.0_dp - 1.0_dp/log(real(n,dp))))
      n2 = n - n1
      if (n1 < 2 .or. n2 < 2) then
         status = fincov_invalid_input
      else
         status = fincov_ok
      end if
   end subroutine cv_settings

   subroutine set_result_metadata(result, n_cv, seed, norm_code)
      type(cv_result), intent(inout) :: result
      integer, intent(in) :: n_cv, seed
      character(len=1), intent(in) :: norm_code
      result%n_cv = n_cv
      result%seed = seed
      result%norm = norm_code
   end subroutine set_result_metadata

   subroutine split_sample(data, n1, seed, part1, part2, status)
      real(dp), intent(in) :: data(:,:)
      integer, intent(in) :: n1, seed
      real(dp), allocatable, intent(out) :: part1(:,:), part2(:,:)
      integer, intent(out) :: status
      integer, allocatable :: index(:)
      logical, allocatable :: selected(:)
      integer :: n, p, i, j

      n = size(data,1)
      p = size(data,2)
      call sample_without_replacement(n, n1, seed, index, status)
      if (status /= fincov_ok) then
         allocate(part1(0,0),part2(0,0))
         return
      end if
      allocate(selected(n), part1(n1,p), part2(n-n1,p))
      selected = .false.
      do i = 1, n1
         selected(index(i)) = .true.
         part1(i,:) = data(index(i),:)
      end do
      j = 0
      do i = 1, n
         if (.not. selected(i)) then
            j = j + 1
            part2(j,:) = data(i,:)
         end if
      end do
      status = fincov_ok
   end subroutine split_sample

   function matrix_loss(difference, norm_code, status) result(loss)
      real(dp), intent(in) :: difference(:,:)
      character(len=1), intent(in) :: norm_code
      integer, intent(out) :: status
      real(dp) :: loss

      if (norm_code == 'f') then
         loss = frobenius_norm_squared(difference)
         status = fincov_ok
      else if (norm_code == 'o') then
         loss = spectral_norm_squared(difference, status)
      else
         loss = huge(1.0_dp)
         status = fincov_invalid_input
      end if
   end function matrix_loss
end module fincov_cv
