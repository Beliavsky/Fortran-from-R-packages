! SPDX-License-Identifier: GPL-2.0-only
module fincov_factor_models
   use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
   use fincov_kinds, only : dp
   use fincov_status, only : fincov_ok, fincov_invalid_input, fincov_size_mismatch, &
      fincov_singular_matrix, fincov_no_convergence
   use fincov_utils, only : lowercase
   use fincov_linalg, only : sample_covariance, column_variances, solve_linear_system, symmetric_eigen_jacobi, add_diagonal
   implicit none
   private

   public :: macro_factor_cov, fundamental_factor_cov, stat_factor_cov

   interface macro_factor_cov
      module procedure macro_factor_cov_vector
      module procedure macro_factor_cov_matrix
   end interface macro_factor_cov
contains
   function macro_factor_cov_vector(assets, factor, status) result(covariance_estimate)
      real(dp), intent(in) :: assets(:,:), factor(:)
      integer, intent(out), optional :: status
      real(dp) :: covariance_estimate(size(assets,2),size(assets,2))
      real(dp), allocatable :: factor_matrix(:,:)
      integer :: local_status

      if (size(factor) /= size(assets,1)) then
         covariance_estimate = ieee_value(0.0_dp, ieee_quiet_nan)
         if (present(status)) status = fincov_size_mismatch
         return
      end if
      allocate(factor_matrix(size(factor),1))
      factor_matrix(:,1) = factor
      covariance_estimate = macro_factor_cov_matrix(assets, factor_matrix, local_status)
      if (present(status)) status = local_status
   end function macro_factor_cov_vector

   function macro_factor_cov_matrix(assets, factors, status) result(covariance_estimate)
      real(dp), intent(in) :: assets(:,:), factors(:,:)
      integer, intent(out), optional :: status
      real(dp) :: covariance_estimate(size(assets,2),size(assets,2))
      real(dp), allocatable :: design(:,:), xtx(:,:), xty(:,:), beta(:,:), residuals(:,:)
      real(dp), allocatable :: factor_covariance(:,:), slopes(:,:), residual_variances(:)
      integer :: n, p, q, df, j, local_status

      n = size(assets,1)
      p = size(assets,2)
      q = size(factors,2)
      covariance_estimate = ieee_value(0.0_dp, ieee_quiet_nan)
      if (n < 2 .or. p < 1 .or. q < 1 .or. size(factors,1) /= n) then
         if (present(status)) status = fincov_size_mismatch
         return
      end if
      df = n - q - 1
      if (df < 1) then
         if (present(status)) status = fincov_invalid_input
         return
      end if

      allocate(design(n,q+1), xtx(q+1,q+1), xty(q+1,p), residuals(n,p), &
         slopes(q,p), residual_variances(p))
      design(:,1) = 1.0_dp
      design(:,2:q+1) = factors
      xtx = matmul(transpose(design), design)
      xty = matmul(transpose(design), assets)
      call solve_linear_system(xtx, xty, beta, local_status)
      if (local_status /= fincov_ok) then
         if (present(status)) status = local_status
         return
      end if
      residuals = assets - matmul(design, beta)
      do j = 1, p
         residual_variances(j) = dot_product(residuals(:,j), residuals(:,j)) / real(df,dp)
      end do
      call sample_covariance(factors, factor_covariance, local_status)
      if (local_status /= fincov_ok) then
         if (present(status)) status = local_status
         return
      end if
      slopes = beta(2:q+1,:)
      covariance_estimate = matmul(transpose(slopes), matmul(factor_covariance, slopes))
      call add_diagonal(covariance_estimate, residual_variances)
      covariance_estimate = 0.5_dp*(covariance_estimate + transpose(covariance_estimate))
      if (present(status)) status = fincov_ok
   end function macro_factor_cov_matrix

   function fundamental_factor_cov(assets, exposure, method, status) result(covariance_estimate)
      real(dp), intent(in) :: assets(:,:), exposure(:,:)
      character(len=*), intent(in), optional :: method
      integer, intent(out), optional :: status
      real(dp) :: covariance_estimate(size(assets,2),size(assets,2))
      character(len=:), allocatable :: method_name
      real(dp), allocatable :: xtx(:,:), rhs(:,:), factor_ols(:,:), residuals(:,:)
      real(dp), allocatable :: residual_variances(:), factor_covariance(:,:)
      real(dp), allocatable :: dinv(:,:), middle(:,:), h_rhs(:,:), hmat(:,:), factor_wls(:,:)
      integer :: n, p, q, i, local_status

      n = size(assets,1)
      p = size(assets,2)
      q = size(exposure,2)
      covariance_estimate = ieee_value(0.0_dp, ieee_quiet_nan)
      method_name = 'wls'
      if (present(method)) method_name = trim(lowercase(method))
      if (method_name /= 'ols' .and. method_name /= 'wls') then
         if (present(status)) status = fincov_invalid_input
         return
      end if
      if (n < 2 .or. p < 1 .or. q < 1 .or. size(exposure,1) /= p) then
         if (present(status)) status = fincov_size_mismatch
         return
      end if

      allocate(xtx(q,q), rhs(q,n), residuals(n,p))
      xtx = matmul(transpose(exposure), exposure)
      rhs = matmul(transpose(exposure), transpose(assets))
      call solve_linear_system(xtx, rhs, factor_ols, local_status)
      if (local_status /= fincov_ok) then
         if (present(status)) status = local_status
         return
      end if
      residuals = assets - transpose(matmul(exposure, factor_ols))
      call column_variances(residuals, residual_variances, local_status)
      if (local_status /= fincov_ok) then
         if (present(status)) status = local_status
         return
      end if

      if (method_name == 'ols') then
         call sample_covariance(transpose(factor_ols), factor_covariance, local_status)
         if (local_status /= fincov_ok) then
            if (present(status)) status = local_status
            return
         end if
      else
         if (any(residual_variances <= 100.0_dp*tiny(1.0_dp))) then
            if (present(status)) status = fincov_singular_matrix
            return
         end if
         allocate(dinv(p,p), middle(q,q), h_rhs(q,p))
         dinv = 0.0_dp
         do i = 1, p
            dinv(i,i) = 1.0_dp/residual_variances(i)
         end do
         middle = matmul(transpose(exposure), matmul(dinv, exposure))
         h_rhs = matmul(transpose(exposure), dinv)
         call solve_linear_system(middle, h_rhs, hmat, local_status)
         if (local_status /= fincov_ok) then
            if (present(status)) status = local_status
            return
         end if
         factor_wls = matmul(assets, transpose(hmat))
         residuals = assets - matmul(factor_wls, transpose(exposure))
         call column_variances(residuals, residual_variances, local_status)
         if (local_status /= fincov_ok) then
            if (present(status)) status = local_status
            return
         end if
         call sample_covariance(factor_wls, factor_covariance, local_status)
         if (local_status /= fincov_ok) then
            if (present(status)) status = local_status
            return
         end if
      end if

      covariance_estimate = matmul(exposure, matmul(factor_covariance, transpose(exposure)))
      call add_diagonal(covariance_estimate, residual_variances)
      covariance_estimate = 0.5_dp*(covariance_estimate + transpose(covariance_estimate))
      if (present(status)) status = fincov_ok
   end function fundamental_factor_cov

   function stat_factor_cov(assets, k, status, selected_k) result(covariance_estimate)
      real(dp), intent(in) :: assets(:,:)
      integer, intent(in), optional :: k
      integer, intent(out), optional :: status
      integer, intent(out), optional :: selected_k
      real(dp) :: covariance_estimate(size(assets,2),size(assets,2))
      real(dp), allocatable :: covariance_sample(:,:), eigenvalues(:), eigenvectors(:,:)
      real(dp), allocatable :: common_covariance(:,:), specific_variances(:)
      real(dp) :: singular_value
      integer :: n, p, k_requested, k_use, k_max, i, local_status

      n = size(assets,1)
      p = size(assets,2)
      covariance_estimate = ieee_value(0.0_dp, ieee_quiet_nan)
      k_requested = 0
      if (present(k)) k_requested = k
      if (n < 2 .or. p < 1 .or. k_requested < 0) then
         if (present(status)) status = fincov_invalid_input
         if (present(selected_k)) selected_k = 0
         return
      end if
      k_max = min(n,p)
      if (k_requested > k_max) then
         if (present(status)) status = fincov_invalid_input
         if (present(selected_k)) selected_k = 0
         return
      end if

      call sample_covariance(assets, covariance_sample, local_status)
      if (local_status /= fincov_ok) then
         if (present(status)) status = local_status
         if (present(selected_k)) selected_k = 0
         return
      end if
      call symmetric_eigen_jacobi(covariance_sample, eigenvalues, eigenvectors, local_status)
      if (local_status /= fincov_ok .and. local_status /= fincov_no_convergence) then
         if (present(status)) status = local_status
         if (present(selected_k)) selected_k = 0
         return
      end if

      if (k_requested == 0) then
         k_use = 0
         do i = 1, k_max
            singular_value = sqrt(max(0.0_dp, eigenvalues(i))*real(n-1,dp))
            if (singular_value > 1.0_dp) k_use = k_use + 1
         end do
      else
         k_use = k_requested
      end if
      if (present(selected_k)) selected_k = k_use

      allocate(common_covariance(p,p), specific_variances(p))
      common_covariance = 0.0_dp
      do i = 1, k_use
         if (eigenvalues(i) > 0.0_dp) then
            common_covariance = common_covariance + eigenvalues(i) * &
               spread(eigenvectors(:,i),2,p) * spread(eigenvectors(:,i),1,p)
         end if
      end do
      do i = 1, p
         specific_variances(i) = covariance_sample(i,i) - common_covariance(i,i)
      end do
      covariance_estimate = common_covariance
      call add_diagonal(covariance_estimate, specific_variances)
      covariance_estimate = 0.5_dp*(covariance_estimate + transpose(covariance_estimate))
      if (present(status)) status = local_status
   end function stat_factor_cov
end module fincov_factor_models
