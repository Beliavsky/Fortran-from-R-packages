! SPDX-License-Identifier: GPL-2.0-or-later
! Upstream mitml 0.4-5 (2023-03-08), authored by Simon Grund,
! Alexander Robitzsch, and Oliver Luedtke; upstream license GPL (>= 2).
! Modern free-form Fortran translation for Fortran-from-R-packages.
! Gaussian likelihood helpers translated from internal-methods-likelihood.R.
module mitml_likelihood
   use r_kinds, only : dp
   use r_linalg, only : signed_log_determinant, solve_spd
   use mitml_types, only : MITML_ERR_ARGUMENT, MITML_ERR_DIMENSION, MITML_ERR_LINALG, MITML_OK
   implicit none
   private

   public :: gaussian_lm_loglik
   public :: gaussian_lmm_loglik

contains

   pure subroutine gaussian_lm_loglik(y, x, beta, sigma2, value, status)
      real(dp), intent(in) :: y(:) !! Response vector with one value per observation.
      real(dp), intent(in) :: x(:, :) !! Fixed-effect design matrix with observations in rows.
      real(dp), intent(in) :: beta(:) !! Fixed-effect coefficient vector matching the columns of `x`.
      real(dp), intent(in) :: sigma2 !! Positive Gaussian residual variance.
      real(dp), intent(out) :: value !! Gaussian log likelihood, including the normalizing constant.
      integer, intent(out) :: status !! `MITML_OK` on success or a dimension/argument error code.
      real(dp), allocatable :: residual(:)
      real(dp) :: pi
      integer :: n

      value = 0.0_dp
      status = MITML_OK
      n = size(y)

      if (size(x, 1) /= n .or. size(x, 2) /= size(beta)) then
         status = MITML_ERR_DIMENSION
         return
      end if
      if (sigma2 <= 0.0_dp) then
         status = MITML_ERR_ARGUMENT
         return
      end if

      allocate(residual(n))
      residual = y - matmul(x, beta)
      pi = acos(-1.0_dp)
      value = -0.5_dp * real(n, dp) * log(2.0_dp * pi * sigma2) &
         - 0.5_dp * dot_product(residual, residual) / sigma2
   end subroutine gaussian_lm_loglik

   pure subroutine gaussian_lmm_loglik(y, x, z, cluster, beta, tau, sigma2, value, status)
      real(dp), intent(in) :: y(:) !! Response vector with one value per observation.
      real(dp), intent(in) :: x(:, :) !! Fixed-effect design matrix with observations in rows.
      real(dp), intent(in) :: z(:, :) !! Random-effect design matrix with observations in rows.
      integer, intent(in) :: cluster(:) !! Integer cluster label for each observation; labels need not be consecutive.
      real(dp), intent(in) :: beta(:) !! Fixed-effect coefficient vector matching the columns of `x`.
      real(dp), intent(in) :: tau(:, :) !! Random-effect covariance matrix matching the columns of `z`.
      real(dp), intent(in) :: sigma2 !! Positive level-1 Gaussian residual variance.
      real(dp), intent(out) :: value !! Upstream mitml LMM log-likelihood kernel, excluding the constant `-n*log(2*pi)/2`.
      integer, intent(out) :: status !! `MITML_OK` on success or a dimension/argument/linear-algebra error code.
      integer, allocatable :: labels(:)
      real(dp), allocatable :: residual(:), rhs(:), v(:, :), xi(:, :), yi(:), zi(:, :)
      real(dp) :: determinant_sign, log_determinant
      integer :: i, info, j, k, n, n_clusters, ni, p, q, row
      logical :: known_label

      value = 0.0_dp
      status = MITML_OK
      n = size(y)
      p = size(beta)
      q = size(tau, 1)

      if (size(x, 1) /= n .or. size(z, 1) /= n .or. size(cluster) /= n) then
         status = MITML_ERR_DIMENSION
         return
      end if
      if (size(x, 2) /= p .or. size(z, 2) /= q .or. size(tau, 2) /= q) then
         status = MITML_ERR_DIMENSION
         return
      end if
      if (sigma2 <= 0.0_dp) then
         status = MITML_ERR_ARGUMENT
         return
      end if
      if (n == 0) return

      allocate(labels(n))
      n_clusters = 0
      do i = 1, n
         known_label = .false.
         do j = 1, n_clusters
            if (cluster(i) == labels(j)) then
               known_label = .true.
               exit
            end if
         end do
         if (.not. known_label) then
            n_clusters = n_clusters + 1
            labels(n_clusters) = cluster(i)
         end if
      end do

      do k = 1, n_clusters
         ni = count(cluster == labels(k))
         allocate(yi(ni), xi(ni, p), zi(ni, q), residual(ni), rhs(ni), v(ni, ni))
         row = 0
         do i = 1, n
            if (cluster(i) /= labels(k)) cycle
            row = row + 1
            yi(row) = y(i)
            if (p > 0) xi(row, :) = x(i, :)
            if (q > 0) zi(row, :) = z(i, :)
         end do

         residual = yi
         if (p > 0) residual = residual - matmul(xi, beta)
         v = 0.0_dp
         do i = 1, ni
            v(i, i) = sigma2
         end do
         if (q > 0) v = v + matmul(matmul(zi, tau), transpose(zi))

         call signed_log_determinant(v, determinant_sign, log_determinant, info)
         if (info /= 0 .or. determinant_sign <= 0.0_dp) then
            status = MITML_ERR_LINALG
            value = 0.0_dp
            return
         end if
         call solve_spd(v, residual, rhs, info)
         if (info /= 0) then
            status = MITML_ERR_LINALG
            value = 0.0_dp
            return
         end if

         value = value - 0.5_dp * (log_determinant + dot_product(residual, rhs))
         deallocate(yi, xi, zi, residual, rhs, v)
      end do
   end subroutine gaussian_lmm_loglik

end module mitml_likelihood
