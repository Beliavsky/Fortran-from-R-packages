! SPDX-License-Identifier: MIT
! Copyright (c) 2026 Dmitriy Mayorov
module vasicekfit_linalg
   use vasicekfit_kinds, only : dp
   implicit none
   private

   public :: invert_matrix, solve_linear_system, ols_regression
   public :: sample_variance, symmetrize_matrix, hac_covariance

contains

   subroutine invert_matrix(a, inverse, ok)
      real(dp), intent(in) :: a(:,:)
      real(dp), allocatable, intent(out) :: inverse(:,:)
      logical, intent(out) :: ok
      real(dp), allocatable :: aug(:,:), row_tmp(:)
      real(dp) :: pivot, factor, scale
      integer :: n, i, j, pivot_row

      n = size(a, 1)
      ok = size(a, 2) == n .and. n > 0
      if (.not. ok) then
         allocate(inverse(0,0))
         return
      end if

      allocate(aug(n, 2*n), row_tmp(2*n), inverse(n,n))
      aug = 0.0_dp
      aug(:, 1:n) = a
      do i = 1, n
         aug(i, n+i) = 1.0_dp
      end do

      scale = max(1.0_dp, maxval(abs(a)))
      do i = 1, n
         pivot_row = i - 1 + maxloc(abs(aug(i:n, i)), dim=1)
         pivot = aug(pivot_row, i)
         if (abs(pivot) <= 100.0_dp * epsilon(1.0_dp) * scale) then
            ok = .false.
            inverse = 0.0_dp
            return
         end if
         if (pivot_row /= i) then
            row_tmp = aug(i,:)
            aug(i,:) = aug(pivot_row,:)
            aug(pivot_row,:) = row_tmp
         end if
         aug(i,:) = aug(i,:) / aug(i,i)
         do j = 1, n
            if (j /= i) then
               factor = aug(j,i)
               aug(j,:) = aug(j,:) - factor * aug(i,:)
            end if
         end do
      end do
      inverse = aug(:, n+1:2*n)
      call symmetrize_matrix(inverse)
   end subroutine invert_matrix

   subroutine solve_linear_system(a, b, x, ok)
      real(dp), intent(in) :: a(:,:), b(:)
      real(dp), allocatable, intent(out) :: x(:)
      logical, intent(out) :: ok
      real(dp), allocatable :: inverse(:,:)

      if (size(a,1) /= size(b)) then
         ok = .false.
         allocate(x(0))
         return
      end if
      call invert_matrix(a, inverse, ok)
      if (ok) then
         x = matmul(inverse, b)
      else
         allocate(x(0))
      end if
   end subroutine solve_linear_system

   subroutine ols_regression(x, y, beta, fitted, residuals, covariance, ok)
      real(dp), intent(in) :: x(:,:), y(:)
      real(dp), allocatable, intent(out) :: beta(:), fitted(:), residuals(:), covariance(:,:)
      logical, intent(out) :: ok
      real(dp), allocatable :: xtx(:,:), xtx_inv(:,:)
      real(dp) :: residual_variance
      integer :: n, k

      n = size(x,1)
      k = size(x,2)
      ok = n == size(y) .and. n > k .and. k > 0
      if (.not. ok) then
         allocate(beta(0), fitted(0), residuals(0), covariance(0,0))
         return
      end if

      xtx = matmul(transpose(x), x)
      call invert_matrix(xtx, xtx_inv, ok)
      if (.not. ok) then
         allocate(beta(0), fitted(0), residuals(0), covariance(0,0))
         return
      end if
      beta = matmul(xtx_inv, matmul(transpose(x), y))
      fitted = matmul(x, beta)
      residuals = y - fitted
      residual_variance = dot_product(residuals, residuals) / real(n - k, dp)
      covariance = residual_variance * xtx_inv
      call symmetrize_matrix(covariance)
   end subroutine ols_regression

   pure real(dp) function sample_variance(x) result(value)
      real(dp), intent(in) :: x(:)
      real(dp) :: mean_x
      if (size(x) < 2) then
         value = 0.0_dp
      else
         mean_x = sum(x) / real(size(x), dp)
         value = sum((x - mean_x)**2) / real(size(x) - 1, dp)
      end if
   end function sample_variance

   subroutine symmetrize_matrix(a)
      real(dp), intent(inout) :: a(:,:)
      integer :: i, j, n
      real(dp) :: value
      n = min(size(a,1), size(a,2))
      do j = 1, n
         do i = j + 1, n
            value = 0.5_dp * (a(i,j) + a(j,i))
            a(i,j) = value
            a(j,i) = value
         end do
      end do
   end subroutine symmetrize_matrix

   function hac_covariance(influence, lag) result(covariance)
      real(dp), intent(in) :: influence(:,:)
      integer, intent(in), optional :: lag
      real(dp), allocatable :: covariance(:,:)
      real(dp), allocatable :: centered(:,:), gamma(:,:)
      real(dp) :: weight
      integer :: n, k, l, max_lag, i, row, col

      n = size(influence,1)
      k = size(influence,2)
      allocate(covariance(k,k))
      covariance = 0.0_dp
      if (n <= 0 .or. k <= 0) return

      allocate(centered(n,k))
      centered = influence
      do i = 1, k
         centered(:,i) = centered(:,i) - sum(centered(:,i)) / real(n,dp)
      end do

      if (present(lag)) then
         max_lag = max(0, min(lag, n-1))
      else
         max_lag = min(n-1, int(floor(4.0_dp * (real(n,dp) / 100.0_dp)**(2.0_dp / 9.0_dp))))
      end if

      covariance = matmul(transpose(centered), centered) / real(n,dp)
      do l = 1, max_lag
         weight = 1.0_dp - real(l,dp) / real(max_lag + 1,dp)
         allocate(gamma(k,k))
         gamma = 0.0_dp
         do i = 1, n - l
            do col = 1, k
               do row = 1, k
                  gamma(row,col) = gamma(row,col) + centered(i+l,row) * centered(i,col)
               end do
            end do
         end do
         gamma = gamma / real(n,dp)
         covariance = covariance + weight * (gamma + transpose(gamma))
         deallocate(gamma)
      end do
      covariance = covariance / real(n,dp)
      call symmetrize_matrix(covariance)
   end function hac_covariance

end module vasicekfit_linalg
