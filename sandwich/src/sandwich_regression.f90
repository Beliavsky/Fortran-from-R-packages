! SPDX-License-Identifier: GPL-2.0-only OR GPL-3.0-only
module sandwich_regression
   use sandwich_kinds, only : dp
   use sandwich_status, only : SANDWICH_SUCCESS, SANDWICH_INVALID_ARGUMENT, &
      SANDWICH_DIMENSION_MISMATCH
   use sandwich_linalg, only : ols_coefficients, inverse_matrix
   implicit none
   private

   type, public :: ols_model
      integer :: nobs = 0
      integer :: ncoef = 0
      integer :: df_residual = 0
      real(dp) :: sigma2 = 0.0_dp
      real(dp), allocatable :: coefficients(:)
      real(dp), allocatable :: fitted(:)
      real(dp), allocatable :: residuals(:)
      real(dp), allocatable :: weights(:)
      real(dp), allocatable :: bread(:, :)
      real(dp), allocatable :: covariance(:, :)
      real(dp), allocatable :: scores(:, :)
      real(dp), allocatable :: hat(:)
   end type ols_model

   public :: fit_ols, ols_scores, ols_bread, ols_hatvalues

contains

   subroutine fit_ols(x, y, model, status, weights, offset)
      real(dp), intent(in) :: x(:, :), y(:)
      type(ols_model), intent(out) :: model
      integer, intent(out), optional :: status
      real(dp), intent(in), optional :: weights(:), offset(:)
      real(dp), allocatable :: w(:), off(:), response(:), xtwx(:, :), xtwx_inv(:, :)
      real(dp) :: weighted_sse
      integer :: n, k, i, j, info

      n = size(x, 1)
      k = size(x, 2)
      if (n <= k .or. k <= 0 .or. size(y) /= n) then
         if (present(status)) status = SANDWICH_DIMENSION_MISMATCH
         return
      end if

      allocate(w(n), off(n), response(n))
      w = 1.0_dp
      off = 0.0_dp
      if (present(weights)) then
         if (size(weights) /= n .or. any(weights < 0.0_dp)) then
            if (present(status)) status = SANDWICH_INVALID_ARGUMENT
            return
         end if
         w = weights
      end if
      if (present(offset)) then
         if (size(offset) /= n) then
            if (present(status)) status = SANDWICH_DIMENSION_MISMATCH
            return
         end if
         off = offset
      end if
      response = y - off

      call ols_coefficients(x, response, model%coefficients, info, w)
      if (info /= SANDWICH_SUCCESS) then
         if (present(status)) status = info
         return
      end if

      model%nobs = n
      model%ncoef = k
      model%df_residual = n - k
      allocate(model%fitted(n), model%residuals(n), model%weights(n))
      allocate(model%scores(n, k), model%hat(n))
      model%weights = w
      model%fitted = matmul(x, model%coefficients) + off
      model%residuals = y - model%fitted

      do i = 1, n
         model%scores(i, :) = model%residuals(i) * w(i) * x(i, :)
      end do

      xtwx = 0.0_dp * matmul(transpose(x), x)
      do i = 1, n
         do j = 1, k
            xtwx(j, :) = xtwx(j, :) + w(i) * x(i, j) * x(i, :)
         end do
      end do
      call inverse_matrix(xtwx, xtwx_inv, info)
      if (info /= SANDWICH_SUCCESS) then
         if (present(status)) status = info
         return
      end if

      model%bread = real(n, dp) * xtwx_inv
      weighted_sse = sum(w * model%residuals**2)
      model%sigma2 = weighted_sse / real(model%df_residual, dp)
      model%covariance = model%sigma2 * xtwx_inv
      do i = 1, n
         model%hat(i) = w(i) * dot_product(x(i, :), matmul(xtwx_inv, x(i, :)))
      end do

      if (present(status)) status = SANDWICH_SUCCESS
   end subroutine fit_ols

   subroutine ols_scores(x, residuals, scores, status, weights)
      real(dp), intent(in) :: x(:, :), residuals(:)
      real(dp), allocatable, intent(out) :: scores(:, :)
      integer, intent(out), optional :: status
      real(dp), intent(in), optional :: weights(:)
      integer :: n, i

      n = size(x, 1)
      if (size(residuals) /= n) then
         allocate(scores(0, 0))
         if (present(status)) status = SANDWICH_DIMENSION_MISMATCH
         return
      end if
      allocate(scores(n, size(x, 2)))
      if (present(weights)) then
         if (size(weights) /= n) then
            deallocate(scores)
            allocate(scores(0, 0))
            if (present(status)) status = SANDWICH_DIMENSION_MISMATCH
            return
         end if
         do i = 1, n
            scores(i, :) = residuals(i) * weights(i) * x(i, :)
         end do
      else
         do i = 1, n
            scores(i, :) = residuals(i) * x(i, :)
         end do
      end if
      if (present(status)) status = SANDWICH_SUCCESS
   end subroutine ols_scores

   subroutine ols_bread(x, bread, status, weights)
      real(dp), intent(in) :: x(:, :)
      real(dp), allocatable, intent(out) :: bread(:, :)
      integer, intent(out), optional :: status
      real(dp), intent(in), optional :: weights(:)
      real(dp), allocatable :: information(:, :)
      integer :: n, k, i, j, info

      n = size(x, 1)
      k = size(x, 2)
      if (n <= 0 .or. k <= 0) then
         allocate(bread(0, 0))
         if (present(status)) status = SANDWICH_INVALID_ARGUMENT
         return
      end if
      allocate(information(k, k))
      information = 0.0_dp
      if (present(weights)) then
         if (size(weights) /= n) then
            deallocate(information)
            allocate(bread(0, 0))
            if (present(status)) status = SANDWICH_DIMENSION_MISMATCH
            return
         end if
         do i = 1, n
            do j = 1, k
               information(j, :) = information(j, :) + weights(i) * x(i, j) * x(i, :)
            end do
         end do
      else
         information = matmul(transpose(x), x)
      end if
      call inverse_matrix(information, bread, info)
      if (info == SANDWICH_SUCCESS) bread = real(n, dp) * bread
      if (present(status)) status = info
   end subroutine ols_bread

   subroutine ols_hatvalues(x, hat, status, weights)
      real(dp), intent(in) :: x(:, :)
      real(dp), allocatable, intent(out) :: hat(:)
      integer, intent(out), optional :: status
      real(dp), intent(in), optional :: weights(:)
      real(dp), allocatable :: bread(:, :), invinfo(:, :)
      real(dp) :: wi
      integer :: n, i, info

      n = size(x, 1)
      if (present(weights)) then
         call ols_bread(x, bread, info, weights)
      else
         call ols_bread(x, bread, info)
      end if
      if (info /= SANDWICH_SUCCESS) then
         allocate(hat(0))
         if (present(status)) status = info
         return
      end if
      invinfo = bread / real(n, dp)
      allocate(hat(n))
      do i = 1, n
         wi = 1.0_dp
         if (present(weights)) wi = weights(i)
         hat(i) = wi * dot_product(x(i, :), matmul(invinfo, x(i, :)))
      end do
      if (present(status)) status = SANDWICH_SUCCESS
   end subroutine ols_hatvalues

end module sandwich_regression
