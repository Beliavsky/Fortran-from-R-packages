! SPDX-License-Identifier: GPL-3.0-only
!
! Experimental modern Fortran translation of computational methods from the
! R package rmgarch, copyright (C) 2008-2025 Alexios Galanos.
! Fortran translation copyright (C) 2026 translation contributors.
! Distributed under the GNU General Public License, version 3 only.
module rmgarch_varx
   use rmgarch_kinds, only : dp
   use rmgarch_math, only : covariance_matrix, inverse_spd, cholesky_lower
   use rmgarch_rng, only : random_normal
   use rmgarch_types, only : varx_fit_result
   implicit none
   private

   public :: fit_varx, fit_varx_robust, filter_varx, forecast_varx, simulate_varx

contains

   function fit_varx(x, order, include_intercept, exogen, ridge) result(fit)
      real(dp), intent(in) :: x(:,:)
      integer, intent(in) :: order
      logical, intent(in), optional :: include_intercept
      real(dp), intent(in), optional :: exogen(:,:)
      real(dp), intent(in), optional :: ridge
      type(varx_fit_result) :: fit
      real(dp), allocatable :: design(:,:), target(:,:), xtx(:,:), xty(:,:), inv(:,:)
      real(dp) :: ridge_value
      integer :: n, m, kexo, cols, rows, t, lag, col, i
      logical :: intercept, ok

      n = size(x,1); m = size(x,2)
      intercept = .true.; if (present(include_intercept)) intercept = include_intercept
      kexo = 0
      if (present(exogen)) then
         if (size(exogen,1) /= n) then
            fit%status = 2; return
         end if
         kexo = size(exogen,2)
      end if
      if (order < 0 .or. n <= order .or. m < 1) then
         fit%status = 2; return
      end if
      rows = n-order
      cols = order*m+merge(1,0,intercept)+kexo
      if (cols == 0) then
         fit%status = 2; return
      end if
      allocate(design(rows,cols),target(rows,m),xtx(cols,cols),xty(cols,m),inv(cols,cols))
      do t = 1, rows
         col = 0
         do lag = 1, order
            design(t,col+1:col+m) = x(order+t-lag,:)
            col = col+m
         end do
         if (intercept) then
            col = col+1; design(t,col) = 1.0_dp
         end if
         if (present(exogen)) design(t,col+1:col+kexo) = exogen(order+t,:)
         target(t,:) = x(order+t,:)
      end do
      xtx = matmul(transpose(design),design)
      ridge_value = 1.0e-10_dp; if (present(ridge)) ridge_value = max(0.0_dp,ridge)
      do i = 1, cols
         xtx(i,i) = xtx(i,i)+ridge_value
      end do
      xty = matmul(transpose(design),target)
      inv = inverse_spd(xtx,ok)
      if (.not. ok) then
         fit%status = 3; return
      end if
      fit%order = order; fit%include_intercept = intercept
      allocate(fit%coefficients(cols,m),fit%residuals(rows,m),fit%sigma(m,m))
      fit%coefficients = matmul(inv,xty)
      fit%residuals = target-matmul(design,fit%coefficients)
      fit%sigma = covariance_matrix(fit%residuals)
      fit%status = 0
   end function fit_varx

   subroutine filter_varx(x, fit, fitted, residuals, exogen)
      real(dp), intent(in) :: x(:,:)
      type(varx_fit_result), intent(in) :: fit
      real(dp), intent(out) :: fitted(size(x,1)-fit%order,size(x,2))
      real(dp), intent(out) :: residuals(size(x,1)-fit%order,size(x,2))
      real(dp), intent(in), optional :: exogen(:,:)
      real(dp) :: row(size(fit%coefficients,1))
      integer :: t, lag, col, m, kexo
      m = size(x,2); kexo = 0
      if (present(exogen)) kexo = size(exogen,2)
      do t = 1, size(fitted,1)
         col = 0
         do lag = 1, fit%order
            row(col+1:col+m) = x(fit%order+t-lag,:); col = col+m
         end do
         if (fit%include_intercept) then
            col = col+1; row(col) = 1.0_dp
         end if
         if (present(exogen)) row(col+1:col+kexo) = exogen(fit%order+t,:)
         fitted(t,:) = matmul(row,fit%coefficients)
         residuals(t,:) = x(fit%order+t,:)-fitted(t,:)
      end do
   end subroutine filter_varx

   subroutine forecast_varx(history, fit, horizons, forecast, exogen_future)
      real(dp), intent(in) :: history(:,:)
      type(varx_fit_result), intent(in) :: fit
      integer, intent(in) :: horizons
      real(dp), intent(out) :: forecast(horizons,size(history,2))
      real(dp), intent(in), optional :: exogen_future(:,:)
      real(dp), allocatable :: work(:,:)
      real(dp) :: row(size(fit%coefficients,1))
      integer :: n, m, h, lag, col, kexo
      n = size(history,1); m = size(history,2); kexo = 0
      if (present(exogen_future)) kexo = size(exogen_future,2)
      allocate(work(n+horizons,m)); work(1:n,:) = history
      do h = 1, horizons
         col = 0
         do lag = 1, fit%order
            row(col+1:col+m) = work(n+h-lag,:); col = col+m
         end do
         if (fit%include_intercept) then
            col = col+1; row(col) = 1.0_dp
         end if
         if (present(exogen_future)) row(col+1:col+kexo) = exogen_future(h,:)
         forecast(h,:) = matmul(row,fit%coefficients)
         work(n+h,:) = forecast(h,:)
      end do
   end subroutine forecast_varx

   subroutine simulate_varx(nobs, fit, initial, simulated, exogen)
      integer, intent(in) :: nobs
      type(varx_fit_result), intent(in) :: fit
      real(dp), intent(in) :: initial(:,:)
      real(dp), intent(out) :: simulated(nobs,size(initial,2))
      real(dp), intent(in), optional :: exogen(:,:)
      real(dp), allocatable :: work(:,:)
      real(dp) :: row(size(fit%coefficients,1)), shock(size(initial,2)), l(size(initial,2),size(initial,2))
      integer :: p, m, t, lag, col, j, kexo
      logical :: ok
      p = fit%order; m = size(initial,2); kexo = 0
      if (present(exogen)) kexo = size(exogen,2)
      allocate(work(p+nobs,m)); work(1:p,:) = initial(size(initial,1)-p+1:size(initial,1),:)
      call cholesky_lower(fit%sigma,l,ok)
      if (.not. ok) l = 0.0_dp
      do t = 1, nobs
         col = 0
         do lag = 1, p
            row(col+1:col+m) = work(p+t-lag,:); col = col+m
         end do
         if (fit%include_intercept) then
            col = col+1; row(col) = 1.0_dp
         end if
         if (present(exogen)) row(col+1:col+kexo) = exogen(t,:)
         do j = 1, m
            shock(j) = random_normal()
         end do
         work(p+t,:) = matmul(row,fit%coefficients)+matmul(l,shock)
         simulated(t,:) = work(p+t,:)
      end do
   end subroutine simulate_varx

   function fit_varx_robust(x, order, include_intercept, exogen, tuning, max_iterations, tolerance, ridge) result(fit)
      !! Huber iteratively reweighted least-squares VARX estimator.
      real(dp), intent(in) :: x(:,:)
      integer, intent(in) :: order
      logical, intent(in), optional :: include_intercept
      real(dp), intent(in), optional :: exogen(:,:)
      real(dp), intent(in), optional :: tuning, tolerance, ridge
      integer, intent(in), optional :: max_iterations
      type(varx_fit_result) :: fit
      real(dp), allocatable :: design(:,:), target(:,:), coefficients(:,:), previous(:,:)
      real(dp), allocatable :: residuals(:,:), weights(:), xtwx(:,:), xtwy(:,:), inv(:,:)
      real(dp), allocatable :: norms(:), work(:)
      real(dp) :: huber, tol, ridge_value, scale, change
      integer :: n, m, kexo, cols, rows, t, lag, col, i, iter, maxit
      logical :: intercept, ok

      n = size(x,1)
      m = size(x,2)
      intercept = .true.
      if (present(include_intercept)) intercept = include_intercept
      kexo = 0
      if (present(exogen)) then
         if (size(exogen,1) /= n) then
            fit%status = 2
            return
         end if
         kexo = size(exogen,2)
      end if
      if (order < 0 .or. n <= order .or. m < 1) then
         fit%status = 2
         return
      end if
      rows = n-order
      cols = order*m+merge(1,0,intercept)+kexo
      if (cols < 1) then
         fit%status = 2
         return
      end if
      huber = 1.5_dp
      if (present(tuning)) huber = max(tuning,0.1_dp)
      tol = 1.0e-7_dp
      if (present(tolerance)) tol = max(tolerance,epsilon(1.0_dp))
      maxit = 100
      if (present(max_iterations)) maxit = max(1,max_iterations)
      ridge_value = 1.0e-10_dp
      if (present(ridge)) ridge_value = max(0.0_dp,ridge)

      allocate(design(rows,cols),target(rows,m),coefficients(cols,m),previous(cols,m))
      allocate(residuals(rows,m),weights(rows),xtwx(cols,cols),xtwy(cols,m),inv(cols,cols))
      allocate(norms(rows),work(rows))
      do t = 1, rows
         col = 0
         do lag = 1, order
            design(t,col+1:col+m) = x(order+t-lag,:)
            col = col+m
         end do
         if (intercept) then
            col = col+1
            design(t,col) = 1.0_dp
         end if
         if (present(exogen)) design(t,col+1:col+kexo) = exogen(order+t,:)
         target(t,:) = x(order+t,:)
      end do

      weights = 1.0_dp
      coefficients = 0.0_dp
      do iter = 1, maxit
         previous = coefficients
         xtwx = 0.0_dp
         xtwy = 0.0_dp
         do t = 1, rows
            do i = 1, cols
               xtwx(i,:) = xtwx(i,:)+weights(t)*design(t,i)*design(t,:)
               xtwy(i,:) = xtwy(i,:)+weights(t)*design(t,i)*target(t,:)
            end do
         end do
         do i = 1, cols
            xtwx(i,i) = xtwx(i,i)+ridge_value
         end do
         inv = inverse_spd(xtwx,ok)
         if (.not. ok) then
            fit%status = 3
            return
         end if
         coefficients = matmul(inv,xtwy)
         residuals = target-matmul(design,coefficients)
         norms = sqrt(sum(residuals*residuals,dim=2))
         work = norms
         call sort_vector(work)
         if (mod(rows,2) == 1) then
            scale = work((rows+1)/2)
         else
            scale = 0.5_dp*(work(rows/2)+work(rows/2+1))
         end if
         scale = max(scale/0.6744897501960817_dp,1.0e-12_dp)
         do t = 1, rows
            if (norms(t) <= huber*scale) then
               weights(t) = 1.0_dp
            else
               weights(t) = huber*scale/max(norms(t),tiny(1.0_dp))
            end if
         end do
         change = maxval(abs(coefficients-previous))/max(1.0_dp,maxval(abs(previous)))
         if (iter > 1 .and. change < tol) exit
      end do

      fit%order = order
      fit%include_intercept = intercept
      allocate(fit%coefficients(cols,m),fit%residuals(rows,m),fit%sigma(m,m))
      fit%coefficients = coefficients
      fit%residuals = residuals
      fit%sigma = covariance_matrix(residuals)
      fit%status = merge(0,1,iter <= maxit)
   end function fit_varx_robust

   subroutine sort_vector(x)
      real(dp), intent(inout) :: x(:)
      real(dp) :: key
      integer :: i, j
      do i = 2, size(x)
         key = x(i)
         j = i-1
         do while (j >= 1)
            if (x(j) <= key) exit
            x(j+1) = x(j)
            j = j-1
         end do
         x(j+1) = key
      end do
   end subroutine sort_vector

end module rmgarch_varx
