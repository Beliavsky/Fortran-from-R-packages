! SPDX-License-Identifier: GPL-2.0-only
! Modern Fortran translation of computational code from the Spillover R package.
module spillover_var
   use spillover_kinds, only : dp
   use spillover_status, only : spillover_success, spillover_invalid_argument, &
      spillover_singular_matrix, set_status
   use spillover_linalg, only : solve_linear_system, identity_matrix, is_finite_matrix
   implicit none
   private

   integer, parameter, public :: var_none = 0
   integer, parameter, public :: var_const = 1
   integer, parameter, public :: var_trend = 2
   integer, parameter, public :: var_const_trend = 3

   type, public :: var_model
      integer :: k = 0
      integer :: p = 0
      integer :: nobs = 0
      integer :: n_effective = 0
      integer :: ncoef = 0
      integer :: deterministic = var_const
      real(dp), allocatable :: ar(:, :, :)
      real(dp), allocatable :: intercept(:)
      real(dp), allocatable :: trend(:)
      real(dp), allocatable :: sigma(:, :)
      real(dp), allocatable :: residuals(:, :)
   end type var_model

   public :: fit_var
   public :: initialize_var_model
   public :: ma_coefficients

contains

   subroutine fit_var(y, p, model, deterministic, info, message)
      real(dp), intent(in) :: y(:, :)
      integer, intent(in) :: p
      type(var_model), intent(out) :: model
      integer, intent(in), optional :: deterministic
      integer, intent(out), optional :: info
      character(len=*), intent(out), optional :: message

      real(dp), allocatable :: x(:, :), target(:, :), xtx(:, :), xty(:, :), beta(:, :)
      real(dp), allocatable :: rhs(:, :)
      integer :: n, k, neff, ndet, ncoef, det, i, j, lag, row, pos, solve_info
      real(dp) :: denominator
      character(len=160) :: solve_message

      call clear_model(model)
      call set_status(info, message, spillover_success, 'success')
      n = size(y, 1)
      k = size(y, 2)
      det = var_const
      if (present(deterministic)) det = deterministic

      if (n < 2 .or. k < 1 .or. p < 1 .or. p >= n) then
         call set_status(info, message, spillover_invalid_argument, &
            'VAR requires n >= 2, k >= 1, and 1 <= p < n')
         return
      end if
      if (det < var_none .or. det > var_const_trend) then
         call set_status(info, message, spillover_invalid_argument, &
            'invalid deterministic specification')
         return
      end if
      if (.not. is_finite_matrix(y)) then
         call set_status(info, message, spillover_invalid_argument, &
            'VAR data contain nonfinite values')
         return
      end if

      ndet = 0
      if (det == var_const .or. det == var_const_trend) ndet = ndet + 1
      if (det == var_trend .or. det == var_const_trend) ndet = ndet + 1
      neff = n - p
      ncoef = k * p + ndet
      if (neff <= ncoef) then
         call set_status(info, message, spillover_invalid_argument, &
            'too few effective observations for the requested VAR')
         return
      end if

      allocate(x(neff, ncoef), target(neff, k))
      x = 0.0_dp
      do row = 1, neff
         i = p + row
         pos = 1
         do lag = 1, p
            x(row, pos:pos + k - 1) = y(i - lag, :)
            pos = pos + k
         end do
         if (det == var_const .or. det == var_const_trend) then
            x(row, pos) = 1.0_dp
            pos = pos + 1
         end if
         if (det == var_trend .or. det == var_const_trend) then
            x(row, pos) = real(i, dp)
         end if
         target(row, :) = y(i, :)
      end do

      xtx = matmul(transpose(x), x)
      xty = matmul(transpose(x), target)
      allocate(rhs(ncoef, k))
      rhs = xty
      call solve_linear_system(xtx, rhs, beta, solve_info, solve_message)
      if (solve_info /= spillover_success) then
         call set_status(info, message, spillover_singular_matrix, trim(solve_message))
         return
      end if

      model%k = k
      model%p = p
      model%nobs = n
      model%n_effective = neff
      model%ncoef = ncoef
      model%deterministic = det
      allocate(model%ar(k, k, p), model%intercept(k), model%trend(k))
      allocate(model%residuals(neff, k), model%sigma(k, k))
      model%ar = 0.0_dp
      model%intercept = 0.0_dp
      model%trend = 0.0_dp

      pos = 1
      do lag = 1, p
         do j = 1, k
            do i = 1, k
               model%ar(i, j, lag) = beta(pos + j - 1, i)
            end do
         end do
         pos = pos + k
      end do
      if (det == var_const .or. det == var_const_trend) then
         model%intercept = beta(pos, :)
         pos = pos + 1
      end if
      if (det == var_trend .or. det == var_const_trend) then
         model%trend = beta(pos, :)
      end if

      model%residuals = target - matmul(x, beta)
      denominator = real(neff - ncoef, dp)
      model%sigma = matmul(transpose(model%residuals), model%residuals) / denominator
      model%sigma = 0.5_dp * (model%sigma + transpose(model%sigma))
   end subroutine fit_var

   subroutine initialize_var_model(ar, sigma, model, intercept, trend, info, message)
      real(dp), intent(in) :: ar(:, :, :)
      real(dp), intent(in) :: sigma(:, :)
      type(var_model), intent(out) :: model
      real(dp), intent(in), optional :: intercept(:)
      real(dp), intent(in), optional :: trend(:)
      integer, intent(out), optional :: info
      character(len=*), intent(out), optional :: message

      integer :: k, p

      call clear_model(model)
      call set_status(info, message, spillover_success, 'success')
      k = size(ar, 1)
      p = size(ar, 3)
      if (k < 1 .or. p < 1 .or. size(ar, 2) /= k .or. &
          size(sigma, 1) /= k .or. size(sigma, 2) /= k) then
         call set_status(info, message, spillover_invalid_argument, &
            'AR and covariance dimensions are inconsistent')
         return
      end if
      if (.not. is_finite_matrix(sigma)) then
         call set_status(info, message, spillover_invalid_argument, &
            'covariance matrix contains nonfinite values')
         return
      end if
      if (present(intercept)) then
         if (size(intercept) /= k) then
            call set_status(info, message, spillover_invalid_argument, &
               'intercept length does not match the VAR dimension')
            return
         end if
      end if
      if (present(trend)) then
         if (size(trend) /= k) then
            call set_status(info, message, spillover_invalid_argument, &
               'trend length does not match the VAR dimension')
            return
         end if
      end if

      model%k = k
      model%p = p
      model%nobs = 0
      model%n_effective = 0
      model%ncoef = 0
      model%deterministic = var_none
      allocate(model%ar(k, k, p), model%sigma(k, k), model%intercept(k), model%trend(k))
      allocate(model%residuals(0, k))
      model%ar = ar
      model%sigma = 0.5_dp * (sigma + transpose(sigma))
      model%intercept = 0.0_dp
      model%trend = 0.0_dp
      if (present(intercept)) then
         model%intercept = intercept
         model%deterministic = var_const
      end if
      if (present(trend)) then
         model%trend = trend
         if (present(intercept)) then
            model%deterministic = var_const_trend
         else
            model%deterministic = var_trend
         end if
      end if
   end subroutine initialize_var_model

   subroutine ma_coefficients(model, horizon, phi, info, message)
      type(var_model), intent(in) :: model
      integer, intent(in) :: horizon
      real(dp), allocatable, intent(out) :: phi(:, :, :)
      integer, intent(out), optional :: info
      character(len=*), intent(out), optional :: message

      integer :: k, step, lag

      call set_status(info, message, spillover_success, 'success')
      k = model%k
      if (k < 1 .or. model%p < 1 .or. .not. allocated(model%ar) .or. horizon < 1) then
         allocate(phi(0, 0, 0))
         call set_status(info, message, spillover_invalid_argument, &
            'invalid VAR model or forecast horizon')
         return
      end if

      allocate(phi(k, k, horizon))
      phi = 0.0_dp
      phi(:, :, 1) = identity_matrix(k)
      do step = 1, horizon - 1
         do lag = 1, min(model%p, step)
            phi(:, :, step + 1) = phi(:, :, step + 1) + &
               matmul(phi(:, :, step - lag + 1), model%ar(:, :, lag))
         end do
      end do
   end subroutine ma_coefficients

   subroutine clear_model(model)
      type(var_model), intent(out) :: model

      model%k = 0
      model%p = 0
      model%nobs = 0
      model%n_effective = 0
      model%ncoef = 0
      model%deterministic = var_const
   end subroutine clear_model

end module spillover_var
