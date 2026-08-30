! SPDX-License-Identifier: GPL-2.0-or-later
! Computational translation of R package vars 1.6-1; see NOTICE.md and UPSTREAM.md.
module vars_bootstrap
   use r_kinds, only : dp
   use r_linalg, only : least_squares_svd
   use vars_types
   use vars_regression, only : rebuild_var_statistics
   use vars_dynamics, only : impulse_response
   use vars_utils, only : quantile_linear, center_columns
   implicit none
   private

   public :: residual_bootstrap_path, bootstrap_irf_indices

contains

   subroutine residual_bootstrap_path(model, indices, ysampled, info)
      type(var_model), intent(in) :: model
      integer, intent(in) :: indices(:)
      real(dp), allocatable, intent(out) :: ysampled(:, :)
      integer, intent(out) :: info
      real(dp), allocatable :: centered_resid(:, :), lagvec(:), reg(:)
      integer :: i, lag, k, p, idx

      if (size(indices) /= model%nobs) then
         info = vars_invalid_argument
         allocate(ysampled(0, 0))
         return
      end if
      if (any(indices < 1) .or. any(indices > model%nobs)) then
         info = vars_invalid_argument
         allocate(ysampled(0, 0))
         return
      end if
      k = model%k
      p = model%p
      allocate(centered_resid(model%nobs, k))
      call center_columns(model%resid, centered_resid)
      allocate(ysampled(model%totobs, k), lagvec(k * p), reg(model%nreg))
      ysampled(1:p, :) = model%y(1:p, :)
      do i = 1, model%nobs
         do lag = 1, p
            lagvec((lag - 1) * k + 1:lag * k) = ysampled(p + i - lag, :)
         end do
         reg(1:k * p) = lagvec
         if (model%nreg > k * p) reg(k * p + 1:model%nreg) = model%x(i, k * p + 1:model%nreg)
         idx = indices(i)
         ysampled(p + i, :) = matmul(model%coef, reg) + centered_resid(idx, :)
      end do
      info = vars_success
   end subroutine residual_bootstrap_path

   subroutine bootstrap_irf_indices(model, indices, n_ahead, orthogonal, cumulative, confidence, &
      lower, upper, info)
      type(var_model), intent(in) :: model
      integer, intent(in) :: indices(:, :)
      integer, intent(in) :: n_ahead
      logical, intent(in) :: orthogonal, cumulative
      real(dp), intent(in) :: confidence
      real(dp), allocatable, intent(out) :: lower(:, :, :), upper(:, :, :)
      integer, intent(out) :: info
      type(var_model) :: boot_model
      real(dp), allocatable :: ysampled(:, :), irf(:, :, :), all_irf(:, :, :, :), sample(:)
      integer :: runs, run, response, shock, h, k
      real(dp) :: alpha

      runs = size(indices, 2)
      if (size(indices, 1) /= model%nobs .or. runs < 1 .or. n_ahead < 0 .or. &
          confidence <= 0.0_dp .or. confidence >= 1.0_dp) then
         info = vars_invalid_argument
         allocate(lower(0, 0, 0), upper(0, 0, 0))
         return
      end if
      k = model%k
      allocate(all_irf(k, k, n_ahead + 1, runs))
      do run = 1, runs
         call residual_bootstrap_path(model, indices(:, run), ysampled, info)
         if (info /= 0) return
         call refit_bootstrap_model(model, ysampled, boot_model, info)
         if (info /= 0) return
         call impulse_response(boot_model, n_ahead, orthogonal, cumulative, irf, info)
         if (info /= 0) return
         all_irf(:, :, :, run) = irf
      end do
      allocate(lower(k, k, n_ahead + 1), upper(k, k, n_ahead + 1), sample(runs))
      alpha = 1.0_dp - confidence
      do h = 1, n_ahead + 1
         do shock = 1, k
            do response = 1, k
               sample = all_irf(response, shock, h, :)
               lower(response, shock, h) = quantile_linear(sample, alpha / 2.0_dp)
               upper(response, shock, h) = quantile_linear(sample, 1.0_dp - alpha / 2.0_dp)
            end do
         end do
      end do
      info = vars_success
   end subroutine bootstrap_irf_indices

   subroutine refit_bootstrap_model(template, y, model, info)
      type(var_model), intent(in) :: template
      real(dp), intent(in) :: y(:, :)
      type(var_model), intent(out) :: model
      integer, intent(out) :: info
      real(dp), allocatable :: x(:, :), response(:, :), beta(:, :)
      integer :: i, j, lag, k, p, rank

      k = template%k
      p = template%p
      allocate(x(template%nobs, template%nreg), response(template%nobs, k), beta(template%nreg, k))
      response = y(p + 1:template%totobs, :)
      do lag = 1, p
         do j = 1, k
            do i = 1, template%nobs
               x(i, (lag - 1) * k + j) = y(p + i - lag, j)
            end do
         end do
      end do
      if (template%nreg > k * p) x(:, k * p + 1:template%nreg) = template%x(:, k * p + 1:template%nreg)
      call least_squares_svd(x, response, beta, rank, info)
      if (info /= 0 .or. rank < template%nreg) then
         info = vars_singular
         return
      end if
      model%p = p
      model%k = k
      model%nobs = template%nobs
      model%totobs = template%totobs
      model%nreg = template%nreg
      model%deterministic = template%deterministic
      model%season = template%season
      model%exogen_cols = template%exogen_cols
      model%y = y
      model%response = response
      model%x = x
      model%coef = transpose(beta)
      allocate(model%active(k, template%nreg), model%df_resid(k))
      model%active = .true.
      model%df_resid = template%nobs - template%nreg
      call rebuild_var_statistics(model, info)
   end subroutine refit_bootstrap_model

end module vars_bootstrap
