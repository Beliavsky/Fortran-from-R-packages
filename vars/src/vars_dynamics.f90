! SPDX-License-Identifier: GPL-2.0-or-later
! Computational translation of R package vars 1.6-1; see NOTICE.md and UPSTREAM.md.
module vars_dynamics
   use r_kinds, only : dp
   use r_distributions, only : r_qnorm
   use r_linalg, only : cholesky_factor, general_real_eigenvalues, inverse_matrix
   use vars_types
   use vars_utils, only : identity_matrix
   implicit none
   private

   public :: phi_from_a, psi_from_a_sigma, var_roots
   public :: forecast_covariance, forecast_var, impulse_response, fevd_var
   public :: structural_impulse_response, structural_fevd
   public :: bq_identification, vec2var_coefficients
   public :: residual_covariance_unbiased

contains

   subroutine phi_from_a(a, nstep, phi, info)
      real(dp), intent(in) :: a(:, :, :)
      integer, intent(in) :: nstep
      real(dp), allocatable, intent(out) :: phi(:, :, :)
      integer, intent(out) :: info
      integer :: h, j, k, p

      k = size(a, 1)
      p = size(a, 3)
      if (size(a, 2) /= k .or. nstep < 0 .or. p < 1) then
         info = vars_invalid_argument
         allocate(phi(0, 0, 0))
         return
      end if
      allocate(phi(k, k, nstep + 1))
      phi = 0.0_dp
      phi(:, :, 1) = identity_matrix(k)
      do h = 1, nstep
         do j = 1, min(p, h)
            phi(:, :, h + 1) = phi(:, :, h + 1) + matmul(phi(:, :, h - j + 1), a(:, :, j))
         end do
      end do
      info = vars_success
   end subroutine phi_from_a

   subroutine psi_from_a_sigma(a, sigma_u, nstep, psi, info)
      real(dp), intent(in) :: a(:, :, :), sigma_u(:, :)
      integer, intent(in) :: nstep
      real(dp), allocatable, intent(out) :: psi(:, :, :)
      integer, intent(out) :: info
      real(dp), allocatable :: phi(:, :, :), lower(:, :)
      integer :: h

      call phi_from_a(a, nstep, phi, info)
      if (info /= 0) then
         allocate(psi(0, 0, 0))
         return
      end if
      call cholesky_factor(sigma_u, lower, info)
      if (info /= 0) then
         allocate(psi(0, 0, 0))
         return
      end if
      allocate(psi(size(phi, 1), size(phi, 2), size(phi, 3)))
      do h = 1, size(phi, 3)
         psi(:, :, h) = matmul(phi(:, :, h), lower)
      end do
   end subroutine psi_from_a_sigma

   subroutine residual_covariance_unbiased(model, sigma, info)
      type(var_model), intent(in) :: model
      real(dp), allocatable, intent(out) :: sigma(:, :)
      integer, intent(out) :: info
      integer :: i, j
      real(dp) :: denom

      if (.not. allocated(model%resid) .or. .not. allocated(model%df_resid)) then
         info = vars_invalid_argument
         allocate(sigma(0, 0))
         return
      end if
      allocate(sigma(model%k, model%k))
      do j = 1, model%k
         do i = 1, model%k
            denom = sqrt(real(model%df_resid(i) * model%df_resid(j), dp))
            sigma(i, j) = dot_product(model%resid(:, i), model%resid(:, j)) / denom
         end do
      end do
      info = vars_success
   end subroutine residual_covariance_unbiased

   subroutine var_roots(a, real_part, imag_part, modulus, info)
      real(dp), intent(in) :: a(:, :, :)
      real(dp), allocatable, intent(out) :: real_part(:), imag_part(:), modulus(:)
      integer, intent(out) :: info
      real(dp), allocatable :: companion(:, :)
      integer :: k, p, kp, i

      k = size(a, 1)
      p = size(a, 3)
      if (size(a, 2) /= k .or. p < 1) then
         info = vars_invalid_argument
         allocate(real_part(0), imag_part(0), modulus(0))
         return
      end if
      kp = k * p
      allocate(companion(kp, kp))
      companion = 0.0_dp
      do i = 1, p
         companion(1:k, (i - 1) * k + 1:i * k) = a(:, :, i)
      end do
      if (p > 1) companion(k + 1:kp, 1:kp - k) = identity_matrix(kp - k)
      call general_real_eigenvalues(companion, real_part, imag_part, info)
      if (info /= 0) then
         allocate(modulus(0))
         return
      end if
      allocate(modulus(size(real_part)))
      modulus = sqrt(real_part ** 2 + imag_part ** 2)
   end subroutine var_roots

   subroutine forecast_covariance(a, sigma_u, n_ahead, covariance, info)
      real(dp), intent(in) :: a(:, :, :), sigma_u(:, :)
      integer, intent(in) :: n_ahead
      real(dp), allocatable, intent(out) :: covariance(:, :, :)
      integer, intent(out) :: info
      real(dp), allocatable :: phi(:, :, :)
      integer :: h, j, k

      if (n_ahead < 1) then
         info = vars_invalid_argument
         allocate(covariance(0, 0, 0))
         return
      end if
      k = size(a, 1)
      call phi_from_a(a, n_ahead - 1, phi, info)
      if (info /= 0) then
         allocate(covariance(0, 0, 0))
         return
      end if
      allocate(covariance(k, k, n_ahead))
      covariance = 0.0_dp
      do h = 1, n_ahead
         do j = 1, h
            covariance(:, :, h) = covariance(:, :, h) + &
               matmul(phi(:, :, j), matmul(sigma_u, transpose(phi(:, :, j))))
         end do
      end do
   end subroutine forecast_covariance

   subroutine forecast_var(model, n_ahead, confidence, result, info, future_nonlag)
      type(var_model), intent(in) :: model
      integer, intent(in) :: n_ahead
      real(dp), intent(in) :: confidence
      type(forecast_result), intent(out) :: result
      integer, intent(out) :: info
      real(dp), intent(in), optional :: future_nonlag(:, :)
      real(dp), allocatable :: nonlag(:, :), lagvec(:), reg(:), sigma(:, :), fcov(:, :, :)
      real(dp) :: critical
      integer :: h, k, p, tail_cols, col, lag, j, cycle, original_t

      if (n_ahead < 1 .or. confidence <= 0.0_dp .or. confidence >= 1.0_dp) then
         info = vars_invalid_argument
         return
      end if
      k = model%k
      p = model%p
      tail_cols = model%nreg - k * p
      allocate(nonlag(n_ahead, tail_cols))
      nonlag = 0.0_dp
      if (present(future_nonlag)) then
         if (size(future_nonlag, 1) /= n_ahead .or. size(future_nonlag, 2) /= tail_cols) then
            info = vars_invalid_argument
            return
         end if
         nonlag = future_nonlag
      else
         if (model%exogen_cols > 0) then
            info = vars_invalid_argument
            return
         end if
         col = 0
         if (model%deterministic == var_const .or. model%deterministic == var_both) then
            col = col + 1
            nonlag(:, col) = 1.0_dp
         end if
         if (model%deterministic == var_trend .or. model%deterministic == var_both) then
            col = col + 1
            do h = 1, n_ahead
               nonlag(h, col) = real(model%totobs + h, dp)
            end do
         end if
         if (model%season > 1) then
            do j = 1, model%season - 1
               col = col + 1
               do h = 1, n_ahead
                  original_t = model%totobs + h
                  cycle = mod(original_t - 1, model%season) + 1
                  nonlag(h, col) = -1.0_dp / real(model%season, dp)
                  if (cycle == j) nonlag(h, col) = nonlag(h, col) + 1.0_dp
               end do
            end do
         end if
      end if

      allocate(result%point(n_ahead, k), result%se(n_ahead, k))
      allocate(result%lower(n_ahead, k), result%upper(n_ahead, k))
      allocate(lagvec(k * p), reg(model%nreg))
      col = 0
      do lag = 1, p
         lagvec((lag - 1) * k + 1:lag * k) = model%y(model%totobs - lag + 1, :)
      end do
      do h = 1, n_ahead
         reg(1:k * p) = lagvec
         if (tail_cols > 0) reg(k * p + 1:model%nreg) = nonlag(h, :)
         result%point(h, :) = matmul(model%coef, reg)
         if (p > 1) lagvec(k + 1:k * p) = lagvec(1:k * (p - 1))
         lagvec(1:k) = result%point(h, :)
      end do
      call residual_covariance_unbiased(model, sigma, info)
      if (info /= 0) return
      call forecast_covariance(model%a, sigma, n_ahead, fcov, info)
      if (info /= 0) return
      critical = -r_qnorm((1.0_dp - confidence) / 2.0_dp)
      do h = 1, n_ahead
         do j = 1, k
            result%se(h, j) = sqrt(max(0.0_dp, fcov(j, j, h)))
         end do
      end do
      result%lower = result%point - critical * result%se
      result%upper = result%point + critical * result%se
      info = vars_success
   end subroutine forecast_var

   subroutine impulse_response(model, n_ahead, orthogonal, cumulative, irf, info)
      type(var_model), intent(in) :: model
      integer, intent(in) :: n_ahead
      logical, intent(in) :: orthogonal, cumulative
      real(dp), allocatable, intent(out) :: irf(:, :, :)
      integer, intent(out) :: info
      real(dp), allocatable :: sigma(:, :)
      integer :: h

      if (orthogonal) then
         call residual_covariance_unbiased(model, sigma, info)
         if (info /= 0) then
            allocate(irf(0, 0, 0))
            return
         end if
         call psi_from_a_sigma(model%a, sigma, n_ahead, irf, info)
      else
         call phi_from_a(model%a, n_ahead, irf, info)
      end if
      if (info /= 0) return
      if (cumulative) then
         do h = 2, size(irf, 3)
            irf(:, :, h) = irf(:, :, h) + irf(:, :, h - 1)
         end do
      end if
   end subroutine impulse_response

   subroutine structural_impulse_response(a_lags, impact, n_ahead, cumulative, irf, info)
      real(dp), intent(in) :: a_lags(:, :, :), impact(:, :)
      integer, intent(in) :: n_ahead
      logical, intent(in) :: cumulative
      real(dp), allocatable, intent(out) :: irf(:, :, :)
      integer, intent(out) :: info
      real(dp), allocatable :: phi(:, :, :)
      integer :: h

      call phi_from_a(a_lags, n_ahead, phi, info)
      if (info /= 0) then
         allocate(irf(0, 0, 0))
         return
      end if
      allocate(irf(size(phi, 1), size(phi, 2), size(phi, 3)))
      do h = 1, size(phi, 3)
         irf(:, :, h) = matmul(phi(:, :, h), impact)
      end do
      if (cumulative) then
         do h = 2, size(irf, 3)
            irf(:, :, h) = irf(:, :, h) + irf(:, :, h - 1)
         end do
      end if
   end subroutine structural_impulse_response

   subroutine fevd_var(model, n_ahead, omega, info)
      type(var_model), intent(in) :: model
      integer, intent(in) :: n_ahead
      real(dp), allocatable, intent(out) :: omega(:, :, :)
      integer, intent(out) :: info
      real(dp), allocatable :: psi(:, :, :), sigma(:, :), fcov(:, :, :)
      real(dp) :: denom
      integer :: h, response, shock, j, k

      k = model%k
      call residual_covariance_unbiased(model, sigma, info)
      if (info /= 0) then
         allocate(omega(0, 0, 0))
         return
      end if
      call psi_from_a_sigma(model%a, sigma, n_ahead - 1, psi, info)
      if (info /= 0) then
         allocate(omega(0, 0, 0))
         return
      end if
      call forecast_covariance(model%a, sigma, n_ahead, fcov, info)
      if (info /= 0) then
         allocate(omega(0, 0, 0))
         return
      end if
      allocate(omega(k, k, n_ahead))
      omega = 0.0_dp
      do h = 1, n_ahead
         do response = 1, k
            denom = fcov(response, response, h)
            do shock = 1, k
               do j = 1, h
                  omega(response, shock, h) = omega(response, shock, h) + psi(response, shock, j) ** 2
               end do
               if (denom > 0.0_dp) omega(response, shock, h) = omega(response, shock, h) / denom
            end do
         end do
      end do
   end subroutine fevd_var

   subroutine structural_fevd(a_lags, impact, n_ahead, omega, info)
      real(dp), intent(in) :: a_lags(:, :, :), impact(:, :)
      integer, intent(in) :: n_ahead
      real(dp), allocatable, intent(out) :: omega(:, :, :)
      integer, intent(out) :: info
      real(dp), allocatable :: irf(:, :, :)
      real(dp) :: denom
      integer :: h, response, shock, j, k

      k = size(a_lags, 1)
      call structural_impulse_response(a_lags, impact, n_ahead - 1, .false., irf, info)
      if (info /= 0) then
         allocate(omega(0, 0, 0))
         return
      end if
      allocate(omega(k, k, n_ahead))
      omega = 0.0_dp
      do h = 1, n_ahead
         do response = 1, k
            denom = 0.0_dp
            do shock = 1, k
               do j = 1, h
                  omega(response, shock, h) = omega(response, shock, h) + irf(response, shock, j) ** 2
               end do
               denom = denom + omega(response, shock, h)
            end do
            if (denom > 0.0_dp) omega(response, :, h) = omega(response, :, h) / denom
         end do
      end do
   end subroutine structural_fevd

   subroutine bq_identification(a_lags, sigma_u, long_run_impact, short_run_impact, info)
      real(dp), intent(in) :: a_lags(:, :, :), sigma_u(:, :)
      real(dp), allocatable, intent(out) :: long_run_impact(:, :), short_run_impact(:, :)
      integer, intent(out) :: info
      real(dp), allocatable :: long_run_matrix(:, :), inv_long(:, :), long_cov(:, :), lower(:, :)
      integer :: j, k

      k = size(a_lags, 1)
      allocate(long_run_matrix(k, k))
      long_run_matrix = identity_matrix(k)
      do j = 1, size(a_lags, 3)
         long_run_matrix = long_run_matrix - a_lags(:, :, j)
      end do
      call inverse_matrix(long_run_matrix, inv_long, info)
      if (info /= 0) then
         allocate(long_run_impact(0, 0), short_run_impact(0, 0))
         return
      end if
      long_cov = matmul(inv_long, matmul(sigma_u, transpose(inv_long)))
      call cholesky_factor(long_cov, lower, info)
      if (info /= 0) then
         allocate(long_run_impact(0, 0), short_run_impact(0, 0))
         return
      end if
      allocate(long_run_impact(k, k), short_run_impact(k, k))
      long_run_impact = lower
      short_run_impact = matmul(long_run_matrix, lower)
   end subroutine bq_identification

   subroutine vec2var_coefficients(pi_matrix, gamma_blocks, specification, a_lags, info)
      real(dp), intent(in) :: pi_matrix(:, :)
      real(dp), intent(in) :: gamma_blocks(:, :, :)
      character(len = *), intent(in) :: specification
      real(dp), allocatable, intent(out) :: a_lags(:, :, :)
      integer, intent(out) :: info
      integer :: k, p, i

      k = size(pi_matrix, 1)
      p = size(gamma_blocks, 3) + 1
      if (size(pi_matrix, 2) /= k .or. size(gamma_blocks, 1) /= k .or. &
          size(gamma_blocks, 2) /= k) then
         info = vars_invalid_argument
         allocate(a_lags(0, 0, 0))
         return
      end if
      allocate(a_lags(k, k, p))
      a_lags = 0.0_dp
      select case (trim(adjustl(specification)))
      case ("transitory")
         a_lags(:, :, 1) = gamma_blocks(:, :, 1) + pi_matrix + identity_matrix(k)
         if (p > 2) then
            do i = 2, p - 1
               a_lags(:, :, i) = gamma_blocks(:, :, i) - gamma_blocks(:, :, i - 1)
            end do
         end if
         a_lags(:, :, p) = -gamma_blocks(:, :, p - 1)
      case ("longrun")
         a_lags(:, :, 1) = gamma_blocks(:, :, 1) + identity_matrix(k)
         if (p > 2) then
            do i = 2, p - 1
               a_lags(:, :, i) = gamma_blocks(:, :, i) - gamma_blocks(:, :, i - 1)
            end do
         end if
         a_lags(:, :, p) = pi_matrix - gamma_blocks(:, :, p - 1)
      case default
         info = vars_invalid_argument
         deallocate(a_lags)
         allocate(a_lags(0, 0, 0))
         return
      end select
      info = vars_success
   end subroutine vec2var_coefficients

end module vars_dynamics
