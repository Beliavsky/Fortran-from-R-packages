! Copyright (C) 1998-2025 T. W. Yee, University of Auckland.
! Modern Fortran computational translation, 2026.
! SPDX-License-Identifier: GPL-3.0-only
module vgam_zoa_beta_models
   use vgam_kinds, only : dp
   use vgam_distributions, only : dbeta_v
   use vgam_optim, only : bfgs_minimize, numerical_hessian
   use vgam_linalg, only : invert_matrix
   implicit none
   private

   type, public :: zoa_beta_result_t
      real(dp), allocatable :: mean_coefficients(:)
      real(dp), allocatable :: precision_coefficients(:)
      real(dp), allocatable :: zero_coefficients(:)
      real(dp), allocatable :: one_coefficients(:)
      real(dp), allocatable :: covariance(:, :)
      real(dp), allocatable :: fitted_beta_mean(:)
      real(dp), allocatable :: fitted_precision(:)
      real(dp), allocatable :: fitted_pzero(:)
      real(dp), allocatable :: fitted_pone(:)
      real(dp), allocatable :: fitted_mean(:)
      real(dp) :: loglik = -huge(1.0_dp)
      real(dp) :: aic = huge(1.0_dp)
      integer :: status = 0
      logical :: converged = .false.
   contains
      procedure :: predict => predict_zoa_beta
   end type zoa_beta_result_t

   public :: fit_zoa_beta_regression

contains

   subroutine fit_zoa_beta_regression(y, x_mean, x_mass, result, x_precision, weights, max_iter, tol)
      real(dp), intent(in) :: y(:), x_mean(:, :), x_mass(:, :)
      type(zoa_beta_result_t), intent(out) :: result
      real(dp), intent(in), optional :: x_precision(:, :), weights(:), tol
      integer, intent(in), optional :: max_iter
      real(dp), allocatable :: xp(:, :), w(:), par(:), h(:, :), cov(:, :)
      real(dp) :: fval, tolerance, p0, p1, mui, phi, a, b
      integer :: n, pm, pp, pz, np, i, stat, stat2, niter

      n = size(y); pm = size(x_mean, 2); pz = size(x_mass, 2)
      if (n < 4 .or. pm < 1 .or. pz < 1 .or. size(x_mean, 1) /= n .or. &
          size(x_mass, 1) /= n .or. any(y < 0.0_dp) .or. any(y > 1.0_dp)) then
         result%status = 1
         return
      end if
      if (present(x_precision)) then
         if (size(x_precision, 1) /= n .or. size(x_precision, 2) < 1) then
            result%status = 2; return
         end if
         pp = size(x_precision, 2); allocate(xp(n, pp)); xp = x_precision
      else
         pp = 1; allocate(xp(n, 1)); xp = 1.0_dp
      end if
      allocate(w(n)); w = 1.0_dp
      if (present(weights)) then
         if (size(weights) /= n .or. any(weights < 0.0_dp)) then; result%status = 3; return; end if
         w = weights
      end if
      np = pm + pp + 2*pz; allocate(par(np)); par = 0.0_dp
      call initialize(y, x_mean, xp, x_mass, w, par, pm, pp, pz)
      niter = 450; if (present(max_iter)) niter = max_iter
      tolerance = 1.0e-7_dp; if (present(tol)) tolerance = tol
      call bfgs_minimize(objective, par, fval, stat, max_iter=niter, tol=tolerance)

      result%mean_coefficients = par(1:pm)
      result%precision_coefficients = par(pm + 1:pm + pp)
      result%zero_coefficients = par(pm + pp + 1:pm + pp + pz)
      result%one_coefficients = par(pm + pp + pz + 1:np)
      result%status = stat; result%converged = stat == 0
      result%loglik = -fval; result%aic = 2.0_dp*fval + 2.0_dp*real(np, dp)
      allocate(result%fitted_beta_mean(n), result%fitted_precision(n), result%fitted_pzero(n), &
         result%fitted_pone(n), result%fitted_mean(n))
      do i = 1, n
         call fitted_row(result, x_mean(i, :), xp(i, :), x_mass(i, :), mui, phi, p0, p1)
         result%fitted_beta_mean(i) = mui; result%fitted_precision(i) = phi
         result%fitted_pzero(i) = p0; result%fitted_pone(i) = p1
         result%fitted_mean(i) = p1 + (1.0_dp - p0 - p1)*mui
      end do
      allocate(h(np, np)); call numerical_hessian(objective, par, h); call invert_matrix(h, cov, stat2)
      if (stat2 == 0) then; result%covariance = cov; else; allocate(result%covariance(0, 0)); end if
   contains
      real(dp) function objective(theta) result(nll)
         real(dp), intent(in) :: theta(:)
         real(dp) :: mu, ph, z0, z1, den, q0, q1, ld
         integer :: row
         nll = 0.0_dp
         do row = 1, n
            mu = logistic(dot_product(x_mean(row, :), theta(1:pm)))
            ph = exp(clamp_eta(dot_product(xp(row, :), theta(pm + 1:pm + pp))))
            z0 = clamp_eta(dot_product(x_mass(row, :), theta(pm + pp + 1:pm + pp + pz)))
            z1 = clamp_eta(dot_product(x_mass(row, :), theta(pm + pp + pz + 1:np)))
            call softmax_masses(z0, z1, q0, q1, den)
            if (y(row) == 0.0_dp) then
               ld = log(max(q0, tiny(1.0_dp)))
            else if (y(row) == 1.0_dp) then
               ld = log(max(q1, tiny(1.0_dp)))
            else
               a = mu*ph; b = (1.0_dp - mu)*ph
               ld = -log(den) + dbeta_v(y(row), a, b, .true.)
            end if
            if (.not. finite_scalar(ld)) then; nll = huge(1.0_dp)/100.0_dp; return; end if
            nll = nll - w(row)*ld
         end do
      end function objective
   end subroutine fit_zoa_beta_regression

   subroutine predict_zoa_beta(self, x_mean, x_mass, beta_mean, precision, pzero, pone, mean, x_precision)
      class(zoa_beta_result_t), intent(in) :: self
      real(dp), intent(in) :: x_mean(:, :), x_mass(:, :)
      real(dp), allocatable, intent(out) :: beta_mean(:), precision(:), pzero(:), pone(:), mean(:)
      real(dp), intent(in), optional :: x_precision(:, :)
      real(dp), allocatable :: xp(:, :)
      integer :: n, pp, i
      n = size(x_mean, 1); pp = size(self%precision_coefficients)
      if (size(x_mean, 2) /= size(self%mean_coefficients) .or. size(x_mass, 1) /= n .or. &
          size(x_mass, 2) /= size(self%zero_coefficients)) then
         allocate(beta_mean(0), precision(0), pzero(0), pone(0), mean(0)); return
      end if
      if (present(x_precision)) then
         if (size(x_precision, 1) /= n .or. size(x_precision, 2) /= pp) then
            allocate(beta_mean(0), precision(0), pzero(0), pone(0), mean(0)); return
         end if
         allocate(xp(n, pp)); xp = x_precision
      else
         if (pp /= 1) then
            allocate(beta_mean(0), precision(0), pzero(0), pone(0), mean(0)); return
         end if
         allocate(xp(n, 1)); xp = 1.0_dp
      end if
      allocate(beta_mean(n), precision(n), pzero(n), pone(n), mean(n))
      do i = 1, n
         call fitted_row(self, x_mean(i, :), xp(i, :), x_mass(i, :), &
            beta_mean(i), precision(i), pzero(i), pone(i))
         mean(i) = pone(i) + (1.0_dp - pzero(i) - pone(i))*beta_mean(i)
      end do
   end subroutine predict_zoa_beta

   subroutine fitted_row(self, xm, xp, xz, mu, phi, p0, p1)
      class(zoa_beta_result_t), intent(in) :: self
      real(dp), intent(in) :: xm(:), xp(:), xz(:)
      real(dp), intent(out) :: mu, phi, p0, p1
      real(dp) :: z0, z1, den
      mu = logistic(dot_product(xm, self%mean_coefficients))
      phi = exp(clamp_eta(dot_product(xp, self%precision_coefficients)))
      z0 = clamp_eta(dot_product(xz, self%zero_coefficients))
      z1 = clamp_eta(dot_product(xz, self%one_coefficients))
      call softmax_masses(z0, z1, p0, p1, den)
   end subroutine fitted_row

   subroutine initialize(y, xm, xp, xz, w, par, pm, pp, pz)
      real(dp), intent(in) :: y(:), xm(:, :), xp(:, :), xz(:, :), w(:)
      real(dp), intent(out) :: par(:)
      integer, intent(in) :: pm, pp, pz
      real(dp) :: sw, p0, p1, pint, mu, vv, phi
      logical, allocatable :: interior(:)
      real(dp), allocatable :: wi(:)
      integer :: ni
      sw = max(sum(w), tiny(1.0_dp)); par = 0.0_dp
      p0 = sum(w*merge(1.0_dp, 0.0_dp, y == 0.0_dp))/sw
      p1 = sum(w*merge(1.0_dp, 0.0_dp, y == 1.0_dp))/sw
      pint = max(1.0e-4_dp, 1.0_dp - p0 - p1)
      allocate(interior(size(y))); interior = y > 0.0_dp .and. y < 1.0_dp
      ni = count(interior); allocate(wi(max(1, ni)))
      if (ni > 0) then
         wi = pack(w, interior)
         mu = sum(wi*pack(y, interior))/max(sum(wi), tiny(1.0_dp))
         vv = sum(wi*(pack(y, interior) - mu)**2)/max(sum(wi), 1.0_dp)
         phi = mu*(1.0_dp - mu)/max(vv, 1.0e-5_dp) - 1.0_dp
         phi = min(1000.0_dp, max(0.5_dp, phi))
      else
         mu = 0.5_dp; phi = 10.0_dp
      end if
      if (is_intercept_design(xm)) par(1) = logit(mu)
      if (is_intercept_design(xp)) par(pm + 1) = log(phi)
      if (is_intercept_design(xz)) then
         par(pm + pp + 1) = log(max(p0, 1.0e-5_dp)/pint)
         par(pm + pp + pz + 1) = log(max(p1, 1.0e-5_dp)/pint)
      end if
   end subroutine initialize

   subroutine softmax_masses(z0, z1, p0, p1, den)
      real(dp), intent(in) :: z0, z1
      real(dp), intent(out) :: p0, p1, den
      real(dp) :: m, e0, e1, eb
      m = max(0.0_dp, max(z0, z1))
      e0 = exp(z0 - m); e1 = exp(z1 - m); eb = exp(-m)
      den = (e0 + e1 + eb)/eb
      p0 = e0/(e0 + e1 + eb); p1 = e1/(e0 + e1 + eb)
   end subroutine softmax_masses

   logical function is_intercept_design(x) result(ok)
      real(dp), intent(in) :: x(:, :)
      ok = size(x, 2) >= 1 .and. maxval(abs(x(:, 1) - 1.0_dp)) < 1.0e-12_dp
   end function is_intercept_design

   elemental real(dp) function logistic(x) result(p)
      real(dp), intent(in) :: x
      if (x >= 0.0_dp) then; p = 1.0_dp/(1.0_dp + exp(-x)); else; p = exp(x)/(1.0_dp + exp(x)); end if
   end function logistic

   elemental real(dp) function logit(p) result(x)
      real(dp), intent(in) :: p
      real(dp) :: q
      q = min(1.0_dp - 1.0e-8_dp, max(1.0e-8_dp, p)); x = log(q/(1.0_dp - q))
   end function logit

   elemental real(dp) function clamp_eta(x) result(y)
      real(dp), intent(in) :: x
      y = min(30.0_dp, max(-30.0_dp, x))
   end function clamp_eta

   elemental logical function finite_scalar(x) result(ok)
      use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
      real(dp), intent(in) :: x
      ok = ieee_is_finite(x)
   end function finite_scalar
end module vgam_zoa_beta_models
