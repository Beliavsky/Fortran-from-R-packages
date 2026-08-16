! Copyright (C) 1998-2025 T. W. Yee, University of Auckland.
! Modern Fortran computational translation, 2026.
! SPDX-License-Identifier: GPL-3.0-only
module vgam_normal_special
   use vgam_kinds, only : dp
   use vgam_distributions, only : dnorm_v, pnorm_v, qnorm_v, rnorm_v
   use vgam_optim, only : bfgs_minimize, numerical_hessian
   use vgam_linalg, only : invert_matrix, weighted_least_squares
   implicit none
   private

   type, public :: tobit_result_t
      real(dp), allocatable :: mean_coefficients(:)
      real(dp), allocatable :: scale_coefficients(:)
      real(dp), allocatable :: covariance(:, :)
      real(dp), allocatable :: fitted_latent_mean(:)
      real(dp), allocatable :: fitted_sd(:)
      real(dp) :: lower = 0.0_dp
      real(dp) :: upper = huge(1.0_dp)
      real(dp) :: loglik = -huge(1.0_dp)
      real(dp) :: aic = huge(1.0_dp)
      integer :: status = 0
      logical :: converged = .false.
   contains
      procedure :: predict => predict_tobit
   end type tobit_result_t

   type, public :: folded_normal_result_t
      real(dp), allocatable :: mean_coefficients(:)
      real(dp), allocatable :: scale_coefficients(:)
      real(dp), allocatable :: covariance(:, :)
      real(dp), allocatable :: fitted_mean_parameter(:)
      real(dp), allocatable :: fitted_sd(:)
      real(dp) :: a1 = 1.0_dp
      real(dp) :: a2 = 1.0_dp
      real(dp) :: loglik = -huge(1.0_dp)
      real(dp) :: aic = huge(1.0_dp)
      integer :: status = 0
      logical :: converged = .false.
   contains
      procedure :: predict => predict_folded_normal
   end type folded_normal_result_t

   public :: dtobit_v, ptobit_v, qtobit_v, rtobit_v, fit_tobit
   public :: dfoldnorm_v, pfoldnorm_v, qfoldnorm_v, rfoldnorm_v, fit_folded_normal

contains

   elemental real(dp) function dtobit_v(x, mean, sd, lower, upper, log_density) result(v)
      real(dp), intent(in) :: x, mean, sd, lower, upper
      logical, intent(in), optional :: log_density
      real(dp) :: ld, p
      logical :: lg
      lg = .false.; if (present(log_density)) lg = log_density
      if (sd <= 0.0_dp .or. lower >= upper) then
         ld = -huge(1.0_dp)
      else if (x < lower .or. x > upper) then
         ld = -huge(1.0_dp)
      else if (x == lower) then
         p = pnorm_v(lower, mean, sd)
         ld = log(max(p, tiny(1.0_dp)))
      else if (x == upper) then
         p = 1.0_dp - pnorm_v(upper, mean, sd)
         ld = log(max(p, tiny(1.0_dp)))
      else
         ld = dnorm_v(x, mean, sd, .true.)
      end if
      if (lg) then
         v = ld
      else if (ld <= -700.0_dp) then
         v = 0.0_dp
      else
         v = exp(ld)
      end if
   end function dtobit_v

   elemental real(dp) function ptobit_v(q, mean, sd, lower, upper) result(p)
      real(dp), intent(in) :: q, mean, sd, lower, upper
      if (sd <= 0.0_dp .or. lower >= upper) then
         p = -1.0_dp
      else if (q < lower) then
         p = 0.0_dp
      else if (q >= upper) then
         p = 1.0_dp
      else
         p = pnorm_v(q, mean, sd)
      end if
   end function ptobit_v

   elemental real(dp) function qtobit_v(p, mean, sd, lower, upper) result(q)
      real(dp), intent(in) :: p, mean, sd, lower, upper
      real(dp) :: pl, pu
      if (p < 0.0_dp .or. p > 1.0_dp .or. sd <= 0.0_dp .or. lower >= upper) then
         q = huge(1.0_dp)
         return
      end if
      pl = pnorm_v(lower, mean, sd)
      pu = pnorm_v(upper, mean, sd)
      if (p <= pl) then
         q = lower
      else if (p >= pu) then
         q = upper
      else
         q = qnorm_v(p, mean, sd)
      end if
   end function qtobit_v

   real(dp) function rtobit_v(mean, sd, lower, upper) result(x)
      real(dp), intent(in) :: mean, sd, lower, upper
      x = rnorm_v(mean, sd)
      x = min(upper, max(lower, x))
   end function rtobit_v

   subroutine fit_tobit(y, x_mean, result, lower, upper, x_scale, weights, max_iter, tol)
      real(dp), intent(in) :: y(:), x_mean(:, :), lower, upper
      type(tobit_result_t), intent(out) :: result
      real(dp), intent(in), optional :: x_scale(:, :), weights(:), tol
      integer, intent(in), optional :: max_iter
      real(dp), allocatable :: xs(:, :), w(:), par(:), h(:, :), cov(:, :), beta0(:), dummy(:, :)
      real(dp) :: fval, tolerance, resid_sd
      integer :: n, pm, ps, np, stat, stat2, niter, i

      n = size(y); pm = size(x_mean, 2)
      if (n < 3 .or. pm < 1 .or. size(x_mean, 1) /= n .or. lower >= upper .or. &
          any(y < lower) .or. any(y > upper)) then
         result%status = 1
         return
      end if
      if (present(x_scale)) then
         if (size(x_scale, 1) /= n .or. size(x_scale, 2) < 1) then
            result%status = 2
            return
         end if
         ps = size(x_scale, 2); allocate(xs(n, ps)); xs = x_scale
      else
         ps = 1; allocate(xs(n, 1)); xs = 1.0_dp
      end if
      allocate(w(n)); w = 1.0_dp
      if (present(weights)) then
         if (size(weights) /= n .or. any(weights < 0.0_dp)) then
            result%status = 3
            return
         end if
         w = weights
      end if
      np = pm + ps; allocate(par(np)); par = 0.0_dp
      call weighted_least_squares(x_mean, y, w, beta0, dummy, stat2)
      if (stat2 == 0) par(1:pm) = beta0
      resid_sd = sqrt(max(sum(w*(y - matmul(x_mean, par(1:pm)))**2)/max(sum(w), 1.0_dp), 1.0e-4_dp))
      if (is_intercept_design(xs)) par(pm + 1) = log(resid_sd)
      niter = 350; if (present(max_iter)) niter = max_iter
      tolerance = 1.0e-7_dp; if (present(tol)) tolerance = tol
      call bfgs_minimize(objective, par, fval, stat, max_iter=niter, tol=tolerance)

      result%mean_coefficients = par(1:pm)
      result%scale_coefficients = par(pm + 1:np)
      result%lower = lower; result%upper = upper
      result%status = stat; result%converged = stat == 0
      result%loglik = -fval; result%aic = 2.0_dp*fval + 2.0_dp*real(np, dp)
      allocate(result%fitted_latent_mean(n), result%fitted_sd(n))
      result%fitted_latent_mean = matmul(x_mean, result%mean_coefficients)
      do i = 1, n
         result%fitted_sd(i) = exp(clamp_eta(dot_product(xs(i, :), result%scale_coefficients)))
      end do
      allocate(h(np, np)); call numerical_hessian(objective, par, h); call invert_matrix(h, cov, stat2)
      if (stat2 == 0) then; result%covariance = cov; else; allocate(result%covariance(0, 0)); end if
   contains
      real(dp) function objective(theta) result(nll)
         real(dp), intent(in) :: theta(:)
         real(dp) :: mu, sdv, ld
         integer :: row
         nll = 0.0_dp
         do row = 1, n
            mu = dot_product(x_mean(row, :), theta(1:pm))
            sdv = exp(clamp_eta(dot_product(xs(row, :), theta(pm + 1:np))))
            ld = dtobit_v(y(row), mu, sdv, lower, upper, .true.)
            if (.not. finite_scalar(ld)) then
               nll = huge(1.0_dp)/100.0_dp; return
            end if
            nll = nll - w(row)*ld
         end do
      end function objective
   end subroutine fit_tobit

   subroutine predict_tobit(self, x_mean, latent_mean, sd, x_scale)
      class(tobit_result_t), intent(in) :: self
      real(dp), intent(in) :: x_mean(:, :)
      real(dp), allocatable, intent(out) :: latent_mean(:), sd(:)
      real(dp), intent(in), optional :: x_scale(:, :)
      real(dp), allocatable :: xs(:, :)
      integer :: n, ps, i
      n = size(x_mean, 1); ps = size(self%scale_coefficients)
      if (size(x_mean, 2) /= size(self%mean_coefficients)) then
         allocate(latent_mean(0), sd(0)); return
      end if
      if (present(x_scale)) then
         if (size(x_scale, 1) /= n .or. size(x_scale, 2) /= ps) then
            allocate(latent_mean(0), sd(0)); return
         end if
         allocate(xs(n, ps)); xs = x_scale
      else
         if (ps /= 1) then; allocate(latent_mean(0), sd(0)); return; end if
         allocate(xs(n, 1)); xs = 1.0_dp
      end if
      allocate(latent_mean(n), sd(n)); latent_mean = matmul(x_mean, self%mean_coefficients)
      do i = 1, n
         sd(i) = exp(clamp_eta(dot_product(xs(i, :), self%scale_coefficients)))
      end do
   end subroutine predict_tobit

   elemental real(dp) function dfoldnorm_v(x, mean, sd, a1, a2, log_density) result(v)
      real(dp), intent(in) :: x, mean, sd, a1, a2
      logical, intent(in), optional :: log_density
      real(dp) :: d
      logical :: lg
      lg = .false.; if (present(log_density)) lg = log_density
      if (x < 0.0_dp .or. sd <= 0.0_dp .or. a1 <= 0.0_dp .or. a2 <= 0.0_dp) then
         d = 0.0_dp
      else
         d = dnorm_v(x/(a1*sd) - mean/sd, 0.0_dp, 1.0_dp)/(a1*sd) + &
             dnorm_v(x/(a2*sd) + mean/sd, 0.0_dp, 1.0_dp)/(a2*sd)
      end if
      if (lg) then
         v = log(max(d, tiny(1.0_dp)))
      else
         v = d
      end if
   end function dfoldnorm_v

   elemental real(dp) function pfoldnorm_v(q, mean, sd, a1, a2) result(p)
      real(dp), intent(in) :: q, mean, sd, a1, a2
      if (sd <= 0.0_dp .or. a1 <= 0.0_dp .or. a2 <= 0.0_dp) then
         p = -1.0_dp
      else if (q <= 0.0_dp) then
         p = 0.0_dp
      else
         p = pnorm_v(q/(a1*sd) - mean/sd, 0.0_dp, 1.0_dp) - &
             pnorm_v(-q/(a2*sd) - mean/sd, 0.0_dp, 1.0_dp)
         p = min(1.0_dp, max(0.0_dp, p))
      end if
   end function pfoldnorm_v

   real(dp) function qfoldnorm_v(p, mean, sd, a1, a2) result(q)
      real(dp), intent(in) :: p, mean, sd, a1, a2
      real(dp) :: lo, hi, mid, pmid
      integer :: it
      if (p < 0.0_dp .or. p > 1.0_dp .or. sd <= 0.0_dp .or. a1 <= 0.0_dp .or. a2 <= 0.0_dp) then
         q = huge(1.0_dp); return
      end if
      if (p == 0.0_dp) then; q = 0.0_dp; return; end if
      lo = 0.0_dp
      hi = max(a1, a2)*(abs(mean) + 8.0_dp*sd + 1.0_dp)
      do while (pfoldnorm_v(hi, mean, sd, a1, a2) < p)
         hi = 2.0_dp*hi
         if (hi > 1.0e100_dp) exit
      end do
      do it = 1, 100
         mid = 0.5_dp*(lo + hi); pmid = pfoldnorm_v(mid, mean, sd, a1, a2)
         if (pmid < p) then; lo = mid; else; hi = mid; end if
      end do
      q = 0.5_dp*(lo + hi)
   end function qfoldnorm_v

   real(dp) function rfoldnorm_v(mean, sd, a1, a2) result(y)
      real(dp), intent(in) :: mean, sd, a1, a2
      real(dp) :: x
      if (sd <= 0.0_dp .or. a1 <= 0.0_dp .or. a2 <= 0.0_dp) then
         y = huge(1.0_dp); return
      end if
      x = rnorm_v(mean, sd)
      y = max(a1*x, -a2*x)
   end function rfoldnorm_v

   subroutine fit_folded_normal(y, x_mean, result, a1, a2, x_scale, weights, max_iter, tol)
      real(dp), intent(in) :: y(:), x_mean(:, :), a1, a2
      type(folded_normal_result_t), intent(out) :: result
      real(dp), intent(in), optional :: x_scale(:, :), weights(:), tol
      integer, intent(in), optional :: max_iter
      real(dp), allocatable :: xs(:, :), w(:), par(:), h(:, :), cov(:, :)
      real(dp) :: fval, tolerance, ysd
      integer :: n, pm, ps, np, stat, stat2, niter, i
      n = size(y); pm = size(x_mean, 2)
      if (n < 3 .or. pm < 1 .or. size(x_mean, 1) /= n .or. any(y < 0.0_dp) .or. &
          a1 <= 0.0_dp .or. a2 <= 0.0_dp) then
         result%status = 1; return
      end if
      if (present(x_scale)) then
         if (size(x_scale, 1) /= n .or. size(x_scale, 2) < 1) then; result%status = 2; return; end if
         ps = size(x_scale, 2); allocate(xs(n, ps)); xs = x_scale
      else
         ps = 1; allocate(xs(n, 1)); xs = 1.0_dp
      end if
      allocate(w(n)); w = 1.0_dp
      if (present(weights)) then
         if (size(weights) /= n .or. any(weights < 0.0_dp)) then; result%status = 3; return; end if
         w = weights
      end if
      np = pm + ps; allocate(par(np)); par = 0.0_dp
      ysd = sqrt(max(sum(w*(y - sum(w*y)/max(sum(w), tiny(1.0_dp)))**2)/max(sum(w), 1.0_dp), 1.0e-4_dp))
      if (is_intercept_design(xs)) par(pm + 1) = log(ysd/max(a1, a2))
      niter = 350; if (present(max_iter)) niter = max_iter
      tolerance = 1.0e-7_dp; if (present(tol)) tolerance = tol
      call bfgs_minimize(objective, par, fval, stat, max_iter=niter, tol=tolerance)
      result%mean_coefficients = par(1:pm); result%scale_coefficients = par(pm + 1:np)
      result%a1 = a1; result%a2 = a2; result%status = stat; result%converged = stat == 0
      result%loglik = -fval; result%aic = 2.0_dp*fval + 2.0_dp*real(np, dp)
      allocate(result%fitted_mean_parameter(n), result%fitted_sd(n))
      result%fitted_mean_parameter = matmul(x_mean, result%mean_coefficients)
      do i = 1, n
         result%fitted_sd(i) = exp(clamp_eta(dot_product(xs(i, :), result%scale_coefficients)))
      end do
      allocate(h(np, np)); call numerical_hessian(objective, par, h); call invert_matrix(h, cov, stat2)
      if (stat2 == 0) then; result%covariance = cov; else; allocate(result%covariance(0, 0)); end if
   contains
      real(dp) function objective(theta) result(nll)
         real(dp), intent(in) :: theta(:)
         real(dp) :: mu, sdv, ld
         integer :: row
         nll = 0.0_dp
         do row = 1, n
            mu = dot_product(x_mean(row, :), theta(1:pm))
            sdv = exp(clamp_eta(dot_product(xs(row, :), theta(pm + 1:np))))
            ld = dfoldnorm_v(y(row), mu, sdv, a1, a2, .true.)
            if (.not. finite_scalar(ld)) then; nll = huge(1.0_dp)/100.0_dp; return; end if
            nll = nll - w(row)*ld
         end do
      end function objective
   end subroutine fit_folded_normal

   subroutine predict_folded_normal(self, x_mean, mean_parameter, sd, x_scale)
      class(folded_normal_result_t), intent(in) :: self
      real(dp), intent(in) :: x_mean(:, :)
      real(dp), allocatable, intent(out) :: mean_parameter(:), sd(:)
      real(dp), intent(in), optional :: x_scale(:, :)
      real(dp), allocatable :: xs(:, :)
      integer :: n, ps, i
      n = size(x_mean, 1); ps = size(self%scale_coefficients)
      if (size(x_mean, 2) /= size(self%mean_coefficients)) then
         allocate(mean_parameter(0), sd(0)); return
      end if
      if (present(x_scale)) then
         if (size(x_scale, 1) /= n .or. size(x_scale, 2) /= ps) then
            allocate(mean_parameter(0), sd(0)); return
         end if
         allocate(xs(n, ps)); xs = x_scale
      else
         if (ps /= 1) then; allocate(mean_parameter(0), sd(0)); return; end if
         allocate(xs(n, 1)); xs = 1.0_dp
      end if
      allocate(mean_parameter(n), sd(n)); mean_parameter = matmul(x_mean, self%mean_coefficients)
      do i = 1, n
         sd(i) = exp(clamp_eta(dot_product(xs(i, :), self%scale_coefficients)))
      end do
   end subroutine predict_folded_normal

   logical function is_intercept_design(x) result(ok)
      real(dp), intent(in) :: x(:, :)
      ok = size(x, 2) >= 1 .and. maxval(abs(x(:, 1) - 1.0_dp)) < 1.0e-12_dp
   end function is_intercept_design

   elemental real(dp) function clamp_eta(x) result(y)
      real(dp), intent(in) :: x
      y = min(30.0_dp, max(-30.0_dp, x))
   end function clamp_eta

   elemental logical function finite_scalar(x) result(ok)
      use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
      real(dp), intent(in) :: x
      ok = ieee_is_finite(x)
   end function finite_scalar
end module vgam_normal_special
