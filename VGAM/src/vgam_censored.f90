! Copyright (C) 1998-2025 T. W. Yee, University of Auckland.
! Modern Fortran computational translation, 2026.
! SPDX-License-Identifier: GPL-3.0-only
module vgam_censored
   use vgam_kinds, only : dp
   use vgam_distributions, only : dnorm_v, pnorm_v, dexp_v, pexp_v, dpois_v, ppois_v, drayleigh, prayleigh
   use vgam_optim, only : bfgs_minimize, numerical_hessian
   use vgam_linalg, only : invert_matrix
   implicit none
   private

   integer, parameter, public :: censor_exact = 0
   integer, parameter, public :: censor_left = 1
   integer, parameter, public :: censor_right = 2
   integer, parameter, public :: censor_interval = 3
   integer, parameter, public :: censored_normal_family = 1
   integer, parameter, public :: censored_poisson_family = 2
   integer, parameter, public :: censored_exponential_family = 3
   integer, parameter, public :: censored_rayleigh_family = 4

   type, public :: censored_regression_result_t
      real(dp), allocatable :: coefficients(:)
      real(dp), allocatable :: scale_coefficients(:)
      real(dp), allocatable :: covariance(:, :)
      real(dp), allocatable :: fitted(:)
      real(dp), allocatable :: fitted_scale(:)
      real(dp) :: loglik = -huge(1.0_dp)
      real(dp) :: aic = huge(1.0_dp)
      integer :: family = 0
      integer :: status = 0
      logical :: converged = .false.
   contains
      procedure :: predict => predict_censored
   end type censored_regression_result_t

   public :: fit_censored_normal, fit_censored_poisson
   public :: fit_censored_exponential, fit_censored_rayleigh
   public :: censored_normal_logprob, censored_poisson_logprob
   public :: censored_exponential_logprob, censored_rayleigh_logprob

contains

   subroutine fit_censored_normal(lower, upper, ctype, x_mean, result, x_sd, weights, max_iter, tol)
      real(dp), intent(in) :: lower(:), upper(:), x_mean(:, :)
      integer, intent(in) :: ctype(:)
      type(censored_regression_result_t), intent(out) :: result
      real(dp), intent(in), optional :: x_sd(:, :), weights(:), tol
      integer, intent(in), optional :: max_iter
      real(dp), allocatable :: xs(:, :)
      if (present(x_sd)) then
         xs = x_sd
      else
         allocate(xs(size(x_mean, 1), 1)); xs = 1.0_dp
      end if
      call fit_censored_common(lower, upper, ctype, x_mean, xs, censored_normal_family, &
         result, weights, max_iter, tol)
   end subroutine fit_censored_normal

   subroutine fit_censored_poisson(lower, upper, ctype, x_mean, result, weights, max_iter, tol)
      real(dp), intent(in) :: lower(:), upper(:), x_mean(:, :)
      integer, intent(in) :: ctype(:)
      type(censored_regression_result_t), intent(out) :: result
      real(dp), intent(in), optional :: weights(:), tol
      integer, intent(in), optional :: max_iter
      real(dp), allocatable :: xs(:, :)
      allocate(xs(size(x_mean, 1), 0))
      call fit_censored_common(lower, upper, ctype, x_mean, xs, censored_poisson_family, &
         result, weights, max_iter, tol)
   end subroutine fit_censored_poisson

   subroutine fit_censored_exponential(lower, upper, ctype, x_rate, result, weights, max_iter, tol)
      real(dp), intent(in) :: lower(:), upper(:), x_rate(:, :)
      integer, intent(in) :: ctype(:)
      type(censored_regression_result_t), intent(out) :: result
      real(dp), intent(in), optional :: weights(:), tol
      integer, intent(in), optional :: max_iter
      real(dp), allocatable :: xs(:, :)
      allocate(xs(size(x_rate, 1), 0))
      call fit_censored_common(lower, upper, ctype, x_rate, xs, censored_exponential_family, &
         result, weights, max_iter, tol)
   end subroutine fit_censored_exponential

   subroutine fit_censored_rayleigh(lower, upper, ctype, x_scale, result, weights, max_iter, tol)
      real(dp), intent(in) :: lower(:), upper(:), x_scale(:, :)
      integer, intent(in) :: ctype(:)
      type(censored_regression_result_t), intent(out) :: result
      real(dp), intent(in), optional :: weights(:), tol
      integer, intent(in), optional :: max_iter
      real(dp), allocatable :: xs(:, :)
      allocate(xs(size(x_scale, 1), 0))
      call fit_censored_common(lower, upper, ctype, x_scale, xs, censored_rayleigh_family, &
         result, weights, max_iter, tol)
   end subroutine fit_censored_rayleigh

   subroutine fit_censored_common(lower, upper, ctype, x1, x2, family, result, weights, max_iter, tol)
      real(dp), intent(in) :: lower(:), upper(:), x1(:, :), x2(:, :)
      integer, intent(in) :: ctype(:), family
      type(censored_regression_result_t), intent(out) :: result
      real(dp), intent(in), optional :: weights(:), tol
      integer, intent(in), optional :: max_iter
      real(dp), allocatable :: w(:), par(:), hess(:, :), cov(:, :)
      real(dp) :: fval, tolerance, m0, sd0
      integer :: n, p1, p2, np, stat, stat2, niter, i

      n = size(lower); p1 = size(x1, 2); p2 = size(x2, 2)
      if (n <= 0 .or. p1 <= 0 .or. size(upper) /= n .or. size(ctype) /= n .or. &
          size(x1, 1) /= n .or. size(x2, 1) /= n .or. &
          any(ctype < censor_exact) .or. any(ctype > censor_interval)) then
         result%status = 1; return
      end if
      do i = 1, n
         if (ctype(i) == censor_interval .and. upper(i) <= lower(i)) then
            result%status = 2; return
         end if
      end do
      if (family /= censored_normal_family .and. p2 /= 0) then
         result%status = 3; return
      end if
      allocate(w(n)); w = 1.0_dp
      if (present(weights)) then
         if (size(weights) /= n .or. any(weights < 0.0_dp)) then
            result%status = 4; return
         end if
         w = weights
      end if
      np = p1 + p2
      allocate(par(np)); par = 0.0_dp
      m0 = sum(w*lower)/max(sum(w), tiny(1.0_dp))
      if (family == censored_normal_family) then
         if (all(abs(x1(:, 1) - 1.0_dp) < 100.0_dp*epsilon(1.0_dp))) par(1) = m0
         sd0 = sqrt(max(sum(w*(lower - m0)**2)/max(sum(w), 1.0_dp), 1.0e-4_dp))
         if (p2 > 0 .and. all(abs(x2(:, 1) - 1.0_dp) < 100.0_dp*epsilon(1.0_dp))) par(p1 + 1) = log(sd0)
      else if (family == censored_poisson_family) then
         if (all(abs(x1(:, 1) - 1.0_dp) < 100.0_dp*epsilon(1.0_dp))) par(1) = log(max(m0, 0.1_dp))
      else if (family == censored_exponential_family) then
         if (all(abs(x1(:, 1) - 1.0_dp) < 100.0_dp*epsilon(1.0_dp))) par(1) = log(1.0_dp/max(m0, 0.1_dp))
      else if (family == censored_rayleigh_family) then
         if (all(abs(x1(:, 1) - 1.0_dp) < 100.0_dp*epsilon(1.0_dp))) par(1) = log(max(m0*sqrt(2.0_dp/acos(-1.0_dp)), 0.1_dp))
      end if
      niter = 350; if (present(max_iter)) niter = max_iter
      tolerance = 1.0e-7_dp; if (present(tol)) tolerance = tol
      call bfgs_minimize(objective, par, fval, stat, max_iter=niter, tol=tolerance)
      result%status = stat
      result%converged = stat == 0
      result%family = family
      result%coefficients = par(1:p1)
      if (p2 > 0) then
         result%scale_coefficients = par(p1 + 1:np)
      else
         allocate(result%scale_coefficients(0))
      end if
      result%loglik = -fval
      result%aic = 2.0_dp*fval + 2.0_dp*real(np, dp)
      allocate(result%fitted(n))
      if (family == censored_normal_family) then
         result%fitted = matmul(x1, result%coefficients)
         allocate(result%fitted_scale(n)); result%fitted_scale = exp(clamp_eta(matmul(x2, result%scale_coefficients)))
      else if (family == censored_exponential_family) then
         result%fitted = exp(-clamp_eta(matmul(x1, result%coefficients)))
         allocate(result%fitted_scale(0))
      else
         result%fitted = exp(clamp_eta(matmul(x1, result%coefficients)))
         allocate(result%fitted_scale(0))
      end if
      allocate(hess(np, np)); call numerical_hessian(objective, par, hess)
      call invert_matrix(hess, cov, stat2)
      if (stat2 == 0) then
         result%covariance = cov
      else
         allocate(result%covariance(0, 0))
      end if

   contains
      real(dp) function objective(theta) result(nll)
         real(dp), intent(in) :: theta(:)
         real(dp) :: eta1, eta2, a1, a2, lp
         integer :: ii
         nll = 0.0_dp
         do ii = 1, n
            eta1 = dot_product(x1(ii, :), theta(1:p1))
            select case (family)
            case (censored_normal_family)
               eta2 = dot_product(x2(ii, :), theta(p1 + 1:np))
               a1 = eta1; a2 = exp(clamp_eta(eta2))
               lp = censored_normal_logprob(lower(ii), upper(ii), ctype(ii), a1, a2)
            case (censored_poisson_family)
               a1 = exp(clamp_eta(eta1))
               lp = censored_poisson_logprob(lower(ii), upper(ii), ctype(ii), a1)
            case (censored_exponential_family)
               a1 = exp(clamp_eta(eta1))
               lp = censored_exponential_logprob(lower(ii), upper(ii), ctype(ii), a1)
            case (censored_rayleigh_family)
               a1 = exp(clamp_eta(eta1))
               lp = censored_rayleigh_logprob(lower(ii), upper(ii), ctype(ii), a1)
            case default
               lp = -huge(1.0_dp)
            end select
            if (.not. finite_scalar(lp) .or. lp < -huge(1.0_dp)/100.0_dp) then
               nll = huge(1.0_dp)/100.0_dp; return
            end if
            nll = nll - w(ii)*lp
         end do
      end function objective
   end subroutine fit_censored_common

   elemental real(dp) function censored_normal_logprob(lower, upper, ctype, mu, sd) result(lp)
      real(dp), intent(in) :: lower, upper, mu, sd
      integer, intent(in) :: ctype
      real(dp) :: p1, p2
      if (sd <= 0.0_dp) then
         lp = -huge(1.0_dp); return
      end if
      select case (ctype)
      case (censor_exact)
         lp = dnorm_v(lower, mu, sd, .true.)
      case (censor_left)
         p1 = pnorm_v(lower, mu, sd); lp = safe_log(p1)
      case (censor_right)
         p1 = pnorm_v(lower, mu, sd); lp = safe_log1m(p1)
      case (censor_interval)
         p1 = pnorm_v(lower, mu, sd); p2 = pnorm_v(upper, mu, sd); lp = safe_log(p2 - p1)
      case default
         lp = -huge(1.0_dp)
      end select
   end function censored_normal_logprob

   elemental real(dp) function censored_poisson_logprob(lower, upper, ctype, mu) result(lp)
      real(dp), intent(in) :: lower, upper, mu
      integer, intent(in) :: ctype
      real(dp) :: p1, p2
      integer :: lo, hi
      if (mu <= 0.0_dp) then
         lp = -huge(1.0_dp); return
      end if
      lo = nint(lower); hi = nint(upper)
      select case (ctype)
      case (censor_exact)
         lp = dpois_v(lo, mu, .true.)
      case (censor_left)
         p1 = ppois_v(lo - 1, mu); lp = safe_log(p1)
      case (censor_right)
         p1 = ppois_v(lo - 1, mu); lp = safe_log1m(p1)
      case (censor_interval)
         p1 = ppois_v(lo, mu); p2 = ppois_v(hi, mu); lp = safe_log(p2 - p1)
      case default
         lp = -huge(1.0_dp)
      end select
   end function censored_poisson_logprob

   elemental real(dp) function censored_exponential_logprob(lower, upper, ctype, rate) result(lp)
      real(dp), intent(in) :: lower, upper, rate
      integer, intent(in) :: ctype
      real(dp) :: p1, p2
      if (rate <= 0.0_dp .or. lower < 0.0_dp) then
         lp = -huge(1.0_dp); return
      end if
      select case (ctype)
      case (censor_exact)
         lp = dexp_v(lower, rate, .true.)
      case (censor_left)
         p1 = pexp_v(lower, rate); lp = safe_log(p1)
      case (censor_right)
         lp = -rate*lower
      case (censor_interval)
         if (upper <= lower) then
            lp = -huge(1.0_dp)
         else
            p1 = pexp_v(lower, rate); p2 = pexp_v(upper, rate); lp = safe_log(p2 - p1)
         end if
      case default
         lp = -huge(1.0_dp)
      end select
   end function censored_exponential_logprob

   elemental real(dp) function censored_rayleigh_logprob(lower, upper, ctype, scale) result(lp)
      real(dp), intent(in) :: lower, upper, scale
      integer, intent(in) :: ctype
      real(dp) :: p1, p2
      if (scale <= 0.0_dp .or. lower < 0.0_dp) then
         lp = -huge(1.0_dp); return
      end if
      select case (ctype)
      case (censor_exact)
         lp = drayleigh(lower, scale, .true.)
      case (censor_left)
         p1 = prayleigh(lower, scale); lp = safe_log(p1)
      case (censor_right)
         p1 = prayleigh(lower, scale); lp = safe_log1m(p1)
      case (censor_interval)
         p1 = prayleigh(lower, scale); p2 = prayleigh(upper, scale); lp = safe_log(p2 - p1)
      case default
         lp = -huge(1.0_dp)
      end select
   end function censored_rayleigh_logprob

   subroutine predict_censored(self, x, fitted, x_scale, fitted_scale)
      class(censored_regression_result_t), intent(in) :: self
      real(dp), intent(in) :: x(:, :)
      real(dp), allocatable, intent(out) :: fitted(:)
      real(dp), intent(in), optional :: x_scale(:, :)
      real(dp), allocatable, intent(out), optional :: fitted_scale(:)
      real(dp), allocatable :: eta(:)
      integer :: n
      n = size(x, 1)
      if (size(x, 2) /= size(self%coefficients)) then
         allocate(fitted(0)); if (present(fitted_scale)) allocate(fitted_scale(0)); return
      end if
      eta = matmul(x, self%coefficients); allocate(fitted(n))
      select case (self%family)
      case (censored_normal_family)
         fitted = eta
         if (present(fitted_scale)) then
            if (.not. present(x_scale)) then
               allocate(fitted_scale(0))
            else if (size(x_scale, 1) /= n .or. size(x_scale, 2) /= size(self%scale_coefficients)) then
               allocate(fitted_scale(0))
            else
               allocate(fitted_scale(n)); fitted_scale = exp(clamp_eta(matmul(x_scale, self%scale_coefficients)))
            end if
         end if
      case (censored_exponential_family)
         fitted = exp(-clamp_eta(eta)); if (present(fitted_scale)) allocate(fitted_scale(0))
      case default
         fitted = exp(clamp_eta(eta)); if (present(fitted_scale)) allocate(fitted_scale(0))
      end select
   end subroutine predict_censored

   elemental real(dp) function safe_log(p) result(v)
      real(dp), intent(in) :: p
      if (p <= 0.0_dp) then
         v = -huge(1.0_dp)
      else
         v = log(p)
      end if
   end function safe_log

   elemental real(dp) function safe_log1m(p) result(v)
      real(dp), intent(in) :: p
      if (p >= 1.0_dp) then
         v = -huge(1.0_dp)
      else if (p <= 0.0_dp) then
         v = 0.0_dp
      else
         v = log(1.0_dp - p)
      end if
   end function safe_log1m

   elemental real(dp) function clamp_eta(x) result(y)
      real(dp), intent(in) :: x
      y = min(35.0_dp, max(-35.0_dp, x))
   end function clamp_eta

   elemental logical function finite_scalar(x) result(ok)
      real(dp), intent(in) :: x
      ok = abs(x) <= huge(x)
   end function finite_scalar

end module vgam_censored
