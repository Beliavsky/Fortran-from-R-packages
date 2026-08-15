! Copyright (C) 1998-2025 T. W. Yee, University of Auckland.
! Modern Fortran computational translation, 2026.
! SPDX-License-Identifier: GPL-3.0-only
module vgam_garma
   use vgam_kinds, only : dp
   use vgam_links, only : link_value, link_inverse, link_identity, link_log, &
      link_logit, link_probit, link_cloglog, link_cauchit, link_reciprocal
   use vgam_vglm, only : vglm_result_t, fit_vglm, family_gaussian, &
      family_poisson, family_binomial
   use vgam_optim, only : bfgs_minimize, numerical_hessian
   use vgam_linalg, only : invert_matrix
   implicit none
   private

   type, public :: garma_result_t
      real(dp), allocatable :: coefficients(:)
      real(dp), allocatable :: ar(:)
      real(dp), allocatable :: covariance(:, :)
      real(dp), allocatable :: fitted(:)
      real(dp), allocatable :: linear_predictor(:)
      real(dp), allocatable :: residuals(:)
      integer :: link = link_identity
      integer :: ar_order = 0
      integer :: iterations = 0
      integer :: status = 0
      logical :: converged = .false.
      real(dp) :: objective = huge(1.0_dp)
      real(dp) :: loglik = -huge(1.0_dp)
      real(dp) :: aic = huge(1.0_dp)
   contains
      procedure :: forecast => forecast_garma
   end type garma_result_t

   public :: fit_garma

contains

   subroutine fit_garma(y, x, ar_order, link_id, result, weights, coefstart, &
                        arstart, max_iter, tol)
      real(dp), intent(in) :: y(:), x(:, :)
      integer, intent(in) :: ar_order, link_id
      type(garma_result_t), intent(out) :: result
      real(dp), intent(in), optional :: weights(:), coefstart(:), arstart(:), tol
      integer, intent(in), optional :: max_iter
      type(vglm_result_t) :: initfit
      real(dp), allocatable :: w(:), par(:), hess(:, :), cov(:, :)
      real(dp) :: fval, tolerance
      integer :: n, p, q, stat, stat2, niter, family

      n = size(y)
      p = size(x, 2)
      q = ar_order
      if (n <= q + 2 .or. p <= 0 .or. size(x, 1) /= n .or. q < 0) then
         result%status = 1
         return
      end if
      if (.not. supported_link(link_id)) then
         result%status = 2
         return
      end if
      allocate(w(n))
      w = 1.0_dp
      if (present(weights)) then
         if (size(weights) /= n .or. any(weights < 0.0_dp)) then
            result%status = 3
            return
         end if
         w = weights
      end if

      select case (link_id)
      case (link_identity)
         family = family_gaussian
      case (link_log, link_reciprocal)
         family = family_poisson
      case default
         family = family_binomial
      end select
      if (family == family_binomial) then
         if (any(y < 0.0_dp) .or. any(y > 1.0_dp)) then
            result%status = 4
            return
         end if
      else if (family == family_poisson) then
         if (any(y < 0.0_dp)) then
            result%status = 5
            return
         end if
      end if

      allocate(par(p + q))
      if (present(coefstart)) then
         if (size(coefstart) /= p) then
            result%status = 6
            return
         end if
         par(1:p) = coefstart
      else
         call fit_vglm(y(q + 1:n), x(q + 1:n, :), family, initfit, &
                       link_id=link_id, weights=w(q + 1:n), max_iter=100)
         if (allocated(initfit%coefficients)) then
            par(1:p) = initfit%coefficients
         else
            par(1:p) = 0.0_dp
         end if
      end if
      if (q > 0) then
         if (present(arstart)) then
            if (size(arstart) /= q) then
               result%status = 7
               return
            end if
            par(p + 1:p + q) = arstart
         else
            par(p + 1:p + q) = 0.0_dp
         end if
      end if

      niter = 300
      if (present(max_iter)) niter = max_iter
      tolerance = 1.0e-8_dp
      if (present(tol)) tolerance = tol
      call bfgs_minimize(objective_fun, par, fval, stat, max_iter=niter, tol=tolerance)
      result%status = stat
      result%converged = stat == 0
      result%coefficients = par(1:p)
      if (q > 0) then
         result%ar = par(p + 1:p + q)
      else
         allocate(result%ar(0))
      end if
      result%link = link_id
      result%ar_order = q
      result%objective = fval
      result%loglik = -fval
      result%aic = 2.0_dp*fval + 2.0_dp*real(p + q, dp)
      call fitted_from_par(par, result%linear_predictor, result%fitted)
      allocate(result%residuals(n - q))
      result%residuals = y(q + 1:n) - result%fitted

      allocate(hess(p + q, p + q))
      call numerical_hessian(objective_fun, par, hess)
      call invert_matrix(hess, cov, stat2)
      if (stat2 == 0) then
         result%covariance = cov
      else
         allocate(result%covariance(p + q, p + q))
         result%covariance = 0.0_dp
      end if

   contains

      real(dp) function objective_fun(theta) result(val)
         real(dp), intent(in) :: theta(:)
         real(dp), allocatable :: eta(:), mu(:)
         real(dp) :: yi, mui, eps
         integer :: t
         call fitted_from_par(theta, eta, mu)
         eps = sqrt(epsilon(1.0_dp))
         val = 0.0_dp
         do t = q + 1, n
            yi = y(t)
            mui = mu(t - q)
            select case (link_id)
            case (link_identity)
               val = val + 0.5_dp*w(t)*(yi - mui)**2
            case (link_log, link_reciprocal)
               mui = max(eps, mui)
               val = val + w(t)*(mui - yi*log(mui))
            case default
               mui = min(1.0_dp - eps, max(eps, mui))
               val = val - w(t)*(yi*log(mui) + (1.0_dp - yi)*log(1.0_dp - mui))
            end select
         end do
         if (.not. finite_scalar(val)) val = huge(1.0_dp)/100.0_dp
      end function objective_fun

      subroutine fitted_from_par(theta, eta, mu)
         real(dp), intent(in) :: theta(:)
         real(dp), allocatable, intent(out) :: eta(:), mu(:)
         real(dp) :: beta(p), phi(q), gy, eps
         integer :: t, lag
         beta = theta(1:p)
         if (q > 0) phi = theta(p + 1:p + q)
         allocate(eta(n - q), mu(n - q))
         eps = sqrt(epsilon(1.0_dp))
         do t = q + 1, n
            eta(t - q) = dot_product(x(t, :), beta)
            do lag = 1, q
               gy = safe_link_y(y(t - lag), link_id, eps)
               eta(t - q) = eta(t - q) + phi(lag)* &
                  (gy - dot_product(x(t - lag, :), beta))
            end do
            mu(t - q) = link_inverse(eta(t - q), link_id)
            call clamp_mu(mu(t - q), link_id, eps)
         end do
      end subroutine fitted_from_par

   end subroutine fit_garma

   subroutine forecast_garma(self, x_history, y_history, x_future, mean)
      class(garma_result_t), intent(in) :: self
      real(dp), intent(in) :: x_history(:, :), y_history(:), x_future(:, :)
      real(dp), allocatable, intent(out) :: mean(:)
      real(dp), allocatable :: xall(:, :), yall(:)
      real(dp) :: eta, gy, eps
      integer :: q, p, h, lag, idx, nh

      q = self%ar_order
      p = size(self%coefficients)
      nh = size(x_future, 1)
      if (size(x_future, 2) /= p .or. size(x_history, 2) /= p .or. &
          size(x_history, 1) < q .or. size(y_history) < q) then
         allocate(mean(0))
         return
      end if
      allocate(mean(nh), xall(q + nh, p), yall(q + nh))
      if (q > 0) then
         xall(1:q, :) = x_history(size(x_history, 1) - q + 1:, :)
         yall(1:q) = y_history(size(y_history) - q + 1:)
      end if
      xall(q + 1:, :) = x_future
      eps = sqrt(epsilon(1.0_dp))
      do h = 1, nh
         idx = q + h
         eta = dot_product(xall(idx, :), self%coefficients)
         do lag = 1, q
            gy = safe_link_y(yall(idx - lag), self%link, eps)
            eta = eta + self%ar(lag)*(gy - &
               dot_product(xall(idx - lag, :), self%coefficients))
         end do
         mean(h) = link_inverse(eta, self%link)
         call clamp_mu(mean(h), self%link, eps)
         yall(idx) = mean(h)
      end do
   end subroutine forecast_garma

   elemental real(dp) function safe_link_y(y, link_id, eps) result(eta)
      real(dp), intent(in) :: y, eps
      integer, intent(in) :: link_id
      real(dp) :: yy
      select case (link_id)
      case (link_log, link_reciprocal)
         yy = max(eps, y)
      case (link_logit, link_probit, link_cloglog, link_cauchit)
         yy = min(1.0_dp - eps, max(eps, y))
      case default
         yy = y
      end select
      eta = link_value(yy, link_id)
   end function safe_link_y

   subroutine clamp_mu(mu, link_id, eps)
      real(dp), intent(inout) :: mu
      integer, intent(in) :: link_id
      real(dp), intent(in) :: eps
      select case (link_id)
      case (link_log, link_reciprocal)
         mu = max(eps, mu)
      case (link_logit, link_probit, link_cloglog, link_cauchit)
         mu = min(1.0_dp - eps, max(eps, mu))
      end select
   end subroutine clamp_mu

   elemental logical function supported_link(link_id) result(ok)
      integer, intent(in) :: link_id
      ok = link_id == link_identity .or. link_id == link_log .or. &
           link_id == link_reciprocal .or. link_id == link_logit .or. &
           link_id == link_probit .or. link_id == link_cloglog .or. &
           link_id == link_cauchit
   end function supported_link

   elemental logical function finite_scalar(x) result(ok)
      real(dp), intent(in) :: x
      ok = abs(x) <= huge(x)
   end function finite_scalar

end module vgam_garma
