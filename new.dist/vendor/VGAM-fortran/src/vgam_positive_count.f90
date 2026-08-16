! Copyright (C) 1998-2025 T. W. Yee, University of Auckland.
! Modern Fortran computational translation, 2026.
! SPDX-License-Identifier: GPL-3.0-only
module vgam_positive_count
   use vgam_kinds, only : dp
   use vgam_distributions, only : dpois_v, ppois_v, qpois_v, rpois_v, &
      dnbinom_v, pnbinom_v, qnbinom_v, rnbinom_v
   use vgam_optim, only : bfgs_minimize, numerical_hessian
   use vgam_linalg, only : invert_matrix
   implicit none
   private

   type, public :: positive_nb_result_t
      real(dp), allocatable :: coefficients(:)
      real(dp), allocatable :: covariance(:, :)
      real(dp), allocatable :: fitted_parent_mean(:)
      real(dp), allocatable :: fitted_truncated_mean(:)
      real(dp) :: size = 1.0_dp
      real(dp) :: loglik = -huge(1.0_dp)
      real(dp) :: aic = huge(1.0_dp)
      integer :: status = 0
      logical :: converged = .false.
   contains
      procedure :: predict => predict_positive_nb
   end type positive_nb_result_t

   public :: dpospois_v, ppospois_v, qpospois_v, rpospois_v
   public :: dposnbinom_v, pposnbinom_v, qposnbinom_v, rposnbinom_v
   public :: fit_positive_negative_binomial

contains

   elemental real(dp) function dpospois_v(y, lambda, log_density) result(v)
      integer, intent(in) :: y
      real(dp), intent(in) :: lambda
      logical, intent(in), optional :: log_density
      real(dp) :: ld, p0
      logical :: lg
      lg = .false.; if (present(log_density)) lg = log_density
      if (lambda <= 0.0_dp .or. y < 1) then
         ld = -huge(1.0_dp)
      else
         p0 = exp(-lambda)
         ld = dpois_v(y, lambda, .true.) - log(max(1.0_dp - p0, tiny(1.0_dp)))
      end if
      if (lg) then
         v = ld
      else if (ld <= -700.0_dp) then
         v = 0.0_dp
      else
         v = exp(ld)
      end if
   end function dpospois_v

   elemental real(dp) function ppospois_v(q, lambda) result(p)
      integer, intent(in) :: q
      real(dp), intent(in) :: lambda
      real(dp) :: p0
      if (lambda <= 0.0_dp) then
         p = -1.0_dp
      else if (q < 1) then
         p = 0.0_dp
      else
         p0 = exp(-lambda)
         p = (ppois_v(q, lambda) - p0)/max(1.0_dp - p0, tiny(1.0_dp))
         p = min(1.0_dp, max(0.0_dp, p))
      end if
   end function ppospois_v

   integer function qpospois_v(p, lambda) result(q)
      real(dp), intent(in) :: p, lambda
      real(dp) :: p0, target
      if (p < 0.0_dp .or. p > 1.0_dp .or. lambda <= 0.0_dp) then
         q = -1
         return
      end if
      if (p <= 0.0_dp) then; q = 1; return; end if
      p0 = exp(-lambda); target = p0 + p*(1.0_dp - p0)
      q = max(1, qpois_v(target, lambda))
   end function qpospois_v

   integer function rpospois_v(lambda) result(y)
      real(dp), intent(in) :: lambda
      if (lambda <= 0.0_dp) then; y = -1; return; end if
      do
         y = rpois_v(lambda)
         if (y > 0) exit
      end do
   end function rpospois_v

   elemental real(dp) function dposnbinom_v(y, mu, size, log_density) result(v)
      integer, intent(in) :: y
      real(dp), intent(in) :: mu, size
      logical, intent(in), optional :: log_density
      real(dp) :: ld, prob, p0
      logical :: lg
      lg = .false.; if (present(log_density)) lg = log_density
      if (mu <= 0.0_dp .or. size <= 0.0_dp .or. y < 1) then
         ld = -huge(1.0_dp)
      else
         prob = size/(size + mu)
         p0 = exp(dnbinom_v(0, size, prob, .true.))
         ld = dnbinom_v(y, size, prob, .true.) - log(max(1.0_dp - p0, tiny(1.0_dp)))
      end if
      if (lg) then
         v = ld
      else if (ld <= -700.0_dp) then
         v = 0.0_dp
      else
         v = exp(ld)
      end if
   end function dposnbinom_v

   elemental real(dp) function pposnbinom_v(q, mu, size) result(p)
      integer, intent(in) :: q
      real(dp), intent(in) :: mu, size
      real(dp) :: prob, p0
      if (mu <= 0.0_dp .or. size <= 0.0_dp) then
         p = -1.0_dp
      else if (q < 1) then
         p = 0.0_dp
      else
         prob = size/(size + mu)
         p0 = exp(dnbinom_v(0, size, prob, .true.))
         p = (pnbinom_v(q, size, prob) - p0)/max(1.0_dp - p0, tiny(1.0_dp))
         p = min(1.0_dp, max(0.0_dp, p))
      end if
   end function pposnbinom_v

   integer function qposnbinom_v(p, mu, size) result(q)
      real(dp), intent(in) :: p, mu, size
      real(dp) :: prob, p0, target
      if (p < 0.0_dp .or. p > 1.0_dp .or. mu <= 0.0_dp .or. size <= 0.0_dp) then
         q = -1
         return
      end if
      if (p <= 0.0_dp) then; q = 1; return; end if
      prob = size/(size + mu); p0 = exp(dnbinom_v(0, size, prob, .true.))
      target = p0 + p*(1.0_dp - p0)
      q = max(1, qnbinom_v(target, size, prob))
   end function qposnbinom_v

   integer function rposnbinom_v(mu, size) result(y)
      real(dp), intent(in) :: mu, size
      real(dp) :: prob
      if (mu <= 0.0_dp .or. size <= 0.0_dp) then; y = -1; return; end if
      prob = size/(size + mu)
      do
         y = rnbinom_v(size, prob)
         if (y > 0) exit
      end do
   end function rposnbinom_v

   subroutine fit_positive_negative_binomial(y, x, result, weights, max_iter, tol)
      integer, intent(in) :: y(:)
      real(dp), intent(in) :: x(:, :)
      type(positive_nb_result_t), intent(out) :: result
      real(dp), intent(in), optional :: weights(:), tol
      integer, intent(in), optional :: max_iter
      real(dp), allocatable :: w(:), par(:), h(:, :), cov(:, :)
      real(dp) :: fval, tolerance, mean_y, sizev, mu, p0
      integer :: n, p, np, i, stat, stat2, niter
      n = size(y); p = size(x, 2); np = p + 1
      if (n < 3 .or. p < 1 .or. size(x, 1) /= n .or. any(y < 1)) then
         result%status = 1; return
      end if
      allocate(w(n)); w = 1.0_dp
      if (present(weights)) then
         if (size(weights) /= n .or. any(weights < 0.0_dp)) then; result%status = 2; return; end if
         w = weights
      end if
      allocate(par(np)); par = 0.0_dp
      mean_y = sum(w*real(y, dp))/max(sum(w), tiny(1.0_dp))
      if (is_intercept_design(x)) par(1) = log(max(0.2_dp, mean_y - 0.5_dp))
      par(np) = log(2.0_dp)
      niter = 350; if (present(max_iter)) niter = max_iter
      tolerance = 1.0e-7_dp; if (present(tol)) tolerance = tol
      call bfgs_minimize(objective, par, fval, stat, max_iter=niter, tol=tolerance)
      result%coefficients = par(1:p); result%size = exp(min(30.0_dp, par(np)))
      result%status = stat; result%converged = stat == 0; result%loglik = -fval
      result%aic = 2.0_dp*fval + 2.0_dp*real(np, dp)
      allocate(result%fitted_parent_mean(n), result%fitted_truncated_mean(n))
      sizev = result%size
      do i = 1, n
         mu = exp(clamp_eta(dot_product(x(i, :), result%coefficients)))
         p0 = (sizev/(sizev + mu))**sizev
         result%fitted_parent_mean(i) = mu
         result%fitted_truncated_mean(i) = mu/max(1.0_dp - p0, tiny(1.0_dp))
      end do
      allocate(h(np, np)); call numerical_hessian(objective, par, h); call invert_matrix(h, cov, stat2)
      if (stat2 == 0) then; result%covariance = cov; else; allocate(result%covariance(0, 0)); end if
   contains
      real(dp) function objective(theta) result(nll)
         real(dp), intent(in) :: theta(:)
         real(dp) :: muv, sv, ld
         integer :: row
         sv = exp(min(30.0_dp, max(-20.0_dp, theta(np)))); nll = 0.0_dp
         do row = 1, n
            muv = exp(clamp_eta(dot_product(x(row, :), theta(1:p))))
            ld = dposnbinom_v(y(row), muv, sv, .true.)
            if (.not. finite_scalar(ld)) then; nll = huge(1.0_dp)/100.0_dp; return; end if
            nll = nll - w(row)*ld
         end do
      end function objective
   end subroutine fit_positive_negative_binomial

   subroutine predict_positive_nb(self, x, parent_mean, truncated_mean)
      class(positive_nb_result_t), intent(in) :: self
      real(dp), intent(in) :: x(:, :)
      real(dp), allocatable, intent(out) :: parent_mean(:), truncated_mean(:)
      real(dp) :: p0
      integer :: i, n
      n = size(x, 1)
      if (size(x, 2) /= size(self%coefficients)) then
         allocate(parent_mean(0), truncated_mean(0)); return
      end if
      allocate(parent_mean(n), truncated_mean(n))
      do i = 1, n
         parent_mean(i) = exp(clamp_eta(dot_product(x(i, :), self%coefficients)))
         p0 = (self%size/(self%size + parent_mean(i)))**self%size
         truncated_mean(i) = parent_mean(i)/max(1.0_dp - p0, tiny(1.0_dp))
      end do
   end subroutine predict_positive_nb

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
end module vgam_positive_count
