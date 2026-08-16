! Copyright (C) 1998-2025 T. W. Yee, University of Auckland.
! Modern Fortran computational translation, 2026.
! SPDX-License-Identifier: GPL-3.0-only
module vgam_count_models
   use vgam_kinds, only : dp
   use vgam_distributions, only : dpois_v, dnbinom_v
   use vgam_optim, only : bfgs_minimize, numerical_hessian
   use vgam_linalg, only : invert_matrix
   implicit none
   private

   type, public :: zero_truncated_poisson_result_t
      real(dp), allocatable :: coefficients(:)
      real(dp), allocatable :: covariance(:, :)
      real(dp), allocatable :: fitted_mean(:)
      real(dp) :: loglik = -huge(1.0_dp)
      real(dp) :: aic = huge(1.0_dp)
      integer :: status = 0
      logical :: converged = .false.
   contains
      procedure :: predict => predict_ztpoisson
   end type zero_truncated_poisson_result_t

   type, public :: hurdle_count_result_t
      real(dp), allocatable :: count_coefficients(:)
      real(dp), allocatable :: zero_coefficients(:)
      real(dp), allocatable :: covariance(:, :)
      real(dp), allocatable :: fitted_mean(:)
      real(dp), allocatable :: zero_probability(:)
      real(dp) :: size = huge(1.0_dp)
      real(dp) :: loglik = -huge(1.0_dp)
      real(dp) :: aic = huge(1.0_dp)
      integer :: status = 0
      logical :: converged = .false.
      logical :: negative_binomial = .false.
   contains
      procedure :: predict => predict_hurdle
   end type hurdle_count_result_t

   type, public :: zinb_result_t
      real(dp), allocatable :: count_coefficients(:)
      real(dp), allocatable :: zero_coefficients(:)
      real(dp), allocatable :: covariance(:, :)
      real(dp), allocatable :: fitted_mean(:)
      real(dp), allocatable :: structural_zero_probability(:)
      real(dp) :: size = 1.0_dp
      real(dp) :: loglik = -huge(1.0_dp)
      real(dp) :: aic = huge(1.0_dp)
      integer :: status = 0
      logical :: converged = .false.
   contains
      procedure :: predict => predict_zinb
   end type zinb_result_t

   public :: dztpois_v, dhurdlepois_v, dhurdlenbinom_v, dzinb_v
   public :: fit_zero_truncated_poisson, fit_hurdle_poisson
   public :: fit_hurdle_negative_binomial, fit_zero_inflated_negative_binomial

contains

   elemental real(dp) function dztpois_v(y, lambda, log_density) result(v)
      integer, intent(in) :: y
      real(dp), intent(in) :: lambda
      logical, intent(in), optional :: log_density
      real(dp) :: ld
      logical :: lg
      lg = .false.
      if (present(log_density)) lg = log_density
      if (lambda <= 0.0_dp) then
         ld = -huge(1.0_dp)
      else if (y <= 0) then
         ld = -huge(1.0_dp)
      else
         ld = dpois_v(y, lambda, .true.) - log1mexp(-lambda)
      end if
      v = merge(ld, exp(ld), lg)
   end function dztpois_v

   elemental real(dp) function dhurdlepois_v(y, lambda, pzero, log_density) result(v)
      integer, intent(in) :: y
      real(dp), intent(in) :: lambda, pzero
      logical, intent(in), optional :: log_density
      real(dp) :: ld
      logical :: lg
      lg = .false.
      if (present(log_density)) lg = log_density
      if (lambda <= 0.0_dp .or. pzero < 0.0_dp .or. pzero > 1.0_dp .or. y < 0) then
         ld = -huge(1.0_dp)
      else if (y == 0) then
         ld = log(max(pzero, tiny(1.0_dp)))
      else
         ld = log(max(1.0_dp - pzero, tiny(1.0_dp))) + dztpois_v(y, lambda, .true.)
      end if
      v = merge(ld, exp(ld), lg)
   end function dhurdlepois_v

   elemental real(dp) function dhurdlenbinom_v(y, mu, size, pzero, log_density) result(v)
      integer, intent(in) :: y
      real(dp), intent(in) :: mu, size, pzero
      logical, intent(in), optional :: log_density
      real(dp) :: ld, prob, p0
      logical :: lg
      lg = .false.
      if (present(log_density)) lg = log_density
      if (mu <= 0.0_dp .or. size <= 0.0_dp .or. pzero < 0.0_dp .or. pzero > 1.0_dp .or. y < 0) then
         ld = -huge(1.0_dp)
      else if (y == 0) then
         ld = log(max(pzero, tiny(1.0_dp)))
      else
         prob = size/(size + mu)
         p0 = exp(dnbinom_v(0, size, prob, .true.))
         ld = log(max(1.0_dp - pzero, tiny(1.0_dp))) + dnbinom_v(y, size, prob, .true.)
         ld = ld - log(max(1.0_dp - p0, tiny(1.0_dp)))
      end if
      v = merge(ld, exp(ld), lg)
   end function dhurdlenbinom_v

   elemental real(dp) function dzinb_v(y, mu, size, pzero, log_density) result(v)
      integer, intent(in) :: y
      real(dp), intent(in) :: mu, size, pzero
      logical, intent(in), optional :: log_density
      real(dp) :: ld, prob, p0, mix
      logical :: lg
      lg = .false.
      if (present(log_density)) lg = log_density
      if (mu <= 0.0_dp .or. size <= 0.0_dp .or. pzero < 0.0_dp .or. pzero > 1.0_dp .or. y < 0) then
         ld = -huge(1.0_dp)
      else
         prob = size/(size + mu)
         if (y == 0) then
            p0 = exp(dnbinom_v(0, size, prob, .true.))
            mix = pzero + (1.0_dp - pzero)*p0
            ld = log(max(mix, tiny(1.0_dp)))
         else
            ld = log(max(1.0_dp - pzero, tiny(1.0_dp))) + dnbinom_v(y, size, prob, .true.)
         end if
      end if
      v = merge(ld, exp(ld), lg)
   end function dzinb_v

   subroutine fit_zero_truncated_poisson(y, x, result, weights, max_iter, tol)
      integer, intent(in) :: y(:)
      real(dp), intent(in) :: x(:, :)
      type(zero_truncated_poisson_result_t), intent(out) :: result
      real(dp), intent(in), optional :: weights(:), tol
      integer, intent(in), optional :: max_iter
      real(dp), allocatable :: w(:), par(:), hess(:, :), cov(:, :), mu(:)
      real(dp) :: fval, tolerance
      integer :: n, p, stat, stat2, niter

      n = size(y)
      p = size(x, 2)
      if (n <= 0 .or. p <= 0 .or. size(x, 1) /= n .or. any(y <= 0)) then
         result%status = 1
         return
      end if
      call get_weights(n, weights, w, stat)
      if (stat /= 0) then
         result%status = 2
         return
      end if
      allocate(par(p))
      par = 0.0_dp
      call initialize_count_intercept(x, real(y, dp), w, par)
      niter = 250
      if (present(max_iter)) niter = max_iter
      tolerance = 1.0e-7_dp
      if (present(tol)) tolerance = tol
      call bfgs_minimize(objective, par, fval, stat, niter, tolerance)
      result%status = stat
      result%converged = stat == 0
      result%coefficients = par
      allocate(mu(n))
      mu = exp(clamp_eta(matmul(x, par)))
      result%fitted_mean = mu/(1.0_dp - exp(-mu))
      result%loglik = -fval
      result%aic = 2.0_dp*fval + 2.0_dp*real(p, dp)
      allocate(hess(p, p))
      call numerical_hessian(objective, par, hess)
      call invert_matrix(hess, cov, stat2)
      if (stat2 == 0) then
         result%covariance = cov
      else
         allocate(result%covariance(p, p))
         result%covariance = 0.0_dp
      end if
   contains
      real(dp) function objective(theta) result(nll)
         real(dp), intent(in) :: theta(:)
         real(dp) :: lambda
         integer :: i
         nll = 0.0_dp
         do i = 1, n
            lambda = exp(clamp_scalar(dot_product(x(i, :), theta)))
            nll = nll - w(i)*dztpois_v(y(i), lambda, .true.)
         end do
      end function objective
   end subroutine fit_zero_truncated_poisson

   subroutine fit_hurdle_poisson(y, x_count, x_zero, result, weights, max_iter, tol)
      integer, intent(in) :: y(:)
      real(dp), intent(in) :: x_count(:, :), x_zero(:, :)
      type(hurdle_count_result_t), intent(out) :: result
      real(dp), intent(in), optional :: weights(:), tol
      integer, intent(in), optional :: max_iter
      call fit_hurdle_common(y, x_count, x_zero, result, .false., weights, max_iter, tol)
   end subroutine fit_hurdle_poisson

   subroutine fit_hurdle_negative_binomial(y, x_count, x_zero, result, weights, max_iter, tol)
      integer, intent(in) :: y(:)
      real(dp), intent(in) :: x_count(:, :), x_zero(:, :)
      type(hurdle_count_result_t), intent(out) :: result
      real(dp), intent(in), optional :: weights(:), tol
      integer, intent(in), optional :: max_iter
      call fit_hurdle_common(y, x_count, x_zero, result, .true., weights, max_iter, tol)
   end subroutine fit_hurdle_negative_binomial

   subroutine fit_hurdle_common(y, xc, xz, result, use_nb, weights, max_iter, tol)
      integer, intent(in) :: y(:)
      real(dp), intent(in) :: xc(:, :), xz(:, :)
      type(hurdle_count_result_t), intent(out) :: result
      logical, intent(in) :: use_nb
      real(dp), intent(in), optional :: weights(:), tol
      integer, intent(in), optional :: max_iter
      real(dp), allocatable :: w(:), par(:), hess(:, :), cov(:, :), jac(:, :), covn(:, :)
      real(dp), allocatable :: mu(:), pzero(:)
      real(dp) :: fval, tolerance, zfrac, sizev
      integer :: n, pc, pz, np, stat, stat2, i, niter

      n = size(y)
      pc = size(xc, 2)
      pz = size(xz, 2)
      if (n <= 0 .or. pc <= 0 .or. pz <= 0 .or. size(xc, 1) /= n .or. &
          size(xz, 1) /= n .or. any(y < 0)) then
         result%status = 1
         return
      end if
      call get_weights(n, weights, w, stat)
      if (stat /= 0) then
         result%status = 2
         return
      end if
      np = pc + pz + merge(1, 0, use_nb)
      allocate(par(np))
      par = 0.0_dp
      call initialize_count_intercept(xc, real(y, dp), w, par(1:pc))
      zfrac = sum(w*merge(1.0_dp, 0.0_dp, y == 0))/max(sum(w), tiny(1.0_dp))
      call initialize_binary_intercept(xz, zfrac, par(pc + 1:pc + pz))
      if (use_nb) par(np) = log(2.0_dp)
      niter = 300
      if (present(max_iter)) niter = max_iter
      tolerance = 1.0e-7_dp
      if (present(tol)) tolerance = tol
      call bfgs_minimize(objective, par, fval, stat, niter, tolerance)
      result%status = stat
      result%converged = stat == 0
      result%negative_binomial = use_nb
      result%count_coefficients = par(1:pc)
      result%zero_coefficients = par(pc + 1:pc + pz)
      sizev = huge(1.0_dp)
      if (use_nb) sizev = exp(min(par(np), 50.0_dp))
      result%size = sizev
      allocate(mu(n), pzero(n))
      mu = exp(clamp_eta(matmul(xc, result%count_coefficients)))
      pzero = logistic_vec(matmul(xz, result%zero_coefficients))
      result%zero_probability = pzero
      if (use_nb) then
         allocate(result%fitted_mean(n))
         do i = 1, n
            result%fitted_mean(i) = (1.0_dp - pzero(i))*truncated_nb_mean(mu(i), sizev)
         end do
      else
         allocate(result%fitted_mean(n))
         result%fitted_mean = (1.0_dp - pzero)*mu/(1.0_dp - exp(-mu))
      end if
      result%loglik = -fval
      result%aic = 2.0_dp*fval + 2.0_dp*real(np, dp)
      allocate(hess(np, np))
      call numerical_hessian(objective, par, hess)
      call invert_matrix(hess, cov, stat2)
      if (stat2 == 0) then
         if (use_nb) then
            allocate(jac(np, np), covn(np, np))
            jac = 0.0_dp
            do i = 1, np - 1
               jac(i, i) = 1.0_dp
            end do
            jac(np, np) = sizev
            covn = matmul(jac, matmul(cov, transpose(jac)))
            result%covariance = covn
         else
            result%covariance = cov
         end if
      else
         allocate(result%covariance(np, np))
         result%covariance = 0.0_dp
      end if
   contains
      real(dp) function objective(theta) result(nll)
         real(dp), intent(in) :: theta(:)
         real(dp) :: muv, pzv, sz
         integer :: ii
         sz = huge(1.0_dp)
         if (use_nb) sz = exp(min(theta(np), 50.0_dp))
         nll = 0.0_dp
         do ii = 1, n
            muv = exp(clamp_scalar(dot_product(xc(ii, :), theta(1:pc))))
            pzv = logistic_scalar(dot_product(xz(ii, :), theta(pc + 1:pc + pz)))
            if (use_nb) then
               nll = nll - w(ii)*dhurdlenbinom_v(y(ii), muv, sz, pzv, .true.)
            else
               nll = nll - w(ii)*dhurdlepois_v(y(ii), muv, pzv, .true.)
            end if
         end do
      end function objective
   end subroutine fit_hurdle_common

   subroutine fit_zero_inflated_negative_binomial(y, x_count, x_zero, result, weights, max_iter, tol)
      integer, intent(in) :: y(:)
      real(dp), intent(in) :: x_count(:, :), x_zero(:, :)
      type(zinb_result_t), intent(out) :: result
      real(dp), intent(in), optional :: weights(:), tol
      integer, intent(in), optional :: max_iter
      real(dp), allocatable :: w(:), par(:), hess(:, :), cov(:, :), jac(:, :), covn(:, :)
      real(dp), allocatable :: mu(:), pzero(:)
      real(dp) :: fval, tolerance, zfrac, sizev
      integer :: n, pc, pz, np, stat, stat2, i, niter

      n = size(y)
      pc = size(x_count, 2)
      pz = size(x_zero, 2)
      np = pc + pz + 1
      if (n <= 0 .or. pc <= 0 .or. pz <= 0 .or. size(x_count, 1) /= n .or. &
          size(x_zero, 1) /= n .or. any(y < 0)) then
         result%status = 1
         return
      end if
      call get_weights(n, weights, w, stat)
      if (stat /= 0) then
         result%status = 2
         return
      end if
      allocate(par(np))
      par = 0.0_dp
      call initialize_count_intercept(x_count, real(y, dp), w, par(1:pc))
      zfrac = sum(w*merge(1.0_dp, 0.0_dp, y == 0))/max(sum(w), tiny(1.0_dp))
      call initialize_binary_intercept(x_zero, max(0.02_dp, 0.5_dp*zfrac), par(pc + 1:pc + pz))
      par(np) = log(2.0_dp)
      niter = 350
      if (present(max_iter)) niter = max_iter
      tolerance = 1.0e-7_dp
      if (present(tol)) tolerance = tol
      call bfgs_minimize(objective, par, fval, stat, niter, tolerance)
      result%status = stat
      result%converged = stat == 0
      result%count_coefficients = par(1:pc)
      result%zero_coefficients = par(pc + 1:pc + pz)
      sizev = exp(min(par(np), 50.0_dp))
      result%size = sizev
      allocate(mu(n), pzero(n), result%fitted_mean(n))
      mu = exp(clamp_eta(matmul(x_count, result%count_coefficients)))
      pzero = logistic_vec(matmul(x_zero, result%zero_coefficients))
      result%structural_zero_probability = pzero
      result%fitted_mean = (1.0_dp - pzero)*mu
      result%loglik = -fval
      result%aic = 2.0_dp*fval + 2.0_dp*real(np, dp)
      allocate(hess(np, np))
      call numerical_hessian(objective, par, hess)
      call invert_matrix(hess, cov, stat2)
      if (stat2 == 0) then
         allocate(jac(np, np), covn(np, np))
         jac = 0.0_dp
         do i = 1, np - 1
            jac(i, i) = 1.0_dp
         end do
         jac(np, np) = sizev
         covn = matmul(jac, matmul(cov, transpose(jac)))
         result%covariance = covn
      else
         allocate(result%covariance(np, np))
         result%covariance = 0.0_dp
      end if
   contains
      real(dp) function objective(theta) result(nll)
         real(dp), intent(in) :: theta(:)
         real(dp) :: muv, pzv, sz
         integer :: ii
         sz = exp(min(theta(np), 50.0_dp))
         nll = 0.0_dp
         do ii = 1, n
            muv = exp(clamp_scalar(dot_product(x_count(ii, :), theta(1:pc))))
            pzv = logistic_scalar(dot_product(x_zero(ii, :), theta(pc + 1:pc + pz)))
            nll = nll - w(ii)*dzinb_v(y(ii), muv, sz, pzv, .true.)
         end do
      end function objective
   end subroutine fit_zero_inflated_negative_binomial

   subroutine predict_ztpoisson(self, x, mean)
      class(zero_truncated_poisson_result_t), intent(in) :: self
      real(dp), intent(in) :: x(:, :)
      real(dp), allocatable, intent(out) :: mean(:)
      real(dp), allocatable :: lambda(:)
      lambda = exp(clamp_eta(matmul(x, self%coefficients)))
      mean = lambda/(1.0_dp - exp(-lambda))
   end subroutine predict_ztpoisson

   subroutine predict_hurdle(self, x_count, x_zero, mean, pzero)
      class(hurdle_count_result_t), intent(in) :: self
      real(dp), intent(in) :: x_count(:, :), x_zero(:, :)
      real(dp), allocatable, intent(out) :: mean(:), pzero(:)
      real(dp), allocatable :: mu(:)
      integer :: i
      mu = exp(clamp_eta(matmul(x_count, self%count_coefficients)))
      pzero = logistic_vec(matmul(x_zero, self%zero_coefficients))
      allocate(mean(size(mu)))
      if (self%negative_binomial) then
         do i = 1, size(mu)
            mean(i) = (1.0_dp - pzero(i))*truncated_nb_mean(mu(i), self%size)
         end do
      else
         mean = (1.0_dp - pzero)*mu/(1.0_dp - exp(-mu))
      end if
   end subroutine predict_hurdle

   subroutine predict_zinb(self, x_count, x_zero, mean, structural_zero)
      class(zinb_result_t), intent(in) :: self
      real(dp), intent(in) :: x_count(:, :), x_zero(:, :)
      real(dp), allocatable, intent(out) :: mean(:), structural_zero(:)
      real(dp), allocatable :: mu(:)
      mu = exp(clamp_eta(matmul(x_count, self%count_coefficients)))
      structural_zero = logistic_vec(matmul(x_zero, self%zero_coefficients))
      mean = (1.0_dp - structural_zero)*mu
   end subroutine predict_zinb

   subroutine get_weights(n, input, w, status)
      integer, intent(in) :: n
      real(dp), intent(in), optional :: input(:)
      real(dp), allocatable, intent(out) :: w(:)
      integer, intent(out) :: status
      status = 0
      allocate(w(n))
      w = 1.0_dp
      if (present(input)) then
         if (size(input) /= n .or. any(input < 0.0_dp)) then
            status = 1
            return
         end if
         w = input
      end if
   end subroutine get_weights

   subroutine initialize_count_intercept(x, y, w, beta)
      real(dp), intent(in) :: x(:, :), y(:), w(:)
      real(dp), intent(inout) :: beta(:)
      real(dp) :: mu
      beta = 0.0_dp
      if (size(beta) == 0) return
      if (all(abs(x(:, 1) - 1.0_dp) <= 100.0_dp*epsilon(1.0_dp))) then
         mu = max(sum(w*y)/max(sum(w), tiny(1.0_dp)), 0.1_dp)
         beta(1) = log(mu)
      end if
   end subroutine initialize_count_intercept

   subroutine initialize_binary_intercept(x, p, beta)
      real(dp), intent(in) :: x(:, :), p
      real(dp), intent(inout) :: beta(:)
      real(dp) :: pp
      beta = 0.0_dp
      if (size(beta) == 0) return
      if (all(abs(x(:, 1) - 1.0_dp) <= 100.0_dp*epsilon(1.0_dp))) then
         pp = min(0.98_dp, max(0.02_dp, p))
         beta(1) = log(pp/(1.0_dp - pp))
      end if
   end subroutine initialize_binary_intercept

   elemental real(dp) function truncated_nb_mean(mu, size) result(v)
      real(dp), intent(in) :: mu, size
      real(dp) :: prob, p0
      prob = size/(size + mu)
      p0 = exp(dnbinom_v(0, size, prob, .true.))
      v = mu/max(1.0_dp - p0, tiny(1.0_dp))
   end function truncated_nb_mean

   elemental real(dp) function logistic_scalar(x) result(v)
      real(dp), intent(in) :: x
      if (x >= 0.0_dp) then
         v = 1.0_dp/(1.0_dp + exp(-min(x, 700.0_dp)))
      else
         v = exp(max(x, -700.0_dp))/(1.0_dp + exp(max(x, -700.0_dp)))
      end if
   end function logistic_scalar

   pure function logistic_vec(x) result(v)
      real(dp), intent(in) :: x(:)
      real(dp) :: v(size(x))
      integer :: i
      do i = 1, size(x)
         v(i) = logistic_scalar(x(i))
      end do
   end function logistic_vec

   elemental real(dp) function clamp_scalar(x) result(v)
      real(dp), intent(in) :: x
      v = min(40.0_dp, max(-40.0_dp, x))
   end function clamp_scalar

   pure function clamp_eta(x) result(v)
      real(dp), intent(in) :: x(:)
      real(dp) :: v(size(x))
      v = min(40.0_dp, max(-40.0_dp, x))
   end function clamp_eta

   elemental real(dp) function log1mexp(a) result(v)
      real(dp), intent(in) :: a
      if (a >= 0.0_dp) then
         v = -huge(1.0_dp)
      else if (a < log(0.5_dp)) then
         v = log(1.0_dp - exp(a))
      else
         v = log(-expm1_local(a))
      end if
   end function log1mexp

   elemental real(dp) function expm1_local(x) result(v)
      real(dp), intent(in) :: x
      if (abs(x) < 1.0e-5_dp) then
         v = x*(1.0_dp + x*(0.5_dp + x*(1.0_dp/6.0_dp + x/24.0_dp)))
      else
         v = exp(x) - 1.0_dp
      end if
   end function expm1_local

end module vgam_count_models
