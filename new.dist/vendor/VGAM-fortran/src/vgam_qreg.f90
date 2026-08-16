! Copyright (C) 1998-2025 T. W. Yee, University of Auckland.
! Modern Fortran computational translation, 2026.
! SPDX-License-Identifier: GPL-3.0-only
module vgam_qreg
   use vgam_kinds, only : dp, pi
   use vgam_special, only : normal_quantile, log1p_v, expm1_v
   use vgam_optim, only : bfgs_minimize, numerical_hessian
   use vgam_linalg, only : invert_matrix
   implicit none
   private

   type, public :: yj_normal_result_t
      real(dp) :: lambda = 1.0_dp
      real(dp) :: sigma = 1.0_dp
      real(dp), allocatable :: coefficients(:)
      real(dp), allocatable :: covariance(:, :)
      real(dp), allocatable :: transformed(:)
      real(dp), allocatable :: fitted_transformed(:)
      real(dp) :: loglik = -huge(1.0_dp)
      real(dp) :: aic = huge(1.0_dp)
      integer :: status = 0
      logical :: converged = .false.
   contains
      procedure :: predict_quantile => predict_yj_quantile
   end type yj_normal_result_t

   type, public :: lms_yj_result_t
      real(dp), allocatable :: lambda_coefficients(:)
      real(dp), allocatable :: mu_coefficients(:)
      real(dp), allocatable :: log_sigma_coefficients(:)
      real(dp), allocatable :: covariance(:, :)
      real(dp), allocatable :: lambda(:)
      real(dp), allocatable :: mu(:)
      real(dp), allocatable :: sigma(:)
      real(dp), allocatable :: transformed(:)
      real(dp) :: loglik = -huge(1.0_dp)
      real(dp) :: aic = huge(1.0_dp)
      integer :: status = 0
      logical :: converged = .false.
   contains
      procedure :: predict_quantile => predict_lms_yj_quantile
   end type lms_yj_result_t

   public :: yeo_johnson, yeo_johnson_inverse, yeo_johnson_lambda_derivative
   public :: dyj_dy, fit_yj_normal, fit_lms_yj

contains

   recursive elemental real(dp) function yeo_johnson(y, lambda, derivative) result(ans)
      real(dp), intent(in) :: y, lambda
      integer, intent(in), optional :: derivative
      integer :: d
      real(dp) :: psi, eps
      d = 0
      if (present(derivative)) d = derivative
      eps = sqrt(epsilon(1.0_dp))
      if (d < 0) then
         ans = huge(1.0_dp)
         return
      end if
      if (d == 0) then
         if (y >= 0.0_dp) then
            if (abs(lambda) > eps) then
               ans = expm1_v(lambda*log1p_v(y))/lambda
            else
               ans = log1p_v(y)
            end if
         else
            if (abs(lambda - 2.0_dp) > eps) then
               ans = -expm1_v((2.0_dp - lambda)*log1p_v(-y))/(2.0_dp - lambda)
            else
               ans = -log1p_v(-y)
            end if
         end if
         return
      end if

      psi = yeo_johnson(y, lambda, d - 1)
      if (y >= 0.0_dp) then
         if (abs(lambda) > eps) then
            ans = ((1.0_dp + y)**lambda*log1p_v(y)**d - real(d, dp)*psi)/lambda
         else
            ans = log1p_v(y)**(d + 1)/real(d + 1, dp)
         end if
      else
         if (abs(lambda - 2.0_dp) > eps) then
            ans = -((1.0_dp - y)**(2.0_dp - lambda)*(-log1p_v(-y))**d - &
                    real(d, dp)*psi)/(2.0_dp - lambda)
         else
            ans = (-log1p_v(-y))**(d + 1)/real(d + 1, dp)
         end if
      end if
   end function yeo_johnson

   elemental real(dp) function yeo_johnson_inverse(psi, lambda) result(y)
      real(dp), intent(in) :: psi, lambda
      real(dp) :: eps
      eps = sqrt(epsilon(1.0_dp))
      if (psi >= 0.0_dp) then
         if (abs(lambda) > eps) then
            y = exp(log(1.0_dp + lambda*psi)/lambda) - 1.0_dp
         else
            y = expm1_v(psi)
         end if
      else
         if (abs(lambda - 2.0_dp) > eps) then
            y = 1.0_dp - exp(log(1.0_dp - (2.0_dp - lambda)*psi)/(2.0_dp - lambda))
         else
            y = -expm1_v(-psi)
         end if
      end if
   end function yeo_johnson_inverse

   elemental real(dp) function yeo_johnson_lambda_derivative(y, lambda) result(v)
      real(dp), intent(in) :: y, lambda
      v = yeo_johnson(y, lambda, 1)
   end function yeo_johnson_lambda_derivative

   elemental real(dp) function dyj_dy(y, lambda) result(v)
      real(dp), intent(in) :: y, lambda
      if (y >= 0.0_dp) then
         v = exp((lambda - 1.0_dp)*log1p_v(y))
      else
         v = exp((1.0_dp - lambda)*log1p_v(-y))
      end if
   end function dyj_dy

   subroutine fit_yj_normal(y, x, result, weights, lambda_start, max_iter, tol)
      real(dp), intent(in) :: y(:), x(:, :)
      type(yj_normal_result_t), intent(out) :: result
      real(dp), intent(in), optional :: weights(:), lambda_start, tol
      integer, intent(in), optional :: max_iter
      real(dp), allocatable :: w(:), par(:), hess(:, :), cov(:, :), jac(:, :), covn(:, :)
      real(dp), allocatable :: psi(:)
      real(dp) :: fval, tolerance, mu0, sd0, lam
      integer :: n, p, np, stat, stat2, i, niter

      n = size(y)
      p = size(x, 2)
      np = p + 2
      if (n <= 0 .or. p <= 0 .or. size(x, 1) /= n) then
         result%status = 1
         return
      end if
      allocate(w(n))
      w = 1.0_dp
      if (present(weights)) then
         if (size(weights) /= n .or. any(weights < 0.0_dp)) then
            result%status = 2
            return
         end if
         w = weights
      end if
      allocate(par(np), psi(n))
      par = 0.0_dp
      par(1) = 1.0_dp
      if (present(lambda_start)) par(1) = lambda_start
      do i = 1, n
         psi(i) = yeo_johnson(y(i), par(1))
      end do
      mu0 = sum(w*psi)/max(sum(w), tiny(1.0_dp))
      sd0 = sqrt(max(sum(w*(psi - mu0)**2)/max(sum(w), tiny(1.0_dp)), 1.0e-6_dp))
      if (all(abs(x(:, 1) - 1.0_dp) <= 100.0_dp*epsilon(1.0_dp))) par(2) = mu0
      par(np) = log(sd0)
      niter = 300
      if (present(max_iter)) niter = max_iter
      tolerance = 1.0e-7_dp
      if (present(tol)) tolerance = tol
      call bfgs_minimize(objective, par, fval, stat, niter, tolerance)
      result%status = stat
      result%converged = stat == 0
      lam = par(1)
      result%lambda = lam
      result%coefficients = par(2:p + 1)
      result%sigma = exp(min(par(np), 50.0_dp))
      allocate(result%transformed(n), result%fitted_transformed(n))
      do i = 1, n
         result%transformed(i) = yeo_johnson(y(i), lam)
      end do
      result%fitted_transformed = matmul(x, result%coefficients)
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
         jac(np, np) = result%sigma
         covn = matmul(jac, matmul(cov, transpose(jac)))
         result%covariance = covn
      else
         allocate(result%covariance(np, np))
         result%covariance = 0.0_dp
      end if
   contains
      real(dp) function objective(theta) result(nll)
         real(dp), intent(in) :: theta(:)
         real(dp) :: sigma, z, tr, logjac
         integer :: ii
         sigma = exp(min(theta(np), 50.0_dp))
         nll = 0.0_dp
         do ii = 1, n
            tr = yeo_johnson(y(ii), theta(1))
            z = (tr - dot_product(x(ii, :), theta(2:p + 1)))/sigma
            if (y(ii) >= 0.0_dp) then
               logjac = (theta(1) - 1.0_dp)*log1p_v(y(ii))
            else
               logjac = (1.0_dp - theta(1))*log1p_v(-y(ii))
            end if
            nll = nll + w(ii)*(0.5_dp*z*z + log(sigma) + 0.5_dp*log(2.0_dp*pi) - logjac)
         end do
      end function objective
   end subroutine fit_yj_normal

   subroutine predict_yj_quantile(self, x, probability, quantile)
      class(yj_normal_result_t), intent(in) :: self
      real(dp), intent(in) :: x(:, :), probability
      real(dp), allocatable, intent(out) :: quantile(:)
      real(dp), allocatable :: eta(:)
      real(dp) :: z
      integer :: i
      if (probability <= 0.0_dp .or. probability >= 1.0_dp) then
         allocate(quantile(0))
         return
      end if
      eta = matmul(x, self%coefficients)
      z = normal_quantile(probability)
      allocate(quantile(size(eta)))
      do i = 1, size(eta)
         quantile(i) = yeo_johnson_inverse(eta(i) + self%sigma*z, self%lambda)
      end do
   end subroutine predict_yj_quantile


   subroutine fit_lms_yj(y, x_lambda, x_mu, x_sigma, result, weights, &
                         max_iter, tol)
      real(dp), intent(in) :: y(:), x_lambda(:, :), x_mu(:, :), x_sigma(:, :)
      type(lms_yj_result_t), intent(out) :: result
      real(dp), intent(in), optional :: weights(:), tol
      integer, intent(in), optional :: max_iter
      real(dp), allocatable :: w(:), par(:), hess(:, :), cov(:, :)
      real(dp) :: fval, tolerance, psi0, sd0
      integer :: n, pl, pm, ps, np, il, im, is, stat, stat2, niter, i

      n = size(y)
      pl = size(x_lambda, 2)
      pm = size(x_mu, 2)
      ps = size(x_sigma, 2)
      if (n <= 0 .or. pl <= 0 .or. pm <= 0 .or. ps <= 0 .or. &
          size(x_lambda, 1) /= n .or. size(x_mu, 1) /= n .or. &
          size(x_sigma, 1) /= n) then
         result%status = 1
         return
      end if
      allocate(w(n))
      w = 1.0_dp
      if (present(weights)) then
         if (size(weights) /= n .or. any(weights < 0.0_dp)) then
            result%status = 2
            return
         end if
         w = weights
      end if
      np = pl + pm + ps
      il = 1
      im = pl + 1
      is = pl + pm + 1
      allocate(par(np))
      par = 0.0_dp
      if (all(abs(x_lambda(:, 1) - 1.0_dp) <= 100.0_dp*epsilon(1.0_dp))) then
         par(il) = 1.0_dp
      end if
      psi0 = sum([(yeo_johnson(y(i), 1.0_dp), i=1,n)])/real(n, dp)
      if (all(abs(x_mu(:, 1) - 1.0_dp) <= 100.0_dp*epsilon(1.0_dp))) then
         par(im) = psi0
      end if
      sd0 = sqrt(max(sum([(yeo_johnson(y(i), 1.0_dp) - psi0, i=1,n)]**2)/ &
                     real(n, dp), 1.0e-6_dp))
      if (all(abs(x_sigma(:, 1) - 1.0_dp) <= 100.0_dp*epsilon(1.0_dp))) then
         par(is) = log(sd0)
      end if
      niter = 400
      if (present(max_iter)) niter = max_iter
      tolerance = 1.0e-7_dp
      if (present(tol)) tolerance = tol
      call bfgs_minimize(objective, par, fval, stat, max_iter=niter, tol=tolerance)
      result%status = stat
      result%converged = stat == 0
      result%lambda_coefficients = par(1:pl)
      result%mu_coefficients = par(pl + 1:pl + pm)
      result%log_sigma_coefficients = par(pl + pm + 1:np)
      result%lambda = matmul(x_lambda, result%lambda_coefficients)
      result%mu = matmul(x_mu, result%mu_coefficients)
      result%sigma = exp(min(matmul(x_sigma, result%log_sigma_coefficients), 50.0_dp))
      allocate(result%transformed(n))
      do i = 1, n
         result%transformed(i) = yeo_johnson(y(i), result%lambda(i))
      end do
      result%loglik = -fval
      result%aic = 2.0_dp*fval + 2.0_dp*real(np, dp)
      allocate(hess(np, np))
      call numerical_hessian(objective, par, hess)
      call invert_matrix(hess, cov, stat2)
      if (stat2 == 0) then
         result%covariance = cov
      else
         allocate(result%covariance(np, np))
         result%covariance = 0.0_dp
      end if

   contains
      real(dp) function objective(theta) result(nll)
         real(dp), intent(in) :: theta(:)
         real(dp) :: lam, muv, sig, tr, z, logjac
         integer :: ii
         nll = 0.0_dp
         do ii = 1, n
            lam = dot_product(x_lambda(ii, :), theta(1:pl))
            muv = dot_product(x_mu(ii, :), theta(pl + 1:pl + pm))
            sig = exp(min(dot_product(x_sigma(ii, :), &
                      theta(pl + pm + 1:np)), 50.0_dp))
            tr = yeo_johnson(y(ii), lam)
            z = (tr - muv)/sig
            if (y(ii) >= 0.0_dp) then
               logjac = (lam - 1.0_dp)*log1p_v(y(ii))
            else
               logjac = (1.0_dp - lam)*log1p_v(-y(ii))
            end if
            nll = nll + w(ii)*(0.5_dp*z*z + log(sig) + &
                  0.5_dp*log(2.0_dp*pi) - logjac)
         end do
      end function objective
   end subroutine fit_lms_yj

   subroutine predict_lms_yj_quantile(self, x_lambda, x_mu, x_sigma, &
                                      probability, quantile)
      class(lms_yj_result_t), intent(in) :: self
      real(dp), intent(in) :: x_lambda(:, :), x_mu(:, :), x_sigma(:, :)
      real(dp), intent(in) :: probability
      real(dp), allocatable, intent(out) :: quantile(:)
      real(dp), allocatable :: lam(:), muv(:), sig(:)
      real(dp) :: z
      integer :: i, n
      n = size(x_mu, 1)
      if (probability <= 0.0_dp .or. probability >= 1.0_dp .or. &
          size(x_lambda, 1) /= n .or. size(x_sigma, 1) /= n .or. &
          size(x_lambda, 2) /= size(self%lambda_coefficients) .or. &
          size(x_mu, 2) /= size(self%mu_coefficients) .or. &
          size(x_sigma, 2) /= size(self%log_sigma_coefficients)) then
         allocate(quantile(0))
         return
      end if
      lam = matmul(x_lambda, self%lambda_coefficients)
      muv = matmul(x_mu, self%mu_coefficients)
      sig = exp(min(matmul(x_sigma, self%log_sigma_coefficients), 50.0_dp))
      z = normal_quantile(probability)
      allocate(quantile(n))
      do i = 1, n
         quantile(i) = yeo_johnson_inverse(muv(i) + sig(i)*z, lam(i))
      end do
   end subroutine predict_lms_yj_quantile

end module vgam_qreg
