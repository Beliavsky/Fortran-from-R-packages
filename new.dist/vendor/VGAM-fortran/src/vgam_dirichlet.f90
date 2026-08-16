! Copyright (C) 1998-2025 T. W. Yee, University of Auckland.
! Modern Fortran computational translation, 2026.
! SPDX-License-Identifier: GPL-3.0-only
module vgam_dirichlet
   use vgam_kinds, only : dp
   use vgam_distributions, only : rgamma_v
   use vgam_special, only : trigamma
   use vgam_optim, only : bfgs_minimize, numerical_hessian
   use vgam_linalg, only : invert_matrix
   implicit none
   private

   type, public :: dirichlet_regression_result_t
      real(dp), allocatable :: coefficients(:, :)
      real(dp), allocatable :: covariance(:, :)
      real(dp), allocatable :: fitted_shape(:, :)
      real(dp), allocatable :: fitted_mean(:, :)
      real(dp) :: loglik = -huge(1.0_dp)
      real(dp) :: aic = huge(1.0_dp)
      integer :: status = 0
      logical :: converged = .false.
      logical :: parallel_slopes = .false.
   contains
      procedure :: predict => predict_dirichlet_regression
   end type dirichlet_regression_result_t

   public :: dirichlet_pdf, dirichlet_logpdf, random_dirichlet
   public :: dirichlet_eim_shape, dirichlet_eim_logshape
   public :: fit_dirichlet_regression

contains

   pure real(dp) function dirichlet_logpdf(y, alpha) result(ld)
      real(dp), intent(in) :: y(:), alpha(:)
      integer :: j
      if (size(y) < 2 .or. size(alpha) /= size(y) .or. any(y <= 0.0_dp) .or. &
          any(alpha <= 0.0_dp) .or. abs(sum(y) - 1.0_dp) > 1.0e-8_dp) then
         ld = -huge(1.0_dp)
         return
      end if
      ld = log_gamma(sum(alpha))
      do j = 1, size(alpha)
         ld = ld - log_gamma(alpha(j)) + (alpha(j) - 1.0_dp)*log(y(j))
      end do
   end function dirichlet_logpdf

   pure real(dp) function dirichlet_pdf(y, alpha) result(d)
      real(dp), intent(in) :: y(:), alpha(:)
      real(dp) :: ld
      ld = dirichlet_logpdf(y, alpha)
      if (ld <= -700.0_dp) then
         d = 0.0_dp
      else
         d = exp(ld)
      end if
   end function dirichlet_pdf

   subroutine random_dirichlet(alpha, y, status)
      real(dp), intent(in) :: alpha(:)
      real(dp), intent(out) :: y(:)
      integer, intent(out), optional :: status
      real(dp) :: s
      integer :: j
      if (size(y) /= size(alpha) .or. size(alpha) < 2 .or. any(alpha <= 0.0_dp)) then
         if (size(y) > 0) y = 0.0_dp
         if (present(status)) status = 1
         return
      end if
      do j = 1, size(alpha)
         y(j) = rgamma_v(alpha(j), 1.0_dp)
      end do
      s = sum(y)
      if (s <= 0.0_dp) then
         y = 1.0_dp/real(size(y), dp)
         if (present(status)) status = 2
      else
         y = y/s
         if (present(status)) status = 0
      end if
   end subroutine random_dirichlet

   subroutine dirichlet_eim_shape(alpha, info)
      real(dp), intent(in) :: alpha(:)
      real(dp), allocatable, intent(out) :: info(:, :)
      real(dp) :: common
      integer :: j, m
      m = size(alpha)
      if (m < 2 .or. any(alpha <= 0.0_dp)) then
         allocate(info(0, 0))
         return
      end if
      allocate(info(m, m))
      common = trigamma(sum(alpha))
      info = -common
      do j = 1, m
         info(j, j) = info(j, j) + trigamma(alpha(j))
      end do
   end subroutine dirichlet_eim_shape

   subroutine dirichlet_eim_logshape(alpha, info)
      real(dp), intent(in) :: alpha(:)
      real(dp), allocatable, intent(out) :: info(:, :)
      real(dp), allocatable :: raw(:, :)
      integer :: j, k, m
      call dirichlet_eim_shape(alpha, raw)
      m = size(alpha)
      if (size(raw, 1) == 0) then
         allocate(info(0, 0))
         return
      end if
      allocate(info(m, m))
      do k = 1, m
         do j = 1, m
            info(j, k) = alpha(j)*raw(j, k)*alpha(k)
         end do
      end do
   end subroutine dirichlet_eim_logshape

   subroutine fit_dirichlet_regression(y, x, result, weights, parallel_slopes, max_iter, tol)
      real(dp), intent(in) :: y(:, :), x(:, :)
      type(dirichlet_regression_result_t), intent(out) :: result
      real(dp), intent(in), optional :: weights(:), tol
      logical, intent(in), optional :: parallel_slopes
      integer, intent(in), optional :: max_iter
      real(dp), allocatable :: w(:), par(:), h(:, :), covred(:, :), map(:, :)
      real(dp), allocatable :: alpha(:, :), means(:, :), beta(:, :), covfull(:, :)
      real(dp) :: fval, tolerance
      logical :: parallel
      integer :: n, p, m, np, i, j, stat, stat2, niter

      n = size(y, 1); m = size(y, 2); p = size(x, 2)
      if (n < 2 .or. m < 2 .or. p < 1 .or. size(x, 1) /= n .or. any(y <= 0.0_dp)) then
         result%status = 1
         return
      end if
      do i = 1, n
         if (abs(sum(y(i, :)) - 1.0_dp) > 1.0e-7_dp) then
            result%status = 2
            return
         end if
      end do
      allocate(w(n)); w = 1.0_dp
      if (present(weights)) then
         if (size(weights) /= n .or. any(weights < 0.0_dp)) then
            result%status = 3
            return
         end if
         w = weights
      end if
      parallel = .false.; if (present(parallel_slopes)) parallel = parallel_slopes
      if (parallel .and. p > 1) then
         np = m + p - 1
      else
         np = p*m
      end if
      allocate(par(np)); par = 0.0_dp
      call initialize_parameters(y, x, w, parallel, par)
      niter = 400; if (present(max_iter)) niter = max_iter
      tolerance = 1.0e-7_dp; if (present(tol)) tolerance = tol
      call bfgs_minimize(objective, par, fval, stat, max_iter=niter, tol=tolerance)

      allocate(beta(p, m)); call unpack_parameters(par, p, m, parallel, beta)
      allocate(alpha(n, m), means(n, m))
      do i = 1, n
         do j = 1, m
            alpha(i, j) = exp(clamp_eta(dot_product(x(i, :), beta(:, j))))
         end do
         means(i, :) = alpha(i, :)/sum(alpha(i, :))
      end do
      result%coefficients = beta
      result%fitted_shape = alpha
      result%fitted_mean = means
      result%parallel_slopes = parallel
      result%status = stat
      result%converged = stat == 0
      result%loglik = -fval
      result%aic = 2.0_dp*fval + 2.0_dp*real(np, dp)

      allocate(h(np, np)); call numerical_hessian(objective, par, h)
      call invert_matrix(h, covred, stat2)
      if (stat2 == 0) then
         call coefficient_map(p, m, parallel, map)
         covfull = matmul(map, matmul(covred, transpose(map)))
         result%covariance = covfull
      else
         allocate(result%covariance(0, 0))
      end if
   contains
      real(dp) function objective(theta) result(nll)
         real(dp), intent(in) :: theta(:)
         real(dp) :: arow(m), brow(p, m), ld
         integer :: row, col
         call unpack_parameters(theta, p, m, parallel, brow)
         nll = 0.0_dp
         do row = 1, n
            do col = 1, m
               arow(col) = exp(clamp_eta(dot_product(x(row, :), brow(:, col))))
            end do
            ld = dirichlet_logpdf(y(row, :), arow)
            if (.not. finite_scalar(ld)) then
               nll = huge(1.0_dp)/100.0_dp
               return
            end if
            nll = nll - w(row)*ld
         end do
      end function objective
   end subroutine fit_dirichlet_regression

   subroutine predict_dirichlet_regression(self, x, shape, mean)
      class(dirichlet_regression_result_t), intent(in) :: self
      real(dp), intent(in) :: x(:, :)
      real(dp), allocatable, intent(out) :: shape(:, :), mean(:, :)
      integer :: i, j, n, m
      n = size(x, 1); m = size(self%coefficients, 2)
      if (size(x, 2) /= size(self%coefficients, 1)) then
         allocate(shape(0, 0), mean(0, 0))
         return
      end if
      allocate(shape(n, m), mean(n, m))
      do i = 1, n
         do j = 1, m
            shape(i, j) = exp(clamp_eta(dot_product(x(i, :), self%coefficients(:, j))))
         end do
         mean(i, :) = shape(i, :)/sum(shape(i, :))
      end do
   end subroutine predict_dirichlet_regression

   subroutine initialize_parameters(y, x, w, parallel, par)
      real(dp), intent(in) :: y(:, :), x(:, :), w(:)
      logical, intent(in) :: parallel
      real(dp), intent(out) :: par(:)
      real(dp), allocatable :: mu(:), vv(:)
      real(dp) :: sw, conc, candidate
      integer :: j, m, p, used
      m = size(y, 2); p = size(x, 2); sw = max(sum(w), tiny(1.0_dp))
      allocate(mu(m), vv(m))
      do j = 1, m
         mu(j) = sum(w*y(:, j))/sw
         vv(j) = sum(w*(y(:, j) - mu(j))**2)/sw
      end do
      conc = 0.0_dp; used = 0
      do j = 1, m
         if (vv(j) > 1.0e-10_dp .and. mu(j) > 0.0_dp .and. mu(j) < 1.0_dp) then
            candidate = mu(j)*(1.0_dp - mu(j))/vv(j) - 1.0_dp
            if (candidate > 0.1_dp .and. candidate < 1.0e6_dp) then
               conc = conc + candidate
               used = used + 1
            end if
         end if
      end do
      if (used > 0) then
         conc = conc/real(used, dp)
      else
         conc = 10.0_dp
      end if
      conc = min(1000.0_dp, max(0.5_dp, conc))
      par = 0.0_dp
      if (is_intercept_design(x)) then
         if (parallel .and. p > 1) then
            do j = 1, m
               par(j) = log(max(mu(j)*conc, 1.0e-3_dp))
            end do
         else
            do j = 1, m
               par((j - 1)*p + 1) = log(max(mu(j)*conc, 1.0e-3_dp))
            end do
         end if
      end if
   end subroutine initialize_parameters

   subroutine unpack_parameters(par, p, m, parallel, beta)
      real(dp), intent(in) :: par(:)
      integer, intent(in) :: p, m
      logical, intent(in) :: parallel
      real(dp), intent(out) :: beta(p, m)
      integer :: j
      beta = 0.0_dp
      if (parallel .and. p > 1) then
         do j = 1, m
            beta(1, j) = par(j)
            beta(2:p, j) = par(m + 1:m + p - 1)
         end do
      else
         do j = 1, m
            beta(:, j) = par((j - 1)*p + 1:j*p)
         end do
      end if
   end subroutine unpack_parameters

   subroutine coefficient_map(p, m, parallel, map)
      integer, intent(in) :: p, m
      logical, intent(in) :: parallel
      real(dp), allocatable, intent(out) :: map(:, :)
      integer :: j, k, row, np
      if (parallel .and. p > 1) then
         np = m + p - 1
         allocate(map(p*m, np)); map = 0.0_dp
         do j = 1, m
            row = (j - 1)*p + 1
            map(row, j) = 1.0_dp
            do k = 2, p
               map(row + k - 1, m + k - 1) = 1.0_dp
            end do
         end do
      else
         np = p*m
         allocate(map(np, np)); map = 0.0_dp
         do j = 1, np
            map(j, j) = 1.0_dp
         end do
      end if
   end subroutine coefficient_map

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
end module vgam_dirichlet
