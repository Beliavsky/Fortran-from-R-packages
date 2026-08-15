! Copyright (C) 1998-2025 T. W. Yee, University of Auckland.
! Modern Fortran computational translation, 2026.
! SPDX-License-Identifier: GPL-3.0-only
module vgam_student_t
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_value, ieee_quiet_nan
   use vgam_kinds, only : dp, pi
   use vgam_special, only : regularized_beta
   use vgam_random, only : random_normal, random_gamma
   use vgam_optim, only : bfgs_minimize, numerical_hessian
   use vgam_linalg, only : invert_matrix
   implicit none
   private

   type, public :: bivariate_student_t_result_t
      real(dp), allocatable :: df_coefficients(:)
      real(dp), allocatable :: rho_coefficients(:)
      real(dp), allocatable :: covariance(:, :)
      real(dp), allocatable :: fitted_df(:)
      real(dp), allocatable :: fitted_rho(:)
      real(dp) :: loglik = -huge(1.0_dp)
      real(dp) :: aic = huge(1.0_dp)
      integer :: status = 0
      logical :: converged = .false.
   end type bivariate_student_t_result_t

   type, public :: student_t_copula_result_t
      real(dp), allocatable :: rho_coefficients(:)
      real(dp), allocatable :: covariance(:, :)
      real(dp), allocatable :: fitted_rho(:)
      real(dp) :: df = 0.0_dp
      real(dp) :: loglik = -huge(1.0_dp)
      real(dp) :: aic = huge(1.0_dp)
      integer :: status = 0
      logical :: converged = .false.
   end type student_t_copula_result_t

   public :: student_t_pdf, student_t_logpdf, student_t_cdf, student_t_quantile, random_student_t
   public :: bivariate_student_t_pdf, bivariate_student_t_logpdf, random_bivariate_student_t
   public :: student_t_copula_pdf, random_student_t_copula
   public :: fit_bivariate_student_t, fit_student_t_copula

contains

   elemental real(dp) function student_t_logpdf(x, df) result(ld)
      real(dp), intent(in) :: x, df
      if (df <= 0.0_dp) then
         ld = ieee_value(0.0_dp, ieee_quiet_nan)
      else if (.not. ieee_is_finite(x)) then
         ld = -huge(1.0_dp)
      else
         ld = log_gamma(0.5_dp*(df + 1.0_dp)) - log_gamma(0.5_dp*df) &
              - 0.5_dp*log(df*pi) - 0.5_dp*(df + 1.0_dp)*log(1.0_dp + x*x/df)
      end if
   end function student_t_logpdf

   elemental real(dp) function student_t_pdf(x, df) result(d)
      real(dp), intent(in) :: x, df
      real(dp) :: ld
      ld = student_t_logpdf(x, df)
      if (.not. ieee_is_finite(ld)) then
         if (ld < 0.0_dp) then
            d = 0.0_dp
         else
            d = ieee_value(0.0_dp, ieee_quiet_nan)
         end if
      else
         d = exp(ld)
      end if
   end function student_t_pdf

   elemental real(dp) function student_t_cdf(x, df) result(p)
      real(dp), intent(in) :: x, df
      real(dp) :: ib, z
      if (df <= 0.0_dp) then
         p = ieee_value(0.0_dp, ieee_quiet_nan)
      else if (x == 0.0_dp) then
         p = 0.5_dp
      else if (.not. ieee_is_finite(x)) then
         p = merge(1.0_dp, 0.0_dp, x > 0.0_dp)
      else
         z = df/(df + x*x)
         ib = regularized_beta(z, 0.5_dp*df, 0.5_dp)
         if (x > 0.0_dp) then
            p = 1.0_dp - 0.5_dp*ib
         else
            p = 0.5_dp*ib
         end if
         p = min(1.0_dp, max(0.0_dp, p))
      end if
   end function student_t_cdf

   real(dp) function student_t_quantile(probability, df) result(q)
      real(dp), intent(in) :: probability, df
      real(dp) :: lo, hi, mid
      integer :: i
      if (df <= 0.0_dp .or. probability < 0.0_dp .or. probability > 1.0_dp) then
         q = ieee_value(0.0_dp, ieee_quiet_nan)
         return
      end if
      if (probability == 0.0_dp) then
         q = -huge(1.0_dp)
         return
      else if (probability == 1.0_dp) then
         q = huge(1.0_dp)
         return
      else if (probability == 0.5_dp) then
         q = 0.0_dp
         return
      end if
      lo = -1.0_dp
      hi = 1.0_dp
      do while (student_t_cdf(lo, df) > probability)
         hi = lo
         lo = 2.0_dp*lo
         if (lo < -1.0e12_dp) exit
      end do
      do while (student_t_cdf(hi, df) < probability)
         lo = hi
         hi = 2.0_dp*hi
         if (hi > 1.0e12_dp) exit
      end do
      do i = 1, 120
         mid = 0.5_dp*(lo + hi)
         if (student_t_cdf(mid, df) < probability) then
            lo = mid
         else
            hi = mid
         end if
      end do
      q = 0.5_dp*(lo + hi)
   end function student_t_quantile

   real(dp) function random_student_t(df) result(x)
      real(dp), intent(in) :: df
      real(dp) :: chi2
      if (df <= 0.0_dp) then
         x = ieee_value(0.0_dp, ieee_quiet_nan)
         return
      end if
      chi2 = random_gamma(0.5_dp*df, 2.0_dp)
      x = random_normal()/sqrt(chi2/df)
   end function random_student_t

   elemental real(dp) function bivariate_student_t_logpdf(x1, x2, df, rho) result(ld)
      real(dp), intent(in) :: x1, x2, df, rho
      real(dp) :: den, quad
      if (df <= 0.0_dp .or. abs(rho) >= 1.0_dp) then
         ld = ieee_value(0.0_dp, ieee_quiet_nan)
      else if (.not. ieee_is_finite(x1) .or. .not. ieee_is_finite(x2)) then
         ld = -huge(1.0_dp)
      else
         den = 1.0_dp - rho*rho
         quad = (x1*x1 + x2*x2 - 2.0_dp*rho*x1*x2)/(df*den)
         ld = -(0.5_dp*df + 1.0_dp)*log(1.0_dp + quad) - log(2.0_dp*pi) - 0.5_dp*log(den)
      end if
   end function bivariate_student_t_logpdf

   elemental real(dp) function bivariate_student_t_pdf(x1, x2, df, rho) result(d)
      real(dp), intent(in) :: x1, x2, df, rho
      real(dp) :: ld
      ld = bivariate_student_t_logpdf(x1, x2, df, rho)
      if (.not. ieee_is_finite(ld)) then
         if (ld < 0.0_dp) then
            d = 0.0_dp
         else
            d = ieee_value(0.0_dp, ieee_quiet_nan)
         end if
      else
         d = exp(ld)
      end if
   end function bivariate_student_t_pdf

   subroutine random_bivariate_student_t(df, rho, x1, x2)
      real(dp), intent(in) :: df, rho
      real(dp), intent(out) :: x1, x2
      real(dp) :: z1, z2, scale, chi2
      if (df <= 0.0_dp .or. abs(rho) >= 1.0_dp) then
         x1 = ieee_value(0.0_dp, ieee_quiet_nan)
         x2 = x1
         return
      end if
      z1 = random_normal()
      z2 = rho*z1 + sqrt(1.0_dp - rho*rho)*random_normal()
      chi2 = random_gamma(0.5_dp*df, 2.0_dp)
      scale = sqrt(df/chi2)
      x1 = z1*scale
      x2 = z2*scale
   end subroutine random_bivariate_student_t

   real(dp) function student_t_copula_pdf(u, v, df, rho) result(d)
      real(dp), intent(in) :: u, v, df, rho
      real(dp) :: x1, x2, ld
      if (u <= 0.0_dp .or. u >= 1.0_dp .or. v <= 0.0_dp .or. v >= 1.0_dp .or. &
          df <= 0.0_dp .or. abs(rho) >= 1.0_dp) then
         if (u <= 0.0_dp .or. u >= 1.0_dp .or. v <= 0.0_dp .or. v >= 1.0_dp) then
            d = 0.0_dp
         else
            d = ieee_value(0.0_dp, ieee_quiet_nan)
         end if
         return
      end if
      x1 = student_t_quantile(u, df)
      x2 = student_t_quantile(v, df)
      ld = bivariate_student_t_logpdf(x1, x2, df, rho) &
           - student_t_logpdf(x1, df) - student_t_logpdf(x2, df)
      d = exp(min(700.0_dp, ld))
   end function student_t_copula_pdf

   subroutine random_student_t_copula(df, rho, u, v)
      real(dp), intent(in) :: df, rho
      real(dp), intent(out) :: u, v
      real(dp) :: x1, x2
      call random_bivariate_student_t(df, rho, x1, x2)
      u = student_t_cdf(x1, df)
      v = student_t_cdf(x2, df)
   end subroutine random_student_t_copula

   subroutine fit_bivariate_student_t(x1, x2, xdf, xrho, result, weights, max_iter, tol)
      real(dp), intent(in) :: x1(:), x2(:), xdf(:, :), xrho(:, :)
      type(bivariate_student_t_result_t), intent(out) :: result
      real(dp), intent(in), optional :: weights(:), tol
      integer, intent(in), optional :: max_iter
      real(dp), allocatable :: w(:), par(:), hess(:, :), cov(:, :)
      real(dp) :: fval, tolerance, r0
      integer :: n, pd, pr, np, stat, stat2, niter, i

      n = size(x1); pd = size(xdf, 2); pr = size(xrho, 2); np = pd + pr
      if (n <= 2 .or. size(x2) /= n .or. size(xdf, 1) /= n .or. size(xrho, 1) /= n .or. &
          pd <= 0 .or. pr <= 0) then
         result%status = 1
         return
      end if
      allocate(w(n)); w = 1.0_dp
      if (present(weights)) then
         if (size(weights) /= n .or. any(weights < 0.0_dp)) then
            result%status = 2
            return
         end if
         w = weights
      end if
      allocate(par(np)); par = 0.0_dp
      if (all(abs(xdf(:, 1) - 1.0_dp) < 100.0_dp*epsilon(1.0_dp))) par(1) = log(9.0_dp)
      if (all(abs(xrho(:, 1) - 1.0_dp) < 100.0_dp*epsilon(1.0_dp))) then
         r0 = sample_corr(x1, x2)
         par(pd + 1) = atanh(min(0.95_dp, max(-0.95_dp, r0)))
      end if
      niter = 350; if (present(max_iter)) niter = max_iter
      tolerance = 1.0e-7_dp; if (present(tol)) tolerance = tol
      call bfgs_minimize(objective, par, fval, stat, max_iter=niter, tol=tolerance)
      result%status = stat; result%converged = stat == 0
      result%df_coefficients = par(1:pd); result%rho_coefficients = par(pd + 1:np)
      allocate(result%fitted_df(n), result%fitted_rho(n))
      do i = 1, n
         result%fitted_df(i) = 1.0_dp + exp(min(30.0_dp, max(-30.0_dp, dot_product(xdf(i, :), par(1:pd)))))
         result%fitted_rho(i) = tanh(dot_product(xrho(i, :), par(pd + 1:np)))
      end do
      result%loglik = -fval; result%aic = 2.0_dp*fval + 2.0_dp*real(np, dp)
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
         real(dp) :: dfv, rhov, ld
         integer :: j
         nll = 0.0_dp
         do j = 1, n
            dfv = 1.0_dp + exp(min(30.0_dp, max(-30.0_dp, dot_product(xdf(j, :), theta(1:pd)))))
            rhov = tanh(dot_product(xrho(j, :), theta(pd + 1:np)))
            ld = bivariate_student_t_logpdf(x1(j), x2(j), dfv, rhov)
            if (.not. ieee_is_finite(ld)) then
               nll = huge(1.0_dp)/100.0_dp
               return
            end if
            nll = nll - w(j)*ld
         end do
      end function objective
   end subroutine fit_bivariate_student_t

   subroutine fit_student_t_copula(u, v, xrho, result, estimate_df, initial_df, weights, max_iter, tol)
      real(dp), intent(in) :: u(:), v(:), xrho(:, :)
      type(student_t_copula_result_t), intent(out) :: result
      logical, intent(in), optional :: estimate_df
      real(dp), intent(in), optional :: initial_df, weights(:), tol
      integer, intent(in), optional :: max_iter
      real(dp), allocatable :: w(:), par(:), hess(:, :), cov(:, :), xu(:), xv(:)
      real(dp) :: fval, tolerance, df0, r0
      integer :: n, pr, np, stat, stat2, niter, i
      logical :: estdf

      n = size(u); pr = size(xrho, 2); estdf = .true.; if (present(estimate_df)) estdf = estimate_df
      np = pr + merge(1, 0, estdf)
      if (n <= 2 .or. size(v) /= n .or. size(xrho, 1) /= n .or. pr <= 0 .or. &
          any(u <= 0.0_dp) .or. any(u >= 1.0_dp) .or. any(v <= 0.0_dp) .or. any(v >= 1.0_dp)) then
         result%status = 1
         return
      end if
      allocate(w(n)); w = 1.0_dp
      if (present(weights)) then
         if (size(weights) /= n .or. any(weights < 0.0_dp)) then
            result%status = 2
            return
         end if
         w = weights
      end if
      df0 = 8.0_dp; if (present(initial_df)) df0 = max(0.2_dp, initial_df)
      allocate(xu(n), xv(n))
      if (.not. estdf) then
         do i = 1, n
            xu(i) = student_t_quantile(u(i), df0)
            xv(i) = student_t_quantile(v(i), df0)
         end do
      else
         xu = 0.0_dp
         xv = 0.0_dp
      end if
      allocate(par(np)); par = 0.0_dp
      r0 = normal_score_corr(u, v, df0)
      if (all(abs(xrho(:, 1) - 1.0_dp) < 100.0_dp*epsilon(1.0_dp))) &
         par(1) = atanh(min(0.95_dp, max(-0.95_dp, r0)))
      if (estdf) par(np) = log(df0)
      niter = 300; if (present(max_iter)) niter = max_iter
      tolerance = 1.0e-7_dp; if (present(tol)) tolerance = tol
      call bfgs_minimize(objective, par, fval, stat, max_iter=niter, tol=tolerance)
      result%status = stat; result%converged = stat == 0
      result%rho_coefficients = par(1:pr)
      result%df = merge(exp(min(20.0_dp, max(-2.0_dp, par(np)))), df0, estdf)
      allocate(result%fitted_rho(n))
      do i = 1, n
         result%fitted_rho(i) = tanh(dot_product(xrho(i, :), par(1:pr)))
      end do
      result%loglik = -fval; result%aic = 2.0_dp*fval + 2.0_dp*real(np, dp)
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
         real(dp) :: dfv, rhov, dens
         integer :: j
         dfv = df0
         if (estdf) dfv = exp(min(20.0_dp, max(-2.0_dp, theta(np))))
         nll = 0.0_dp
         do j = 1, n
            rhov = tanh(dot_product(xrho(j, :), theta(1:pr)))
            if (estdf) then
               dens = student_t_copula_pdf(u(j), v(j), dfv, rhov)
            else
               dens = exp(min(700.0_dp, bivariate_student_t_logpdf(xu(j), xv(j), dfv, rhov) &
                  - student_t_logpdf(xu(j), dfv) - student_t_logpdf(xv(j), dfv)))
            end if
            if (dens <= tiny(1.0_dp) .or. .not. ieee_is_finite(dens)) then
               nll = huge(1.0_dp)/100.0_dp
               return
            end if
            nll = nll - w(j)*log(dens)
         end do
      end function objective
   end subroutine fit_student_t_copula

   real(dp) function sample_corr(x, y) result(r)
      real(dp), intent(in) :: x(:), y(:)
      real(dp) :: mx, my, sx, sy
      mx = sum(x)/real(size(x), dp); my = sum(y)/real(size(y), dp)
      sx = sum((x - mx)**2); sy = sum((y - my)**2)
      if (sx <= tiny(1.0_dp) .or. sy <= tiny(1.0_dp)) then
         r = 0.0_dp
      else
         r = sum((x - mx)*(y - my))/sqrt(sx*sy)
      end if
   end function sample_corr

   real(dp) function normal_score_corr(u, v, df) result(r)
      real(dp), intent(in) :: u(:), v(:), df
      real(dp), allocatable :: x(:), y(:)
      integer :: i
      allocate(x(size(u)), y(size(v)))
      do i = 1, size(u)
         x(i) = student_t_quantile(u(i), df)
         y(i) = student_t_quantile(v(i), df)
      end do
      r = sample_corr(x, y)
   end function normal_score_corr

end module vgam_student_t
