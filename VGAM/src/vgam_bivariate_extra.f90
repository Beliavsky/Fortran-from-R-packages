! Copyright (C) 1998-2025 T. W. Yee, University of Auckland.
! Modern Fortran computational translation, 2026.
! SPDX-License-Identifier: GPL-3.0-only
module vgam_bivariate_extra
   use vgam_kinds, only : dp, pi
   use vgam_distributions, only : rnorm_v, rexp_v
   use vgam_optim, only : bfgs_minimize, numerical_hessian
   use vgam_linalg, only : invert_matrix
   implicit none
   private

   type, public :: bivariate_normal_result_t
      real(dp), allocatable :: mean1_coefficients(:), mean2_coefficients(:)
      real(dp), allocatable :: sd1_coefficients(:), sd2_coefficients(:), rho_coefficients(:)
      real(dp), allocatable :: covariance(:, :)
      real(dp), allocatable :: fitted_mean1(:), fitted_mean2(:), fitted_sd1(:), fitted_sd2(:), fitted_rho(:)
      real(dp) :: loglik = -huge(1.0_dp), aic = huge(1.0_dp)
      integer :: status = 0
      logical :: converged = .false.
   end type bivariate_normal_result_t

   type, public :: bivariate_logistic_result_t
      real(dp), allocatable :: location1_coefficients(:), scale1_coefficients(:)
      real(dp), allocatable :: location2_coefficients(:), scale2_coefficients(:)
      real(dp), allocatable :: covariance(:, :)
      real(dp), allocatable :: fitted_location1(:), fitted_scale1(:), fitted_location2(:), fitted_scale2(:)
      real(dp) :: loglik = -huge(1.0_dp), aic = huge(1.0_dp)
      integer :: status = 0
      logical :: converged = .false.
   end type bivariate_logistic_result_t

   type, public :: freund61_result_t
      real(dp), allocatable :: coefficients(:, :)
      real(dp), allocatable :: covariance(:, :)
      real(dp), allocatable :: fitted_parameters(:, :)
      real(dp) :: loglik = -huge(1.0_dp), aic = huge(1.0_dp)
      integer :: status = 0
      logical :: converged = .false.
   end type freund61_result_t

   public :: bivariate_normal_pdf, bivariate_normal_logpdf, random_bivariate_normal
   public :: fit_bivariate_normal
   public :: bivariate_logistic_pdf, bivariate_logistic_logpdf, bivariate_logistic_cdf
   public :: random_bivariate_logistic, fit_bivariate_logistic
   public :: freund61_pdf, freund61_logpdf, random_freund61, fit_freund61

contains

   elemental real(dp) function bivariate_normal_logpdf(x1, x2, mean1, mean2, sd1, sd2, rho) result(ld)
      real(dp), intent(in) :: x1, x2, mean1, mean2, sd1, sd2, rho
      real(dp) :: z1, z2, den
      if (sd1 <= 0.0_dp .or. sd2 <= 0.0_dp .or. abs(rho) >= 1.0_dp) then
         ld = -huge(1.0_dp); return
      end if
      z1 = (x1 - mean1)/sd1; z2 = (x2 - mean2)/sd2; den = 1.0_dp - rho*rho
      ld = -log(2.0_dp*pi) - log(sd1) - log(sd2) - 0.5_dp*log(den) &
           - 0.5_dp*(z1*z1 - 2.0_dp*rho*z1*z2 + z2*z2)/den
   end function bivariate_normal_logpdf

   elemental real(dp) function bivariate_normal_pdf(x1, x2, mean1, mean2, sd1, sd2, rho) result(d)
      real(dp), intent(in) :: x1, x2, mean1, mean2, sd1, sd2, rho
      real(dp) :: ld
      ld = bivariate_normal_logpdf(x1, x2, mean1, mean2, sd1, sd2, rho)
      if (ld < -700.0_dp) then
         d = 0.0_dp
      else
         d = exp(ld)
      end if
   end function bivariate_normal_pdf

   subroutine random_bivariate_normal(mean1, mean2, sd1, sd2, rho, x1, x2)
      real(dp), intent(in) :: mean1, mean2, sd1, sd2, rho
      real(dp), intent(out) :: x1, x2
      real(dp) :: z1, z2
      if (sd1 <= 0.0_dp .or. sd2 <= 0.0_dp .or. abs(rho) >= 1.0_dp) then
         x1 = huge(1.0_dp); x2 = huge(1.0_dp); return
      end if
      z1 = rnorm_v(0.0_dp, 1.0_dp); z2 = rnorm_v(0.0_dp, 1.0_dp)
      x1 = mean1 + sd1*z1
      x2 = mean2 + sd2*(rho*z1 + sqrt(1.0_dp - rho*rho)*z2)
   end subroutine random_bivariate_normal

   subroutine fit_bivariate_normal(y1, y2, xm1, xm2, xs1, xs2, xr, result, weights, max_iter, tol)
      real(dp), intent(in) :: y1(:), y2(:), xm1(:, :), xm2(:, :), xs1(:, :), xs2(:, :), xr(:, :)
      type(bivariate_normal_result_t), intent(out) :: result
      real(dp), intent(in), optional :: weights(:), tol
      integer, intent(in), optional :: max_iter
      real(dp), allocatable :: w(:), par(:), h(:, :), cov(:, :)
      real(dp) :: fval, tolerance, m1, m2, s1, s2
      integer :: n, p1, p2, p3, p4, p5, np, stat, stat2, niter, k
      n = size(y1); p1 = size(xm1, 2); p2 = size(xm2, 2); p3 = size(xs1, 2)
      p4 = size(xs2, 2); p5 = size(xr, 2); np = p1 + p2 + p3 + p4 + p5
      if (n <= 1 .or. size(y2) /= n .or. min(p1, p2, p3, p4, p5) <= 0 .or. &
          size(xm1, 1) /= n .or. size(xm2, 1) /= n .or. size(xs1, 1) /= n .or. &
          size(xs2, 1) /= n .or. size(xr, 1) /= n) then
         result%status = 1; return
      end if
      allocate(w(n)); w = 1.0_dp
      if (present(weights)) then
         if (size(weights) /= n .or. any(weights < 0.0_dp)) then
            result%status = 2; return
         end if
         w = weights
      end if
      allocate(par(np)); par = 0.0_dp
      m1 = sum(w*y1)/max(sum(w), tiny(1.0_dp)); m2 = sum(w*y2)/max(sum(w), tiny(1.0_dp))
      s1 = sqrt(max(sum(w*(y1 - m1)**2)/max(sum(w), 1.0_dp), 1.0e-6_dp))
      s2 = sqrt(max(sum(w*(y2 - m2)**2)/max(sum(w), 1.0_dp), 1.0e-6_dp))
      k = 0
      if (is_intercept(xm1)) par(k + 1) = m1; k = k + p1
      if (is_intercept(xm2)) par(k + 1) = m2; k = k + p2
      if (is_intercept(xs1)) par(k + 1) = log(s1); k = k + p3
      if (is_intercept(xs2)) par(k + 1) = log(s2); k = k + p4
      if (is_intercept(xr)) par(k + 1) = atanh(max(-0.95_dp, min(0.95_dp, weighted_corr(y1, y2, w))))
      niter = 350; if (present(max_iter)) niter = max_iter
      tolerance = 1.0e-7_dp; if (present(tol)) tolerance = tol
      call bfgs_minimize(objective, par, fval, stat, max_iter=niter, tol=tolerance)
      result%status = stat; result%converged = stat == 0; result%loglik = -fval
      result%aic = 2.0_dp*fval + 2.0_dp*real(np, dp)
      k = 0; result%mean1_coefficients = par(k + 1:k + p1); k = k + p1
      result%mean2_coefficients = par(k + 1:k + p2); k = k + p2
      result%sd1_coefficients = par(k + 1:k + p3); k = k + p3
      result%sd2_coefficients = par(k + 1:k + p4); k = k + p4
      result%rho_coefficients = par(k + 1:k + p5)
      result%fitted_mean1 = matmul(xm1, result%mean1_coefficients)
      result%fitted_mean2 = matmul(xm2, result%mean2_coefficients)
      result%fitted_sd1 = exp(clamp_eta(matmul(xs1, result%sd1_coefficients)))
      result%fitted_sd2 = exp(clamp_eta(matmul(xs2, result%sd2_coefficients)))
      result%fitted_rho = tanh(matmul(xr, result%rho_coefficients))
      allocate(h(np, np)); call numerical_hessian(objective, par, h); call invert_matrix(h, cov, stat2)
      if (stat2 == 0) then; result%covariance = cov; else; allocate(result%covariance(0, 0)); end if
   contains
      real(dp) function objective(theta) result(nll)
         real(dp), intent(in) :: theta(:)
         real(dp) :: a, b, c, d, r, ld
         integer :: row, pos
         nll = 0.0_dp
         do row = 1, n
            pos = 0; a = dot_product(xm1(row, :), theta(pos + 1:pos + p1)); pos = pos + p1
            b = dot_product(xm2(row, :), theta(pos + 1:pos + p2)); pos = pos + p2
            c = exp(clamp_eta(dot_product(xs1(row, :), theta(pos + 1:pos + p3)))); pos = pos + p3
            d = exp(clamp_eta(dot_product(xs2(row, :), theta(pos + 1:pos + p4)))); pos = pos + p4
            r = tanh(dot_product(xr(row, :), theta(pos + 1:pos + p5)))
            ld = bivariate_normal_logpdf(y1(row), y2(row), a, b, c, d, r)
            if (.not. finite_scalar(ld)) then; nll = huge(1.0_dp)/100.0_dp; return; end if
            nll = nll - w(row)*ld
         end do
      end function objective
   end subroutine fit_bivariate_normal

   elemental real(dp) function bivariate_logistic_logpdf(x1, x2, loc1, scale1, loc2, scale2) result(ld)
      real(dp), intent(in) :: x1, x2, loc1, scale1, loc2, scale2
      real(dp) :: z1, z2, vmax, lden
      if (scale1 <= 0.0_dp .or. scale2 <= 0.0_dp) then
         ld = -huge(1.0_dp); return
      end if
      z1 = (x1 - loc1)/scale1; z2 = (x2 - loc2)/scale2
      vmax = max(0.0_dp, max(-z1, -z2))
      lden = vmax + log(exp(-vmax) + exp(-z1 - vmax) + exp(-z2 - vmax))
      ld = log(2.0_dp) - z1 - z2 - log(scale1) - log(scale2) - 3.0_dp*lden
   end function bivariate_logistic_logpdf

   elemental real(dp) function bivariate_logistic_pdf(x1, x2, loc1, scale1, loc2, scale2) result(d)
      real(dp), intent(in) :: x1, x2, loc1, scale1, loc2, scale2
      real(dp) :: ld
      ld = bivariate_logistic_logpdf(x1, x2, loc1, scale1, loc2, scale2)
      d = merge(exp(ld), 0.0_dp, ld > -700.0_dp)
   end function bivariate_logistic_pdf

   elemental real(dp) function bivariate_logistic_cdf(x1, x2, loc1, scale1, loc2, scale2) result(p)
      real(dp), intent(in) :: x1, x2, loc1, scale1, loc2, scale2
      if (scale1 <= 0.0_dp .or. scale2 <= 0.0_dp) then
         p = -1.0_dp
      else
         p = 1.0_dp/(1.0_dp + exp(min(700.0_dp, -(x1 - loc1)/scale1)) + &
            exp(min(700.0_dp, -(x2 - loc2)/scale2)))
      end if
   end function bivariate_logistic_cdf

   subroutine random_bivariate_logistic(loc1, scale1, loc2, scale2, x1, x2)
      real(dp), intent(in) :: loc1, scale1, loc2, scale2
      real(dp), intent(out) :: x1, x2
      real(dp) :: u1, u2, ez1
      if (scale1 <= 0.0_dp .or. scale2 <= 0.0_dp) then; x1 = huge(1.0_dp); x2 = huge(1.0_dp); return; end if
      call random_number(u1); call random_number(u2)
      u1 = min(1.0_dp - epsilon(1.0_dp), max(epsilon(1.0_dp), u1))
      u2 = min(1.0_dp - epsilon(1.0_dp), max(epsilon(1.0_dp), u2))
      x1 = loc1 + scale1*log(u1/(1.0_dp - u1)); ez1 = exp(-(x1 - loc1)/scale1)
      x2 = loc2 - scale2*log(1.0_dp/sqrt(u2/(1.0_dp + ez1)**2) - 1.0_dp - ez1)
   end subroutine random_bivariate_logistic

   subroutine fit_bivariate_logistic(y1, y2, x, result, weights, max_iter, tol)
      real(dp), intent(in) :: y1(:), y2(:), x(:, :)
      type(bivariate_logistic_result_t), intent(out) :: result
      real(dp), intent(in), optional :: weights(:), tol
      integer, intent(in), optional :: max_iter
      real(dp), allocatable :: w(:), par(:), h(:, :), cov(:, :)
      real(dp) :: fval, tolerance, m1, m2, s1, s2
      integer :: n, p, np, stat, stat2, niter
      n = size(y1); p = size(x, 2); np = 4*p
      if (n <= 1 .or. size(y2) /= n .or. size(x, 1) /= n .or. p <= 0) then; result%status = 1; return; end if
      allocate(w(n)); w = 1.0_dp
      if (present(weights)) then
         if (size(weights) /= n .or. any(weights < 0.0_dp)) then; result%status = 2; return; end if
         w = weights
      end if
      allocate(par(np)); par = 0.0_dp
      m1 = sum(w*y1)/sum(w); m2 = sum(w*y2)/sum(w)
      s1 = sqrt(max(3.0_dp*sum(w*(y1 - m1)**2)/sum(w), 1.0e-6_dp))/pi
      s2 = sqrt(max(3.0_dp*sum(w*(y2 - m2)**2)/sum(w), 1.0e-6_dp))/pi
      if (is_intercept(x)) then; par(1) = m1; par(p + 1) = log(s1); par(2*p + 1) = m2; par(3*p + 1) = log(s2); end if
      niter = 350; if (present(max_iter)) niter = max_iter
      tolerance = 1.0e-7_dp; if (present(tol)) tolerance = tol
      call bfgs_minimize(objective, par, fval, stat, max_iter=niter, tol=tolerance)
      result%status = stat; result%converged = stat == 0; result%loglik = -fval
      result%aic = 2.0_dp*fval + 2.0_dp*real(np, dp)
      result%location1_coefficients = par(1:p); result%scale1_coefficients = par(p + 1:2*p)
      result%location2_coefficients = par(2*p + 1:3*p); result%scale2_coefficients = par(3*p + 1:4*p)
      result%fitted_location1 = matmul(x, result%location1_coefficients)
      result%fitted_scale1 = exp(clamp_eta(matmul(x, result%scale1_coefficients)))
      result%fitted_location2 = matmul(x, result%location2_coefficients)
      result%fitted_scale2 = exp(clamp_eta(matmul(x, result%scale2_coefficients)))
      allocate(h(np, np)); call numerical_hessian(objective, par, h); call invert_matrix(h, cov, stat2)
      if (stat2 == 0) then; result%covariance = cov; else; allocate(result%covariance(0, 0)); end if
   contains
      real(dp) function objective(theta) result(nll)
         real(dp), intent(in) :: theta(:)
         real(dp) :: l1, l2, sc1, sc2, ld
         integer :: row
         nll = 0.0_dp
         do row = 1, n
            l1 = dot_product(x(row, :), theta(1:p)); sc1 = exp(clamp_eta(dot_product(x(row, :), theta(p + 1:2*p))))
            l2 = dot_product(x(row, :), theta(2*p + 1:3*p))
            sc2 = exp(clamp_eta(dot_product(x(row, :), theta(3*p + 1:4*p))))
            ld = bivariate_logistic_logpdf(y1(row), y2(row), l1, sc1, l2, sc2)
            if (.not. finite_scalar(ld)) then; nll = huge(1.0_dp)/100.0_dp; return; end if
            nll = nll - w(row)*ld
         end do
      end function objective
   end subroutine fit_bivariate_logistic

   elemental real(dp) function freund61_logpdf(x1, x2, a, ap, b, bp) result(ld)
      real(dp), intent(in) :: x1, x2, a, ap, b, bp
      if (min(a, ap, b, bp) <= 0.0_dp .or. min(x1, x2) < 0.0_dp) then
         ld = -huge(1.0_dp)
      else if (x1 < x2) then
         ld = log(a) + log(bp) - bp*x2 - (a + b - bp)*x1
      else
         ld = log(b) + log(ap) - ap*x1 - (a + b - ap)*x2
      end if
   end function freund61_logpdf

   elemental real(dp) function freund61_pdf(x1, x2, a, ap, b, bp) result(d)
      real(dp), intent(in) :: x1, x2, a, ap, b, bp
      real(dp) :: ld
      ld = freund61_logpdf(x1, x2, a, ap, b, bp); d = merge(exp(ld), 0.0_dp, ld > -700.0_dp)
   end function freund61_pdf

   subroutine random_freund61(a, ap, b, bp, x1, x2)
      real(dp), intent(in) :: a, ap, b, bp
      real(dp), intent(out) :: x1, x2
      real(dp) :: t1, t2
      if (min(a, ap, b, bp) <= 0.0_dp) then; x1 = huge(1.0_dp); x2 = huge(1.0_dp); return; end if
      t1 = rexp_v(a); t2 = rexp_v(b)
      if (t1 < t2) then
         x1 = t1; x2 = t1 + rexp_v(bp)
      else
         x2 = t2; x1 = t2 + rexp_v(ap)
      end if
   end subroutine random_freund61

   subroutine fit_freund61(y1, y2, x, result, weights, independent, max_iter, tol)
      real(dp), intent(in) :: y1(:), y2(:), x(:, :)
      type(freund61_result_t), intent(out) :: result
      real(dp), intent(in), optional :: weights(:), tol
      logical, intent(in), optional :: independent
      integer, intent(in), optional :: max_iter
      real(dp), allocatable :: w(:), par(:), h(:, :), cov(:, :)
      real(dp) :: fval, tolerance, mean1, mean2
      logical :: indep
      integer :: n, p, q, np, stat, stat2, niter, row
      n = size(y1); p = size(x, 2); indep = .false.; if (present(independent)) indep = independent
      q = merge(2, 4, indep); np = p*q
      if (n <= 2 .or. size(y2) /= n .or. size(x, 1) /= n .or. p <= 0 .or. any(y1 < 0.0_dp) .or. any(y2 < 0.0_dp)) then
         result%status = 1; return
      end if
      if (all(y1 < y2) .or. all(y2 <= y1)) then; result%status = 2; return; end if
      allocate(w(n)); w = 1.0_dp
      if (present(weights)) then
         if (size(weights) /= n .or. any(weights < 0.0_dp)) then; result%status = 3; return; end if
         w = weights
      end if
      allocate(par(np)); par = 0.0_dp
      mean1 = sum(w*y1)/sum(w); mean2 = sum(w*y2)/sum(w)
      if (is_intercept(x)) then
         par(1) = log(1.0_dp/max(mean1, 0.1_dp)); par(p + 1) = log(1.0_dp/max(mean2, 0.1_dp))
         if (.not. indep) then; par(2*p + 1) = par(1); par(3*p + 1) = par(p + 1); end if
      end if
      niter = 400; if (present(max_iter)) niter = max_iter
      tolerance = 1.0e-7_dp; if (present(tol)) tolerance = tol
      call bfgs_minimize(objective, par, fval, stat, max_iter=niter, tol=tolerance)
      result%status = stat; result%converged = stat == 0; result%loglik = -fval
      result%aic = 2.0_dp*fval + 2.0_dp*real(np, dp)
      allocate(result%coefficients(p, 4)); result%coefficients = 0.0_dp
      if (indep) then
         result%coefficients(:, 1) = par(1:p); result%coefficients(:, 2) = par(1:p)
         result%coefficients(:, 3) = par(p + 1:2*p); result%coefficients(:, 4) = par(p + 1:2*p)
      else
         result%coefficients = reshape(par, [p, 4])
      end if
      allocate(result%fitted_parameters(n, 4))
      do row = 1, n
         result%fitted_parameters(row, :) = exp(clamp_eta(matmul(transpose(result%coefficients), x(row, :))))
      end do
      allocate(h(np, np)); call numerical_hessian(objective, par, h); call invert_matrix(h, cov, stat2)
      if (stat2 == 0) then; result%covariance = cov; else; allocate(result%covariance(0, 0)); end if
   contains
      real(dp) function objective(theta) result(nll)
         real(dp), intent(in) :: theta(:)
         real(dp) :: av, apv, bv, bpv, ld
         integer :: r
         nll = 0.0_dp
         do r = 1, n
            if (indep) then
               av = exp(clamp_eta(dot_product(x(r, :), theta(1:p)))); apv = av
               bv = exp(clamp_eta(dot_product(x(r, :), theta(p + 1:2*p)))); bpv = bv
            else
               av = exp(clamp_eta(dot_product(x(r, :), theta(1:p))))
               apv = exp(clamp_eta(dot_product(x(r, :), theta(p + 1:2*p))))
               bv = exp(clamp_eta(dot_product(x(r, :), theta(2*p + 1:3*p))))
               bpv = exp(clamp_eta(dot_product(x(r, :), theta(3*p + 1:4*p))))
            end if
            ld = freund61_logpdf(y1(r), y2(r), av, apv, bv, bpv)
            if (.not. finite_scalar(ld)) then; nll = huge(1.0_dp)/100.0_dp; return; end if
            nll = nll - w(r)*ld
         end do
      end function objective
   end subroutine fit_freund61

   logical function is_intercept(x) result(ans)
      real(dp), intent(in) :: x(:, :)
      ans = size(x, 2) >= 1 .and. all(abs(x(:, 1) - 1.0_dp) < 100.0_dp*epsilon(1.0_dp))
   end function is_intercept

   real(dp) function weighted_corr(x, y, w) result(r)
      real(dp), intent(in) :: x(:), y(:), w(:)
      real(dp) :: sw, mx, my, vx, vy, cv
      sw = sum(w); mx = sum(w*x)/sw; my = sum(w*y)/sw
      vx = sum(w*(x - mx)**2); vy = sum(w*(y - my)**2); cv = sum(w*(x - mx)*(y - my))
      r = cv/sqrt(max(vx*vy, tiny(1.0_dp)))
   end function weighted_corr

   elemental real(dp) function clamp_eta(x) result(y)
      real(dp), intent(in) :: x
      y = min(30.0_dp, max(-30.0_dp, x))
   end function clamp_eta

   elemental logical function finite_scalar(x) result(ok)
      real(dp), intent(in) :: x
      ok = abs(x) <= huge(x)
   end function finite_scalar

end module vgam_bivariate_extra
