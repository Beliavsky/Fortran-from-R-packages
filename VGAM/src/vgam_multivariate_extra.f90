! Copyright (C) 1998-2025 T. W. Yee, University of Auckland.
! Modern Fortran computational translation, 2026.
! SPDX-License-Identifier: GPL-3.0-only
module vgam_multivariate_extra
   use vgam_kinds, only : dp, log2pi
   use vgam_distributions, only : rnorm_v
   use vgam_copulas, only : random_fgm_copula
   use vgam_optim, only : bfgs_minimize, numerical_hessian
   use vgam_linalg, only : invert_matrix
   implicit none
   private

   type, public :: trivariate_normal_result_t
      real(dp) :: mean(3) = 0.0_dp
      real(dp) :: sd(3) = 1.0_dp
      real(dp) :: rho12 = 0.0_dp, rho13 = 0.0_dp, rho23 = 0.0_dp
      real(dp), allocatable :: covariance(:, :)
      real(dp) :: loglik = -huge(1.0_dp), aic = huge(1.0_dp)
      integer :: status = 0
      logical :: converged = .false.
   end type trivariate_normal_result_t

   type, public :: bifgm_exponential_result_t
      real(dp), allocatable :: coefficients(:)
      real(dp), allocatable :: covariance(:, :)
      real(dp), allocatable :: fitted_alpha(:)
      real(dp) :: loglik = -huge(1.0_dp), aic = huge(1.0_dp)
      integer :: status = 0
      logical :: converged = .false.
   end type bifgm_exponential_result_t

   public :: trivariate_normal_pdf, trivariate_normal_logpdf, random_trivariate_normal
   public :: fit_trivariate_normal
   public :: bifgm_exponential_pdf, bifgm_exponential_logpdf, bifgm_exponential_cdf
   public :: random_bifgm_exponential, fit_bifgm_exponential
   public :: kendall_tau_v

contains

   elemental real(dp) function trivariate_normal_logpdf(x1, x2, x3, mean1, mean2, mean3, &
      sd1, sd2, sd3, rho12, rho13, rho23) result(ld)
      real(dp), intent(in) :: x1, x2, x3, mean1, mean2, mean3, sd1, sd2, sd3
      real(dp), intent(in) :: rho12, rho13, rho23
      real(dp) :: z1, z2, z3, det, quad
      if (min(sd1, sd2, sd3) <= 0.0_dp .or. max(abs(rho12), abs(rho13), abs(rho23)) >= 1.0_dp) then
         ld = -huge(1.0_dp); return
      end if
      det = correlation_determinant(rho12, rho13, rho23)
      if (det <= 0.0_dp) then; ld = -huge(1.0_dp); return; end if
      z1 = (x1 - mean1)/sd1; z2 = (x2 - mean2)/sd2; z3 = (x3 - mean3)/sd3
      quad = (1.0_dp - rho23*rho23)*z1*z1 + (1.0_dp - rho13*rho13)*z2*z2 + &
         (1.0_dp - rho12*rho12)*z3*z3 + 2.0_dp*(rho13*rho23 - rho12)*z1*z2 + &
         2.0_dp*(rho12*rho23 - rho13)*z1*z3 + 2.0_dp*(rho12*rho13 - rho23)*z2*z3
      quad = quad/det
      ld = -1.5_dp*log2pi - log(sd1*sd2*sd3) - 0.5_dp*log(det) - 0.5_dp*quad
   end function trivariate_normal_logpdf

   elemental real(dp) function trivariate_normal_pdf(x1, x2, x3, mean1, mean2, mean3, &
      sd1, sd2, sd3, rho12, rho13, rho23) result(d)
      real(dp), intent(in) :: x1, x2, x3, mean1, mean2, mean3, sd1, sd2, sd3
      real(dp), intent(in) :: rho12, rho13, rho23
      real(dp) :: ld
      ld = trivariate_normal_logpdf(x1, x2, x3, mean1, mean2, mean3, sd1, sd2, sd3, rho12, rho13, rho23)
      d = merge(exp(ld), 0.0_dp, ld > -700.0_dp)
   end function trivariate_normal_pdf

   subroutine random_trivariate_normal(mean, sd, rho12, rho13, rho23, x, status)
      real(dp), intent(in) :: mean(3), sd(3), rho12, rho13, rho23
      real(dp), intent(out) :: x(3)
      integer, intent(out), optional :: status
      real(dp) :: z(3), l22, l32, l33
      if (minval(sd) <= 0.0_dp .or. correlation_determinant(rho12, rho13, rho23) <= 0.0_dp .or. &
          max(abs(rho12), abs(rho13), abs(rho23)) >= 1.0_dp) then
         x = huge(1.0_dp); if (present(status)) status = 1; return
      end if
      z = [rnorm_v(0.0_dp, 1.0_dp), rnorm_v(0.0_dp, 1.0_dp), rnorm_v(0.0_dp, 1.0_dp)]
      l22 = sqrt(1.0_dp - rho12*rho12)
      l32 = (rho23 - rho12*rho13)/l22
      l33 = sqrt(max(0.0_dp, 1.0_dp - rho13*rho13 - l32*l32))
      x(1) = mean(1) + sd(1)*z(1)
      x(2) = mean(2) + sd(2)*(rho12*z(1) + l22*z(2))
      x(3) = mean(3) + sd(3)*(rho13*z(1) + l32*z(2) + l33*z(3))
      if (present(status)) status = 0
   end subroutine random_trivariate_normal

   subroutine fit_trivariate_normal(y, result, weights, max_iter, tol)
      real(dp), intent(in) :: y(:, :)
      type(trivariate_normal_result_t), intent(out) :: result
      real(dp), intent(in), optional :: weights(:), tol
      integer, intent(in), optional :: max_iter
      real(dp), allocatable :: w(:), par(:), h(:, :), cov(:, :)
      real(dp) :: fval, tolerance, m(3), s(3), r12, r13, r23
      integer :: n, j, stat, stat2, niter
      n = size(y, 1)
      if (n < 4 .or. size(y, 2) /= 3) then; result%status = 1; return; end if
      allocate(w(n)); w = 1.0_dp
      if (present(weights)) then
         if (size(weights) /= n .or. any(weights < 0.0_dp)) then; result%status = 2; return; end if
         w = weights
      end if
      do j = 1, 3
         m(j) = sum(w*y(:, j))/max(sum(w), tiny(1.0_dp))
         s(j) = sqrt(max(sum(w*(y(:, j) - m(j))**2)/max(sum(w), 1.0_dp), 1.0e-8_dp))
      end do
      r12 = weighted_corr(y(:, 1), y(:, 2), w); r13 = weighted_corr(y(:, 1), y(:, 3), w)
      r23 = weighted_corr(y(:, 2), y(:, 3), w)
      call shrink_to_positive_definite(r12, r13, r23)
      allocate(par(9)); par(1:3) = m; par(4:6) = log(s)
      par(7:9) = atanh([r12, r13, r23])
      niter = 400; if (present(max_iter)) niter = max_iter
      tolerance = 1.0e-7_dp; if (present(tol)) tolerance = tol
      call bfgs_minimize(objective, par, fval, stat, max_iter=niter, tol=tolerance)
      result%status = stat; result%converged = stat == 0; result%loglik = -fval
      result%aic = 2.0_dp*fval + 18.0_dp; result%mean = par(1:3); result%sd = exp(par(4:6))
      result%rho12 = tanh(par(7)); result%rho13 = tanh(par(8)); result%rho23 = tanh(par(9))
      allocate(h(9, 9)); call numerical_hessian(objective, par, h); call invert_matrix(h, cov, stat2)
      if (stat2 == 0) then; result%covariance = cov; else; allocate(result%covariance(0, 0)); end if
   contains
      real(dp) function objective(theta) result(nll)
         real(dp), intent(in) :: theta(:)
         real(dp) :: mm(3), ss(3), rr(3), det, ld
         integer :: i
         mm = theta(1:3); ss = exp(min(30.0_dp, max(-30.0_dp, theta(4:6))))
         rr = tanh(theta(7:9)); det = correlation_determinant(rr(1), rr(2), rr(3))
         if (det <= 1.0e-8_dp) then; nll = huge(1.0_dp)/100.0_dp; return; end if
         nll = 0.0_dp
         do i = 1, n
            ld = trivariate_normal_logpdf(y(i, 1), y(i, 2), y(i, 3), mm(1), mm(2), mm(3), &
               ss(1), ss(2), ss(3), rr(1), rr(2), rr(3))
            if (.not. finite_scalar(ld)) then; nll = huge(1.0_dp)/100.0_dp; return; end if
            nll = nll - w(i)*ld
         end do
      end function objective
   end subroutine fit_trivariate_normal

   elemental real(dp) function bifgm_exponential_logpdf(x1, x2, alpha) result(ld)
      real(dp), intent(in) :: x1, x2, alpha
      real(dp) :: den
      if (x1 < 0.0_dp .or. x2 < 0.0_dp .or. abs(alpha) >= 1.0_dp) then
         ld = -huge(1.0_dp); return
      end if
      den = 1.0_dp + alpha - 2.0_dp*alpha*(exp(-x1) + exp(-x2)) + 4.0_dp*alpha*exp(-x1 - x2)
      if (den <= 0.0_dp) then; ld = -huge(1.0_dp); else; ld = -x1 - x2 + log(den); end if
   end function bifgm_exponential_logpdf

   elemental real(dp) function bifgm_exponential_pdf(x1, x2, alpha) result(d)
      real(dp), intent(in) :: x1, x2, alpha
      real(dp) :: ld
      ld = bifgm_exponential_logpdf(x1, x2, alpha); d = merge(exp(ld), 0.0_dp, ld > -700.0_dp)
   end function bifgm_exponential_pdf

   elemental real(dp) function bifgm_exponential_cdf(x1, x2, alpha) result(p)
      real(dp), intent(in) :: x1, x2, alpha
      real(dp) :: f1, f2
      if (x1 <= 0.0_dp .or. x2 <= 0.0_dp) then; p = 0.0_dp; return; end if
      if (abs(alpha) >= 1.0_dp) then; p = -1.0_dp; return; end if
      f1 = 1.0_dp - exp(-x1); f2 = 1.0_dp - exp(-x2)
      p = f1*f2*(1.0_dp + alpha*(1.0_dp - f1)*(1.0_dp - f2))
   end function bifgm_exponential_cdf

   subroutine random_bifgm_exponential(alpha, x1, x2)
      real(dp), intent(in) :: alpha
      real(dp), intent(out) :: x1, x2
      real(dp) :: u, v
      call random_fgm_copula(alpha, u, v)
      if (u <= 0.0_dp .or. v <= 0.0_dp) then; x1 = huge(1.0_dp); x2 = huge(1.0_dp); return; end if
      x1 = -log(1.0_dp - u); x2 = -log(1.0_dp - v)
   end subroutine random_bifgm_exponential

   subroutine fit_bifgm_exponential(y1, y2, x, result, weights, max_iter, tol)
      real(dp), intent(in) :: y1(:), y2(:), x(:, :)
      type(bifgm_exponential_result_t), intent(out) :: result
      real(dp), intent(in), optional :: weights(:), tol
      integer, intent(in), optional :: max_iter
      real(dp), allocatable :: w(:), par(:), h(:, :), cov(:, :)
      real(dp) :: fval, tolerance, tau0
      integer :: n, p, stat, stat2, niter
      n = size(y1); p = size(x, 2)
      if (n <= 1 .or. size(y2) /= n .or. size(x, 1) /= n .or. p <= 0 .or. any(y1 < 0.0_dp) .or. &
          any(y2 < 0.0_dp)) then; result%status = 1; return; end if
      allocate(w(n)); w = 1.0_dp
      if (present(weights)) then
         if (size(weights) /= n .or. any(weights < 0.0_dp)) then; result%status = 2; return; end if
         w = weights
      end if
      allocate(par(p)); par = 0.0_dp
      if (is_intercept(x)) then
         tau0 = kendall_tau_v(y1, y2); par(1) = atanh(max(-0.85_dp, min(0.85_dp, 4.5_dp*tau0)))
      end if
      niter = 300; if (present(max_iter)) niter = max_iter
      tolerance = 1.0e-7_dp; if (present(tol)) tolerance = tol
      call bfgs_minimize(objective, par, fval, stat, max_iter=niter, tol=tolerance)
      result%status = stat; result%converged = stat == 0; result%loglik = -fval
      result%aic = 2.0_dp*fval + 2.0_dp*real(p, dp); result%coefficients = par
      result%fitted_alpha = tanh(matmul(x, par))
      allocate(h(p, p)); call numerical_hessian(objective, par, h); call invert_matrix(h, cov, stat2)
      if (stat2 == 0) then; result%covariance = cov; else; allocate(result%covariance(0, 0)); end if
   contains
      real(dp) function objective(theta) result(nll)
         real(dp), intent(in) :: theta(:)
         real(dp) :: a, ld
         integer :: i
         nll = 0.0_dp
         do i = 1, n
            a = tanh(dot_product(x(i, :), theta)); ld = bifgm_exponential_logpdf(y1(i), y2(i), a)
            if (.not. finite_scalar(ld)) then; nll = huge(1.0_dp)/100.0_dp; return; end if
            nll = nll - w(i)*ld
         end do
      end function objective
   end subroutine fit_bifgm_exponential

   real(dp) function kendall_tau_v(x, y) result(tau)
      real(dp), intent(in) :: x(:), y(:)
      real(dp) :: con, dis, prod
      integer :: n, i, j
      n = size(x)
      if (size(y) /= n .or. n < 2) then; tau = 0.0_dp; return; end if
      con = 0.0_dp; dis = 0.0_dp
      do i = 1, n - 1
         do j = i + 1, n
            prod = (x(i) - x(j))*(y(i) - y(j))
            if (prod > 0.0_dp) then; con = con + 1.0_dp
            else if (prod < 0.0_dp) then; dis = dis + 1.0_dp
            else; con = con + 0.5_dp; dis = dis + 0.5_dp; end if
         end do
      end do
      if (con + dis <= 0.0_dp) then; tau = 0.0_dp; else; tau = (con - dis)/(con + dis); end if
   end function kendall_tau_v

   elemental real(dp) function correlation_determinant(r12, r13, r23) result(det)
      real(dp), intent(in) :: r12, r13, r23
      det = 1.0_dp + 2.0_dp*r12*r13*r23 - r12*r12 - r13*r13 - r23*r23
   end function correlation_determinant

   subroutine shrink_to_positive_definite(r12, r13, r23)
      real(dp), intent(inout) :: r12, r13, r23
      integer :: k
      r12 = max(-0.98_dp, min(0.98_dp, r12)); r13 = max(-0.98_dp, min(0.98_dp, r13))
      r23 = max(-0.98_dp, min(0.98_dp, r23))
      do k = 1, 100
         if (correlation_determinant(r12, r13, r23) > 1.0e-4_dp) exit
         r12 = 0.97_dp*r12; r13 = 0.97_dp*r13; r23 = 0.97_dp*r23
      end do
   end subroutine shrink_to_positive_definite

   real(dp) function weighted_corr(x, y, w) result(r)
      real(dp), intent(in) :: x(:), y(:), w(:)
      real(dp) :: mx, my, vx, vy, cv, sw
      sw = max(sum(w), tiny(1.0_dp)); mx = sum(w*x)/sw; my = sum(w*y)/sw
      vx = sum(w*(x - mx)**2); vy = sum(w*(y - my)**2); cv = sum(w*(x - mx)*(y - my))
      if (vx <= 0.0_dp .or. vy <= 0.0_dp) then; r = 0.0_dp; else; r = cv/sqrt(vx*vy); end if
      r = max(-0.98_dp, min(0.98_dp, r))
   end function weighted_corr

   logical function is_intercept(x) result(ok)
      real(dp), intent(in) :: x(:, :)
      ok = size(x, 2) >= 1 .and. all(abs(x(:, 1) - 1.0_dp) < 100.0_dp*epsilon(1.0_dp))
   end function is_intercept

   elemental logical function finite_scalar(x) result(ok)
      real(dp), intent(in) :: x
      ok = abs(x) <= huge(x)
   end function finite_scalar
end module vgam_multivariate_extra
