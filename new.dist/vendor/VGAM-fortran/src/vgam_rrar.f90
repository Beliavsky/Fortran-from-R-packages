! Copyright (C) 1998-2025 T. W. Yee, University of Auckland.
! Modern Fortran computational translation, 2026.
! SPDX-License-Identifier: GPL-3.0-only
module vgam_rrar
   use vgam_kinds, only : dp, log2pi
   use vgam_linalg, only : weighted_least_squares, cholesky_factor, invert_matrix
   use vgam_optim, only : bfgs_minimize, numerical_hessian
   implicit none
   private

   type, public :: rrar_result_t
      integer, allocatable :: ranks(:)
      real(dp), allocatable :: common_basis(:, :)
      real(dp), allocatable :: right_factors(:, :, :)
      real(dp), allocatable :: phi(:, :, :)
      real(dp), allocatable :: omega(:, :)
      real(dp), allocatable :: fitted(:, :)
      real(dp), allocatable :: residuals(:, :)
      real(dp), allocatable :: coefficients(:)
      real(dp), allocatable :: covariance(:, :)
      real(dp), allocatable :: transformed_series(:, :)
      real(dp) :: loglik = -huge(1.0_dp)
      real(dp) :: aic = huge(1.0_dp)
      integer :: lag = 0
      integer :: rank_max = 0
      integer :: status = 0
      logical :: converged = .false.
   contains
      procedure :: forecast => forecast_rrar
   end type rrar_result_t

   public :: fit_rrar

contains

   subroutine fit_rrar(y, ranks, result, max_iter, tol, compute_covariance)
      real(dp), intent(in) :: y(:, :)
      integer, intent(in) :: ranks(:)
      type(rrar_result_t), intent(out) :: result
      integer, intent(in), optional :: max_iter
      real(dp), intent(in), optional :: tol
      logical, intent(in), optional :: compute_covariance
      real(dp), allocatable :: xlag(:, :), bvar(:, :), beta(:), cov(:, :), ones(:)
      real(dp), allocatable :: a(:, :), c(:, :, :), phi(:, :, :), theta(:)
      real(dp), allocatable :: resid(:, :), omega(:, :), fitted(:, :), hess(:, :), cinv(:, :)
      real(dp) :: fval, tolerance
      integer :: n, m, lags, rmax, neff, pvar, j, l, r, pos, npar, stat, stat2
      integer :: niter, nfree_a, t, k
      logical :: do_cov

      n = size(y, 1)
      m = size(y, 2)
      lags = size(ranks)
      if (n <= lags .or. m <= 1 .or. lags <= 0) then
         result%status = 1
         return
      end if
      if (any(ranks < 1) .or. any(ranks > m)) then
         result%status = 2
         return
      end if
      if (lags > 1) then
         if (any(ranks(2:) > ranks(:lags - 1))) then
            result%status = 3
            return
         end if
      end if
      rmax = ranks(1)
      neff = n - lags
      pvar = m*lags
      allocate(xlag(neff, pvar), bvar(pvar, m), ones(neff))
      ones = 1.0_dp
      do t = 1, neff
         pos = 0
         do l = 1, lags
            xlag(t, pos + 1:pos + m) = y(lags + t - l, :)
            pos = pos + m
         end do
      end do
      do j = 1, m
         call weighted_least_squares(xlag, y(lags + 1:n, j), ones, beta, cov, stat)
         if (stat /= 0) then
            result%status = 10 + stat
            return
         end if
         bvar(:, j) = beta
      end do

      allocate(a(m, rmax), c(m, rmax, lags), phi(m, m, lags))
      a = 0.0_dp
      do j = 1, rmax
         a(j, j) = 1.0_dp
      end do
      c = 0.0_dp
      pos = 0
      do l = 1, lags
         r = ranks(l)
         ! With the corner identification A(1:rmax,:) = I, the first r
         ! rows of Phi supply a natural right-factor starting value.
         c(:, 1:r, l) = bvar(pos + 1:pos + m, 1:r)
         pos = pos + m
      end do

      nfree_a = max(0, m - rmax)*rmax
      npar = nfree_a + m*sum(ranks)
      allocate(theta(npar))
      pos = 0
      if (nfree_a > 0) then
         theta(1:nfree_a) = reshape(a(rmax + 1:m, :), [nfree_a])
         pos = nfree_a
      end if
      do l = 1, lags
         r = ranks(l)
         theta(pos + 1:pos + m*r) = reshape(c(:, 1:r, l), [m*r])
         pos = pos + m*r
      end do

      niter = 300
      if (present(max_iter)) niter = max_iter
      tolerance = 1.0e-7_dp
      if (present(tol)) tolerance = tol
      call bfgs_minimize(objective, theta, fval, stat, max_iter=niter, tol=tolerance)
      result%status = stat
      result%converged = stat == 0
      call unpack(theta, a, c, phi)
      call residual_matrices(phi, resid, fitted, omega)
      result%ranks = ranks
      result%common_basis = a
      result%right_factors = c
      result%phi = phi
      result%omega = omega
      result%fitted = fitted
      result%residuals = resid
      result%coefficients = theta
      result%lag = lags
      result%rank_max = rmax
      result%loglik = -concentrated_nll(omega, neff, m)
      result%aic = -2.0_dp*result%loglik + 2.0_dp*real(npar + m*(m + 1)/2, dp)

      do_cov = npar <= 60
      if (present(compute_covariance)) do_cov = compute_covariance
      if (do_cov) then
         allocate(hess(npar, npar))
         call numerical_hessian(objective, theta, hess)
         call invert_matrix(hess, cov, stat2)
         if (stat2 == 0) then
            result%covariance = cov
         else
            allocate(result%covariance(0, 0))
         end if
      else
         allocate(result%covariance(0, 0))
      end if

      if (rmax == m) then
         call invert_matrix(a, cinv, stat2)
         if (stat2 == 0) result%transformed_series = matmul(y, transpose(cinv))
      end if

   contains

      real(dp) function objective(par) result(val)
         real(dp), intent(in) :: par(:)
         real(dp) :: aa(m, rmax), cc(m, rmax, lags), pp(m, m, lags)
         real(dp), allocatable :: rr(:, :), ff(:, :), oo(:, :)
         call unpack(par, aa, cc, pp)
         call residual_matrices(pp, rr, ff, oo)
         val = concentrated_nll(oo, neff, m)
         val = val + 1.0e-12_dp*sum(par*par)
         if (.not. finite_scalar(val)) val = huge(1.0_dp)/100.0_dp
      end function objective

      subroutine unpack(par, aa, cc, pp)
         real(dp), intent(in) :: par(:)
         real(dp), intent(out) :: aa(:, :), cc(:, :, :), pp(:, :, :)
         integer :: ll, rr, p0, jj
         aa = 0.0_dp
         do jj = 1, rmax
            aa(jj, jj) = 1.0_dp
         end do
         cc = 0.0_dp
         p0 = 0
         if (nfree_a > 0) then
            aa(rmax + 1:m, :) = reshape(par(1:nfree_a), [m - rmax, rmax])
            p0 = nfree_a
         end if
         do ll = 1, lags
            rr = ranks(ll)
            cc(:, 1:rr, ll) = reshape(par(p0 + 1:p0 + m*rr), [m, rr])
            pp(:, :, ll) = matmul(aa(:, 1:rr), transpose(cc(:, 1:rr, ll)))
            p0 = p0 + m*rr
         end do
      end subroutine unpack

      subroutine residual_matrices(pp, rr, ff, oo)
         real(dp), intent(in) :: pp(:, :, :)
         real(dp), allocatable, intent(out) :: rr(:, :), ff(:, :), oo(:, :)
         integer :: tt, ll
         allocate(rr(neff, m), ff(neff, m), oo(m, m))
         ff = 0.0_dp
         do tt = 1, neff
            do ll = 1, lags
               ff(tt, :) = ff(tt, :) + matmul(pp(:, :, ll), y(lags + tt - ll, :))
            end do
         end do
         rr = y(lags + 1:n, :) - ff
         oo = matmul(transpose(rr), rr)/real(neff, dp)
         do k = 1, m
            oo(k, k) = oo(k, k) + 1.0e-10_dp*max(1.0_dp, sum(abs(oo(k, :))))
         end do
      end subroutine residual_matrices

   end subroutine fit_rrar

   subroutine forecast_rrar(self, history, steps, forecast)
      class(rrar_result_t), intent(in) :: self
      real(dp), intent(in) :: history(:, :)
      integer, intent(in) :: steps
      real(dp), allocatable, intent(out) :: forecast(:, :)
      real(dp), allocatable :: work(:, :)
      integer :: n, m, h, l, row

      n = size(history, 1)
      m = size(history, 2)
      if (.not. allocated(self%phi) .or. steps < 0 .or. m /= size(self%phi, 1) .or. n < self%lag) then
         allocate(forecast(0, 0))
         return
      end if
      allocate(work(n + steps, m), forecast(steps, m))
      work(1:n, :) = history
      do h = 1, steps
         row = n + h
         work(row, :) = 0.0_dp
         do l = 1, self%lag
            work(row, :) = work(row, :) + matmul(self%phi(:, :, l), work(row - l, :))
         end do
         forecast(h, :) = work(row, :)
      end do
   end subroutine forecast_rrar

   real(dp) function concentrated_nll(omega, nobs, m) result(v)
      real(dp), intent(in) :: omega(:, :)
      integer, intent(in) :: nobs, m
      real(dp), allocatable :: l(:, :)
      real(dp) :: logdet
      integer :: stat, j
      call cholesky_factor(omega, l, stat)
      if (stat /= 0) then
         v = huge(1.0_dp)/100.0_dp
         return
      end if
      logdet = 0.0_dp
      do j = 1, m
         logdet = logdet + 2.0_dp*log(max(l(j, j), tiny(1.0_dp)))
      end do
      v = 0.5_dp*real(nobs, dp)*(real(m, dp)*(1.0_dp + log2pi) + logdet)
   end function concentrated_nll

   elemental logical function finite_scalar(x) result(ok)
      real(dp), intent(in) :: x
      ok = abs(x) <= huge(x)
   end function finite_scalar

end module vgam_rrar
