! Copyright (C) 1998-2025 T. W. Yee, University of Auckland.
! Modern Fortran computational translation, 2026.
! SPDX-License-Identifier: GPL-3.0-only
module vgam_gaitd_nb_dispersion
   use vgam_kinds, only : dp
   use vgam_distributions, only : dnbinom_v
   use vgam_optim, only : bfgs_minimize, numerical_hessian
   use vgam_linalg, only : invert_matrix
   implicit none
   private

   type, public :: gaitd_mix_nb_dispersion_result_t
      real(dp), allocatable :: parent_coefficients(:)
      real(dp), allocatable :: mass_coefficients(:, :)
      real(dp), allocatable :: outer_mean_coefficients(:, :)
      real(dp), allocatable :: size_coefficients(:, :)
      real(dp), allocatable :: covariance(:, :)
      real(dp), allocatable :: fitted_mean(:)
      real(dp), allocatable :: fitted_mass(:, :)
      real(dp), allocatable :: fitted_outer_mean(:, :)
      real(dp), allocatable :: fitted_size(:, :)
      integer, allocatable :: a_mix(:), i_mix(:), d_mix(:), truncate(:)
      real(dp) :: loglik = -huge(1.0_dp)
      real(dp) :: aic = huge(1.0_dp)
      integer :: status = 0
      logical :: converged = .false.
   contains
      procedure :: predict => predict_gaitd_mix_nb_dispersion
   end type gaitd_mix_nb_dispersion_result_t

   public :: fit_gaitd_mix_nb_dispersion_regression

contains

   subroutine fit_gaitd_mix_nb_dispersion_regression(y, xp, xmass, xouter, xsize, result, &
      a_mix, i_mix, d_mix, truncate, weights, max_iter, tol)
      integer, intent(in) :: y(:)
      real(dp), intent(in) :: xp(:, :), xmass(:, :), xouter(:, :), xsize(:, :)
      type(gaitd_mix_nb_dispersion_result_t), intent(out) :: result
      integer, intent(in), optional :: a_mix(:), i_mix(:), d_mix(:), truncate(:), max_iter
      real(dp), intent(in), optional :: weights(:), tol
      integer, allocatable :: aa(:), ii(:), dd(:), tr(:)
      real(dp), allocatable :: w(:), par(:), hess(:, :), cov(:, :)
      real(dp) :: fval, tolerance, mean_y, mass(3), om(3), sizes(4), mu0
      integer :: n, pp, pm, po, ps, np, niter, stat, stat2, r, idx

      call copy_int(a_mix, aa); call copy_int(i_mix, ii)
      call copy_int(d_mix, dd); call copy_int(truncate, tr)
      n = size(y); pp = size(xp, 2); pm = size(xmass, 2); po = size(xouter, 2); ps = size(xsize, 2)
      if (n <= 0 .or. min(pp, pm, po, ps) <= 0 .or. any(y < 0) .or. size(xp, 1) /= n .or. &
          size(xmass, 1) /= n .or. size(xouter, 1) /= n .or. size(xsize, 1) /= n) then
         result%status = 1; return
      end if
      if (size(aa) + size(ii) + size(dd) == 0) then
         result%status = 2; return
      end if
      if (any([aa, ii, dd, tr] < 0) .or. any_duplicate([aa, ii, dd, tr])) then
         result%status = 3; return
      end if
      allocate(w(n)); w = 1.0_dp
      if (present(weights)) then
         if (size(weights) /= n .or. any(weights < 0.0_dp)) then
            result%status = 4; return
         end if
         w = weights
      end if

      np = pp + 3*pm + 3*po + 4*ps
      allocate(par(np)); par = 0.0_dp
      mean_y = sum(w*real(y, dp))/max(sum(w), tiny(1.0_dp))
      if (is_intercept(xp)) par(1) = log(max(mean_y, 0.2_dp))
      idx = pp + 3*pm
      if (is_intercept(xouter)) then
         par(idx + 1) = log(max(mean_y, 0.2_dp))
         par(idx + po + 1) = log(max(mean_y, 0.2_dp))
         par(idx + 2*po + 1) = log(max(mean_y, 0.2_dp))
      end if
      idx = pp + 3*pm + 3*po
      if (is_intercept(xsize)) then
         par(idx + 1) = log(2.0_dp)
         par(idx + ps + 1) = log(2.0_dp)
         par(idx + 2*ps + 1) = log(2.0_dp)
         par(idx + 3*ps + 1) = log(2.0_dp)
      end if
      niter = 600; if (present(max_iter)) niter = max_iter
      tolerance = 2.0e-7_dp; if (present(tol)) tolerance = tol
      call bfgs_minimize(objective, par, fval, stat, max_iter=niter, tol=tolerance)

      result%status = stat; result%converged = stat == 0
      result%parent_coefficients = par(1:pp)
      result%mass_coefficients = reshape(par(pp + 1:pp + 3*pm), [pm, 3])
      idx = pp + 3*pm
      result%outer_mean_coefficients = reshape(par(idx + 1:idx + 3*po), [po, 3])
      idx = idx + 3*po
      result%size_coefficients = reshape(par(idx + 1:idx + 4*ps), [ps, 4])
      result%a_mix = aa; result%i_mix = ii; result%d_mix = dd; result%truncate = tr
      result%loglik = -fval; result%aic = 2.0_dp*fval + 2.0_dp*real(np, dp)
      allocate(result%fitted_mean(n), result%fitted_mass(n, 3), result%fitted_outer_mean(n, 3))
      allocate(result%fitted_size(n, 4))
      do r = 1, n
         mu0 = exp(clamp_eta(dot_product(xp(r, :), result%parent_coefficients)))
         call mass_probabilities(xmass(r, :), result%mass_coefficients, aa, ii, dd, mass)
         om = exp(clamp_eta(matmul(transpose(result%outer_mean_coefficients), xouter(r, :))))
         sizes = exp(clamp_eta(matmul(transpose(result%size_coefficients), xsize(r, :))))
         result%fitted_mass(r, :) = mass; result%fitted_outer_mean(r, :) = om
         result%fitted_size(r, :) = sizes
         result%fitted_mean(r) = mix_mean(mu0, sizes, aa, ii, dd, tr, mass, om)
      end do
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
         real(dp) :: beta(pp), gm(pm, 3), bo(po, 3), bs(ps, 4)
         real(dp) :: mu, mm(3), outer_mu(3), sz(4), py
         integer :: row, k
         beta = theta(1:pp)
         gm = reshape(theta(pp + 1:pp + 3*pm), [pm, 3])
         k = pp + 3*pm
         bo = reshape(theta(k + 1:k + 3*po), [po, 3]); k = k + 3*po
         bs = reshape(theta(k + 1:k + 4*ps), [ps, 4])
         nll = 0.0_dp
         do row = 1, n
            mu = exp(clamp_eta(dot_product(xp(row, :), beta)))
            call mass_probabilities(xmass(row, :), gm, aa, ii, dd, mm)
            outer_mu = exp(clamp_eta(matmul(transpose(bo), xouter(row, :))))
            sz = exp(clamp_eta(matmul(transpose(bs), xsize(row, :))))
            py = mix_probability(y(row), mu, sz, aa, ii, dd, tr, mm, outer_mu)
            if (py <= tiny(1.0_dp) .or. .not. finite_scalar(py)) then
               nll = huge(1.0_dp)/100.0_dp; return
            end if
            nll = nll - w(row)*log(py)
         end do
      end function objective
   end subroutine fit_gaitd_mix_nb_dispersion_regression

   subroutine predict_gaitd_mix_nb_dispersion(self, xp, xmass, xouter, xsize, fitted_mean, &
      fitted_mass, fitted_outer_mean, fitted_size)
      class(gaitd_mix_nb_dispersion_result_t), intent(in) :: self
      real(dp), intent(in) :: xp(:, :), xmass(:, :), xouter(:, :), xsize(:, :)
      real(dp), allocatable, intent(out) :: fitted_mean(:)
      real(dp), allocatable, intent(out), optional :: fitted_mass(:, :), fitted_outer_mean(:, :), fitted_size(:, :)
      real(dp) :: mu, mass(3), om(3), sz(4)
      integer :: n, r
      n = size(xp, 1)
      if (size(xmass, 1) /= n .or. size(xouter, 1) /= n .or. size(xsize, 1) /= n .or. &
          size(xp, 2) /= size(self%parent_coefficients) .or. &
          size(xmass, 2) /= size(self%mass_coefficients, 1) .or. &
          size(xouter, 2) /= size(self%outer_mean_coefficients, 1) .or. &
          size(xsize, 2) /= size(self%size_coefficients, 1)) then
         allocate(fitted_mean(0)); if (present(fitted_mass)) allocate(fitted_mass(0, 0))
         if (present(fitted_outer_mean)) allocate(fitted_outer_mean(0, 0))
         if (present(fitted_size)) allocate(fitted_size(0, 0)); return
      end if
      allocate(fitted_mean(n)); if (present(fitted_mass)) allocate(fitted_mass(n, 3))
      if (present(fitted_outer_mean)) allocate(fitted_outer_mean(n, 3))
      if (present(fitted_size)) allocate(fitted_size(n, 4))
      do r = 1, n
         mu = exp(clamp_eta(dot_product(xp(r, :), self%parent_coefficients)))
         call mass_probabilities(xmass(r, :), self%mass_coefficients, self%a_mix, self%i_mix, self%d_mix, mass)
         om = exp(clamp_eta(matmul(transpose(self%outer_mean_coefficients), xouter(r, :))))
         sz = exp(clamp_eta(matmul(transpose(self%size_coefficients), xsize(r, :))))
         fitted_mean(r) = mix_mean(mu, sz, self%a_mix, self%i_mix, self%d_mix, self%truncate, mass, om)
         if (present(fitted_mass)) fitted_mass(r, :) = mass
         if (present(fitted_outer_mean)) fitted_outer_mean(r, :) = om
         if (present(fitted_size)) fitted_size(r, :) = sz
      end do
   end subroutine predict_gaitd_mix_nb_dispersion

   subroutine mass_probabilities(xrow, gamma, aa, ii, dd, mass)
      real(dp), intent(in) :: xrow(:), gamma(:, :)
      integer, intent(in) :: aa(:), ii(:), dd(:)
      real(dp), intent(out) :: mass(3)
      real(dp) :: eta(3), e(3), den, vmax
      logical :: active(3)
      active = [size(aa) > 0, size(ii) > 0, size(dd) > 0]
      eta = matmul(transpose(gamma), xrow); vmax = max(0.0_dp, maxval(eta, mask=active))
      e = 0.0_dp; where (active) e = exp(eta - vmax)
      den = exp(-vmax) + sum(e); mass = e/den
   end subroutine mass_probabilities

   real(dp) function mix_probability(y, mu, sizes, aa, ii, dd, tr, mass, om) result(py)
      integer, intent(in) :: y, aa(:), ii(:), dd(:), tr(:)
      real(dp), intent(in) :: mu, sizes(4), mass(3), om(3)
      real(dp) :: den, delta
      den = 1.0_dp - point_mass(mu, sizes(1), tr) - point_mass(mu, sizes(1), aa)
      if (den <= tiny(1.0_dp)) then; py = 0.0_dp; return; end if
      delta = (1.0_dp - mass(1) - mass(2) + mass(3))/den
      py = 0.0_dp
      if (.not. any(tr == y) .and. .not. any(aa == y)) py = delta*nb_pmf(y, mu, sizes(1))
      py = py + mass(1)*restricted_weight(y, aa, om(1), sizes(2))
      py = py + mass(2)*restricted_weight(y, ii, om(2), sizes(3))
      py = py - mass(3)*restricted_weight(y, dd, om(3), sizes(4))
   end function mix_probability

   real(dp) function mix_mean(mu, sizes, aa, ii, dd, tr, mass, om) result(ans)
      real(dp), intent(in) :: mu, sizes(4), mass(3), om(3)
      integer, intent(in) :: aa(:), ii(:), dd(:), tr(:)
      real(dp) :: den, delta, removed
      den = 1.0_dp - point_mass(mu, sizes(1), tr) - point_mass(mu, sizes(1), aa)
      if (den <= tiny(1.0_dp)) then; ans = huge(1.0_dp); return; end if
      delta = (1.0_dp - mass(1) - mass(2) + mass(3))/den
      removed = point_first_moment(mu, sizes(1), tr) + point_first_moment(mu, sizes(1), aa)
      ans = delta*(mu - removed)
      ans = ans + mass(1)*restricted_mean(aa, om(1), sizes(2))
      ans = ans + mass(2)*restricted_mean(ii, om(2), sizes(3))
      ans = ans - mass(3)*restricted_mean(dd, om(3), sizes(4))
   end function mix_mean

   real(dp) function restricted_weight(y, points, mu, sizev) result(w)
      integer, intent(in) :: y, points(:)
      real(dp), intent(in) :: mu, sizev
      real(dp) :: den
      if (size(points) == 0 .or. .not. any(points == y)) then; w = 0.0_dp; return; end if
      den = point_mass(mu, sizev, points)
      if (den <= tiny(1.0_dp)) then; w = 0.0_dp; else; w = nb_pmf(y, mu, sizev)/den; end if
   end function restricted_weight

   real(dp) function restricted_mean(points, mu, sizev) result(ans)
      integer, intent(in) :: points(:)
      real(dp), intent(in) :: mu, sizev
      real(dp) :: den
      if (size(points) == 0) then; ans = 0.0_dp; return; end if
      den = point_mass(mu, sizev, points)
      if (den <= tiny(1.0_dp)) then
         ans = 0.0_dp
      else
         ans = point_first_moment(mu, sizev, points)/den
      end if
   end function restricted_mean

   real(dp) function point_mass(mu, sizev, points) result(ans)
      real(dp), intent(in) :: mu, sizev
      integer, intent(in) :: points(:)
      integer :: j
      ans = 0.0_dp
      do j = 1, size(points); ans = ans + nb_pmf(points(j), mu, sizev); end do
   end function point_mass

   real(dp) function point_first_moment(mu, sizev, points) result(ans)
      real(dp), intent(in) :: mu, sizev
      integer, intent(in) :: points(:)
      integer :: j
      ans = 0.0_dp
      do j = 1, size(points)
         ans = ans + real(points(j), dp)*nb_pmf(points(j), mu, sizev)
      end do
   end function point_first_moment

   elemental real(dp) function nb_pmf(y, mu, sizev) result(p)
      integer, intent(in) :: y
      real(dp), intent(in) :: mu, sizev
      real(dp) :: prob
      if (y < 0 .or. mu <= 0.0_dp .or. sizev <= 0.0_dp) then
         p = 0.0_dp
      else
         prob = sizev/(sizev + mu); p = dnbinom_v(y, sizev, prob)
      end if
   end function nb_pmf

   subroutine copy_int(src, dst)
      integer, intent(in), optional :: src(:)
      integer, allocatable, intent(out) :: dst(:)
      if (present(src)) then; dst = src; else; allocate(dst(0)); end if
   end subroutine copy_int

   logical function any_duplicate(x) result(ans)
      integer, intent(in) :: x(:)
      integer :: j, k
      ans = .false.
      do j = 1, size(x); do k = j + 1, size(x)
         if (x(j) == x(k)) then; ans = .true.; return; end if
      end do; end do
   end function any_duplicate

   logical function is_intercept(x) result(ok)
      real(dp), intent(in) :: x(:, :)
      ok = size(x, 2) >= 1 .and. all(abs(x(:, 1) - 1.0_dp) < 100.0_dp*epsilon(1.0_dp))
   end function is_intercept

   elemental real(dp) function clamp_eta(x) result(y)
      real(dp), intent(in) :: x
      y = min(30.0_dp, max(-30.0_dp, x))
   end function clamp_eta

   elemental logical function finite_scalar(x) result(ok)
      real(dp), intent(in) :: x
      ok = abs(x) <= huge(x)
   end function finite_scalar
end module vgam_gaitd_nb_dispersion
