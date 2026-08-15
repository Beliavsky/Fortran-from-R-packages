! Copyright (C) 1998-2025 T. W. Yee, University of Auckland.
! Modern Fortran computational translation, 2026.
! SPDX-License-Identifier: GPL-3.0-only
module vgam_gaitd_mix_regression
   use vgam_kinds, only : dp
   use vgam_distributions, only : dpois_v, dnbinom_v
   use vgam_optim, only : bfgs_minimize, numerical_hessian
   use vgam_linalg, only : invert_matrix
   implicit none
   private

   integer, parameter, public :: gaitd_mix_base_poisson = 1
   integer, parameter, public :: gaitd_mix_base_negbinomial = 2

   type, public :: gaitd_mix_regression_result_t
      real(dp), allocatable :: parent_coefficients(:)
      real(dp), allocatable :: mass_coefficients(:, :)
      real(dp), allocatable :: outer_coefficients(:, :)
      real(dp), allocatable :: covariance(:, :)
      real(dp), allocatable :: fitted_mean(:)
      real(dp), allocatable :: fitted_mass(:, :)
      real(dp), allocatable :: fitted_outer_mean(:, :)
      integer, allocatable :: a_mix(:), i_mix(:), d_mix(:), truncate(:)
      real(dp) :: size = huge(1.0_dp)
      real(dp) :: loglik = -huge(1.0_dp)
      real(dp) :: aic = huge(1.0_dp)
      integer :: base_family = gaitd_mix_base_poisson
      integer :: status = 0
      logical :: converged = .false.
   contains
      procedure :: predict => predict_gaitd_mix_regression
   end type gaitd_mix_regression_result_t

   public :: fit_gaitd_mix_poisson_regression, fit_gaitd_mix_nb_regression

contains

   subroutine fit_gaitd_mix_poisson_regression(y, x_parent, x_mass, x_outer, result, &
      a_mix, i_mix, d_mix, truncate, weights, max_iter, tol)
      integer, intent(in) :: y(:)
      real(dp), intent(in) :: x_parent(:, :), x_mass(:, :), x_outer(:, :)
      type(gaitd_mix_regression_result_t), intent(out) :: result
      integer, intent(in), optional :: a_mix(:), i_mix(:), d_mix(:), truncate(:), max_iter
      real(dp), intent(in), optional :: weights(:), tol
      call fit_mix_common(y, x_parent, x_mass, x_outer, gaitd_mix_base_poisson, result, &
         a_mix, i_mix, d_mix, truncate, weights, max_iter, tol)
   end subroutine fit_gaitd_mix_poisson_regression

   subroutine fit_gaitd_mix_nb_regression(y, x_parent, x_mass, x_outer, result, &
      a_mix, i_mix, d_mix, truncate, weights, max_iter, tol)
      integer, intent(in) :: y(:)
      real(dp), intent(in) :: x_parent(:, :), x_mass(:, :), x_outer(:, :)
      type(gaitd_mix_regression_result_t), intent(out) :: result
      integer, intent(in), optional :: a_mix(:), i_mix(:), d_mix(:), truncate(:), max_iter
      real(dp), intent(in), optional :: weights(:), tol
      call fit_mix_common(y, x_parent, x_mass, x_outer, gaitd_mix_base_negbinomial, result, &
         a_mix, i_mix, d_mix, truncate, weights, max_iter, tol)
   end subroutine fit_gaitd_mix_nb_regression

   subroutine fit_mix_common(y, xp, xmass, xouter, family, result, &
      a_mix_in, i_mix_in, d_mix_in, truncate_in, weights, max_iter, tol)
      integer, intent(in) :: y(:), family
      real(dp), intent(in) :: xp(:, :), xmass(:, :), xouter(:, :)
      type(gaitd_mix_regression_result_t), intent(out) :: result
      integer, intent(in), optional :: a_mix_in(:), i_mix_in(:), d_mix_in(:), truncate_in(:), max_iter
      real(dp), intent(in), optional :: weights(:), tol
      integer, allocatable :: aa(:), ii(:), dd(:), tr(:)
      real(dp), allocatable :: w(:), par(:), hess(:, :), cov(:, :)
      real(dp) :: fval, tolerance, mean_y, muv, sizev, masses(3), om(3)
      integer :: n, pp, pm, po, np, niter, stat, stat2, row, j, idx

      call copy_int(a_mix_in, aa); call copy_int(i_mix_in, ii)
      call copy_int(d_mix_in, dd); call copy_int(truncate_in, tr)
      n = size(y); pp = size(xp, 2); pm = size(xmass, 2); po = size(xouter, 2)
      if (n <= 0 .or. pp <= 0 .or. pm <= 0 .or. po <= 0 .or. any(y < 0) .or. &
          size(xp, 1) /= n .or. size(xmass, 1) /= n .or. size(xouter, 1) /= n) then
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

      np = pp + 3*pm + 3*po + merge(1, 0, family == gaitd_mix_base_negbinomial)
      allocate(par(np)); par = 0.0_dp
      mean_y = sum(w*real(y, dp))/max(sum(w), tiny(1.0_dp))
      if (all(abs(xp(:, 1) - 1.0_dp) < 100.0_dp*epsilon(1.0_dp))) par(1) = log(max(mean_y, 0.2_dp))
      idx = pp + 3*pm
      do j = 1, 3
         if (all(abs(xouter(:, 1) - 1.0_dp) < 100.0_dp*epsilon(1.0_dp))) then
            par(idx + (j - 1)*po + 1) = log(max(mean_y, 0.2_dp))
         end if
      end do
      if (family == gaitd_mix_base_negbinomial) par(np) = log(2.0_dp)
      niter = 450; if (present(max_iter)) niter = max_iter
      tolerance = 2.0e-7_dp; if (present(tol)) tolerance = tol
      call bfgs_minimize(objective, par, fval, stat, max_iter=niter, tol=tolerance)

      result%status = stat; result%converged = stat == 0; result%base_family = family
      result%parent_coefficients = par(1:pp)
      result%mass_coefficients = reshape(par(pp + 1:pp + 3*pm), [pm, 3])
      idx = pp + 3*pm
      result%outer_coefficients = reshape(par(idx + 1:idx + 3*po), [po, 3])
      sizev = huge(1.0_dp)
      if (family == gaitd_mix_base_negbinomial) sizev = exp(min(25.0_dp, par(np)))
      result%size = sizev
      result%a_mix = aa; result%i_mix = ii; result%d_mix = dd; result%truncate = tr
      result%loglik = -fval; result%aic = 2.0_dp*fval + 2.0_dp*real(np, dp)
      allocate(result%fitted_mean(n), result%fitted_mass(n, 3), result%fitted_outer_mean(n, 3))
      do row = 1, n
         muv = exp(clamp_eta(dot_product(xp(row, :), result%parent_coefficients)))
         call mass_probabilities(xmass(row, :), result%mass_coefficients, aa, ii, dd, masses)
         om = exp(clamp_eta(matmul(transpose(result%outer_coefficients), xouter(row, :))))
         result%fitted_mass(row, :) = masses
         result%fitted_outer_mean(row, :) = om
         result%fitted_mean(row) = mix_mean(muv, sizev, family, aa, ii, dd, tr, masses, om)
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
         real(dp) :: beta(pp), gm(pm, 3), bo(po, 3), mu0, sz, mass(3), outer_mu(3), py
         integer :: r, kk
         beta = theta(1:pp)
         gm = reshape(theta(pp + 1:pp + 3*pm), [pm, 3])
         kk = pp + 3*pm
         bo = reshape(theta(kk + 1:kk + 3*po), [po, 3])
         sz = huge(1.0_dp)
         if (family == gaitd_mix_base_negbinomial) sz = exp(min(25.0_dp, theta(np)))
         nll = 0.0_dp
         do r = 1, n
            mu0 = exp(clamp_eta(dot_product(xp(r, :), beta)))
            call mass_probabilities(xmass(r, :), gm, aa, ii, dd, mass)
            outer_mu = exp(clamp_eta(matmul(transpose(bo), xouter(r, :))))
            py = mix_probability(y(r), mu0, sz, family, aa, ii, dd, tr, mass, outer_mu)
            if (py <= tiny(1.0_dp) .or. .not. finite_scalar(py)) then
               nll = huge(1.0_dp)/100.0_dp; return
            end if
            nll = nll - w(r)*log(py)
         end do
      end function objective
   end subroutine fit_mix_common

   subroutine predict_gaitd_mix_regression(self, x_parent, x_mass, x_outer, fitted_mean, &
      fitted_mass, fitted_outer_mean)
      class(gaitd_mix_regression_result_t), intent(in) :: self
      real(dp), intent(in) :: x_parent(:, :), x_mass(:, :), x_outer(:, :)
      real(dp), allocatable, intent(out) :: fitted_mean(:)
      real(dp), allocatable, intent(out), optional :: fitted_mass(:, :), fitted_outer_mean(:, :)
      real(dp) :: mu0, mass(3), om(3)
      integer :: n, r
      n = size(x_parent, 1)
      if (size(x_mass, 1) /= n .or. size(x_outer, 1) /= n .or. &
          size(x_parent, 2) /= size(self%parent_coefficients) .or. &
          size(x_mass, 2) /= size(self%mass_coefficients, 1) .or. &
          size(x_outer, 2) /= size(self%outer_coefficients, 1)) then
         allocate(fitted_mean(0)); if (present(fitted_mass)) allocate(fitted_mass(0, 0))
         if (present(fitted_outer_mean)) allocate(fitted_outer_mean(0, 0)); return
      end if
      allocate(fitted_mean(n)); if (present(fitted_mass)) allocate(fitted_mass(n, 3))
      if (present(fitted_outer_mean)) allocate(fitted_outer_mean(n, 3))
      do r = 1, n
         mu0 = exp(clamp_eta(dot_product(x_parent(r, :), self%parent_coefficients)))
         call mass_probabilities(x_mass(r, :), self%mass_coefficients, self%a_mix, self%i_mix, self%d_mix, mass)
         om = exp(clamp_eta(matmul(transpose(self%outer_coefficients), x_outer(r, :))))
         fitted_mean(r) = mix_mean(mu0, self%size, self%base_family, self%a_mix, self%i_mix, &
            self%d_mix, self%truncate, mass, om)
         if (present(fitted_mass)) fitted_mass(r, :) = mass
         if (present(fitted_outer_mean)) fitted_outer_mean(r, :) = om
      end do
   end subroutine predict_gaitd_mix_regression

   subroutine mass_probabilities(xrow, gamma, aa, ii, dd, mass)
      real(dp), intent(in) :: xrow(:), gamma(:, :)
      integer, intent(in) :: aa(:), ii(:), dd(:)
      real(dp), intent(out) :: mass(3)
      real(dp) :: eta(3), e(3), den, vmax
      logical :: active(3)
      active = [size(aa) > 0, size(ii) > 0, size(dd) > 0]
      eta = matmul(transpose(gamma), xrow)
      vmax = max(0.0_dp, maxval(eta, mask=active))
      e = 0.0_dp
      where (active) e = exp(eta - vmax)
      den = exp(-vmax) + sum(e)
      mass = e/den
   end subroutine mass_probabilities

   real(dp) function mix_probability(y, mu0, sizev, family, aa, ii, dd, tr, mass, om) result(py)
      integer, intent(in) :: y, family, aa(:), ii(:), dd(:), tr(:)
      real(dp), intent(in) :: mu0, sizev, mass(3), om(3)
      real(dp) :: den, delta, wa, wi, wd
      den = 1.0_dp - point_mass(mu0, sizev, family, tr) - point_mass(mu0, sizev, family, aa)
      if (den <= tiny(1.0_dp)) then
         py = 0.0_dp; return
      end if
      delta = (1.0_dp - mass(1) - mass(2) + mass(3))/den
      py = 0.0_dp
      if (.not. any(tr == y) .and. .not. any(aa == y)) py = delta*base_pmf(y, mu0, sizev, family)
      wa = restricted_weight(y, aa, om(1), sizev, family)
      wi = restricted_weight(y, ii, om(2), sizev, family)
      wd = restricted_weight(y, dd, om(3), sizev, family)
      py = py + mass(1)*wa + mass(2)*wi - mass(3)*wd
   end function mix_probability

   real(dp) function mix_mean(mu0, sizev, family, aa, ii, dd, tr, mass, om) result(ans)
      real(dp), intent(in) :: mu0, sizev, mass(3), om(3)
      integer, intent(in) :: family, aa(:), ii(:), dd(:), tr(:)
      real(dp) :: den, delta, removed
      den = 1.0_dp - point_mass(mu0, sizev, family, tr) - point_mass(mu0, sizev, family, aa)
      if (den <= tiny(1.0_dp)) then
         ans = huge(1.0_dp); return
      end if
      delta = (1.0_dp - mass(1) - mass(2) + mass(3))/den
      removed = point_first_moment(mu0, sizev, family, tr) + point_first_moment(mu0, sizev, family, aa)
      ans = delta*(mu0 - removed)
      ans = ans + mass(1)*restricted_mean(aa, om(1), sizev, family)
      ans = ans + mass(2)*restricted_mean(ii, om(2), sizev, family)
      ans = ans - mass(3)*restricted_mean(dd, om(3), sizev, family)
   end function mix_mean

   real(dp) function restricted_weight(y, points, mu, sizev, family) result(w)
      integer, intent(in) :: y, points(:), family
      real(dp), intent(in) :: mu, sizev
      real(dp) :: den
      if (size(points) == 0 .or. .not. any(points == y)) then
         w = 0.0_dp; return
      end if
      den = point_mass(mu, sizev, family, points)
      if (den <= tiny(1.0_dp)) then
         w = 0.0_dp
      else
         w = base_pmf(y, mu, sizev, family)/den
      end if
   end function restricted_weight

   real(dp) function restricted_mean(points, mu, sizev, family) result(ans)
      integer, intent(in) :: points(:), family
      real(dp), intent(in) :: mu, sizev
      real(dp) :: den
      if (size(points) == 0) then
         ans = 0.0_dp; return
      end if
      den = point_mass(mu, sizev, family, points)
      if (den <= tiny(1.0_dp)) then
         ans = 0.0_dp
      else
         ans = point_first_moment(mu, sizev, family, points)/den
      end if
   end function restricted_mean

   real(dp) function point_mass(mu, sizev, family, points) result(ans)
      real(dp), intent(in) :: mu, sizev
      integer, intent(in) :: family, points(:)
      integer :: j
      ans = 0.0_dp
      do j = 1, size(points)
         ans = ans + base_pmf(points(j), mu, sizev, family)
      end do
   end function point_mass

   real(dp) function point_first_moment(mu, sizev, family, points) result(ans)
      real(dp), intent(in) :: mu, sizev
      integer, intent(in) :: family, points(:)
      integer :: j
      ans = 0.0_dp
      do j = 1, size(points)
         ans = ans + real(points(j), dp)*base_pmf(points(j), mu, sizev, family)
      end do
   end function point_first_moment

   elemental real(dp) function base_pmf(y, mu, sizev, family) result(p)
      integer, intent(in) :: y, family
      real(dp), intent(in) :: mu, sizev
      real(dp) :: prob
      if (y < 0 .or. mu <= 0.0_dp) then
         p = 0.0_dp
      else if (family == gaitd_mix_base_poisson) then
         p = dpois_v(y, mu)
      else
         prob = sizev/(sizev + mu); p = dnbinom_v(y, sizev, prob)
      end if
   end function base_pmf

   subroutine copy_int(src, dst)
      integer, intent(in), optional :: src(:)
      integer, allocatable, intent(out) :: dst(:)
      if (present(src)) then
         dst = src
      else
         allocate(dst(0))
      end if
   end subroutine copy_int

   logical function any_duplicate(x) result(ans)
      integer, intent(in) :: x(:)
      integer :: j, k
      ans = .false.
      do j = 1, size(x)
         do k = j + 1, size(x)
            if (x(j) == x(k)) then
               ans = .true.; return
            end if
         end do
      end do
   end function any_duplicate

   elemental real(dp) function clamp_eta(x) result(y)
      real(dp), intent(in) :: x
      y = min(30.0_dp, max(-30.0_dp, x))
   end function clamp_eta

   elemental logical function finite_scalar(x) result(ok)
      real(dp), intent(in) :: x
      ok = abs(x) <= huge(x)
   end function finite_scalar

end module vgam_gaitd_mix_regression
