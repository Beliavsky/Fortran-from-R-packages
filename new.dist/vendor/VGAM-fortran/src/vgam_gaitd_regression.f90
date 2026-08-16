! Copyright (C) 1998-2025 T. W. Yee, University of Auckland.
! Modern Fortran computational translation, 2026.
! SPDX-License-Identifier: GPL-3.0-only
module vgam_gaitd_regression
   use vgam_kinds, only : dp
   use vgam_distributions, only : dpois_v, dnbinom_v
   use vgam_optim, only : bfgs_minimize, numerical_hessian
   use vgam_linalg, only : invert_matrix
   implicit none
   private

   integer, parameter, public :: gaitd_altered = 1
   integer, parameter, public :: gaitd_inflated = 2
   integer, parameter, public :: gaitd_base_poisson = 1
   integer, parameter, public :: gaitd_base_negbinomial = 2

   type, public :: gaitd_regression_result_t
      real(dp), allocatable :: mean_coefficients(:)
      real(dp), allocatable :: mass_coefficients(:, :)
      real(dp), allocatable :: covariance(:, :)
      real(dp), allocatable :: fitted_mean(:)
      real(dp), allocatable :: baseline_probability(:)
      real(dp), allocatable :: special_probabilities(:, :)
      integer, allocatable :: special_points(:)
      integer, allocatable :: special_modes(:)
      integer, allocatable :: truncate(:)
      integer :: base_family = gaitd_base_poisson
      real(dp) :: size = huge(1.0_dp)
      real(dp) :: loglik = -huge(1.0_dp)
      real(dp) :: aic = huge(1.0_dp)
      integer :: status = 0
      logical :: converged = .false.
   contains
      procedure :: predict => predict_gaitd_regression
   end type gaitd_regression_result_t

   public :: fit_gaitd_poisson_regression, fit_gaitd_nb_regression

contains

   subroutine fit_gaitd_poisson_regression(y, x_mean, x_mass, special_points, &
                                            special_modes, result, truncate, weights, &
                                            max_iter, tol)
      integer, intent(in) :: y(:), special_points(:), special_modes(:)
      real(dp), intent(in) :: x_mean(:, :), x_mass(:, :)
      type(gaitd_regression_result_t), intent(out) :: result
      integer, intent(in), optional :: truncate(:), max_iter
      real(dp), intent(in), optional :: weights(:), tol
      call fit_gaitd_common(y, x_mean, x_mass, special_points, special_modes, &
                            gaitd_base_poisson, result, truncate, weights, max_iter, tol)
   end subroutine fit_gaitd_poisson_regression

   subroutine fit_gaitd_nb_regression(y, x_mean, x_mass, special_points, &
                                      special_modes, result, truncate, weights, &
                                      max_iter, tol)
      integer, intent(in) :: y(:), special_points(:), special_modes(:)
      real(dp), intent(in) :: x_mean(:, :), x_mass(:, :)
      type(gaitd_regression_result_t), intent(out) :: result
      integer, intent(in), optional :: truncate(:), max_iter
      real(dp), intent(in), optional :: weights(:), tol
      call fit_gaitd_common(y, x_mean, x_mass, special_points, special_modes, &
                            gaitd_base_negbinomial, result, truncate, weights, max_iter, tol)
   end subroutine fit_gaitd_nb_regression

   subroutine fit_gaitd_common(y, xm, xz, points, modes, base_family, result, &
                               truncate, weights, max_iter, tol)
      integer, intent(in) :: y(:), points(:), modes(:), base_family
      real(dp), intent(in) :: xm(:, :), xz(:, :)
      type(gaitd_regression_result_t), intent(out) :: result
      integer, intent(in), optional :: truncate(:), max_iter
      real(dp), intent(in), optional :: weights(:), tol
      real(dp), allocatable :: w(:), par(:), hess(:, :), cov(:, :), probs(:)
      real(dp), allocatable :: spec(:, :), pbase(:), fmean(:)
      integer, allocatable :: tr(:)
      real(dp) :: fval, tolerance, mean_y, frac, sizev
      integer :: n, pm, pz, q, np, stat, stat2, niter, i, j, pos

      n = size(y)
      pm = size(xm, 2)
      pz = size(xz, 2)
      q = size(points)
      if (n <= 0 .or. pm <= 0 .or. pz <= 0 .or. q <= 0 .or. &
          size(xm, 1) /= n .or. size(xz, 1) /= n .or. size(modes) /= q .or. &
          any(y < 0) .or. any(points < 0) .or. any((modes /= gaitd_altered) .and. &
          (modes /= gaitd_inflated))) then
         result%status = 1
         return
      end if
      if (any_duplicate(points)) then
         result%status = 2
         return
      end if
      if (present(truncate)) then
         if (any(truncate < 0) .or. any_duplicate(truncate)) then
            result%status = 3
            return
         end if
         do i = 1, size(truncate)
            if (any(points == truncate(i))) then
               result%status = 4
               return
            end if
         end do
         tr = truncate
      else
         allocate(tr(0))
      end if
      allocate(w(n))
      w = 1.0_dp
      if (present(weights)) then
         if (size(weights) /= n .or. any(weights < 0.0_dp)) then
            result%status = 5
            return
         end if
         w = weights
      end if

      np = pm + pz*q + merge(1, 0, base_family == gaitd_base_negbinomial)
      allocate(par(np))
      par = 0.0_dp
      mean_y = sum(w*real(y, dp))/max(sum(w), tiny(1.0_dp))
      if (all(abs(xm(:, 1) - 1.0_dp) < 100.0_dp*epsilon(1.0_dp))) then
         par(1) = log(max(mean_y, 0.1_dp))
      end if
      pos = pm
      do j = 1, q
         frac = sum(w*merge(1.0_dp, 0.0_dp, y == points(j)))/max(sum(w), tiny(1.0_dp))
         if (all(abs(xz(:, 1) - 1.0_dp) < 100.0_dp*epsilon(1.0_dp))) then
            par(pos + 1) = log(max(frac, 1.0e-3_dp)/max(1.0_dp - frac, 1.0e-3_dp))
         end if
         pos = pos + pz
      end do
      if (base_family == gaitd_base_negbinomial) par(np) = log(2.0_dp)

      niter = 400
      if (present(max_iter)) niter = max_iter
      tolerance = 1.0e-7_dp
      if (present(tol)) tolerance = tol
      call bfgs_minimize(objective, par, fval, stat, max_iter=niter, tol=tolerance)
      result%status = stat
      result%converged = stat == 0
      result%base_family = base_family
      result%mean_coefficients = par(1:pm)
      allocate(result%mass_coefficients(pz, q))
      result%mass_coefficients = reshape(par(pm + 1:pm + pz*q), [pz, q])
      sizev = huge(1.0_dp)
      if (base_family == gaitd_base_negbinomial) sizev = exp(min(par(np), 50.0_dp))
      result%size = sizev
      result%special_points = points
      result%special_modes = modes
      result%truncate = tr
      result%loglik = -fval
      result%aic = 2.0_dp*fval + 2.0_dp*real(np, dp)

      allocate(spec(n, q), pbase(n), fmean(n), probs(q + 1))
      do i = 1, n
         call mixing_probabilities(xz(i, :), result%mass_coefficients, probs)
         pbase(i) = probs(1)
         spec(i, :) = probs(2:)
         fmean(i) = observation_mean(xm(i, :), result%mean_coefficients, &
                                     sizev, base_family, probs)
      end do
      result%baseline_probability = pbase
      result%special_probabilities = spec
      result%fitted_mean = fmean

      allocate(hess(np, np))
      call numerical_hessian(objective, par, hess)
      call invert_matrix(hess, cov, stat2)
      if (stat2 == 0) then
         result%covariance = cov
      else
         allocate(result%covariance(0, 0))
      end if

   contains

      real(dp) function objective(theta) result(nll)
         real(dp), intent(in) :: theta(:)
         real(dp) :: beta(pm), gamma(pz, q), mix(q + 1), muv, sz, py
         integer :: ii
         beta = theta(1:pm)
         gamma = reshape(theta(pm + 1:pm + pz*q), [pz, q])
         sz = huge(1.0_dp)
         if (base_family == gaitd_base_negbinomial) sz = exp(min(theta(np), 50.0_dp))
         nll = 0.0_dp
         do ii = 1, n
            muv = exp(clamp_eta(dot_product(xm(ii, :), beta)))
            call mixing_probabilities(xz(ii, :), gamma, mix)
            py = gaitd_probability(y(ii), muv, sz, base_family, mix)
            if (py <= 0.0_dp .or. .not. finite_scalar(py)) then
               nll = huge(1.0_dp)/100.0_dp
               return
            end if
            nll = nll - w(ii)*log(py)
         end do
      end function objective

      real(dp) function gaitd_probability(yy, muv, sz, fam, mix) result(py)
         integer, intent(in) :: yy, fam
         real(dp), intent(in) :: muv, sz, mix(:)
         real(dp) :: norm, basep
         integer :: jj
         norm = base_normalizer(muv, sz, fam)
         if (norm <= tiny(1.0_dp)) then
            py = 0.0_dp
            return
         end if
         py = 0.0_dp
         if (.not. any(tr == yy) .and. .not. altered_point(yy)) then
            basep = base_pmf(yy, muv, sz, fam)/norm
            py = mix(1)*basep
         end if
         do jj = 1, q
            if (yy == points(jj)) py = py + mix(jj + 1)
         end do
      end function gaitd_probability

      real(dp) function observation_mean(xrow, beta, sz, fam, mix) result(ans)
         real(dp), intent(in) :: xrow(:), beta(:), sz, mix(:)
         integer, intent(in) :: fam
         real(dp) :: muv, norm, removed_mean
         integer :: jj
         muv = exp(clamp_eta(dot_product(xrow, beta)))
         norm = base_normalizer(muv, sz, fam)
         removed_mean = 0.0_dp
         do jj = 1, size(tr)
            removed_mean = removed_mean + real(tr(jj), dp)*base_pmf(tr(jj), muv, sz, fam)
         end do
         do jj = 1, q
            if (modes(jj) == gaitd_altered) then
               removed_mean = removed_mean + real(points(jj), dp)*base_pmf(points(jj), muv, sz, fam)
            end if
         end do
         ans = mix(1)*(muv - removed_mean)/max(norm, tiny(1.0_dp))
         do jj = 1, q
            ans = ans + mix(jj + 1)*real(points(jj), dp)
         end do
      end function observation_mean

      real(dp) function base_normalizer(muv, sz, fam) result(norm)
         real(dp), intent(in) :: muv, sz
         integer, intent(in) :: fam
         integer :: jj
         norm = 1.0_dp
         do jj = 1, size(tr)
            norm = norm - base_pmf(tr(jj), muv, sz, fam)
         end do
         do jj = 1, q
            if (modes(jj) == gaitd_altered) norm = norm - base_pmf(points(jj), muv, sz, fam)
         end do
      end function base_normalizer

      logical function altered_point(yy) result(ans)
         integer, intent(in) :: yy
         integer :: jj
         ans = .false.
         do jj = 1, q
            if (modes(jj) == gaitd_altered .and. points(jj) == yy) then
               ans = .true.
               return
            end if
         end do
      end function altered_point

   end subroutine fit_gaitd_common

   subroutine predict_gaitd_regression(self, x_mean, x_mass, fitted_mean, &
                                       baseline_probability, special_probabilities)
      class(gaitd_regression_result_t), intent(in) :: self
      real(dp), intent(in) :: x_mean(:, :), x_mass(:, :)
      real(dp), allocatable, intent(out) :: fitted_mean(:)
      real(dp), allocatable, intent(out), optional :: baseline_probability(:)
      real(dp), allocatable, intent(out), optional :: special_probabilities(:, :)
      real(dp), allocatable :: mix(:)
      real(dp) :: muv, norm, removed_mean
      integer :: n, q, i, j

      n = size(x_mean, 1)
      q = size(self%special_points)
      if (size(x_mass, 1) /= n .or. size(x_mean, 2) /= size(self%mean_coefficients) .or. &
          size(x_mass, 2) /= size(self%mass_coefficients, 1)) then
         allocate(fitted_mean(0))
         if (present(baseline_probability)) allocate(baseline_probability(0))
         if (present(special_probabilities)) allocate(special_probabilities(0, 0))
         return
      end if
      allocate(fitted_mean(n), mix(q + 1))
      if (present(baseline_probability)) allocate(baseline_probability(n))
      if (present(special_probabilities)) allocate(special_probabilities(n, q))
      do i = 1, n
         call mixing_probabilities(x_mass(i, :), self%mass_coefficients, mix)
         muv = exp(clamp_eta(dot_product(x_mean(i, :), self%mean_coefficients)))
         norm = 1.0_dp
         removed_mean = 0.0_dp
         do j = 1, size(self%truncate)
            norm = norm - base_pmf(self%truncate(j), muv, self%size, self%base_family)
            removed_mean = removed_mean + real(self%truncate(j), dp)* &
               base_pmf(self%truncate(j), muv, self%size, self%base_family)
         end do
         do j = 1, q
            if (self%special_modes(j) == gaitd_altered) then
               norm = norm - base_pmf(self%special_points(j), muv, self%size, self%base_family)
               removed_mean = removed_mean + real(self%special_points(j), dp)* &
                  base_pmf(self%special_points(j), muv, self%size, self%base_family)
            end if
         end do
         fitted_mean(i) = mix(1)*(muv - removed_mean)/max(norm, tiny(1.0_dp))
         do j = 1, q
            fitted_mean(i) = fitted_mean(i) + mix(j + 1)*real(self%special_points(j), dp)
         end do
         if (present(baseline_probability)) baseline_probability(i) = mix(1)
         if (present(special_probabilities)) special_probabilities(i, :) = mix(2:)
      end do
   end subroutine predict_gaitd_regression

   subroutine mixing_probabilities(xrow, gamma, probs)
      real(dp), intent(in) :: xrow(:), gamma(:, :)
      real(dp), intent(out) :: probs(:)
      real(dp), allocatable :: eta(:)
      real(dp) :: vmax, den
      integer :: q
      q = size(gamma, 2)
      allocate(eta(q))
      eta = matmul(transpose(gamma), xrow)
      vmax = max(0.0_dp, maxval(eta))
      den = exp(-vmax) + sum(exp(eta - vmax))
      probs(1) = exp(-vmax)/den
      probs(2:q + 1) = exp(eta - vmax)/den
   end subroutine mixing_probabilities

   elemental real(dp) function base_pmf(y, mu, sizev, family) result(p)
      integer, intent(in) :: y, family
      real(dp), intent(in) :: mu, sizev
      real(dp) :: prob
      if (y < 0 .or. mu <= 0.0_dp) then
         p = 0.0_dp
      else if (family == gaitd_base_poisson) then
         p = dpois_v(y, mu)
      else
         prob = sizev/(sizev + mu)
         p = dnbinom_v(y, sizev, prob)
      end if
   end function base_pmf

   elemental real(dp) function clamp_eta(x) result(y)
      real(dp), intent(in) :: x
      y = min(35.0_dp, max(-35.0_dp, x))
   end function clamp_eta

   logical function any_duplicate(x) result(ans)
      integer, intent(in) :: x(:)
      integer :: i, j
      ans = .false.
      do i = 1, size(x)
         do j = i + 1, size(x)
            if (x(i) == x(j)) then
               ans = .true.
               return
            end if
         end do
      end do
   end function any_duplicate

   elemental logical function finite_scalar(x) result(ok)
      real(dp), intent(in) :: x
      ok = abs(x) <= huge(x)
   end function finite_scalar

end module vgam_gaitd_regression
