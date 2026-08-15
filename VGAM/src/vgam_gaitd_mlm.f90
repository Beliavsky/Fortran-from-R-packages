! Copyright (C) 1998-2025 T. W. Yee, University of Auckland.
! Modern Fortran computational translation, 2026.
! SPDX-License-Identifier: GPL-3.0-only
module vgam_gaitd_mlm
   use vgam_kinds, only : dp
   use vgam_distributions, only : dpois_v, dnbinom_v
   use vgam_gaitd, only : gaitd_distribution_t
   use vgam_gaitd_regression, only : gaitd_altered, gaitd_inflated, &
      gaitd_base_poisson, gaitd_base_negbinomial
   use vgam_optim, only : bfgs_minimize, numerical_hessian
   use vgam_linalg, only : invert_matrix
   implicit none
   private

   integer, parameter, public :: gaitd_deflated = 3

   type, public :: gaitd_mlm_regression_result_t
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
      procedure :: predict => predict_gaitd_mlm_regression
   end type gaitd_mlm_regression_result_t

   public :: gaitd_mlm_poisson, gaitd_mlm_negative_binomial
   public :: fit_gaitd_mlm_poisson_regression, fit_gaitd_mlm_nb_regression

contains

   subroutine gaitd_mlm_poisson(lambda, max_support, dist, truncate, altered_points, &
                                 altered_probabilities, inflated_points, inflation_probabilities, &
                                 deflated_points, deflation_probabilities)
      real(dp), intent(in) :: lambda
      integer, intent(in) :: max_support
      type(gaitd_distribution_t), intent(out) :: dist
      integer, intent(in), optional :: truncate(:), altered_points(:), inflated_points(:), deflated_points(:)
      real(dp), intent(in), optional :: altered_probabilities(:), inflation_probabilities(:)
      real(dp), intent(in), optional :: deflation_probabilities(:)
      real(dp), allocatable :: base(:)
      integer :: k
      if (lambda <= 0.0_dp .or. max_support < 0) then
         dist%status = 1
         return
      end if
      allocate(base(max_support + 1))
      do k = 0, max_support
         base(k + 1) = dpois_v(k, lambda)
      end do
      call build_mlm_distribution(base, dist, truncate, altered_points, altered_probabilities, &
                                  inflated_points, inflation_probabilities, deflated_points, &
                                  deflation_probabilities)
   end subroutine gaitd_mlm_poisson

   subroutine gaitd_mlm_negative_binomial(mu, sizev, max_support, dist, truncate, altered_points, &
                                           altered_probabilities, inflated_points, inflation_probabilities, &
                                           deflated_points, deflation_probabilities)
      real(dp), intent(in) :: mu, sizev
      integer, intent(in) :: max_support
      type(gaitd_distribution_t), intent(out) :: dist
      integer, intent(in), optional :: truncate(:), altered_points(:), inflated_points(:), deflated_points(:)
      real(dp), intent(in), optional :: altered_probabilities(:), inflation_probabilities(:)
      real(dp), intent(in), optional :: deflation_probabilities(:)
      real(dp), allocatable :: base(:)
      real(dp) :: prob
      integer :: k
      if (mu <= 0.0_dp .or. sizev <= 0.0_dp .or. max_support < 0) then
         dist%status = 1
         return
      end if
      prob = sizev/(sizev + mu)
      allocate(base(max_support + 1))
      do k = 0, max_support
         base(k + 1) = dnbinom_v(k, sizev, prob)
      end do
      call build_mlm_distribution(base, dist, truncate, altered_points, altered_probabilities, &
                                  inflated_points, inflation_probabilities, deflated_points, &
                                  deflation_probabilities)
   end subroutine gaitd_mlm_negative_binomial

   subroutine build_mlm_distribution(base, dist, truncate, altered_points, altered_probabilities, &
                                      inflated_points, inflation_probabilities, deflated_points, &
                                      deflation_probabilities)
      real(dp), intent(in) :: base(:)
      type(gaitd_distribution_t), intent(out) :: dist
      integer, intent(in), optional :: truncate(:), altered_points(:), inflated_points(:), deflated_points(:)
      real(dp), intent(in), optional :: altered_probabilities(:), inflation_probabilities(:)
      real(dp), intent(in), optional :: deflation_probabilities(:)
      integer, allocatable :: tr(:), ap(:), ip(:), dpnt(:)
      real(dp), allocatable :: pa(:), pi(:), pd(:), work(:)
      real(dp) :: parent_total, sumt, suma_parent, delta, tmp6, ex2, x
      integer :: k, idx

      call unpack_integer(truncate, tr)
      call unpack_integer(altered_points, ap)
      call unpack_integer(inflated_points, ip)
      call unpack_integer(deflated_points, dpnt)
      call unpack_real(altered_probabilities, pa)
      call unpack_real(inflation_probabilities, pi)
      call unpack_real(deflation_probabilities, pd)
      if (size(ap) /= size(pa) .or. size(ip) /= size(pi) .or. size(dpnt) /= size(pd)) then
         dist%status = 2
         return
      end if
      if (any(pa < 0.0_dp) .or. any(pi < 0.0_dp) .or. any(pd < 0.0_dp) .or. &
          any(pa >= 1.0_dp) .or. any(pi >= 1.0_dp) .or. any(pd >= 1.0_dp)) then
         dist%status = 3
         return
      end if
      if (.not. valid_point_sets(size(base), tr, ap, ip, dpnt)) then
         dist%status = 4
         return
      end if
      parent_total = sum(base)
      if (parent_total <= 0.0_dp) then
         dist%status = 5
         return
      end if
      sumt = point_mass_sum(base, tr)
      suma_parent = point_mass_sum(base, ap)
      tmp6 = 1.0_dp - sum(pa) - sum(pi) + sum(pd)
      if (tmp6 <= 0.0_dp .or. parent_total - sumt - suma_parent <= tiny(1.0_dp)) then
         dist%status = 6
         return
      end if
      delta = tmp6/(parent_total - sumt - suma_parent)
      allocate(work(size(base)))
      work = delta*base
      do k = 1, size(tr)
         work(tr(k) + 1) = 0.0_dp
      end do
      do k = 1, size(ap)
         work(ap(k) + 1) = pa(k)
      end do
      do k = 1, size(ip)
         work(ip(k) + 1) = work(ip(k) + 1) + pi(k)
      end do
      do k = 1, size(dpnt)
         idx = dpnt(k) + 1
         work(idx) = work(idx) - pd(k)
      end do
      if (any(work < -1.0e-12_dp) .or. any(work > 1.0_dp + 1.0e-12_dp)) then
         dist%status = 7
         return
      end if
      where (abs(work) < 1.0e-14_dp) work = 0.0_dp
      work = work/sum(work)
      dist%min_support = 0
      dist%pmf = work
      allocate(dist%cdf(size(work)))
      dist%cdf(1) = work(1)
      do k = 2, size(work)
         dist%cdf(k) = min(1.0_dp, dist%cdf(k - 1) + work(k))
      end do
      dist%cdf(size(work)) = 1.0_dp
      dist%mean = 0.0_dp
      ex2 = 0.0_dp
      do k = 1, size(work)
         x = real(k - 1, dp)
         dist%mean = dist%mean + x*work(k)
         ex2 = ex2 + x*x*work(k)
      end do
      dist%variance = max(0.0_dp, ex2 - dist%mean*dist%mean)
      dist%status = 0
   end subroutine build_mlm_distribution

   subroutine fit_gaitd_mlm_poisson_regression(y, x_mean, x_mass, special_points, special_modes, &
                                                result, truncate, weights, max_iter, tol)
      integer, intent(in) :: y(:), special_points(:), special_modes(:)
      real(dp), intent(in) :: x_mean(:, :), x_mass(:, :)
      type(gaitd_mlm_regression_result_t), intent(out) :: result
      integer, intent(in), optional :: truncate(:), max_iter
      real(dp), intent(in), optional :: weights(:), tol
      call fit_gaitd_mlm_common(y, x_mean, x_mass, special_points, special_modes, &
                                gaitd_base_poisson, result, truncate, weights, max_iter, tol)
   end subroutine fit_gaitd_mlm_poisson_regression

   subroutine fit_gaitd_mlm_nb_regression(y, x_mean, x_mass, special_points, special_modes, &
                                           result, truncate, weights, max_iter, tol)
      integer, intent(in) :: y(:), special_points(:), special_modes(:)
      real(dp), intent(in) :: x_mean(:, :), x_mass(:, :)
      type(gaitd_mlm_regression_result_t), intent(out) :: result
      integer, intent(in), optional :: truncate(:), max_iter
      real(dp), intent(in), optional :: weights(:), tol
      call fit_gaitd_mlm_common(y, x_mean, x_mass, special_points, special_modes, &
                                gaitd_base_negbinomial, result, truncate, weights, max_iter, tol)
   end subroutine fit_gaitd_mlm_nb_regression

   subroutine fit_gaitd_mlm_common(y, xm, xz, points, modes, base_family, result, &
                                    truncate, weights, max_iter, tol)
      integer, intent(in) :: y(:), points(:), modes(:), base_family
      real(dp), intent(in) :: xm(:, :), xz(:, :)
      type(gaitd_mlm_regression_result_t), intent(out) :: result
      integer, intent(in), optional :: truncate(:), max_iter
      real(dp), intent(in), optional :: weights(:), tol
      integer, allocatable :: tr(:)
      real(dp), allocatable :: w(:), par(:), hess(:, :), cov(:, :), probs(:)
      real(dp), allocatable :: spec(:, :), pbase(:), fmean(:)
      real(dp) :: fval, tolerance, mean_y, frac, sizev
      integer :: n, pm, pz, q, np, stat, stat2, niter, i, j, pos

      n = size(y)
      pm = size(xm, 2)
      pz = size(xz, 2)
      q = size(points)
      if (n <= 0 .or. pm <= 0 .or. pz <= 0 .or. q <= 0 .or. size(xm, 1) /= n .or. &
          size(xz, 1) /= n .or. size(modes) /= q .or. any(y < 0) .or. any(points < 0) .or. &
          any((modes /= gaitd_altered) .and. (modes /= gaitd_inflated) .and. (modes /= gaitd_deflated))) then
         result%status = 1
         return
      end if
      if (any_duplicate(points)) then
         result%status = 2
         return
      end if
      call unpack_integer(truncate, tr)
      if (any(tr < 0) .or. any_duplicate(tr)) then
         result%status = 3
         return
      end if
      do i = 1, size(tr)
         if (any(points == tr(i))) then
            result%status = 4
            return
         end if
      end do
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
      if (all(abs(xm(:, 1) - 1.0_dp) < 100.0_dp*epsilon(1.0_dp))) par(1) = log(max(mean_y, 0.1_dp))
      pos = pm
      do j = 1, q
         if (modes(j) == gaitd_deflated) then
            frac = 1.0e-3_dp
         else
            frac = sum(w*merge(1.0_dp, 0.0_dp, y == points(j)))/max(sum(w), tiny(1.0_dp))
            frac = max(1.0e-3_dp, min(0.15_dp, 0.5_dp*frac))
         end if
         if (all(abs(xz(:, 1) - 1.0_dp) < 100.0_dp*epsilon(1.0_dp))) par(pos + 1) = log(frac)
         pos = pos + pz
      end do
      if (base_family == gaitd_base_negbinomial) par(np) = log(2.0_dp)
      niter = 450
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
         fmean(i) = mlm_mean(xm(i, :), result%mean_coefficients, sizev, base_family, probs)
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
            if (.not. mlm_valid(muv, sz, base_family, mix)) then
               nll = huge(1.0_dp)/100.0_dp
               return
            end if
            py = mlm_probability(y(ii), muv, sz, base_family, mix)
            if (py <= tiny(1.0_dp) .or. .not. finite_scalar(py)) then
               nll = huge(1.0_dp)/100.0_dp
               return
            end if
            nll = nll - w(ii)*log(py)
         end do
      end function objective

      logical function mlm_valid(muv, sz, fam, mix) result(ok)
         real(dp), intent(in) :: muv, sz, mix(:)
         integer, intent(in) :: fam
         real(dp) :: delta, denom, tmp6, pval
         integer :: jj
         call mlm_scale(muv, sz, fam, mix, delta, denom, tmp6)
         ok = denom > tiny(1.0_dp) .and. tmp6 > tiny(1.0_dp) .and. finite_scalar(delta)
         if (.not. ok) return
         do jj = 1, q
            if (modes(jj) == gaitd_deflated) then
               pval = delta*base_pmf(points(jj), muv, sz, fam) - mix(jj + 1)
               if (pval < -1.0e-12_dp) then
                  ok = .false.
                  return
               end if
            else if (modes(jj) == gaitd_inflated) then
               pval = delta*base_pmf(points(jj), muv, sz, fam) + mix(jj + 1)
               if (pval > 1.0_dp + 1.0e-12_dp) then
                  ok = .false.
                  return
               end if
            end if
         end do
      end function mlm_valid

      real(dp) function mlm_probability(yy, muv, sz, fam, mix) result(py)
         integer, intent(in) :: yy, fam
         real(dp), intent(in) :: muv, sz, mix(:)
         real(dp) :: delta, denom, tmp6
         integer :: jj
         if (any(tr == yy)) then
            py = 0.0_dp
            return
         end if
         call mlm_scale(muv, sz, fam, mix, delta, denom, tmp6)
         py = delta*base_pmf(yy, muv, sz, fam)
         do jj = 1, q
            if (yy /= points(jj)) cycle
            select case (modes(jj))
            case (gaitd_altered)
               py = mix(jj + 1)
            case (gaitd_inflated)
               py = py + mix(jj + 1)
            case (gaitd_deflated)
               py = py - mix(jj + 1)
            end select
         end do
      end function mlm_probability

      real(dp) function mlm_mean(xrow, beta, sz, fam, mix) result(ans)
         real(dp), intent(in) :: xrow(:), beta(:), sz, mix(:)
         integer, intent(in) :: fam
         real(dp) :: muv, delta, denom, tmp6, removed
         integer :: jj
         muv = exp(clamp_eta(dot_product(xrow, beta)))
         call mlm_scale(muv, sz, fam, mix, delta, denom, tmp6)
         removed = 0.0_dp
         do jj = 1, size(tr)
            removed = removed + real(tr(jj), dp)*base_pmf(tr(jj), muv, sz, fam)
         end do
         do jj = 1, q
            if (modes(jj) == gaitd_altered) then
               removed = removed + real(points(jj), dp)*base_pmf(points(jj), muv, sz, fam)
            end if
         end do
         ans = delta*(muv - removed)
         do jj = 1, q
            select case (modes(jj))
            case (gaitd_altered, gaitd_inflated)
               ans = ans + real(points(jj), dp)*mix(jj + 1)
            case (gaitd_deflated)
               ans = ans - real(points(jj), dp)*mix(jj + 1)
            end select
         end do
      end function mlm_mean

      subroutine mlm_scale(muv, sz, fam, mix, delta, denom, tmp6)
         real(dp), intent(in) :: muv, sz, mix(:)
         integer, intent(in) :: fam
         real(dp), intent(out) :: delta, denom, tmp6
         real(dp) :: sumt, suma, suma_prob, sumi_prob, sumd_prob
         integer :: jj
         sumt = 0.0_dp
         suma = 0.0_dp
         suma_prob = 0.0_dp
         sumi_prob = 0.0_dp
         sumd_prob = 0.0_dp
         do jj = 1, size(tr)
            sumt = sumt + base_pmf(tr(jj), muv, sz, fam)
         end do
         do jj = 1, q
            select case (modes(jj))
            case (gaitd_altered)
               suma = suma + base_pmf(points(jj), muv, sz, fam)
               suma_prob = suma_prob + mix(jj + 1)
            case (gaitd_inflated)
               sumi_prob = sumi_prob + mix(jj + 1)
            case (gaitd_deflated)
               sumd_prob = sumd_prob + mix(jj + 1)
            end select
         end do
         denom = 1.0_dp - sumt - suma
         tmp6 = 1.0_dp - suma_prob - sumi_prob + sumd_prob
         delta = tmp6/max(denom, tiny(1.0_dp))
      end subroutine mlm_scale

   end subroutine fit_gaitd_mlm_common

   subroutine predict_gaitd_mlm_regression(self, x_mean, x_mass, fitted_mean, &
                                            baseline_probability, special_probabilities)
      class(gaitd_mlm_regression_result_t), intent(in) :: self
      real(dp), intent(in) :: x_mean(:, :), x_mass(:, :)
      real(dp), allocatable, intent(out) :: fitted_mean(:)
      real(dp), allocatable, intent(out), optional :: baseline_probability(:)
      real(dp), allocatable, intent(out), optional :: special_probabilities(:, :)
      real(dp), allocatable :: mix(:)
      real(dp) :: muv, delta, denom, tmp6, removed, suma, sumt
      real(dp) :: suma_prob, sumi_prob, sumd_prob
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
         sumt = 0.0_dp
         suma = 0.0_dp
         removed = 0.0_dp
         suma_prob = 0.0_dp
         sumi_prob = 0.0_dp
         sumd_prob = 0.0_dp
         do j = 1, size(self%truncate)
            sumt = sumt + base_pmf(self%truncate(j), muv, self%size, self%base_family)
            removed = removed + real(self%truncate(j), dp)* &
                      base_pmf(self%truncate(j), muv, self%size, self%base_family)
         end do
         do j = 1, q
            select case (self%special_modes(j))
            case (gaitd_altered)
               suma = suma + base_pmf(self%special_points(j), muv, self%size, self%base_family)
               removed = removed + real(self%special_points(j), dp)* &
                         base_pmf(self%special_points(j), muv, self%size, self%base_family)
               suma_prob = suma_prob + mix(j + 1)
            case (gaitd_inflated)
               sumi_prob = sumi_prob + mix(j + 1)
            case (gaitd_deflated)
               sumd_prob = sumd_prob + mix(j + 1)
            end select
         end do
         denom = 1.0_dp - sumt - suma
         tmp6 = 1.0_dp - suma_prob - sumi_prob + sumd_prob
         delta = tmp6/max(denom, tiny(1.0_dp))
         fitted_mean(i) = delta*(muv - removed)
         do j = 1, q
            if (self%special_modes(j) == gaitd_deflated) then
               fitted_mean(i) = fitted_mean(i) - real(self%special_points(j), dp)*mix(j + 1)
            else
               fitted_mean(i) = fitted_mean(i) + real(self%special_points(j), dp)*mix(j + 1)
            end if
         end do
         if (present(baseline_probability)) baseline_probability(i) = mix(1)
         if (present(special_probabilities)) special_probabilities(i, :) = mix(2:)
      end do
   end subroutine predict_gaitd_mlm_regression

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

   pure real(dp) function point_mass_sum(base, points) result(ans)
      real(dp), intent(in) :: base(:)
      integer, intent(in) :: points(:)
      integer :: i
      ans = 0.0_dp
      do i = 1, size(points)
         ans = ans + base(points(i) + 1)
      end do
   end function point_mass_sum

   logical function valid_point_sets(nbase, tr, ap, ip, dpnt) result(ok)
      integer, intent(in) :: nbase, tr(:), ap(:), ip(:), dpnt(:)
      integer, allocatable :: allp(:)
      integer :: n, pos
      n = size(tr) + size(ap) + size(ip) + size(dpnt)
      allocate(allp(n))
      pos = 0
      call append_points(allp, pos, tr)
      call append_points(allp, pos, ap)
      call append_points(allp, pos, ip)
      call append_points(allp, pos, dpnt)
      ok = all(allp >= 0) .and. all(allp < nbase) .and. .not. any_duplicate(allp)
   end function valid_point_sets

   subroutine append_points(dest, pos, src)
      integer, intent(inout) :: dest(:), pos
      integer, intent(in) :: src(:)
      if (size(src) > 0) then
         dest(pos + 1:pos + size(src)) = src
         pos = pos + size(src)
      end if
   end subroutine append_points

   subroutine unpack_integer(x, y)
      integer, intent(in), optional :: x(:)
      integer, allocatable, intent(out) :: y(:)
      if (present(x)) then
         y = x
      else
         allocate(y(0))
      end if
   end subroutine unpack_integer

   subroutine unpack_real(x, y)
      real(dp), intent(in), optional :: x(:)
      real(dp), allocatable, intent(out) :: y(:)
      if (present(x)) then
         y = x
      else
         allocate(y(0))
      end if
   end subroutine unpack_real

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

   elemental real(dp) function clamp_eta(x) result(y)
      real(dp), intent(in) :: x
      y = min(35.0_dp, max(-35.0_dp, x))
   end function clamp_eta

   elemental logical function finite_scalar(x) result(ok)
      real(dp), intent(in) :: x
      ok = abs(x) <= huge(x)
   end function finite_scalar

end module vgam_gaitd_mlm
