! SPDX-License-Identifier: GPL-2.0-or-later
! Copyright (C) 2008-2021 David Ardia
! Modern Fortran translation of computational routines from bayesGARCH.
!
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 2 of the License, or
! (at your option) any later version.
module bayesgarch_sampler
   use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
   use bayesgarch_garch, only : filter_alpha, filter_w, garch11_filter, quasi_difference
   use bayesgarch_kinds, only : dp
   use bayesgarch_math, only : digamma, draw_mvn2, inverse_2x2, log_mvn2_density, &
      log_normal_density
   use bayesgarch_rng, only : random_exponential, random_gamma, random_normal, random_uniform
   implicit none
   private

   type, public :: bayesgarch_priors
      real(dp) :: mu_alpha(2) = [0.0_dp, 0.0_dp]
      real(dp) :: sigma_alpha(2, 2) = reshape([1000.0_dp, 0.0_dp, 0.0_dp, 1000.0_dp], [2, 2])
      real(dp) :: mu_beta = 0.0_dp
      real(dp) :: sigma_beta = 1000.0_dp
      real(dp) :: lambda = 0.01_dp
      real(dp) :: delta = 2.0_dp
   end type bayesgarch_priors

   type, public :: bayesgarch_control
      integer :: n_chains = 1
      integer :: n_iter = 10000
      integer :: refresh = 0
      logical :: enforce_stationarity = .false.
      real(dp) :: start(4) = [0.01_dp, 0.1_dp, 0.7_dp, 100.0_dp]
   end type bayesgarch_control

   type, public :: bayesgarch_result
      real(dp), allocatable :: draws(:, :, :)
      integer, allocatable :: alpha_accept(:)
      integer, allocatable :: beta_accept(:)
      integer, allocatable :: nu_updates(:)
      integer, allocatable :: constraint_reject(:)
      integer, allocatable :: nu_root_fail(:)
   end type bayesgarch_result

   abstract interface
      logical function parameter_constraint(psi)
         import dp
         real(dp), intent(in) :: psi(4)
      end function parameter_constraint
   end interface

   public :: run_bayesgarch
   public :: draw_latent_scales
   public :: augmented_log_posterior
   public :: nu_root_function
   public :: nu_log_acceptance
   public :: solve_nu_proposal_rate
   public :: regression_posterior_2
   public :: regression_posterior_1
   public :: parameter_constraint

contains

   subroutine run_bayesgarch(y, result, priors, control, start_values, constraint)
      real(dp), intent(in) :: y(:)
      type(bayesgarch_result), intent(out) :: result
      type(bayesgarch_priors), intent(in), optional :: priors
      type(bayesgarch_control), intent(in), optional :: control
      real(dp), intent(in), optional :: start_values(:, :)
      procedure(parameter_constraint), optional :: constraint
      type(bayesgarch_priors) :: prior
      type(bayesgarch_control) :: ctrl
      real(dp) :: prior_precision(2, 2)
      real(dp) :: determinant
      real(dp) :: current(4)
      real(dp) :: proposed(4)
      real(dp), allocatable :: w(:)
      integer :: chain
      integer :: iteration
      logical :: accepted_alpha
      logical :: accepted_beta
      logical :: changed_nu
      logical :: root_failed
      logical :: valid

      if (size(y) < 2) error stop "run_bayesgarch: y must contain at least two observations"
      if (any(ieee_is_nan(y))) error stop "run_bayesgarch: y contains NaN"

      prior = bayesgarch_priors()
      if (present(priors)) prior = priors
      ctrl = bayesgarch_control()
      if (present(control)) ctrl = control
      call validate_inputs(prior, ctrl, start_values)
      call inverse_2x2(prior%sigma_alpha, prior_precision, determinant, valid)
      if (.not. valid) error stop "run_bayesgarch: sigma_alpha is not positive definite"

      allocate(result%draws(ctrl%n_iter, 4, ctrl%n_chains))
      allocate(result%alpha_accept(ctrl%n_chains), result%beta_accept(ctrl%n_chains))
      allocate(result%nu_updates(ctrl%n_chains), result%constraint_reject(ctrl%n_chains))
      allocate(result%nu_root_fail(ctrl%n_chains))
      result%alpha_accept = 0
      result%beta_accept = 0
      result%nu_updates = 0
      result%constraint_reject = 0
      result%nu_root_fail = 0
      allocate(w(size(y)))

      do chain = 1, ctrl%n_chains
         if (present(start_values)) then
            current = start_values(:, chain)
         else
            current = ctrl%start
         end if
         result%draws(1, :, chain) = current

         do iteration = 2, ctrl%n_iter
            call draw_latent_scales(y, current(1:2), current(3), current(4), w)
            proposed = current
            call draw_alpha_mh(y, current(1:2), current(3), current(4), w, prior, prior_precision, &
               proposed(1:2), accepted_alpha)
            call draw_beta_mh(y, proposed(1:2), current(3), current(4), w, prior, prior_precision, &
               proposed(3), accepted_beta)
            call draw_nu_conditional(w, current(4), prior%lambda, prior%delta, proposed(4), &
               changed_nu, root_failed)

            valid = all(proposed > 0.0_dp) .and. proposed(4) > prior%delta
            if (ctrl%enforce_stationarity) valid = valid .and. proposed(2) + proposed(3) < 1.0_dp
            if (present(constraint)) then
               if (.not. constraint(proposed)) valid = .false.
            end if

            if (valid) then
               current = proposed
               if (accepted_alpha) result%alpha_accept(chain) = result%alpha_accept(chain) + 1
               if (accepted_beta) result%beta_accept(chain) = result%beta_accept(chain) + 1
               if (changed_nu) result%nu_updates(chain) = result%nu_updates(chain) + 1
            else
               result%constraint_reject(chain) = result%constraint_reject(chain) + 1
            end if
            if (root_failed) result%nu_root_fail(chain) = result%nu_root_fail(chain) + 1
            result%draws(iteration, :, chain) = current

            if (ctrl%refresh > 0) then
               if (mod(iteration, ctrl%refresh) == 0) then
                  write(*, '(a,i0,a,i0,a,4(1x,es13.5))') "chain ", chain, " iteration ", iteration, &
                     " parameters", current
               end if
            end if
         end do
      end do
   end subroutine run_bayesgarch

   subroutine validate_inputs(prior, control, start_values)
      type(bayesgarch_priors), intent(in) :: prior
      type(bayesgarch_control), intent(in) :: control
      real(dp), intent(in), optional :: start_values(:, :)
      real(dp) :: inverse(2, 2)
      real(dp) :: determinant
      logical :: ok

      if (control%n_chains < 1) error stop "run_bayesgarch: n_chains must be positive"
      if (control%n_iter < 2) error stop "run_bayesgarch: n_iter must be at least two"
      if (prior%sigma_beta <= 0.0_dp) error stop "run_bayesgarch: sigma_beta must be positive"
      if (prior%lambda <= 0.0_dp) error stop "run_bayesgarch: lambda must be positive"
      if (prior%delta < 2.0_dp) error stop "run_bayesgarch: delta must be at least two"
      if (abs(prior%sigma_alpha(1, 2) - prior%sigma_alpha(2, 1)) > 100.0_dp * epsilon(1.0_dp)) &
         error stop "run_bayesgarch: sigma_alpha must be symmetric"
      call inverse_2x2(prior%sigma_alpha, inverse, determinant, ok)
      if (.not. ok .or. prior%sigma_alpha(1, 1) <= 0.0_dp) &
         error stop "run_bayesgarch: sigma_alpha must be positive definite"
      if (any(control%start <= 0.0_dp)) error stop "run_bayesgarch: start values must be positive"
      if (control%start(4) <= prior%delta) error stop "run_bayesgarch: starting nu must exceed delta"
      if (present(start_values)) then
         if (size(start_values, 1) /= 4 .or. size(start_values, 2) /= control%n_chains) &
            error stop "run_bayesgarch: start_values must have shape (4,n_chains)"
         if (any(start_values <= 0.0_dp)) error stop "run_bayesgarch: start_values must be positive"
         if (any(start_values(4, :) <= prior%delta)) &
            error stop "run_bayesgarch: all starting nu values must exceed delta"
      end if
   end subroutine validate_inputs

   subroutine draw_latent_scales(y, alpha, beta, nu, w)
      real(dp), intent(in) :: y(:)
      real(dp), intent(in) :: alpha(2)
      real(dp), intent(in) :: beta
      real(dp), intent(in) :: nu
      real(dp), intent(out) :: w(size(y))
      real(dp) :: h(size(y) + 1)
      real(dp) :: shape
      real(dp) :: rate
      integer :: i

      call garch11_filter(y, [alpha, beta], h)
      shape = 0.5_dp * (nu + 1.0_dp)
      do i = 1, size(y)
         rate = 0.5_dp * (y(i)**2 / h(i) + nu - 2.0_dp)
         w(i) = 1.0_dp / random_gamma(shape, rate)
      end do
   end subroutine draw_latent_scales

   function augmented_log_posterior(y, augmented_variance, alpha, beta, nu, prior, prior_precision) result(value)
      real(dp), intent(in) :: y(:)
      real(dp), intent(in) :: augmented_variance(:)
      real(dp), intent(in) :: alpha(2)
      real(dp), intent(in) :: beta
      real(dp), intent(in) :: nu
      type(bayesgarch_priors), intent(in) :: prior
      real(dp), intent(in), optional :: prior_precision(2, 2)
      real(dp) :: value
      real(dp) :: precision(2, 2)
      real(dp) :: inverse(2, 2)
      real(dp) :: determinant
      real(dp) :: difference(2)
      logical :: ok

      if (size(augmented_variance) /= size(y) .or. any(augmented_variance <= 0.0_dp)) then
         value = -huge(1.0_dp)
         return
      end if
      if (present(prior_precision)) then
         precision = prior_precision
      else
         call inverse_2x2(prior%sigma_alpha, inverse, determinant, ok)
         if (.not. ok) then
            value = -huge(1.0_dp)
            return
         end if
         precision = inverse
      end if
      difference = alpha - prior%mu_alpha
      value = -0.5_dp * sum(log(augmented_variance) + y**2 / augmented_variance)
      value = value - 0.5_dp * dot_product(difference, matmul(precision, difference))
      value = value - 0.5_dp * (beta - prior%mu_beta)**2 / prior%sigma_beta
      value = value - prior%lambda * nu
   end function augmented_log_posterior

   subroutine draw_alpha_mh(y, alpha_old, beta, nu, w, prior, prior_precision, alpha_new, accepted)
      real(dp), intent(in) :: y(:)
      real(dp), intent(in) :: alpha_old(2)
      real(dp), intent(in) :: beta
      real(dp), intent(in) :: nu
      real(dp), intent(in) :: w(:)
      type(bayesgarch_priors), intent(in) :: prior
      real(dp), intent(in) :: prior_precision(2, 2)
      real(dp), intent(out) :: alpha_new(2)
      logical, intent(out) :: accepted
      real(dp) :: candidate(2)
      real(dp) :: covariance(2, 2)
      real(dp) :: covariance_new(2, 2)
      real(dp) :: h(size(y) + 1)
      real(dp) :: h_new(size(y) + 1)
      real(dp) :: mean(2)
      real(dp) :: mean_new(2)
      real(dp) :: post_old
      real(dp) :: post_new
      real(dp) :: proposal_new
      real(dp) :: proposal_old
      real(dp) :: tau(size(y))
      real(dp) :: tau_new(size(y))
      real(dp) :: x(size(y), 2)
      logical :: ok
      real(dp) :: ratio

      alpha_new = alpha_old
      accepted = .false.
      call garch11_filter(y, [alpha_old, beta], h)
      post_old = augmented_log_posterior(y, w * h(:size(y)), alpha_old, beta, nu, prior, prior_precision)
      tau = 2.0_dp * h(:size(y))**2
      call filter_alpha(y, beta, x)
      call regression_posterior_2(y**2 / w, x, tau, prior%mu_alpha, prior_precision, covariance, mean, ok)
      if (.not. ok) return
      candidate = draw_mvn2(mean, covariance)
      if (any(candidate <= 0.0_dp)) return

      call garch11_filter(y, [candidate, beta], h_new)
      tau_new = 2.0_dp * h_new(:size(y))**2
      post_new = augmented_log_posterior(y, w * h_new(:size(y)), candidate, beta, nu, prior, prior_precision)
      proposal_new = log_mvn2_density(candidate, mean, covariance)
      call regression_posterior_2(y**2 / w, x, tau_new, prior%mu_alpha, prior_precision, &
         covariance_new, mean_new, ok)
      if (.not. ok) return
      proposal_old = log_mvn2_density(alpha_old, mean_new, covariance_new)
      ratio = post_new - post_old + proposal_old - proposal_new
      if (log(random_uniform()) < min(0.0_dp, ratio)) then
         alpha_new = candidate
         accepted = .true.
      end if
   end subroutine draw_alpha_mh

   subroutine draw_beta_mh(y, alpha, beta_old, nu, w, prior, prior_precision, beta_new, accepted)
      real(dp), intent(in) :: y(:)
      real(dp), intent(in) :: alpha(2)
      real(dp), intent(in) :: beta_old
      real(dp), intent(in) :: nu
      real(dp), intent(in) :: w(:)
      type(bayesgarch_priors), intent(in) :: prior
      real(dp), intent(in) :: prior_precision(2, 2)
      real(dp), intent(out) :: beta_new
      logical, intent(out) :: accepted
      real(dp) :: candidate
      real(dp) :: variance
      real(dp) :: variance_new
      real(dp) :: mean
      real(dp) :: mean_new
      real(dp) :: h(size(y) + 1)
      real(dp) :: h_new(size(y) + 1)
      real(dp) :: tau(size(y))
      real(dp) :: tau_new(size(y))
      real(dp) :: q(size(y))
      real(dp) :: z(size(y))
      real(dp) :: post_old
      real(dp) :: post_new
      real(dp) :: proposal_new
      real(dp) :: proposal_old
      real(dp) :: ratio
      logical :: ok

      beta_new = beta_old
      accepted = .false.
      call garch11_filter(y, [alpha, beta_old], h)
      post_old = augmented_log_posterior(y, w * h(:size(y)), alpha, beta_old, nu, prior, prior_precision)
      tau = 2.0_dp * h(:size(y))**2
      call beta_regression_data(y, alpha, beta_old, w, z, q)
      call regression_posterior_1(z, q, tau, prior%mu_beta, 1.0_dp / prior%sigma_beta, &
         variance, mean, ok)
      if (.not. ok) return
      candidate = mean + sqrt(variance) * random_normal()
      if (candidate <= 0.0_dp) return

      call garch11_filter(y, [alpha, candidate], h_new)
      tau_new = 2.0_dp * h_new(:size(y))**2
      post_new = augmented_log_posterior(y, w * h_new(:size(y)), alpha, candidate, nu, prior, prior_precision)
      proposal_new = log_normal_density(candidate, mean, variance)
      call beta_regression_data(y, alpha, candidate, w, z, q)
      call regression_posterior_1(z, q, tau_new, prior%mu_beta, 1.0_dp / prior%sigma_beta, &
         variance_new, mean_new, ok)
      if (.not. ok) return
      proposal_old = log_normal_density(beta_old, mean_new, variance_new)
      ratio = post_new - post_old + proposal_old - proposal_new
      if (log(random_uniform()) < min(0.0_dp, ratio)) then
         beta_new = candidate
         accepted = .true.
      end if
   end subroutine draw_beta_mh

   subroutine beta_regression_data(y, alpha, beta, latent_w, z, q)
      real(dp), intent(in) :: y(:)
      real(dp), intent(in) :: alpha(2)
      real(dp), intent(in) :: beta
      real(dp), intent(in) :: latent_w(:)
      real(dp), intent(out) :: z(size(y))
      real(dp), intent(out) :: q(size(y))
      real(dp) :: residual(size(y))

      call filter_w(y, alpha, beta, residual)
      call quasi_difference(y**2 - residual, -beta, q)
      z = residual + q * beta + y**2 * (1.0_dp / latent_w - 1.0_dp)
   end subroutine beta_regression_data

   subroutine regression_posterior_2(y, x, variance, prior_mean, prior_precision, covariance, mean, ok)
      real(dp), intent(in) :: y(:)
      real(dp), intent(in) :: x(:, :)
      real(dp), intent(in) :: variance(:)
      real(dp), intent(in) :: prior_mean(2)
      real(dp), intent(in) :: prior_precision(2, 2)
      real(dp), intent(out) :: covariance(2, 2)
      real(dp), intent(out) :: mean(2)
      logical, intent(out) :: ok
      real(dp) :: precision(2, 2)
      real(dp) :: d(2)
      real(dp) :: determinant
      integer :: i

      if (size(x, 1) /= size(y) .or. size(x, 2) /= 2 .or. size(variance) /= size(y) .or. &
          any(variance <= 0.0_dp)) then
         covariance = 0.0_dp
         mean = 0.0_dp
         ok = .false.
         return
      end if
      precision = prior_precision
      d = matmul(prior_precision, prior_mean)
      do i = 1, size(y)
         precision(1, 1) = precision(1, 1) + x(i, 1) * x(i, 1) / variance(i)
         precision(1, 2) = precision(1, 2) + x(i, 1) * x(i, 2) / variance(i)
         precision(2, 1) = precision(2, 1) + x(i, 2) * x(i, 1) / variance(i)
         precision(2, 2) = precision(2, 2) + x(i, 2) * x(i, 2) / variance(i)
         d = d + x(i, :) * y(i) / variance(i)
      end do
      call inverse_2x2(precision, covariance, determinant, ok)
      if (ok) mean = matmul(covariance, d)
   end subroutine regression_posterior_2

   subroutine regression_posterior_1(y, x, variance_data, prior_mean, prior_precision, &
      posterior_variance, posterior_mean, ok)
      real(dp), intent(in) :: y(:)
      real(dp), intent(in) :: x(:)
      real(dp), intent(in) :: variance_data(:)
      real(dp), intent(in) :: prior_mean
      real(dp), intent(in) :: prior_precision
      real(dp), intent(out) :: posterior_variance
      real(dp), intent(out) :: posterior_mean
      logical, intent(out) :: ok
      real(dp) :: precision

      if (size(x) /= size(y) .or. size(variance_data) /= size(y) .or. any(variance_data <= 0.0_dp)) then
         posterior_variance = 0.0_dp
         posterior_mean = 0.0_dp
         ok = .false.
         return
      end if
      precision = sum(x * x / variance_data) + prior_precision
      ok = precision > 0.0_dp
      if (.not. ok) then
         posterior_variance = 0.0_dp
         posterior_mean = 0.0_dp
         return
      end if
      posterior_variance = 1.0_dp / precision
      posterior_mean = posterior_variance * (sum(x * y / variance_data) + prior_precision * prior_mean)
   end subroutine regression_posterior_1

   pure function nu_root_function(rate, n, delta, phi) result(value)
      real(dp), intent(in) :: rate
      integer, intent(in) :: n
      real(dp), intent(in) :: delta
      real(dp), intent(in) :: phi
      real(dp) :: value
      real(dp) :: numerator
      real(dp) :: denominator

      if (rate <= 0.0_dp) then
         value = huge(1.0_dp)
         return
      end if
      numerator = 1.0_dp + rate * delta
      denominator = 1.0_dp + rate * (delta - 2.0_dp)
      value = 0.5_dp * real(n, dp) * (log(denominator / (2.0_dp * rate)) + &
         numerator / denominator - digamma(numerator / (2.0_dp * rate))) + rate - phi
   end function nu_root_function

   pure function nu_log_acceptance(n, rate, delta, nu, phi) result(value)
      integer, intent(in) :: n
      real(dp), intent(in) :: rate
      real(dp), intent(in) :: delta
      real(dp), intent(in) :: nu
      real(dp), intent(in) :: phi
      real(dp) :: value
      real(dp) :: u
      real(dp) :: v
      real(dp) :: w

      if (rate <= 0.0_dp .or. nu <= 2.0_dp) then
         value = -huge(1.0_dp)
         return
      end if
      u = (1.0_dp + rate * delta) / (2.0_dp * rate)
      v = 0.5_dp * nu
      w = (1.0_dp + rate * (delta - 2.0_dp)) / (2.0_dp * rate)
      value = real(n, dp) * (log_gamma(u) - log_gamma(v))
      value = value + real(n, dp) * v * log(0.5_dp * (nu - 2.0_dp))
      value = value - real(n, dp) * u * log(w)
      value = value + (nu - delta) * (rate - phi) + phi / rate - 1.0_dp
   end function nu_log_acceptance

   subroutine solve_nu_proposal_rate(n, delta, phi, rate, ok, lower, upper)
      integer, intent(in) :: n
      real(dp), intent(in) :: delta
      real(dp), intent(in) :: phi
      real(dp), intent(out) :: rate
      logical, intent(out) :: ok
      real(dp), intent(in), optional :: lower
      real(dp), intent(in), optional :: upper
      real(dp) :: a
      real(dp) :: b
      real(dp) :: c
      real(dp) :: fa
      real(dp) :: fb
      real(dp) :: fc
      integer :: iteration

      a = 1.0e-5_dp
      b = 500.0_dp
      if (present(lower)) a = lower
      if (present(upper)) b = upper
      fa = nu_root_function(a, n, delta, phi)
      fb = nu_root_function(b, n, delta, phi)
      ok = fa * fb <= 0.0_dp
      if (.not. ok) then
         rate = 0.0_dp
         return
      end if

      do iteration = 1, 200
         c = 0.5_dp * (a + b)
         fc = nu_root_function(c, n, delta, phi)
         if (abs(fc) <= 1.0e-12_dp .or. abs(b - a) <= 1.0e-12_dp * max(1.0_dp, abs(c))) exit
         if (fa * fc <= 0.0_dp) then
            b = c
            fb = fc
         else
            a = c
            fa = fc
         end if
      end do
      rate = 0.5_dp * (a + b)
      ok = rate > a .and. rate < b
      if (.not. ok) ok = rate > 0.0_dp
   end subroutine solve_nu_proposal_rate

   subroutine draw_nu_conditional(w, nu_old, lambda, delta, nu_new, changed, root_failed)
      real(dp), intent(in) :: w(:)
      real(dp), intent(in) :: nu_old
      real(dp), intent(in) :: lambda
      real(dp), intent(in) :: delta
      real(dp), intent(out) :: nu_new
      logical, intent(out) :: changed
      logical, intent(out) :: root_failed
      real(dp) :: phi
      real(dp) :: rate
      real(dp) :: candidate
      real(dp) :: log_probability
      logical :: ok
      integer :: k

      nu_new = nu_old
      changed = .false.
      root_failed = .false.
      phi = 0.5_dp * sum(log(w) + 1.0_dp / w) + lambda
      call solve_nu_proposal_rate(size(w), delta, phi, rate, ok)
      if (.not. ok) then
         root_failed = .true.
         return
      end if
      do k = 1, 50000
         candidate = delta + random_exponential(rate)
         log_probability = nu_log_acceptance(size(w), rate, delta, candidate, phi)
         if (log(random_uniform()) < min(0.0_dp, log_probability)) then
            nu_new = candidate
            changed = abs(nu_new - nu_old) > 0.0_dp
            return
         end if
      end do
   end subroutine draw_nu_conditional

end module bayesgarch_sampler
