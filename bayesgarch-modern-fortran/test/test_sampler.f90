! SPDX-License-Identifier: GPL-2.0-or-later
! Copyright (C) 2008-2021 David Ardia
! Modern Fortran translation of computational routines from bayesGARCH.
!
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 2 of the License, or
! (at your option) any later version.
program test_sampler
   use bayesgarch_garch, only : filter_alpha, garch11_filter, simulate_garch11_student
   use bayesgarch_kinds, only : dp
   use bayesgarch_rng, only : seed_rng
   use bayesgarch_sample, only : form_posterior_sample, posterior_mean, posterior_sd
   use bayesgarch_sampler, only : augmented_log_posterior, bayesgarch_control, &
      bayesgarch_priors, bayesgarch_result, nu_log_acceptance, nu_root_function, &
      regression_posterior_1, regression_posterior_2, run_bayesgarch, &
      solve_nu_proposal_rate
   use test_support, only : assert_close, assert_true, test_custom_constraint
   implicit none

   real(dp), parameter :: y_small(5) = [0.2_dp, -0.1_dp, 0.3_dp, -0.4_dp, 0.25_dp]
   real(dp), parameter :: theta_small(3) = [0.05_dp, 0.12_dp, 0.8_dp]
   real(dp), parameter :: latent_w(5) = [1.2_dp, 0.9_dp, 1.1_dp, 0.8_dp, 1.3_dp]
   real(dp) :: h(6)
   real(dp) :: x(5, 2)
   real(dp) :: tau(5)
   real(dp) :: target(5)
   real(dp) :: covariance(2, 2)
   real(dp) :: mean2(2)
   real(dp) :: variance1
   real(dp) :: mean1
   real(dp) :: rate
   real(dp) :: phi
   real(dp) :: posterior
   real(dp) :: prior_precision(2, 2)
   real(dp) :: q(5)
   real(dp) :: z(5)
   real(dp) :: y(160)
   real(dp) :: hsim(161)
   real(dp), allocatable :: sample(:, :)
   real(dp) :: means(4)
   real(dp) :: sds(4)
   type(bayesgarch_priors) :: prior
   type(bayesgarch_control) :: control
   type(bayesgarch_result) :: result
   type(bayesgarch_result) :: constrained_result
   type(bayesgarch_control) :: constrained_control
   logical :: ok
   integer :: i

   prior = bayesgarch_priors()
   prior_precision = 0.0_dp
   prior_precision(1, 1) = 0.001_dp
   prior_precision(2, 2) = 0.001_dp
   call garch11_filter(y_small, theta_small, h)
   posterior = augmented_log_posterior(y_small, latent_w * h(:5), theta_small(:2), &
      theta_small(3), 8.0_dp, prior, prior_precision)
   call assert_close(posterior, 3.786714886825154_dp, 1.0e-12_dp, "augmented posterior")

   call filter_alpha(y_small, theta_small(3), x)
   tau = 2.0_dp * h(:5)**2
   target = y_small**2 / latent_w
   call regression_posterior_2(target, x, tau, [0.0_dp, 0.0_dp], prior_precision, &
      covariance, mean2, ok)
   call assert_true(ok, "two-parameter regression posterior exists")
   call assert_close(covariance, reshape([0.0026240998150997_dp, -0.0514244950886562_dp, &
      -0.0514244950886562_dp, 1.7895557842938894_dp], [2, 2]), 1.0e-12_dp, &
      "two-parameter posterior covariance")
   call assert_close(mean2, [0.0324441370878436_dp, -0.0443330523860012_dp], &
      1.0e-12_dp, "two-parameter posterior mean")

   q = [0.0_dp, 0.05_dp, 0.1348_dp, 0.23488_dp, 0.350336_dp]
   z = [-0.0166666666666667_dp, -0.0436888888888889_dp, 0.0626181818181818_dp, &
      0.225472_dp, 0.129200123076923_dp]
   call regression_posterior_1(z, q, tau, 0.0_dp, 0.001_dp, variance1, mean1, ok)
   call assert_true(ok, "one-parameter regression posterior exists")
   call assert_close(variance1, 0.303406387021123_dp, 1.0e-12_dp, &
      "one-parameter posterior variance")
   call assert_close(mean1, 0.520111732996687_dp, 1.0e-12_dp, &
      "one-parameter posterior mean")

   phi = 0.5_dp * sum(log(latent_w) + 1.0_dp / latent_w) + 0.01_dp
   call solve_nu_proposal_rate(5, 2.0_dp, phi, rate, ok)
   call assert_true(ok, "nu proposal root exists")
   call assert_close(rate, 0.014841567001108_dp, 1.0e-11_dp, "nu proposal rate")
   call assert_close(nu_root_function(rate, 5, 2.0_dp, phi), 0.0_dp, 2.0e-11_dp, &
      "nu root equation")
   call assert_close(nu_log_acceptance(5, rate, 2.0_dp, 8.0_dp, phi), &
      -3.883789865725476_dp, 2.0e-10_dp, "nu log acceptance")

   call seed_rng(192837)
   call simulate_garch11_student(size(y), [0.03_dp, 0.08_dp, 0.85_dp], 8.0_dp, y, hsim, burn=300)
   control = bayesgarch_control()
   control%n_chains = 2
   control%n_iter = 220
   control%start = [0.04_dp, 0.10_dp, 0.75_dp, 10.0_dp]
   control%enforce_stationarity = .true.
   call run_bayesgarch(y, result, control=control)

   call assert_true(all(result%draws(:, 1:3, :) > 0.0_dp), "sampled GARCH parameters are positive")
   call assert_true(all(result%draws(:, 4, :) > 2.0_dp), "sampled nu exceeds two")
   call assert_true(all(result%draws(:, 2, :) + result%draws(:, 3, :) < 1.0_dp), &
      "stationarity constraint is enforced")
   call assert_true(all(result%alpha_accept > 0), "alpha block moves in every chain")
   call assert_true(all(result%beta_accept > 0), "beta block moves in every chain")
   call assert_true(all(result%nu_updates > 0), "nu block moves in every chain")
   call assert_true(all(result%nu_root_fail == 0), "nu proposal roots are found")

   call form_posterior_sample(result, burn=60, thin=4, sample=sample)
   call assert_true(size(sample, 1) == 80, "posterior sample row count")
   call assert_true(size(sample, 2) == 4, "posterior sample column count")
   call assert_close(sample(1, :), result%draws(61, :, 1), 0.0_dp, &
      "posterior sample first retained draw")
   call assert_close(sample(41, :), result%draws(61, :, 2), 0.0_dp, &
      "posterior sample chain concatenation")
   means = posterior_mean(sample)
   sds = posterior_sd(sample)
   call assert_true(all(means > 0.0_dp), "posterior means are positive")
   call assert_true(all(sds > 0.0_dp), "posterior standard deviations are positive")

   do i = 1, 4
      call assert_true(means(i) < huge(1.0_dp), "posterior means are finite")
   end do

   constrained_control = bayesgarch_control()
   constrained_control%n_iter = 80
   constrained_control%start = [0.04_dp, 0.10_dp, 0.75_dp, 10.0_dp]
   call run_bayesgarch(y, constrained_result, control=constrained_control, constraint=test_custom_constraint)
   call assert_true(all(constrained_result%draws(:, 1, 1) < 0.06_dp), &
      "custom prior constraint is enforced")
   call assert_true(all(constrained_result%draws(:, 4, 1) < 50.0_dp), &
      "custom nu constraint is enforced")
   call assert_true(constrained_result%constraint_reject(1) > 0, &
      "custom prior constraint rejects at least one block")

   write(*, '(a)') "Bayesian sampler and posterior-sample tests passed."

end program test_sampler
