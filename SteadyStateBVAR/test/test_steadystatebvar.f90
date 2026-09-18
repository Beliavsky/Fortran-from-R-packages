! SPDX-License-Identifier: GPL-3.0-or-later
! Derived from SteadyStateBVAR 0.2.0, Copyright (c) 2026 Mark Becker.
program test_steadystatebvar
   use steadystatebvar, only : dp, bvar_model, forecast_result, irf_result, bvar_summary_result
   use steadystatebvar, only : bvar_create, setup_model, set_priors, restrict_beta
   use steadystatebvar, only : fit_model, forecast_model, conditional_forecast_model, irf_model
   use steadystatebvar, only : ppi, bvar_summary, ssbvar_success, ssbvar_invalid_input
   use steadystatebvar, only : ssbvar_unsupported
   implicit none

   integer, parameter :: n = 48, k = 2
   real(dp) :: data(n, k), restriction(k, k), mean_value, variance
   real(dp) :: own_mean(k), theta_psi(k), omega_psi(k, k), dummy(n)
   real(dp) :: theta_a(1), omega_a(1, 1), theta_log_lambda_1(k), omega_log_lambda_1(k, k)
   real(dp) :: alpha_phi(k), beta_phi(k), theta_gamma_0(k), omega_gamma_0(k, k)
   real(dp) :: theta_gamma_1(k), omega_gamma_1(k, k), v_phi(k, k)
   real(dp) :: expected_d(3, 2)
   integer :: t, info
   type(bvar_model) :: model, model_repeat, model_trend, model_sv, model_sv_ar1
   type(forecast_result) :: fcst, fcst_growth, cfcst
   type(irf_result) :: irf, girf
   type(bvar_summary_result) :: summary

   do t = 1, n
      data(t, 1) = 0.6_dp * sin(0.19_dp * real(t, dp)) + 0.02_dp * real(t, dp) + &
         0.03_dp * sin(0.017_dp * real(t * t, dp))
      data(t, 2) = 0.4_dp * cos(0.13_dp * real(t, dp)) + 0.25_dp * data(t, 1) + &
         0.02_dp * cos(0.011_dp * real(t * t, dp))
      dummy(t) = merge(1.0_dp, 0.0_dp, t > n / 2)
   end do

   call bvar_create(data, model, info)
   call check(info == ssbvar_success, "bvar_create")
   call setup_model(model, 1, "constant", info=info)
   call check(info == ssbvar_success, "setup_model")
   call check(model%setup%n == n - 1, "setup N")
   call check(model%setup%q == 1, "setup q")
   call check(all(shape(model%setup%beta_ols) == [k, k]), "beta OLS shape")

   own_mean = [0.7_dp, 0.7_dp]
   theta_psi = 0.0_dp
   omega_psi = 0.0_dp
   omega_psi(1, 1) = 0.5_dp
   omega_psi(2, 2) = 0.5_dp
   call set_priors(model, lambda_1=0.25_dp, lambda_2=0.5_dp, lambda_3=0.0_dp, info=info)
   call check(info == ssbvar_invalid_input, "lambda_3 must be positive")
   call set_priors(model, lambda_1=0.25_dp, lambda_2=0.5_dp, lambda_3=1.0_dp, &
      first_own_lag_prior_mean=own_mean, theta_psi=theta_psi, omega_psi=omega_psi, &
      jeffreys=.false., info=info)
   call check(info == ssbvar_success, "set_priors")
   call check(model%priors%m == k + 2, "inverse-Wishart df")

   restriction = 1.0_dp
   restriction(1, 2) = 0.0_dp
   call restrict_beta(model, restriction, info)
   call check(info == ssbvar_success, "restrict_beta")
   call check(abs(model%priors%omega_beta(3, 3) - 1.0e-7_dp) < 1.0e-15_dp, "restriction variance")

   call ppi(1.7_dp, 2.3_dp, interval=0.95_dp, annualized_growthrate=.true., freq=4, &
      mean_value=mean_value, variance=variance, info=info)
   call check(info == ssbvar_success, "ppi")
   call check(abs(mean_value - 0.5_dp) < 1.0e-12_dp, "ppi mean")
   call check(variance > 0.0_dp, "ppi variance")

   model_repeat = model
   call fit_model(model, h=4, n_draws=40, burnin=30, seed=20260915, info=info)
   call check(info == ssbvar_success, "fit_model")
   call check(all(shape(model%fit%beta) == [k, k, 40]), "beta draws shape")
   call check(all(shape(model%fit%y_pred) == [40, 4, k]), "predictive draws shape")
   call fit_model(model_repeat, h=4, n_draws=40, burnin=30, seed=20260915, info=info)
   call check(info == ssbvar_success, "fit_model reproducibility run")
   call check(maxval(abs(model%fit%beta - model_repeat%fit%beta)) < 1.0e-15_dp, "deterministic beta RNG")
   call check(maxval(abs(model%fit%y_pred - model_repeat%fit%y_pred)) < 1.0e-15_dp, "deterministic forecast RNG")

   call forecast_model(model, pi=0.9_dp, result=fcst, info=info)
   call check(info == ssbvar_success, "forecast_model")
   call check(all(shape(fcst%forecast) == [4, k]), "forecast shape")
   call check(all(fcst%lower <= fcst%upper), "forecast interval ordering")
   call forecast_model(model, pi=0.9_dp, use_median=.true., growth_rate_idx=[1], freq=4, &
      result=fcst_growth, info=info)
   call check(info == ssbvar_success, "forecast annual growth")
   call check(all(fcst_growth%lower <= fcst_growth%upper), "annual forecast interval ordering")

   call conditional_forecast_model(model, [1], [2], [0.15_dp], pi=0.9_dp, result=cfcst, info=info)
   call check(info == ssbvar_success, "conditional_forecast_model")
   call check(abs(cfcst%forecast(2, 1) - 0.15_dp) < 5.0e-8_dp, "conditional equality")

   call irf_model(model, h=5, orthogonal=.true., use_median=.false., ci=0.9_dp, result=irf, info=info)
   call check(info == ssbvar_success, "irf_model")
   call check(all(shape(irf%center) == [k, k, 6]), "IRF shape")
   call check(all(irf%lower <= irf%upper), "IRF interval ordering")
   call irf_model(model, h=5, orthogonal=.false., use_median=.true., ci=0.9_dp, &
      growth_rate_idx=[1], freq=4, result=girf, info=info)
   call check(info == ssbvar_success, "generalized annualized IRF")
   call check(all(girf%lower <= girf%upper), "GIRF interval ordering")

   call bvar_summary(model, result=summary, info=info)
   call check(info == ssbvar_success, "bvar_summary")
   call check(all(shape(summary%sigma_u) == [k, k]), "summary covariance shape")
   call bvar_summary(model, use_median=.true., result=summary, info=info)
   call check(info == ssbvar_success, "bvar_summary median")

   call bvar_create(data, model_trend, info)
   call check(info == ssbvar_success, "trend bvar_create")
   call setup_model(model_trend, 2, "constant_and_trend", info=info)
   call check(info == ssbvar_success, "trend setup")
   call check(model_trend%setup%q == 2, "trend q")
   call check(abs(model_trend%setup%dt(n, 2) - real(n, dp)) < 1.0e-12_dp, "trend endpoint")
   call set_priors(model_trend, jeffreys=.true., info=info)
   call check(info == ssbvar_success, "trend priors")
   call fit_model(model_trend, h=3, n_draws=8, burnin=6, seed=101, info=info)
   call check(info == ssbvar_success, "trend Jeffreys fit")
   expected_d(:, 1) = 1.0_dp
   expected_d(:, 2) = [real(n + 1, dp), real(n + 2, dp), real(n + 3, dp)]
   call check(maxval(abs(model_trend%fit%d_pred - expected_d)) < 1.0e-12_dp, "future trend construction")

   theta_a = 0.0_dp
   omega_a = 0.0_dp
   omega_a(1, 1) = 4.0_dp
   theta_log_lambda_1 = -1.0_dp
   omega_log_lambda_1 = 0.0_dp
   omega_log_lambda_1(1, 1) = 1.0_dp
   omega_log_lambda_1(2, 2) = 1.0_dp
   alpha_phi = 4.0_dp
   beta_phi = 0.08_dp

   call bvar_create(data, model_sv, info)
   call check(info == ssbvar_success, "SV bvar_create")
   call setup_model(model_sv, 1, "constant_and_dummy", dummy=dummy, info=info)
   call check(info == ssbvar_success, "dummy setup")
   call check(maxval(abs(model_sv%setup%dt(:, 2) - dummy)) < 1.0e-15_dp, "dummy deterministic column")
   call set_priors(model_sv, sv=.true., sv_type="BAD", info=info)
   call check(info == ssbvar_invalid_input, "invalid SV type")
   call set_priors(model_sv, sv=.true., sv_type="RW", info=info)
   call check(info == ssbvar_invalid_input, "RW SV priors are required")
   call set_priors(model_sv, sv=.true., sv_type="RW", theta_a=theta_a, omega_a=omega_a, &
      theta_log_lambda_1=theta_log_lambda_1, omega_log_lambda_1=omega_log_lambda_1, &
      alpha_phi=alpha_phi, beta_phi=beta_phi, info=info)
   call check(info == ssbvar_success, "RW SV prior surface")
   call check(.not. model_sv%priors%jeffreys, "SV disables Jeffreys covariance prior")
   call fit_model(model_sv, h=2, n_draws=6, burnin=6, seed=3, sv_state_scale=0.12_dp, info=info)
   call check(info == ssbvar_success, "RW SV fit")
   call check(.not. model_sv%fit%homoscedastic, "RW SV fit marker")
   call check(all(shape(model_sv%fit%phi) == [k, 6]), "RW phi draw shape")
   call check(all(shape(model_sv%fit%sigma_u_time) == [k, k, n - 1, 6]), "RW covariance draw shape")
   call check(all(shape(model_sv%fit%sigma_u_pred) == [k, k, 2, 6]), "RW forecast covariance shape")
   call check(all(shape(model_sv%fit%sigma_time_mean) == [k, k, n - 1]), "RW covariance mean path shape")
   call check(all(model_sv%fit%phi > 0.0_dp), "RW phi positivity")
   call forecast_model(model_sv, pi=0.8_dp, result=fcst, info=info)
   call check(info == ssbvar_success, "RW SV forecast")
   call irf_model(model_sv, h=3, orthogonal=.true., ci=0.8_dp, result=irf, info=info, t=n)
   call check(info == ssbvar_success, "RW SV time-specific IRF")
   call bvar_summary(model_sv, result=summary, info=info, t=n)
   call check(info == ssbvar_success, "RW SV summary")
   call check(allocated(summary%a), "RW summary A")
   call check(allocated(summary%phi), "RW summary phi")
   call check(all(shape(summary%sigma_u) == [k, k]), "RW summary covariance shape")
   call conditional_forecast_model(model_sv, [1], [1], [0.0_dp], result=cfcst, info=info)
   call check(info == ssbvar_unsupported, "RW conditional forecast matches upstream rejection")

   theta_gamma_0 = -0.1_dp
   theta_gamma_1 = 0.8_dp
   omega_gamma_0 = 0.0_dp
   omega_gamma_1 = 0.0_dp
   v_phi = 0.0_dp
   omega_gamma_0(1, 1) = 0.5_dp
   omega_gamma_0(2, 2) = 0.5_dp
   omega_gamma_1(1, 1) = 0.05_dp
   omega_gamma_1(2, 2) = 0.05_dp
   v_phi(1, 1) = 0.10_dp
   v_phi(2, 2) = 0.10_dp

   call bvar_create(data, model_sv_ar1, info)
   call check(info == ssbvar_success, "AR1 SV bvar_create")
   call setup_model(model_sv_ar1, 1, "constant", info=info)
   call check(info == ssbvar_success, "AR1 SV setup")
   call set_priors(model_sv_ar1, sv=.true., sv_type="AR1", theta_a=theta_a, omega_a=omega_a, &
      theta_log_lambda_1=theta_log_lambda_1, omega_log_lambda_1=omega_log_lambda_1, &
      theta_gamma_0=theta_gamma_0, omega_gamma_0=omega_gamma_0, theta_gamma_1=theta_gamma_1, &
      omega_gamma_1=omega_gamma_1, m_phi=k + 3, v_phi=v_phi, info=info)
   call check(info == ssbvar_success, "AR1 SV prior surface")
   call fit_model(model_sv_ar1, h=2, n_draws=5, burnin=5, seed=17, sv_state_scale=0.10_dp, &
      sv_parameter_scale=0.02_dp, info=info)
   call check(info == ssbvar_success, "AR1 SV fit")
   call check(all(shape(model_sv_ar1%fit%gamma_0) == [k, 5]), "AR1 gamma_0 draw shape")
   call check(all(shape(model_sv_ar1%fit%gamma_1) == [k, 5]), "AR1 gamma_1 draw shape")
   call check(all(abs(model_sv_ar1%fit%gamma_1) < 1.0_dp), "AR1 stationarity bounds")
   call check(all(shape(model_sv_ar1%fit%phi_cov) == [k, k, 5]), "AR1 Phi draw shape")
   call check(all(shape(model_sv_ar1%fit%sigma_time_median) == [k, k, n - 1]), &
      "AR1 covariance median path shape")
   call forecast_model(model_sv_ar1, pi=0.8_dp, result=fcst, info=info)
   call check(info == ssbvar_success, "AR1 SV forecast")
   call irf_model(model_sv_ar1, h=2, orthogonal=.false., ci=0.8_dp, result=girf, info=info, t=n)
   call check(info == ssbvar_success, "AR1 SV generalized IRF")
   call bvar_summary(model_sv_ar1, use_median=.true., result=summary, info=info, t=n)
   call check(info == ssbvar_success, "AR1 SV summary")
   call check(allocated(summary%a), "AR1 summary A")
   call check(allocated(summary%gamma_0), "AR1 summary gamma_0")
   call check(allocated(summary%gamma_1), "AR1 summary gamma_1")
   call check(allocated(summary%phi_cov), "AR1 summary Phi")

   print '(a)', "All SteadyStateBVAR deterministic tests passed."

contains

   subroutine check(condition, label)
      logical, intent(in) :: condition !! Test predicate that must evaluate true.
      character(len=*), intent(in) :: label !! Short diagnostic label printed before stopping on a failed test.
      if (.not. condition) then
         print '(a)', "FAILED: " // trim(label)
         error stop 1
      end if
   end subroutine check

end program test_steadystatebvar
