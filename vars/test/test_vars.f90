! SPDX-License-Identifier: GPL-2.0-or-later
! Computational translation of R package vars 1.6-1; see NOTICE.md and UPSTREAM.md.
program test_vars
   use r_kinds, only : dp
   use vars, only : var_model, var_selection_result, forecast_result, vars_test_result, svar_result, svec_result
   use vars, only : vars_success, var_const, fit_var, var_select, phi_from_a, var_roots
   use vars, only : forecast_var, residual_covariance_unbiased, bq_identification, fevd_var
   use vars, only : arch_test_univariate, arch_test_multivariate
   use vars, only : jarque_bera_univariate, jarque_bera_multivariate
   use vars, only : portmanteau_tests, bg_serial_tests, granger_causality, instantaneous_causality
   use vars, only : svar_fit_scoring, svar_negloglik, vec2var_coefficients, var_loglik
   use vars, only : residual_bootstrap_path, svec_long_run_matrix, svec_fit_scoring
   implicit none

   type(var_model) :: model
   type(var_selection_result) :: selection
   type(forecast_result) :: forecast
   type(vars_test_result) :: arch1, arch_multi, jb1, jb_multi, jb_skew, jb_kurt
   type(vars_test_result) :: pt1, pt2, bg_lm, bg_es, granger, instant
   type(svar_result) :: svar
   type(svec_result) :: svec
   real(dp) :: y(20, 2), innov(19, 2), atrue(2, 2), ctrue(2)
   real(dp), allocatable :: phi(:, :, :), omega(:, :, :), sigma(:, :)
   real(dp), allocatable :: lr(:, :), sr(:, :), roots_r(:), roots_i(:), roots_m(:)
   real(dp), allocatable :: vec_a(:, :, :), ysampled(:, :), xi(:, :)
   real(dp) :: sigma_s(2, 2), a_fixed(2, 2), b_fixed(2, 2), nll, loglik
   real(dp) :: pi_matrix(2, 2), gamma_blocks(2, 2, 2), gamma_zero(2, 2, 1)
   real(dp) :: alpha(2, 1), beta(2, 1)
   logical :: free_a(2, 2), free_b(2, 2), free_lr(2, 2), free_sr(2, 2), cause(2)
   integer :: info, i, nparams, indices(19)

   atrue = reshape([0.5_dp, -0.2_dp, 0.1_dp, 0.4_dp], [2, 2])
   ctrue = [1.0_dp, -0.5_dp]
   innov(:, 1) = [0.10_dp, -0.08_dp, 0.05_dp, -0.02_dp, 0.06_dp, -0.04_dp, &
      0.03_dp, 0.00_dp, -0.05_dp, 0.04_dp, -0.01_dp, 0.02_dp, -0.03_dp, &
      0.05_dp, -0.02_dp, 0.01_dp, -0.04_dp, 0.03_dp, -0.01_dp]
   innov(:, 2) = [-0.05_dp, 0.04_dp, 0.02_dp, -0.03_dp, 0.01_dp, 0.05_dp, &
      -0.02_dp, 0.03_dp, -0.01_dp, 0.02_dp, -0.04_dp, 0.01_dp, 0.00_dp, &
      -0.02_dp, 0.04_dp, -0.01_dp, 0.02_dp, 0.00_dp, -0.03_dp]
   y = 0.0_dp
   y(1, :) = [0.4_dp, -0.3_dp]
   do i = 2, 20
      y(i, :) = ctrue + matmul(atrue, y(i - 1, :)) + innov(i - 1, :)
   end do

   call fit_var(y, 1, var_const, model, info)
   call check_int("fit_var info", info, vars_success)
   call check_close("coef 11", model%coef(1, 1), 0.2230426009_dp, 1.0e-8_dp)
   call check_close("coef 12", model%coef(1, 2), -0.1477675927_dp, 1.0e-8_dp)
   call check_close("coef c1", model%coef(1, 3), 1.1377985098_dp, 1.0e-8_dp)
   call check_close("coef 21", model%coef(2, 1), -0.0915708467_dp, 1.0e-8_dp)
   call check_close("coef 22", model%coef(2, 2), 0.4913892744_dp, 1.0e-8_dp)
   call check_close("coef c2", model%coef(2, 3), -0.5576431812_dp, 1.0e-8_dp)

   call var_select(y, 3, var_const, selection, info)
   call check_int("var_select info", info, vars_success)
   do i = 1, 4
      call check_int("var_select lag", selection%selection(i), 2)
   end do
   call check_close("AIC lag2", selection%criteria(1, 2), -14.5262526_dp, 2.0e-7_dp)

   call phi_from_a(model%a, 2, phi, info)
   call check_int("phi info", info, vars_success)
   call check_close("phi2 11", phi(1, 1, 3), 0.06327920_dp, 2.0e-8_dp)
   call check_close("phi2 22", phi(2, 2, 3), 0.25499462_dp, 2.0e-8_dp)

   call var_roots(model%a, roots_r, roots_i, roots_m, info)
   call check_int("roots info", info, vars_success)
   call check_close("root min", minval(roots_m), 0.17963866_dp, 2.0e-8_dp)
   call check_close("root max", maxval(roots_m), 0.53479320_dp, 2.0e-8_dp)

   call forecast_var(model, 2, 0.95_dp, forecast, info)
   call check_int("forecast info", info, vars_success)
   call check_close("forecast 1 1", forecast%point(1, 1), 1.73270410_dp, 2.0e-8_dp)
   call check_close("forecast 2 2", forecast%point(2, 2), -1.41376320_dp, 2.0e-8_dp)
   call check_close("forecast se 2 1", forecast%se(2, 1), sqrt(0.00146685_dp), 2.0e-7_dp)

   call var_loglik(model, loglik, nparams, info)
   call check_int("loglik info", info, vars_success)
   call check_int("loglik parameters", nparams, 6)
   call check_close("VAR loglik", loglik, 80.47169710_dp, 3.0e-7_dp)

   call residual_covariance_unbiased(model, sigma, info)
   call check_int("sigma info", info, vars_success)
   call check_close("sigma 11", sigma(1, 1), 0.00136820_dp, 2.0e-8_dp)
   call bq_identification(model%a, sigma, lr, sr, info)
   call check_int("BQ info", info, vars_success)
   call check_close("BQ LR 11", lr(1, 1), 0.05254320_dp, 2.0e-8_dp)
   call check_close("BQ SR 22", sr(2, 2), 0.02552624_dp, 2.0e-8_dp)

   call fevd_var(model, 2, omega, info)
   call check_int("FEVD info", info, vars_success)
   call check_close("FEVD h1 21", omega(2, 1, 1), 0.04108473_dp, 2.0e-8_dp)
   call check_close("FEVD h2 12", omega(1, 2, 2), 0.01101946_dp, 2.0e-8_dp)

   call arch_test_univariate(model%resid(:, 1), 2, arch1, info)
   call check_int("ARCH info", info, vars_success)
   call check_close("ARCH statistic", arch1%statistic, 0.41108721_dp, 2.0e-7_dp)
   call arch_test_multivariate(model%resid, 1, arch_multi, info)
   call check_int("multivariate ARCH info", info, vars_success)
   call check_close("multivariate ARCH", arch_multi%statistic, 9.69723432_dp, 3.0e-7_dp)
   call jarque_bera_univariate(model%resid(:, 1), jb1, info)
   call check_int("JB info", info, vars_success)
   call check_close("JB statistic", jb1%statistic, 0.51373454_dp, 2.0e-7_dp)
   call jarque_bera_multivariate(model%resid, jb_multi, jb_skew, jb_kurt, info)
   call check_int("multivariate JB info", info, vars_success)
   call check_close("multivariate skewness", jb_skew%statistic, 0.73863931_dp, 2.0e-7_dp)
   call check_close("multivariate kurtosis", jb_kurt%statistic, 0.95399718_dp, 2.0e-7_dp)
   call check_close("multivariate JB", jb_multi%statistic, 1.69263650_dp, 3.0e-7_dp)

   call portmanteau_tests(model, 2, pt1, pt2, info)
   call check_int("PT info", info, vars_success)
   call check_close("PT statistic", pt1%statistic, 12.38196057_dp, 2.0e-7_dp)
   call check_close("PT adjusted", pt2%statistic, 13.35631156_dp, 2.0e-7_dp)
   call bg_serial_tests(model, 1, bg_lm, bg_es, info)
   call check_int("BG info", info, vars_success)
   call check_close("BG statistic", bg_lm%statistic, 9.62442830_dp, 3.0e-7_dp)
   call check_close("ES statistic", bg_es%statistic, 2.38118994_dp, 3.0e-7_dp)

   cause = [.true., .false.]
   call granger_causality(model, cause, granger, info)
   call check_int("Granger info", info, vars_success)
   call check_close("Granger statistic", granger%statistic, 1.63205569_dp, 2.0e-6_dp)
   call instantaneous_causality(model, cause, instant, info)
   call check_int("instant info", info, vars_success)
   call check_close("instant statistic", instant%statistic, 0.74980433_dp, 2.0e-7_dp)

   do i = 1, 19
      indices(i) = i
   end do
   call residual_bootstrap_path(model, indices, ysampled, info)
   call check_int("bootstrap path info", info, vars_success)
   call check_close("bootstrap identity path", maxval(abs(ysampled - y)), 0.0_dp, 1.0e-10_dp)

   pi_matrix = reshape([-0.2_dp, 0.05_dp, 0.1_dp, -0.1_dp], [2, 2])
   gamma_blocks(:, :, 1) = reshape([0.3_dp, 0.2_dp, 0.1_dp, 0.25_dp], [2, 2])
   gamma_blocks(:, :, 2) = reshape([0.1_dp, 0.15_dp, -0.05_dp, 0.05_dp], [2, 2])
   call vec2var_coefficients(pi_matrix, gamma_blocks, "transitory", vec_a, info)
   call check_int("vec2var transitory info", info, vars_success)
   call check_close("vec2var transitory A1(1,1)", vec_a(1, 1, 1), 1.1_dp, 1.0e-12_dp)
   call check_close("vec2var transitory A2(1,2)", vec_a(1, 2, 2), -0.15_dp, 1.0e-12_dp)
   call check_close("vec2var transitory A3(2,1)", vec_a(2, 1, 3), -0.15_dp, 1.0e-12_dp)
   call vec2var_coefficients(pi_matrix, gamma_blocks, "longrun", vec_a, info)
   call check_int("vec2var longrun info", info, vars_success)
   call check_close("vec2var longrun A1(1,1)", vec_a(1, 1, 1), 1.3_dp, 1.0e-12_dp)
   call check_close("vec2var longrun A3(1,2)", vec_a(1, 2, 3), 0.15_dp, 1.0e-12_dp)

   alpha(:, 1) = [-0.2_dp, 0.1_dp]
   beta(:, 1) = [1.0_dp, -1.0_dp]
   gamma_zero = 0.0_dp
   call svec_long_run_matrix(alpha, beta, gamma_zero, xi, info)
   call check_int("SVEC Xi info", info, vars_success)
   call check_close("SVEC Xi 11", xi(1, 1), 1.0_dp / 3.0_dp, 1.0e-10_dp)
   call check_close("SVEC Xi 12", xi(1, 2), 2.0_dp / 3.0_dp, 1.0e-10_dp)

   sigma_s = reshape([1.0_dp, 0.3_dp, 0.3_dp, 0.8_dp], [2, 2])
   free_lr = .true.
   free_sr = .true.
   free_sr(1, 2) = .false.
   call svec_fit_scoring(alpha, beta, gamma_zero, sigma_s, 100, free_lr, free_sr, svec, info, &
      max_iter = 100, conv_crit = 1.0e-10_dp)
   call check_int("SVEC scoring info", info, vars_success)
   call check_close("SVEC SR11", svec%sr(1, 1), 1.0_dp, 2.0e-7_dp)
   call check_close("SVEC SR21", svec%sr(2, 1), 0.3_dp, 2.0e-7_dp)
   call check_close("SVEC SR22", svec%sr(2, 2), sqrt(0.71_dp), 2.0e-7_dp)
   call check_close("SVEC LR11", svec%lr(1, 1), 0.5333333333_dp, 2.0e-7_dp)

   a_fixed = 0.0_dp
   a_fixed(1, 1) = 1.0_dp
   a_fixed(2, 2) = 1.0_dp
   b_fixed = 0.0_dp
   free_a = .false.
   free_b = .false.
   free_b(1, 1) = .true.
   free_b(2, 1) = .true.
   free_b(2, 2) = .true.
   call svar_fit_scoring(sigma_s, 100, a_fixed, b_fixed, free_a, free_b, svar, info, &
      start = [1.0_dp, 0.3_dp, 0.84_dp], max_iter = 100, conv_crit = 1.0e-10_dp)
   call check_int("SVAR info", info, vars_success)
   call check_close("SVAR B11", svar%b(1, 1), 1.0_dp, 2.0e-7_dp)
   call check_close("SVAR B21", svar%b(2, 1), 0.3_dp, 2.0e-7_dp)
   call check_close("SVAR B22", svar%b(2, 2), sqrt(0.71_dp), 2.0e-7_dp)
   call svar_negloglik(svar%a, svar%b, sigma_s, 100, nll, info)
   call check_int("SVAR nll info", info, vars_success)
   call check_close("SVAR nll", nll, 266.66319119_dp, 2.0e-6_dp)

   write (*, '(a)') "All vars tests passed."

contains

   subroutine check_close(name, actual, expected, tolerance)
      character(len = *), intent(in) :: name
      real(dp), intent(in) :: actual, expected, tolerance

      if (abs(actual - expected) > tolerance) then
         write (*, '(a,2es24.14)') trim(name) // " failed: ", actual, expected
         error stop 1
      end if
   end subroutine check_close

   subroutine check_int(name, actual, expected)
      character(len = *), intent(in) :: name
      integer, intent(in) :: actual, expected

      if (actual /= expected) then
         write (*, '(a,2i12)') trim(name) // " failed: ", actual, expected
         error stop 1
      end if
   end subroutine check_int

end program test_vars
