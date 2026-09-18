! SPDX-License-Identifier: MIT
! SPDX-FileComment: Deterministic nonlinear least-squares parity tests against R stats.
program test_nls
   use stats_test_assertions, only: assert_close, assert_matrix_close, assert_true, &
                                      assert_vector_close
   use r_stats, only: dp, nls_confint_profile, nls_fit, nls_fit_plinear, nls_fit_t
   use r_stats, only: ss_asymp, ss_gompertz, ss_logis, ss_micmen, ss_weibull
   use r_stats, only: ss_micmen_jacobian, ss_micmen_model
   implicit none

   real(dp), allocatable :: intervals(:, :)
   real(dp) :: x_exp(9, 1), x_mm(10, 1), x_plinear(11, 1)
   real(dp) :: y_exp(9), y_mm(10), y_plinear(11)
   type(nls_fit_t) :: fit
   integer :: i

   do i = 1, size(x_exp, 1)
      x_exp(i, 1) = real(i - 1, dp)
   end do
   y_exp = [2.1_dp, 3.0_dp, 4.7_dp, 6.8_dp, 10.2_dp, 14.9_dp, 22.3_dp, 32.8_dp, 49.1_dp]
   fit = nls_fit(x_exp, y_exp, [2.0_dp, 0.35_dp], exponential_model, &
                 jacobian=exponential_jacobian, tolerance=1.0e-10_dp, max_iterations=200)
   call assert_true(fit%converged, "exponential NLS convergence")
   call assert_true(fit%analytic_jacobian, "exponential analytic Jacobian")
   call assert_vector_close(fit%coefficients, &
      [2.07549872916724_dp, 0.395210292036399_dp], 2.0e-7_dp, &
      "exponential NLS coefficients")
   call assert_close(fit%rss, 0.0979256731223441_dp, 2.0e-10_dp, "exponential NLS RSS")
   call assert_close(fit%residual_standard_error, 0.118276717381827_dp, 2.0e-10_dp, &
                     "exponential NLS residual standard error")
   call assert_matrix_close(fit%covariance, reshape( &
      [4.99907845207183e-4_dp, -3.26390551427007e-5_dp, &
       -3.26390551427007e-5_dp, 2.19115372964706e-6_dp], [2, 2]), &
      2.0e-9_dp, "exponential NLS covariance")
   intervals = nls_confint_profile(x_exp, y_exp, fit, exponential_model, &
                                   jacobian=exponential_jacobian)
   call assert_matrix_close(intervals, reshape( &
      [2.02295447760017_dp, 0.391718852015299_dp, &
       2.12886285632919_dp, 0.398730984623080_dp], [2, 2]), &
      2.0e-6_dp, "exponential NLS profile confidence intervals")

   do i = 1, size(x_plinear, 1)
      x_plinear(i, 1) = real(i - 1, dp)
   end do
   y_plinear = 1.5_dp + 4.0_dp*exp(-0.35_dp*x_plinear(:, 1)) + &
                [0.02_dp, -0.03_dp, 0.01_dp, 0.04_dp, -0.02_dp, 0.01_dp, &
                 -0.01_dp, 0.02_dp, -0.02_dp, 0.01_dp, -0.01_dp]
   fit = nls_fit_plinear(x_plinear, y_plinear, [0.2_dp], 2, exponential_basis, &
                         tolerance=1.0e-9_dp, max_iterations=200)
   call assert_true(fit%converged, "partially linear NLS convergence")
   call assert_vector_close(fit%coefficients, &
      [0.348990910556766_dp, 1.49650487146996_dp, 4.00892899398785_dp], &
      2.0e-6_dp, "partially linear NLS coefficients")
   call assert_close(fit%rss, 0.00448142773706529_dp, 2.0e-10_dp, &
                     "partially linear NLS RSS")
   call assert_matrix_close(fit%covariance, reshape( &
      [3.30358351763113e-5_dp, 9.41223949159293e-5_dp, -2.66430251817925e-5_dp, &
       9.41223949159293e-5_dp, 3.70854606938973e-4_dp, -2.47348066668415e-4_dp, &
       -2.66430251817925e-5_dp, -2.47348066668415e-4_dp, 5.89271932674426e-4_dp], &
      [3, 3]), 3.0e-7_dp, "partially linear NLS covariance")

   do i = 1, size(x_mm, 1)
      x_mm(i, 1) = real(i, dp)
   end do
   y_mm = [1.7_dp, 2.8_dp, 3.6_dp, 4.2_dp, 4.8_dp, 5.1_dp, 5.5_dp, 5.7_dp, 6.0_dp, 6.2_dp]
   fit = nls_fit(x_mm, y_mm, [8.0_dp, 3.0_dp], ss_micmen_model, &
                 jacobian=ss_micmen_jacobian, tolerance=1.0e-10_dp, max_iterations=200)
   call assert_true(fit%converged, "Michaelis-Menten NLS convergence")
   call assert_true(fit%analytic_jacobian, "self-start analytic Jacobian")
   call assert_vector_close(fit%coefficients, &
      [8.87364392883630_dp, 4.35469904036509_dp], 2.0e-6_dp, &
      "Michaelis-Menten NLS coefficients")
   call assert_close(fit%rss, 0.0133394478856795_dp, 2.0e-10_dp, &
                     "Michaelis-Menten NLS RSS")
   call assert_matrix_close(fit%covariance, reshape( &
      [0.00901457932736784_dp, 0.0101638172516955_dp, &
       0.0101638172516955_dp, 0.0122432356589129_dp], [2, 2]), &
      2.0e-7_dp, "Michaelis-Menten NLS covariance")

   fit = nls_fit(x_mm, y_mm, [8.0_dp, 3.0_dp], ss_micmen_model, &
                 jacobian=ss_micmen_jacobian, tolerance=1.0e-10_dp, max_iterations=200, &
                 lower=[0.0_dp, 0.0_dp], upper=[8.5_dp, huge(1.0_dp)])
   call assert_true(fit%converged, "bounded NLS convergence")
   call assert_true(fit%bounded, "bounded NLS metadata")
   call assert_vector_close(fit%coefficients, &
      [8.5_dp, 3.93757084911233_dp], 2.0e-6_dp, "bounded NLS coefficients")
   call assert_close(fit%rss, 0.0416227206694995_dp, 2.0e-9_dp, "bounded NLS RSS")
   call assert_close(fit%residual_standard_error, 0.0721307152583935_dp, 2.0e-9_dp, &
                     "bounded NLS residual standard error")
   call assert_vector_close(fit%fitted_values, &
      [1.72149428529782_dp, 2.86312373056424_dp, 3.67563813827757_dp, &
       4.28342633361216_dp, 4.75520705989380_dp, 5.13203888297869_dp, &
       5.43996476190405_dp, 5.69630127096221_dp, 5.91301109707537_dp, &
       6.09862370711562_dp], 2.0e-6_dp, "bounded NLS fitted values")
   call assert_matrix_close(fit%covariance, reshape( &
      [0.0231006999846804_dp, 0.0256781124699358_dp, &
       0.0256781124699358_dp, 0.0307559514935225_dp], [2, 2]), &
      5.0e-7_dp, "bounded NLS covariance")

   fit = nls_fit(x_mm, y_mm, [1.0_dp, 1.0_dp], redundant_linear_model, &
                 jacobian=redundant_linear_jacobian)
   call assert_true(.not. fit%converged, "rank-deficient NLS convergence flag")
   call assert_true(fit%status /= 0, "rank-deficient NLS status")
   call assert_true(fit%rank == 1, "rank-deficient NLS Jacobian rank")

   call assert_close(ss_asymp(2.0_dp, 10.0_dp, 2.0_dp, log(0.4_dp)), &
                     6.40536828706223_dp, 2.0e-14_dp, "SSasymp value")
   call assert_close(ss_gompertz(2.0_dp, 10.0_dp, 2.0_dp, 0.7_dp), &
                     3.75311098851399_dp, 2.0e-14_dp, "SSgompertz value")
   call assert_close(ss_logis(2.0_dp, 10.0_dp, 3.0_dp, 0.8_dp), &
                     2.22700138825309_dp, 2.0e-14_dp, "SSlogis value")
   call assert_close(ss_micmen(2.0_dp, 10.0_dp, 3.0_dp), 4.0_dp, &
                     2.0e-14_dp, "SSmicmen value")
   call assert_close(ss_weibull(2.0_dp, 10.0_dp, 8.0_dp, log(0.4_dp), 1.5_dp), &
                     7.41927416237091_dp, 2.0e-14_dp, "SSweibull value")

   print *, "test_nls: PASS"

contains

   pure subroutine exponential_model(x, parameters, values)
      !! Evaluates `a*exp(b*x)` for the exponential-growth fixture.
      real(dp), intent(in) :: x(:, :) !! Predictor matrix with one column.
      real(dp), intent(in) :: parameters(:) !! Parameters `[a,b]`.
      real(dp), intent(out) :: values(size(x, 1)) !! Model values.

      values = parameters(1)*exp(parameters(2)*x(:, 1))
   end subroutine exponential_model

   pure subroutine exponential_jacobian(x, parameters, jacobian)
      !! Evaluates the analytic Jacobian of the exponential-growth model.
      real(dp), intent(in) :: x(:, :) !! Predictor matrix with one column.
      real(dp), intent(in) :: parameters(:) !! Parameters `[a,b]`.
      real(dp), intent(out) :: jacobian(size(x, 1), size(parameters)) !! Model Jacobian.

      jacobian(:, 1) = exp(parameters(2)*x(:, 1))
      jacobian(:, 2) = parameters(1)*x(:, 1)*jacobian(:, 1)
   end subroutine exponential_jacobian

   pure subroutine exponential_basis(x, nonlinear_parameters, basis_values)
      !! Evaluates constant and exponential columns for the partially linear fixture.
      real(dp), intent(in) :: x(:, :) !! Predictor matrix with one column.
      real(dp), intent(in) :: nonlinear_parameters(:) !! Decay-rate parameter `[k]`.
      real(dp), intent(out) :: basis_values(:, :) !! Basis columns `[1,exp(-k*x)]`.

      basis_values(:, 1) = 1.0_dp
      basis_values(:, 2) = exp(-nonlinear_parameters(1)*x(:, 1))
   end subroutine exponential_basis

   pure subroutine redundant_linear_model(x, parameters, values)
      !! Evaluates a deliberately unidentified two-parameter linear model.
      real(dp), intent(in) :: x(:, :) !! Predictor matrix with one column.
      real(dp), intent(in) :: parameters(:) !! Two parameters entering only through their sum.
      real(dp), intent(out) :: values(size(x, 1)) !! Model values.

      values = (parameters(1) + parameters(2))*x(:, 1)
   end subroutine redundant_linear_model

   pure subroutine redundant_linear_jacobian(x, parameters, jacobian)
      !! Evaluates the rank-one Jacobian of the redundant linear model.
      real(dp), intent(in) :: x(:, :) !! Predictor matrix with one column.
      real(dp), intent(in) :: parameters(:) !! Unused two-parameter vector.
      real(dp), intent(out) :: jacobian(size(x, 1), size(parameters)) !! Rank-one Jacobian.

      jacobian(:, 1) = x(:, 1)
      jacobian(:, 2) = x(:, 1)
   end subroutine redundant_linear_jacobian



end program test_nls
