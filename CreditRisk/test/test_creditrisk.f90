program test_creditrisk
  use creditrisk
  use creditrisk_calibration
  use creditrisk_reference
  use r_distributions, only: r_pnorm
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_is_nan, ieee_value, ieee_quiet_nan
  implicit none
  real(dp), parameter :: time(4) = [0.5_dp, 1.0_dp, 2.0_dp, 4.0_dp]
  real(dp), parameter :: rate(4) = [0.015_dp, 0.017_dp, 0.022_dp, 0.028_dp]
  real(dp), parameter :: sigma(4) = [0.19_dp, 0.23_dp, 0.21_dp, 0.26_dp]
  type(credit_curve_t) :: curve, other, first, second
  type(calibration_result_t) :: fit
  type(lbfgsb_control_t) :: control
  real(dp) :: actual(3), expected(3), x, textbook, source, dt, hazard, i0, i1, accrued
  real(dp) :: nan, path_time(4), innovations(4,3), sumsq, initial_error
  integer :: i, checks
  checks = 0
  nan = ieee_value(0.0_dp, ieee_quiet_nan)
  call close_vector(cum_normal_density([-8.0_dp, -2.0_dp, 0.0_dp, 0.5_dp, 3.0_dp, 8.0_dp]), &
    ref_normal(:,1), 1.0e-14_dp, 'R Hastings CDF')
  actual = generalized_black_scholes('c', 100.0_dp, 105.0_dp, [0.25_dp, 1.0_dp, 3.0_dp], &
    0.03_dp, 0.01_dp, 0.22_dp)
  call close_vector(actual, ref_call(:,1), 2.0e-12_dp, 'R calls')
  actual = generalized_black_scholes('p', 100.0_dp, 105.0_dp, [0.25_dp, 1.0_dp, 3.0_dp], &
    0.03_dp, 0.01_dp, 0.22_dp)
  call close_vector(actual, ref_put(:,1), 2.0e-12_dp, 'R puts')
  curve = merton(0.55_dp, 1.0_dp, 0.25_dp, 0.02_dp, time)
  call close_curve(curve, ref_merton, 2.0e-12_dp, 'R Merton')
  call close_vector(curve%values(:,3) + curve%values(:,4), curve%values(:,2), 1.0e-14_dp, 'equity plus debt')
  curve = black_cox(0.55_dp, 1.0_dp, 0.25_dp, 0.02_dp, 0.04_dp, time, 0.4_dp)
  call close_curve(curve, ref_black_cox, 1.0e-13_dp, 'R BlackCox')
  curve = at1p(1.0_dp, 0.65_dp, 0.4_dp, sigma, rate, time)
  call close_curve(curve, ref_at1p, 1.0e-13_dp, 'R AT1P')
  curve = sbtv(1.0_dp, [0.4_dp, 0.7_dp], [0.8_dp, 0.2_dp], 0.15_dp, sigma, rate, time)
  call close_curve(curve, ref_sbtv, 1.0e-13_dp, 'R SBTV')
  curve = cds(time, [0.01_dp, 0.03_dp, 0.025_dp, 0.045_dp], rate)
  call close_curve(curve, ref_cds_full, 1.0e-13_dp, 'R CDS continuous')
  curve = cds(time, [0.01_dp, 0.03_dp, 0.025_dp, 0.045_dp], rate, 0.012_dp, 0.25_dp, .true.)
  call close_curve(curve, ref_cds_simple, 1.0e-13_dp, 'R CDS discrete')
  curve = cds2(time, [0.25_dp, 0.75_dp, 1.5_dp, 3.0_dp, 5.0_dp], [0.25_dp, 2.0_dp, 5.0_dp], &
    [0.01_dp, 0.025_dp, 0.02_dp], [0.1_dp, 1.5_dp, 4.5_dp], [0.015_dp, 0.04_dp, 0.05_dp], 0.009_dp, 0.35_dp)
  call close_curve(curve, ref_cds2, 1.0e-13_dp, 'R CDS interpolation')
  path_time = [0.0_dp, 0.25_dp, 1.0_dp, 2.0_dp]
  curve = merton_sim_from_shocks(20.0_dp, 0.05_dp, 0.2_dp, path_time, ref_shocks)
  call close_curve(curve, ref_simulation, 2.0e-13_dp, 'R simulation with paired innovations')
  curve = merton_sim(20.0_dp, 0.05_dp, 0.2_dp, path_time, 3, 1729)
  other = merton_sim(20.0_dp, 0.05_dp, 0.2_dp, path_time, 3, 1729)
  call close_matrix(curve%values, other%values, 0.0_dp, 'seed reproducibility')
  call require(all(curve%values > 0.0_dp), 'simulated prices positive')
  call close_vector(curve%values(1,:), [20.0_dp, 20.0_dp, 20.0_dp], 0.0_dp, 'initial date simulation')

  ! Independent identities and limits, including preserved source formula quirks.
  actual = generalized_black_scholes('c', 100.0_dp, 105.0_dp, [0.25_dp, 1.0_dp, 3.0_dp], &
    0.03_dp, 0.01_dp, 0.22_dp) - generalized_black_scholes('p', 100.0_dp, 105.0_dp, &
    [0.25_dp, 1.0_dp, 3.0_dp], 0.03_dp, 0.01_dp, 0.22_dp)
  expected = 100.0_dp * exp(-0.02_dp * [0.25_dp, 1.0_dp, 3.0_dp]) &
    - 105.0_dp * exp(-0.03_dp * [0.25_dp, 1.0_dp, 3.0_dp])
  call close_vector(actual, expected, 1.0e-13_dp, 'independent put-call parity')
  x = generalized_black_scholes('c', 100.0_dp, 100.0_dp, 0.0_dp, 0.03_dp, 0.01_dp, 0.2_dp)
  call require(x == 0.0_dp, 'ATM expiry extended exact limit')
  x = generalized_black_scholes('p', 100.0_dp, 105.0_dp, 2.0_dp, 0.03_dp, 0.01_dp, 0.0_dp)
  call require(abs(x - max(105.0_dp * exp(-0.06_dp) - 100.0_dp * exp(-0.04_dp), 0.0_dp)) &
    < 1.0e-13_dp, 'zero volatility option limit')
  curve = merton(0.8_dp, 1.0_dp, 0.3_dp, 0.08_dp, [2.0_dp])
  source = r_pnorm((log(1.0_dp / 0.8_dp) + 0.08_dp - 0.5_dp * 0.3_dp ** 2 * 2.0_dp) / (0.3_dp * sqrt(2.0_dp)))
  textbook = r_pnorm((log(1.0_dp / 0.8_dp) + (0.08_dp - 0.5_dp * 0.3_dp ** 2) * 2.0_dp) / (0.3_dp * sqrt(2.0_dp)))
  call require(abs(curve%values(1,5) - source) < 1.0e-14_dp, 'preserved upstream Merton rate term')
  call require(abs(curve%values(1,5) - textbook) > 0.03_dp, 'Merton source differs from textbook explicitly')
  curve = at1p(1.0_dp, 0.5_dp, 0.5_dp, [0.2_dp, 0.2_dp], [0.03_dp, 0.04_dp], [1.0_dp, 2.0_dp])
  expected(1:2) = 2.0_dp * r_pnorm(log(2.0_dp) / sqrt([0.04_dp, 0.08_dp])) - 1.0_dp
  call close_vector(curve%values(:,4), expected(1:2), 1.0e-14_dp, 'AT1P B=half reflection formula')
  curve = at1p(1.0_dp, 0.6_dp, 0.4_dp, [0.0_dp, 0.0_dp], [0.01_dp, 0.02_dp], [1.0_dp, 2.0_dp])
  call close_vector(curve%values(:,4), [1.0_dp, 1.0_dp], 0.0_dp, 'zero variance survival')
  call close_vector(curve%values(:,5), [0.0_dp, 0.0_dp], 0.0_dp, 'zero variance hazards')
  first = at1p(1.0_dp, 0.4_dp, 0.15_dp, sigma, rate, time)
  second = at1p(1.0_dp, 0.7_dp, 0.15_dp, sigma, rate, time)
  curve = sbtv(1.0_dp, [0.4_dp, 0.7_dp], [0.8_dp, 0.2_dp], 0.15_dp, sigma, rate, time)
  call close_vector(curve%values(:,2), 0.8_dp * first%values(:,4) + 0.2_dp * second%values(:,4), &
    1.0e-15_dp, 'scenario survival mixture')
  other = cds(time, curve%values(:,3), rate)
  call close_vector(other%values(:,2), curve%values(:,2), 1.0e-14_dp, 'hazard reconstructs survival')
  curve = cds([1.0_dp, 2.0_dp], [0.0_dp, 0.0_dp], [0.0_dp])
  call close_vector(curve%values(:,3), [0.005_dp, 0.01_dp], 1.0e-16_dp, 'zero hazard and zero rate premium')
  call close_vector(curve%values(:,5), [0.0_dp, 0.0_dp], 0.0_dp, 'zero hazard rates')
  hazard = 0.1_dp
  dt = 1.5_dp
  curve = cds([dt], [hazard], [-hazard])
  call require(abs(curve%values(1,3) - 0.005_dp * (dt - hazard * dt ** 2 / 2.0_dp)) &
    < 1.0e-15_dp, 'cancelling hazard and rate stable premium limit')
  call require(abs(curve%values(1,4) - 0.6_dp * hazard * dt) < 1.0e-15_dp, 'cancelling rate protection limit')
  curve = cds([dt], [hazard], [0.03_dp])
  i0 = (1.0_dp - exp(-0.13_dp * dt)) / 0.13_dp
  i1 = (1.0_dp - (1.0_dp + 0.13_dp * dt) * exp(-0.13_dp * dt)) / 0.13_dp ** 2
  accrued = hazard * i1
  source = 0.005_dp * (dt * exp(-0.13_dp * dt) - accrued)
  textbook = 0.005_dp * (dt * exp(-0.13_dp * dt) + accrued)
  call require(abs(curve%values(1,3) - source) < 1.0e-15_dp, 'preserved source accrued-premium sign')
  call require(abs(curve%values(1,3) - textbook) > 1.0e-4_dp, 'source premium differs from conventional accrued leg')
  call require(abs(curve%values(1,4) - 0.6_dp * hazard * i0) < 1.0e-15_dp, 'analytic continuous protection integral')
  curve = cds([dt], [hazard], [0.03_dp], simplified=.true.)
  call require(abs(curve%values(1,5) - 0.6_dp * (exp(hazard * dt) - 1.0_dp) / dt) &
    < 1.0e-14_dp, 'discrete one-period CDS identity')
  innovations = 0.0_dp
  curve = merton_sim_from_shocks(20.0_dp, 0.05_dp, 0.0_dp, path_time, innovations)
  call close_vector(curve%values(:,1), 20.0_dp * exp(0.05_dp * path_time), 1.0e-13_dp, 'deterministic simulation')
  call close_vector(linear_constant([1.0_dp, 3.0_dp], [2.0_dp, 6.0_dp], [-1.0_dp, 2.0_dp, 4.0_dp]), &
    [2.0_dp, 4.0_dp, 6.0_dp], 0.0_dp, 'interpolation and constant extrapolation')

  ! Invalid domains and shapes must return a status rather than indexing out of bounds.
  curve = cds([1.0_dp, 1.0_dp], [0.1_dp, 0.2_dp], [0.01_dp])
  call require(curve%status == credit_invalid, 'duplicate payment dates rejected')
  curve = at1p(1.0_dp, 0.5_dp, 0.4_dp, sigma(:3), rate, time)
  call require(curve%status == credit_invalid, 'volatility shape rejected')
  curve = black_cox(0.4_dp, 1.0_dp, 0.2_dp, 0.02_dp, 0.01_dp, time, 0.6_dp)
  call require(curve%status == credit_invalid, 'barrier above debt rejected')
  curve = sbtv(1.0_dp, [0.4_dp, 0.7_dp], [0.8_dp, 0.3_dp], 0.15_dp, sigma, rate, time)
  call require(curve%status == credit_invalid, 'bad probability sum rejected')
  curve = cds([1.0_dp], [nan], [0.01_dp])
  call require(curve%status == credit_invalid, 'NaN input rejected')
  x = generalized_black_scholes('x', 100.0_dp, 105.0_dp, 1.0_dp, 0.03_dp, 0.01_dp, 0.2_dp)
  call require(ieee_is_nan(x), 'invalid option flag returns NaN')
  curve = merton_sim(20.0_dp, 0.05_dp, 0.2_dp, path_time, 0)
  call require(curve%status == credit_invalid, 'zero simulation paths rejected')

  ! Stable large-tail/ratio handling and explicit nonfinite-output status.
  call close_vector(cum_normal_density([-1.0e200_dp, 1.0e200_dp]), [0.0_dp, 1.0_dp], 0.0_dp, 'large finite CDF tails')
  x = generalized_black_scholes('c', 1.0e300_dp, 1.0e-300_dp, 1.0_dp, 0.03_dp, 0.03_dp, 0.2_dp)
  call require(ieee_is_finite(x), 'extreme spot-strike ratio finite call')
  call require(abs(x / 1.0e300_dp - 1.0_dp) < 1.0e-14_dp, 'extreme ratio call approaches spot')
  x = generalized_black_scholes('p', 1.0e-300_dp, 1.0e300_dp, 1.0_dp, 0.03_dp, 0.03_dp, 0.2_dp)
  call require(abs(x / 1.0e300_dp - exp(-0.03_dp)) < 1.0e-14_dp, 'extreme ratio put approaches discounted strike')
  x = generalized_black_scholes('c', 1.0e300_dp, 100.0_dp, 1.0_dp, 0.0_dp, 1000.0_dp, 0.2_dp)
  call require(ieee_is_nan(x), 'unrepresentable option returns NaN')
  curve = at1p(1.0_dp, 0.5_dp, 0.4_dp, [0.2_dp], [1000.0_dp], [1.0_dp])
  call require(curve%status == credit_invalid, 'overflowing firm curve has invalid status')
  curve = cds([1.0_dp], [0.1_dp], [-1000.0_dp])
  call require(curve%status == credit_invalid, 'overflowing CDS legs have invalid status')
  curve = cds([real(dp) ::], [real(dp) ::], [0.02_dp])
  call require(curve%status == credit_invalid, 'empty payment grid rejected')
  curve = sbtv(1.0_dp, [0.4_dp], [1.0_dp], 0.15_dp, sigma, rate, time)
  call close_vector(curve%values(:,2), first%values(:,4), 1.0e-15_dp, 'singleton scenario equals AT1P')
  curve = black_cox(0.55_dp, 1.0_dp, 0.25_dp, 0.02_dp, 0.04_dp, [1.0_dp], 0.4_dp)
  call require(curve%status == credit_ok, 'singleton BlackCox grid supported')
  call require(curve%values(1,3) == 0.55_dp, 'singleton BlackCox final debt barrier')

  ! Genuine R optim fixtures, followed by independent objective reconstruction.
  fit = calibrate_black_cox(1.0_dp, ref_quote_black_cox(:,1), 0.02_dp, time)
  call check_fit(fit, ref_calibrate_black_cox_parameters(:,1), ref_calibrate_black_cox_rates(:,1), &
    ref_calibrate_black_cox_error(1,1), ref_quote_black_cox(:,1), 'BlackCox calibration')
  curve = black_cox(fit%parameters(1), 1.0_dp, fit%parameters(3), 0.02_dp, fit%parameters(4), time, fit%parameters(2))
  other = cds(time, curve%values(:,5), [0.02_dp])
  call close_vector(other%values(:,5), fit%rates, 1.0e-13_dp, 'BlackCox returned fitted L/K reconstruct objective')
  call require(abs(fit%parameters(1) - ref_calibrate_black_cox_reported(1,2)) > 1.0e-4_dp, &
    'corrected fitted L differs from upstream initial L')
  fit = calibrate_at1p(1.0_dp, ref_quote_at1p(:,1), rate, time)
  call check_fit(fit, ref_calibrate_at1p_parameters(:,1), ref_calibrate_at1p_rates(:,1), &
    ref_calibrate_at1p_error(1,1), ref_quote_at1p(:,1), 'AT1P calibration')
  curve = at1p(1.0_dp, fit%parameters(1), fit%parameters(2), fit%parameters(3:), rate, time)
  other = cds(time, curve%values(:,5), rate)
  call close_vector(other%values(:,5), fit%rates, 1.0e-13_dp, 'AT1P objective reconstruction')
  fit = calibrate_sbtv(1.0_dp, [0.8_dp, 0.2_dp], ref_quote_sbtv(:,1), rate, time)
  call check_fit(fit, ref_calibrate_sbtv_parameters(:,1), ref_calibrate_sbtv_rates(:,1), &
    ref_calibrate_sbtv_error(1,1), ref_quote_sbtv(:,1), 'SBTV calibration')
  curve = sbtv(1.0_dp, fit%parameters(6:), [0.8_dp, 0.2_dp], fit%parameters(1), fit%parameters(2:5), rate, time)
  other = cds(time, curve%values(:,3), rate)
  call close_vector(other%values(:,5), fit%rates, 1.0e-13_dp, 'SBTV objective reconstruction')
  fit = calibrate_cds(rate, time, [0.75_dp, 1.5_dp, 3.0_dp], [0.008_dp, 0.013_dp, 0.02_dp])
  call check_fit(fit, ref_calibrate_cds_parameters(:,1), ref_calibrate_cds_rates(:,1), &
    ref_calibrate_cds_error(1,1), [0.008_dp, 0.013_dp, 0.02_dp], 'CDS calibration')
  curve = cds(time, linear_constant([0.75_dp, 1.5_dp, 3.0_dp], fit%parameters, time), rate)
  call close_vector(linear_constant(time, curve%values(:,5), [0.75_dp, 1.5_dp, 3.0_dp]), &
    fit%rates, 1.0e-14_dp, 'CDS calibration interpolates rates not legs')
  curve = cds(time, [0.1_dp, 0.1_dp, 0.1_dp, 0.1_dp], rate)
  actual = linear_constant(time, curve%values(:,5), [0.75_dp, 1.5_dp, 3.0_dp])
  initial_error = sum((actual - [0.008_dp, 0.013_dp, 0.02_dp]) ** 2) / 3.0_dp
  call require(fit%error < initial_error * 1.0e-6_dp, 'CDS fit improves initial objective')
  control%max_iterations = 1
  fit = calibrate_cds(rate, time, [0.75_dp, 1.5_dp, 3.0_dp], [0.008_dp, 0.013_dp, 0.02_dp], control=control)
  call require(.not. fit%optimizer%success, 'iteration limit diagnostics exposed')
  call require(fit%optimizer%iterations <= 1, 'iteration cap honored')
  call require(len_trim(fit%warning) > 0, 'nonconvergence warning exposed')
  fit = calibrate_sbtv(1.0_dp, [0.7_dp, 0.3_dp], ref_sbtv_stress_quotes(:,1), &
    [0.011_dp, 0.015_dp, 0.024_dp, 0.035_dp], [0.5_dp, 1.5_dp, 3.0_dp, 7.0_dp])
  call require(fit%status == credit_ok .and. fit%optimizer%success, 'SBTV invalid-trial regression converges')
  call require(fit%error < 1.0e-8_dp, 'SBTV regression obtains a meaningful market fit')
  call require(fit%error <= ref_sbtv_stress_error(1,1) + 1.0e-9_dp, 'SBTV regression R objective quality')
  call close_vector(fit%rates, ref_sbtv_stress_rates(:,1), 5.0e-5_dp, 'SBTV regression R fitted rates')
  call require(fit%max_abs_residual < 0.0001_dp, 'SBTV regression residual below one basis point')
  call require(fit%error < fit%initial_error * 1.0e-5_dp, 'SBTV regression improves initial fit substantially')
  call require(abs(fit%rms_error ** 2 - fit%error) < 1.0e-20_dp, 'RMS error diagnostic consistent')
  call require(len_trim(fit%warning) == 0, 'good calibrated fit has no warning')
  curve = sbtv(1.0_dp, fit%parameters(6:), [0.7_dp, 0.3_dp], fit%parameters(1), fit%parameters(2:5), &
    [0.011_dp, 0.015_dp, 0.024_dp, 0.035_dp], [0.5_dp, 1.5_dp, 3.0_dp, 7.0_dp])
  other = cds([0.5_dp, 1.5_dp, 3.0_dp, 7.0_dp], curve%values(:,3), [0.011_dp, 0.015_dp, 0.024_dp, 0.035_dp])
  call close_vector(other%values(:,5), fit%rates, 1.0e-13_dp, 'SBTV regression objective reconstructed')
  fit = calibrate_cds([0.02_dp], [1.0_dp, 2.0_dp], [1.0_dp, 2.0_dp], [-0.01_dp, -0.01_dp])
  call require(fit%optimizer%success .and. fit%status == credit_ok, 'incompatible quotes can reach optimizer convergence')
  call require(fit%max_abs_residual >= 0.01_dp, 'incompatible quotes expose actual residual')
  call require(len_trim(fit%warning) > 0, 'poor market fit warns despite optimizer convergence')
  fit = calibrate_cds([0.02_dp], [1.0_dp, 2.0_dp], [1.0_dp, 2.0_dp], [1.0e300_dp, 1.0e300_dp])
  call require(fit%status == credit_invalid, 'unrepresentable objective rejected before optimization')
  call require(fit%optimizer%function_evaluations == 0, 'huge-quote objective does not enter optimizer')
  print '(a,i0,a)', 'CreditRisk: ', checks, ' checks passed.'  
contains
  subroutine require(condition, label)
    logical, intent(in) :: condition !! Assertion expected to hold for this test case.
    character(len=*), intent(in) :: label !! Human-readable description printed on failure.
    checks = checks + 1
    if (.not. condition) then
      print '(a)', 'FAIL: ' // label
      error stop 1
    end if
  end subroutine require
  subroutine close_vector(observed, target, tolerance, label)
    real(dp), intent(in) :: observed(:) !! Computed vector to compare with a trusted reference.
    real(dp), intent(in) :: target(:) !! Reference vector of equal length.
    real(dp), intent(in) :: tolerance !! Nonnegative maximum absolute error allowed.
    character(len=*), intent(in) :: label !! Description identifying the numerical assertion.
    call require(size(observed) == size(target), label // ' shape')
    call require(all(ieee_is_finite(observed)), label // ' finite')
    if (maxval(abs(observed - target)) > tolerance) then
      print *, 'observed ', observed
      print *, 'target   ', target
      print *, 'tolerance ', tolerance
    end if
    call require(maxval(abs(observed - target)) <= tolerance, label)
  end subroutine close_vector
  subroutine close_matrix(observed, target, tolerance, label)
    real(dp), intent(in) :: observed(:,:) !! Computed table to compare with a trusted reference.
    real(dp), intent(in) :: target(:,:) !! Reference table with the same row and column counts.
    real(dp), intent(in) :: tolerance !! Maximum absolute elementwise error allowed.
    character(len=*), intent(in) :: label !! Description identifying the numerical assertion.
    call require(all(shape(observed) == shape(target)), label // ' shape')
    do i = 1, size(target,2)
      call close_vector(observed(:,i), target(:,i), tolerance, label)
    end do
  end subroutine close_matrix
  subroutine close_curve(result, target, tolerance, label)
    type(credit_curve_t), intent(in) :: result !! Model output with an expected successful status.
    real(dp), intent(in) :: target(:,:) !! R numeric reference table.
    real(dp), intent(in) :: tolerance !! Maximum absolute elementwise error allowed.
    character(len=*), intent(in) :: label !! Description identifying the curve test.
    call require(result%status == credit_ok, label // ' status')
    call close_matrix(result%values, target, tolerance, label)
  end subroutine close_curve
  subroutine check_fit(result, parameters, rates, objective, quotes, label)
    type(calibration_result_t), intent(in) :: result !! Computed calibration including convergence diagnostics.
    real(dp), intent(in) :: parameters(:) !! Optimized parameter vector recorded from actual R optim.
    real(dp), intent(in) :: rates(:) !! Fitted CDS rates recorded from the R optimizer output.
    real(dp), intent(in) :: objective !! Mean squared error reported by R optim.
    real(dp), intent(in) :: quotes(:) !! Input market quotes used for independent error reconstruction.
    character(len=*), intent(in) :: label !! Calibrator name used in test diagnostics.
    call require(result%status == credit_ok, label // ' valid output')
    call require(result%optimizer%success, label // ' optimizer converged')
    call close_vector(result%parameters, parameters, 1.0e-5_dp, label // ' R parameters')
    call close_vector(result%rates, rates, 2.0e-7_dp, label // ' R rates')
    call require(abs(result%error - objective) < 1.0e-10_dp, label // ' R objective')
    sumsq = sum((result%rates - quotes) ** 2) / real(size(quotes), dp)
    call require(abs(result%error - sumsq) < 1.0e-20_dp, label // ' reconstructed mean-square error')
    call require(result%optimizer%function_evaluations > 1, label // ' actual optimization performed')
  end subroutine check_fit
end program test_creditrisk
