program test_burstfin
   use, intrinsic :: ieee_arithmetic, only: ieee_is_nan, ieee_value, ieee_quiet_nan, ieee_positive_inf
   use burstfin
   implicit none
   include 'r_reference.inc'
   type(statistical_factor_model) :: model, scaled_model
   type(alpha_surface) :: surface
   real(dp), allocatable :: cov(:, :), systematic(:, :), specific(:, :), result(:, :), weights(:)
   real(dp), allocatable :: cube(:, :, :), output_cube(:, :, :), relative(:, :)
   real(dp) :: missing(15, 6), scaled(9, 4), intensity, reference_intensity, nan, inf
   real(dp) :: rank_data(6, 3), constant_data(9, 5), w15(15)
   integer :: status, i, j, tests
   real(dp), parameter :: heterogeneous_scales(4) = [1.0e90_dp, 1.0e-90_dp, 1.0e40_dp, 1.0e-40_dp]
   tests = 0
   nan = ieee_value(0.0_dp, ieee_quiet_nan)
   inf = ieee_value(0.0_dp, ieee_positive_inf)

   call factor_model_stat(x_reference, model, status)
   call check(status == 0, 'factor default status')
   call factor_covariance(model, cov, status)
   call near_matrix(cov, factor_default_cov, 3.0e-12_dp, 'R factor default covariance')
   call factor_covariance(model, systematic, status, 'systematic')
   call near_matrix(systematic, factor_default_systematic, 3.0e-12_dp, 'R systematic output')
   call factor_covariance(model, specific, status, 'specific')
   call near_matrix(specific, factor_default_specific, 3.0e-12_dp, 'R specific output')
   call near_matrix(cov, systematic + specific, 1.0e-14_dp, 'covariance components sum')
   call near_vector(model%uniquenesses, factor_default_unique, 3.0e-12_dp, 'R uniquenesses')
   call near_vector(model%sdev, factor_default_sd, 3.0e-12_dp, 'R standard deviations')
   call near_vector(model%cumulative_variance_fraction, factor_default_cumulative, 3.0e-12_dp, 'R cumulative variance')
   call factor_model_stat(x_reference, model, status, weights = [1.0_dp], center = .false., frac_var = 0.8_dp, &
      iter_max = 0, range_factors = [1, 2], specific_floor = 0.25_dp, floor_type = 'fraction')
   call check(status == 0, 'factor explicit options status')
   call factor_covariance(model, cov, status)
   call near_matrix(cov, factor_noiter_cov, 3.0e-12_dp, 'R factor no iteration and fraction floor')

   call var_shrink_eqcor(x_reference, cov, intensity, status)
   call check(status == 0, 'shrink default status')
   call near_matrix(cov, shrink_default_cov, 5.0e-12_dp, 'R shrink default covariance')
   call check(abs(intensity - shrink_default_intensity(1)) < 5.0e-12_dp, 'R automatic intensity')
   call var_shrink_eqcor(x_reference, cov, intensity, status, weights = [0.0_dp, 1.0_dp, 2.0_dp, &
      3.0_dp, 0.0_dp, 2.0_dp, 1.0_dp, 3.0_dp, 2.0_dp], shrink = 0.37_dp, &
      center_values = [0.002_dp, -0.003_dp, 0.001_dp, 0.005_dp], vol_shrink = 0.65_dp, &
      compatible = .true., tol = 0.03_dp)
   call check(status == 0, 'shrink explicit options status')
   call near_matrix(cov, shrink_explicit_cov, 5.0e-12_dp, 'R all explicit shrink options')
   call check(intensity == 0.37_dp, 'explicit intensity retained')

   missing = missing_base
   missing(15, :) = nan
   missing(2, 4) = nan
   missing([1, 3], 5) = nan
   missing(3:14, 6) = nan
   do i = 1, 15
      w15(i) = 0.2_dp + real(i - 1, dp) * 0.1_dp
   end do
   call factor_model_stat(missing, model, status, weights = w15, full_min = 2, reg_min = 3, &
      sd_min = 4, nfac_miss = [0, 1, 2], iter_max = 3, range_factors = [2, 2])
   call check(status == 0, 'factor missing status')
   call factor_covariance(model, cov, status)
   call near_matrix(cov, factor_missing_cov, 2.0e-11_dp, 'R missing factor covariance')
   call near_vector(model%sdev, factor_missing_sd, 5.0e-12_dp, 'R missing volatility replacement')
   call near_vector(model%uniquenesses, factor_missing_unique, 2.0e-11_dp, 'R missing regression uniqueness')
   call near_vector(model%cumulative_variance_fraction, factor_missing_cumulative, 5.0e-12_dp, 'R missing cumulative')
   call var_shrink_eqcor(missing, cov, intensity, status, weights = [1.0_dp], sd_min = 4, &
      vol_shrink = 0.3_dp, tol = 0.005_dp)
   call check(status == 0, 'shrink missing status')
   call near_matrix(cov, shrink_missing_cov, 2.0e-11_dp, 'R missing automatic shrink covariance')
   call check(abs(intensity - shrink_missing_intensity(1)) < 2.0e-11_dp, 'R missing automatic intensity')
   call var_shrink_eqcor(missing, cov, intensity, status, weights = [1.0_dp], sd_min = 4, &
      shrink = 0.2_dp, center = .false., tol = 0.0_dp)
   call near_matrix(cov, shrink_missing_fixed_cov, 2.0e-11_dp, 'R missing explicit shrink covariance')

   missing(4, 5) = nan
   call factor_model_stat(missing, model, status, weights = [1.0_dp], center_values = spread(0.001_dp, 1, 6), &
      zero_load = .true., full_min = 2, reg_min = 3, sd_min = 4, quan_sd = 0.6_dp, &
      nfac_miss = [0, 1, 2], iter_max = 8, tol = 1.0e-9_dp, range_factors = [2, 2])
   call check(status == 0, 'multiple missing-score regressions status')
   call factor_covariance(model, cov, status)
   call near_matrix(cov, factor_multireg_cov, 2.0e-11_dp, 'R multivariate scores, explicit centers and zero fallback')
   call near_vector(model%uniquenesses, factor_multireg_unique, 2.0e-11_dp, 'R multiple-score regression uniqueness')
   call near_vector(model%sdev, factor_multireg_sd, 5.0e-12_dp, 'R explicit volatility replacement quantile')

   call alpha_proxy_grid([0.8_dp, 0.2_dp, 0.5_dp], [0.24_dp], [0.15_dp], [0.18_dp], &
      [0.5_dp, -0.3_dp], [0.2_dp], surface, status, transpose_output = .true.)
   call check(status == 0, 'alpha grid status')
   call near_matrix(surface%z, alpha_grid_expected, 1.0e-13_dp, 'R sorted transposed alpha grid')
   call near_vector(surface%x, [-0.3_dp, 0.5_dp], 1.0e-14_dp, 'alpha first sorted axis')
   call near_vector(surface%y, [0.2_dp, 0.5_dp, 0.8_dp], 1.0e-14_dp, 'alpha second sorted axis')
   call check(alpha_proxy() == 0.0_dp, 'identical manager and benchmark zero alpha')
   call check(abs(alpha_proxy(0.5_dp, 0.3_dp, 0.2_dp, 0.1_dp, 0.0_dp, 0.0_dp) + 250.0_dp) < &
      1.0e-12_dp, 'independent alpha analytic value')
   call check(ieee_is_nan(alpha_proxy(weight = 0.0_dp)), 'alpha invalid boundary NaN')
   call check(ieee_is_nan(alpha_proxy(cor_man = nan)), 'alpha missing NaN')

   call slide_weight(8, weights, status, fractions = [0.1875_dp, 0.6875_dp])
   call near_vector(weights, slide_fraction_expected, 1.0e-14_dp, 'R weights half-even fractions')
   call slide_weight(8, weights, status, observations = [9.0_dp, 2.0_dp])
   call near_vector(weights, slide_observations_expected, 1.0e-14_dp, 'R weights ages outside interval')
   call slide_weight(4, weights, status)
   call near_vector(weights, [0.2_dp, 0.4_dp, 0.6_dp, 0.8_dp], 1.0e-14_dp, 'analytic default weights')
   call slide_weight(4, weights, status, locations = [1.0_dp, 1.0_dp])
   call near_vector(weights, [0.0_dp, 1.0_dp, 1.0_dp, 1.0_dp], 1.0e-14_dp, 'equal endpoint ramp')
   call slide_weight(4, weights, status, locations = [4.0_dp, 5.0_dp])
   call check(status /= 0, 'all-zero ramp rejected')
   call slide_weight(0, weights, status)
   call check(status /= 0, 'empty sample rejected')
   call slide_weight(4, weights, status, observations = [0.0_dp, 1.0_dp], locations = [0.0_dp, 1.0_dp])
   call check(status /= 0, 'conflicting ramp coordinates rejected')

   call var_add_benchmark(benchmark_input, [-2.0_dp, 3.0_dp], result, status, asset_indices = [3, 1])
   call check(status == 0, 'benchmark sparse reordered weights status')
   call near_matrix(result, benchmark_added_expected, 1.0e-13_dp, 'R appended benchmark')
   call var_relative_benchmark(result, 4, relative, status)
   call near_matrix(relative, benchmark_relative_expected, 1.0e-13_dp, 'R relative benchmark')
   call var_add_benchmark(benchmark_input, [1.0_dp, 0.0_dp, 0.0_dp], result, status)
   call check(result(4, 4) == benchmark_input(1, 1), 'unit benchmark variance')
   call near_vector(result(:3, 4), benchmark_input(:, 1), 1.0e-14_dp, 'unit benchmark cross covariance')
   call var_add_benchmark(benchmark_input, [0.0_dp, 0.0_dp, 0.0_dp], result, status)
   call check(status /= 0, 'zero normalization rejected')
   call var_add_benchmark(benchmark_input, [0.0_dp, 0.0_dp, 0.0_dp], result, status, sum_to_one = .false.)
   call check(status == 0 .and. result(4, 4) == 0.0_dp, 'zero unnormalized benchmark')
   allocate(cube(3, 3, 2))
   cube(:, :, 1) = benchmark_input
   cube(:, :, 2) = 2.0_dp * benchmark_input
   call var_add_benchmark(cube, [-2.0_dp, 3.0_dp], output_cube, status, asset_indices = [3, 1])
   call near_matrix(output_cube(:, :, 2), 2.0_dp * benchmark_added_expected, 1.0e-13_dp, 'benchmark cube')
   cube = output_cube
   call var_relative_benchmark(cube, 4, output_cube, status)
   call near_matrix(output_cube(:, :, 2), 2.0_dp * benchmark_relative_expected, 1.0e-13_dp, 'relative cube')
   call var_relative_benchmark(reshape([2.0_dp], [1, 1]), 1, result, status)
   call check(status == 0 .and. size(result) == 0, 'single benchmark empty covariance')

   ! Independent rank-one covariance: population variance of 1:6 is 35/12.
   rank_data(:, 1) = [(real(i, dp), i = 1, 6)]
   rank_data(:, 2) = 2.0_dp * rank_data(:, 1)
   rank_data(:, 3) = -rank_data(:, 1)
   call factor_model_stat(rank_data, model, status, weights = [1.0_dp], specific_floor = 0.0_dp)
   call check(status == 0, 'rank-deficient factor accepted')
   call factor_covariance(model, cov, status)
   result = reshape([1.0_dp, 2.0_dp, -1.0_dp, 2.0_dp, 4.0_dp, -2.0_dp, &
      -1.0_dp, -2.0_dp, 1.0_dp], [3, 3]) * (35.0_dp / 12.0_dp)
   call near_matrix(cov, result, 5.0e-12_dp, 'rank-one analytic covariance')
   call factor_model_stat(x_reference, model, status, range_factors = [0, 0], specific_floor = 0.0_dp)
   call check(status == 0 .and. size(model%loadings, 2) == 0, 'zero-factor extension')
   call factor_covariance(model, cov, status, 'systematic')
   call check(maxval(abs(cov)) == 0.0_dp, 'zero-factor systematic covariance')
   constant_data(:, :4) = x_reference
   constant_data(:, 5) = 0.0_dp
   call factor_model_stat(constant_data, model, status)
   call check(status == 3, 'constant asset default rejection')
   call factor_model_stat(constant_data, model, status, constant_returns_okay = .true.)
   call check(status == 0 .and. model%constant_assets(5), 'constant asset permitted flag')
   call factor_covariance(model, cov, status)
   call check(maxval(abs(cov(5, :))) == 0.0_dp, 'constant asset zero covariance')

   do j = 1, 4
      scaled(:, j) = x_reference(:, j) * heterogeneous_scales(j)
   end do
   call var_shrink_eqcor(scaled, cov, intensity, status, shrink = 0.5_dp, tol = 0.0_dp, sd_min = 2)
   call check(status == 0, 'explicit shrink heterogeneous dynamic range accepted')
   do j = 1, 4
      do i = 1, 4
         cov(i, j) = (cov(i, j) / heterogeneous_scales(i)) / heterogeneous_scales(j)
      end do
   end do
   call near_matrix(cov, shrink_heterogeneous_scaled, 1.0e-11_dp, 'R heterogeneous covariance every entry rescaled')

   ! Independent homogeneity and weight-scale invariance, over 280 orders of magnitude.
   call var_shrink_eqcor(x_reference, result, reference_intensity, status)
   scaled = x_reference * 1.0e140_dp
   call var_shrink_eqcor(scaled, cov, intensity, status)
   call check(status == 0, 'large-scale shrink status')
   call near_matrix(cov / 1.0e280_dp, result, 5.0e-12_dp, 'large-scale shrink homogeneity')
   call check(abs(intensity - reference_intensity) < 5.0e-12_dp, 'scale invariant shrink intensity')
   scaled = x_reference * 1.0e-140_dp
   call var_shrink_eqcor(scaled, cov, intensity, status)
   call check(status == 0, 'small-scale shrink status')
   call near_matrix(cov / 1.0e-280_dp, result, 5.0e-12_dp, 'small-scale shrink homogeneity')
   call factor_model_stat(x_reference, model, status)
   call factor_covariance(model, result, status)
   call factor_model_stat(x_reference * 1.0e140_dp, scaled_model, status)
   call check(status == 0, 'large-scale factor status')
   call factor_covariance(scaled_model, cov, status)
   call near_matrix(cov / 1.0e280_dp, result, 5.0e-12_dp, 'large-scale factor homogeneity')
   call factor_model_stat(x_reference * 1.0e-140_dp, scaled_model, status)
   call factor_covariance(scaled_model, cov, status)
   call near_matrix(cov / 1.0e-280_dp, result, 5.0e-12_dp, 'small-scale factor homogeneity')
   call factor_model_stat(x_reference, model, status, weights = [(real(i, dp), i = 1, 9)])
   call factor_model_stat(x_reference, scaled_model, status, weights = [(real(i, dp) * 1.0e300_dp, i = 1, 9)])
   call near_vector(model%sdev, scaled_model%sdev, 5.0e-12_dp, 'huge weight normalization')

   call factor_model_stat(missing, model, status)
   call check(status == 2, 'insufficient complete columns')
   call factor_model_stat(x_reference, model, status, weights = [0.0_dp])
   call check(status == 1, 'invalid all-zero factor weights')
   call factor_model_stat(x_reference, model, status, frac_var = 1.1_dp)
   call check(status == 1, 'invalid factor fraction')
   call var_shrink_eqcor(x_reference, cov, intensity, status, weights = [-1.0_dp])
   call check(status == 1, 'negative shrink weights rejected')
   scaled = x_reference
   scaled(1, 1) = inf
   call var_shrink_eqcor(scaled, cov, intensity, status)
   call check(status == 1, 'infinite shrink input rejected')
   call factor_model_stat(scaled, model, status, full_min = 2, sd_min = 2, reg_min = 3)
   call check(status == 0, 'factor infinite input treated as missing')
   call var_shrink_eqcor(x_reference * 0.0_dp, cov, intensity, status)
   call check(status == 3, 'all-zero shrink data rejected')
   call var_shrink_eqcor(x_reference, cov, intensity, status, shrink = 2.0_dp)
   call check(status == 0 .and. intensity == 1.0_dp, 'explicit shrink upper clipping')
   print '(a,i0,a)', 'BurStFin: ', tests, ' deterministic checks passed.'
contains
   subroutine check(condition, message)
      logical, intent(in) :: condition !! Assertion that must hold for the current numerical test.
      character(*), intent(in) :: message !! Description printed if the assertion fails.
      if (.not. condition) then
         print '(a)', 'FAILED: ' // message
         error stop 1
      end if
      tests = tests + 1
   end subroutine check

   subroutine near_matrix(actual, expected, tolerance, message)
      real(dp), intent(in) :: actual(:, :) !! Computed finite matrix to compare with its reference.
      real(dp), intent(in) :: expected(:, :) !! Oracle or analytically expected matrix of matching shape.
      real(dp), intent(in) :: tolerance !! Relative maximum-norm error tolerance with an absolute tiny floor.
      character(*), intent(in) :: message !! Context for a failed matrix comparison.
      real(dp) :: error, scale
      call check(all(shape(actual) == shape(expected)), message // ' shape')
      scale = max(maxval(abs(expected)), tiny(1.0_dp))
      error = maxval(abs(actual - expected)) / scale
      if (error > tolerance) print '(a,es16.7)', 'relative matrix error: ', error
      call check(error <= tolerance, message)
   end subroutine near_matrix

   subroutine near_vector(actual, expected, tolerance, message)
      real(dp), intent(in) :: actual(:) !! Computed finite vector to compare with its reference.
      real(dp), intent(in) :: expected(:) !! Oracle or analytic vector of matching length.
      real(dp), intent(in) :: tolerance !! Relative maximum-norm tolerance.
      character(*), intent(in) :: message !! Context for a failed vector comparison.
      real(dp) :: error, scale
      call check(size(actual) == size(expected), message // ' size')
      scale = max(maxval(abs(expected)), tiny(1.0_dp))
      error = maxval(abs(actual - expected)) / scale
      if (error > tolerance) print '(a,es16.7)', 'relative vector error: ', error
      call check(error <= tolerance, message)
   end subroutine near_vector
end program test_burstfin
