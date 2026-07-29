! SPDX-License-Identifier: GPL-3.0-only
! SPDX-FileComment: Diagnostic tests translated from the R rugarch package.
module rugarch_diagnostics_mod
   !! Stability, distributional, duration, moment, and PIT diagnostics.
   use kind_mod, only: dp
   use distribution_mod, only: standardized_cdf
   use linalg_mod, only: symmetric_pseudoinverse
   use special_functions_mod, only: regularized_gamma_q
   use stats_mod, only: variance
   use rugarch_mod, only: rugarch_fit_t
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use, intrinsic :: iso_fortran_env, only: output_unit
   implicit none
   private

   type, public :: rugarch_nyblom_test_t
      !! Hansen-Nyblom individual and joint parameter-stability statistics.
      real(dp), allocatable :: individual_statistic(:)
      real(dp) :: joint_statistic = 0.0_dp
      real(dp) :: individual_critical(3) = [0.353_dp, 0.470_dp, 0.748_dp]
      real(dp) :: joint_critical(3) = 0.0_dp
      integer :: info = 0
   end type rugarch_nyblom_test_t

   type, public :: rugarch_gof_test_t
      !! Pearson probability-integral-transform goodness-of-fit tests.
      integer, allocatable :: groups(:)
      real(dp), allocatable :: statistic(:)
      real(dp), allocatable :: p_value(:)
      integer :: info = 0
   end type rugarch_gof_test_t

   type, public :: rugarch_var_duration_test_t
      !! Weibull duration test for independence of value-at-risk exceedances.
      real(dp) :: shape = 1.0_dp
      real(dp) :: unrestricted_log_likelihood = -huge(1.0_dp)
      real(dp) :: restricted_log_likelihood = -huge(1.0_dp)
      real(dp) :: likelihood_ratio = 0.0_dp
      real(dp) :: p_value = 1.0_dp
      integer :: exceedances = 0
      integer :: info = 0
   end type rugarch_var_duration_test_t

   type, public :: rugarch_gmm_test_t
      !! Orthogonality and serial moment tests for standardized innovations.
      real(dp) :: moment_mean(4) = 0.0_dp
      real(dp) :: moment_variance(4) = 0.0_dp
      real(dp) :: moment_t(4) = 0.0_dp
      real(dp) :: serial_statistic(3) = 0.0_dp
      real(dp) :: serial_p_value(3) = 1.0_dp
      real(dp) :: joint_statistic = 0.0_dp
      real(dp) :: joint_p_value = 1.0_dp
      integer :: lags = 0
      integer :: info = 0
   end type rugarch_gmm_test_t

   type, public :: rugarch_hong_li_test_t
      !! Hong-Li moment and nonparametric PIT specification statistics.
      real(dp) :: statistic(7) = 0.0_dp
      real(dp) :: p_value(7) = 1.0_dp
      integer :: lags = 0
      integer :: info = 0
   end type rugarch_hong_li_test_t

   interface display
      module procedure display_rugarch_nyblom_test
      module procedure display_rugarch_gof_test
      module procedure display_rugarch_var_duration_test
      module procedure display_rugarch_gmm_test
      module procedure display_rugarch_hong_li_test
   end interface display

   interface rugarch_nyblom_test
      module procedure rugarch_nyblom_test_scores
      module procedure rugarch_nyblom_test_fit
   end interface rugarch_nyblom_test

   public :: rugarch_nyblom_test, rugarch_gof_test
   public :: rugarch_var_duration_test, rugarch_gmm_test
   public :: rugarch_hong_li_test, display

contains

   pure function rugarch_nyblom_test_scores(scores) result(out)
      !! Test parameter constancy from per-observation likelihood scores.
      real(dp), intent(in) :: scores(:, :) !! Score by observation and parameter.
      type(rugarch_nyblom_test_t) :: out
      real(dp), allocatable :: information(:, :), inverse(:, :)
      real(dp), allocatable :: cumulative(:, :), cross_product(:, :)
      integer :: observations, parameters, time, inverse_info

      observations = size(scores, 1)
      parameters = size(scores, 2)
      if (observations < 2 .or. parameters < 1 .or. parameters > 20 .or. &
         .not. all(ieee_is_finite(scores))) then
         out%info = 1
         return
      end if
      allocate(information(parameters, parameters), &
         inverse(parameters, parameters), cumulative(observations, parameters), &
         cross_product(parameters, parameters), &
         out%individual_statistic(parameters))
      information = matmul(transpose(scores), scores)
      call symmetric_pseudoinverse(information, inverse, inverse_info)
      if (inverse_info /= 0) then
         out%info = 2
         return
      end if
      cumulative(1, :) = scores(1, :)
      do time = 2, observations
         cumulative(time, :) = cumulative(time - 1, :) + scores(time, :)
      end do
      cross_product = matmul(transpose(cumulative), cumulative)
      out%joint_statistic = sum_diagonal(matmul(cross_product, inverse))/ &
         real(observations, dp)
      out%individual_statistic = diagonal(cross_product)/ &
         max(diagonal(information)*real(observations, dp), tiny(1.0_dp))
      out%joint_critical = nyblom_critical(parameters)
   end function rugarch_nyblom_test_scores

   pure function rugarch_nyblom_test_fit(fit) result(out)
      !! Test parameter constancy directly from a fitted rugarch model.
      type(rugarch_fit_t), intent(in) :: fit !! Fitted model containing likelihood scores.
      type(rugarch_nyblom_test_t) :: out

      if (fit%info /= 0 .or. .not. allocated(fit%scores)) then
         out%info = 1
         return
      end if
      out = rugarch_nyblom_test_scores(fit%scores)
   end function rugarch_nyblom_test_fit

   pure function rugarch_gof_test(standardized_residuals, distribution, &
      groups, shape, skew, lambda) result(out)
      !! Apply grouped Pearson tests to innovation probability transforms.
      real(dp), intent(in) :: standardized_residuals(:) !! Standardized innovations.
      integer, intent(in) :: distribution !! Innovation distribution code.
      integer, intent(in) :: groups(:) !! Requested numbers of equiprobable cells.
      real(dp), intent(in), optional :: shape !! Distribution shape parameter.
      real(dp), intent(in), optional :: skew !! Distribution skew parameter.
      real(dp), intent(in), optional :: lambda !! Generalized-hyperbolic lambda.
      type(rugarch_gof_test_t) :: out
      real(dp), allocatable :: probability(:)
      integer, allocatable :: counts(:)
      real(dp) :: selected_shape, selected_skew, selected_lambda, expected
      integer :: test, observation, cell

      if (size(standardized_residuals) < 2 .or. size(groups) < 1 .or. &
         any(groups < 2) .or. .not. all(ieee_is_finite( &
         standardized_residuals))) then
         out%info = 1
         return
      end if
      selected_shape = 5.0_dp
      selected_skew = 1.0_dp
      selected_lambda = -0.5_dp
      if (present(shape)) selected_shape = shape
      if (present(skew)) selected_skew = skew
      if (present(lambda)) selected_lambda = lambda
      probability = standardized_cdf(standardized_residuals, distribution, &
         selected_shape, selected_skew, selected_lambda)
      allocate(out%groups(size(groups)), out%statistic(size(groups)), &
         out%p_value(size(groups)))
      out%groups = groups
      do test = 1, size(groups)
         allocate(counts(groups(test)))
         counts = 0
         do observation = 1, size(probability)
            cell = min(groups(test), max(1, 1 + int(probability(observation)* &
               real(groups(test), dp))))
            counts(cell) = counts(cell) + 1
         end do
         expected = real(size(probability), dp)/real(groups(test), dp)
         out%statistic(test) = sum((real(counts, dp) - expected)**2/expected)
         out%p_value(test) = regularized_gamma_q(0.5_dp* &
            real(groups(test) - 1, dp), 0.5_dp*out%statistic(test))
         deallocate(counts)
      end do
   end function rugarch_gof_test

   pure function rugarch_var_duration_test(alpha, actual, value_at_risk) &
      result(out)
      !! Test whether VaR exceedance durations have exponential memorylessness.
      real(dp), intent(in) :: alpha !! Nominal lower-tail probability.
      real(dp), intent(in) :: actual(:) !! Realized observations.
      real(dp), intent(in) :: value_at_risk(:) !! Matching VaR forecasts.
      type(rugarch_var_duration_test_t) :: out
      integer, allocatable :: hit_index(:)
      real(dp), allocatable :: duration(:)
      logical, allocatable :: censored(:), hit(:)
      real(dp) :: lower, upper, first, second, first_value, second_value
      integer :: observations, hits, duration_count, index, iteration

      observations = size(actual)
      if (alpha <= 0.0_dp .or. alpha >= 1.0_dp .or. observations < 3 .or. &
         size(value_at_risk) /= observations .or. &
         .not. all(ieee_is_finite(actual)) .or. &
         .not. all(ieee_is_finite(value_at_risk))) then
         out%info = 1
         return
      end if
      hit = actual < value_at_risk
      hits = count(hit)
      out%exceedances = hits
      if (hits < 2) then
         out%info = 2
         return
      end if
      hit_index = pack([(index, index=1, observations)], hit)
      duration_count = hits - 1
      if (.not. hit(1)) duration_count = duration_count + 1
      if (.not. hit(observations)) duration_count = duration_count + 1
      allocate(duration(duration_count), censored(duration_count))
      index = 0
      if (.not. hit(1)) then
         index = index + 1
         duration(index) = real(hit_index(1), dp)
         censored(index) = .true.
      end if
      duration(index + 1:index + hits - 1) = real(hit_index(2:) - &
         hit_index(:hits - 1), dp)
      censored(index + 1:index + hits - 1) = .false.
      index = index + hits - 1
      if (.not. hit(observations)) then
         duration(duration_count) = real(observations - hit_index(hits), dp)
         censored(duration_count) = .true.
      end if
      lower = 0.001_dp
      upper = 10.0_dp
      do iteration = 1, 100
         first = lower + 0.3819660112501051_dp*(upper - lower)
         second = upper - 0.3819660112501051_dp*(upper - lower)
         first_value = duration_negative_log_likelihood(first, duration, &
            censored)
         second_value = duration_negative_log_likelihood(second, duration, &
            censored)
         if (first_value < second_value) then
            upper = second
         else
            lower = first
         end if
      end do
      out%shape = 0.5_dp*(lower + upper)
      out%unrestricted_log_likelihood = -duration_negative_log_likelihood( &
         out%shape, duration, censored)
      out%restricted_log_likelihood = -duration_negative_log_likelihood( &
         1.0_dp, duration, censored)
      out%likelihood_ratio = max(0.0_dp, 2.0_dp*( &
         out%unrestricted_log_likelihood - out%restricted_log_likelihood))
      out%p_value = regularized_gamma_q(0.5_dp, &
         0.5_dp*out%likelihood_ratio)
   end function rugarch_var_duration_test

   pure function rugarch_gmm_test(standardized_residuals, lags, skewness, &
      kurtosis) result(out)
      !! Test unconditional and serial orthogonality moment restrictions.
      real(dp), intent(in) :: standardized_residuals(:) !! Standardized innovations.
      integer, intent(in), optional :: lags !! Positive serial-moment lag count.
      real(dp), intent(in), optional :: skewness !! Null third raw moment.
      real(dp), intent(in), optional :: kurtosis !! Null fourth raw moment.
      type(rugarch_gmm_test_t) :: out
      real(dp), allocatable :: moments(:, :), covariance(:, :), inverse(:, :)
      real(dp), allocatable :: average(:), serial(:, :), serial_average(:)
      real(dp) :: selected_skewness, selected_kurtosis
      integer :: selected_lags, observations, effective, condition, power
      integer :: lag, time, inverse_info, offset

      selected_lags = 1
      if (present(lags)) selected_lags = lags
      selected_skewness = 0.0_dp
      if (present(skewness)) selected_skewness = skewness
      selected_kurtosis = 3.0_dp
      if (present(kurtosis)) selected_kurtosis = kurtosis
      observations = size(standardized_residuals)
      effective = observations - selected_lags
      if (selected_lags < 1 .or. effective < 8 .or. &
         .not. all(ieee_is_finite(standardized_residuals))) then
         out%info = 1
         return
      end if
      out%lags = selected_lags
      allocate(moments(4 + 3*selected_lags, effective), &
         average(4 + 3*selected_lags))
      do time = 1, effective
         moments(1, time) = standardized_residuals(time + selected_lags)
         moments(2, time) = moments(1, time)**2 - 1.0_dp
         moments(3, time) = moments(1, time)**3 - selected_skewness
         moments(4, time) = moments(1, time)**4 - selected_kurtosis
      end do
      do power = 2, 4
         offset = 4 + (power - 2)*selected_lags
         do lag = 1, selected_lags
            do time = 1, effective
               moments(offset + lag, time) = &
                  (standardized_residuals(time + selected_lags)**power - &
                  merge(1.0_dp, merge(selected_skewness, selected_kurtosis, &
                  power == 3), power == 2))* &
                  (standardized_residuals(time + selected_lags - lag)**power - &
                  merge(1.0_dp, merge(selected_skewness, selected_kurtosis, &
                  power == 3), power == 2))
            end do
         end do
      end do
      average = sum(moments, dim=2)/real(effective, dp)
      do condition = 1, 4
         out%moment_mean(condition) = average(condition)
         out%moment_variance(condition) = sum(moments(condition, :)**2)/ &
            real(effective, dp)
         out%moment_t(condition) = sqrt(real(effective, dp))* &
            average(condition)/sqrt(max(out%moment_variance(condition), &
            tiny(1.0_dp)))
      end do
      do power = 1, 3
         offset = 4 + (power - 1)*selected_lags
         serial = moments(offset + 1:offset + selected_lags, :)
         serial_average = average(offset + 1:offset + selected_lags)
         covariance = matmul(serial, transpose(serial))/real(effective, dp)
         allocate(inverse(selected_lags, selected_lags))
         call symmetric_pseudoinverse(covariance, inverse, inverse_info)
         if (inverse_info /= 0) then
            out%info = 2
            return
         end if
         out%serial_statistic(power) = real(effective, dp)*dot_product( &
            serial_average, matmul(inverse, serial_average))
         out%serial_p_value(power) = regularized_gamma_q(0.5_dp* &
            real(selected_lags, dp), 0.5_dp*out%serial_statistic(power))
         deallocate(inverse)
      end do
      covariance = matmul(moments, transpose(moments))/real(effective, dp)
      allocate(inverse(size(covariance, 1), size(covariance, 2)))
      call symmetric_pseudoinverse(covariance, inverse, inverse_info)
      if (inverse_info /= 0) then
         out%info = 2
         return
      end if
      out%joint_statistic = real(effective, dp)*dot_product(average, &
         matmul(inverse, average))
      out%joint_p_value = regularized_gamma_q(0.5_dp*real(size(average), dp), &
         0.5_dp*out%joint_statistic)
   end function rugarch_gmm_test

   pure function rugarch_hong_li_test(probability_transform, lags) result(out)
      !! Apply Hong-Li moment and quartic-kernel specification tests to PITs.
      real(dp), intent(in) :: probability_transform(:) !! PIT values in the unit interval.
      integer, intent(in), optional :: lags !! Positive lag count.
      type(rugarch_hong_li_test_t) :: out
      real(dp), parameter :: variance_constant = 2.0_dp*(50.0_dp/49.0_dp - &
         300.0_dp/294.0_dp + 1950.0_dp/1960.0_dp - 900.0_dp/1568.0_dp + &
         450.0_dp/2304.0_dp)**2
      real(dp) :: bandwidth, a_constant, q_sum, statistic
      integer :: selected_lags, lag, index

      selected_lags = 4
      if (present(lags)) selected_lags = lags
      if (selected_lags < 1 .or. size(probability_transform) <= &
         selected_lags + 4 .or. any(probability_transform < 0.0_dp) .or. &
         any(probability_transform > 1.0_dp) .or. &
         .not. all(ieee_is_finite(probability_transform))) then
         out%info = 1
         return
      end if
      out%lags = selected_lags
      bandwidth = sqrt(max(variance(probability_transform), tiny(1.0_dp)))* &
         real(size(probability_transform), dp)**(-1.0_dp/6.0_dp)
      a_constant = hong_li_a_constant(bandwidth)
      q_sum = 0.0_dp
      do lag = 1, selected_lags
         statistic = hong_li_density_integral(probability_transform, lag, &
            bandwidth)
         q_sum = q_sum + ((real(size(probability_transform) - lag, dp)* &
            bandwidth*statistic - bandwidth*a_constant)/ &
            sqrt(variance_constant))
      end do
      out%statistic(7) = q_sum/sqrt(real(selected_lags, dp))
      out%statistic(1) = hong_li_moment_statistic(1, 1, selected_lags, &
         probability_transform)
      out%statistic(2) = hong_li_moment_statistic(2, 2, selected_lags, &
         probability_transform)
      out%statistic(3) = hong_li_moment_statistic(3, 3, selected_lags, &
         probability_transform)
      out%statistic(4) = hong_li_moment_statistic(4, 4, selected_lags, &
         probability_transform)
      out%statistic(5) = hong_li_moment_statistic(1, 2, selected_lags, &
         probability_transform)
      out%statistic(6) = hong_li_moment_statistic(2, 1, selected_lags, &
         probability_transform)
      do index = 1, 7
         out%p_value(index) = 0.5_dp*erfc(out%statistic(index)/sqrt(2.0_dp))
      end do
   end function rugarch_hong_li_test

   pure real(dp) function duration_negative_log_likelihood(shape, duration, &
      censored) result(value)
      !! Evaluate the concentrated censored-Weibull duration likelihood.
      real(dp), intent(in) :: shape !! Weibull shape parameter.
      real(dp), intent(in) :: duration(:) !! Positive exceedance durations.
      logical, intent(in) :: censored(:) !! Left- or right-censoring flags.
      real(dp) :: scale, log_likelihood
      integer :: observation, uncensored

      uncensored = size(duration) - merge(1, 0, censored(1)) - &
         merge(1, 0, censored(size(censored)))
      scale = (real(max(1, uncensored), dp)/sum(duration**shape))** &
         (1.0_dp/shape)
      log_likelihood = 0.0_dp
      do observation = 1, size(duration)
         if (censored(observation)) then
            log_likelihood = log_likelihood - &
               (scale*duration(observation))**shape
         else
            log_likelihood = log_likelihood + shape*log(scale) + log(shape) + &
               (shape - 1.0_dp)*log(duration(observation)) - &
               (scale*duration(observation))**shape
         end if
      end do
      value = -log_likelihood
      if (.not. ieee_is_finite(value)) value = huge(1.0_dp)
   end function duration_negative_log_likelihood

   pure real(dp) function hong_li_moment_statistic(first_power, second_power, &
      lags, probability_transform) result(value)
      !! Compute one Hong-Li weighted PIT cross-moment statistic.
      integer, intent(in) :: first_power !! Power of the leading PIT.
      integer, intent(in) :: second_power !! Power of the lagged PIT.
      integer, intent(in) :: lags !! Maximum lag.
      real(dp), intent(in) :: probability_transform(:) !! PIT observations.
      real(dp) :: weight, numerator, centering, fourth_sum
      integer :: lag

      numerator = 0.0_dp
      centering = 0.0_dp
      fourth_sum = 0.0_dp
      do lag = 1, lags
         weight = 1.0_dp - real(lag, dp)/real(lags, dp)
         numerator = numerator + weight**2* &
            real(size(probability_transform) - lag, dp)* &
            (pit_cross_covariance(lag, first_power, second_power, &
            probability_transform)/max(pit_cross_covariance(0, first_power, &
            second_power, probability_transform), tiny(1.0_dp)))**2
         centering = centering + weight**2
         fourth_sum = fourth_sum + weight**4
      end do
      value = (numerator - centering)/max(fourth_sum, tiny(1.0_dp))
   end function hong_li_moment_statistic

   pure real(dp) function pit_cross_covariance(lag, first_power, second_power, &
      probability_transform) result(value)
      !! Calculate the finite-sample PIT power cross covariance at one lag.
      integer, intent(in) :: lag !! Nonnegative lag.
      integer, intent(in) :: first_power !! Leading PIT power.
      integer, intent(in) :: second_power !! Lagged PIT power.
      real(dp), intent(in) :: probability_transform(:) !! PIT observations.
      integer :: observations
      real(dp) :: first_mean, second_mean

      observations = size(probability_transform)
      first_mean = sum(probability_transform(lag + 1:)**first_power)/ &
         real(observations, dp)
      second_mean = sum(probability_transform(:observations - lag)** &
         second_power)/real(observations, dp)
      value = sum(probability_transform(lag + 1:)**first_power* &
         probability_transform(:observations - lag)**second_power)/ &
         real(observations, dp) - first_mean*second_mean
   end function pit_cross_covariance

   pure real(dp) function hong_li_density_integral(probability_transform, lag, &
      bandwidth) result(value)
      !! Integrate the squared boundary-corrected bivariate PIT density error.
      real(dp), intent(in) :: probability_transform(:) !! PIT observations.
      integer, intent(in) :: lag !! Dependence lag.
      real(dp), intent(in) :: bandwidth !! Kernel bandwidth.
      real(dp), parameter :: nodes(12) = [-0.9815606_dp, -0.9041173_dp, &
         -0.7699027_dp, -0.5873180_dp, -0.3678315_dp, -0.1252334_dp, &
         0.1252334_dp, 0.3678315_dp, 0.5873180_dp, 0.7699027_dp, &
         0.9041173_dp, 0.9815606_dp]
      real(dp), parameter :: weights(12) = [0.04717534_dp, 0.10693933_dp, &
         0.16007833_dp, 0.20316743_dp, 0.23349254_dp, 0.24914705_dp, &
         0.24914705_dp, 0.23349254_dp, 0.20316743_dp, 0.16007833_dp, &
         0.10693933_dp, 0.04717534_dp]
      real(dp) :: first, second, density
      integer :: row, column, time, observations

      observations = size(probability_transform) - lag
      value = 0.0_dp
      do row = 1, 12
         first = 0.5_dp*(nodes(row) + 1.0_dp)
         do column = 1, 12
            second = 0.5_dp*(nodes(column) + 1.0_dp)
            density = 0.0_dp
            do time = 1, observations
               density = density + boundary_kernel(first, &
                  probability_transform(time + lag), bandwidth)* &
                  boundary_kernel(second, probability_transform(time), &
                  bandwidth)
            end do
            density = density/real(observations, dp)
            value = value + weights(row)*weights(column)*(density - 1.0_dp)**2
         end do
      end do
      value = 0.25_dp*value
   end function hong_li_density_integral

   pure elemental real(dp) function boundary_kernel(point, observation, &
      bandwidth) result(value)
      !! Evaluate the quartic kernel with unit-interval boundary correction.
      real(dp), intent(in) :: point !! Evaluation point in the unit interval.
      real(dp), intent(in) :: observation !! PIT observation.
      real(dp), intent(in) :: bandwidth !! Positive bandwidth.
      real(dp) :: normalization

      value = quartic_kernel((point - observation)/bandwidth)/bandwidth
      if (point < bandwidth) then
         normalization = quartic_integral(-point/bandwidth, 1.0_dp)
         value = value/max(normalization, tiny(1.0_dp))
      else if (point > 1.0_dp - bandwidth) then
         normalization = quartic_integral(-1.0_dp, &
            (1.0_dp - point)/bandwidth)
         value = value/max(normalization, tiny(1.0_dp))
      end if
   end function boundary_kernel

   pure elemental real(dp) function quartic_kernel(value) result(kernel)
      !! Evaluate the compact quartic kernel.
      real(dp), intent(in) :: value !! Scaled kernel argument.

      if (abs(value) <= 1.0_dp) then
         kernel = 15.0_dp/16.0_dp*(1.0_dp - value**2)**2
      else
         kernel = 0.0_dp
      end if
   end function quartic_kernel

   pure elemental real(dp) function quartic_antiderivative(value) &
      result(integral)
      !! Evaluate an antiderivative of the quartic kernel on its support.
      real(dp), intent(in) :: value !! Argument clipped to the kernel support.
      real(dp) :: clipped

      clipped = max(-1.0_dp, min(1.0_dp, value))
      integral = 15.0_dp/16.0_dp*(clipped - 2.0_dp*clipped**3/3.0_dp + &
         clipped**5/5.0_dp)
   end function quartic_antiderivative

   pure elemental real(dp) function quartic_integral(lower, upper) &
      result(integral)
      !! Integrate the quartic kernel between two bounds.
      real(dp), intent(in) :: lower !! Lower integration bound.
      real(dp), intent(in) :: upper !! Upper integration bound.

      integral = quartic_antiderivative(upper) - quartic_antiderivative(lower)
   end function quartic_integral

   pure real(dp) function hong_li_a_constant(bandwidth) result(value)
      !! Compute the Hong-Li boundary adjustment constant.
      real(dp), intent(in) :: bandwidth !! Kernel bandwidth.
      integer, parameter :: intervals = 1024
      real(dp) :: point, width, first, second, integral
      integer :: index

      width = 1.0_dp/real(intervals, dp)
      integral = 0.0_dp
      do index = 1, intervals
         point = (real(index, dp) - 0.5_dp)*width
         first = (8.0_dp/15.0_dp + point - 2.0_dp*point**3/3.0_dp + &
            point**5/5.0_dp)**(-2)
         second = point*(1.0_dp - point**2)**4 + 128.0_dp/315.0_dp + &
            8.0_dp*point**3/3.0_dp - 24.0_dp*point**5/5.0_dp + &
            24.0_dp*point**7/7.0_dp - 8.0_dp*point**9/9.0_dp
         integral = integral + first*second
      end do
      integral = integral*width
      value = ((1.0_dp/bandwidth - 2.0_dp)*5.0_dp/7.0_dp + &
         2.0_dp*integral)**2 - 1.0_dp
   end function hong_li_a_constant

   pure function diagonal(matrix) result(values)
      !! Extract the diagonal of a square real matrix.
      real(dp), intent(in) :: matrix(:, :) !! Square input matrix.
      real(dp) :: values(min(size(matrix, 1), size(matrix, 2)))
      integer :: index

      do index = 1, size(values)
         values(index) = matrix(index, index)
      end do
   end function diagonal

   pure real(dp) function sum_diagonal(matrix) result(value)
      !! Sum the diagonal of a square real matrix.
      real(dp), intent(in) :: matrix(:, :) !! Square input matrix.

      value = sum(diagonal(matrix))
   end function sum_diagonal

   pure function nyblom_critical(parameters) result(values)
      !! Return Hansen's 10%, 5%, and 1% joint critical values.
      integer, intent(in) :: parameters !! Parameter count from one through twenty.
      real(dp) :: values(3)
      real(dp), parameter :: table(3, 20) = reshape([ &
         0.353_dp, 0.470_dp, 0.748_dp, 0.610_dp, 0.749_dp, 1.07_dp, &
         0.846_dp, 1.01_dp, 1.35_dp, 1.07_dp, 1.24_dp, 1.60_dp, &
         1.28_dp, 1.47_dp, 1.88_dp, 1.49_dp, 1.68_dp, 2.12_dp, &
         1.69_dp, 1.90_dp, 2.35_dp, 1.89_dp, 2.11_dp, 2.59_dp, &
         2.10_dp, 2.32_dp, 2.82_dp, 2.29_dp, 2.54_dp, 3.05_dp, &
         2.49_dp, 2.75_dp, 3.27_dp, 2.69_dp, 2.96_dp, 3.51_dp, &
         2.89_dp, 3.15_dp, 3.69_dp, 3.08_dp, 3.34_dp, 3.90_dp, &
         3.26_dp, 3.54_dp, 4.07_dp, 3.46_dp, 3.75_dp, 4.30_dp, &
         3.64_dp, 3.95_dp, 4.51_dp, 3.83_dp, 4.14_dp, 4.73_dp, &
         4.03_dp, 4.33_dp, 4.92_dp, 4.22_dp, 4.52_dp, 5.13_dp], [3, 20])

      values = table(:, max(1, min(20, parameters)))
   end function nyblom_critical

   subroutine display_rugarch_nyblom_test(value, unit)
      !! Display Hansen-Nyblom stability statistics.
      type(rugarch_nyblom_test_t), intent(in) :: value !! Stability-test result.
      integer, intent(in), optional :: unit !! Output unit.
      integer :: output

      output = output_unit
      if (present(unit)) output = unit
      write(output, '(a)') 'rugarch Nyblom stability test'
      write(output, '(a,f12.5)') 'joint statistic: ', value%joint_statistic
      write(output, '(a,3(1x,f8.3))') 'joint critical values:', &
         value%joint_critical
   end subroutine display_rugarch_nyblom_test

   subroutine display_rugarch_gof_test(value, unit)
      !! Display grouped probability-transform goodness-of-fit tests.
      type(rugarch_gof_test_t), intent(in) :: value !! Goodness-of-fit result.
      integer, intent(in), optional :: unit !! Output unit.
      integer :: output, test

      output = output_unit
      if (present(unit)) output = unit
      write(output, '(a)') 'rugarch grouped goodness-of-fit test'
      do test = 1, size(value%groups)
         write(output, '(i6,2(1x,f12.5))') value%groups(test), &
            value%statistic(test), value%p_value(test)
      end do
   end subroutine display_rugarch_gof_test

   subroutine display_rugarch_var_duration_test(value, unit)
      !! Display value-at-risk duration test results.
      type(rugarch_var_duration_test_t), intent(in) :: value !! Duration-test result.
      integer, intent(in), optional :: unit !! Output unit.
      integer :: output

      output = output_unit
      if (present(unit)) output = unit
      write(output, '(a)') 'rugarch VaR duration test'
      write(output, '(a,f12.5)') 'Weibull shape: ', value%shape
      write(output, '(a,f12.5)') 'likelihood ratio: ', value%likelihood_ratio
      write(output, '(a,f12.5)') 'p-value: ', value%p_value
   end subroutine display_rugarch_var_duration_test

   subroutine display_rugarch_gmm_test(value, unit)
      !! Display GMM moment-test results.
      type(rugarch_gmm_test_t), intent(in) :: value !! GMM-test result.
      integer, intent(in), optional :: unit !! Output unit.
      integer :: output

      output = output_unit
      if (present(unit)) output = unit
      write(output, '(a)') 'rugarch GMM test'
      write(output, '(a,4(1x,f10.4))') 'moment t values:', value%moment_t
      write(output, '(a,f12.5)') 'joint statistic: ', value%joint_statistic
      write(output, '(a,f12.5)') 'joint p-value: ', value%joint_p_value
   end subroutine display_rugarch_gmm_test

   subroutine display_rugarch_hong_li_test(value, unit)
      !! Display Hong-Li PIT specification tests.
      type(rugarch_hong_li_test_t), intent(in) :: value !! Hong-Li result.
      integer, intent(in), optional :: unit !! Output unit.
      integer :: output

      output = output_unit
      if (present(unit)) output = unit
      write(output, '(a)') 'rugarch Hong-Li test'
      write(output, '(a,7(1x,f10.4))') 'statistics:', value%statistic
      write(output, '(a,7(1x,f10.4))') 'p-values:', value%p_value
   end subroutine display_rugarch_hong_li_test

end module rugarch_diagnostics_mod
