program compare_rfortran_core
   use r_kinds, only : dp
   use r_descriptive, only : r_correlation, r_covariance
   use r_descriptive, only : r_weighted_correlation, r_weighted_covariance
   use r_descriptive, only : r_weighted_mean, r_weighted_variance
   use r_rolling, only : r_roll_correlation_right, r_roll_covariance_valid, r_roll_mean_right
   use r_rolling, only : r_roll_max_valid, r_roll_mean_valid, r_roll_min_right
   use r_rolling, only : r_roll_sd_right, r_roll_sum_valid, r_roll_variance_valid
   use r_quantiles, only : r_weighted_quantile_ecdf, r_weighted_quantile_frequency_type7
   use r_quantiles, only : r_weighted_quantile_isotone, r_weighted_quantile_linear_cdf
   use r_quantiles, only : r_weighted_quantile_survey
   use r_quantiles, only : r_median, r_quantile_type7
   use r_robust, only : r_mad
   use r_distributions, only : r_df, r_dchisq, r_dt, r_pf, r_pchisq, r_pt
   use r_distributions, only : r_qf, r_qchisq, r_qnorm, r_qt
   use r_ordering, only : r_order
   use r_special, only : r_digamma, r_log_beta, r_log_choose, r_log_factorial, r_trigamma
   use r_special, only : r_regularized_beta, r_regularized_gamma_p, r_regularized_gamma_q
   use r_stability, only : r_log_mean_exp, r_log_sum_exp
   use r_transforms, only : r_expm1, r_log1mexp, r_log1p, r_log1pexp, r_logistic, r_logit
   implicit none

   real(dp), parameter :: points(6) = [0.1_dp, 0.5_dp, 1.0_dp, 2.5_dp, 10.0_dp, 100.0_dp]
   real(dp), parameter :: log1p_points(6) = [-0.9_dp, -1.0e-12_dp, 0.0_dp, 1.0e-12_dp, 1.0_dp, 1.0e6_dp]
   real(dp), parameter :: expm1_points(6) = [-100.0_dp, -1.0_dp, -1.0e-12_dp, 1.0e-12_dp, 1.0_dp, 10.0_dp]
   real(dp), parameter :: log1mexp_points(4) = [-100.0_dp, -10.0_dp, -1.0_dp, -1.0e-12_dp]
   real(dp), parameter :: softplus_points(5) = [-100.0_dp, -1.0_dp, 0.0_dp, 1.0_dp, 100.0_dp]
   real(dp), parameter :: log_weight_points(6) = [-1000.0_dp, -100.0_dp, -2.0_dp, 0.0_dp, 1.0_dp, 100.0_dp]
   real(dp), parameter :: probability_points(6) = [1.0e-12_dp, 0.01_dp, 0.25_dp, &
      0.5_dp, 0.99_dp, 1.0_dp - 1.0e-12_dp]
   real(dp), parameter :: log_probability_points(4) = [-1.0_dp, -10.0_dp, -100.0_dp, -1000.0_dp]
   real(dp), parameter :: quantile_x(6) = [8.0_dp, 1.0_dp, 5.0_dp, 3.0_dp, 9.0_dp, 2.0_dp]
   real(dp), parameter :: order_x(6) = [3.0_dp, 1.0_dp, 3.0_dp, 2.0_dp, 1.0_dp, 3.0_dp]
   real(dp), parameter :: gamma_shapes(5) = [0.5_dp, 1.0_dp, 2.5_dp, 10.0_dp, 100.0_dp]
   real(dp), parameter :: gamma_x(5) = [0.1_dp, 1.0_dp, 3.0_dp, 8.0_dp, 110.0_dp]
   real(dp), parameter :: beta_x(5) = [0.01_dp, 0.1_dp, 0.5_dp, 0.8_dp, 0.99_dp]
   real(dp), parameter :: beta_shape1(5) = [0.5_dp, 1.0_dp, 2.0_dp, 5.0_dp, 20.0_dp]
   real(dp), parameter :: beta_shape2(5) = [2.0_dp, 4.0_dp, 3.0_dp, 1.5_dp, 10.0_dp]
   real(dp), parameter :: distribution_x(6) = [0.1_dp, 0.5_dp, 1.0_dp, 2.0_dp, 5.0_dp, 10.0_dp]
   real(dp), parameter :: student_x(6) = [-3.0_dp, -1.0_dp, 0.0_dp, 0.5_dp, 2.0_dp, 5.0_dp]
   real(dp), parameter :: distribution_probs(6) = [0.01_dp, 0.1_dp, 0.25_dp, 0.5_dp, 0.9_dp, 0.99_dp]
   real(dp), parameter :: distribution_df(6) = [1.0_dp, 2.0_dp, 5.0_dp, 10.0_dp, 30.0_dp, 100.0_dp]
   real(dp), parameter :: f_df1(6) = [1.0_dp, 2.0_dp, 5.0_dp, 10.0_dp, 20.0_dp, 50.0_dp]
   real(dp), parameter :: f_df2(6) = [2.0_dp, 5.0_dp, 8.0_dp, 12.0_dp, 30.0_dp, 100.0_dp]
   real(dp), parameter :: quantile_weights(6) = [1.0_dp, 2.0_dp, 4.0_dp, 1.0_dp, 3.0_dp, 2.0_dp]
   real(dp), parameter :: quantile_probs(6) = [0.0_dp, 0.1_dp, 0.25_dp, 0.5_dp, 0.9_dp, 1.0_dp]
   integer, parameter :: factorial_points(5) = [0, 1, 5, 20, 1000]
   integer, parameter :: choose_n(4) = [5, 20, 1000, 1000000]
   integer, parameter :: choose_k(4) = [2, 17, 500, 2]
   real(dp) :: pair_y(100), weights(100), x(100), value, t0, t1
   real(dp), allocatable :: rolling(:)
   integer, allocatable :: order_indices(:)
   integer :: i, j, reps, unit
   character(512) :: output

   call get_command_argument(1, output)
   if (len_trim(output) == 0) output = 'fortran_results.csv'
   x = [(real(i, dp), i=1,100)]
   pair_y = 0.5_dp*x + sin(x/7.0_dp)
   weights = [(real(1 + mod(i,7), dp), i=1,100)]
   reps = 20000

   open(newunit=unit, file=trim(output), status='replace')
   write(unit, '(a)') 'case,value,seconds,abs_tol,rel_tol'

   call cpu_time(t0)
   do i = 1, reps
      value = sum(r_digamma(points)*[(real(j, dp), j=1,size(points))])
   end do
   call cpu_time(t1)
   call emit('digamma_checksum', value, t1 - t0)

   call cpu_time(t0)
   do i = 1, reps
      value = sum(r_trigamma(points)*[(real(j, dp), j=1,size(points))])
   end do
   call cpu_time(t1)
   call emit('trigamma_checksum', value, t1 - t0)

   call cpu_time(t0)
   do i = 1, reps
      value = sum(r_log_beta(points, points(size(points):1:-1))*[(real(j, dp), j=1,size(points))])
   end do
   call cpu_time(t1)
   call emit('log_beta_checksum', value, t1 - t0)

   call cpu_time(t0)
   do i = 1, reps
      call r_roll_mean_valid(x, 5, rolling)
      value = sum(rolling*[(real(j, dp), j=1,size(rolling))])
   end do
   call cpu_time(t1)
   call emit('rolling_mean_valid', value, t1 - t0)

   call cpu_time(t0)
   do i = 1, reps
      call r_roll_mean_right(x, 5, rolling)
      value = sum(rolling(5:)*[(real(j, dp), j=5,100)])
   end do
   call cpu_time(t1)
   call emit('rolling_mean_right', value, t1 - t0)

   call cpu_time(t0)
   do i = 1, reps
      value = sum(r_log1p(log1p_points)*[(real(j, dp), j=1,size(log1p_points))])
   end do
   call cpu_time(t1)
   call emit('log1p_checksum', value, t1 - t0)

   call cpu_time(t0)
   do i = 1, reps
      value = sum(r_expm1(expm1_points)*[(real(j, dp), j=1,size(expm1_points))])
   end do
   call cpu_time(t1)
   call emit('expm1_checksum', value, t1 - t0)

   call cpu_time(t0)
   do i = 1, reps
      value = sum(r_log1mexp(log1mexp_points)*[(real(j, dp), j=1,size(log1mexp_points))])
   end do
   call cpu_time(t1)
   call emit('log1mexp_checksum', value, t1 - t0)

   call cpu_time(t0)
   do i = 1, reps
      value = sum(r_log1pexp(softplus_points)*[(real(j, dp), j=1,size(softplus_points))])
   end do
   call cpu_time(t1)
   call emit('log1pexp_checksum', value, t1 - t0)

   call cpu_time(t0)
   do i = 1, reps
      value = r_log_sum_exp(log_weight_points)
   end do
   call cpu_time(t1)
   call emit('log_sum_exp', value, t1 - t0)

   call cpu_time(t0)
   do i = 1, reps
      value = r_log_mean_exp(log_weight_points)
   end do
   call cpu_time(t1)
   call emit('log_mean_exp', value, t1 - t0)

   call cpu_time(t0)
   do i = 1, reps
      value = sum(r_logistic(softplus_points)*[(real(j, dp), j=1,size(softplus_points))])
   end do
   call cpu_time(t1)
   call emit('logistic_checksum', value, t1 - t0)

   call cpu_time(t0)
   do i = 1, reps
      value = sum(r_logit(probability_points)*[(real(j, dp), j=1,size(probability_points))])
   end do
   call cpu_time(t1)
   call emit('logit_checksum', value, t1 - t0)

   call cpu_time(t0)
   do i = 1, reps
      value = sum(r_qnorm(probability_points)*[(real(j, dp), j=1,size(probability_points))])
   end do
   call cpu_time(t1)
   call emit('qnorm_checksum', value, t1 - t0)

   call cpu_time(t0)
   do i = 1, reps
      value = sum(r_qnorm(log_probability_points, log_probability=.true.)* &
         [(real(j, dp), j=1,size(log_probability_points))])
   end do
   call cpu_time(t1)
   call emit('qnorm_log_probability_checksum', value, t1 - t0)

   reps = 5000
   call cpu_time(t0)
   do i = 1, reps
      value = sum(r_dt(student_x,distribution_df)*[(real(j,dp),j=1,size(student_x))])
   end do
   call cpu_time(t1)
   call emit('student_t_density_checksum',value,t1-t0)

   call cpu_time(t0)
   do i = 1, reps
      value = sum(r_pt(student_x,distribution_df)*[(real(j,dp),j=1,size(student_x))])
   end do
   call cpu_time(t1)
   call emit('student_t_cdf_checksum',value,t1-t0)

   call cpu_time(t0)
   do i = 1, reps
      value = sum(r_pt(abs(student_x),distribution_df,lower_tail=.false., &
         log_probability=.true.)*[(real(j,dp),j=1,size(student_x))])
   end do
   call cpu_time(t1)
   call emit('student_t_log_survival_checksum',value,t1-t0)

   reps = 100
   call cpu_time(t0)
   do i = 1, reps
      value = sum(r_qt(distribution_probs,distribution_df)*[(real(j,dp),j=1,size(distribution_probs))])
   end do
   call cpu_time(t1)
   call emit('student_t_quantile_checksum',value,t1-t0)

   reps = 5000
   call cpu_time(t0)
   do i = 1, reps
      value = sum(r_dchisq(distribution_x,distribution_df)*[(real(j,dp),j=1,size(distribution_x))])
   end do
   call cpu_time(t1)
   call emit('chi_square_density_checksum',value,t1-t0)

   call cpu_time(t0)
   do i = 1, reps
      value = sum(r_pchisq(distribution_x,distribution_df)*[(real(j,dp),j=1,size(distribution_x))])
   end do
   call cpu_time(t1)
   call emit('chi_square_cdf_checksum',value,t1-t0)

   call cpu_time(t0)
   do i = 1, reps
      value = sum(r_pchisq(distribution_x,distribution_df,lower_tail=.false., &
         log_probability=.true.)*[(real(j,dp),j=1,size(distribution_x))])
   end do
   call cpu_time(t1)
   call emit('chi_square_log_survival_checksum',value,t1-t0)

   reps = 100
   call cpu_time(t0)
   do i = 1, reps
      value = sum(r_qchisq(distribution_probs,distribution_df)*[(real(j,dp),j=1,size(distribution_probs))])
   end do
   call cpu_time(t1)
   call emit('chi_square_quantile_checksum',value,t1-t0)

   reps = 5000
   call cpu_time(t0)
   do i = 1, reps
      value = sum(r_df(distribution_x,f_df1,f_df2)*[(real(j,dp),j=1,size(distribution_x))])
   end do
   call cpu_time(t1)
   call emit('f_density_checksum',value,t1-t0)

   call cpu_time(t0)
   do i = 1, reps
      value = sum(r_pf(distribution_x,f_df1,f_df2)*[(real(j,dp),j=1,size(distribution_x))])
   end do
   call cpu_time(t1)
   call emit('f_cdf_checksum',value,t1-t0)

   call cpu_time(t0)
   do i = 1, reps
      value = sum(r_pf(distribution_x,f_df1,f_df2,lower_tail=.false., &
         log_probability=.true.)*[(real(j,dp),j=1,size(distribution_x))])
   end do
   call cpu_time(t1)
   call emit('f_log_survival_checksum',value,t1-t0)

   reps = 100
   call cpu_time(t0)
   do i = 1, reps
      value = sum(r_qf(distribution_probs,f_df1,f_df2)*[(real(j,dp),j=1,size(distribution_probs))])
   end do
   call cpu_time(t1)
   call emit('f_quantile_checksum',value,t1-t0)

   reps = 20000
   call cpu_time(t0)
   do i = 1, reps
      value = r_quantile_type7(quantile_x, 0.37_dp)
   end do
   call cpu_time(t1)
   call emit('quantile_type7', value, t1 - t0)

   call cpu_time(t0)
   do i = 1, reps
      value = r_median(quantile_x)
   end do
   call cpu_time(t1)
   call emit('median', value, t1 - t0)

   call cpu_time(t0)
   do i = 1, reps
      value = r_mad(quantile_x)
   end do
   call cpu_time(t1)
   call emit('mad', value, t1 - t0)

   call cpu_time(t0)
   do i = 1, reps
      value = sum(r_regularized_gamma_p(gamma_shapes, gamma_x)* &
         [(real(j, dp), j=1,size(gamma_shapes))])
   end do
   call cpu_time(t1)
   call emit('regularized_gamma_p_checksum', value, t1 - t0)

   call cpu_time(t0)
   do i = 1, reps
      value = sum(r_regularized_gamma_q(gamma_shapes, gamma_x)* &
         [(real(j, dp), j=1,size(gamma_shapes))])
   end do
   call cpu_time(t1)
   call emit('regularized_gamma_q_checksum', value, t1 - t0)

   call cpu_time(t0)
   do i = 1, reps
      value = sum(r_regularized_beta(beta_x, beta_shape1, beta_shape2)* &
         [(real(j, dp), j=1,size(beta_x))])
   end do
   call cpu_time(t1)
   call emit('regularized_beta_checksum', value, t1 - t0)

   reps = 5000
   call cpu_time(t0)
   do i = 1, reps
      call r_order(order_x, order_indices)
      value = sum(real(order_indices, dp)*[(real(j, dp), j=1,size(order_indices))])
   end do
   call cpu_time(t1)
   call emit('stable_order_checksum', value, t1 - t0)

   reps = 1000
   call cpu_time(t0)
   do i = 1, reps
      call r_roll_variance_valid(x, 5, rolling)
      value = sum(rolling*[(real(j, dp), j=1,size(rolling))])
   end do
   call cpu_time(t1)
   call emit('rolling_variance_valid', value, t1 - t0)

   call cpu_time(t0)
   do i = 1, reps
      call r_roll_sd_right(x, 5, rolling)
      value = sum(rolling(5:)*[(real(j, dp), j=5,100)])
   end do
   call cpu_time(t1)
   call emit('rolling_sd_right', value, t1 - t0)

   call cpu_time(t0)
   do i = 1, reps
      call r_roll_sum_valid(x, 5, rolling)
      value = sum(rolling*[(real(j, dp), j=1,size(rolling))])
   end do
   call cpu_time(t1)
   call emit('rolling_sum_valid', value, t1 - t0)

   call cpu_time(t0)
   do i = 1, reps
      call r_roll_min_right(x, 5, rolling)
      value = sum(rolling(5:)*[(real(j, dp), j=5,100)])
   end do
   call cpu_time(t1)
   call emit('rolling_min_right', value, t1 - t0)

   call cpu_time(t0)
   do i = 1, reps
      call r_roll_max_valid(x, 5, rolling)
      value = sum(rolling*[(real(j, dp), j=1,size(rolling))])
   end do
   call cpu_time(t1)
   call emit('rolling_max_valid', value, t1 - t0)

   reps = 20000
   call cpu_time(t0)
   do i = 1, reps
      value = r_covariance(x, pair_y)
   end do
   call cpu_time(t1)
   call emit('covariance', value, t1 - t0)

   call cpu_time(t0)
   do i = 1, reps
      value = r_correlation(x, pair_y)
   end do
   call cpu_time(t1)
   call emit('correlation', value, t1 - t0)

   reps = 5000
   call cpu_time(t0)
   do i = 1, reps
      value = r_weighted_mean(x, weights)
   end do
   call cpu_time(t1)
   call emit('weighted_mean', value, t1 - t0)

   call cpu_time(t0)
   do i = 1, reps
      value = r_weighted_variance(x, weights)
   end do
   call cpu_time(t1)
   call emit('weighted_variance_ml', value, t1 - t0)

   call cpu_time(t0)
   do i = 1, reps
      value = r_weighted_covariance(x, pair_y, weights, unbiased=.true.)
   end do
   call cpu_time(t1)
   call emit('weighted_covariance_unbiased', value, t1 - t0)

   call cpu_time(t0)
   do i = 1, reps
      value = r_weighted_correlation(x, pair_y, weights)
   end do
   call cpu_time(t1)
   call emit('weighted_correlation', value, t1 - t0)

   call cpu_time(t0)
   do i = 1, reps
      value = 0.0_dp
      do j = 1, size(quantile_probs)
         value = value + real(j,dp)*r_weighted_quantile_ecdf(quantile_x,quantile_weights,quantile_probs(j))
      end do
   end do
   call cpu_time(t1)
   call emit('weighted_quantile_ecdf', value, t1 - t0)

   call cpu_time(t0)
   do i = 1, reps
      value = 0.0_dp
      do j = 1, size(quantile_probs)
         value = value + real(j,dp)*r_weighted_quantile_linear_cdf(quantile_x,quantile_weights,quantile_probs(j))
      end do
   end do
   call cpu_time(t1)
   call emit('weighted_quantile_linear_cdf', value, t1 - t0)

   call cpu_time(t0)
   do i = 1, reps
      value = 0.0_dp
      do j = 1, size(quantile_probs)
         value = value + real(j,dp)*r_weighted_quantile_frequency_type7(quantile_x,quantile_weights,quantile_probs(j))
      end do
   end do
   call cpu_time(t1)
   call emit('weighted_quantile_frequency_type7', value, t1 - t0)

   call cpu_time(t0)
   do i = 1, reps
      value = 0.0_dp
      do j = 1, size(quantile_probs)
         value = value + real(j,dp)*r_weighted_quantile_isotone(quantile_x,quantile_weights,quantile_probs(j))
      end do
   end do
   call cpu_time(t1)
   call emit('weighted_quantile_isotone', value, t1 - t0)

   call cpu_time(t0)
   do i = 1, reps
      value = 0.0_dp
      do j = 1, 12
         value = value + real(j,dp)*r_weighted_quantile_survey(quantile_x,quantile_weights,0.37_dp,j)
      end do
   end do
   call cpu_time(t1)
   call emit('weighted_quantile_survey_rules', value, t1 - t0)

   reps = 1000
   call cpu_time(t0)
   do i = 1, reps
      call r_roll_covariance_valid(x, pair_y, 5, rolling)
      value = sum(rolling*[(real(j, dp), j=1,size(rolling))])
   end do
   call cpu_time(t1)
   call emit('rolling_covariance_valid', value, t1 - t0)

   call cpu_time(t0)
   do i = 1, reps
      call r_roll_correlation_right(x, pair_y, 5, rolling)
      value = sum(rolling(5:)*[(real(j, dp), j=5,100)])
   end do
   call cpu_time(t1)
   call emit('rolling_correlation_right', value, t1 - t0)

   reps = 20000
   call cpu_time(t0)
   do i = 1, reps
      value = sum(r_log_factorial(factorial_points)*[(real(j, dp), j=1,size(factorial_points))])
   end do
   call cpu_time(t1)
   call emit('log_factorial_checksum', value, t1 - t0)

   call cpu_time(t0)
   do i = 1, reps
      value = sum(r_log_choose(choose_n, choose_k)*[(real(j, dp), j=1,size(choose_n))])
   end do
   call cpu_time(t1)
   call emit('log_choose_checksum', value, t1 - t0)

   close(unit)

contains

   subroutine emit(name, result, seconds)
      character(*), intent(in) :: name
      real(dp), intent(in) :: result, seconds

      write(unit, '(a,",",es25.16e3,",",es16.8,",",es12.4,",",es12.4)') &
         trim(name), result, seconds, 1.0e-11_dp, 1.0e-11_dp
   end subroutine emit

end program compare_rfortran_core
