program parity_driver
   use flexmix, only : dp, flexmix_result, flexmix_control
   use flexmix, only : fit_gaussian_regression_component, fit_poisson_regression_component
   use flexmix, only : fit_binomial_regression_component, fit_mvnormal_component
   use flexmix, only : flexmix_gamma, flexmix_poisson, flexmix_multinomial
   use flexmix, only : flexmix_glmfix_gaussian, fit_factor_analysis_component
   use flexmix, only : flexmix_glmnet_gaussian, flexmix_lmm, flexmix_lmmc, flexmix_lmc
   use flexmix, only : flexmix_mgcv_gaussian, flexmix_mgcv_poisson, flexmix_mgcv_binomial
   use flexmix, only : fit_conditional_logit_component
   implicit none
   real(dp) :: x(5,2), yg(5), yp(5), success(5), trials(5), weights(5)
   real(dp) :: beta(2), sigma, center(2), covariance(2,2), ym(5,2)
   real(dp) :: x_gamma(4,1), y_gamma(4), x_offset(5,1), y_offset(5), offset(5)
   real(dp) :: x_multi(10,1)
   integer :: y_class(10)
   real(dp) :: xfix(20,3,2), yfix(20), tfix
   real(dp) :: y_factor(60,4), w_factor(60), center_factor(4), covariance_factor(4,4)
   real(dp) :: marginal_factor(4), loadings_factor(4,1), uniqueness_factor(4), z1, z2
   real(dp) :: x_cond(36,1), y_cond(36), w_cond(36), beta_cond(1)
   integer :: strata_cond(36), selected, r, s
   logical :: design_fix(2,3)
   integer :: initial_fix(20)
   type(flexmix_result) :: gamma_fit, offset_fit, multi_fit, fixed_fit
   type(flexmix_control) :: ctl
   real(dp) :: x_pen(20,4), y_pen(20), penalty_grid(1)
   real(dp) :: x_lmm(48,2), z_lmm(48,1), y_lmm(48), y_lmmc(48), random_effect
   real(dp) :: x_multi_cens(12,2), z_multi_cens(12,1), y_multi_cens(12)
   integer :: group_lmm(48), group_multi_cens(12), g, rr
   logical :: censored_lmm(48), censored_multi(12)
   type(flexmix_result) :: glmnet_fit, lmm_fit, lmmc_fit, lmmc_multi_fit, lmc_fit
   real(dp) :: x_smooth(10,3), y_smooth(10), y_smooth_pois(10), smooth_success(10), smooth_failure(10)
   real(dp) :: smooth_penalty(3,3), smooth_lambda(1), smooth_t
   type(flexmix_result) :: smooth_gfit, smooth_pfit, smooth_bfit
   integer :: df, info, unit, i

   x(:,1) = 1.0_dp
   x(:,2) = [-2.0_dp, -1.0_dp, 0.0_dp, 1.0_dp, 2.0_dp]
   yg = [1.2_dp, 1.8_dp, 3.1_dp, 4.2_dp, 5.1_dp]
   yp = [1.0_dp, 1.0_dp, 2.0_dp, 4.0_dp, 7.0_dp]
   success = [1.0_dp, 2.0_dp, 4.0_dp, 7.0_dp, 9.0_dp]
   trials = 10.0_dp
   weights = 1.0_dp
   ym(:,1) = [-1.0_dp, 0.0_dp, 1.0_dp, 2.0_dp, 3.0_dp]
   ym(:,2) = [2.0_dp, 1.0_dp, 4.0_dp, 3.0_dp, 5.0_dp]

   open(newunit=unit, file='fortran_parity.csv', status='replace', action='write')
   write(unit,'(a)') 'metric,value'

   call fit_gaussian_regression_component(x, yg, weights, beta, sigma, df, info)
   call emit(unit, 'gaussian_beta0', beta(1))
   call emit(unit, 'gaussian_beta1', beta(2))
   call emit(unit, 'gaussian_sigma', sigma)

   call fit_poisson_regression_component(x, yp, weights, beta, df, info)
   call emit(unit, 'poisson_beta0', beta(1))
   call emit(unit, 'poisson_beta1', beta(2))

   call fit_binomial_regression_component(x, success, trials, weights, beta, df, info)
   call emit(unit, 'binomial_beta0', beta(1))
   call emit(unit, 'binomial_beta1', beta(2))

   call fit_mvnormal_component(ym, weights, center, covariance, .false., df, info)
   do i = 1, 2
      write(unit,'(a,i0,a,es24.16)') 'mvnorm_center', i, ',', center(i)
   end do
   call emit(unit, 'mvnorm_cov11', covariance(1,1))
   call emit(unit, 'mvnorm_cov12', covariance(1,2))
   call emit(unit, 'mvnorm_cov22', covariance(2,2))

   x_gamma = 1.0_dp
   y_gamma = [1.0_dp, 2.0_dp, 4.0_dp, 8.0_dp]
   call flexmix_gamma(x_gamma, y_gamma, 1, gamma_fit)
   call emit(unit, 'gamma_inverse_beta0', gamma_fit%beta(1,1))
   call emit(unit, 'gamma_shape', gamma_fit%shape(1))

   x_offset = 1.0_dp
   y_offset = [1.0_dp, 2.0_dp, 4.0_dp, 8.0_dp, 16.0_dp]
   offset = log([1.0_dp, 1.0_dp, 2.0_dp, 2.0_dp, 4.0_dp])
   call flexmix_poisson(x_offset, y_offset, 1, offset_fit, offset=offset)
   call emit(unit, 'poisson_offset_beta0', offset_fit%beta(1,1))

   x_multi = 1.0_dp
   y_class = [1,1,2,2,2,3,3,3,3,3]
   call flexmix_multinomial(x_multi, y_class, 1, multi_fit)
   call emit(unit, 'multinomial_beta_class2', multi_fit%multinomial_coef(1,1,1))
   call emit(unit, 'multinomial_beta_class3', multi_fit%multinomial_coef(1,2,1))

   xfix = 0.0_dp
   design_fix(1,:) = [.true., .false., .true.]
   design_fix(2,:) = [.false., .true., .true.]
   do i = 1, 20
      tfix = real(mod(i-1,10),dp) / 9.0_dp
      xfix(i,1,1) = 1.0_dp
      xfix(i,3,1) = tfix
      xfix(i,2,2) = 1.0_dp
      xfix(i,3,2) = tfix
      if (i <= 10) then
         yfix(i) = 1.0_dp + 2.0_dp*tfix + 0.02_dp*real((-1)**i,dp)
         initial_fix(i) = 1
      else
         yfix(i) = 8.0_dp + 2.0_dp*tfix + 0.02_dp*real((-1)**i,dp)
         initial_fix(i) = 2
      end if
   end do
   ctl%minprior = 0.0_dp
   ctl%tolerance = 1.0e-10_dp
   call flexmix_glmfix_gaussian(xfix, yfix, design_fix, fixed_fit, ctl, initial_cluster=initial_fix)
   call emit(unit, 'glmfix_intercept1', fixed_fit%beta(1,1))
   call emit(unit, 'glmfix_intercept2', fixed_fit%beta(2,2))
   call emit(unit, 'glmfix_shared_slope', fixed_fit%beta(3,1))

   do s = 1, 12
      selected = 1
      if (s > 3 .and. s <= 7) selected = 2
      if (s > 7) selected = 3
      do r = 1, 3
         strata_cond((s-1)*3+r) = s
         x_cond((s-1)*3+r,1) = real(r-2,dp)
         y_cond((s-1)*3+r) = merge(1.0_dp,0.0_dp,r == selected)
      end do
   end do
   w_cond = 1.0_dp
   call fit_conditional_logit_component(x_cond, y_cond, strata_cond, w_cond, beta_cond, df, info)
   call emit(unit, 'conditional_logit_beta1', beta_cond(1))

   do i = 1, 60
      z1 = sin(0.37_dp * real(i,dp)) + 0.3_dp * cos(0.11_dp * real(i,dp))
      z2 = cos(0.29_dp * real(i,dp))
      y_factor(i,1) = 2.0_dp + z1 + 0.18_dp * sin(1.7_dp * real(i,dp))
      y_factor(i,2) = -1.0_dp + 0.85_dp * z1 + 0.25_dp * cos(1.3_dp * real(i,dp))
      y_factor(i,3) = 0.5_dp + 0.65_dp * z1 + 0.35_dp * sin(1.1_dp * real(i,dp))
      y_factor(i,4) = 3.0_dp + 0.45_dp * z1 + 0.55_dp * z2
   end do
   w_factor = 1.0_dp
   call fit_factor_analysis_component(y_factor, w_factor, 1, center_factor, covariance_factor, marginal_factor, &
                                      loadings_factor, uniqueness_factor, df, info)
   do i = 1, 4
      write(unit,'(a,i0,a,es24.16)') 'factanal_uniqueness', i, ',', uniqueness_factor(i)
   end do
   do i = 1, 20
      x_pen(i,1) = 1.0_dp
      x_pen(i,2) = (real(i,dp) - 10.5_dp) / 5.0_dp
      x_pen(i,3) = sin(real(i,dp))
      x_pen(i,4) = cos(0.3_dp*real(i,dp))
      y_pen(i) = 2.0_dp + 3.0_dp*x_pen(i,2) + 0.02_dp*sin(2.0_dp*real(i,dp))
   end do
   penalty_grid = 0.05_dp
   call flexmix_glmnet_gaussian(x_pen, y_pen, 1, glmnet_fit, adaptive=.false., alpha=1.0_dp, &
                                lambda_grid=penalty_grid, nfolds=2)
   call emit(unit, 'glmnet_gaussian_intercept', glmnet_fit%beta(1,1))
   call emit(unit, 'glmnet_gaussian_slope1', glmnet_fit%beta(2,1))
   call emit(unit, 'glmnet_gaussian_slope2', glmnet_fit%beta(3,1))
   call emit(unit, 'glmnet_gaussian_slope3', glmnet_fit%beta(4,1))

   i = 0
   do g = 1, 12
      random_effect = 0.8_dp*sin(0.7_dp*real(g,dp))
      do rr = 1, 4
         i = i + 1
         group_lmm(i) = g
         x_lmm(i,1) = 1.0_dp
         x_lmm(i,2) = real(rr-1,dp) / 3.0_dp
         z_lmm(i,1) = 1.0_dp
         y_lmm(i) = 2.0_dp + 1.5_dp*x_lmm(i,2) + random_effect + 0.15_dp*cos(real(3*g+rr,dp))
      end do
   end do
   ctl%iter_max = 300
   ctl%tolerance = 1.0e-9_dp
   ctl%minprior = 0.0_dp
   call flexmix_lmm(x_lmm, y_lmm, z_lmm, group_lmm, 1, lmm_fit, control=ctl)
   call emit(unit, 'lmm_beta0', lmm_fit%beta(1,1))
   call emit(unit, 'lmm_beta1', lmm_fit%beta(2,1))
   call emit(unit, 'lmm_random_var', lmm_fit%random_covariance(1,1,1))
   call emit(unit, 'lmm_residual_var', lmm_fit%residual_variance(1))
   call emit(unit, 'lmm_loglik', lmm_fit%loglik)

   y_lmmc = y_lmm
   censored_lmm = .false.
   do i = 1, 48
      if (mod(i-1,4) /= 0) cycle
      censored_lmm(i) = .true.
      y_lmmc(i) = y_lmm(i) + 0.22_dp
   end do
   ctl%iter_max = 700
   ctl%tolerance = 1.0e-8_dp
   call flexmix_lmmc(x_lmm, y_lmmc, z_lmm, group_lmm, censored_lmm, 1, lmmc_fit, control=ctl)
   call emit(unit, 'lmmc_beta0', lmmc_fit%beta(1,1))
   call emit(unit, 'lmmc_beta1', lmmc_fit%beta(2,1))
   call emit(unit, 'lmmc_random_var', lmmc_fit%random_covariance(1,1,1))
   call emit(unit, 'lmmc_residual_var', lmmc_fit%residual_variance(1))
   call emit(unit, 'lmmc_loglik', lmmc_fit%loglik)

   ctl%iter_max = 3000
   ctl%tolerance = 1.0e-13_dp
   call flexmix_lmc(x_lmm, y_lmmc, group_lmm, censored_lmm, 1, lmc_fit, control=ctl)
   call emit(unit, 'lmc_beta0', lmc_fit%beta(1,1))
   call emit(unit, 'lmc_beta1', lmc_fit%beta(2,1))
   call emit(unit, 'lmc_residual_var', lmc_fit%residual_variance(1))
   call emit(unit, 'lmc_loglik', lmc_fit%loglik)

   i = 0
   do g = 1, 4
      random_effect = 0.6_dp*sin(0.5_dp*real(g,dp))
      do rr = 1, 3
         i = i + 1
         group_multi_cens(i) = g
         x_multi_cens(i,1) = 1.0_dp
         x_multi_cens(i,2) = real(rr-1,dp) / 2.0_dp
         z_multi_cens(i,1) = 1.0_dp
         y_multi_cens(i) = 1.3_dp + 1.1_dp*x_multi_cens(i,2) + random_effect + &
                            0.08_dp*cos(real(2*g+rr,dp))
         censored_multi(i) = rr <= 2
         if (censored_multi(i)) y_multi_cens(i) = y_multi_cens(i) + 0.15_dp
      end do
   end do
   ctl%iter_max = 12
   ctl%tolerance = 1.0e-12_dp
   call flexmix_lmmc(x_multi_cens, y_multi_cens, z_multi_cens, group_multi_cens, censored_multi, 1, &
                     lmmc_multi_fit, control=ctl)
   call emit(unit, 'lmmc_multi_beta0', lmmc_multi_fit%beta(1,1))
   call emit(unit, 'lmmc_multi_beta1', lmmc_multi_fit%beta(2,1))
   call emit(unit, 'lmmc_multi_random_var', lmmc_multi_fit%random_covariance(1,1,1))
   call emit(unit, 'lmmc_multi_residual_var', lmmc_multi_fit%residual_variance(1))
   call emit(unit, 'lmmc_multi_loglik', lmmc_multi_fit%loglik)

   smooth_penalty = 0.0_dp
   smooth_penalty(3,3) = 1.0_dp
   smooth_lambda = 5.0_dp
   do i = 1, 10
      smooth_t = -1.0_dp + 2.0_dp*real(i-1,dp)/9.0_dp
      x_smooth(i,:) = [1.0_dp,smooth_t,smooth_t*smooth_t]
      y_smooth(i) = 1.0_dp + 2.0_dp*smooth_t + 0.5_dp*smooth_t*smooth_t + 0.05_dp*sin(real(i,dp))
      y_smooth_pois(i) = real(mod(i,4),dp) + 1.0_dp
      smooth_success(i) = real(mod(i,3),dp)
      smooth_failure(i) = 2.0_dp - smooth_success(i)
   end do
   call flexmix_mgcv_gaussian(x_smooth,y_smooth,smooth_penalty,1,smooth_gfit,lambda_grid=smooth_lambda)
   call emit(unit, 'mgcv_gaussian_beta0', smooth_gfit%beta(1,1))
   call emit(unit, 'mgcv_gaussian_beta1', smooth_gfit%beta(2,1))
   call emit(unit, 'mgcv_gaussian_beta2', smooth_gfit%beta(3,1))
   call emit(unit, 'mgcv_gaussian_sigma', smooth_gfit%sigma(1))
   call emit(unit, 'mgcv_gaussian_edf', smooth_gfit%component_effective_df(1))
   call flexmix_mgcv_poisson(x_smooth,y_smooth_pois,smooth_penalty,1,smooth_pfit,lambda_grid=smooth_lambda)
   call emit(unit, 'mgcv_poisson_beta0', smooth_pfit%beta(1,1))
   call emit(unit, 'mgcv_poisson_beta1', smooth_pfit%beta(2,1))
   call emit(unit, 'mgcv_poisson_beta2', smooth_pfit%beta(3,1))
   call flexmix_mgcv_binomial(x_smooth,smooth_success,smooth_failure,smooth_penalty,1,smooth_bfit, &
                             lambda_grid=smooth_lambda)
   call emit(unit, 'mgcv_binomial_beta0', smooth_bfit%beta(1,1))
   call emit(unit, 'mgcv_binomial_beta1', smooth_bfit%beta(2,1))
   call emit(unit, 'mgcv_binomial_beta2', smooth_bfit%beta(3,1))

   close(unit)

contains

   subroutine emit(unit, metric, value)
      integer, intent(in) :: unit !! Connected output unit receiving one CSV metric row.
      character(len=*), intent(in) :: metric !! Metric name written in the first CSV column.
      real(dp), intent(in) :: value !! Numeric metric value written in the second CSV column.
      write(unit,'(a,a,es24.16)') trim(metric), ',', value
   end subroutine emit

end program parity_driver
