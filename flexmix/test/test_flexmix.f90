program test_flexmix
   use flexmix, only : dp, flexmix_result, flexmix_control, flexmix_step_result, flexmix_boot_result
   use flexmix, only : flexmix_gaussian, flexmix_poisson, flexmix_binomial, flexmix_gamma, flexmix_multinomial
   use flexmix, only : flexmix_glmnet_gaussian, flexmix_glmnet_poisson, flexmix_glmnet_binomial
   use flexmix, only : flexmix_lmm, flexmix_lmer, flexmix_lmmc, flexmix_lmc
   use flexmix, only : flexmix_mgcv_gaussian, flexmix_mgcv_poisson, flexmix_mgcv_binomial
   use flexmix, only : flexmix_conditional_logit
   use flexmix, only : flexmix_ziglm_poisson, flexmix_ziglm_binomial, flexmix_robust_gaussian, flexmix_robust_poisson
   use flexmix, only : flexmix_mvnorm, flexmix_factanal, flexmix_mvbinary, flexmix_mvpois, flexmix_mvcombi
   use flexmix, only : flexmix_dist1, flexmix_loglik, flexmix_nobs, flexmix_aic, flexmix_bic
   use flexmix, only : flexmix_cloglik, flexmix_icl, flexmix_eic
   use flexmix, only : flexmix_get_k, flexmix_get_obs, flexmix_get_prior, flexmix_get_posterior
   use flexmix, only : flexmix_get_clusters, flexmix_get_parameters, flexmix_predict_regression, flexmix_predict_multinomial
   use flexmix, only : flexmix_model_poisson, flexmix_role_structural_zero, flexmix_role_robust_background
   use flexmix, only : flexmix_relabel, kl_divergence_matrix, kl_divergence_regression, kl_divergence_mvnorm
   use flexmix, only : step_flexmix_gaussian, step_flexmix_unique, get_model
   use flexmix, only : flexmix_check_result, flexmix_remove_component, flexmix_refit_gaussian
   use flexmix, only : flexmix_simulate_regression, flexmix_make_distribution, ex_linear, ex_npreg, ex_nclus
   use flexmix, only : flexmix_get_design, flexmix_parameter_vector, flexmix_replace_parameters
   use flexmix, only : flexmix_exist_gradient, flexmix_loglik_regression_at, flexmix_gradient_regression
   use flexmix, only : flexmix_refit_optim_regression
   use flexmix, only : flexmix_boot_regression, flexmix_lr_test_regression
   use flexmix, only : flexmix_glmfix_gaussian, flexmix_glmfix_poisson, flexmix_glmfix_binomial, flexmix_glmfix_gamma
   implicit none
   integer :: failures

   failures = 0
   call test_gaussian_regression(failures)
   call test_poisson_and_binomial(failures)
   call test_extended_glm_drivers(failures)
   call test_penalized_and_mixed_drivers(failures)
   call test_smooth_drivers(failures)
   call test_multivariate_models(failures)
   call test_factor_analysis(failures)
   call test_conditional_logit(failures)
   call test_dist1(failures)
   call test_accessors_information_and_kl(failures)
   call test_grouping_and_concomitant(failures)
   call test_step_selection(failures)
   call test_bootstrap_and_unique(failures)
   call test_refit_simulate_and_remove(failures)
   call test_parameter_refit_tools(failures)
   call test_example_generators(failures)
   call test_glmfix(failures)
   if (failures /= 0) then
      write(*,'(a,i0)') 'FAILED tests: ', failures
      error stop 1
   end if
   write(*,'(a)') 'All flexmix tests passed.'

contains

   subroutine test_gaussian_regression(failures)
      integer, intent(inout) :: failures !! Running count of failed assertions updated by this test case.
      integer, parameter :: n = 20
      real(dp) :: x(n,2), y(n)
      type(flexmix_result) :: fit
      type(flexmix_control) :: ctl
      integer :: i, info
      ctl%minprior = 0.01_dp
      ctl%tolerance = 1.0e-10_dp
      do i = 1, n
         x(i,1) = 1.0_dp
         x(i,2) = real(mod(i-1,10),dp) / 9.0_dp
         if (i <= 10) then
            y(i) = 1.0_dp + 2.0_dp * x(i,2) + 0.02_dp * real((-1)**i,dp)
         else
            y(i) = 8.0_dp - 1.5_dp * x(i,2) + 0.02_dp * real((-1)**i,dp)
         end if
      end do
      call flexmix_gaussian(x, y, 2, fit, ctl)
      call assert_true(fit%status == 0, 'Gaussian EM status', failures)
      call assert_true(fit%converged, 'Gaussian EM convergence', failures)
      call assert_true(fit%k == 2, 'Gaussian component count', failures)
      call sort_by_intercept(fit, info)
      call assert_true(info == 0, 'Gaussian relabel', failures)
      call assert_close(fit%prior(1), 0.5_dp, 1.0e-10_dp, 'Gaussian prior 1', failures)
      call assert_close(fit%prior(2), 0.5_dp, 1.0e-10_dp, 'Gaussian prior 2', failures)
      call assert_close(fit%beta(1,1), 0.9945454545454546_dp, 2.0e-10_dp, 'Gaussian intercept 1', failures)
      call assert_close(fit%beta(2,1), 2.0109090909090903_dp, 2.0e-10_dp, 'Gaussian slope 1', failures)
      call assert_close(fit%beta(1,2), 7.9945454545454542_dp, 2.0e-10_dp, 'Gaussian intercept 2', failures)
      call assert_close(fit%beta(2,2), -1.4890909090909117_dp, 2.0e-10_dp, 'Gaussian slope 2', failures)
      call assert_true(all(fit%cluster(1:10) == 1), 'Gaussian first component labels', failures)
      call assert_true(all(fit%cluster(11:20) == 2), 'Gaussian second component labels', failures)
   end subroutine test_gaussian_regression

   subroutine test_poisson_and_binomial(failures)
      integer, intent(inout) :: failures !! Running count of failed assertions updated by this test case.
      real(dp) :: x(5,1), y(5), xb(4,1), success(4), failure(4)
      type(flexmix_result) :: pfit, bfit
      x = 1.0_dp
      y = [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 5.0_dp]
      call flexmix_poisson(x, y, 1, pfit)
      call assert_true(pfit%status == 0, 'Poisson status', failures)
      call assert_close(pfit%beta(1,1), log(3.0_dp), 2.0e-9_dp, 'Poisson intercept MLE', failures)
      xb = 1.0_dp
      success = [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp]
      failure = [4.0_dp, 3.0_dp, 2.0_dp, 1.0_dp]
      call flexmix_binomial(xb, success, failure, 1, bfit)
      call assert_true(bfit%status == 0, 'Binomial status', failures)
      call assert_close(bfit%beta(1,1), 0.0_dp, 2.0e-9_dp, 'Binomial intercept MLE', failures)
   end subroutine test_poisson_and_binomial

   subroutine test_extended_glm_drivers(failures)
      integer, intent(inout) :: failures !! Running count of failed assertions updated by this test case.
      real(dp) :: xg(4,1), yg(4), xp(5,1), yp(5), offsetp(5), xm(10,1), xs(6,1)
      real(dp) :: ys(6), success(6), failure(6), yrob(6)
      real(dp), allocatable :: prediction(:,:), class_probability(:,:,:), parameters(:,:)
      integer :: y_class(10), initial_cluster(6), info
      type(flexmix_result) :: gamma_fit, poisson_fit, multi_fit, zigp_fit, zigb_fit, robustg_fit, robustp_fit

      xg = 1.0_dp
      yg = [1.0_dp, 2.0_dp, 4.0_dp, 8.0_dp]
      call flexmix_gamma(xg, yg, 1, gamma_fit)
      call assert_true(gamma_fit%status == 0, 'Gamma regression status', failures)
      call assert_close(gamma_fit%beta(1,1), 1.0_dp / 3.75_dp, 1.0e-10_dp, &
                        'Gamma inverse-link intercept', failures)
      call assert_close(gamma_fit%shape(1), 1.7728291787272257_dp, 2.0e-10_dp, &
                        'Gamma deviance shape', failures)
      call flexmix_get_parameters(gamma_fit, parameters, info)
      call assert_true(info == 0 .and. size(parameters,1) == 2, 'Gamma parameter packing', failures)

      xp = 1.0_dp
      yp = [1.0_dp, 2.0_dp, 4.0_dp, 8.0_dp, 16.0_dp]
      offsetp = log([1.0_dp, 1.0_dp, 2.0_dp, 2.0_dp, 4.0_dp])
      call flexmix_poisson(xp, yp, 1, poisson_fit, offset=offsetp)
      call assert_true(poisson_fit%status == 0, 'Poisson offset status', failures)
      call assert_close(poisson_fit%beta(1,1), log(3.1_dp), 2.0e-9_dp, 'Poisson offset intercept', failures)

      xm = 1.0_dp
      y_class = [1,1,2,2,2,3,3,3,3,3]
      call flexmix_multinomial(xm, y_class, 1, multi_fit)
      call assert_true(multi_fit%status == 0 .and. multi_fit%nclass == 3, 'Multinomial status/classes', failures)
      call assert_close(multi_fit%multinomial_coef(1,1,1), log(1.5_dp), 2.0e-9_dp, &
                        'Multinomial class-2 coefficient', failures)
      call assert_close(multi_fit%multinomial_coef(1,2,1), log(2.5_dp), 2.0e-9_dp, &
                        'Multinomial class-3 coefficient', failures)
      call flexmix_predict_multinomial(multi_fit, xm(1:1,:), class_probability, info)
      call assert_true(info == 0, 'Multinomial prediction status', failures)
      call assert_close(class_probability(1,1,1), 0.2_dp, 2.0e-10_dp, 'Multinomial probability 1', failures)
      call assert_close(class_probability(1,2,1), 0.3_dp, 2.0e-10_dp, 'Multinomial probability 2', failures)
      call assert_close(class_probability(1,3,1), 0.5_dp, 2.0e-10_dp, 'Multinomial probability 3', failures)

      xs = 1.0_dp
      initial_cluster = [1,1,1,2,2,2]
      ys = [0.0_dp, 0.0_dp, 0.0_dp, 2.0_dp, 3.0_dp, 4.0_dp]
      call flexmix_ziglm_poisson(xs, ys, 2, zigp_fit, initial_cluster=initial_cluster)
      call assert_true(zigp_fit%status == 0 .and. zigp_fit%k == 2, 'ZIGLM Poisson status', failures)
      call assert_true(zigp_fit%component_role(1) == flexmix_role_structural_zero, &
                       'ZIGLM Poisson structural role', failures)
      call flexmix_predict_regression(zigp_fit, xs(1:1,:), prediction, info)
      call assert_true(info == 0, 'ZIGLM Poisson prediction status', failures)
      call assert_close(prediction(1,1), 0.0_dp, 0.0_dp, 'ZIGLM structural prediction', failures)
      call assert_true(zigp_fit%df == 2, 'ZIGLM Poisson degrees of freedom', failures)

      success = [0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp, 2.0_dp, 3.0_dp]
      failure = [4.0_dp, 4.0_dp, 4.0_dp, 3.0_dp, 2.0_dp, 1.0_dp]
      call flexmix_ziglm_binomial(xs, success, failure, 2, zigb_fit, initial_cluster=initial_cluster)
      call assert_true(zigb_fit%status == 0 .and. zigb_fit%component_role(1) == flexmix_role_structural_zero, &
                       'ZIGLM binomial structural component', failures)
      call flexmix_predict_regression(zigb_fit, xs(1:1,:), prediction, info)
      call assert_close(prediction(1,1), 0.0_dp, 0.0_dp, 'ZIGLM binomial structural prediction', failures)

      yrob = [0.0_dp, 1.0_dp, 2.0_dp, 10.0_dp, 11.0_dp, 12.0_dp]
      call flexmix_robust_gaussian(xs, yrob, 2, robustg_fit, initial_cluster=initial_cluster)
      call assert_true(robustg_fit%status == 0 .and. &
                       robustg_fit%component_role(1) == flexmix_role_robust_background, &
                       'Robust Gaussian background role', failures)
      call assert_close(robustg_fit%beta(1,1), 6.0_dp, 1.0e-12_dp, 'Robust Gaussian background mean', failures)
      call assert_close(robustg_fit%sigma(1), sqrt(30.8_dp), 2.0e-10_dp, 'Robust Gaussian background scale', failures)
      call assert_true(robustg_fit%df == 3, 'Robust Gaussian degrees of freedom', failures)

      call flexmix_robust_poisson(xs, yrob, 2, robustp_fit, initial_cluster=initial_cluster)
      call assert_true(robustp_fit%status == 0 .and. &
                       robustp_fit%component_role(1) == flexmix_role_robust_background, &
                       'Robust Poisson background role', failures)
      call assert_close(robustp_fit%beta(1,1), log(18.0_dp), 2.0e-10_dp, &
                        'Robust Poisson background intercept', failures)
      call assert_true(robustp_fit%df == 2, 'Robust Poisson degrees of freedom', failures)

      call flexmix_remove_component(robustp_fit, 1, info)
      call assert_true(info == 0 .and. robustp_fit%model_kind == flexmix_model_poisson, &
                       'Robust component removal converts family', failures)
   end subroutine test_extended_glm_drivers

   subroutine test_penalized_and_mixed_drivers(failures)
      integer, intent(inout) :: failures !! Running count of failed assertions updated by this test case.
      integer, parameter :: npen = 20, ng = 12, m = 4, nmix = ng*m
      real(dp) :: xp(npen,4), yg(npen), yp(npen), success(npen), failure(npen), lambda_grid(1)
      real(dp) :: xm(nmix,2), zm(nmix,1), ym(nmix), yc(nmix), prediction(1,1)
      integer :: group(nmix), i, g, r, info
      logical :: censored(nmix), no_censor(nmix)
      type(flexmix_result) :: gfit, pfit, bfit, lmmfit, lmerfit, lmmcfit, lmmc0fit, lmcfit
      type(flexmix_control) :: ctl
      real(dp), allocatable :: pred(:,:)
      real(dp) :: random_effect

      do i = 1, npen
         xp(i,1) = 1.0_dp
         xp(i,2) = (real(i,dp) - 10.5_dp) / 5.0_dp
         xp(i,3) = sin(real(i,dp))
         xp(i,4) = cos(0.3_dp * real(i,dp))
         yg(i) = 2.0_dp + 3.0_dp*xp(i,2) + 0.02_dp*sin(2.0_dp*real(i,dp))
         yp(i) = max(0.0_dp,real(nint(exp(0.2_dp + 0.4_dp*xp(i,2))),dp))
         success(i) = real(mod(i,5),dp)
         failure(i) = 4.0_dp - success(i)
      end do
      lambda_grid = 0.05_dp
      ctl%iter_max = 20
      ctl%minprior = 0.0_dp
      ctl%tolerance = 1.0e-10_dp
      call flexmix_glmnet_gaussian(xp, yg, 1, gfit, control=ctl, adaptive=.false., alpha=1.0_dp, &
                                   lambda_grid=lambda_grid, nfolds=2)
      call assert_true(gfit%status == 0 .and. allocated(gfit%penalty_lambda), 'FLXMRglmnet Gaussian status', failures)
      call assert_close(gfit%penalty_lambda(1), 0.05_dp, 0.0_dp, 'FLXMRglmnet selected lambda', failures)
      call assert_close(gfit%beta(1,1), 2.0009077208905435_dp, 2.0e-8_dp, 'FLXMRglmnet Gaussian intercept', failures)
      call assert_close(gfit%beta(2,1), 2.96277677_dp, 1.0e-6_dp, 'FLXMRglmnet Gaussian selected slope', failures)
      call assert_close(gfit%beta(3,1), 0.0_dp, 1.0e-10_dp, 'FLXMRglmnet Gaussian zero term 1', failures)
      call assert_close(gfit%beta(4,1), 0.0_dp, 1.0e-10_dp, 'FLXMRglmnet Gaussian zero term 2', failures)

      call flexmix_glmnet_poisson(xp, yp, 1, pfit, control=ctl, adaptive=.false., alpha=1.0_dp, &
                                  lambda_grid=lambda_grid, nfolds=2)
      call assert_true(pfit%status == 0, 'FLXMRglmnet Poisson status', failures)
      call assert_close(pfit%beta(1,1), 0.2910207155_dp, 2.0e-6_dp, 'FLXMRglmnet Poisson intercept', failures)
      call flexmix_glmnet_binomial(xp, success, failure, 1, bfit, control=ctl, adaptive=.false., alpha=1.0_dp, &
                                   lambda_grid=lambda_grid, nfolds=2)
      call assert_true(bfit%status == 0, 'FLXMRglmnet binomial status', failures)
      call assert_close(bfit%beta(3,1), 0.19529733_dp, 2.0e-6_dp, 'FLXMRglmnet binomial selected coefficient', failures)

      i = 0
      do g = 1, ng
         random_effect = 0.8_dp * sin(0.7_dp*real(g,dp))
         do r = 1, m
            i = i + 1
            group(i) = g
            xm(i,1) = 1.0_dp
            xm(i,2) = real(r-1,dp) / 3.0_dp
            zm(i,1) = 1.0_dp
            ym(i) = 2.0_dp + 1.5_dp*xm(i,2) + random_effect + 0.15_dp*cos(real(3*g+r,dp))
         end do
      end do
      ctl%iter_max = 300
      ctl%tolerance = 1.0e-9_dp
      call flexmix_lmm(xm, ym, zm, group, 1, lmmfit, control=ctl)
      call assert_true(lmmfit%status == 0 .and. lmmfit%converged, 'FLXMRlmm convergence', failures)
      call assert_close(lmmfit%beta(1,1), 2.1717724871_dp, 2.0e-8_dp, 'FLXMRlmm intercept', failures)
      call assert_close(lmmfit%beta(2,1), 1.4998264762_dp, 2.0e-8_dp, 'FLXMRlmm slope', failures)
      call assert_close(lmmfit%random_covariance(1,1,1), 0.3250235869_dp, 2.0e-8_dp, 'FLXMRlmm random variance', failures)
      call assert_close(lmmfit%residual_variance(1), 0.01557732169_dp, 2.0e-9_dp, 'FLXMRlmm residual variance', failures)
      call flexmix_predict_regression(lmmfit, reshape([1.0_dp,0.5_dp],[1,2]), pred, info)
      call assert_true(info == 0, 'FLXMRlmm fixed-effect prediction status', failures)
      prediction = pred
      call assert_close(prediction(1,1), lmmfit%beta(1,1) + 0.5_dp*lmmfit%beta(2,1), 1.0e-12_dp, &
                        'FLXMRlmm fixed-effect prediction', failures)
      call flexmix_lmer(xm, ym, zm, group, 1, lmerfit, control=ctl)
      call assert_true(lmerfit%status == 0, 'FLXMRlmer status', failures)
      call assert_close(lmerfit%loglik, lmmfit%loglik, 1.0e-10_dp, 'FLXMRlmer marginal likelihood alias', failures)

      no_censor = .false.
      call flexmix_lmmc(xm, ym, zm, group, no_censor, 1, lmmc0fit, control=ctl)
      call assert_true(lmmc0fit%status == 0 .and. lmmc0fit%converged, 'FLXMRlmmc uncensored status', failures)
      call assert_close(lmmc0fit%beta(1,1), lmmfit%beta(1,1), 2.0e-7_dp, &
                        'FLXMRlmmc uncensored intercept parity', failures)
      call assert_close(lmmc0fit%beta(2,1), lmmfit%beta(2,1), 2.0e-10_dp, &
                        'FLXMRlmmc uncensored slope parity', failures)
      call assert_close(lmmc0fit%random_covariance(1,1,1), lmmfit%random_covariance(1,1,1), 2.0e-7_dp, &
                        'FLXMRlmmc uncensored random variance parity', failures)
      call assert_close(lmmc0fit%loglik, lmmfit%loglik, 2.0e-9_dp, &
                        'FLXMRlmmc uncensored likelihood parity', failures)

      yc = ym
      censored = .false.
      do i = 1, nmix
         if (mod(i-1,m) /= 0) cycle
         censored(i) = .true.
         yc(i) = ym(i) + 0.22_dp
      end do
      ctl%iter_max = 700
      ctl%tolerance = 1.0e-8_dp
      call flexmix_lmmc(xm, yc, zm, group, censored, 1, lmmcfit, control=ctl)
      call assert_true(lmmcfit%status == 0 .and. lmmcfit%converged, 'FLXMRlmmc censored convergence', failures)
      call assert_close(lmmcfit%beta(1,1), 2.13142006_dp, 8.0e-4_dp, 'FLXMRlmmc censored intercept', failures)
      call assert_close(lmmcfit%beta(2,1), 1.55170806_dp, 2.0e-5_dp, 'FLXMRlmmc censored slope', failures)
      call assert_close(lmmcfit%random_covariance(1,1,1), 0.32832676_dp, 2.0e-5_dp, &
                        'FLXMRlmmc censored random variance', failures)
      call assert_close(lmmcfit%residual_variance(1), 0.01072646_dp, 2.0e-6_dp, &
                        'FLXMRlmmc censored residual variance', failures)
      call assert_close(lmmcfit%loglik, 0.59033647_dp, 2.0e-5_dp, 'FLXMRlmmc censored likelihood', failures)

      ctl%iter_max = 3000
      ctl%tolerance = 1.0e-13_dp
      call flexmix_lmc(xm, yc, group, censored, 1, lmcfit, control=ctl)
      call assert_true(lmcfit%status == 0 .and. lmcfit%converged, 'FLXMRlmc censored convergence', failures)
      call assert_close(lmcfit%beta(1,1), 1.67736853_dp, 2.0e-6_dp, 'FLXMRlmc censored intercept', failures)
      call assert_close(lmcfit%beta(2,1), 2.13548870_dp, 2.0e-6_dp, 'FLXMRlmc censored slope', failures)
      call assert_close(lmcfit%residual_variance(1), 0.368141256_dp, 2.0e-7_dp, &
                        'FLXMRlmc censored residual variance', failures)
      call assert_close(lmcfit%loglik, -37.0240551054_dp, 2.0e-8_dp, 'FLXMRlmc censored likelihood', failures)
   end subroutine test_penalized_and_mixed_drivers

   subroutine test_smooth_drivers(failures)
      integer, intent(inout) :: failures !! Running count of failed assertions updated by this test case.
      integer, parameter :: n = 10
      real(dp) :: x(n,3), yg(n), yp(n), success(n), failure(n), penalty(3,3), lambda_grid(1), t
      type(flexmix_result) :: gfit, pfit, bfit
      integer :: i

      penalty = 0.0_dp
      penalty(3,3) = 1.0_dp
      lambda_grid = 5.0_dp
      do i = 1, n
         t = -1.0_dp + 2.0_dp*real(i-1,dp)/real(n-1,dp)
         x(i,:) = [1.0_dp,t,t*t]
         yg(i) = 1.0_dp + 2.0_dp*t + 0.5_dp*t*t + 0.05_dp*sin(real(i,dp))
         yp(i) = real(mod(i,4),dp) + 1.0_dp
         success(i) = real(mod(i,3),dp)
         failure(i) = 2.0_dp - success(i)
      end do
      call flexmix_mgcv_gaussian(x,yg,penalty,1,gfit,lambda_grid=lambda_grid)
      call assert_true(gfit%status == 0 .and. gfit%converged, 'FLXMRmgcv Gaussian status', failures)
      call assert_close(gfit%beta(1,1), 1.1663084047886485_dp, 2.0e-11_dp, &
                        'FLXMRmgcv Gaussian intercept', failures)
      call assert_close(gfit%beta(2,1), 1.9907441624531801_dp, 2.0e-11_dp, &
                        'FLXMRmgcv Gaussian linear term', failures)
      call assert_close(gfit%beta(3,1), 0.10910759098371964_dp, 2.0e-11_dp, &
                        'FLXMRmgcv Gaussian penalized term', failures)
      call assert_close(gfit%sigma(1), 0.17597916333946315_dp, 2.0e-11_dp, &
                        'FLXMRmgcv Gaussian residual scale', failures)
      call assert_close(gfit%component_effective_df(1), 3.2047851065377064_dp, 2.0e-11_dp, &
                        'FLXMRmgcv Gaussian effective df', failures)
      call assert_close(gfit%smoothing_lambda(1), 5.0_dp, 0.0_dp, 'FLXMRmgcv Gaussian lambda', failures)
      call assert_close(flexmix_aic(gfit), -2.0_dp*gfit%loglik + 2.0_dp*gfit%effective_df, 2.0e-12_dp, &
                        'FLXMRmgcv effective-df AIC', failures)
      call assert_close(flexmix_bic(gfit), &
                        -2.0_dp*gfit%loglik + gfit%effective_df*log(real(n,dp)), 2.0e-12_dp, &
                        'FLXMRmgcv effective-df BIC', failures)

      call flexmix_mgcv_poisson(x,yp,penalty,1,pfit,lambda_grid=lambda_grid)
      call assert_true(pfit%status == 0 .and. pfit%converged, 'FLXMRmgcv Poisson status', failures)
      call assert_close(pfit%beta(1,1), 0.9161579684057697_dp, 2.0e-9_dp, &
                        'FLXMRmgcv Poisson intercept', failures)
      call assert_close(pfit%beta(2,1), -0.03273235717192185_dp, 2.0e-9_dp, &
                        'FLXMRmgcv Poisson linear term', failures)
      call assert_close(pfit%beta(3,1), -0.00020977722351206652_dp, 2.0e-9_dp, &
                        'FLXMRmgcv Poisson penalized term', failures)

      call flexmix_mgcv_binomial(x,success,failure,penalty,1,bfit,lambda_grid=lambda_grid)
      call assert_true(bfit%status == 0 .and. bfit%converged, 'FLXMRmgcv binomial status', failures)
      call assert_close(bfit%beta(1,1), 0.0_dp, 2.0e-9_dp, 'FLXMRmgcv binomial intercept', failures)
      call assert_close(bfit%beta(2,1), -0.3294079895020593_dp, 2.0e-9_dp, &
                        'FLXMRmgcv binomial linear term', failures)
   end subroutine test_smooth_drivers

   subroutine test_multivariate_models(failures)
      integer, intent(inout) :: failures !! Running count of failed assertions updated by this test case.
      real(dp) :: yn(12,2), yb(6,2), yp(5,2), yc(6,2)
      logical :: binary(2)
      type(flexmix_result) :: nfit, bfit, pfit, cfit
      integer :: i, info
      do i = 1, 6
         yn(i,1) = -4.0_dp + 0.1_dp * real(i-3,dp)
         yn(i,2) = -1.0_dp + 0.05_dp * real((-1)**i,dp)
         yn(i+6,1) = 5.0_dp + 0.1_dp * real(i-3,dp)
         yn(i+6,2) = 3.0_dp + 0.05_dp * real((-1)**i,dp)
      end do
      call flexmix_mvnorm(yn, 2, nfit, diagonal=.true.)
      call assert_true(nfit%status == 0, 'MV normal status', failures)
      call assert_true(nfit%k == 2, 'MV normal component count', failures)
      call sort_by_center(nfit, info)
      call assert_true(info == 0, 'MV normal relabel', failures)
      call assert_close(nfit%center(1,1), -3.95_dp, 1.0e-10_dp, 'MV normal center 1', failures)
      call assert_close(nfit%center(1,2), 5.05_dp, 1.0e-10_dp, 'MV normal center 2', failures)
      yb = reshape([0.0_dp,0.0_dp,0.0_dp,1.0_dp,1.0_dp,1.0_dp, &
                    0.0_dp,1.0_dp,0.0_dp,1.0_dp,1.0_dp,1.0_dp], [6,2])
      call flexmix_mvbinary(yb, 1, bfit)
      call assert_close(bfit%probability(1,1), 0.5_dp, 1.0e-12_dp, 'MV binary probability 1', failures)
      call assert_close(bfit%probability(2,1), 4.0_dp/6.0_dp, 1.0e-12_dp, 'MV binary probability 2', failures)
      yp = reshape([1.0_dp,2.0_dp,3.0_dp,4.0_dp,5.0_dp, 2.0_dp,4.0_dp,6.0_dp,8.0_dp,10.0_dp], [5,2])
      call flexmix_mvpois(yp, 1, pfit)
      call assert_close(pfit%lambda(1,1), 3.0_dp, 1.0e-12_dp, 'MV Poisson mean 1', failures)
      call assert_close(pfit%lambda(2,1), 6.0_dp, 1.0e-12_dp, 'MV Poisson mean 2', failures)
      yc(:,1) = [0.0_dp,0.0_dp,1.0_dp,1.0_dp,1.0_dp,0.0_dp]
      yc(:,2) = [1.0_dp,2.0_dp,3.0_dp,4.0_dp,5.0_dp,6.0_dp]
      binary = [.true., .false.]
      call flexmix_mvcombi(yc, binary, 1, cfit)
      call assert_true(cfit%status == 0, 'MV combined status', failures)
      call assert_close(cfit%center(1,1), 0.5_dp, 1.0e-12_dp, 'MV combined Bernoulli center', failures)
      call assert_close(cfit%center(2,1), 3.5_dp, 1.0e-12_dp, 'MV combined Gaussian center', failures)
   end subroutine test_multivariate_models

   subroutine test_factor_analysis(failures)
      integer, intent(inout) :: failures !! Running count of failed assertions updated by this test case.
      integer, parameter :: n = 60, d = 4
      real(dp) :: y(n,d), ymix(80,d), z1, z2
      integer :: initial_cluster(80), i, info
      type(flexmix_result) :: fit, mixture
      type(flexmix_control) :: ctl

      do i = 1, n
         z1 = sin(0.37_dp*real(i,dp)) + 0.3_dp*cos(0.11_dp*real(i,dp))
         z2 = cos(0.29_dp*real(i,dp))
         y(i,1) = 2.0_dp + z1 + 0.18_dp*sin(1.7_dp*real(i,dp))
         y(i,2) = -1.0_dp + 0.85_dp*z1 + 0.25_dp*cos(1.3_dp*real(i,dp))
         y(i,3) = 0.5_dp + 0.65_dp*z1 + 0.35_dp*sin(1.1_dp*real(i,dp))
         y(i,4) = 3.0_dp + 0.45_dp*z1 + 0.55_dp*z2
      end do
      ctl%minprior = 0.001_dp
      ctl%tolerance = 1.0e-10_dp
      call flexmix_factanal(y, 1, 1, fit, ctl)
      call assert_true(fit%status == 0 .and. fit%factors == 1, 'factor analyzer one-component status', failures)
      call assert_true(all(shape(fit%factor_loadings) == [d,1,1]), 'factor analyzer loading shape', failures)
      call assert_close(fit%uniqueness(1,1), 0.02588148_dp, 2.0e-6_dp, 'factor analyzer uniqueness 1', failures)
      call assert_close(fit%uniqueness(2,1), 0.06880303_dp, 2.0e-6_dp, 'factor analyzer uniqueness 2', failures)
      call assert_close(fit%uniqueness(3,1), 0.19983414_dp, 2.0e-6_dp, 'factor analyzer uniqueness 3', failures)
      call assert_close(fit%uniqueness(4,1), 0.45035140_dp, 2.0e-6_dp, 'factor analyzer uniqueness 4', failures)
      call assert_true(fit%df == 12, 'factor analyzer upstream component df', failures)

      ymix(1:40,:) = y(1:40,:)
      ymix(41:80,:) = y(1:40,:)
      ymix(41:80,1) = ymix(41:80,1) + 9.0_dp
      ymix(41:80,2) = ymix(41:80,2) + 7.0_dp
      initial_cluster(1:40) = 1
      initial_cluster(41:80) = 2
      call flexmix_factanal(ymix, 2, 1, mixture, ctl, initial_cluster=initial_cluster)
      call assert_true(mixture%status == 0 .and. mixture%k == 2, 'factor analyzer mixture status', failures)
      call sort_by_center(mixture, info)
      call assert_true(info == 0, 'factor analyzer mixture relabel', failures)
      call assert_true(all(mixture%cluster(1:40) == 1) .and. all(mixture%cluster(41:80) == 2), &
                       'factor analyzer mixture labels', failures)
   end subroutine test_factor_analysis

   subroutine test_conditional_logit(failures)
      integer, intent(inout) :: failures !! Running count of failed assertions updated by this test case.
      integer, parameter :: nalt = 3, nstrata = 12, n = nalt*nstrata
      integer, parameter :: nmix_strata = 24, nmix = nalt*nmix_strata
      real(dp) :: x(n,1), y(n), xmix(nmix,1), ymix(nmix)
      integer :: strata(n), strata_mix(nmix), initial_cluster(nmix)
      integer :: i, r, selected
      type(flexmix_result) :: fit, mixture
      type(flexmix_control) :: ctl

      do i = 1, nstrata
         selected = 1
         if (i > 3 .and. i <= 7) selected = 2
         if (i > 7) selected = 3
         do r = 1, nalt
            strata((i-1)*nalt+r) = i
            x((i-1)*nalt+r,1) = real(r-2,dp)
            y((i-1)*nalt+r) = merge(1.0_dp,0.0_dp,r == selected)
         end do
      end do
      call flexmix_conditional_logit(x, y, strata, 1, fit)
      call assert_true(fit%status == 0 .and. fit%k == 1, 'conditional logit one-component status', failures)
      call assert_close(fit%beta(1,1), 0.2526512522141187_dp, 2.0e-9_dp, &
                        'conditional logit coefficient', failures)
      call assert_true(fit%df == 1, 'conditional logit one-component df', failures)

      do i = 1, nmix_strata
         if (i <= 12) then
            if (i <= 7) then
               selected = 1
            else if (i <= 11) then
               selected = 2
            else
               selected = 3
            end if
         else
            if (i <= 13) then
               selected = 1
            else if (i <= 17) then
               selected = 2
            else
               selected = 3
            end if
         end if
         do r = 1, nalt
            strata_mix((i-1)*nalt+r) = i
            xmix((i-1)*nalt+r,1) = real(r-2,dp)
            ymix((i-1)*nalt+r) = merge(1.0_dp,0.0_dp,r == selected)
            initial_cluster((i-1)*nalt+r) = merge(1,2,i <= 12)
         end do
      end do
      ctl%minprior = 0.001_dp
      ctl%tolerance = 1.0e-6_dp
      call flexmix_conditional_logit(xmix, ymix, strata_mix, 2, mixture, ctl, initial_cluster=initial_cluster)
      call assert_true(mixture%status == 0 .and. mixture%k == 2, 'conditional logit mixture status', failures)
      call assert_true(mixture%beta(1,1) < 0.0_dp .and. mixture%beta(1,2) > 0.0_dp, &
                       'conditional logit mixture coefficient signs', failures)
      call assert_true(mixture%df == 3, 'conditional logit mixture df', failures)
   end subroutine test_conditional_logit

   subroutine test_dist1(failures)
      integer, intent(inout) :: failures !! Running count of failed assertions updated by this test case.
      real(dp) :: y(5), yl(4)
      type(flexmix_result) :: fit
      y = [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 5.0_dp]
      call flexmix_dist1(y, 1, 'exp', fit)
      call assert_true(fit%status == 0, 'Exponential mixture status', failures)
      call assert_close(fit%rate(1), 1.0_dp/3.0_dp, 1.0e-12_dp, 'Exponential rate', failures)
      yl = exp([0.0_dp, 1.0_dp, 2.0_dp, 3.0_dp])
      call flexmix_dist1(yl, 1, 'lnorm', fit)
      call assert_close(fit%center(1,1), 1.5_dp, 1.0e-12_dp, 'Lognormal meanlog', failures)
      call assert_close(fit%sigma(1), sqrt(1.25_dp), 1.0e-12_dp, 'Lognormal sdlog', failures)
      call flexmix_dist1(y, 1, 'gamma', fit)
      call assert_true(fit%shape(1) > 0.0_dp .and. fit%rate(1) > 0.0_dp, 'Gamma positive parameters', failures)
      call flexmix_dist1(y, 1, 'weibull', fit)
      call assert_true(fit%shape(1) > 0.0_dp .and. fit%scale(1) > 0.0_dp, 'Weibull positive parameters', failures)
      call flexmix_dist1(y, 1, 'invgauss', fit)
      call assert_true(fit%center(1,1) > 0.0_dp .and. fit%lambda(1,1) > 0.0_dp, 'Inverse Gaussian parameters', failures)
   end subroutine test_dist1

   subroutine test_accessors_information_and_kl(failures)
      integer, intent(inout) :: failures !! Running count of failed assertions updated by this test case.
      real(dp) :: x(20,2), y(20), object(3,2), div(2,2)
      real(dp), allocatable :: prior(:), post(:,:), pars(:,:), pred(:,:), kld(:,:)
      integer, allocatable :: cluster(:)
      type(flexmix_result) :: fit, nfit
      type(flexmix_control) :: ctl
      real(dp) :: yn(12,2)
      integer :: i, info
      ctl%minprior = 0.01_dp
      do i = 1, 20
         x(i,1) = 1.0_dp
         x(i,2) = real(mod(i-1,10),dp) / 9.0_dp
         if (i <= 10) then
            y(i) = 1.0_dp + 2.0_dp*x(i,2) + 0.02_dp*real((-1)**i,dp)
         else
            y(i) = 8.0_dp - 1.5_dp*x(i,2) + 0.02_dp*real((-1)**i,dp)
         end if
      end do
      call flexmix_gaussian(x,y,2,fit,ctl)
      call sort_by_intercept(fit, info)
      call flexmix_get_prior(fit, prior)
      call flexmix_get_posterior(fit, post)
      call flexmix_get_clusters(fit, cluster)
      call flexmix_get_parameters(fit, pars, info)
      call flexmix_predict_regression(fit, x, pred, info)
      call assert_true(flexmix_get_k(fit) == 2, 'get K', failures)
      call assert_true(flexmix_get_obs(fit) == 20, 'get obs', failures)
      call assert_close(flexmix_nobs(fit), 20.0_dp, 1.0e-12_dp, 'nobs', failures)
      call assert_close(flexmix_loglik(fit), fit%loglik, 0.0_dp, 'logLik accessor', failures)
      call assert_close(flexmix_aic(fit), -2.0_dp*fit%loglik + 2.0_dp*real(fit%df,dp), 1.0e-12_dp, 'AIC', failures)
      call assert_close(flexmix_bic(fit), -2.0_dp*fit%loglik + real(fit%df,dp)*log(20.0_dp), 1.0e-12_dp, 'BIC', failures)
      call assert_true(flexmix_icl(fit) >= flexmix_bic(fit) - 1.0e-8_dp, 'ICL finite classification penalty', failures)
      call assert_true(flexmix_cloglik(fit) <= fit%loglik + 1.0e-8_dp, 'classification logLik', failures)
      call assert_true(flexmix_eic(fit) >= -1.0e-10_dp .and. flexmix_eic(fit) <= 1.0_dp + 1.0e-10_dp, 'EIC range', failures)
      call assert_true(size(prior) == 2 .and. size(post,2) == 2 .and. size(cluster) == 20, 'accessor shapes', failures)
      call assert_true(size(pars,1) == 3 .and. size(pred,2) == 2, 'parameter/predict shapes', failures)
      object = reshape([0.8_dp,0.1_dp,0.1_dp, 0.1_dp,0.2_dp,0.7_dp], [3,2])
      call kl_divergence_matrix(object, div)
      call assert_true(div(1,2) > 0.0_dp .and. div(2,1) > 0.0_dp, 'discrete KL positive', failures)
      call assert_close(div(1,1), 0.0_dp, 0.0_dp, 'discrete KL diagonal', failures)
      call kl_divergence_regression(fit, x, kld, info)
      call assert_true(info == 0 .and. kld(1,2) > 0.0_dp .and. kld(2,1) > 0.0_dp, 'regression KL', failures)
      do i = 1, 6
         yn(i,1) = -4.0_dp + 0.1_dp * real(i-3,dp)
         yn(i,2) = -1.0_dp + 0.05_dp * real((-1)**i,dp)
         yn(i+6,1) = 5.0_dp + 0.1_dp * real(i-3,dp)
         yn(i+6,2) = 3.0_dp + 0.05_dp * real((-1)**i,dp)
      end do
      call flexmix_mvnorm(yn, 2, nfit, diagonal=.true.)
      call kl_divergence_mvnorm(nfit, kld, info)
      call assert_true(info == 0 .and. kld(1,2) > 0.0_dp .and. kld(2,1) > 0.0_dp, 'MV normal KL', failures)
   end subroutine test_accessors_information_and_kl

   subroutine test_grouping_and_concomitant(failures)
      integer, intent(inout) :: failures !! Running count of failed assertions updated by this test case.
      real(dp) :: x(12,1), y(12), concomitant(12,2)
      integer :: group(12), i
      type(flexmix_result) :: fit
      type(flexmix_control) :: ctl
      ctl%minprior = 0.001_dp
      ctl%tolerance = 1.0e-8_dp
      do i = 1, 12
         x(i,1) = 1.0_dp
         group(i) = 1 + (i-1)/2
         concomitant(i,1) = 1.0_dp
         concomitant(i,2) = -1.0_dp + 2.0_dp*real(group(i)-1,dp)/5.0_dp
         if (group(i) <= 3) then
            y(i) = -2.0_dp + 0.05_dp*real((-1)**i,dp)
         else
            y(i) = 3.0_dp + 0.05_dp*real((-1)**i,dp)
         end if
      end do
      call flexmix_gaussian(x, y, 2, fit, ctl, group=group, concomitant_x=concomitant)
      call assert_true(fit%status == 0, 'group/concomitant status', failures)
      call assert_true(allocated(fit%concomitant_coef), 'concomitant coefficients stored', failures)
      do i = 1, 12, 2
         call assert_close(fit%posterior(i,1), fit%posterior(i+1,1), 1.0e-12_dp, 'group posterior equality', failures)
      end do
   end subroutine test_grouping_and_concomitant

   subroutine test_step_selection(failures)
      integer, intent(inout) :: failures !! Running count of failed assertions updated by this test case.
      real(dp) :: x(20,2), y(20)
      integer :: k_values(3), i, info
      type(flexmix_step_result) :: step
      type(flexmix_result) :: selected
      type(flexmix_control) :: ctl
      ctl%minprior = 0.01_dp
      do i = 1, 20
         x(i,1) = 1.0_dp
         x(i,2) = real(mod(i-1,10),dp) / 9.0_dp
         if (i <= 10) then
            y(i) = 1.0_dp + 2.0_dp*x(i,2) + 0.02_dp*real((-1)**i,dp)
         else
            y(i) = 8.0_dp - 1.5_dp*x(i,2) + 0.02_dp*real((-1)**i,dp)
         end if
      end do
      k_values = [1,2,3]
      call step_flexmix_gaussian(x, y, k_values, 2, step, ctl)
      call assert_true(all(shape(step%logliks) == [3,2]), 'step loglik shape', failures)
      call get_model(step, 'BIC', selected, info)
      call assert_true(info == 0, 'getModel status', failures)
      call assert_true(selected%k >= 1 .and. selected%k <= 3, 'getModel component range', failures)
      call get_model(step, 'logLik', selected, info)
      call assert_true(info == 0 .and. selected%loglik >= maxval(step%logliks) - 1.0e-8_dp, 'getModel logLik', failures)
   end subroutine test_step_selection


   subroutine test_refit_simulate_and_remove(failures)
      integer, intent(inout) :: failures !! Running count of failed assertions updated by this test case.
      integer, parameter :: n = 20
      real(dp) :: x(n,2), y(n), ysim(n), beta(2,2)
      real(dp), allocatable :: xgen(:,:), ygen(:)
      integer, allocatable :: cgen(:)
      integer :: csim(n), ncomp(2), i, info
      logical :: valid
      type(flexmix_result) :: fit, refit
      type(flexmix_control) :: ctl
      ctl%minprior = 0.01_dp
      do i = 1, n
         x(i,1) = 1.0_dp
         x(i,2) = real(mod(i-1,10),dp) / 9.0_dp
         if (i <= 10) then
            y(i) = 1.0_dp + 2.0_dp*x(i,2) + 0.02_dp*real((-1)**i,dp)
         else
            y(i) = 8.0_dp - 1.5_dp*x(i,2) + 0.02_dp*real((-1)**i,dp)
         end if
      end do
      call flexmix_gaussian(x,y,2,fit,ctl)
      call flexmix_check_result(fit, valid)
      call assert_true(valid, 'fit structure check', failures)
      call flexmix_refit_gaussian(fit, x, y, refit, ctl)
      call assert_true(refit%status == 0 .and. refit%k == 2, 'Gaussian refit', failures)
      call set_test_seed()
      call flexmix_simulate_regression(fit, x, ysim, csim, info)
      call assert_true(info == 0 .and. all(csim >= 1) .and. all(csim <= 2), 'regression simulation', failures)
      call flexmix_remove_component(fit, 2, info)
      call assert_true(info == 0 .and. fit%k == 1 .and. size(fit%prior) == 1, 'remove component', failures)
      call flexmix_check_result(fit, valid)
      call assert_true(valid, 'removed fit structure check', failures)
      beta(:,1) = [1.0_dp, 2.0_dp]
      beta(:,2) = [4.0_dp, -1.0_dp]
      ncomp = [5,7]
      call ex_linear(beta, ncomp, 'gaussian', xgen, ygen, cgen, info)
      call assert_true(info == 0 .and. size(xgen,1) == 12 .and. size(ygen) == 12, 'ExLinear generator shape', failures)
      call assert_true(count(cgen == 1) == 5 .and. count(cgen == 2) == 7, 'ExLinear component labels', failures)
   end subroutine test_refit_simulate_and_remove

   subroutine test_bootstrap_and_unique(failures)
      integer, intent(inout) :: failures !! Running count of failed assertions updated by this test case.
      integer, parameter :: n = 24
      real(dp) :: x(n,2), y(n), statistic, p_value
      real(dp), allocatable :: lr_statistics(:)
      integer :: k_values(2), i, info
      type(flexmix_result) :: null_fit
      type(flexmix_control) :: ctl
      type(flexmix_boot_result) :: boot, boot_lr
      type(flexmix_step_result) :: step, unique_step

      do i = 1, n
         x(i,1) = 1.0_dp
         x(i,2) = real(i-1,dp) / real(n-1,dp)
         y(i) = 2.0_dp + 0.75_dp * x(i,2) + 0.15_dp * sin(real(i,dp))
      end do
      ctl%minprior = 0.001_dp
      ctl%tolerance = 1.0e-9_dp
      call flexmix_gaussian(x, y, 1, null_fit, ctl)
      call set_test_seed()
      call flexmix_boot_regression(null_fit, x, y, 4, [1], 'parametric', boot, info, ctl, initialize_solution=.true.)
      call assert_true(info == 0 .and. all(shape(boot%loglik) == [4,1]), 'bootstrap result shape', failures)
      call assert_close(boot%loglik(1,1), null_fit%loglik, 1.0e-12_dp, 'bootstrap observed row', failures)
      call assert_true(all(boot%fitted_k(:,1) == 1), 'bootstrap retained component count', failures)

      call set_test_seed()
      call flexmix_lr_test_regression(null_fit, x, y, 5, 'greater', statistic, p_value, lr_statistics, boot_lr, info, ctl)
      call assert_true(info == 0, 'bootstrap LR test status', failures)
      call assert_true(size(lr_statistics) >= 1, 'bootstrap LR statistic count', failures)
      if (size(lr_statistics) >= 1) then
         call assert_close(lr_statistics(1), statistic, 1.0e-12_dp, 'bootstrap LR statistic retention', failures)
      end if
      call assert_true(p_value >= 0.0_dp .and. p_value <= 1.0_dp, 'bootstrap LR p-value range', failures)

      k_values = [1,1]
      call step_flexmix_gaussian(x, y, k_values, 2, step, ctl)
      call step_flexmix_unique(step, unique_step)
      call assert_true(size(unique_step%models) == 1 .and. unique_step%k(1) == 1, 'step unique retained count', failures)
      call assert_close(unique_step%models(1)%loglik, maxval([(step%models(i)%loglik,i=1,size(step%models))]), &
                        1.0e-12_dp, 'step unique best likelihood', failures)
   end subroutine test_bootstrap_and_unique

   subroutine test_parameter_refit_tools(failures)
      integer, intent(inout) :: failures !! Running count of failed assertions updated by this test case.
      integer, parameter :: n = 20
      real(dp) :: x(n,2), y(n), loglik0, loglik_plus, loglik_minus, h
      real(dp) :: beta_fixed(2,2), prior_fixed(2), sigma_fixed(2)
      real(dp), allocatable :: parameters(:), perturbed(:), gradient(:), ysim(:), estimate(:), se(:), vcov(:,:)
      logical, allocatable :: design(:,:)
      integer, allocatable :: csim(:)
      type(flexmix_result) :: fit, changed, distribution
      type(flexmix_control) :: ctl
      integer :: i, info

      ctl%minprior = 0.01_dp
      ctl%tolerance = 1.0e-11_dp
      do i = 1, n
         x(i,1) = 1.0_dp
         x(i,2) = real(mod(i-1,10),dp) / 9.0_dp
         if (i <= 10) then
            y(i) = 1.0_dp + 2.0_dp*x(i,2) + 0.02_dp*real((-1)**i,dp)
         else
            y(i) = 8.0_dp - 1.5_dp*x(i,2) + 0.02_dp*real((-1)**i,dp)
         end if
      end do
      call flexmix_gaussian(x, y, 2, fit, ctl)
      call sort_by_intercept(fit, info)
      call flexmix_get_design(fit, design, info)
      call assert_true(info == 0 .and. all(shape(design) == [2,6]), 'FLXgetDesign shape', failures)
      call assert_true(all(design(1,1:3)) .and. .not. any(design(1,4:6)), 'FLXgetDesign component 1', failures)
      call assert_true(all(design(2,4:6)) .and. .not. any(design(2,1:3)), 'FLXgetDesign component 2', failures)
      call flexmix_parameter_vector(fit, parameters, info)
      call assert_true(info == 0 .and. size(parameters) == 7, 'parameter vector size', failures)
      call assert_close(parameters(3), log(fit%sigma(1)), 1.0e-12_dp, 'log sigma parameterization', failures)
      call flexmix_replace_parameters(fit, parameters, changed, info)
      call assert_true(info == 0, 'replace parameters status', failures)
      call assert_close(maxval(abs(changed%beta-fit%beta)), 0.0_dp, 1.0e-12_dp, 'replace beta roundtrip', failures)
      call assert_close(maxval(abs(changed%sigma-fit%sigma)), 0.0_dp, 1.0e-12_dp, 'replace sigma roundtrip', failures)
      call flexmix_loglik_regression_at(fit, x, y, parameters, loglik0, info)
      call assert_true(info == 0, 'parameter loglik status', failures)
      call assert_close(loglik0, fit%loglik, 2.0e-8_dp, 'parameter loglik at fit', failures)
      call assert_true(flexmix_exist_gradient(fit), 'analytic gradient availability', failures)
      perturbed = parameters
      perturbed(1) = perturbed(1) + 0.07_dp
      perturbed(4) = perturbed(4) - 0.05_dp
      perturbed(7) = perturbed(7) + 0.12_dp
      call flexmix_gradient_regression(fit, x, y, perturbed, gradient, info)
      call assert_true(info == 0 .and. size(gradient) == size(perturbed), 'analytic gradient status', failures)
      h = 1.0e-6_dp
      do i = 1, size(perturbed)
         parameters = perturbed
         parameters(i) = parameters(i) + h
         call flexmix_loglik_regression_at(fit, x, y, parameters, loglik_plus, info)
         parameters = perturbed
         parameters(i) = parameters(i) - h
         call flexmix_loglik_regression_at(fit, x, y, parameters, loglik_minus, info)
         call assert_close(gradient(i), (loglik_plus-loglik_minus)/(2.0_dp*h), 3.0e-5_dp, &
                           'analytic versus finite-difference score', failures)
      end do

      beta_fixed(:,1) = [1.0_dp, 2.0_dp]
      beta_fixed(:,2) = [4.0_dp, -1.0_dp]
      prior_fixed = [1.0_dp, 3.0_dp]
      sigma_fixed = [0.5_dp, 1.5_dp]
      call flexmix_make_distribution('gaussian', beta_fixed, prior_fixed, distribution, info, sigma=sigma_fixed)
      call assert_true(info == 0 .and. distribution%k == 2, 'FLXdist numerical constructor', failures)
      call assert_close(distribution%prior(1), 0.25_dp, 1.0e-12_dp, 'FLXdist prior normalization', failures)
      call flexmix_refit_optim_regression(fit, x, y, estimate, se, vcov, info)
      call assert_true(info == 0 .and. size(estimate) == 7 .and. all(se >= 0.0_dp), 'refit_optim covariance', failures)
      call assert_true(all(shape(vcov) == [7,7]), 'refit_optim covariance shape', failures)
      allocate(ysim(n), csim(n))
      call set_test_seed()
      call flexmix_simulate_regression(distribution, x, ysim, csim, info)
      call assert_true(info == 0 .and. all(csim >= 1) .and. all(csim <= 2), 'FLXdist simulation use', failures)
   end subroutine test_parameter_refit_tools

   subroutine test_example_generators(failures)
      integer, intent(inout) :: failures !! Running count of failed assertions updated by this test case.
      real(dp), allocatable :: x(:), yn(:), yp(:), yb(:), yclus(:,:)
      integer, allocatable :: class_label(:), id1(:), id2(:), cluster(:)
      integer :: i, info

      call set_test_seed()
      call ex_npreg(10, x, yn, yp, yb, class_label, id1, id2, info)
      call assert_true(info == 0 .and. size(x) == 20 .and. size(yn) == 20, 'ExNPreg output shape', failures)
      call assert_true(all(class_label(1:10) == 1) .and. all(class_label(11:20) == 2), 'ExNPreg classes', failures)
      call assert_true(all(id1(1:20:2) == [(i,i=1,10)]) .and. all(id1(2:20:2) == [(i,i=1,10)]), &
                       'ExNPreg id1 grouping', failures)
      call assert_true(all(id2(1:4) == 1) .and. all(id2(5:8) == 2) .and. all(id2(9:12) == 3) .and. &
                       all(id2(13:16) == 4) .and. all(id2(17:20) == 5), 'ExNPreg id2 grouping', failures)
      call assert_true(all(abs(yb-real(nint(yb),dp)) < 1.0e-12_dp) .and. all(yb >= 0.0_dp) .and. &
                       all(yb <= 1.0_dp) .and. all(yp >= 0.0_dp), 'ExNPreg response support', failures)
      call ex_npreg(9, x, yn, yp, yb, class_label, id1, id2, info)
      call assert_true(info /= 0, 'ExNPreg odd n rejection', failures)

      call set_test_seed()
      call ex_nclus(10, yclus, cluster, info)
      call assert_true(info == 0 .and. all(shape(yclus) == [55,2]), 'ExNclus output shape', failures)
      call assert_true(count(cluster == 1) == 10 .and. count(cluster == 2) == 10, 'ExNclus first counts', failures)
      call assert_true(count(cluster == 3) == 15 .and. count(cluster == 4) == 20, 'ExNclus last counts', failures)
   end subroutine test_example_generators

   subroutine test_glmfix(failures)
      integer, intent(inout) :: failures !! Running count of failed assertions updated by this test case.
      integer, parameter :: n = 20
      real(dp) :: x(n,3,2), y(n), t, xp(5,1,1), yp(5), xb(4,1,1), success(4), failure(4)
      real(dp) :: xg(4,1,1), yg(4)
      logical :: design(2,3), design1(1,1)
      integer :: initial_cluster(n), i, info
      logical, allocatable :: returned_design(:,:)
      real(dp), allocatable :: parameters(:)
      type(flexmix_result) :: fit, pfit, bfit, gfit
      type(flexmix_control) :: ctl

      x = 0.0_dp
      design(1,:) = [.true., .false., .true.]
      design(2,:) = [.false., .true., .true.]
      do i = 1, n
         t = real(mod(i-1,10),dp) / 9.0_dp
         x(i,1,1) = 1.0_dp
         x(i,3,1) = t
         x(i,2,2) = 1.0_dp
         x(i,3,2) = t
         if (i <= 10) then
            y(i) = 1.0_dp + 2.0_dp*t + 0.02_dp*real((-1)**i,dp)
            initial_cluster(i) = 1
         else
            y(i) = 8.0_dp + 2.0_dp*t + 0.02_dp*real((-1)**i,dp)
            initial_cluster(i) = 2
         end if
      end do
      ctl%minprior = 0.0_dp
      ctl%tolerance = 1.0e-10_dp
      call flexmix_glmfix_gaussian(x, y, design, fit, ctl, initial_cluster=initial_cluster)
      call assert_true(fit%status == 0 .and. fit%k == 2, 'FLXMRglmfix Gaussian status', failures)
      call assert_close(fit%beta(1,1), 0.9945454545454546_dp, 3.0e-8_dp, 'FLXMRglmfix intercept 1', failures)
      call assert_close(fit%beta(2,2), 7.994545454545454_dp, 3.0e-8_dp, 'FLXMRglmfix intercept 2', failures)
      call assert_close(fit%beta(3,1), fit%beta(3,2), 1.0e-12_dp, 'FLXMRglmfix shared slope equality', failures)
      call assert_close(fit%beta(3,1), 2.01090909090909_dp, 3.0e-8_dp, 'FLXMRglmfix shared slope', failures)
      call flexmix_get_design(fit, returned_design, info)
      call assert_true(info == 0 .and. all(shape(returned_design) == [2,5]), 'FLXMRglmfix refit design shape', failures)
      call assert_true(all(returned_design(:,3)), 'FLXMRglmfix shared design column', failures)
      call flexmix_parameter_vector(fit, parameters, info)
      call assert_true(info == 0 .and. size(parameters) == 6, 'FLXMRglmfix unique parameter vector', failures)

      design1 = .true.
      xp = 1.0_dp
      yp = [1.0_dp,2.0_dp,3.0_dp,4.0_dp,5.0_dp]
      call flexmix_glmfix_poisson(xp, yp, design1, pfit)
      call assert_true(pfit%status == 0, 'FLXMRglmfix Poisson status', failures)
      call assert_close(pfit%beta(1,1), log(3.0_dp), 2.0e-9_dp, 'FLXMRglmfix Poisson intercept', failures)

      xb = 1.0_dp
      success = [1.0_dp,2.0_dp,3.0_dp,4.0_dp]
      failure = [4.0_dp,3.0_dp,2.0_dp,1.0_dp]
      call flexmix_glmfix_binomial(xb, success, failure, design1, bfit)
      call assert_true(bfit%status == 0, 'FLXMRglmfix binomial status', failures)
      call assert_close(bfit%beta(1,1), 0.0_dp, 2.0e-9_dp, 'FLXMRglmfix binomial intercept', failures)

      xg = 1.0_dp
      yg = [1.0_dp,2.0_dp,4.0_dp,8.0_dp]
      call flexmix_glmfix_gamma(xg, yg, design1, gfit)
      call assert_true(gfit%status == 0, 'FLXMRglmfix Gamma status', failures)
      call assert_close(gfit%beta(1,1), 1.0_dp/3.75_dp, 2.0e-9_dp, 'FLXMRglmfix Gamma intercept', failures)
   end subroutine test_glmfix

   subroutine sort_by_intercept(fit, info)
      type(flexmix_result), intent(inout) :: fit !! Gaussian regression fit reordered so component intercepts increase.
      integer, intent(out) :: info !! Relabeling status returned by the public relabel operation.
      integer :: permutation(2)
      if (fit%k /= 2) then
         info = 1
         return
      end if
      if (fit%beta(1,1) <= fit%beta(1,2)) then
         permutation = [1,2]
      else
         permutation = [2,1]
      end if
      call flexmix_relabel(fit, permutation, info)
   end subroutine sort_by_intercept

   subroutine sort_by_center(fit, info)
      type(flexmix_result), intent(inout) :: fit !! Two-component normal mixture reordered by its first coordinate.
      integer, intent(out) :: info !! Relabeling status returned by the public relabel operation.
      integer :: permutation(2)
      if (fit%k /= 2) then
         info = 1
         return
      end if
      if (fit%center(1,1) <= fit%center(1,2)) then
         permutation = [1,2]
      else
         permutation = [2,1]
      end if
      call flexmix_relabel(fit, permutation, info)
   end subroutine sort_by_center


   subroutine set_test_seed()
      integer, allocatable :: seed(:)
      integer :: i, nseed
      call random_seed(size=nseed)
      allocate(seed(nseed))
      do i = 1, nseed
         seed(i) = 1729 + 97 * i
      end do
      call random_seed(put=seed)
   end subroutine set_test_seed

   subroutine assert_close(actual, expected, tolerance, label, failures)
      real(dp), intent(in) :: actual !! Computed scalar value being checked.
      real(dp), intent(in) :: expected !! Reference scalar value expected from the deterministic test.
      real(dp), intent(in) :: tolerance !! Maximum allowed absolute difference.
      character(len=*), intent(in) :: label !! Human-readable assertion label printed on failure.
      integer, intent(inout) :: failures !! Running failure count incremented when the assertion fails.
      if (abs(actual - expected) > tolerance) then
         write(*,'(a,2(1x,es24.16))') 'FAIL '//trim(label)//':', actual, expected
         failures = failures + 1
      end if
   end subroutine assert_close

   subroutine assert_true(condition, label, failures)
      logical, intent(in) :: condition !! Boolean condition that must hold for the test to pass.
      character(len=*), intent(in) :: label !! Human-readable assertion label printed on failure.
      integer, intent(inout) :: failures !! Running failure count incremented when the assertion fails.
      if (.not. condition) then
         write(*,'(a)') 'FAIL '//trim(label)
         failures = failures + 1
      end if
   end subroutine assert_true

end program test_flexmix
