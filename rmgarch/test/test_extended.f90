! SPDX-License-Identifier: GPL-3.0-only
!
! Experimental modern Fortran translation of computational methods from the
! R package rmgarch, copyright (C) 2008-2025 Alexios Galanos.
! Fortran translation copyright (C) 2026 translation contributors.
! Distributed under the GNU General Public License, version 3 only.
program test_extended
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
   use rmgarch
   implicit none

   integer, parameter :: n = 180, m = 2, nsim = 120, horizon = 2
   real(dp) :: identity(m,m), corr(m,m), mean0(m), x0(m), params(4)
   real(dp) :: z(n,m), q(m,m,n), r(m,m,n), qbar(m,m), ll
   real(dp) :: qf(m,m,3), rf(m,m,3), us(n,m), zs(n,m)
   real(dp) :: scores(n,m), uniforms(n,m), back(n,m)
   real(dp), allocatable :: fa(:,:), fb(:,:), fc(:,:)
   real(dp) :: fz(n,m), fq(m,m,n), fr(m,m,n)
   real(dp) :: ffq(m,m,2), ffr(m,m,2)
   real(dp) :: scenarios(horizon,m,nsim), scenario_cor(m,m,horizon,nsim)
   real(dp) :: scenario_mean(horizon,m), scenario_cov(m,m,horizon)
   real(dp) :: scenario_skew(m,m,m,horizon), scenario_kurt(m,m,m,m,horizon)
   real(dp) :: portfolio(horizon,nsim), weights1(1,m)
   real(dp) :: source(n,m), mixed(n,m), reconstructed(n,m)
   real(dp) :: mixing(m,m), component_mean(n,m), component_sigma(n,m), asset_returns(n,m)
   real(dp) :: factor_mean(m), factor_variance(m), factor_third(m), factor_exkurt(m)
   real(dp) :: asset_mean(m), covariance(m,m), coskew(m,m,m), cokurt(m,m,m,m)
   real(dp) :: beta2(m), beta3(m), beta4(m), port_mean, port_var, port_skew, port_kurt
   real(dp) :: gx_mean, gx_var, gx_skew, gx_kurt, value, p
   real(dp) :: margin_cov(m,m,3), margin_mean(1,m), margin_parameters(3,4)
   real(dp) :: vx(n,m), exogen(n,1)
   real(dp) :: gg_sigma(3,m), gg_cov(m,m,3), gg_cor(m,m,3), gg_sim(30,m), gg_factors(30,m)
   real(dp) :: gg_component(n,m), gg_filter_sigma(n,m), gg_filter_std(n,m)
   real(dp) :: gg_filter_cov(m,m,n), gg_filter_cor(m,m,n)
   real(dp) :: cg_sigma(n,m), cg_std(n,m), cg_u(n,m), cg_scores(n,m)
   real(dp) :: cg_q(m,m,n), cg_r(m,m,n), cg_sim(30,m)
   integer :: groups(m), i, t
   logical :: ok
   type(dcc_spec) :: spec, laplace_spec
   type(dcc_fit_result) :: dfit
   type(fdcc_fit_result) :: ffit
   type(copula_fit_result) :: cfit
   type(ica_result) :: ifit, rfit
   type(grid_distribution) :: gd1, gd2, gconv
   type(rolling_dcc_result) :: rolling, raw_rolling
   type(rolling_gogarch_result) :: go_rolling
   type(varx_fit_result) :: robust_fit
   type(gogarch_fit_result) :: ggfit
   type(multivariate_garch_fit) :: twostep
   type(copula_garch_fit_result) :: cgfit

   call seed_rng(24680)
   identity = 0.0_dp
   identity(1,1) = 1.0_dp
   identity(2,2) = 1.0_dp
   corr = reshape([1.0_dp,0.4_dp,0.4_dp,1.0_dp],[m,m])
   mean0 = 0.0_dp

   value = multivariate_normal_logpdf(mean0,mean0,identity,ok)
   call assert_true(ok .and. ieee_is_finite(value),'multivariate Normal log density')
   call assert_close(value,-log(2.0_dp*pi),1.0e-10_dp,'multivariate Normal reference')
   value = multivariate_student_logpdf(mean0,mean0,identity,8.0_dp,ok)
   call assert_true(ok .and. ieee_is_finite(value),'multivariate Student log density')
   value = multivariate_laplace_logpdf(mean0,mean0,identity,ok)
   call assert_true(ok .and. ieee_is_finite(value),'multivariate Laplace log density')
   call random_multivariate_laplace(mean0,corr,x0,ok)
   call assert_true(ok .and. all(ieee_is_finite(x0)),'multivariate Laplace simulation')
   call weighted_margin(dist_student,[0.6_dp,0.4_dp],[0.1_dp,-0.1_dp],corr,params,8.0_dp,ok)
   call assert_true(ok,'weighted margin')
   call assert_close(params(1),0.02_dp,1.0e-12_dp,'weighted mean')
   call assert_true(params(2) > 0.0_dp .and. abs(params(4)-8.0_dp) < 1.0e-12_dp,'weighted scale and shape')
   weights1 = reshape([0.6_dp,0.4_dp],[1,m])
   margin_mean = reshape([0.1_dp,-0.1_dp],[1,m])
   margin_cov(:,:,1) = corr
   margin_cov(:,:,2) = 1.2_dp*corr
   margin_cov(:,:,3) = 0.8_dp*corr
   call weighted_margin_path(dist_laplace,weights1,margin_mean,margin_cov, &
      margin_parameters,valid=ok)
   call assert_true(ok .and. all(margin_parameters(:,2) > 0.0_dp), &
      'weighted margin path')
   p = standardized_student_cdf(0.7_dp,8.0_dp)
   call assert_close(standardized_student_quantile(p,8.0_dp),0.7_dp,2.0e-6_dp, &
      'Student CDF quantile inversion')
   p = standardized_laplace_cdf(-0.4_dp)
   call assert_close(standardized_laplace_quantile(p),-0.4_dp,1.0e-12_dp, &
      'Laplace CDF quantile inversion')
   value = modified_bessel_k_half_integer_order(0.5_dp,1.3_dp,ok)
   call assert_true(ok,'Bessel K evaluation')
   call assert_close(value,sqrt(pi/(2.0_dp*1.3_dp))*exp(-1.3_dp),1.0e-12_dp, &
      'Bessel K half order reference')

   qbar = corr
   spec = make_dcc_spec([0.03_dp,0.02_dp],[0.92_dp])
   call simulate_dcc(n,spec,qbar,z,q,r,burn=80)
   ll = dcc_log_likelihood(z,spec,ok)
   call assert_true(ok .and. ieee_is_finite(ll),'DCC(2,1) likelihood')
   dfit = fit_dcc(z,p=2,q=1,max_iterations=140)
   call assert_true(dfit%status == 0 .or. dfit%status == 1,'generic DCC fit')
   call assert_true(size(dfit%spec%alpha) == 2,'generic DCC order retained')
   call dcc_forecast_history(dfit%spec,dfit%qbar,dfit%nbar,dfit%q, &
      dfit%standardized_residuals,3,qf,rf)
   call assert_close(rf(1,1,3),1.0_dp,1.0e-10_dp,'DCC history forecast')

   laplace_spec = make_dcc_spec([0.05_dp],[0.90_dp],distribution=dist_laplace)
   call simulate_dcc(n,laplace_spec,qbar,z,q,r,burn=60)
   ll = dcc_log_likelihood(z,laplace_spec,ok)
   call assert_true(ok .and. ieee_is_finite(ll),'Laplace DCC likelihood and simulation')
   dfit = fit_dcc(z,p=1,q=1,distribution=dist_laplace,max_iterations=100)
   call assert_true(dfit%status == 0 .or. dfit%status == 1,'Laplace DCC fit')
   spec = make_dcc_spec([0.05_dp],[0.90_dp],distribution=dist_student,shape=7.0_dp)
   call simulate_dcc(n,spec,qbar,z,q,r,burn=60)
   dfit = fit_dcc(z,p=1,q=1,distribution=dist_student,shape=8.0_dp, &
      estimate_shape=.true.,max_iterations=120)
   call assert_true(dfit%status == 0 .or. dfit%status == 1,'Student DCC shape fit')
   call assert_true(dfit%spec%shape > 2.01_dp .and. dfit%spec%shape < 50.0_dp, &
      'Student DCC estimated shape bounds')

   groups = [1,2]
   call fdcc_parameter_matrices(groups,[0.18_dp,0.14_dp],[0.82_dp,0.84_dp],fa,fb,fc,ok)
   call assert_true(ok,'FDCC parameter matrices')
   call simulate_fdcc(n,fa,fb,fc,qbar,fz,fq,fr,burn=60)
   ffit = fit_fdcc11(fz,groups,max_iterations=140)
   call assert_true(ffit%status == 0 .or. ffit%status == 1,'FDCC grouped fit')
   call fdcc_forecast(ffit%a,ffit%b,ffit%c,ffit%qbar,ffit%q(:,:,n),fz(n,:),2,ffq,ffr)
   call assert_close(ffr(2,2,2),1.0_dp,1.0e-10_dp,'FDCC forecast')

   call simulate_static_copula(n,corr,dist_gaussian,scores,uniforms,valid=ok)
   call assert_true(ok .and. minval(uniforms) > 0.0_dp .and. maxval(uniforms) < 1.0_dp, &
      'static Gaussian copula simulation')
   cfit = fit_copula(scores,distribution=dist_gaussian,time_varying=.false.)
   call assert_true(cfit%status == 0,'static copula fit')
   call copula_score_transform(uniforms,dist_gaussian,back)
   call assert_true(maxval(abs(back-scores)) < 1.0e-6_dp,'Gaussian copula transform inversion')
   call score_to_uniform_transform(scores,dist_gaussian,back)
   call assert_true(maxval(abs(back-uniforms)) < 1.0e-12_dp,'score to uniform transform')
   call empirical_uniform_transform(scores,us)
   call normal_score_transform(us,zs)
   call assert_true(all(ieee_is_finite(zs)),'empirical and normal score transforms')
   cfit = fit_copula(scores,distribution=dist_gaussian,time_varying=.true., &
      max_iterations=120)
   call assert_true(cfit%status == 0 .or. cfit%status == 1,'dynamic copula fit')
   call simulate_static_copula(n,corr,dist_student,scores,uniforms,shape=8.0_dp,valid=ok)
   call assert_true(ok,'static Student copula simulation')
   value = student_copula_log_density(scores(1,:),corr,8.0_dp,ok)
   call assert_true(ok .and. ieee_is_finite(value),'Student copula density')
   cfit = fit_copula(scores,distribution=dist_student,time_varying=.false.,shape=8.0_dp)
   call assert_true(cfit%status == 0,'static Student copula fit')
   call copula_score_transform(uniforms,dist_student,back,8.0_dp)
   call assert_true(maxval(abs(back-scores)) < 3.0e-6_dp,'Student copula transform inversion')
   spec = make_dcc_spec([0.05_dp],[0.90_dp],distribution=dist_student,shape=8.0_dp)
   value = dynamic_student_copula_log_likelihood(scores,spec,8.0_dp,ok)
   call assert_true(ok .and. ieee_is_finite(value),'dynamic Student copula likelihood')
   call simulate_dynamic_copula(n,spec,qbar,scores,uniforms,q,r,burn=40)
   call assert_true(minval(uniforms) > 0.0_dp .and. maxval(uniforms) < 1.0_dp, &
      'dynamic copula simulation')
   cfit = fit_copula(scores,distribution=dist_student,time_varying=.true., &
      shape=8.0_dp,max_iterations=100)
   call assert_true(cfit%status == 0 .or. cfit%status == 1,'dynamic Student copula fit')

   do t = 1, n
      call random_number(source(t,1))
      call random_number(source(t,2))
      source(t,1) = -log(max(source(t,1),tiny(1.0_dp)))-1.0_dp
      source(t,2) = 2.0_dp*source(t,2)-1.0_dp
   end do
   mixing = reshape([1.0_dp,0.45_dp,0.25_dp,1.1_dp],[m,m])
   mixed = matmul(source,transpose(mixing))
   ifit = fastica(mixed,max_iterations=700,nonlinearity='tanh')
   call assert_true(ifit%status == 0 .or. ifit%status == 1,'FastICA extended result')
   reconstructed = matmul(ifit%sources,transpose(ifit%mixing))+spread(ifit%center,1,n)
   call assert_true(maxval(abs(reconstructed-mixed)) < 1.0e-6_dp,'FastICA reconstruction')
   rfit = radical(mixed,max_sweeps=4,angle_points=25)
   call assert_true(rfit%status == 0,'RADICAL pairwise ICA')
   reconstructed = matmul(rfit%sources,transpose(rfit%mixing))+spread(rfit%center,1,n)
   call assert_true(maxval(abs(reconstructed-mixed)) < 1.0e-6_dp,'RADICAL reconstruction')

   factor_mean = [0.01_dp,-0.02_dp]
   factor_variance = [1.0_dp,0.6_dp]
   factor_third = [0.3_dp,-0.1_dp]
   factor_exkurt = [1.5_dp,0.5_dp]
   call gogarch_moments_at(mixing,factor_mean,factor_variance,factor_third, &
      factor_exkurt,asset_mean,covariance,coskew,cokurt)
   weights1 = reshape([0.6_dp,0.4_dp],[1,m])
   beta2 = portfolio_covariance_beta(weights1(1,:),covariance)
   beta3 = portfolio_coskew_beta(weights1(1,:),coskew)
   beta4 = portfolio_cokurt_beta(weights1(1,:),cokurt)
   call assert_close(dot_product(weights1(1,:),beta2),1.0_dp,1.0e-10_dp,'covariance beta')
   call assert_close(dot_product(weights1(1,:),beta3),1.0_dp,1.0e-10_dp,'coskew beta')
   call assert_close(dot_product(weights1(1,:),beta4),1.0_dp,1.0e-10_dp,'cokurt beta')
   call portfolio_factor_moments(weights1(1,:),mixing,factor_mean,factor_variance, &
      factor_third,factor_exkurt,port_mean,port_var,port_skew,port_kurt)
   call assert_true(port_var > 0.0_dp .and. port_kurt > 0.0_dp,'portfolio factor moments')
   component_mean = spread(factor_mean,1,n)
   component_sigma = spread(sqrt(factor_variance),1,n)
   call simulate_gogarch(source,component_mean,component_sigma,mixing,asset_returns)
   call assert_true(all(ieee_is_finite(asset_returns)),'GO-GARCH factor simulation')
   ggfit = fit_gogarch11(mixed,max_iterations=140)
   call assert_true(ggfit%status == 0 .or. ggfit%status == 1,'GO-GARCH model fit')
   call forecast_gogarch11(ggfit,3,gg_sigma,gg_cov,gg_cor)
   call assert_close(gg_cor(1,1,3),1.0_dp,1.0e-10_dp,'GO-GARCH forecast')
   call simulate_fitted_gogarch11(30,ggfit,gg_sim,gg_factors,burn=30)
   call assert_true(all(ieee_is_finite(gg_sim)),'GO-GARCH fitted simulation')
   call filter_gogarch11(mixed,ggfit,gg_component,gg_filter_sigma,gg_filter_std, &
      gg_filter_cov,gg_filter_cor,ok)
   call assert_true(ok .and. all(ieee_is_finite(gg_filter_cov)),'GO-GARCH filtering')
   call assert_true(maxval(abs(gg_component-ggfit%ica%sources)) < 1.0e-10_dp, &
      'GO-GARCH filtered components')
   call assert_close(gg_filter_cor(1,1,n),1.0_dp,1.0e-10_dp,'GO-GARCH filtered correlation')

   twostep = fit_two_step_dcc_general(mixed,1,1,0,max_iterations=120)
   call assert_true(twostep%status == 0 .or. twostep%status == 1,'general two-step DCC')
   call assert_close(twostep%dcc%r(1,1,n),1.0_dp,1.0e-10_dp,'two-step DCC covariance')

   cgfit = fit_copula_garch11(mixed,copula_distribution=dist_gaussian, &
      time_varying=.true.,max_iterations=120)
   call assert_true(cgfit%status == 0 .or. cgfit%status == 1,'copula-GARCH model fit')
   call filter_copula_garch11(mixed,cgfit,cg_sigma,cg_std,cg_u,cg_scores,cg_q,cg_r,ok)
   call assert_true(ok .and. all(ieee_is_finite(cg_scores)),'copula-GARCH filtering')
   call assert_close(cg_r(1,1,n),1.0_dp,1.0e-10_dp,'copula-GARCH correlation')
   call simulate_fitted_copula_garch11(30,cgfit,cg_sim,burn=30,valid=ok)
   call assert_true(ok .and. all(ieee_is_finite(cg_sim)),'copula-GARCH simulation')

   gd1 = normal_grid_distribution(0.0_dp,1.0_dp,-6.0_dp,6.0_dp,257)
   gd2 = normal_grid_distribution(0.0_dp,1.0_dp,-6.0_dp,6.0_dp,257)
   gconv = convolve_grid_distributions(gd1,gd2)
   call assert_true(gconv%status == 0,'FFT grid convolution')
   call grid_moments(gconv,gx_mean,gx_var,gx_skew,gx_kurt)
   call assert_close(gx_mean,0.0_dp,2.0e-2_dp,'FFT convolution mean')
   call assert_close(gx_var,2.0_dp,5.0e-2_dp,'FFT convolution variance')
   call assert_close(grid_cdf(gconv,grid_quantile(gconv,0.25_dp)),0.25_dp,2.0e-2_dp, &
      'grid CDF quantile')
   call assert_true(grid_density(gconv,0.0_dp) > 0.0_dp,'grid density')

   dfit = fit_dcc(scores,p=1,q=1,max_iterations=120)
   call simulate_dcc_scenarios(dfit,horizon,nsim,scenarios,scenario_cor,valid=ok)
   call assert_true(ok,'conditional DCC scenarios')
   call scenario_moments(scenarios,scenario_mean,scenario_cov,scenario_skew,scenario_kurt)
   call assert_true(all([(scenario_cov(i,i,1) > 0.0_dp,i=1,m)]),'scenario moments')
   call assert_true(all(ieee_is_finite(scenario_skew)) .and. &
      all(ieee_is_finite(scenario_kurt)),'scenario higher moments')
   call portfolio_scenarios(scenarios,weights1,portfolio)
   call assert_true(all(ieee_is_finite(portfolio)),'portfolio scenarios')

   rolling = roll_dcc(scores,window=90,refit_every=90,horizon=2,max_iterations=100)
   call assert_true(size(rolling%origin) == 2,'rolling DCC fit count')
   call assert_close(rolling%correlation(1,1,1,1),1.0_dp,1.0e-10_dp,'rolling DCC forecast')
   raw_rolling = roll_two_step_dcc(mixed,window=90,refit_every=90,horizon=2, &
      max_iterations=80)
   call assert_true(size(raw_rolling%origin) == 2,'raw rolling DCC fit count')
   call assert_true(raw_rolling%covariance(1,1,1,1) > 0.0_dp,'raw rolling DCC covariance')
   call assert_close(raw_rolling%correlation(1,1,1,1),1.0_dp,1.0e-10_dp, &
      'raw rolling DCC correlation')
   go_rolling = roll_gogarch11(mixed,window=90,refit_every=90,horizon=2, &
      max_iterations=80)
   call assert_true(size(go_rolling%origin) == 2,'rolling GO-GARCH fit count')
   call assert_true(go_rolling%covariance(1,1,1,1) > 0.0_dp,'rolling GO-GARCH covariance')
   call assert_close(go_rolling%correlation(1,1,1,1),1.0_dp,1.0e-10_dp, &
      'rolling GO-GARCH correlation')

   vx(1,:) = 0.0_dp
   exogen(:,1) = [(sin(0.05_dp*real(t,dp)),t=1,n)]
   do t = 2, n
      vx(t,1) = 0.1_dp+0.5_dp*vx(t-1,1)+0.1_dp*vx(t-1,2)+0.2_dp*exogen(t,1)+0.05_dp*source(t,1)
      vx(t,2) =-0.1_dp+0.2_dp*vx(t-1,1)+0.4_dp*vx(t-1,2)-0.1_dp*exogen(t,1)+0.05_dp*source(t,2)
   end do
   vx(n/2,:) = vx(n/2,:)+20.0_dp
   robust_fit = fit_varx_robust(vx,1,exogen=exogen,max_iterations=80)
   call assert_true(robust_fit%status == 0 .or. robust_fit%status == 1,'robust VARX fit')
   call assert_true(all(ieee_is_finite(robust_fit%coefficients)),'robust VARX coefficients')

   print '(a)', 'Extended computational tests passed.'

contains

   subroutine assert_true(condition,message)
      logical, intent(in) :: condition
      character(len=*), intent(in) :: message
      if (.not. condition) error stop message
   end subroutine assert_true

   subroutine assert_close(actual,expected,tolerance,message)
      real(dp), intent(in) :: actual, expected, tolerance
      character(len=*), intent(in) :: message
      if (abs(actual-expected) > tolerance) then
         print '(a,3(1x,es14.6))', trim(message),actual,expected,tolerance
         error stop message
      end if
   end subroutine assert_close

end program test_extended
