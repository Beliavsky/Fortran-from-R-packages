! SPDX-License-Identifier: GPL-3.0-only
program test_regression
  use mass
  use test_support
  implicit none
  integer, parameter :: n=60
  real(dp) :: x(n,2), raw(n,3), y(n), y_out(n), lambdas(3), xmean(1)
  real(dp), allocatable :: ridge_coef(:,:), profile(:)
  type(regression_result) :: fit, robust_fit, lts_fit
  type(ridge_result) :: ridge
  type(model_selection_result) :: selected
  integer :: i, status

  do i=1,n
    x(i,1)=1.0_dp
    x(i,2)=-1.0_dp+2.0_dp*real(i-1,dp)/real(n-1,dp)
    raw(i,1)=x(i,2)
    raw(i,2)=sin(0.37_dp*real(i,dp))
    raw(i,3)=cos(0.23_dp*real(i,dp))
    y(i)=1.0_dp+2.0_dp*x(i,2)+0.02_dp*sin(real(i,dp))
  end do
  call linear_model_fit(x,y,fit)
  call assert_true(fit%status == mass_success, 'linear model status')
  call assert_close(fit%coefficients(1),1.0_dp,0.01_dp,'linear intercept')
  call assert_close(fit%coefficients(2),2.0_dp,0.01_dp,'linear slope')

  y_out=y
  y_out(n)=y_out(n)+50.0_dp
  call rlm_fit(x,y_out,robust_fit,psi='bisquare')
  call assert_true(robust_fit%status == mass_success, 'rlm status')
  call assert_close(robust_fit%coefficients(2),2.0_dp,0.08_dp,'robust slope')
  call ltsreg(raw(:,1:1),y_out,lts_fit,nsamp=300,seed=345)
  call assert_true(lts_fit%status == mass_success, 'lts status')
  call assert_close(lts_fit%coefficients(size(lts_fit%coefficients)),2.0_dp,0.12_dp,'lts slope')

  lambdas=[0.0_dp,0.1_dp,1.0_dp]
  call lm_ridge(raw(:,1:1),y,lambdas,ridge)
  call assert_true(ridge%status == mass_success, 'ridge status')
  xmean=[sum(raw(:,1))/real(n,dp)]
  call ridge_coefficients(ridge,xmean,ridge_coef)
  call assert_close(ridge_coef(2,1),2.0_dp,0.02_dp,'ridge lambda zero')

  call step_aic_linear(raw,y,selected,direction='both')
  call assert_true(selected%status == mass_success, 'step AIC status')
  call assert_true(selected%selected(1), 'step AIC selects signal')
  call assert_true(.not. selected%selected(2) .and. .not. selected%selected(3), &
    'step AIC excludes noise')

  call boxcox_profile(x, y-minval(y)+1.0_dp, [-1.0_dp,0.0_dp,1.0_dp], profile, status)
  call assert_true(status == mass_success, 'boxcox status')
  call assert_all_finite(profile, 'boxcox finite')
  write(*,'(a)') 'test_regression: PASS'
end program test_regression
