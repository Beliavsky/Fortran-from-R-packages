! SPDX-License-Identifier: GPL-2.0-only
program test_factor_models
   use fincovregularization
   use fincov_linalg, only : sample_covariance
   implicit none
   real(dp) :: factors(8,2), factor_vector(6), assets_one(6,2), assets(8,3)
   real(dp) :: exposure(4,2), factor_returns(8,2), residuals(8,4), fundamental_assets(8,4)
   real(dp) :: expected_one(2,2), cov_est3(3,3)
   real(dp) :: cov_est4(4,4), cov_est4_wls(4,4), stat_assets(7,3), stat_cov(3,3)
   real(dp), allocatable :: expected3(:,:), sample_cov(:,:)
   integer :: i, status, selected_k
   real(dp), parameter :: tol = 2.0e-9_dp

   factor_vector = [-2.0_dp,-1.0_dp,0.0_dp,1.0_dp,2.0_dp,3.0_dp]
   assets_one(:,1) = 1.0_dp + 2.0_dp*factor_vector
   assets_one(:,2) = -1.0_dp - 0.5_dp*factor_vector
   expected_one = reshape([14.0_dp,-3.5_dp,-3.5_dp,0.875_dp],[2,2])
   call assert_matrix_close(macro_factor_cov(assets_one, factor_vector, status), expected_one, tol, 'one-factor macro')
   call assert_true(status == fincov_ok, 'one-factor status')

   do i = 1, 8
      factors(i,1) = real(i-4,dp)
      factors(i,2) = real(mod(3*i,7)-3,dp)
   end do
   assets(:,1) = 0.3_dp + 1.2_dp*factors(:,1) - 0.4_dp*factors(:,2)
   assets(:,2) = -0.2_dp + 0.5_dp*factors(:,1) + 0.8_dp*factors(:,2)
   assets(:,3) = 0.7_dp - 0.3_dp*factors(:,1) + 1.1_dp*factors(:,2)
   call sample_covariance(assets, expected3, status)
   cov_est3 = macro_factor_cov(assets, factors, status)
   call assert_true(status == fincov_ok, 'multi-factor macro status')
   call assert_matrix_close(cov_est3, expected3, 5.0e-9_dp, 'multi-factor macro')

   exposure = reshape([&
      1.0_dp,1.0_dp,0.0_dp,0.0_dp, &
      0.0_dp,0.0_dp,1.0_dp,1.0_dp], [4,2])
   do i = 1, 8
      factor_returns(i,1) = 0.01_dp*real(i-4,dp)
      factor_returns(i,2) = 0.015_dp*real(mod(2*i,5)-2,dp)
      residuals(i,1) = 0.002_dp*real((-1)**i,dp)
      residuals(i,2) = -residuals(i,1)
      residuals(i,3) = 0.003_dp*real((-1)**(i+1),dp)
      residuals(i,4) = -residuals(i,3)
   end do
   fundamental_assets = matmul(factor_returns,transpose(exposure)) + residuals
   cov_est4 = fundamental_factor_cov(fundamental_assets, exposure, 'OLS', status)
   call assert_true(status == fincov_ok, 'fundamental OLS status')
   cov_est4_wls = fundamental_factor_cov(fundamental_assets, exposure, 'WLS', status)
   call assert_true(status == fincov_ok, 'fundamental WLS status')
   call assert_matrix_close(cov_est4_wls, cov_est4, 5.0e-10_dp, 'fundamental WLS vs OLS')
   call assert_true(minval([(cov_est4(i,i),i=1,4)]) > 0.0_dp, 'fundamental positive diagonal')

   do i = 1, 7
      stat_assets(i,1) = real(i,dp)
      stat_assets(i,2) = 0.5_dp*real(i,dp) + real(mod(i,2),dp)
      stat_assets(i,3) = -0.2_dp*real(i,dp) + real(mod(i,3),dp)
   end do
   call sample_covariance(stat_assets, sample_cov, status)
   stat_cov = stat_factor_cov(stat_assets, 3, status, selected_k)
   call assert_true(status == fincov_ok, 'stat factor status')
   call assert_true(selected_k == 3, 'stat factor selected k')
   call assert_matrix_close(stat_cov, sample_cov, 2.0e-8_dp, 'stat factor full rank')

   stat_cov = stat_factor_cov(stat_assets, status=status, selected_k=selected_k)
   call assert_true(status == fincov_ok, 'stat factor automatic status')
   call assert_true(selected_k >= 0 .and. selected_k <= 3, 'stat factor automatic k')

   print '(a)', 'test_factor_models: PASS'
contains
   subroutine assert_matrix_close(actual_matrix, expected_matrix, tolerance, label)
      real(dp), intent(in) :: actual_matrix(:,:), expected_matrix(:,:), tolerance
      character(len=*), intent(in) :: label
      if (maxval(abs(actual_matrix-expected_matrix)) > tolerance) then
         print '(a,es24.14)', trim(label)//' failed, max error: ', maxval(abs(actual_matrix-expected_matrix))
         error stop 1
      end if
   end subroutine assert_matrix_close

   subroutine assert_true(condition, label)
      logical, intent(in) :: condition
      character(len=*), intent(in) :: label
      if (.not. condition) then
         print '(a)', trim(label)//' failed'
         error stop 1
      end if
   end subroutine assert_true
end program test_factor_models
