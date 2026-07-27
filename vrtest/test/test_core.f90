! SPDX-License-Identifier: GPL-2.0-only
! Derived from vrtest 1.2 by Jae H. Kim.
program test_core
   use vrtest, only : dp, normal_quantile, chi_square_quantile, seed_random, &
      ar1_result, lm_result, auto_vr_result, chow_denning_result, wald_result, &
      wright_result, vr_curve_result, ar1_fit, lm_statistic, &
      automatic_variance_ratio, chow_denning, wald_test, wright_tests, &
      variance_ratio_curve
   implicit none
   real(dp), parameter :: y(20) = [ &
      0.012_dp,-0.004_dp,0.009_dp,0.003_dp,-0.011_dp,0.006_dp,0.015_dp,-0.008_dp, &
      0.002_dp,0.010_dp,-0.006_dp,0.004_dp,-0.013_dp,0.007_dp,0.005_dp,-0.002_dp, &
      0.011_dp,-0.009_dp,0.008_dp,-0.001_dp ]
   integer, parameter :: kvec(3) = [2,5,10]
   type(ar1_result) :: ar
   type(lm_result) :: lm
   type(auto_vr_result) :: av
   type(chow_denning_result) :: cd
   type(wald_result) :: wt
   type(wright_result) :: wr
   type(vr_curve_result) :: curve

   call assert_close(normal_quantile(0.975_dp),1.95996398454005_dp,2.0e-9_dp,'normal quantile')
   call assert_close(chi_square_quantile(0.95_dp,1.0_dp),3.84145882069412_dp,2.0e-9_dp,'chi-square quantile')

   ar = ar1_fit(y)
   call assert_close(ar%coefficient,-0.42857142857142855_dp,2.0e-12_dp,'AR(1) coefficient')
   lm = lm_statistic(y,2)
   call assert_close(lm%variance_ratio,0.46272994280414270_dp,2.0e-12_dp,'VR(2)')
   call assert_close(lm%homoskedastic,-2.402744740330274_dp,2.0e-11_dp,'Lo-Mac M1')
   call assert_close(lm%heteroskedastic,-2.9526687685759758_dp,2.0e-11_dp,'Lo-Mac M2')

   av = automatic_variance_ratio(y)
   call assert_close(av%ar1_coefficient,ar%coefficient,1.0e-14_dp,'Auto.VR AR coefficient')
   call assert_true(av%bandwidth > 0.0_dp,'positive automatic bandwidth')
   call assert_true(abs(av%statistic) < 100.0_dp,'finite automatic statistic')

   cd = chow_denning(y,kvec)
   call assert_true(cd%cd_homoskedastic >= abs(lm%homoskedastic),'CD maximum')
   call assert_true(all(cd%critical_values > 0.0_dp),'CD critical values')
   wt = wald_test(y,kvec)
   call assert_true(wt%solve_info == 0,'Wald covariance solve')
   call assert_true(wt%statistic >= 0.0_dp,'nonnegative Wald statistic')

   wr = wright_tests(y,kvec)
   call assert_true(all(shape(wr%statistics) == [3,3]),'Wright matrix shape')
   call assert_true(all(abs(wr%statistics) < 100.0_dp),'finite Wright statistics')

   curve = variance_ratio_curve(y,10)
   call assert_true(size(curve%holding_periods) == 9,'curve length')
   call assert_close(curve%variance_ratios(1),lm%variance_ratio,1.0e-14_dp,'curve VR(2)')

   print '(a)', 'test_core: PASS'
contains
   subroutine assert_close(actual, expected, tolerance, label)
      real(dp), intent(in) :: actual, expected, tolerance
      character(len=*), intent(in) :: label
      if (abs(actual-expected) > tolerance) then
         print '(a,2(1x,es24.16))', trim(label)//' failed:',actual,expected
         error stop 1
      end if
   end subroutine assert_close

   subroutine assert_true(condition, label)
      logical, intent(in) :: condition
      character(len=*), intent(in) :: label
      if (.not. condition) then
         print '(a)', trim(label)//' failed'
         error stop 1
      end if
   end subroutine assert_true
end program test_core
