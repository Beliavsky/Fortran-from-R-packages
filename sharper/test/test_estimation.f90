! SPDX-License-Identifier: LGPL-3.0-or-later
program test_estimation
   use sharper, only: dp, sr_result, sropt_result, moment_vcov_result
   use sharper, only: fit_sr, fit_sropt, sr_standard_error, sr_confint
   use sharper, only: sr_vcov, sric, sample_standardized_cumulants
   implicit none
   real(dp), parameter :: x(8) = [0.01_dp,0.02_dp,-0.01_dp,0.03_dp, &
      0.0_dp,0.015_dp,0.005_dp,0.025_dp]
   real(dp), parameter :: y(8) = [0.005_dp,0.018_dp,-0.005_dp,0.022_dp, &
      0.002_dp,0.012_dp,0.007_dp,0.019_dp]
   real(dp) :: data(8,2), cumulants(4)
   real(dp), allocatable :: se(:), ci(:, :)
   type(sr_result) :: z
   type(sropt_result) :: zo
   type(moment_vcov_result) :: svc

   data(:,1) = x
   data(:,2) = y
   z = fit_sr(data,higher_order=.true.)
   call assert_close(z%value(1),0.889756521002609_dp,2.0e-12_dp)
   call assert_close(z%value(2),1.066003581778052_dp,2.0e-12_dp)
   if (any(z%df /= 7)) error stop 1
   se = sr_standard_error(z,'t')
   if (any(se <= 0.0_dp)) error stop 1
   ci = sr_confint(z,0.95_dp,'t')
   if (any(ci(:,1) >= ci(:,2))) error stop 1

   zo = fit_sropt(data)
   call assert_close(zo%value,1.3398973267965195_dp,2.0e-11_dp)
   call assert_close(zo%t2,14.362598770851674_dp,2.0e-10_dp)
   if (.not. (sric(zo) < zo%value)) error stop 1

   svc = sr_vcov(data)
   call assert_close(svc%mean(1),0.951189731211341_dp,2.0e-12_dp)
   call assert_close(svc%mean(2),1.139605764596379_dp,2.0e-12_dp)
   call assert_close(svc%covariance(1,1),0.2081393961093209_dp,3.0e-12_dp)
   call assert_close(svc%covariance(1,2),0.206986742311975_dp,1.0e-8_dp)

   cumulants = sample_standardized_cumulants(x)
   if (any(.not. (abs(cumulants) < huge(1.0_dp)))) error stop 1

   print '(a)', 'test_estimation: PASS'
contains
   subroutine assert_close(actual, expected, tolerance)
      real(dp), intent(in) :: actual, expected, tolerance
      if (abs(actual-expected) > tolerance) then
         print '(a,3es24.16)', 'assert_close failed: ',actual,expected,tolerance
         error stop 1
      end if
   end subroutine assert_close
end program test_estimation
