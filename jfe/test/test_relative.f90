! SPDX-License-Identifier: GPL-2.0-or-later
program test_relative
   use jfe
   implicit none

   real(dp), parameter :: rb(5) = [0.01_dp, -0.005_dp, 0.015_dp, -0.01_dp, 0.005_dp]
   real(dp), parameter :: ra(5) = [0.018_dp, -0.0065_dp, 0.025_dp, -0.0135_dp, 0.0095_dp]
   real(dp) :: alpha, te, ap, info

   alpha = capm_jensen_alpha(ra, rb, scale=12.0_dp)
   te = tracking_error(ra, rb, 12.0_dp)
   ap = active_premium(ra, rb, 12.0_dp)
   info = information_ratio(ra, rb, 12.0_dp)

   call check_close(alpha, 0.023092466587270342_dp, 1.0e-12_dp, 'Jensen alpha')
   call check_close(te, 0.020310096011589902_dp, 1.0e-12_dp, 'tracking error')
   call check_close(ap, 0.04343307000553631_dp, 1.0e-12_dp, 'active premium')
   call check_close(info, 2.1384965379164798_dp, 1.0e-12_dp, 'information ratio')

   call check_close(treynor_ratio(ra, rb, scale=12.0_dp), &
      return_annualized(ra, 12.0_dp)/1.563953488372093_dp, 1.0e-12_dp, 'Treynor ratio')
   call check_close(appraisal_ratio(ra, rb, scale=12.0_dp, method=appraisal_modified), &
      alpha/1.563953488372093_dp, 1.0e-12_dp, 'modified appraisal ratio')
   call check_close(appraisal_ratio(ra, rb, scale=12.0_dp, method=appraisal_standard), &
      alpha/(residual_population_sd(ra, rb)*sqrt(12.0_dp)), 1.0e-12_dp, &
      'standard appraisal ratio')
   call check_close(appraisal_ratio(ra, rb, scale=12.0_dp, method=appraisal_alternative), &
      alpha/(1.563953488372093_dp*sample_sd(rb)*sqrt(12.0_dp)), 1.0e-12_dp, &
      'alternative appraisal ratio')
   call check_close(m2_sortino(ra, rb, scale=12.0_dp), &
      return_annualized(ra, 12.0_dp) + sortino_ratio(ra)*sqrt(12.0_dp)* &
      (downside_deviation(rb) - downside_deviation(ra)), 1.0e-12_dp, 'M2 Sortino')

   print '(a)', 'test_relative: PASS'

contains

   pure real(dp) function sample_sd(x) result(s)
      real(dp), intent(in) :: x(:)
      real(dp) :: mu
      mu = sum(x)/real(size(x), dp)
      s = sqrt(sum((x - mu)**2)/real(size(x) - 1, dp))
   end function sample_sd

   pure real(dp) function residual_population_sd(y, x) result(s)
      real(dp), intent(in) :: y(:), x(:)
      real(dp) :: mx, my, beta, intercept
      real(dp) :: residuals(size(y))
      mx = sum(x)/real(size(x), dp)
      my = sum(y)/real(size(y), dp)
      beta = sum((x - mx)*(y - my))/sum((x - mx)**2)
      intercept = my - beta*mx
      residuals = y - intercept - beta*x
      s = sqrt(sum((residuals - sum(residuals)/real(size(residuals), dp))**2)/ &
         real(size(residuals), dp))
   end function residual_population_sd

   subroutine check_close(actual, expected, tolerance, label)
      real(dp), intent(in) :: actual, expected, tolerance
      character(*), intent(in) :: label
      if (abs(actual - expected) > tolerance) then
         print '(a,2(1x,es24.16))', 'FAIL: '//trim(label), actual, expected
         error stop 1
      end if
   end subroutine check_close

end program test_relative
