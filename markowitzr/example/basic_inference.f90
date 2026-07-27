! SPDX-License-Identifier: LGPL-3.0-or-later
! Based on MarkowitzR, copyright 2014-2020 Steven E. Pav.
program basic_inference
   use markowitzr, only: dp, theta_result, theta_vcov, itheta_vcov
   use markowitzr, only: covariance_normal
   implicit none
   real(dp) :: x(6,2)
   type(theta_result) :: second_moment, inverse_moment

   x = reshape([1.0_dp,2.0_dp,3.0_dp,4.0_dp,5.0_dp,6.0_dp, &
                2.0_dp,0.0_dp,1.0_dp,-1.0_dp,3.0_dp,2.0_dp],[6,2])

   second_moment = theta_vcov(x,covariance_method=covariance_normal)
   inverse_moment = itheta_vcov(x,covariance_method=covariance_normal)
   if (second_moment%status /= 0 .or. inverse_moment%status /= 0) error stop 1

   print '(a,*(1x,es13.5))', 'vech(Theta):',second_moment%mu
   print '(a,*(1x,es13.5))', 'vech(inv(Theta)):',inverse_moment%mu
end program basic_inference
