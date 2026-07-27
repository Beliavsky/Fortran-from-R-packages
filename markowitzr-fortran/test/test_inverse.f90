! SPDX-License-Identifier: LGPL-3.0-or-later
program test_inverse
   use markowitzr, only: dp, theta_result, itheta_vcov
   use markowitzr, only: covariance_empirical, covariance_normal
   implicit none
   real(dp) :: x(6,2)
   real(dp), parameter :: expected(6) = [ &
      5.328947368421055_dp,-1.144736842105264_dp,-0.276315789473684_dp, &
      0.366541353383459_dp,-0.118421052631579_dp,0.592105263157895_dp]
   type(theta_result) :: empirical, normal
   integer :: i

   x = reshape([1.0_dp,2.0_dp,3.0_dp,4.0_dp,5.0_dp,6.0_dp, &
                2.0_dp,0.0_dp,1.0_dp,-1.0_dp,3.0_dp,2.0_dp],[6,2])

   empirical = itheta_vcov(x,covariance_method=covariance_empirical)
   if (empirical%status /= 0) error stop 1
   if (any(shape(empirical%covariance) /= [6,6])) error stop 1
   if (maxval(abs(empirical%mu-expected)) > 2.0e-13_dp) error stop 1
   if (maxval(abs(empirical%covariance-transpose(empirical%covariance))) > 0.0_dp) error stop 1

   normal = itheta_vcov(x,covariance_method=covariance_normal)
   if (normal%status /= 0) error stop 1
   if (maxval(abs(normal%mu-expected)) > 2.0e-13_dp) error stop 1
   if (any([(normal%covariance(i,i) < -1.0e-10_dp,i=1,size(normal%mu))])) error stop 1

   print '(a)', 'test_inverse: PASS'
end program test_inverse
