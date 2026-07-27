! SPDX-License-Identifier: LGPL-3.0-or-later
program test_theta
   use markowitzr, only: dp, theta_result, theta_vcov
   use markowitzr, only: covariance_empirical, covariance_normal, covariance_hac
   implicit none
   real(dp) :: x(6,2)
   real(dp), parameter :: expected_mu(6) = [ &
      1.0_dp,3.5_dp,1.166666666666667_dp,15.16666666666667_dp, &
      4.666666666666667_dp,3.166666666666667_dp]
   type(theta_result) :: empirical, normal, no_intercept, hac
   integer :: i

   x = reshape([1.0_dp,2.0_dp,3.0_dp,4.0_dp,5.0_dp,6.0_dp, &
                2.0_dp,0.0_dp,1.0_dp,-1.0_dp,3.0_dp,2.0_dp],[6,2])

   empirical = theta_vcov(x,covariance_method=covariance_empirical)
   if (empirical%status /= 0) error stop 1
   if (empirical%n /= 6 .or. empirical%pp /= 3) error stop 1
   call assert_close(empirical%mu,expected_mu,2.0e-14_dp)
   call assert_scalar(empirical%covariance(2,2),0.5833333333333333_dp,2.0e-14_dp)
   call assert_scalar(empirical%covariance(4,4),29.82777777777778_dp,2.0e-13_dp)
   call assert_scalar(empirical%covariance(5,6),3.377777777777778_dp,2.0e-13_dp)
   call assert_symmetric(empirical%covariance)
   if (maxval(abs(empirical%covariance(1,:))) > 0.0_dp) error stop 1

   normal = theta_vcov(x,covariance_method=covariance_normal)
   if (normal%status /= 0) error stop 1
   call assert_close(normal%mu,expected_mu,2.0e-14_dp)
   call assert_scalar(normal%covariance(4,4),32.66666666666667_dp,2.0e-13_dp)
   call assert_scalar(normal%covariance(5,6),3.772222222222223_dp,2.0e-13_dp)
   call assert_symmetric(normal%covariance)

   no_intercept = theta_vcov(x,fit_intercept=.false., &
      covariance_method=covariance_empirical)
   if (no_intercept%status /= 0) error stop 1
   call assert_close(no_intercept%mu,expected_mu(4:6),2.0e-14_dp)

   hac = theta_vcov(x,covariance_method=covariance_hac,hac_lags=2)
   if (hac%status /= 0) error stop 1
   call assert_symmetric(hac%covariance)
   if (any([(hac%covariance(i,i) < -1.0e-12_dp,i=1,size(hac%mu))])) error stop 1

   print '(a)', 'test_theta: PASS'

contains

   subroutine assert_close(actual, reference, tolerance)
      real(dp), intent(in) :: actual(:), reference(:), tolerance
      if (size(actual) /= size(reference)) error stop 1
      if (maxval(abs(actual-reference)) > tolerance) error stop 1
   end subroutine assert_close

   subroutine assert_scalar(actual, reference, tolerance)
      real(dp), intent(in) :: actual, reference, tolerance
      if (abs(actual-reference) > tolerance) error stop 1
   end subroutine assert_scalar

   subroutine assert_symmetric(a)
      real(dp), intent(in) :: a(:, :)
      if (maxval(abs(a-transpose(a))) > 0.0_dp) error stop 1
   end subroutine assert_symmetric

end program test_theta
