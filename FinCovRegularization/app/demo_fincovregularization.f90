! SPDX-License-Identifier: GPL-2.0-only
program demo_fincovregularization
   use fincovregularization
   use fincov_linalg, only : sample_covariance
   implicit none
   real(dp) :: returns(40,5), weights(5)
   real(dp), allocatable :: covariance(:,:)
   type(cv_result) :: selection
   integer :: i, j, status

   do i = 1, size(returns,1)
      do j = 1, size(returns,2)
         returns(i,j) = 0.012_dp*sin(0.13_dp*real(i*j,dp)) + &
            0.005_dp*cos(0.09_dp*real(i*(j+2),dp))
      end do
   end do

   call sample_covariance(returns, covariance, status)
   if (status /= fincov_ok) error stop fincov_status_message(status)

   selection = banding_cv(returns, n_cv=5, norm='F', seed=142857)
   if (selection%status /= fincov_ok) error stop fincov_status_message(selection%status)

   covariance = banding(covariance, int(selection%parameter_opt), status)
   weights = gmvp(covariance, allow_short=.false., status=status)
   if (status /= fincov_ok) error stop fincov_status_message(status)

   print '(a,f8.3)', 'Selected banding parameter: ', selection%parameter_opt
   print '(a,5f11.6)', 'Long-only GMVP weights: ', weights
end program demo_fincovregularization
