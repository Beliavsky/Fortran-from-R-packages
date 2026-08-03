! SPDX-License-Identifier: GPL-2.0-only
program regularization_example
   use fincovregularization
   implicit none
   real(dp) :: sigma(4,4), regularized(4,4)
   integer :: status, i

   sigma = reshape([ &
      0.040_dp, 0.018_dp, 0.012_dp, 0.006_dp, &
      0.018_dp, 0.090_dp, 0.015_dp, 0.008_dp, &
      0.012_dp, 0.015_dp, 0.160_dp, 0.020_dp, &
      0.006_dp, 0.008_dp, 0.020_dp, 0.250_dp], [4,4])

   regularized = soft_thresholding(sigma, 0.01_dp, status)
   if (status /= fincov_ok) error stop fincov_status_message(status)

   print '(a)', 'Soft-thresholded covariance matrix:'
   do i = 1, 4
      print '(4f12.6)', regularized(i,:)
   end do
end program regularization_example
