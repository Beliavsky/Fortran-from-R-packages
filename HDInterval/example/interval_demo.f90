! SPDX-License-Identifier: GPL-3.0-only
module interval_demo_callbacks
   use hdinterval, only: dp
   implicit none
contains
   pure function logistic_quantile(p) result(q)
      real(dp), intent(in) :: p !! Lower-tail probability strictly between zero and one.
      real(dp) :: q
      q = log(p / (1.0_dp - p))
   end function logistic_quantile

   pure function logistic_cdf(x) result(p)
      real(dp), intent(in) :: x !! Standard logistic coordinate on the real line.
      real(dp) :: p, e
      e = exp(-abs(x))
      if (x < 0.0_dp) then
         p = e / (1.0_dp + e)
      else
         p = 1.0_dp / (1.0_dp + e)
      end if
   end function logistic_cdf
end module interval_demo_callbacks

program interval_demo
   use hdinterval
   use interval_demo_callbacks
   implicit none
   real(dp) :: interval(2), height, q(3)
   real(dp), allocatable :: intervals(:, :)
   integer :: code, i
   procedure(scalar_function_interface), pointer :: callback

   call hdi_samples([1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 5.0_dp], interval, 0.6_dp, code)
   if (code /= hdi_ok) error stop 'Sample HDI failed'
   print '(a,2f12.6)', 'Sample 60% interval: ', interval

   call hdi_density([-3.0_dp, -2.0_dp, -1.0_dp, 0.0_dp, 1.0_dp, 2.0_dp, 3.0_dp], &
      [3.0_dp, 4.0_dp, 0.0_dp, 0.0_dp, 2.0_dp, 4.0_dp, 1.0_dp], intervals, height, &
      cred_mass = 0.6_dp, allow_split = .true., status = code)
   if (code /= hdi_ok) error stop 'Density HDI failed'
   print '(a,f12.6)', 'Density threshold: ', height
   do i = 1, size(intervals, 2)
      print '(a,2f12.6)', '  Modal interval: ', intervals(:, i)
   end do

   callback => logistic_quantile
   call hdi_function(callback, interval, cred_mass = 0.8_dp, status = code)
   if (code /= hdi_ok) error stop 'Quantile HDI failed'
   print '(a,2f12.6)', 'Logistic 80% HDI: ', interval

   callback => logistic_cdf
   call inverse_cdf([0.1_dp, 0.5_dp, 0.9_dp], callback, q, tol = 1e-10_dp, status = code)
   if (code /= hdi_ok) error stop 'CDF inversion failed'
   print '(a,3f12.6)', 'Logistic quantiles: ', q
end program interval_demo
