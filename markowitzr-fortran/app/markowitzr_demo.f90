! SPDX-License-Identifier: LGPL-3.0-or-later
! Based on MarkowitzR, copyright 2014-2020 Steven E. Pav.
program markowitzr_demo
   use markowitzr, only: dp, markowitz_result, mp_vcov, covariance_empirical
   implicit none
   real(dp) :: x(8,3)
   type(markowitz_result) :: fit
   integer :: i

   x = reshape([ &
      0.010_dp,0.018_dp,-0.006_dp,0.012_dp,0.021_dp,-0.004_dp,0.016_dp,0.009_dp, &
      0.004_dp,-0.003_dp,0.009_dp,0.006_dp,0.011_dp,0.002_dp,-0.005_dp,0.007_dp, &
      -0.002_dp,0.005_dp,0.003_dp,-0.004_dp,0.008_dp,0.006_dp,0.001_dp,0.004_dp], &
      [8,3])

   fit = mp_vcov(x,covariance_method=covariance_empirical)
   if (fit%status /= 0) then
      print '(a)', trim(fit%message)
      error stop 1
   end if

   print '(a)', 'Unscaled Markowitz portfolio:'
   do i = 1, fit%p
      print '(i3,2x,es14.6,2x,a,es14.6)', i,fit%w(i,1),'SE ', &
         sqrt(max(0.0_dp,fit%w_covariance(i,i)))
   end do
end program markowitzr_demo
