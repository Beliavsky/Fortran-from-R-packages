! SPDX-License-Identifier: LGPL-3.0-or-later
program test_unified
   use sharper, only: dp, moment_vcov_result, sm_vcov, ism_vcov
   implicit none
   real(dp) :: x(60,3)
   integer :: i
   type(moment_vcov_result) :: sm, ism

   do i = 1, 60
      x(i,1) = 0.01_dp+0.02_dp*sin(0.2_dp*real(i,dp))
      x(i,2) = 0.005_dp+0.7_dp*x(i,1)+0.01_dp*cos(0.31_dp*real(i,dp))
      x(i,3) = -0.002_dp+0.2_dp*x(i,1)+0.015_dp*sin(0.47_dp*real(i,dp))
   end do
   sm = sm_vcov(x,normal_model=.false.)
   if (size(sm%mean) /= 9) error stop 1
   if (any(shape(sm%covariance) /= [9,9])) error stop 1
   if (maxval(abs(sm%covariance-transpose(sm%covariance))) > 1.0e-12_dp) error stop 1

   ism = ism_vcov(x,normal_model=.true.)
   if (ism%status /= 0) error stop 1
   if (size(ism%mean) /= 9) error stop 1
   if (any(shape(ism%covariance) /= [9,9])) error stop 1
   if (maxval(abs(ism%covariance-transpose(ism%covariance))) > 1.0e-7_dp) error stop 1

   print '(a)', 'test_unified: PASS'
end program test_unified
