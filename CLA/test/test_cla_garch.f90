! SPDX-License-Identifier: GPL-3.0-or-later
program test_cla_garch
   use kind_mod, only: dp
   use cla, only: cla_garch_result_t, mu_sigma_garch, cla_distribution_normal, &
      cla_distribution_student
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   implicit none
   integer, parameter :: n = 90
   real(dp) :: prices(n,2), r1, r2
   type(cla_garch_result_t) :: result, student_result
   integer :: t

   prices(1,:) = [100.0_dp,80.0_dp]
   do t=2,n
      r1 = 0.0004_dp + 0.009_dp*sin(0.31_dp*real(t,dp)) + &
         0.004_dp*cos(0.07_dp*real(t*t,dp))
      r2 = 0.0002_dp + 0.006_dp*sin(0.31_dp*real(t,dp)+0.4_dp) - &
         0.005_dp*cos(0.11_dp*real(t*t,dp))
      prices(t,1) = prices(t-1,1)*exp(r1)
      prices(t,2) = prices(t-1,2)*exp(r2)
   end do

   result = mu_sigma_garch(prices,distribution=cla_distribution_normal, &
      max_iterations=120,tolerance=2.0e-5_dp)
   if(result%info/=0)error stop 'mu_sigma_garch failed'
   if(.not.all(ieee_is_finite(result%mu)))error stop 'nonfinite means'
   if(.not.all(ieee_is_finite(result%covariance)))error stop 'nonfinite covariance'
   if(maxval(abs(result%covariance-transpose(result%covariance)))>1.0e-12_dp) &
      error stop 'asymmetric covariance'
   if(maxval(abs(diagonal(result%covariance)-result%forecast_sigma**2))>1.0e-10_dp) &
      error stop 'covariance diagonal mismatch'
   if(any(result%forecast_sigma<=0.0_dp))error stop 'nonpositive forecast sigma'
   student_result = mu_sigma_garch(prices,distribution=cla_distribution_student, &
      max_iterations=140,tolerance=3.0e-5_dp)
   if(student_result%info/=0)error stop 'Student GARCH failed'
   if(any(student_result%shape<=2.0_dp))error stop 'invalid Student shape'
   if(.not.all(ieee_is_finite(student_result%covariance)))error stop 'Student covariance'
   print '(a)', 'test_cla_garch: PASS'
contains
   pure function diagonal(a) result(d)
      real(dp),intent(in)::a(:,:)
      real(dp)::d(min(size(a,1),size(a,2)))
      integer::i
      do i=1,size(d)
         d(i)=a(i,i)
      end do
   end function diagonal
end program test_cla_garch
