! SPDX-License-Identifier: GPL-3.0-only
!
! Experimental modern Fortran translation of computational methods from the
! R package rmgarch, copyright (C) 2008-2025 Alexios Galanos.
! Fortran translation copyright (C) 2026 translation contributors.
! Distributed under the GNU General Public License, version 3 only.
program test_core
   use rmgarch
   implicit none
   integer, parameter :: n = 300, m = 2
   real(dp) :: z(n,m), rident(m,m), value, stat, pval
   real(dp) :: a(m,1), b(m,1), c(m,m), q(m,m,n), r(m,m,n), ll(n)
   real(dp) :: factors(n,m), mixing(m,m), covariance(m,m,n), sigma(n,m)
   real(dp) :: x(n,m), fitted(n-1,m), residuals(n-1,m), forecast(4,m), simulated(20,m)
   real(dp) :: skew_stat, skew_p, kurt_stat, kurt_p
   type(varx_fit_result) :: varfit
   type(ica_result) :: icafit
   logical :: ok
   integer :: t

   call seed_rng(9876)
   rident = 0.0_dp; rident(1,1)=1.0_dp; rident(2,2)=1.0_dp
   z = 0.0_dp
   do t = 1, n
      z(t,1) = random_normal(); z(t,2) = random_normal()
   end do
   value = gaussian_copula_log_density(z(1,:),rident,ok)
   call assert_close(value,0.0_dp,1.0e-10_dp,'identity Gaussian copula')

   a(:,1) = [0.20_dp,0.15_dp]
   b(:,1) = [0.80_dp,0.82_dp]
   c = 0.0_dp; c(1,1)=0.20_dp; c(2,2)=0.20_dp; c(1,2)=0.10_dp; c(2,1)=0.10_dp
   call fdcc_filter(z,a,b,c,q,r,ll,valid=ok)
   call assert_true(ok,'FDCC filter')
   call assert_close(r(1,1,n),1.0_dp,1.0e-10_dp,'FDCC diagonal')

   factors(:,1) = 0.8_dp+0.1_dp*abs(z(:,1))
   factors(:,2) = 0.6_dp+0.1_dp*abs(z(:,2))
   mixing = reshape([1.0_dp,0.4_dp,0.2_dp,1.1_dp],[m,m])
   call gogarch_covariance(factors,mixing,covariance)
   call gogarch_sigma(factors,mixing,sigma)
   call assert_true(all(sigma > 0.0_dp),'GO-GARCH sigma')

   x(1,:) = [0.0_dp,0.0_dp]
   do t = 2, n
      x(t,1) = 0.10_dp+0.55_dp*x(t-1,1)+0.10_dp*x(t-1,2)+0.10_dp*z(t,1)
      x(t,2) =-0.05_dp+0.15_dp*x(t-1,1)+0.45_dp*x(t-1,2)+0.10_dp*z(t,2)
   end do
   varfit = fit_varx(x,1)
   call assert_true(varfit%status == 0,'VAR fit')
   call filter_varx(x,varfit,fitted,residuals)
   call forecast_varx(x,varfit,4,forecast)
   call assert_true(all(abs(forecast) < 10.0_dp),'VAR forecast')
   call simulate_varx(20,varfit,x,simulated)
   call assert_true(all(abs(simulated) < huge(1.0_dp)),'VAR simulation')

   icafit = fastica(matmul(z,transpose(mixing)),max_iterations=800)
   call assert_true(icafit%status == 0 .or. icafit%status == 1,'FastICA execution')
   call dcc_constancy_test(z,1,stat,pval)
   call assert_true(pval >= 0.0_dp .and. pval <= 1.0_dp,'DCC test p-value')
   call mardia_test(z,skew_stat,skew_p,kurt_stat,kurt_p)
   call assert_true(skew_p >= 0.0_dp .and. skew_p <= 1.0_dp,'Mardia skew p-value')
   call assert_true(kurt_p >= 0.0_dp .and. kurt_p <= 1.0_dp,'Mardia kurtosis p-value')
   print '(a)', 'Core tests passed.'
contains
   subroutine assert_true(condition,message)
      logical, intent(in) :: condition
      character(len=*), intent(in) :: message
      if (.not. condition) error stop message
   end subroutine assert_true
   subroutine assert_close(actual,expected,tolerance,message)
      real(dp), intent(in) :: actual, expected, tolerance
      character(len=*), intent(in) :: message
      if (abs(actual-expected) > tolerance) error stop message
   end subroutine assert_close
end program test_core
