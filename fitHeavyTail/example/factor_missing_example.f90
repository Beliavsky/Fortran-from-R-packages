! SPDX-License-Identifier: GPL-3.0-only
program factor_missing_example
   use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
   use fitheavytail
   implicit none
   real(dp) :: x(200,4), nan_value
   type(heavy_tail_fit) :: fitted

   call random_mvt_identity(200,4,7.0_dp,x,202)
   x(:,3)=0.5_dp*x(:,1)+0.7_dp*x(:,3)
   x(:,4)=-0.3_dp*x(:,2)+0.8_dp*x(:,4)
   nan_value=ieee_value(0.0_dp,ieee_quiet_nan)
   x(1,1)=nan_value
   x(2,2)=nan_value

   call fit_mvt(x,fitted,fixed_nu=7.0_dp,na_rm=.false., &
      factors=2,max_iter=100)

   print '(a,l1)', 'converged: ',fitted%converged
   print '(a)', 'factor loadings for covariance:'
   print '(*(f12.6,1x))',fitted%loadings
   print '(a)', 'idiosyncratic variances:'
   print '(*(f12.6,1x))',fitted%psi
end program factor_missing_example
