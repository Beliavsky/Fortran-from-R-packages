! SPDX-License-Identifier: GPL-3.0-only
program skew_t_example
   use fitheavytail
   implicit none
   real(dp) :: x(120,2), initial_gamma(2)
   type(heavy_tail_fit) :: fitted

   call random_mvt_identity(120,2,8.0_dp,x,303)
   x(:,1)=x(:,1)+0.15_dp*abs(x(:,2))
   initial_gamma=[0.1_dp,0.0_dp]

   call fit_mvst(x,fitted,fixed_nu=8.0_dp, &
      initial_gamma=initial_gamma,max_iter=100)

   print '(a)', 'estimated skewness vector:'
   print '(*(f12.6,1x))',fitted%gamma
   print '(a)', 'estimated mean:'
   print '(*(f12.6,1x))',fitted%mean
   print '(a)', 'estimated covariance:'
   print '(*(f12.6,1x))',fitted%covariance
end program skew_t_example
