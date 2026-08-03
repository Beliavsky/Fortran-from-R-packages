! SPDX-License-Identifier: GPL-3.0-only
program demo_fitheavytail
   use fitheavytail
   implicit none
   real(dp) :: x(300,3)
   type(heavy_tail_fit) :: fitted

   call random_mvt_identity(300,3,6.0_dp,x,20260730)
   x(:,2)=0.35_dp*x(:,1)+0.9_dp*x(:,2)
   x(:,3)=-0.2_dp*x(:,1)+0.25_dp*x(:,2)+0.8_dp*x(:,3)

   call fit_mvt(x,fitted,nu_method='iterative', &
      nu_iterative_method='POP',max_iter=150)

   print '(a,f10.4)', 'estimated nu: ',fitted%nu
   print '(a)', 'estimated location:'
   print '(*(f12.6,1x))',fitted%mu
   print '(a)', 'estimated covariance:'
   print '(*(f12.6,1x))',fitted%covariance
end program demo_fitheavytail
