! SPDX-License-Identifier: GPL-3.0-only
program tyler_cauchy_example
   use fitheavytail
   implicit none
   real(dp) :: x(250,2)
   type(heavy_tail_fit) :: tyler, cauchy

   call random_mvt_identity(250,2,5.0_dp,x,101)
   x(:,2)=0.4_dp*x(:,1)+0.8_dp*x(:,2)

   call fit_tyler(x,tyler)
   call fit_cauchy(x,cauchy)

   print '(a,f8.3)', 'Tyler implied nu: ',tyler%nu
   print '(a,f8.3)', 'Cauchy implied nu:',cauchy%nu
   print '(a)', 'Tyler covariance:'
   print '(*(f12.6,1x))',tyler%covariance
end program tyler_cauchy_example
