! SPDX-License-Identifier: GPL-3.0-only
! Derived from the GPL-3 R package poilog by Vidar Grotan and Steinar Engen.
program test_bivariate_mle
   use poilog, only : dp,bipoilog_fit,bipoilog_mle_fit,dbipoilog
   implicit none
   integer :: n1(6)=[1,2,1,3,2,4]
   integer :: n2(6)=[1,1,2,2,3,4]
   type(bipoilog_fit) :: fit
   real(dp) :: p

   p=dbipoilog(2,3,0.0_dp,0.0_dp,1.0_dp,1.0_dp,1.0_dp)
   if(p<0.0_dp .or. p>1.0_dp) error stop 1

   fit=bipoilog_mle_fit(n1,n2,start=[0.5_dp,0.5_dp,1.0_dp,1.0_dp,0.2_dp], &
      ztrunc=.true.,method='Nelder-Mead',maxit=120)
   if(fit%convergence/=0) error stop 1
   if(fit%sig1<=0.0_dp .or. fit%sig2<=0.0_dp) error stop 1
   if(fit%rho < -1.0_dp .or. fit%rho > 1.0_dp) error stop 1
   if(.not.(fit%loglik>-huge(1.0_dp))) error stop 1
   print *, 'test_bivariate_mle: PASS',fit%loglik,fit%convergence
end program test_bivariate_mle
