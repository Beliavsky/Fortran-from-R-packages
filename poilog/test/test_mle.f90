! SPDX-License-Identifier: GPL-3.0-only
! Derived from the GPL-3 R package poilog by Vidar Grotan and Steinar Engen.
program test_mle
   use poilog, only : dp,poilog_fit,poilog_mle_fit
   implicit none
   integer :: x(20)=[1,1,1,2,2,2,2,3,3,4,4,5,6,7,8,9,10,12,15,20]
   type(poilog_fit) :: fit
   fit=poilog_mle_fit(x,start_mu=1.0_dp,start_sig=1.5_dp,ztrunc=.true.,maxit=80)
   if(.not.(fit%sig>0.0_dp)) error stop 1
   if(.not.(fit%loglik>-huge(1.0_dp))) error stop 1
   print *, 'test_mle: PASS',fit%mu,fit%sig,fit%loglik,fit%convergence
end program test_mle
