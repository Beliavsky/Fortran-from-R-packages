! SPDX-License-Identifier: GPL-2.0-or-later
program distributions_and_risk
   use ghyp
   implicit none
   type(ghyp_model_type) :: nig, vg, student
   nig = nig_uv(chi=2.0_dp,psi=2.0_dp,mu=0.0_dp,sigma=1.0_dp,gamma=-0.25_dp)
   vg = vg_uv(1.5_dp,psi=3.0_dp,mu=0.0_dp,sigma=1.0_dp,gamma=0.2_dp)
   student = student_t_uv(7.0_dp,mu=0.0_dp,sigma=1.0_dp)
   print '(a,3f12.6)', 'densities at zero: ',dghyp(0.0_dp,nig),dghyp(0.0_dp,vg),dghyp(0.0_dp,student)
   print '(a,3f12.6)', '99% quantiles: ',qghyp(0.99_dp,nig),qghyp(0.99_dp,vg),qghyp(0.99_dp,student)
   print '(a,f12.6)', 'NIG lower-tail ES at 1%: ',esghyp(0.01_dp,nig,loss=.false.)
end program distributions_and_risk
