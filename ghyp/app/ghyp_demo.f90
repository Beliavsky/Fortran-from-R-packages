! SPDX-License-Identifier: GPL-2.0-or-later
program ghyp_demo
   use ghyp
   implicit none
   type(ghyp_model_type) :: model
   type(moments_result) :: moments
   real(dp), allocatable :: sample(:,:)
   logical :: ok

   model = ghyp_ad(0.7_dp,1.8_dp,1.2_dp,[0.3_dp],[0.2_dp],reshape([1.0_dp],[1,1]))
   if (.not. model%ok) error stop trim(model%message)
   moments = ghyp_moments(model)
   print '(a,a)', 'family: ',trim(ghyp_family_name(model))
   print '(a,f12.6)', 'mean: ',moments%mean(1)
   print '(a,f12.6)', 'variance: ',moments%covariance(1,1)
   print '(a,f12.6)', '95% quantile: ',qghyp(0.95_dp,model)
   print '(a,f12.6)', '95% expected shortfall: ',esghyp(0.95_dp,model,loss=.true.)
   call rghyp(5,model,sample,ok,20260801_i8)
   if (.not. ok) error stop 'simulation failed'
   print '(a,5f11.5)', 'sample: ',sample(:,1)
end program ghyp_demo
