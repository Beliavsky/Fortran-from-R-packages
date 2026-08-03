! SPDX-License-Identifier: GPL-3.0-or-later
program test_mpin
   use pinstimation
   implicit none
   type(mpin_parameters) :: truth
   type(mpin_result) :: ml, ecm
   type(trade_counts) :: data
   integer, allocatable :: states(:)
   real(dp), allocatable :: posterior(:,:)
   integer :: i, layers

   allocate(truth%alpha(2),truth%delta(2),truth%mu(2))
   truth%alpha=[0.16_dp,0.24_dp]
   truth%delta=[0.30_dp,0.65_dp]
   truth%mu=[9.0_dp,23.0_dp]
   truth%eps_b=17.0_dp
   truth%eps_s=19.0_dp
   call simulate_mpin(420,truth,data,states,seed=991)
   if (.not. finite_number(mpin_loglik(data,truth))) error stop 'nonfinite MPIN likelihood'
   call mpin_posteriors(data,truth,posterior)
   do i=1,size(posterior,1)
      if(abs(sum(posterior(i,:))-1.0_dp)>1.0e-11_dp) error stop 'MPIN posteriors do not normalize'
   end do
   call fit_mpin_ml(data,2,ml,max_iterations=1800,tolerance=2.0e-7_dp)
   if(.not.finite_number(ml%log_likelihood)) error stop 'MPIN ML failed'
   if(abs(ml%mpin-mpin_value(truth))>0.10_dp) error stop 'MPIN estimate too far from truth'
   call fit_mpin_ecm(data,2,ecm,max_iterations=45,tolerance=2.0e-5_dp)
   if(.not.finite_number(ecm%log_likelihood)) error stop 'MPIN ECM failed'
   layers=detectlayers_eg(data,max_layers=4)
   if(layers<1.or.layers>4) error stop 'invalid detected layer count'
   print '(a)', 'test_mpin: PASS'
end program test_mpin
