! SPDX-License-Identifier: GPL-3.0-or-later
program test_adjpin
   use pinstimation
   implicit none
   type(adjpin_parameters) :: truth
   type(adjpin_restrictions) :: restrictions
   type(adjpin_result) :: ml, ecm
   type(trade_counts) :: data
   integer,allocatable :: states(:)
   real(dp) :: w(6), a, s

   truth=adjpin_parameters(0.32_dp,0.44_dp,0.18_dp,0.18_dp,18.0_dp,18.0_dp,13.0_dp,13.0_dp,5.0_dp,5.0_dp)
   restrictions%equal_theta=.true.
   restrictions%equal_eps=.true.
   restrictions%equal_mu=.true.
   restrictions%equal_d=.true.
   w=adjpin_distribution(truth)
   if(abs(sum(w)-1.0_dp)>1.0e-13_dp.or.minval(w)<0.0_dp) error stop 'invalid AdjPIN distribution'
   call adjpin_values(truth,a,s)
   if(a<=0.0_dp.or.s<=0.0_dp.or.a+s>=1.0_dp) error stop 'invalid AdjPIN values'
   call simulate_adjpin(360,truth,data,states,seed=141)
   call fit_adjpin_ml(data,ml,restrictions=restrictions,max_iterations=1600,tolerance=3.0e-7_dp)
   if(.not.finite_number(ml%log_likelihood)) error stop 'AdjPIN ML failed'
   if(abs(ml%parameters%theta-ml%parameters%theta_p)>1.0e-12_dp) error stop 'theta restriction failed'
   if(abs(ml%parameters%eps_b-ml%parameters%eps_s)>1.0e-12_dp) error stop 'eps restriction failed'
   if(ml%adjpin<0.0_dp.or.ml%psos<0.0_dp.or.ml%adjpin+ml%psos>=1.0_dp) error stop 'invalid estimated measures'
   call fit_adjpin_ecm(data,ecm,restrictions=restrictions,max_iterations=18,tolerance=1.0e-4_dp)
   if(.not.finite_number(ecm%log_likelihood)) error stop 'AdjPIN ECM failed'
   print '(a)', 'test_adjpin: PASS'
end program test_adjpin
