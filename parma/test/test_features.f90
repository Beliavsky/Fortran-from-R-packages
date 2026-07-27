! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2012-2014 Alexios Galanos and Bernhard Pfaff.
! Modern Fortran translation Copyright (C) 2026 OpenAI.
! This derivative work is distributed under GPL-3.0-or-later.
program test_features
   use parma
   implicit none
   type(parma_spec) :: spec
   type(parma_port) :: port
   type(parma_options) :: options
   real(dp) :: data(4,2),w(2),r(4),centered(4),cov(2,2),mu(2),lb(2),ub(2)
   real(dp) :: x3(3),lb3(3),ub3(3),repaired3(3),weights3(3),forecast3(3)
   real(dp) :: benchmark_cov(3),expected
   integer :: info

   data = reshape([0.01_dp,-0.02_dp,0.03_dp,0.00_dp, &
                   0.02_dp,0.01_dp,-0.01_dp,0.04_dp],[4,2])
   w = [0.7_dp,0.3_dp]
   r = matmul(data,w)
   centered = r-sum(r)/4.0_dp
   expected = sum(max(-centered,0.0_dp))/4.0_dp
   call assert_close(lpm_risk(w,data,999.0_dp,1.0_dp),expected,1.0e-12_dp,'threshold 999')
   call assert_close(riskfun(w,data,'lpm',moment=1.0_dp,threshold=999.0_dp), &
      expected,1.0e-12_dp,'riskfun dispatch')

   cov = reshape([1.0_dp,0.0_dp,0.0_dp,4.0_dp],[2,2])
   benchmark_cov = [0.5_dp,0.2_dp,0.1_dp]
   call assert_close(benchmark_variance([0.8_dp,0.2_dp],cov,benchmark_cov), &
      0.94_dp,1.0e-12_dp,'benchmark variance')

   mu = [0.01_dp,0.02_dp]
   lb = 0.0_dp
   ub = 1.0_dp
   call parmaspec(spec,risk=risk_ev,objective=solve_min_risk,mu=mu,covariance=cov, &
      lb=lb,ub=ub,info=info)
   options%max_iter = 2000
   options%tol = 1.0e-10_dp
   call parmasolve(spec,port,options)
   if (port%status /= 0) error stop 'covariance-only portfolio failed'
   call assert_close(port%weights(1),0.8_dp,2.0e-4_dp,'covariance-only weight')

   x3 = [0.9_dp,-0.6_dp,0.2_dp]
   lb3 = -1.0_dp
   ub3 = 1.0_dp
   call parmaspec(spec,risk=risk_ev,objective=solve_min_risk,mu=[0.01_dp,0.02_dp,0.03_dp], &
      covariance=eye(3),lb=lb3,ub=ub3,leverage=1.5_dp,max_positions=2,info=info)
   call repair_weights(x3,spec,repaired3,info)
   call assert_close(sum(abs(repaired3)),1.5_dp,1.0e-10_dp,'leverage repair')
   if (count(abs(repaired3)>1.0e-8_dp) > 2) error stop 'cardinality repair failed'

   forecast3 = [0.02_dp,0.01_dp,-0.01_dp]
   call simweights(weights3,[0.0_dp,0.0_dp,0.0_dp],[1.0_dp,1.0_dp,1.0_dp], &
      1.0_dp,forecast3,.true.,13579,info)
   if (info /= 0) error stop 'simweights failed'
   call assert_close(sum(weights3),1.0_dp,1.0e-12_dp,'simweights budget')
   if (dot_product(weights3,forecast3) < 0.0_dp) error stop 'simweights forecast failed'

   print '(a)', 'test_features: PASS'

contains

   subroutine assert_close(actual,wanted,epsilonx,name)
      real(dp), intent(in) :: actual,wanted,epsilonx
      character(len=*), intent(in) :: name
      if (abs(actual-wanted) > epsilonx*max(1.0_dp,abs(wanted))) then
         write(*,'(a,2es24.14)') trim(name)//' failed: ',actual,wanted
         error stop 1
      end if
   end subroutine assert_close

end program test_features
