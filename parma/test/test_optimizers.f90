! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2012-2014 Alexios Galanos and Bernhard Pfaff.
! Modern Fortran translation Copyright (C) 2026 OpenAI.
! This derivative work is distributed under GPL-3.0-or-later.
program test_optimizers
   use parma
   implicit none
   type(qp_result) :: qres
   type(lp_result) :: lres
   type(milp_result) :: mres
   type(cmaes_result) :: cres
   type(parma_spec) :: spec
   type(parma_port) :: port
   type(parma_options) :: opt
   real(dp) :: h(2,2),g(2),lb(2),ub(2),a_lp(3,2),b_lp(3),c_lp(2)
   real(dp) :: a_bin(1,2),lower_bin(1),upper_bin(1),data(4,2),cov(2,2),mu(2)
   integer :: info

   h = reshape([2.0_dp,0.0_dp,0.0_dp,8.0_dp],[2,2])
   g = 0.0_dp
   lb = 0.0_dp
   ub = 1.0_dp
   call qp_box_budget(h,g,lb,ub,1.0_dp,qres)
   call assert_close(qres%x(1),0.8_dp,1.0e-5_dp,'QP weight 1')
   call assert_close(qres%x(2),0.2_dp,1.0e-5_dp,'QP weight 2')

   c_lp = [3.0_dp,2.0_dp]
   a_lp = reshape([1.0_dp,1.0_dp,0.0_dp,1.0_dp,0.0_dp,1.0_dp],[3,2])
   b_lp = [4.0_dp,2.0_dp,3.0_dp]
   call lp_simplex(c_lp,a_lp,b_lp,lres,maximize=.true.)
   call assert_close(lres%objective,10.0_dp,1.0e-10_dp,'LP objective')

   a_bin = reshape([1.0_dp,1.0_dp],[1,2])
   lower_bin = [-huge(1.0_dp)]
   upper_bin = [1.0_dp]
   call milp_binary_solve([5.0_dp,3.0_dp],mres,a_bin,lower_bin,upper_bin,maximize=.true.)
   if (mres%status /= 0 .or. any(mres%x /= [1,0])) error stop 'MILP failed'

   opt%max_iter = 800
   opt%seed = 7788
   opt%sigma0 = 0.3_dp
   call cmaes_minimize(sphere,[3.0_dp,-2.0_dp],cres,opt,[-5.0_dp,-5.0_dp],[5.0_dp,5.0_dp])
   if (cres%value > 1.0e-8_dp) then
      write(*,'(a,es14.6)') 'CMA-ES sphere failed: ',cres%value
      error stop 1
   end if

   data = reshape([0.01_dp,0.02_dp,-0.01_dp,0.00_dp, &
                   0.02_dp,-0.02_dp,0.01_dp,0.00_dp],[4,2])
   cov = reshape([1.0_dp,0.0_dp,0.0_dp,4.0_dp],[2,2])
   mu = [0.01_dp,0.02_dp]
   call parmaspec(spec,data=data,risk=risk_ev,objective=solve_min_risk,mu=mu, &
      covariance=cov,lb=lb,ub=ub,info=info)
   opt%max_iter = 3000
   opt%tol = 1.0e-10_dp
   call parmasolve(spec,port,opt)
   call assert_close(port%weights(1),0.8_dp,2.0e-4_dp,'portfolio weight 1')
   call assert_close(port%weights(2),0.2_dp,2.0e-4_dp,'portfolio weight 2')
   if (port%status /= 0) error stop 'portfolio status failed'

   print '(a)', 'test_optimizers: PASS'

contains

   function sphere(x) result(value)
      real(dp), intent(in) :: x(:)
      real(dp) :: value
      value = sum(x*x)
   end function sphere

   subroutine assert_close(actual,wanted,epsilonx,name)
      real(dp), intent(in) :: actual,wanted,epsilonx
      character(len=*), intent(in) :: name
      if (abs(actual-wanted) > epsilonx*max(1.0_dp,abs(wanted))) then
         write(*,'(a,2es24.14)') trim(name)//' failed: ',actual,wanted
         error stop 1
      end if
   end subroutine assert_close

end program test_optimizers
