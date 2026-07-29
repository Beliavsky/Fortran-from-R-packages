! SPDX-License-Identifier: GPL-3.0-or-later
program test_linear_quadratic
   use cccp
   implicit none
   real(dp) :: p(2,2), q2(2), a(1,2), b(1), ql(2), g(4,2), h(4)
   type(cone_constraint) :: c(1)
   type(cccp_solution) :: sol
   type(dqp_problem) :: qprob
   integer :: fails
   fails=0

   p=2.0_dp*reshape([2.0_dp,0.5_dp,0.5_dp,1.0_dp],[2,2])
   q2=[1.0_dp,1.0_dp];a=reshape([1.0_dp,1.0_dp],[1,2]);b=1.0_dp
   call solve_qp(p,q2,a,b,sol=sol)
   call check(trim(sol%status)=='optimal','QP status',fails)
   call check(maxval(abs(sol%x-[0.25_dp,0.75_dp]))<5.0e-5_dp,'QP solution',fails)
   call check(abs(sol%state%pobj-1.875_dp)<5.0e-5_dp,'QP objective',fails)
   qprob=dqp(p,q2,a,b)
   sol=cps(qprob)
   call check(maxval(abs(getx(sol)-[0.25_dp,0.75_dp]))<5.0e-5_dp,'DQP/CPS compatibility',fails)

   ql=[-4.0_dp,-5.0_dp]
   g=reshape([2.0_dp,1.0_dp,-1.0_dp,0.0_dp, 1.0_dp,2.0_dp,0.0_dp,-1.0_dp],[4,2])
   h=[3.0_dp,3.0_dp,0.0_dp,0.0_dp]
   c(1)=nnoc(g,h)
   call solve_lp(ql,cones=c,sol=sol)
   call check(trim(sol%status)=='optimal','LP status',fails)
   call check(maxval(abs(matmul(g,sol%x)-h))<10.0_dp .or. all(matmul(g,sol%x)<=h+1e-5_dp),'LP feasible',fails)
   call check(abs(sol%state%pobj+9.0_dp)<2.0e-4_dp,'LP objective',fails)

   if(fails>0)error stop 1
   print '(a)','test_linear_quadratic: PASS'
contains
   subroutine check(ok,name,fails)
      logical,intent(in)::ok
      character(len=*),intent(in)::name
      integer,intent(inout)::fails
      if(.not.ok)then
         print '(a,a)','FAIL: ',name;fails=fails+1
      end if
   end subroutine check
end program test_linear_quadratic
