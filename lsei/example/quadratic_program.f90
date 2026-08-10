program quadratic_program
   use lsei
   implicit none
   real(dp) :: q(2,2),p(2),c(1,2),d(1),lo(2)
   type(ls_result) :: r
   q=0.0_dp; q(1,1)=1.0_dp; q(2,2)=1.0_dp
   p=[-2.0_dp,-3.0_dp]; c(1,:)=[1.0_dp,1.0_dp]; d=[4.0_dp]; lo=0.0_dp
   call qp_solve(q,p,r,c=c,d=d,lower=lo)
   print '(a,2f14.8)', 'QP solution: ',r%x
   print '(a,f14.8)', 'objective: ',r%objective
end program
