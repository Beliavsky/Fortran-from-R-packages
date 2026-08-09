program qprog_example
   use coneproj
   implicit none
   real(dp) :: q(2,2),c(2),a(2,2),b(2)
   type(qprog_result) :: ans
   q=0.0_dp; q(1,1)=2.0_dp; q(2,2)=1.0_dp
   c=[-1.0_dp,3.0_dp]
   a=0.0_dp; a(1,1)=1.0_dp; a(2,2)=1.0_dp
   b=0.0_dp
   call qprog(q,c,a,b,ans)
   print '(a,2(f12.6,1x))','solution: ',ans%theta
   print '(a,es14.6)','objective: ',ans%objective
end program qprog_example
