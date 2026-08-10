program test_qp_pnnqp
   use lsei
   implicit none
   real(dp) :: q(2,2),p(2),c(1,2),d(1),lo(2)
   real(dp) :: a3(3,3),b3(3),q3(3,3),p3(3)
   type(ls_result) :: r1,r2
   q=0.0_dp; q(1,1)=1.0_dp; q(2,2)=1.0_dp
   p=[-2.0_dp,-3.0_dp]; c(1,:)=[1.0_dp,1.0_dp]; d=[4.0_dp]; lo=0.0_dp
   call qp_solve(q,p,r1,c=c,d=d,lower=lo)
   if (.not.r1%succeeded()) error stop 1
   if (maxval(abs(r1%x-[1.5_dp,2.5_dp]))>2e-10_dp) error stop 2
   if (abs(r1%objective+6.25_dp)>2e-10_dp) error stop 3
   a3=0.0_dp; a3(1,1)=1.0_dp; a3(2,2)=1.0_dp; a3(3,3)=1.0_dp
   b3=[.2_dp,.3_dp,.7_dp]; q3=matmul(transpose(a3),a3); p3=-matmul(transpose(a3),b3)
   call pnnls_solve(a3,b3,1,r1,sum_value=1.0_dp)
   call pnnqp_solve(q3,p3,1,r2,sum_value=1.0_dp)
   if (maxval(abs(r1%x-[.4_dp,.1_dp,.9_dp]))>2e-12_dp) error stop 4
   if (maxval(abs(r1%x-r2%x))>2e-12_dp) error stop 5
   if (abs(sum(r1%x(2:))-1.0_dp)>2e-12_dp) error stop 6
   print *, 'PASS test_qp_pnnqp'
end program
