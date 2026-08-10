program test_bounds_status
   use lsei
   implicit none
   real(dp) :: a(2,2),b(2),lo(2),up(2),e(1,2),f(1)
   type(ls_result) :: r
   a=0.0_dp; a(1,1)=1.0_dp; a(2,2)=1.0_dp; b=[-2.0_dp,3.0_dp]
   lo=[0.0_dp,0.0_dp]; up=[1.0_dp,2.0_dp]
   call lsi_solve(a,b,res=r,lower=lo,upper=up)
   if (.not.r%succeeded()) error stop 1
   if (maxval(abs(r%x-[0.0_dp,2.0_dp]))>2e-10_dp) error stop 2
   e(1,:)=[0.0_dp,0.0_dp]; f=[1.0_dp]
   call ldp_solve(e,f,r)
   if (r%mode/=LSEI_INFEASIBLE) error stop 3
   print *, 'PASS test_bounds_status'
end program
