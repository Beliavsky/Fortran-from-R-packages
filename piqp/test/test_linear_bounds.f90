program test_linear_bounds
   use piqp
   implicit none
   real(dp) :: c(1), xl(1), xu(1)
   type(piqp_result_type) :: r
   c=-1.0_dp; xl=0.0_dp; xu=1.0_dp
   call solve_piqp(c=c,result=r,x_l=xl,x_u=xu)
   if(r%info%status/=PIQP_SOLVED) error stop 'linear bounds status'
   if(abs(r%x(1)-1.0_dp)>1.0e-6_dp) error stop 'linear bounds x'
   if(abs(r%z_bu(1)-1.0_dp)>2.0e-5_dp) error stop 'linear bounds dual'
   print *, 'test_linear_bounds: PASS'
end program test_linear_bounds
