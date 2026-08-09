program test_status_settings
   use piqp
   implicit none
   real(dp) :: p(1,1), c(1), xl(1), xu(1)
   type(piqp_result_type) :: r
   type(piqp_settings_type) :: s
   p(1,1)=1.0_dp; c=0.0_dp
   xl=2.0_dp; xu=1.0_dp
   call solve_piqp(p,c,r,x_l=xl,x_u=xu)
   if(r%info%status/=PIQP_PRIMAL_INFEASIBLE) error stop 'infeasible bounds status'
   s%tau=2.0_dp; xl=-1.0_dp; xu=1.0_dp
   call solve_piqp(p,c,r,x_l=xl,x_u=xu,settings=s)
   if(r%info%status/=PIQP_INVALID_SETTINGS) error stop 'invalid settings status'
   if(index(status_description(PIQP_SOLVED),'Solver solved')/=1) error stop 'status description'
   print *, 'test_status_settings: PASS'
end program test_status_settings
