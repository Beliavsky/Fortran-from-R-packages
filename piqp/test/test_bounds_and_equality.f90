program test_bounds_and_equality
   use piqp
   implicit none
   real(dp) :: p1(1,1), c1(1), g1(1,1), hl1(1), hu1(1)
   real(dp) :: p2(2,2), c2(2), a2(1,2), b2(1)
   type(piqp_result_type) :: r
   p1(1,1)=2.0_dp; c1=0.0_dp; g1(1,1)=1.0_dp; hl1=1.0_dp; hu1=2.0_dp
   call solve_piqp(p1,c1,r,gmat=g1,h_l=hl1,h_u=hu1)
   if(r%info%status/=PIQP_SOLVED .or. abs(r%x(1)-1.0_dp)>1.0e-6_dp) error stop 'two sided bound'
   if(abs(r%z_l(1)-2.0_dp)>2.0e-5_dp) error stop 'lower dual'
   p2=0.0_dp; p2(1,1)=1.0_dp; p2(2,2)=1.0_dp
   c2=[-2.0_dp,-5.0_dp]; a2=reshape([1.0_dp,1.0_dp],[1,2]); b2=3.0_dp
   call solve_piqp(p2,c2,r,amat=a2,b=b2)
   if(r%info%status/=PIQP_SOLVED) error stop 'equality status'
   if(maxval(abs(r%x-[0.0_dp,3.0_dp]))>2.0e-6_dp) error stop 'equality x'
   if(abs(r%y(1)-2.0_dp)>2.0e-6_dp) error stop 'equality y'
   print *, 'test_bounds_and_equality: PASS'
end program test_bounds_and_equality
