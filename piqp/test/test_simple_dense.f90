program test_simple_dense
   use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_positive_inf
   use piqp
   implicit none
   real(dp) :: p(2,2), c(2), a(1,2), b(1), g(2,2), hu(2), xl(2), xu(2), inf
   type(piqp_result_type) :: r
   inf=ieee_value(0.0_dp,ieee_positive_inf)
   p=reshape([6.0_dp,0.0_dp,0.0_dp,4.0_dp],[2,2]); c=[-1.0_dp,-4.0_dp]
   a=reshape([1.0_dp,-2.0_dp],[1,2]); b=0.0_dp
   g=reshape([1.0_dp,-1.0_dp,0.0_dp,0.0_dp],[2,2]); hu=[1.0_dp,1.0_dp]
   xl=[-inf,-1.0_dp]; xu=[inf,1.0_dp]
   call solve_piqp(p,c,r,a,b,g,h_u=hu,x_l=xl,x_u=xu)
   if(r%info%status/=PIQP_SOLVED) error stop 'simple dense status'
   if(maxval(abs(r%x-[0.4285714_dp,0.2142857_dp]))>1.0e-6_dp) error stop 'simple dense x'
   if(abs(r%y(1)+1.5714286_dp)>1.0e-6_dp) error stop 'simple dense y'
   if(maxval(abs(r%z_l))>1.0e-6_dp .or. maxval(abs(r%z_bl))>1.0e-6_dp) error stop 'simple dense duals'
   print *, 'test_simple_dense: PASS'
end program test_simple_dense
