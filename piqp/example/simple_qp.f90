program simple_qp
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
   print '(a,i0)', 'status = ',r%info%status
   print '(a,*(f12.8,1x))','x = ',r%x
   print '(a,*(f12.8,1x))','y = ',r%y
end program simple_qp
