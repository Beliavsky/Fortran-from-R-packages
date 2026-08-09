program test_model_update
   use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_positive_inf
   use piqp
   implicit none
   real(dp) :: p(2,2), c(2), a(1,2), b(1), g(2,2), hu(2), xl(2), xu(2), inf
   type(piqp_model_type) :: model
   integer :: info
   inf=ieee_value(0.0_dp,ieee_positive_inf)
   p=reshape([6.0_dp,0.0_dp,0.0_dp,4.0_dp],[2,2]); c=[-1.0_dp,-4.0_dp]
   a=reshape([1.0_dp,-2.0_dp],[1,2]); b=0.0_dp
   g=reshape([1.0_dp,-1.0_dp,0.0_dp,0.0_dp],[2,2]); hu=[1.0_dp,1.0_dp]
   xl=[-inf,-1.0_dp]; xu=[inf,1.0_dp]
   call model%setup(p,c,a,b,g,h_u=hu,x_l=xl,x_u=xu)
   call model%solve()
   p(1,1)=8.0_dp; a=reshape([1.0_dp,-3.0_dp],[1,2]); hu=[2.0_dp,1.0_dp]; xu=[inf,2.0_dp]
   call model%update(pmat=p,amat=a,h_u=hu,x_u=xu,info=info)
   if(info/=0) error stop 'model update info'
   call model%solve()
   if(model%result%info%status/=PIQP_SOLVED) error stop 'model update status'
   if(maxval(abs(model%result%x-[0.2763157_dp,0.0921052_dp]))>1.0e-6_dp) error stop 'model update x'
   if(abs(model%result%y(1)+1.2105263_dp)>1.0e-6_dp) error stop 'model update y'
   print *, 'test_model_update: PASS'
end program test_model_update
