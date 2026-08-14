program test_medde
   use rebayes_kinds, only : dp
   use rebayes_medde
   implicit none
   real(dp)::x(8),grid(21),q(3),out(3)
   type(medde_result)::r
   integer::i
   x=[-1.2_dp,-0.8_dp,-0.5_dp,-0.1_dp,0.1_dp,0.5_dp,0.8_dp,1.2_dp]
   do i=1,21;grid(i)=-1.5_dp+3.0_dp*real(i-1,dp)/20.0_dp;end do
   call medde_fit(x,grid,0.2_dp,1.0_dp,1,r)
   if(any(r%y< -1.0e-10_dp))error stop "medde nonnegative"
   q=[0.25_dp,0.5_dp,0.75_dp];call medde_quantiles(r,q,out)
   if(.not.(out(1)<=out(2).and.out(2)<=out(3)))error stop "medde quantiles"
   print *,"test_medde: PASS"
end program test_medde
