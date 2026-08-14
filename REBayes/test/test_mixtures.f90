program test_mixtures
   use rebayes_kinds, only : dp
   use rebayes_mixtures
   implicit none
   real(dp)::x(8),v(3),sig(1)
   integer::bx(6),bk(6)
   real(dp)::bv(3)
   type(mixture_fit)::g,b
   x=[-1.2_dp,-0.9_dp,-1.1_dp,-0.8_dp,1.0_dp,1.1_dp,0.9_dp,1.2_dp]
   v=[-1.0_dp,0.0_dp,1.0_dp];sig=0.25_dp
   call glmix(x,v,sig,g)
   if(g%mass(2)>1.0e-3_dp)error stop "glmix middle mass"
   if(abs(g%mass(1)-0.5_dp)>2.0e-3_dp)error stop "glmix left mass"
   bx=[0,1,1,3,4,4];bk=4;bv=[0.1_dp,0.5_dp,0.9_dp]
   call bmix(bx,bk,bv,b)
   if(abs(sum(b%mass)-1.0_dp)>1.0e-12_dp)error stop "bmix mass"
   if(any(b%g<=0.0_dp))error stop "bmix g"
   print *,"test_mixtures: PASS"
end program test_mixtures
