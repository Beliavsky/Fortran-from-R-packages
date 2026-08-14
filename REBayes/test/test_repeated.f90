program test_repeated
   use rebayes_kinds, only : dp
   use rebayes_mixtures, only : mixture_fit, bivariate_mixture_fit
   use rebayes_repeated
   implicit none
   real(dp)::y(12),w(12),u(3),v(3)
   integer::id(12)
   type(mixture_fit)::gv
   type(bivariate_mixture_fit)::glv
   integer::i
   y=[-1.1_dp,-0.9_dp,-1.0_dp, 1.0_dp,1.2_dp,0.8_dp, -0.8_dp,-1.2_dp,-1.0_dp, 0.9_dp,1.1_dp,1.0_dp]
   id=[1,1,1,2,2,2,3,3,3,4,4,4];w=1.0_dp;u=[-1.0_dp,0.0_dp,1.0_dp];v=[0.02_dp,0.1_dp,0.5_dp]
   call wgvmix(y,id,w,v,gv)
   if(abs(sum(gv%mass)-1.0_dp)>1.0e-10_dp)error stop "wgvmix mass"
   call wglvmix(y,id,w,u,v,glv)
   if(abs(sum(glv%mass)-1.0_dp)>1.0e-10_dp)error stop "wglvmix mass"
   if(any(glv%g<=0.0_dp))error stop "wglvmix g"
   do i=1,size(glv%post_u)
      if(abs(glv%post_u(i))>1.5_dp)error stop "wglvmix posterior"
   end do
   print *,"test_repeated: PASS"
end program test_repeated
