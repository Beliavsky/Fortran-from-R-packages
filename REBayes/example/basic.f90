program basic
   use rebayes_kinds, only : dp
   use rebayes_mixtures, only : mixture_fit, glmix
   implicit none
   real(dp)::x(10),grid(5),sigma(1)
   type(mixture_fit)::fit
   x=[-1.2_dp,-1.0_dp,-0.8_dp,-1.1_dp,-0.9_dp,0.8_dp,1.0_dp,1.2_dp,0.9_dp,1.1_dp]
   grid=[-1.0_dp,-0.5_dp,0.0_dp,0.5_dp,1.0_dp];sigma=0.25_dp
   call glmix(x,grid,sigma,fit)
   print '(a,5f10.6)',"mixing weights: ",fit%mass
   print '(a,f12.6)',"log likelihood: ",fit%loglik
end program basic
