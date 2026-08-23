program test_special_families
   use lmomco, only : dp, lmomco_params, make_params, lmomco_pdf, lmomco_cdf, lmomco_quantile
   implicit none
   type(lmomco_params) :: p
   real(dp) :: x, f
   p=make_params('rice',[1.0_dp,0.8_dp])
   if(lmomco_pdf(1.2_dp,p)<=0.0_dp) error stop 1
   f=lmomco_cdf(1.5_dp,p)
   if(f<=0.0_dp .or. f>=1.0_dp) error stop 1
   x=lmomco_quantile(0.5_dp,p)
   if(abs(lmomco_cdf(x,p)-0.5_dp)>3.0e-4_dp) error stop 1

   p=make_params('kmu',[0.8_dp,1.4_dp])
   if(lmomco_pdf(1.0_dp,p)<=0.0_dp) error stop 1
   if(lmomco_cdf(1.0_dp,p)<=0.0_dp) error stop 1
   print '(a)', 'test_special_families: PASS'
end program test_special_families
