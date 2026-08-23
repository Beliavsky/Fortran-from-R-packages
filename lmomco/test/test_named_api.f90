program test_named_api
   use lmomco, only : dp, lmomco_params, make_params, pdfgev, cdfgev, quagev, pargev
   implicit none
   type(lmomco_params) :: p, fit
   real(dp) :: x, f, dens, l(3)
   p=make_params('gev',[0.0_dp,1.0_dp,0.1_dp])
   x=quagev(0.73_dp,p)
   f=cdfgev(x,p)
   dens=pdfgev(x,p)
   if(abs(f-0.73_dp)>1.0e-11_dp) error stop 1
   if(dens<=0.0_dp) error stop 1
   l=[0.2_dp,0.8_dp,0.12_dp]
   fit=pargev(l)
   if(fit%npar/=3 .or. fit%p(2)<=0.0_dp) error stop 1
   print '(a)', 'test_named_api: PASS'
end program test_named_api
