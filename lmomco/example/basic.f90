program basic
   use lmomco, only : dp, lmomco_params, make_params, lmomco_pdf, lmomco_cdf, lmomco_quantile
   implicit none
   type(lmomco_params) :: par
   par=make_params('gev',[10.0_dp,2.0_dp,0.1_dp])
   print '(a,f12.6)', 'pdf(11) = ',lmomco_pdf(11.0_dp,par)
   print '(a,f12.6)', 'cdf(11) = ',lmomco_cdf(11.0_dp,par)
   print '(a,f12.6)', 'q(.99)  = ',lmomco_quantile(0.99_dp,par)
end program basic
