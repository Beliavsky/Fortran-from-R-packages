program basic
   use ccd, only : dp, i8, dcc, pcc, qcc, cc_fit_result, cc_mle0
   implicit none
   real(dp) :: y(9)
   type(cc_fit_result) :: fit
   y = real([-4,-2,-1,0,0,1,1,3,5], dp)
   print '(a,f12.8)', 'P(X=0), lambda=1.5: ', dcc(0.0_dp, 0.0_dp, 1.5_dp)
   print '(a,f12.8)', 'P(X<=2):            ', pcc(2_i8, 0.0_dp, 1.5_dp)
   print '(a,i0)', 'median quantile:       ', qcc(0.5_dp, 0.0_dp, 1.5_dp)
   fit = cc_mle0(y)
   print '(a,f12.6)', 'zero-location MLE lambda: ', fit%lambda
end program basic
