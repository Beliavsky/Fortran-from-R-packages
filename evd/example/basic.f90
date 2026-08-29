program basic
   use evd, only : dp, pgev, qgev, abvlog
   implicit none
   real(dp) :: p, q
   p=0.99_dp
   q=qgev(p,0.0_dp,1.0_dp,0.1_dp)
   print '(a,f10.5)', 'GEV 0.99 quantile: ',q
   print '(a,f10.5)', 'CDF at quantile:  ',pgev(q,0.0_dp,1.0_dp,0.1_dp)
   print '(a,f10.5)', 'A(0.5), dep=.7:  ',abvlog(0.5_dp,0.7_dp)
end program
