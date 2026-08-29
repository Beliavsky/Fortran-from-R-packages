program dpq_demo
   use dpq, only: dp, pnchisq, qnchisq, pnt_r, qnt_r
   implicit none
   real(dp) :: p, q

   p = pnchisq(7.25_dp, 4.5_dp, 3.2_dp)
   q = qnchisq(0.73_dp, 4.5_dp, 3.2_dp)
   print '(a,f18.12)', 'P[noncentral chi-square <= 7.25] = ', p
   print '(a,f18.12)', 'noncentral chi-square q(0.73)    = ', q

   p = pnt_r(1.2_dp, 7.5_dp, 0.8_dp)
   q = qnt_r(0.73_dp, 7.5_dp, 0.8_dp)
   print '(a,f18.12)', 'P[noncentral t <= 1.2]           = ', p
   print '(a,f18.12)', 'noncentral t q(0.73)             = ', q
end program dpq_demo
