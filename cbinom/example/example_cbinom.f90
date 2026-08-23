program example_cbinom
   use cbinom, only : dp, dcbinom, pcbinom, qcbinom, rcbinom
   implicit none
   real(dp) :: x(5)

   print '(a,f12.8)', 'density:  ', dcbinom(4.0_dp, 20.0_dp, 0.2_dp)
   print '(a,f12.8)', 'cdf:      ', pcbinom(4.0_dp, 20.0_dp, 0.2_dp)
   print '(a,f12.8)', 'quantile: ', qcbinom(0.5_dp, 20.0_dp, 0.2_dp)
   call rcbinom(x, 20.0_dp, 0.2_dp)
   print '(a,5f10.5)', 'random:   ', x
end program example_cbinom
