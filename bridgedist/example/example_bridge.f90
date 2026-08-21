program example_bridge
   use bridgedist, only : dp, dbridge, pbridge, qbridge, bridge_variance
   implicit none
   real(dp) :: phi

   phi = 0.5_dp
   print '(a,f12.8)', 'density at 1 = ', dbridge(1.0_dp, phi)
   print '(a,f12.8)', 'cdf at 1     = ', pbridge(1.0_dp, phi)
   print '(a,f12.8)', 'q(.90)       = ', qbridge(0.9_dp, phi)
   print '(a,f12.8)', 'variance     = ', bridge_variance(phi)
end program example_bridge
