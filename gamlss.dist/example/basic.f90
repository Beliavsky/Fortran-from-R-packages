program basic
   use gamlss_dist
   implicit none
   real(dp) :: x, probability

   x = qBCT(0.95_dp, mu=10.0_dp, sigma=0.2_dp, nu=0.5_dp, tau=5.0_dp)
   probability = pNBI(4.0_dp, mu=2.7_dp, sigma=0.45_dp)

   print '(a,f10.5)', 'BCT 95% quantile: ', x
   print '(a,f10.5)', 'NBI P(Y <= 4):    ', probability
end program basic
