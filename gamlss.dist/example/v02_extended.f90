program v02_extended
   use gamlss_dist
   implicit none
   real(dp) :: x,p

   x=qST3(0.90_dp,mu=0.0_dp,sigma=1.0_dp,nu=1.4_dp,tau=4.5_dp)
   p=dSICHEL(4.0_dp,mu=2.3_dp,sigma=0.6_dp,nu=-0.4_dp)

   print '(a,f12.6)', 'ST3 90% quantile:       ',x
   print '(a,f12.6)', 'SICHEL P(Y = 4):        ',p
   print '(a,f12.6)', 'SHASH CDF at x = 0.7:   ', &
      pSHASH(0.7_dp,mu=0.2_dp,sigma=1.3_dp,nu=0.8_dp,tau=1.2_dp)
end program v02_extended
