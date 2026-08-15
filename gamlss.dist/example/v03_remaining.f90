program v03_remaining
   use gamlss_dist
   implicit none
   real(dp) :: q

   q=qGT(0.90_dp,mu=0.0_dp,sigma=1.0_dp,nu=4.0_dp,tau=1.5_dp)
   print '(a,f12.6)', 'GT 90% quantile:             ',q
   print '(a,f12.6)', 'Double-binomial P(Y=3):     ',dDBI(3.0_dp,0.35_dp,0.6_dp,8.0_dp)
   print '(a,f12.6)', 'Zero-adjusted PIG P(Y=0):   ',dZAPIG(0.0_dp,2.3_dp,0.6_dp,0.2_dp)
   print '(a,f12.6)', 'Pareto-II 75% quantile:     ',qPARETO2(0.75_dp,1.4_dp,0.6_dp)
end program v03_remaining
