program demo_suppdists
   use suppdists
   implicit none
   type(johnson_parms) :: jp
   real(dp), allocatable :: scores(:)
   jp%gamma=0.3_dp; jp%delta=1.2_dp; jp%xi=0.0_dp; jp%lambda=1.0_dp; jp%family=johnson_su
   print '(a,f12.8)', 'Inverse Gaussian CDF at 1.5: ', pinvgauss(1.5_dp,2.0_dp,3.0_dp)
   print '(a,f12.8)', 'Kendall P(tau <= 0), N=8:   ', pkendall(0.0_dp,8)
   print '(a,f12.8)', 'Johnson SU 90% quantile:      ', qjohnson(0.9_dp,jp)
   print '(a,f12.8)', 'Classic hypergeom P(X=3):     ', dghyper(3,5.0_dp,7.0_dp,20.0_dp)
   call norm_order(7,scores)
   print '(a,7f9.4)', 'Expected normal scores: ', scores
end program demo_suppdists
