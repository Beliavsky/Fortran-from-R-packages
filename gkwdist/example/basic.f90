program basic
   use gkwdist
   implicit none
   real(dp) :: x,p
   real(dp),allocatable :: r(:)
   x=0.4_dp
   p=pgkw(x,1.7_dp,2.4_dp,1.3_dp,0.8_dp,1.2_dp)
   print '(a,f12.8)','GKw CDF at x=0.4: ',p
   print '(a,f12.8)','Recovered quantile: ',qgkw(p,1.7_dp,2.4_dp,1.3_dp,0.8_dp,1.2_dp)
   call seed_rng(123)
   r=rgkw(5,1.7_dp,2.4_dp,1.3_dp,0.8_dp,1.2_dp)
   print '(a,5f10.6)','Five draws: ',r
end program basic
