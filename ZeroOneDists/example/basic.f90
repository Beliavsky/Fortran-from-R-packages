program basic
   use zero_one_dists
   implicit none
   real(dp) :: x
   x=0.37_dp
   print '(a,f12.8)','BER density = ',dber(x,0.42_dp,3.7_dp,0.18_dp)
   print '(a,f12.8)','BER cdf     = ',pber(x,0.42_dp,3.7_dp,0.18_dp)
   print '(a,f12.8)','BER q(.73)  = ',qber(0.73_dp,0.42_dp,3.7_dp,0.18_dp)
end program basic
