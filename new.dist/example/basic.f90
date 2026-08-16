program basic
   use new_dist
   implicit none
   print '(a,f12.8)', 'Lindley density at 1: ', dLd(1.0_dp,2.0_dp)
   print '(a,f12.8)', 'Maxwell CDF at 1:   ', pmd(1.0_dp,2.0_dp)
   print '(a,f12.8)', 'Kumaraswamy q(.8):  ', qkd(0.8_dp,2.0_dp,3.0_dp)
   print '(a,i0)',    'Weighted geom q(.9):', qwgd(0.9_dp,0.2_dp,3.0_dp)
end program basic
