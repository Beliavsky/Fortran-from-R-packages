program basic_matrixdist
   use r_compat, only: dp
   use matrixdist, only: ph_density, ph_cdf, ph_mean, rphasetype
   implicit none
   real(dp) :: alpha(2), s(2,2)
   real(dp), allocatable :: x(:)
   alpha=[0.7_dp,0.3_dp]
   s=reshape([-3.0_dp,0.5_dp,1.0_dp,-2.0_dp],[2,2])
   print '(a,f12.8)', 'density at 0.8 = ', ph_density(0.8_dp,alpha,s)
   print '(a,f12.8)', 'cdf at 0.8     = ', ph_cdf(0.8_dp,alpha,s)
   print '(a,f12.8)', 'mean           = ', ph_mean(alpha,s)
   x=rphasetype(5,alpha,s)
   print '(a,5f10.5)', 'five draws: ', x
end program basic_matrixdist
