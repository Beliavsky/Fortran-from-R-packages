program demo
   use SpatialExtremes, only: dp,pgev,qgev,COV_POWEREXP,covariance_matrix,extremal_coefficient_schlather
   implicit none
   real(dp) :: coord(3,2),c(3,3),p,x
   coord=reshape([0.0_dp,0.0_dp, 1.0_dp,0.0_dp, 0.0_dp,1.0_dp],[3,2],order=[2,1])
   p=pgev(1.0_dp,0.0_dp,1.0_dp,0.1_dp)
   x=qgev(p,0.0_dp,1.0_dp,0.1_dp)
   c=covariance_matrix(coord,COV_POWEREXP,0.1_dp,0.9_dp,2.0_dp,1.5_dp)
   print '(a,f10.6)', 'p/q inversion: ',x
   print '(a,f10.6)', 'theta(h=1): ',extremal_coefficient_schlather(c(1,2))
end program demo
