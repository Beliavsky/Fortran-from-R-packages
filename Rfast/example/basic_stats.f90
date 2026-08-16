program basic_stats
   use rfast
   implicit none
   real(dp) :: x(8), design(8,2), y(8)
   type(mle_result) :: fit
   type(regression_result) :: reg
   integer :: i

   x=[1.0_dp,2.0_dp,2.5_dp,3.0_dp,4.0_dp,5.0_dp,6.0_dp,7.0_dp]
   print '(a,f10.5)', 'mean:   ', mean_r(x)
   print '(a,f10.5)', 'median: ', median_r(x)
   print '(a,f10.5)', 'sd:     ', sqrt(variance_r(x))

   fit=normal_mle(x)
   print '(a,2f12.6)', 'normal MLE mean/variance: ', fit%param

   do i=1,8
      design(i,:)=[1.0_dp,real(i,dp)]
   end do
   y=2.0_dp+0.5_dp*design(:,2)
   reg=lmfit(design,y)
   print '(a,2f12.6)', 'linear-model coefficients: ', reg%beta
end program basic_stats
