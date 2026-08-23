program basic
   use rfast2
   implicit none
   real(dp) :: x(5),a(5,2)
   type(pca_result) :: pc
   type(scalar_test_result) :: jb
   x=[1.0_dp,2.0_dp,3.0_dp,4.0_dp,5.0_dp]
   a(:,1)=x
   a(:,2)=[2.0_dp,1.0_dp,2.0_dp,1.0_dp,2.0_dp]
   print '(a,f8.4)', '25% quantile = ',quantile_rfast2(x,0.25_dp)
   jb=jarque_bera(x)
   print '(a,f8.4)', 'JB statistic  = ',jb%statistic
   pc=pca(a,k=2,vectors=.true.)
   print '(a,2f10.5)', 'PCA values    = ',pc%values
end program basic
