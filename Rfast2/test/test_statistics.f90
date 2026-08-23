program test_statistics
   use rfast2
   implicit none
   type(km_result) :: kmr
   real(dp) :: meta(7)
   type(silhouette_result) :: sr
   type(scalar_test_result) :: tr
   real(dp) :: x(5),y(5),mat(4,2)
   integer :: status(5),cl(4)

   x = [1.0_dp,2.0_dp,3.0_dp,4.0_dp,5.0_dp]
   y = 2.0_dp*x+1.0_dp
   tr = cor_test_pearson(x,y)
   if (abs(tr%estimate-1.0_dp) > 1.0e-12_dp) error stop 1
   status = [1,0,1,1,0]
   kmr = kaplan_meier([1.0_dp,2.0_dp,3.0_dp,4.0_dp,5.0_dp],status)
   if (size(kmr%time) /= 3) error stop 2
   if (abs(kmr%survival(1)-0.8_dp) > 1.0e-12_dp) error stop 3
   meta = wls_meta([0.2_dp,0.4_dp,0.1_dp],[0.01_dp,0.04_dp,0.01_dp])
   if (abs(meta(1)-0.1777777777777778_dp) > 1.0e-12_dp) error stop 4
   mat = reshape([0.0_dp,0.1_dp,5.0_dp,5.1_dp,0.0_dp,0.1_dp,5.0_dp,5.1_dp],[4,2])
   cl = [1,1,2,2]
   sr = silhouette_euclidean(mat,cl)
   if (minval(sr%value) < 0.95_dp) error stop 5
   print '(a)', 'test_statistics: PASS'
end program test_statistics
