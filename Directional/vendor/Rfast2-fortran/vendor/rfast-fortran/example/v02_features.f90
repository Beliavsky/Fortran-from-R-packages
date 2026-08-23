program v02_features
   use rfast
   implicit none
   integer, parameter :: n = 12
   real(dp) :: x(n,2), y(n), z(7)
   integer :: i
   type(selection_result) :: sel
   type(el_result) :: el

   do i = 1, n
      x(i,1) = (real(i,dp) - 6.5_dp)/3.0_dp
      x(i,2) = sin(real(i,dp))
      y(i) = 1.0_dp + 2.0_dp*x(i,1) + 0.05_dp*cos(real(i,dp))
   end do

   sel = ompr(y,x,OMP_BIC,2.0_dp,.true.)
   print '(a,i0)', 'first OMP/BIC predictor: ', sel%selected(1)

   z = [-3.0_dp,-2.0_dp,-1.0_dp,0.0_dp,1.0_dp,2.0_dp,3.0_dp]
   el = el_test1(z,0.0_dp)
   print '(a,es14.6)', 'EL statistic at the sample mean: ', el%statistic
end program v02_features
