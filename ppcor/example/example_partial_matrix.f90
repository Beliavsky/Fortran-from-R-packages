program example_partial_matrix
   use ppcor, only : dp, ppcor_result, pcor, ppcor_pearson
   implicit none
   real(dp) :: x(25,4)
   type(ppcor_result) :: result
   integer :: i
   real(dp) :: t

   do i = 1, size(x,1)
      t = real(i,dp)
      x(i,1) = sin(0.23_dp*t) + 0.4_dp*cos(0.71_dp*t)
      x(i,2) = 0.5_dp*x(i,1) + cos(0.19_dp*t) + 0.2_dp*sin(0.83_dp*t)
      x(i,3) = sin(0.17_dp*t)
      x(i,4) = cos(0.29_dp*t) + 0.01_dp*t
   end do

   call pcor(x, result, ppcor_pearson)
   print '(a)', 'Pearson partial-correlation matrix:'
   do i = 1, size(x,2)
      print '(*(f11.6,1x))', result%estimate(i,:)
   end do
end program example_partial_matrix
