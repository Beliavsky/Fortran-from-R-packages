program example_semi_partial
   use ppcor, only : dp, ppcor_result, spcor, ppcor_pearson
   implicit none
   real(dp) :: x(30,3)
   type(ppcor_result) :: result
   integer :: i
   real(dp) :: t

   do i = 1, size(x,1)
      t = real(i,dp)
      x(i,3) = sin(0.2_dp*t)
      x(i,1) = 0.8_dp*x(i,3) + cos(0.73_dp*t)
      x(i,2) = 0.4_dp*x(i,1) - 0.3_dp*x(i,3) + sin(0.51_dp*t)
   end do

   call spcor(x, result, ppcor_pearson)
   print '(a,f10.6)', 'semi-partial r(1,2 | others from variable 1) = ', &
         result%estimate(1,2)
   print '(a,f10.6)', 'semi-partial r(2,1 | others from variable 2) = ', &
         result%estimate(2,1)
end program example_semi_partial
