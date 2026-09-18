program interpolate_sine
   use stinepack_api, only : dp, stinterp, stinterp_result
   implicit none

   integer :: i
   type(stinterp_result) :: fit
   real(dp), parameter :: pi = acos(-1.0_dp)
   real(dp) :: x(13)
   real(dp) :: xout(25)
   real(dp) :: y(13)

   do i = 1, size(x)
      x(i) = real(i - 1, dp) * pi / 6.0_dp
      y(i) = sin(x(i))
   end do
   do i = 1, size(xout)
      xout(i) = real(i - 1, dp) * pi / 12.0_dp
   end do

   fit = stinterp(x, y, xout)
   if (fit%status /= 0) error stop fit%message

   print '(a)', '       x              sin(x)          stinterp(x)'
   do i = 1, size(xout)
      print '(3f16.8)', xout(i), sin(xout(i)), fit%y(i)
   end do
end program interpolate_sine
