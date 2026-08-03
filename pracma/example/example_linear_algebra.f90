! SPDX-License-Identifier: GPL-3.0-or-later
program example_linear_algebra
   use pracma
   implicit none
   real(dp) :: a(3,3), b(3), x(3)
   complex(dp), allocatable :: r(:)
   integer :: status

   a=reshape([4.0_dp,1.0_dp,1.0_dp, &
              1.0_dp,3.0_dp,0.0_dp, &
              1.0_dp,0.0_dp,2.0_dp],[3,3])
   b=[1.0_dp,2.0_dp,3.0_dp]
   call solve_linear(a,b,x,status)
   r=roots([1.0_dp,-6.0_dp,11.0_dp,-6.0_dp])

   print '(a,3f12.6)','solution: ',x
   print '(a,f12.6)','determinant: ',determinant(a)
   print '(a,3f12.6)','polynomial roots: ',real(r,dp)
end program example_linear_algebra
