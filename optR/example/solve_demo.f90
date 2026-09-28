program solve_demo
   use optr, only : dp, optr_success, optr_solve
   implicit none

   real(dp), allocatable :: x(:)
   real(dp) :: a(3, 3)
   real(dp) :: b(3)
   integer :: status

   a = reshape([6.0_dp, -4.0_dp, 1.0_dp, -4.0_dp, 6.0_dp, -4.0_dp, &
      1.0_dp, -4.0_dp, 6.0_dp], [3, 3], order=[2, 1])
   b = [-14.0_dp, 36.0_dp, 6.0_dp]

   call optr_solve(a, b, 'LU', x, status)
   if (status /= optr_success) error stop 'optR solve failed'

   print '(a,3(1x,f10.6))', 'solution:', x
end program solve_demo
