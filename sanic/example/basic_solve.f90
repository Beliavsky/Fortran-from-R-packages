program basic_solve
   use sanic, only: dp, solve2_dense
   implicit none

   real(dp) :: a(2, 2), b(2, 1)
   real(dp), allocatable :: x(:, :)
   integer :: info

   a = reshape([4.0_dp, 1.0_dp, 1.0_dp, 3.0_dp], [2, 2])
   b(:, 1) = [6.0_dp, 7.0_dp]
   call solve2_dense(a, x, info, b=b)
   if (info /= 0) error stop 'solve failed'
   write(*, '(a,2(1x,f10.6))') 'x =', x(:, 1)
end program basic_solve
