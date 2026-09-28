program basic
   use rmumps_api, only : dp, rmumps_matrix, rmumps_init_triplet, rmumps_solve
   implicit none

   type(rmumps_matrix) :: a
   integer :: info
   integer :: i(4)
   integer :: j(4)
   real(dp) :: v(4)
   real(dp) :: b(2)
   real(dp) :: x(2)

   i = [1, 1, 2, 2]
   j = [1, 2, 1, 2]
   v = [3.0_dp, 1.0_dp, 1.0_dp, 2.0_dp]
   b = [5.0_dp, 5.0_dp]

   call rmumps_init_triplet(a, i, j, v, 2, info=info)
   if (info /= 0) error stop 'rmumps_init_triplet failed'
   call rmumps_solve(a, b, x, info)
   if (info /= 0) error stop 'rmumps_solve failed'

   print '(a,2(1x,f8.4))', 'solution:', x
end program basic
