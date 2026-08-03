program example_second_order_cone
   use clarabel
   implicit none
   real(dp) :: ad(5,3), q(3), b(5)
   type(csc_matrix) :: p, a
   type(clarabel_cone) :: cones(2)
   type(clarabel_solution) :: sol
   integer :: i, code
   character(len=:), allocatable :: message

   ! Minimize t subject to x1=1, x2=2, and sqrt(x1^2+x2^2) <= t.
   ad = 0.0_dp
   ad(1,2) = 1.0_dp
   ad(2,3) = 1.0_dp
   do i = 1, 3
      ad(i+2,i) = -1.0_dp
   end do
   b = [1.0_dp, 2.0_dp, 0.0_dp, 0.0_dp, 0.0_dp]
   q = [1.0_dp, 0.0_dp, 0.0_dp]
   p = csc_empty(3,3)
   a = csc_from_dense(ad)
   cones = [zero_cone(2), second_order_cone(3)]
   call clarabel_solve_problem(p, q, a, b, cones, sol, code=code, message=message)
   if (code /= 0) error stop message
   print '(a,3f12.6)', '[t,x1,x2] = ', sol%x
end program example_second_order_cone
