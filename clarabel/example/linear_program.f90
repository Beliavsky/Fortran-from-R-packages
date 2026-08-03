program linear_program
   use clarabel
   implicit none
   real(dp) :: ad(6,3), q(3), b(6)
   type(csc_matrix) :: p, a
   type(clarabel_cone) :: cones(2)
   type(clarabel_solution) :: sol
   integer :: i, code
   character(len=:), allocatable :: message

   ad = 0.0_dp
   do i = 1, 3
      ad(i,i) = 2.0_dp
      ad(i+3,i) = -2.0_dp
   end do
   q = [3.0_dp, -2.0_dp, 1.0_dp]
   b = 1.0_dp
   p = csc_empty(3,3)
   a = csc_from_dense(ad)
   cones = [nonnegative_cone(3), nonnegative_cone(3)]
   call clarabel_solve_problem(p, q, a, b, cones, sol, code=code, message=message)
   if (code /= 0) error stop message
   print '(a,3f12.6)', 'x = ', sol%x
   print '(a,f12.6)', 'objective = ', sol%obj_val
end program linear_program
