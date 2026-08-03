program equality_qp
   use clarabel
   implicit none
   real(dp) :: pd(2,2), ad(1,2), q(2), b(1)
   type(csc_matrix) :: p, a
   type(clarabel_cone) :: cones(1)
   type(clarabel_solution) :: sol
   integer :: code
   character(len=:), allocatable :: message

   pd = 0.0_dp
   pd(1,1) = 1.0_dp; pd(2,2) = 1.0_dp
   ad = reshape([1.0_dp, 1.0_dp], shape(ad))
   q = 0.0_dp; b = 1.0_dp
   p = csc_from_symmetric_upper(pd)
   a = csc_from_dense(ad)
   cones(1) = zero_cone(1)
   call clarabel_solve_problem(p, q, a, b, cones, sol, code=code, message=message)
   if (code /= 0) error stop message
   print '(a,2f12.6)', 'x = ', sol%x
   print '(a,f12.6)', 'objective = ', sol%obj_val
end program equality_qp
