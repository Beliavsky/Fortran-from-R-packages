program unconstrained_qp
   use clarabel
   implicit none
   real(dp) :: pd(3,3), ad(0,3), q(3), b(0)
   type(csc_matrix) :: p, a
   type(clarabel_cone), allocatable :: cones(:)
   type(clarabel_solution) :: sol
   integer :: code
   character(len=:), allocatable :: message

   pd = 0.0_dp
   pd(1,1) = 1.0_dp; pd(2,2) = 1.0_dp; pd(3,3) = 1.0_dp
   q = [1.0_dp, 2.0_dp, -3.0_dp]
   p = csc_from_symmetric_upper(pd)
   a = csc_from_dense(ad)
   allocate(cones(0))
   call clarabel_solve_problem(p, q, a, b, cones, sol, code=code, message=message)
   if (code /= 0) error stop message
   print '(a,3f12.6)', 'x = ', sol%x
   print '(a,a)', 'status = ', status_name(sol%status)
end program unconstrained_qp
