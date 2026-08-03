program test_validation
   use clarabel
   implicit none
   real(dp) :: pd(2,2), ad(0,2), q(2), b(0), aeq_dense(1,2), anew_dense(1,2)
   type(csc_matrix) :: p, a, aeq, anew
   type(clarabel_cone), allocatable :: cones(:)
   type(clarabel_settings) :: settings
   type(clarabel_solution) :: sol
   type(clarabel_solver_type) :: solver
   integer :: code
   character(len=:), allocatable :: message

   pd = 0.0_dp
   pd(1,1) = 1.0_dp
   pd(2,2) = 1.0_dp
   p = csc_from_symmetric_upper(pd)
   a = csc_from_dense(ad)
   q = 0.0_dp
   allocate(cones(1))
   cones(1) = nonnegative_cone(1)
   call clarabel_solve_problem(p, q, a, b, cones, sol, code=code, message=message)
   if (code /= -16) error stop "cone-row mismatch was not rejected"

   deallocate(cones)
   allocate(cones(1))
   cones(1) = zero_cone(1)
   aeq_dense = reshape([1.0_dp, 0.0_dp], shape(aeq_dense))
   anew_dense = reshape([0.0_dp, 1.0_dp], shape(anew_dense))
   aeq = csc_from_dense(aeq_dense)
   anew = csc_from_dense(anew_dense)
   settings = default_clarabel_settings()
   settings%verbose = .false.
   settings%presolve_enable = .false.
   call solver%initialize(p, q, aeq, [0.0_dp], cones, settings, code, message)
   if (code /= 0) error stop message
   call solver%update(a=anew, code=code, message=message)
   if (code /= -3) error stop "changed A sparsity pattern was not rejected"
   call solver%update(q=[1.0_dp], code=code, message=message)
   if (code /= -4) error stop "wrong q length was not rejected"
   print *, "test_validation: PASS"
end program test_validation
