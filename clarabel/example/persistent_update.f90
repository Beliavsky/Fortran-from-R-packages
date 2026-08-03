program persistent_update
   use clarabel
   implicit none
   real(dp) :: pd(2,2), ad(0,2), q(2), b(0)
   type(csc_matrix) :: p, a
   type(clarabel_cone), allocatable :: cones(:)
   type(clarabel_settings) :: settings
   type(clarabel_solution) :: sol
   type(clarabel_solver_type) :: solver
   integer :: code
   character(len=:), allocatable :: message

   pd = 0.0_dp
   pd(1,1) = 1.0_dp; pd(2,2) = 1.0_dp
   q = [1.0_dp, -2.0_dp]
   p = csc_from_symmetric_upper(pd)
   a = csc_from_dense(ad)
   allocate(cones(0))
   settings = default_clarabel_settings()
   settings%verbose = .false.
   settings%presolve_enable = .false.
   settings%chordal_decomposition_enable = .false.
   settings%input_sparse_dropzeros = .false.
   call solver%initialize(p, q, a, b, cones, settings, code, message)
   if (code /= 0) error stop message
   call solver%solve(sol, code, message)
   if (code /= 0) error stop message
   print '(a,2f10.4)', 'first x = ', sol%x
   call solver%update(q=[-3.0_dp, 4.0_dp], code=code, message=message)
   if (code /= 0) error stop message
   call solver%solve(sol, code, message)
   if (code /= 0) error stop message
   print '(a,2f10.4)', 'updated x = ', sol%x
end program persistent_update
