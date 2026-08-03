program test_backend_qp
   use clarabel
   implicit none
   real(dp) :: pd(2,2), ad(6,2), q(2), b(6), target(2)
   type(csc_matrix) :: p, a
   type(clarabel_cone) :: cones(2)
   type(clarabel_settings) :: settings
   type(clarabel_solution) :: sol
   integer :: code
   character(len=:), allocatable :: message

   pd = reshape([4.0_dp, 1.0_dp, 1.0_dp, 2.0_dp], shape(pd))
   ad = reshape([-1.0_dp, -1.0_dp, 0.0_dp, 1.0_dp, 1.0_dp, 0.0_dp, &
                 -1.0_dp, 0.0_dp, -1.0_dp, 1.0_dp, 0.0_dp, 1.0_dp], shape(ad))
   q = [1.0_dp, 1.0_dp]
   b = [-1.0_dp, 0.0_dp, 0.0_dp, 1.0_dp, 0.7_dp, 0.7_dp]
   target = [0.3_dp, 0.7_dp]
   p = csc_from_symmetric_upper(pd)
   a = csc_from_dense(ad)
   cones = [nonnegative_cone(3), nonnegative_cone(3)]
   settings = default_clarabel_settings(); settings%verbose = .false.
   call clarabel_solve_problem(p, q, a, b, cones, sol, settings, code, message)
   if (code /= 0) error stop message
   if (.not. sol%solved()) error stop "QP status"
   if (maxval(abs(sol%x - target)) > 1.0e-6_dp) error stop "QP solution"
   if (abs(sol%obj_val - 1.88_dp) > 1.0e-6_dp) error stop "QP objective"
   print *, "test_backend_qp: PASS"
end program test_backend_qp
