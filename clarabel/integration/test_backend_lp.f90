program test_backend_lp
   use clarabel
   implicit none
   real(dp) :: ad(6,3), q(3), b(6), target(3)
   type(csc_matrix) :: p, a
   type(clarabel_cone) :: cones(2)
   type(clarabel_settings) :: settings
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
   target = [-0.5_dp, 0.5_dp, -0.5_dp]
   p = csc_empty(3,3)
   a = csc_from_dense(ad)
   cones = [nonnegative_cone(3), nonnegative_cone(3)]
   settings = default_clarabel_settings(); settings%verbose = .false.
   call clarabel_solve_problem(p, q, a, b, cones, sol, settings, code, message)
   if (code /= 0) error stop message
   if (.not. sol%solved()) error stop "LP status"
   if (maxval(abs(sol%x - target)) > 1.0e-7_dp) error stop "LP solution"
   if (abs(sol%obj_val + 3.0_dp) > 1.0e-7_dp) error stop "LP objective"
   print *, "test_backend_lp: PASS"
end program test_backend_lp
