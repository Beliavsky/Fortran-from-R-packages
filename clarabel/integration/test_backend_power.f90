program test_backend_power
   use clarabel
   implicit none
   real(dp) :: ad(8,6), q(6), b(8)
   type(csc_matrix) :: p, a
   type(clarabel_cone) :: cones(3), gp_cones(3)
   type(clarabel_settings) :: settings
   type(clarabel_solution) :: sol
   integer :: i, code
   character(len=:), allocatable :: message

   ad = 0.0_dp
   do i = 1, 6
      ad(i,i) = -1.0_dp
   end do
   ad(7,1) = 1.0_dp; ad(7,2) = 2.0_dp; ad(7,4) = 3.0_dp
   ad(8,5) = 1.0_dp
   q = [0.0_dp, 0.0_dp, -1.0_dp, 0.0_dp, 0.0_dp, -1.0_dp]
   b = [0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 3.0_dp, 1.0_dp]
   p = csc_empty(6,6)
   a = csc_from_dense(ad)
   cones = [power_cone(0.6_dp), power_cone(0.1_dp), zero_cone(2)]
   settings = default_clarabel_settings(); settings%verbose = .false.
   call clarabel_solve_problem(p, q, a, b, cones, sol, settings, code, message)
   if (code /= 0) error stop message
   if (.not. sol%solved()) error stop "power status"
   if (abs(sol%obj_val + 1.8458_dp) > 1.0e-3_dp) error stop "power objective"

   gp_cones(1) = generalized_power_cone([0.6_dp, 0.4_dp], 1)
   gp_cones(2) = generalized_power_cone([0.1_dp, 0.9_dp], 1)
   gp_cones(3) = zero_cone(2)
   call clarabel_solve_problem(p, q, a, b, gp_cones, sol, settings, code, message)
   if (code /= 0) error stop message
   if (.not. sol%solved()) error stop "generalized-power status"
   if (abs(sol%obj_val + 1.8458_dp) > 1.0e-3_dp) error stop "generalized-power objective"
   print *, "test_backend_power: PASS"
end program test_backend_power
