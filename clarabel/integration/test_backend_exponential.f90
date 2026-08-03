program test_backend_exponential
   use clarabel
   implicit none
   real(dp) :: ad(5,3), q(3), b(5), target(3)
   type(csc_matrix) :: p, a
   type(clarabel_cone) :: cones(2)
   type(clarabel_settings) :: settings
   type(clarabel_solution) :: sol
   integer :: i, code
   character(len=:), allocatable :: message

   ad = 0.0_dp
   do i = 1, 3
      ad(i,i) = -1.0_dp
   end do
   ad(4,2) = 1.0_dp
   ad(5,3) = 1.0_dp
   q = [-1.0_dp, 0.0_dp, 0.0_dp]
   b = [0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp, exp(5.0_dp)]
   target = [5.0_dp, 1.0_dp, exp(5.0_dp)]
   p = csc_empty(3,3)
   a = csc_from_dense(ad)
   cones = [exponential_cone(), zero_cone(2)]
   settings = default_clarabel_settings(); settings%verbose = .false.
   call clarabel_solve_problem(p, q, a, b, cones, sol, settings, code, message)
   if (code /= 0) error stop message
   if (.not. sol%solved()) error stop "exponential status"
   if (maxval(abs(sol%x - target)) > 1.0e-5_dp) error stop "exponential solution"
   if (abs(sol%obj_val + 5.0_dp) > 1.0e-6_dp) error stop "exponential objective"
   print *, "test_backend_exponential: PASS"
end program test_backend_exponential
