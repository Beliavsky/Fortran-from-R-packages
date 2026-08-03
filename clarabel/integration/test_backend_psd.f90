program test_backend_psd
   use clarabel
   implicit none
   real(dp) :: pd(6,6), ad(6,6), q(6), b(6), target(6)
   type(csc_matrix) :: p, a
   type(clarabel_cone) :: cones(1)
   type(clarabel_settings) :: settings
   type(clarabel_solution) :: sol
   integer :: i, code
   character(len=:), allocatable :: message

   pd = 0.0_dp; ad = 0.0_dp
   do i = 1, 6
      pd(i,i) = 1.0_dp
      ad(i,i) = 1.0_dp
   end do
   q = 0.0_dp
   b = [-3.0_dp, 1.0_dp, 4.0_dp, 1.0_dp, 2.0_dp, 5.0_dp]
   target = [-3.0729833267361095_dp, 0.3696004167288786_dp, -0.022226685581313674_dp, &
              0.31441213129613066_dp, -0.026739700851545107_dp, -0.016084530571308823_dp]
   p = csc_from_symmetric_upper(pd)
   a = csc_from_dense(ad)
   cones(1) = psd_triangle_cone(3)
   settings = default_clarabel_settings(); settings%verbose = .false.
   call clarabel_solve_problem(p, q, a, b, cones, sol, settings, code, message)
   if (code /= 0) error stop message
   if (.not. sol%solved()) error stop "PSD status"
   if (maxval(abs(sol%x - target)) > 1.0e-6_dp) error stop "PSD solution"
   if (abs(sol%obj_val - 4.840076866013861_dp) > 1.0e-6_dp) error stop "PSD objective"
   print *, "test_backend_psd: PASS"
end program test_backend_psd
