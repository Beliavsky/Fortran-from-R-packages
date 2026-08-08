program test_scalar_bounds
   use pso, only : dp, pso_control, pso_result, psoptim, seed_random, parabola
   implicit none
   real(dp) :: par(3)
   type(pso_control) :: con
   type(pso_result) :: res

   call seed_random(42)
   par = 0.0_dp
   con%maxit = 5
   con%abstol = 1.0e-14_dp
   call psoptim(par, parabola, -5.0_dp, 5.0_dp, res, con)
   if (res%value > 1.0e-14_dp) error stop "scalar bounds/start point failed"
   if (maxval(abs(res%par)) > 1.0e-14_dp) error stop "initial point was not retained"
   print *, "test_scalar_bounds: PASS"
end program test_scalar_bounds
