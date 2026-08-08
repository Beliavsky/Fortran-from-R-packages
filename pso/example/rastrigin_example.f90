program rastrigin_example
   use pso, only : dp, pso_control, pso_result, psoptim, seed_random, rastrigin
   use, intrinsic :: ieee_arithmetic, only : ieee_quiet_nan, ieee_value
   implicit none
   real(dp) :: par(2), lower(2), upper(2)
   type(pso_control) :: control
   type(pso_result) :: result

   call seed_random(1)
   par = ieee_value(0.0_dp, ieee_quiet_nan)
   lower = -5.0_dp
   upper = 5.0_dp
   control%abstol = 1.0e-8_dp
   call psoptim(par, rastrigin, lower, upper, result, control)
   print '(a,2f14.8)', "par   = ", result%par
   print '(a,es14.6)', "value = ", result%value
   print '(a,i0)', "feval = ", result%function_evaluations
end program rastrigin_example
