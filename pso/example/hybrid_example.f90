program hybrid_example
   use pso, only : dp, pso_control, pso_result, psoptim, seed_random, &
      pso_hybrid_improved, parabola, parabola_grad
   use, intrinsic :: ieee_arithmetic, only : ieee_quiet_nan, ieee_value
   implicit none
   real(dp) :: par(4), lower(4), upper(4)
   type(pso_control) :: control
   type(pso_result) :: result

   call seed_random(2)
   par = ieee_value(0.0_dp, ieee_quiet_nan)
   lower = -100.0_dp
   upper = 100.0_dp
   control%hybrid = pso_hybrid_improved
   control%abstol = 1.0e-12_dp
   call psoptim(par, parabola, lower, upper, result, control, parabola_grad)
   print '(a,4es14.6)', "par   = ", result%par
   print '(a,es14.6)', "value = ", result%value
end program hybrid_example
