program test_hybrid
   use pso, only : dp, pso_control, pso_result, psoptim, seed_random, &
      pso_hybrid_improved, pso_hybrid_on, parabola, parabola_grad
   use, intrinsic :: ieee_arithmetic, only : ieee_quiet_nan, ieee_value
   implicit none
   real(dp) :: par(5), lower(5), upper(5)
   type(pso_control) :: con
   type(pso_result) :: res

   call seed_random(777)
   par = ieee_value(0.0_dp, ieee_quiet_nan)
   lower = -100.0_dp
   upper = 100.0_dp
   con%swarm_size = 15
   con%maxit = 100
   con%hybrid = pso_hybrid_improved
   con%hybrid_maxit = 30
   call psoptim(par, parabola, lower, upper, res, con, parabola_grad)
   if (res%value > 1.0e-12_dp) error stop "hybrid local refinement failed"

   call seed_random(778)
   con = pso_control()
   con%swarm_size = 8
   con%maxit = 5
   con%hybrid = pso_hybrid_on
   con%hybrid_maxit = 20
   call psoptim(par, parabola, lower, upper, res, con)
   if (res%value > 1.0e-8_dp) error stop "numerical-gradient hybrid failed"
   print *, "test_hybrid: PASS", res%value
end program test_hybrid
