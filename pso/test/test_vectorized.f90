program test_vectorized
   use pso, only : dp, pso_control, pso_result, psoptim, seed_random, &
      pso_spso2011, parabola
   use, intrinsic :: ieee_arithmetic, only : ieee_quiet_nan, ieee_value
   implicit none
   real(dp) :: par(3), lower(3), upper(3)
   type(pso_control) :: con
   type(pso_result) :: res

   call seed_random(314159)
   par = ieee_value(0.0_dp, ieee_quiet_nan)
   lower = -10.0_dp
   upper = 10.0_dp
   con%pso_type = pso_spso2011
   con%vectorize = .true.
   con%swarm_size = 40
   con%maxit = 600
   con%v_max = 0.5_dp
   call psoptim(par, parabola, lower, upper, res, con)
   if (res%value > 1.0e-6_dp) error stop "vectorized SPSO2011 failed"

   call seed_random(271828)
   con = pso_control()
   con%vectorize = .true.
   con%swarm_size = 30
   con%maxit = 500
   call psoptim(par, parabola, lower, upper, res, con)
   if (res%value > 1.0e-6_dp) error stop "vectorized SPSO2007 failed"
   print *, "test_vectorized: PASS", res%value
end program test_vectorized
