program test_spso2007
   use pso, only : dp, pso_control, pso_result, psoptim, seed_random, parabola
   use, intrinsic :: ieee_arithmetic, only : ieee_quiet_nan, ieee_value
   implicit none
   real(dp) :: par(4), lower(4), upper(4)
   type(pso_control) :: con
   type(pso_result) :: res

   call seed_random(12345)
   par = ieee_value(0.0_dp, ieee_quiet_nan)
   lower = -5.0_dp
   upper = 5.0_dp
   con%swarm_size = 30
   con%maxit = 500
   con%abstol = 1.0e-10_dp
   call psoptim(par, parabola, lower, upper, res, con)
   if (res%value > 1.0e-6_dp) error stop "SPSO2007 did not minimize parabola"
   if (any(res%par < lower) .or. any(res%par > upper)) error stop "bounds violated"
   print *, "test_spso2007: PASS", res%value
end program test_spso2007
