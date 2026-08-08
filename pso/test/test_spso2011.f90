program test_spso2011
   use pso, only : dp, pso_control, pso_result, psoptim, seed_random, &
      pso_spso2011, rastrigin
   use, intrinsic :: ieee_arithmetic, only : ieee_quiet_nan, ieee_value
   implicit none
   real(dp) :: par(2), lower(2), upper(2)
   type(pso_control) :: con
   type(pso_result) :: res

   call seed_random(24680)
   par = ieee_value(0.0_dp, ieee_quiet_nan)
   lower = -5.12_dp
   upper = 5.12_dp
   con%pso_type = pso_spso2011
   con%swarm_size = 40
   con%maxit = 1000
   con%abstol = 1.0e-10_dp
   call psoptim(par, rastrigin, lower, upper, res, con)
   if (res%value > 1.0e-5_dp) error stop "SPSO2011 did not minimize Rastrigin"
   print *, "test_spso2011: PASS", res%value
end program test_spso2011
