module test_controls_functions
   use pso, only : dp
   implicit none
contains
   function concave(x) result(f)
      real(dp), intent(in) :: x(:)
      real(dp) :: f
      f = -(x(1) - 2.0_dp)**2
   end function concave

   function constant_one(x) result(f)
      real(dp), intent(in) :: x(:)
      real(dp) :: f
      f = 1.0_dp + 0.0_dp * sum(x)
   end function constant_one
end module test_controls_functions

program test_controls
   use pso, only : dp, pso_control, pso_result, psoptim, seed_random
   use test_controls_functions, only : concave, constant_one
   use, intrinsic :: ieee_arithmetic, only : ieee_quiet_nan, ieee_value
   implicit none
   real(dp) :: par(1), lower(1), upper(1)
   type(pso_control) :: con
   type(pso_result) :: res

   call seed_random(90210)
   par = ieee_value(0.0_dp, ieee_quiet_nan)
   lower = -5.0_dp
   upper = 5.0_dp
   con%fnscale = -1.0_dp
   con%swarm_size = 20
   con%maxit = 300
   call psoptim(par, concave, lower, upper, res, con)
   if (abs(res%par(1) - 2.0_dp) > 1.0e-3_dp) error stop "fnscale maximization failed"

   call seed_random(111)
   con = pso_control()
   con%swarm_size = 10
   con%maxit = 20
   con%reltol = 2.0_dp
   con%max_restart = 1
   call psoptim(par, constant_one, lower, upper, res, con)
   if (res%convergence /= 3) error stop "restart termination code mismatch"
   if (res%restarts /= 1) error stop "restart count mismatch"
   print *, "test_controls: PASS"
end program test_controls
