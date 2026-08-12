program rastrigin_example
   use psoptim, only : dp, ps_control, ps_result, ps_optimize
   implicit none

   type(ps_control) :: control
   type(ps_result) :: result
   real(dp) :: xmin(2), xmax(2), vmax(2)

   xmin = [-5.12_dp, -5.12_dp]
   xmax = [ 5.12_dp,  5.12_dp]
   vmax = [4.0_dp, 4.0_dp]

   control%n = 80
   control%max_loop = 250
   control%w = 0.8_dp
   control%c1 = 0.5_dp
   control%c2 = 0.5_dp
   control%seed = 5

   call ps_optimize(rastrigin_fitness, xmin, xmax, vmax, result, control)

   print '(a,2f14.8)', 'solution: ', result%sol
   print '(a,es14.6)', 'fitness:  ', result%val

contains

   function rastrigin_fitness(x) result(f)
      real(dp), intent(in) :: x(:)
      real(dp) :: f
      real(dp), parameter :: pi = acos(-1.0_dp)

      f = -(20.0_dp + sum(x*x) - 10.0_dp*sum(cos(2.0_dp*pi*x)))
   end function rastrigin_fitness

end program rastrigin_example
