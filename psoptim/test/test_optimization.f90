program test_optimization
   use psoptim, only : dp, ps_control, ps_result, ps_optimize
   implicit none
   type(ps_control) :: control
   type(ps_result) :: result
   real(dp) :: xmin(2), xmax(2), vmax(2)

   xmin = [-5.12_dp, -5.12_dp]
   xmax = [ 5.12_dp,  5.12_dp]
   vmax = [4.0_dp, 4.0_dp]
   control%n = 100
   control%max_loop = 300
   control%w = 0.75_dp
   control%c1 = 0.6_dp
   control%c2 = 0.6_dp
   control%seed = 5

   call ps_optimize(fitness, xmin, xmax, vmax, result, control)
   if (result%val < -1.0e-5_dp) error stop "Rastrigin maximum not reached"
   if (maxval(abs(result%sol)) > 2.0e-3_dp) error stop "Rastrigin solution inaccurate"
   if (result%nfe /= int(control%n*(control%max_loop+1), kind(result%nfe))) &
      error stop "unexpected objective evaluation count"
   print *, "test_optimization: PASS"
contains
   function fitness(x) result(f)
      real(dp), intent(in) :: x(:)
      real(dp) :: f
      real(dp), parameter :: pi = acos(-1.0_dp)
      f = -(20.0_dp + sum(x*x) - 10.0_dp*sum(cos(2.0_dp*pi*x)))
   end function fitness
end program test_optimization
