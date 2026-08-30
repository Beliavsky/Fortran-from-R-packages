program rosenbrock_example
   use ucminf, only : dp, ucminf_options, ucminf_result, ucminf_minimize
   implicit none

   type(ucminf_options) :: options
   type(ucminf_result) :: result
   real(dp) :: start(2)

   start = [2.0_dp, 0.5_dp]
   options%trace = .true.
   call ucminf_minimize(start, rosenbrock, result, rosenbrock_gradient, options)

   write(*,'(a,2f16.9)') "minimum at: ", result%par
   write(*,'(a,es16.8)') "objective:  ", result%value
   write(*,'(a,i0)') "status:     ", result%convergence

contains

   function rosenbrock(x) result(f)
      real(dp), intent(in) :: x(:)
      real(dp) :: f
      f = (1.0_dp - x(1))**2 + 100.0_dp * (x(2) - x(1)**2)**2
   end function rosenbrock

   subroutine rosenbrock_gradient(x, g)
      real(dp), intent(in) :: x(:)
      real(dp), intent(out) :: g(:)
      g(1) = -400.0_dp * x(1) * (x(2) - x(1)**2) - 2.0_dp * (1.0_dp - x(1))
      g(2) = 200.0_dp * (x(2) - x(1)**2)
   end subroutine rosenbrock_gradient

end program rosenbrock_example
