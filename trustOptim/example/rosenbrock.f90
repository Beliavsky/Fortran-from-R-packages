program rosenbrock_example
   use trustoptim
   implicit none
   real(dp) :: x0(2)
   type(trustoptim_control) :: con
   type(trustoptim_result) :: res

   x0 = [-1.2_dp, 1.0_dp]
   con%maxit = 1000
   con%prec = 1.0e-8_dp
   con%report_level = 0
   call trust_optim_bfgs(x0, rosen, grad, res, con)

   write(*,'(a,2f14.8)') 'solution: ', res%solution
   write(*,'(a,es14.6)') 'fval:     ', res%fval
   write(*,'(a,a)') 'status:   ', trim(res%status_message())
contains
   function rosen(x) result(f)
      real(dp), intent(in) :: x(:)
      real(dp) :: f
      f = 100.0_dp*(x(1)**2-x(2))**2 + (x(1)-1.0_dp)**2
   end function rosen
   subroutine grad(x,g)
      real(dp), intent(in) :: x(:)
      real(dp), intent(out) :: g(:)
      g(1) = 400.0_dp*(x(1)**2-x(2))*x(1) + 2.0_dp*(x(1)-1.0_dp)
      g(2) = -200.0_dp*(x(1)**2-x(2))
   end subroutine grad
end program rosenbrock_example
