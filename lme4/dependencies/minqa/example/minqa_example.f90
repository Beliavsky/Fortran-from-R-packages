program minqa_example
   use minqa_module, only : dp, minqa_control_t, minqa_result_t, bobyqa, newuoa, uobyqa
   implicit none

   real(dp) :: x(2)
   type(minqa_control_t) :: control
   type(minqa_result_t) :: result

   control%rhobeg = 0.25_dp
   control%rhoend = 1.0e-7_dp
   control%maxfun = 5000

   x = [-1.2_dp, 1.0_dp]
   call bobyqa(rosenbrock, x, result, [-2.0_dp, -1.0_dp], [2.0_dp, 3.0_dp], control)
   call print_result('BOBYQA', result)

   x = [-1.2_dp, 1.0_dp]
   call newuoa(rosenbrock, x, result, control)
   call print_result('NEWUOA', result)

   x = [-1.2_dp, 1.0_dp]
   call uobyqa(rosenbrock, x, result, control)
   call print_result('UOBYQA', result)

contains

   function rosenbrock(x) result(f)
      real(dp), intent(in) :: x(:)
      real(dp) :: f

      f = 100.0_dp * (x(2) - x(1)**2)**2 + (1.0_dp - x(1))**2
   end function rosenbrock

   subroutine print_result(name, result)
      character(len=*), intent(in) :: name
      type(minqa_result_t), intent(in) :: result

      write(*, '(/,a)') name
      write(*, '("  status: ",i0," (",a,")")') result%status, result%message
      write(*, '("  evaluations: ",i0)') result%evaluations
      write(*, '("  objective: ",es14.6)') result%fval
      write(*, '("  x:",*(1x,f12.8))') result%x
   end subroutine print_result

end program minqa_example
