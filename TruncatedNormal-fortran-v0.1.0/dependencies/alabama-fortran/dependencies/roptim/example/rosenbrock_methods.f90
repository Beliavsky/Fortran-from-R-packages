program rosenbrock_methods
  use roptim_mod, only : dp, roptim_control_t, roptim_result_t, &
       roptim_minimize, method_nelder_mead, method_bfgs, method_cg, method_lbfgsb
  implicit none

  character(len=16), parameter :: methods(4) = [character(len=16) :: &
       method_nelder_mead, method_bfgs, method_cg, method_lbfgsb]
  real(dp) :: x(2), lower(2), upper(2)
  integer :: i
  type(roptim_control_t) :: control
  type(roptim_result_t) :: result

  lower = [-2.0_dp, -1.0_dp]
  upper = [ 2.0_dp,  3.0_dp]

  do i = 1, size(methods)
    x = [-1.2_dp, 1.0_dp]
    control = roptim_control_t()
    control%max_iterations = 1000
    control%reltol = 1.0e-9_dp
    control%compute_hessian = methods(i) == method_bfgs

    if (methods(i) == method_nelder_mead) then
      call roptim_minimize(x, rosenbrock, result, methods(i), control=control)
    else if (methods(i) == method_lbfgsb) then
      call roptim_minimize(x, rosenbrock, result, methods(i), &
           gradient=rosenbrock_gradient, lower=lower, upper=upper, &
           control=control)
    else
      call roptim_minimize(x, rosenbrock, result, methods(i), &
           gradient=rosenbrock_gradient, control=control)
    end if

    write(*, '(a16,2x,a,l1,2x,a,2f12.7,2x,a,es13.5)') trim(methods(i)), &
         'success=', result%success, 'x=', x, 'f=', result%value
  end do

contains

  function rosenbrock(x, user_data) result(f)
    real(dp), intent(in) :: x(:)
    class(*), intent(inout), optional :: user_data
    real(dp) :: f
    f = 100.0_dp*(x(2)-x(1)*x(1))**2 + (1.0_dp-x(1))**2
  end function rosenbrock

  subroutine rosenbrock_gradient(x, g, user_data)
    real(dp), intent(in) :: x(:)
    real(dp), intent(out) :: g(:)
    class(*), intent(inout), optional :: user_data
    g(1) = -400.0_dp*x(1)*(x(2)-x(1)*x(1)) - 2.0_dp*(1.0_dp-x(1))
    g(2) = 200.0_dp*(x(2)-x(1)*x(1))
  end subroutine rosenbrock_gradient

end program rosenbrock_methods
