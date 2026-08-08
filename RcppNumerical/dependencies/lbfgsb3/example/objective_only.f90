program objective_only
  use lbfgsb3_mod, only : dp, lbfgsb_control_t, lbfgsb_result_t, &
                         lbfgsb_minimize_fd
  implicit none

  real(dp) :: x(2)
  type(lbfgsb_control_t) :: control
  type(lbfgsb_result_t) :: result

  x = [-1.2_dp, 1.0_dp]
  control%pgtol = 1.0e-7_dp
  control%reltol = 0.0_dp
  control%max_evaluations = 20000

  call lbfgsb_minimize_fd(x, rosenbrock_value, result, &
                          lower=[-2.0_dp], upper=[2.0_dp], control=control)

  write(*,'(a,2(1x,f16.10))') 'x:', x
  write(*,'(a,1x,es16.8)') 'f:', result%value
  write(*,'(a,1x,i0)') 'objective calls:', result%function_evaluations

contains

  function rosenbrock_value(z, user_data) result(f)
    real(dp), intent(in) :: z(:)
    class(*), intent(inout), optional :: user_data
    real(dp) :: f

    f = 100.0_dp*(z(2)-z(1)**2)**2 + (1.0_dp-z(1))**2
  end function rosenbrock_value

end program objective_only
