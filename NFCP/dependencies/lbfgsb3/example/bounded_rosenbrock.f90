program bounded_rosenbrock
  use lbfgsb3_mod, only : dp, lbfgsb_control_t, lbfgsb_result_t, &
                         lbfgsb_minimize
  implicit none

  real(dp) :: x(2), lower(2), upper(2)
  type(lbfgsb_control_t) :: control
  type(lbfgsb_result_t) :: result

  x = [-1.2_dp, 1.0_dp]
  lower = -2.0_dp
  upper =  2.0_dp
  control%pgtol = 1.0e-9_dp
  control%reltol = 0.0_dp

  call lbfgsb_minimize(x, rosenbrock, result, lower, upper, control)

  write(*,'(a,2(1x,f16.10))') 'x:', x
  write(*,'(a,1x,es16.8)') 'f:', result%value
  write(*,'(a,1x,a)') 'status:', result%message
  write(*,'(a,1x,i0)') 'iterations:', result%iterations

contains

  subroutine rosenbrock(z, f, g, user_data)
    real(dp), intent(in) :: z(:)
    real(dp), intent(out) :: f, g(:)
    class(*), intent(inout), optional :: user_data

    f = 100.0_dp*(z(2)-z(1)**2)**2 + (1.0_dp-z(1))**2
    g(1) = -400.0_dp*z(1)*(z(2)-z(1)**2) - 2.0_dp*(1.0_dp-z(1))
    g(2) =  200.0_dp*(z(2)-z(1)**2)
  end subroutine rosenbrock

end program bounded_rosenbrock
