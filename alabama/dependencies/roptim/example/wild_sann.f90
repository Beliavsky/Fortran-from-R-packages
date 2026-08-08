program wild_sann
  use roptim_mod, only : dp, roptim_control_t, roptim_result_t, &
       roptim_minimize, method_sann, method_bfgs
  implicit none

  real(dp) :: x(1)
  type(roptim_control_t) :: anneal_control, local_control
  type(roptim_result_t) :: anneal_result, local_result

  x = [50.0_dp]
  anneal_control%max_iterations = 20000
  anneal_control%temperature = 20.0_dp
  anneal_control%seed = 123
  allocate(anneal_control%parscale(1), source=20.0_dp)

  call roptim_minimize(x, wild_function, anneal_result, method_sann, &
       control=anneal_control)

  local_control%max_iterations = 200
  call roptim_minimize(x, wild_function, local_result, method_bfgs, &
       control=local_control)

  write(*, '(a,f14.7,a,es14.6)') 'x = ', x(1), ', f = ', local_result%value

contains

  function wild_function(x, user_data) result(f)
    real(dp), intent(in) :: x(:)
    class(*), intent(inout), optional :: user_data
    real(dp) :: f
    f = 10.0_dp*sin(0.3_dp*x(1))*sin(1.3_dp*x(1)**2) + &
         0.00001_dp*x(1)**4 + 0.2_dp*x(1) + 80.0_dp
  end function wild_function

end program wild_sann
