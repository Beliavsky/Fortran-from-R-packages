program rosenbrock_example
  use ao
  implicit none
  type(ao_result) :: result
  real(dp) :: initial(2)
  initial = [2.0_dp, 2.0_dp]
  call ao_optimize(rosenbrock, initial, result)
  print '(a,2f14.8)', 'estimate: ', result%estimate
  print '(a,es14.6)', 'value:    ', result%value
  print '(a,a)', 'stop:     ', trim(result%stopping_reason)
contains
  function rosenbrock(x) result(value)
    real(dp), intent(in) :: x(:)
    real(dp) :: value
    value = (1.0_dp-x(1))**2 + (x(2)-x(1)**2)**2
  end function rosenbrock
end program rosenbrock_example
