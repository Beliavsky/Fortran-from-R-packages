program test_ao_basic
  use ao
  implicit none
  type(ao_result) :: r
  type(ao_options) :: o
  real(dp) :: x0(2)
  x0=[2.0_dp,2.0_dp]
  o%iteration_limit=100
  o%base_max_iterations=10
  call ao_optimize(rosenbrock,x0,r,o)
  if (r%value > 1.0e-5_dp) error stop 'Rosenbrock did not converge'
  if (maxval(abs(r%estimate-[1.0_dp,1.0_dp])) > 5.0e-3_dp) error stop 'bad estimate'
  print *, 'PASS test_ao_basic', r%estimate, r%value, trim(r%stopping_reason)
contains
  function rosenbrock(x) result(f)
    real(dp),intent(in)::x(:)
    real(dp)::f
    f=(1.0_dp-x(1))**2+(x(2)-x(1)**2)**2
  end function
end program
