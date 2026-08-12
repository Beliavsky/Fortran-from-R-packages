program test_binary_ga
  use ga
  implicit none
  type(ga_control_type)::c
  type(ga_int_result)::r
  c%pop_size=60;c%max_iter=150;c%run=60;c%seed=4567;c%elitism=3
  c%pcrossover=0.9_dp;c%pmutation=0.2_dp
  call ga_binary(onemax,30,c,r)
  if(nint(r%fitness_value)/=30) error stop "binary GA did not solve OneMax"
  if(any(r%solution/=1)) error stop "binary solution"
  print *, "test_binary_ga: PASS",r%fitness_value
contains
  function onemax(x) result(f)
    integer,intent(in)::x(:);real(dp)::f
    f=real(sum(x),dp)
  end function onemax
end program test_binary_ga
