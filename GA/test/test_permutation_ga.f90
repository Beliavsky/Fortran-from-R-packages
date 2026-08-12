program test_permutation_ga
  use ga
  implicit none
  type(ga_control_type)::c
  type(ga_int_result)::r
  c%pop_size=120
  c%max_iter=400
  c%run=120
  c%seed=9876
  c%elitism=6
  c%pcrossover=0.9_dp
  c%pmutation=0.3_dp
  call ga_permutation(target_fit,1,9,c,r)
  if(abs(r%fitness_value-9.0_dp)>1.0e-12_dp) error stop "permutation GA did not hit target"
  if(any(r%solution/=[2,7,1,9,4,6,3,8,5])) error stop "permutation solution"
  print *, "test_permutation_ga: PASS",r%fitness_value
contains
  function target_fit(x) result(f)
    integer,intent(in)::x(:);real(dp)::f
    integer,parameter::target(9)=[2,7,1,9,4,6,3,8,5]
    f=real(count(x==target),dp)
  end function target_fit
end program test_permutation_ga
