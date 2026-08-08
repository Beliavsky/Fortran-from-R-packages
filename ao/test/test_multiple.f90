program test_multiple
  use ao
  implicit none
  type(ao_multi_result)::r
  type(ao_options)::o
  real(dp)::initials(2,2)
  initials(:,1)=[2.0_dp,2.0_dp]; initials(:,2)=[-1.0_dp,1.0_dp]
  o%iteration_limit=80
  call ao_optimize_multiple(rosenbrock,initials, &
       [AO_PARTITION_SEQUENTIAL,AO_PARTITION_NONE], &
       [AO_BASE_BFGS,AO_BASE_NELDER_MEAD],r,o)
  if(size(r%process)/=8) error stop 'wrong process count'
  if(r%value>1.0e-5_dp) error stop 'multi process failed'
  if(r%best_process<1 .or. r%best_process>8) error stop 'bad best process'
  print *, 'PASS test_multiple',r%best_process,r%value
contains
  function rosenbrock(x) result(f)
    real(dp),intent(in)::x(:); real(dp)::f
    f=(1.0_dp-x(1))**2+(x(2)-x(1)**2)**2
  end function
end program
