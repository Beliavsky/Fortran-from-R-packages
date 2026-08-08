program test_maximize
  use ao
  implicit none
  type(ao_result)::r
  type(ao_options)::o
  real(dp)::x0(3)
  x0=[-2.0_dp,4.0_dp,0.0_dp]
  o%minimize=.false.; o%partition=AO_PARTITION_NONE; o%iteration_limit=30
  call ao_optimize(obj,x0,r,o)
  if(abs(r%value)>1.0e-8_dp) error stop 'maximization failed'
  if(maxval(abs(r%estimate-1.0_dp))>1.0e-4_dp) error stop 'maximizer wrong'
  print *, 'PASS test_maximize',r%value
contains
  function obj(x) result(f)
    real(dp),intent(in)::x(:); real(dp)::f
    f=-sum((x-1.0_dp)**2)
  end function
end program
