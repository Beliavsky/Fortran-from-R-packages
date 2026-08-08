program test_gradient_random
  use ao
  implicit none
  type(ao_result) :: r
  type(ao_options) :: o
  real(dp) :: x0(2)
  x0=[0.0_dp,0.0_dp]
  o%partition=AO_PARTITION_RANDOM
  o%minimum_block_number=1
  o%new_block_probability=0.5_dp
  o%iteration_limit=80
  call ao_seed(12345)
  call ao_optimize(himmelblau,x0,r,o,gradient=himmel_grad)
  if(r%value>1.0e-6_dp) error stop 'random gradient AO failed'
  print *, 'PASS test_gradient_random',r%value
contains
  function himmelblau(x) result(f)
    real(dp),intent(in)::x(:); real(dp)::f
    f=(x(1)**2+x(2)-11.0_dp)**2+(x(1)+x(2)**2-7.0_dp)**2
  end function
  subroutine himmel_grad(x,g)
    real(dp),intent(in)::x(:); real(dp),intent(out)::g(:)
    g(1)=4.0_dp*x(1)*(x(1)**2+x(2)-11.0_dp)+2.0_dp*(x(1)+x(2)**2-7.0_dp)
    g(2)=2.0_dp*(x(1)**2+x(2)-11.0_dp)+4.0_dp*x(2)*(x(1)+x(2)**2-7.0_dp)
  end subroutine
end program
