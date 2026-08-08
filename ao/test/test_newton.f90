program test_newton
  use ao
  implicit none
  type(ao_result)::r
  type(ao_options)::o
  real(dp)::x0(2)
  x0=[10.0_dp,-7.0_dp]
  o%base_optimizer=AO_BASE_NEWTON; o%partition=AO_PARTITION_NONE
  o%iteration_limit=5; o%base_max_iterations=5
  call ao_optimize(obj,x0,r,o,gradient=grad,hessian=hess)
  if(r%value>1.0e-20_dp) error stop 'Newton failed'
  if(maxval(abs(r%estimate-[2.0_dp,-1.0_dp]))>1.0e-10_dp) error stop 'Newton estimate wrong'
  print *, 'PASS test_newton',r%value
contains
  function obj(x) result(f)
    real(dp),intent(in)::x(:); real(dp)::f
    f=3.0_dp*(x(1)-2.0_dp)**2+0.5_dp*(x(2)+1.0_dp)**2
  end function
  subroutine grad(x,g)
    real(dp),intent(in)::x(:); real(dp),intent(out)::g(:)
    g=[6.0_dp*(x(1)-2.0_dp),x(2)+1.0_dp]
  end subroutine
  subroutine hess(x,h)
    real(dp),intent(in)::x(:); real(dp),intent(out)::h(:,:)
    if(size(x)<1) error stop
    h=0.0_dp; h(1,1)=6.0_dp; h(2,2)=1.0_dp
  end subroutine
end program
