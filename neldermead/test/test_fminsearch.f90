program test_fminsearch
  use neldermead
  implicit none
  type(nm_result)::r
  real(dp)::x0(2)
  x0=[-1.2_dp,1.0_dp]
  call fminsearch(rosen,x0,r)
  if(maxval(abs(r%x-[1.0_dp,1.0_dp]))>2.0e-3_dp) error stop 'fminsearch Rosenbrock failed'
  if(r%f>1.0e-6_dp) error stop 'fminsearch objective too large'
contains
  function rosen(x) result(f)
    real(dp),intent(in)::x(:);real(dp)::f
    f=100.0_dp*(x(2)-x(1)**2)**2+(1.0_dp-x(1))**2
  end function rosen
end program test_fminsearch
