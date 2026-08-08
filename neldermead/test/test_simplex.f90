program test_simplex
  use neldermead
  implicit none
  type(nm_simplex)::s;integer::ne
  real(dp)::x0(3)
  x0=0.0_dp;ne=0
  call simplex_build_spendley(sphere,x0,1.0_dp,s,ne)
  if(size(s%x,2)/=4) error stop 'Spendley vertex count failed'
  if(ne/=4) error stop 'Spendley eval count failed'
  if(simplex_size(s)<=0.0_dp) error stop 'simplex size failed'
contains
  function sphere(x) result(f)
    real(dp),intent(in)::x(:);real(dp)::f;f=sum(x*x)
  end function sphere
end program test_simplex
