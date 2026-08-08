program test_grid
  use neldermead
  implicit none
  type(grid_result)::g;real(dp)::x0(2),lo(2),hi(2)
  x0=[0.0_dp,0.0_dp];lo=[-1.0_dp,-1.0_dp];hi=[1.0_dp,1.0_dp]
  call fmin_gridsearch(sphere,x0,5,[10.0_dp],g,lo,hi)
  if(size(g%f)/=25) error stop 'grid point count failed'
  if(abs(g%f(1))>epsilon(1.0_dp)) error stop 'grid sorting/minimum failed'
contains
  function sphere(x) result(f)
    real(dp),intent(in)::x(:);real(dp)::f;f=sum(x*x)
  end function sphere
end program test_grid
