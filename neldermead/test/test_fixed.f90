program test_fixed
  use neldermead
  implicit none
  type(nm_options)::o;type(nm_result)::r;real(dp)::x0(2)
  x0=[2.0_dp,-3.0_dp];o%method='fixed';o%simplex0_method='spendley';o%simplex0_length=1.0_dp
  o%max_iter=1500;o%max_fun_evals=5000;o%tol_x_method=.false.;o%tol_fun_method=.false.
  o%tol_simplex_size_method=.true.;o%tol_simplex_size_absolute=1e-7_dp;o%tol_simplex_size_relative=0.0_dp
  call neldermead_search(sphere,x0,o,r)
  if(sqrt(sum(r%x*r%x))>2.0e-3_dp) error stop 'fixed simplex failed'
contains
  function sphere(x) result(f)
    real(dp),intent(in)::x(:);real(dp)::f;f=sum(x*x)
  end function sphere
end program test_fixed
