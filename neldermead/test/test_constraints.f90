program test_constraints
  use neldermead
  implicit none
  type(nm_result)::r;type(nm_options)::o
  real(dp)::x0(2),lo(2),hi(2),c(2)
  x0=[15.0_dp,4.99_dp];lo=[13.0_dp,0.0_dp];hi=[20.0_dp,10.0_dp]
  o=fminbnd_options(2);o%seed=2;o%max_iter=1000;o%max_fun_evals=5000;o%box_npoints=4
  o%box_tol_f=1e-7_dp;o%box_nb_match=5
  call fminbnd(g6,x0,lo,hi,r,o,cons,2)
  call cons(r%x,c)
  if(any(c < -1e-8_dp)) error stop 'nonlinear constraint violated'
  if(r%f > -6900.0_dp) error stop 'G6 objective not sufficiently optimized'
contains
  function g6(x) result(f)
    real(dp),intent(in)::x(:);real(dp)::f
    f=(x(1)-10.0_dp)**3+(x(2)-20.0_dp)**3
  end function g6
  subroutine cons(x,c)
    real(dp),intent(in)::x(:);real(dp),intent(out)::c(:)
    c(1)=(x(1)-5.0_dp)**2+(x(2)-5.0_dp)**2-100.0_dp
    c(2)=82.81_dp-(x(1)-6.0_dp)**2-(x(2)-5.0_dp)**2
  end subroutine cons
end program test_constraints
