program test_box
  use neldermead
  implicit none
  type(nm_result)::r;type(nm_options)::o
  real(dp)::x0(2),lo(2),hi(2)
  x0=[1.2_dp,1.9_dp];lo=[1.0_dp,1.0_dp];hi=[2.0_dp,2.0_dp]
  o=fminbnd_options(2);o%seed=777;o%max_iter=1200;o%max_fun_evals=5000;o%box_tol_f=1e-8_dp;o%box_nb_match=8
  call fminbnd(sphere,x0,lo,hi,r,o)
  if(maxval(abs(r%x-[1.0_dp,1.0_dp]))>2.0e-2_dp) error stop 'Box bounded optimum failed'
  if(any(r%x<lo).or.any(r%x>hi)) error stop 'Box returned infeasible bound point'
contains
  function sphere(x) result(f)
    real(dp),intent(in)::x(:);real(dp)::f;f=sum(x*x)
  end function sphere
end program test_box
