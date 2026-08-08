program box_g6
  use neldermead
  implicit none
  type(nm_result)::r;type(nm_options)::o
  o=fminbnd_options(2);o%max_iter=1000;o%max_fun_evals=5000;o%box_npoints=4;o%seed=2
  call fminbnd(g6,[15.0_dp,4.99_dp],[13.0_dp,0.0_dp],[20.0_dp,10.0_dp],r,o,cons,2)
  print '(a,2(1x,es16.8))','x =',r%x
  print '(a,1x,es16.8)','f =',r%f
  print '(a,1x,a)','status =',trim(r%status)
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
end program box_g6
