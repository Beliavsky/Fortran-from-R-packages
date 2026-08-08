program fminsearch_rosenbrock
  use neldermead
  implicit none
  type(nm_result)::r
  call fminsearch(rosen,[-1.2_dp,1.0_dp],r)
  print '(a,2(1x,es16.8))','x =',r%x
  print '(a,1x,es16.8)','f =',r%f
  print '(a,1x,a)','status =',trim(r%status)
contains
  function rosen(x) result(f)
    real(dp),intent(in)::x(:);real(dp)::f
    f=100.0_dp*(x(2)-x(1)**2)**2+(1.0_dp-x(1))**2
  end function rosen
end program fminsearch_rosenbrock
