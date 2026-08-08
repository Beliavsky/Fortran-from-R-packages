program test_bounds
  use optimflex
  implicit none
  type(optim_result) :: r
  type(optim_control) :: c
  real(dp) :: x0(2), lo(2), hi(2)
  x0=[0.0_dp,0.0_dp]; lo=[-1.0_dp,-3.0_dp]; hi=[1.0_dp,3.0_dp]
  c=lbfgsb_default_control(); c%max_iter=1000; c%use_posdef=.false.
  call l_bfgs_b(x0,obj,r,lo,hi,control=c)
  if(.not.r%converged) error stop 'bounded L-BFGS failed'
  if(abs(r%par(1)-1.0_dp)>1.0e-5_dp .or. abs(r%par(2)+2.0_dp)>1.0e-5_dp) error stop 'bounds mismatch'
  print *, 'test_bounds: PASS'
contains
  real(dp) function obj(x) result(f)
    real(dp), intent(in) :: x(:)
    f=(x(1)-2.0_dp)**2+(x(2)+2.0_dp)**2
  end function obj
end program test_bounds
