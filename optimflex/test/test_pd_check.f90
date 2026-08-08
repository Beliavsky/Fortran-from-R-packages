program test_pd_check
  use optimflex
  implicit none
  type(optim_result) :: r
  type(optim_control) :: c
  real(dp) :: x0(2)
  x0=[0.0_dp,0.0_dp]
  c=newton_default_control(); c%max_iter=2; c%use_posdef=.true.
  call newton_raphson(x0,saddle,r,gradient=grad,hessian=hess,control=c)
  if(r%converged) error stop 'saddle incorrectly accepted as minimum'
  if(r%hess_is_pd) error stop 'saddle Hessian marked PD'
  print *, 'test_pd_check: PASS'
contains
  real(dp) function saddle(x) result(f)
    real(dp), intent(in) :: x(:)
    f=x(1)*x(1)-x(2)*x(2)
  end function saddle
  subroutine grad(x,g)
    real(dp), intent(in) :: x(:)
    real(dp), intent(out) :: g(:)
    g=[2.0_dp*x(1),-2.0_dp*x(2)]
  end subroutine grad
  subroutine hess(x,h)
    real(dp), intent(in) :: x(:)
    real(dp), intent(out) :: h(:,:)
    h=0.0_dp*x(1); h(1,1)=2.0_dp; h(2,2)=-2.0_dp
  end subroutine hess
end program test_pd_check
