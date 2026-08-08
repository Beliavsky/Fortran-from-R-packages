program test_gn_hessian
  use optimflex
  implicit none
  type(optim_result) :: r
  type(optim_control) :: c
  real(dp) :: x0(2), lo(2), hi(2)
  x0=[4.0_dp,-4.0_dp]; lo=-10.0_dp; hi=10.0_dp
  c=lm_default_control(); c%max_iter=200
  call levenberg_marquardt(x0,obj,r,lo,hi,gradient=grad,hessian=hess,gn_hessian=gn,control=c)
  call check(r,'lm gn_hessian')
  c=double_dogleg_default_control(); c%max_iter=200
  call double_dogleg(x0,obj,r,lo,hi,gradient=grad,hessian=hess,gn_hessian=gn,control=c)
  call check(r,'double dogleg gn_hessian')
  print *, 'test_gn_hessian: PASS'
contains
  real(dp) function obj(x) result(f)
    real(dp),intent(in)::x(:)
    f=(x(1)-2.0_dp)**2+2.0_dp*(x(2)+1.0_dp)**2
  end function obj
  subroutine grad(x,g)
    real(dp),intent(in)::x(:); real(dp),intent(out)::g(:)
    g=[2.0_dp*(x(1)-2.0_dp),4.0_dp*(x(2)+1.0_dp)]
  end subroutine grad
  subroutine hess(x,h)
    real(dp),intent(in)::x(:); real(dp),intent(out)::h(:,:)
    h=0.0_dp*x(1); h(1,1)=2.0_dp; h(2,2)=4.0_dp
  end subroutine hess
  subroutine gn(x,h)
    real(dp),intent(in)::x(:); real(dp),intent(out)::h(:,:)
    h=0.0_dp*x(1); h(1,1)=2.0_dp; h(2,2)=4.0_dp
  end subroutine gn
  subroutine check(res,name)
    type(optim_result),intent(in)::res; character(len=*),intent(in)::name
    if(.not.res%converged) then
      print *,trim(name),trim(res%status),res%par
      error stop 'gn_hessian convergence failed'
    end if
    if(maxval(abs(res%par-[2.0_dp,-1.0_dp]))>1.0e-5_dp) error stop 'gn_hessian mismatch'
  end subroutine check
end program test_gn_hessian
