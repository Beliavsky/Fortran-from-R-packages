program test_least_squares
  use optimflex
  implicit none
  type(optim_result) :: r
  type(optim_control) :: c
  real(dp) :: x0(2), lo(2), hi(2)
  x0 = [0.0_dp,0.0_dp]; lo=-10.0_dp; hi=10.0_dp
  c=gauss_newton_default_control(); c%max_iter=100
  call gauss_newton(x0,obj,r,residual=resid,jacobian=jac,control=c); call check(r,'gauss_newton')
  c=lm_default_control(); c%max_iter=100
  call levenberg_marquardt(x0,obj,r,lo,hi,residual=resid,jacobian=jac,control=c); call check(r,'lm_ls')
  c=dogleg_default_control(); c%max_iter=100
  call dogleg(x0,obj,r,lo,hi,residual=resid,jacobian=jac,control=c); call check(r,'dogleg_ls')
  c=double_dogleg_default_control(); c%max_iter=100
  call double_dogleg(x0,obj,r,lo,hi,residual=resid,jacobian=jac,control=c); call check(r,'double_dogleg_ls')
  print *, 'test_least_squares: PASS'
contains
  real(dp) function obj(x) result(f)
    real(dp), intent(in) :: x(:)
    real(dp), allocatable :: rr(:)
    rr=resid(x); f=sum(rr*rr)
  end function obj
  function resid(x) result(rr)
    real(dp), intent(in) :: x(:)
    real(dp), allocatable :: rr(:)
    allocate(rr(3)); rr=[x(1)-2.0_dp, x(2)+1.0_dp, x(1)+x(2)-1.0_dp]
  end function resid
  function jac(x) result(j)
    real(dp), intent(in) :: x(:)
    real(dp), allocatable :: j(:,:)
    allocate(j(3,2)); j=0.0_dp*x(1)
    j(1,1)=1.0_dp; j(2,2)=1.0_dp; j(3,:)=[1.0_dp,1.0_dp]
  end function jac
  subroutine check(res,name)
    type(optim_result), intent(in) :: res
    character(len=*), intent(in) :: name
    if (.not.res%converged) then
      print *, trim(name), trim(res%status), res%par, res%objective
      error stop 'least-squares convergence failed'
    end if
    if(maxval(abs(res%par-[2.0_dp,-1.0_dp]))>1.0e-4_dp) error stop 'least-squares mismatch'
  end subroutine check
end program test_least_squares
