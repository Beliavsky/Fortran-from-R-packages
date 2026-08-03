! SPDX-License-Identifier: Artistic-2.0
program test_fitting
  use ecd_api
  implicit none
  type(optimization_result) :: r
  type(ecd_model) :: d
  real(dp) :: x0(2),lo(2),hi(2)
  x0=[-3.0_dp,4.0_dp]; lo=[-10.0_dp,-10.0_dp]; hi=[10.0_dp,10.0_dp]
  call nelder_mead(quadratic,x0,r,lo,hi,max_iter=500,tolerance=1.0e-11_dp)
  call check(r%status==ecd_ok,'nelder-mead convergence')
  call close(r%parameters(1),1.25_dp,2.0e-5_dp,'nelder-mead x1')
  call close(r%parameters(2),-0.75_dp,2.0e-5_dp,'nelder-mead x2')
  call check(r%objective<1.0e-10_dp,'nelder-mead objective')
  d=ecd_new(alpha=0.0_dp,gamma=0.0_dp,sigma=1.0_dp,beta=0.0_dp,mu=0.0_dp)
  call close(ecd_estimate_const(d),2.0_dp*sqrt(acos(-1.0_dp)/2.0_dp)*sqrt(63.0_dp/8.0_dp), &
    2.0e-14_dp,'constant estimate formula')
  print '(a)', 'test_fitting: PASS'
contains
  function quadratic(x) result(v)
    real(dp),intent(in)::x(:)
    real(dp)::v
    v=(x(1)-1.25_dp)**2+3.0_dp*(x(2)+0.75_dp)**2
  end function quadratic
  subroutine close(a,b,tol,msg)
    real(dp),intent(in)::a,b,tol
    character(len=*),intent(in)::msg
    if(abs(a-b)>tol*max(1.0_dp,abs(b)))then
      write(*,*)trim(msg),a,b; error stop 1
    end if
  end subroutine close
  subroutine check(ok,msg)
    logical,intent(in)::ok
    character(len=*),intent(in)::msg
    if(.not.ok)then;write(*,*)trim(msg);error stop 1;end if
  end subroutine check
end program test_fitting
