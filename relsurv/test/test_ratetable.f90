program test_ratetable
  use relsurv, only : dp, ratetable_type, make_ratetable, expected_survival, rstrans_times
  implicit none
  type(ratetable_type) :: tab
  integer :: dims(1), factor(1), ncuts(1)
  real(dp) :: cuts(1,1), rates(1), x(3,1), t(3), s(3), tt(3), lambda
  dims=1; factor=1; ncuts=0; cuts=0.0_dp; lambda=0.01_dp; rates=lambda
  tab=make_ratetable(dims,factor,cuts,ncuts,rates)
  x=1.0_dp; t=[0.0_dp,2.0_dp,10.0_dp]
  s=expected_survival(tab,x,t)
  call assert_close(maxval(abs(s-exp(-lambda*t))),0.0_dp,1.0e-12_dp,'expected survival')
  call rstrans_times(tab,x,t,tt)
  call assert_close(maxval(abs(tt-(1.0_dp-exp(-lambda*t)))),0.0_dp,1.0e-12_dp,'transformed time')
  print *, 'test_ratetable: PASS'
contains
  subroutine assert_close(a,b,tol,msg)
    real(dp),intent(in)::a,b,tol; character(len=*),intent(in)::msg
    if(abs(a-b)>tol) then; print *, 'FAIL ',msg,a,b; error stop 1; end if
  end subroutine
end program
