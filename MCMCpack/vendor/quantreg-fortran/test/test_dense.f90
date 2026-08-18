program test_dense
  use quantreg, only : dp, rq_result, rq_multi_result, rq_fit_fnb, rq_fit_fnc, rq_fit_qfnb
  implicit none
  real(dp) :: x(5,2), y(5), rmat(1,2), rvec(1)
  type(rq_result) :: fit, cfit
  type(rq_multi_result) :: mf

  x(:,1) = 1.0_dp
  x(:,2) = [0.0_dp,1.0_dp,2.0_dp,3.0_dp,4.0_dp]
  y = [1.0_dp,2.0_dp,4.0_dp,5.0_dp,7.0_dp]
  call rq_fit_fnb(x,y,0.5_dp,fit)
  call assert_true(fit%info == 0, 'fnb info')
  call assert_close(fit%coefficients(1),1.0_dp,1.0e-6_dp,'fnb intercept')
  call assert_close(fit%coefficients(2),1.5_dp,1.0e-6_dp,'fnb slope')

  rmat = 0.0_dp
  rmat(1,2) = 1.0_dp
  rvec(1) = 2.0_dp
  call rq_fit_fnc(x,y,rmat,rvec,0.5_dp,cfit)
  call assert_true(cfit%info == 0, 'fnc info')
  call assert_close(cfit%coefficients(1),0.0_dp,2.0e-4_dp,'fnc intercept')
  call assert_close(cfit%coefficients(2),2.0_dp,2.0e-4_dp,'fnc slope')

  call rq_fit_qfnb(x,y,[0.25_dp,0.5_dp,0.75_dp],mf)
  call assert_true(mf%info == 0, 'qfnb info')
  call assert_close(mf%coefficients(1,2),1.0_dp,1.0e-6_dp,'qfnb median intercept')
  print *, 'test_dense: PASS'
contains
  subroutine assert_true(ok,msg)
    logical, intent(in) :: ok
    character(*), intent(in) :: msg
    if (.not. ok) error stop msg
  end subroutine
  subroutine assert_close(a,b,tol,msg)
    real(dp), intent(in) :: a,b,tol
    character(*), intent(in) :: msg
    if (abs(a-b) > tol) then
      print *, msg, a, b
      error stop 1
    end if
  end subroutine
end program
