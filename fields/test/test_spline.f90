program test_spline
use fields, only: dp,spline_fit,smoothing_spline,spline_predict,spline_df_to_lambda
implicit none
integer,parameter::n=12
real(dp)::x(n),y(n),xn(3),lam
real(dp),allocatable::yp(:),dy(:)
type(spline_fit)::fit
integer::i
do i=1,n;x(i)=real(i-1,dp)/real(n-1,dp);y(i)=2.0_dp+3.0_dp*x(i);end do
fit=smoothing_spline(x,y,1.0_dp)
call check(fit%ierr==0,'css ierr')
call check(maxval(abs(fit%fitted-y))<2e-9_dp,'linear null space')
xn=[0.1_dp,0.45_dp,0.9_dp]
yp=spline_predict(fit,xn);dy=spline_predict(fit,xn,1)
call check(maxval(abs(yp-(2.0_dp+3.0_dp*xn)))<2e-8_dp,'spline prediction')
call check(maxval(abs(dy-3.0_dp))<2e-7_dp,'spline derivative')
lam=spline_df_to_lambda(6.0_dp,x)
call check(lam>0.0_dp .and. lam<huge(1.0_dp),'df to lambda')
print *,'test_spline: PASS'
contains
subroutine check(ok,msg)
logical,intent(in)::ok;character(len=*),intent(in)::msg
if(.not.ok)then;print *,'FAIL: ',msg;error stop 1;end if
end subroutine
end program
