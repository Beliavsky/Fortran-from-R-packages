program test_quantile
use fields
implicit none
integer,parameter::n=15
real(dp)::x(n),y(n),xx(n,2),pred(3),r
integer::i
type(quantile_spline_fit)::q
type(quantile_tps_fit)::qt
do i=1,n
 x(i)=real(i-1,dp)/real(n-1,dp)
 y(i)=sin(4.0_dp*x(i)); if(i==8)y(i)=y(i)+3.0_dp
 xx(i,:)=[x(i),mod(real(3*i,dp),real(n,dp))/real(n,dp)]
end do
call check(abs(qsreg_psi(2.0_dp,alpha=0.5_dp,c=1.0_dp)-1.0_dp)<1e-14_dp,'psi cap')
call check(abs(qsreg_sigma(2.0_dp,alpha=0.5_dp,c=1.0_dp)-1.5_dp)<1e-14_dp,'sigma tail')
q=quantile_smoothing_spline(x,y,lambda=1.0e-3_dp,maxiter=20)
pred=quantile_spline_predict(q,[0.2_dp,0.5_dp,0.8_dp])
call check(all(abs(pred)<10.0_dp),'quantile spline prediction')
call check(q%iterations>=1 .and. q%fit%ierr==0,'quantile spline fit')
qt=quantile_thin_plate_spline(xx,y,lambda=0.1_dp,maxiter=10)
call check(qt%fit%info==0,'quantile tps fit')
r=sum(abs(quantile_tps_predict(qt,xx)-qt%fit%fitted))
call check(r<1e-8_dp,'quantile tps training prediction')
print *,'test_quantile: PASS'
contains
subroutine check(ok,msg)
logical,intent(in)::ok;character(*),intent(in)::msg
if(.not.ok)then;print *,trim(msg);error stop;end if
end subroutine
end program
