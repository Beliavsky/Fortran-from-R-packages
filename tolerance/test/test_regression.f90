program test_regression
  use tolerance
  implicit none
  real(dp)::x(6,2),y(6),xp(3,2),beta(2),mse,cov(2,2)
  type(regression_band)::band
  integer::i,df,info,fail
  fail=0
  do i=1,6;x(i,:)=[1._dp,real(i,dp)];y(i)=2._dp+3._dp*i+0.1_dp*(-1._dp)**i;end do
  xp(1,:)=[1._dp,1.5_dp];xp(2,:)=[1._dp,3.5_dp];xp(3,:)=[1._dp,6.5_dp]
  call linear_regression_fit(x,y,beta,mse,df,cov,info)
  if(info/=0 .or. abs(beta(2)-3._dp)>0.1_dp)call bad('ols')
  band=regtol_int(x,y,xp,0.1_dp,0.9_dp,1)
  if(any(band%upper<=band%lower))call bad('reg tol')
  band=npregtol_int(y,matmul(x,beta),0.1_dp,0.8_dp,1,'WILKS')
  if(any(band%upper<band%lower))call bad('np reg')
  if(fail==0)then;print '(a)','test_regression: PASS';else;error stop 1;end if
contains
  subroutine bad(nm);character(len=*),intent(in)::nm;print *,trim(nm);fail=fail+1;end subroutine
end program
