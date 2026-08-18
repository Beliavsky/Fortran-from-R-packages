program test_nonlinear
  use tolerance
  implicit none
  real(dp)::x(8,1),y(8),xp(3,1),b(2)
  type(regression_band)::band
  integer::i,info,fail
  fail=0
  do i=1,8;x(i,1)=0.5_dp*i;y(i)=1.5_dp+2.2_dp*x(i,1)+0.02_dp*(-1._dp)**i;end do
  xp(:,1)=[0.75_dp,2.25_dp,4.25_dp];b=[1._dp,2._dp]
  call nlregtol_int(model,x,y,b,xp,band,0.1_dp,0.9_dp,1,1000,info)
  if(info/=0 .or. abs(b(1)-1.5_dp)>0.1_dp .or. abs(b(2)-2.2_dp)>0.1_dp)call bad('nl fit')
  if(any(band%upper<=band%lower))call bad('nl band')
  if(fail==0)then;print '(a)','test_nonlinear: PASS';else;error stop 1;end if
contains
  subroutine model(xx,beta,yy)
    real(dp),intent(in)::xx(:,:),beta(:);real(dp),intent(out)::yy(:)
    yy=beta(1)+beta(2)*xx(:,1)
  end subroutine
  subroutine bad(nm);character(len=*),intent(in)::nm;print *,trim(nm);fail=fail+1;end subroutine
end program
