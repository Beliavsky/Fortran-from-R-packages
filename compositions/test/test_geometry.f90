program test_geometry
  use compositions
  implicit none
  real(dp) :: x(3),y(3),z(2),c(3),v3(3,2),sclr(3,3),silr(2,2),t(3,3),sback(3,3)
  real(dp), allocatable :: pb(:,:)
  real(dp) :: dat(6,4)
  integer :: i
  x=[2.0_dp,3.0_dp,5.0_dp]
  y=closure(x)
  call assert_close(sum(y),1.0_dp,1e-12_dp,'closure')
  c=clr(y); call assert_close(sum(c),0.0_dp,1e-12_dp,'clr zero sum')
  v3=ilr_base(3); z=ilr(y,v3); call assert_vec(ilr_inv(z,v3),y,2e-12_dp,'ilr round trip')
  z=alr(y); call assert_vec(alr_inv(z),y,2e-12_dp,'alr round trip')
  z=apt(y); call assert_vec(apt_inv(z),y,2e-12_dp,'apt round trip')
  z=ipt(y,v3); call assert_vec(ipt_inv(z,v3),y,2e-12_dp,'ipt round trip')
  call assert_vec(ilt_inv(ilt(y)),y,2e-12_dp,'ilt round trip')
  silr=reshape([0.4_dp,0.1_dp,0.1_dp,0.25_dp],[2,2])
  sclr=ilrvar_to_clr(silr,v3); t=clrvar_to_variation(sclr); sback=variation_to_clrvar(t)
  call assert_mat(sback,sclr,1e-10_dp,'variation covariance round trip')
  dat=reshape([ &
    0.55_dp,0.10_dp,0.15_dp,0.20_dp, &
    0.50_dp,0.12_dp,0.18_dp,0.20_dp, &
    0.60_dp,0.08_dp,0.12_dp,0.20_dp, &
    0.15_dp,0.55_dp,0.10_dp,0.20_dp, &
    0.12_dp,0.58_dp,0.10_dp,0.20_dp, &
    0.18_dp,0.52_dp,0.10_dp,0.20_dp],[6,4],order=[2,1])
  pb=principal_balance_maxvar(dat)
  if(any(shape(pb)/=[4,3])) error stop 'principal balance shape'
  call assert_mat(matmul(transpose(pb),pb),identity(3),1e-9_dp,'principal balance orthonormal')
  print *, 'test_geometry: PASS'
contains
  function identity(n) result(a)
    integer,intent(in)::n; real(dp)::a(n,n); integer::j
    a=0; do j=1,n; a(j,j)=1; end do
  end function
  subroutine assert_close(a,b,tol,msg)
    real(dp),intent(in)::a,b,tol; character(*),intent(in)::msg
    if(abs(a-b)>tol) then; print *,trim(msg),a,b; error stop 1; end if
  end subroutine
  subroutine assert_vec(a,b,tol,msg)
    real(dp),intent(in)::a(:),b(:),tol; character(*),intent(in)::msg
    if(maxval(abs(a-b))>tol) then; print *,trim(msg),maxval(abs(a-b)); error stop 1; end if
  end subroutine
  subroutine assert_mat(a,b,tol,msg)
    real(dp),intent(in)::a(:,:),b(:,:),tol; character(*),intent(in)::msg
    if(maxval(abs(a-b))>tol) then; print *,trim(msg),maxval(abs(a-b)); error stop 1; end if
  end subroutine
end program
