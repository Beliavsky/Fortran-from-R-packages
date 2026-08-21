program test_diagnostics
  use relsurv, only : dp, rs_br_result, rs_zph_result, rs_br, rs_zph
  implicit none
  real(dp) :: sr(2,1), vv(1,1,2), et(2), nr(2), fv(1,1), beta(1)
  type(rs_br_result) :: br
  type(rs_zph_result) :: zph
  sr(:,1)=[1.0_dp,-1.0_dp];vv=1.0_dp;et=[1.0_dp,2.0_dp];nr=[10.0_dp,5.0_dp]
  call rs_br(sr,vv,et,nr,0.0_dp,'max',br,global=.false.)
  call assert_close(br%statistic(1),sqrt(0.5_dp),1.0e-12_dp,'bridge statistic')
  if(br%p_value(1)<=0.0_dp.or.br%p_value(1)>=1.0_dp)error stop 'bridge p value'
  fv(1,1)=0.25_dp;beta=2.0_dp
  call rs_zph(sr,vv,fv,beta,et,zph,'identity','sum')
  call assert_close(zph%y(1,1),2.5_dp,1.0e-12_dp,'zph first')
  call assert_close(zph%y(2,1),1.5_dp,1.0e-12_dp,'zph second')
  print *, 'test_diagnostics: PASS'
contains
  subroutine assert_close(a,b,tol,msg)
    real(dp),intent(in)::a,b,tol;character(len=*),intent(in)::msg
    if(abs(a-b)>tol)then;print *,'FAIL ',msg,a,b;error stop 1;end if
  end subroutine assert_close
end program test_diagnostics
