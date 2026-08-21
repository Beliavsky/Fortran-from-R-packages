program test_epa
  use relsurv, only : dp, epanechnikov_boundary_matrix, epa_smooth
  implicit none
  real(dp) :: source(4),target(1),bw(1),k(1,4),lam(4),v(1)
  integer :: cuts(2)
  source=[0.0_dp,1.0_dp,2.0_dp,3.0_dp];target=1.5_dp;bw=1.0_dp;cuts=[1,4]
  call epanechnikov_boundary_matrix(target,source,bw,cuts,k)
  call assert_close(k(1,2),0.5625_dp,1.0e-12_dp,'interior left')
  call assert_close(k(1,3),0.5625_dp,1.0e-12_dp,'interior right')
  lam=[0.0_dp,1.0_dp,1.0_dp,0.0_dp]
  call epa_smooth(source,lam,target,v,bwin=25.0_dp,n_bwin=1)
  if(v(1)<0.0_dp)error stop 'epa nonnegative'
  print *, 'test_epa: PASS'
contains
  subroutine assert_close(a,b,tol,msg)
    real(dp),intent(in)::a,b,tol;character(len=*),intent(in)::msg
    if(abs(a-b)>tol)then;print *,'FAIL ',msg,a,b;error stop 1;end if
  end subroutine assert_close
end program test_epa
