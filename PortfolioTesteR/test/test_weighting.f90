! SPDX-License-Identifier: MIT
! PortfolioTesteR modern Fortran translation
program test_weighting
  use portfolio_tester
  implicit none
  real(dp)::cov(3,3),marg(3),rc(3),pv
  real(dp),allocatable::w(:),wh(:),wm(:),ret(:,:)
  integer::st,i
  cov=reshape([0.04_dp,0.006_dp,0.004_dp, &
               0.006_dp,0.09_dp,0.008_dp, &
               0.004_dp,0.008_dp,0.16_dp],[3,3])
  call calculate_erc_weights(cov,w,st)
  call assert_true(st==0,'erc status')
  marg=matmul(cov,w);pv=dot_product(w,marg);rc=w*marg/pv
  call assert_true(maxval(rc)-minval(rc)<1.0e-6_dp,'equal risk contributions')
  call calculate_max_div_weights(cov,wm,st)
  call assert_true(st==0,'max div status')
  call assert_close(sum(wm),1.0_dp,1.0e-12_dp,'max div sum')
  allocate(ret(200,3))
  do i=1,200
    ret(i,1)=0.01_dp*sin(0.13_dp*real(i,dp))
    ret(i,2)=0.015_dp*cos(0.09_dp*real(i,dp))+0.2_dp*ret(i,1)
    ret(i,3)=0.02_dp*sin(0.05_dp*real(i,dp)+1.0_dp)+0.1_dp*ret(i,1)
  end do
  call calculate_hrp_weights(ret,wh,st)
  call assert_true(st==0,'hrp status')
  call assert_close(sum(wh),1.0_dp,1.0e-12_dp,'hrp sum')
  call assert_true(all(wh>=0.0_dp),'hrp nonnegative')
  print '(a)','test_weighting: PASS'
contains
  subroutine assert_true(ok,msg)
    logical,intent(in)::ok;character(len=*),intent(in)::msg
    if(.not.ok)then;write(*,'(a)')'FAIL: '//msg;error stop 1;end if
  end subroutine
  subroutine assert_close(x,y,tol,msg)
    real(dp),intent(in)::x,y,tol;character(len=*),intent(in)::msg
    call assert_true(abs(x-y)<=tol,msg)
  end subroutine
end program test_weighting
