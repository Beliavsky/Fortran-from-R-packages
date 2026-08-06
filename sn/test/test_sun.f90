! SPDX-License-Identifier: GPL-2.0-only OR GPL-3.0-only
! Modern Fortran translation of computational code from the R package sn.
program test_sun
  use iso_fortran_env, only : int64
  use sn, only : dp, sn_mv_params, sun_params, sn_to_sun, dmsn, pmsn, dsun, psun, &
                 rsun, sn_rng_state, sn_ok, marginal_sun, affine_transform_sun
  implicit none
  type(sn_mv_params) :: sp
  type(sun_params) :: su,mar,aff
  type(sn_rng_state) :: rng
  real(dp),allocatable :: xdraw(:,:)
  real(dp) :: x(2),a(1),t(1,2)
  integer :: info

  sp%xi=[0.1_dp,-0.2_dp]
  sp%omega=reshape([1.0_dp,0.25_dp,0.25_dp,1.44_dp],[2,2])
  sp%alpha=[1.2_dp,-0.7_dp]
  sp%tau=0.0_dp
  call sn_to_sun(sp,su,info)
  call assert_true(info==sn_ok,'sn_to_sun')
  x=[0.3_dp,0.4_dp]
  call assert_close(dsun(x,su,info=info),dmsn(x,sp,info=info),2.0e-10_dp,'SUN density')
  call assert_close(psun(x,su,info=info,samples=32768),pmsn(x,sp,info=info,samples=32768),5.0e-5_dp,'SUN cdf')
  call marginal_sun(su,[1],mar,info)
  call assert_true(info==sn_ok .and. mar%dimension()==1,'SUN marginal')
  a=[0.5_dp]; t=reshape([1.0_dp,2.0_dp],[1,2])
  call affine_transform_sun(su,a,t,aff,info)
  call assert_true(info==sn_ok .and. aff%dimension()==1,'SUN affine')
  call rng%seed(1234_int64)
  call rsun(rng,1000,su,xdraw,info)
  call assert_true(info==sn_ok .and. size(xdraw,1)==1000 .and. size(xdraw,2)==2,'SUN random')
  print '(a)','test_sun: PASS'
contains
  subroutine assert_close(a,b,tol,msg)
    real(dp),intent(in)::a,b,tol
    character(len=*),intent(in)::msg
    if(abs(a-b)>tol*max(1.0_dp,abs(a),abs(b))) then
      write(*,*) trim(msg),a,b
      error stop 1
    end if
  end subroutine
  subroutine assert_true(ok,msg)
    logical,intent(in)::ok
    character(len=*),intent(in)::msg
    if(.not.ok) then
      write(*,*) trim(msg)
      error stop 1
    end if
  end subroutine
end program test_sun
