! SPDX-License-Identifier: GPL-2.0-only OR GPL-3.0-only
! Modern Fortran translation of computational code from the R package sn.
program test_multivariate
  use iso_fortran_env, only : int64
  use sn, only : dp, sn_mv_params, st_mv_params, dmsn, pmsn, rmsn, dsn, psn, &
                 dmst, dst, sn_rng_state, sn_ok, mv_operational_params, &
                 dp_to_op_mv, op_to_dp_mv
  implicit none
  type(sn_mv_params) :: p,pback
  type(mv_operational_params) :: op
  type(st_mv_params) :: tp
  type(sn_rng_state) :: rng
  real(dp),allocatable :: draws(:,:)
  real(dp) :: x1(1),v1,v2,mu
  integer :: info

  p%xi=[0.2_dp]
  p%omega=reshape([1.44_dp],[1,1])
  p%alpha=[1.5_dp]
  p%tau=0.0_dp
  x1=[0.7_dp]
  v1=dmsn(x1,p,info=info)
  v2=dsn(x1(1),0.2_dp,1.2_dp,1.5_dp)
  call assert_true(info==sn_ok,'dmsn status')
  call assert_close(v1,v2,2.0e-12_dp,'dmsn 1d')
  call assert_close(pmsn(x1,p,info=info),psn(x1(1),0.2_dp,1.2_dp,1.5_dp),2.0e-8_dp,'pmsn 1d')
  call rng%seed(445566_int64)
  call rmsn(rng,25000,p,draws,info)
  call assert_true(info==sn_ok,'rmsn status')
  mu=sum(draws(:,1))/real(size(draws,1),dp)
  call assert_close(mu,0.2_dp+1.2_dp*(1.5_dp/sqrt(3.25_dp))*sqrt(2.0_dp/acos(-1.0_dp)),0.035_dp,'rmsn mean')

  call dp_to_op_mv(p,op,info)
  call assert_true(info==sn_ok,'mv dp2op')
  call op_to_dp_mv(op,pback,info)
  call assert_true(info==sn_ok,'mv op2dp')
  call assert_close(maxval(abs(pback%omega-p%omega)),0.0_dp,2.0e-10_dp,'mv op omega')
  call assert_close(maxval(abs(pback%alpha-p%alpha)),0.0_dp,2.0e-10_dp,'mv op alpha')

  tp%xi=p%xi; tp%omega=p%omega; tp%alpha=p%alpha; tp%nu=9.0_dp
  call assert_close(dmst(x1,tp,info=info),dst(x1(1),0.2_dp,1.2_dp,1.5_dp,9.0_dp),2.0e-12_dp,'dmst 1d')
  print '(a)','test_multivariate: PASS'
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
end program test_multivariate
