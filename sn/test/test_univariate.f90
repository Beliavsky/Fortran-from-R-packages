! SPDX-License-Identifier: GPL-2.0-only OR GPL-3.0-only
! Modern Fortran translation of computational code from the R package sn.
program test_univariate
  use sn, only : dp, pi, dsn, psn, qsn, dst, psc, qsc, student_t_pdf, &
                 sn_uv_params, uv_operational_params, dp_to_cp_sn, cp_to_dp_sn, &
                 dp_to_op_uv, op_to_dp_uv, sn_ok
  implicit none
  real(dp) :: x,p,q,mean,sd,skew
  type(sn_uv_params) :: par,back
  type(uv_operational_params) :: op
  integer :: info

  call assert_close(dsn(0.0_dp,alpha=0.0_dp),1.0_dp/sqrt(2.0_dp*pi),2.0e-12_dp,'dsn normal')
  call assert_close(psn(0.0_dp,alpha=2.0_dp),0.5_dp-atan(2.0_dp)/pi,3.0e-10_dp,'psn zero')
  x=0.73_dp
  p=psn(x,alpha=1.7_dp)
  q=qsn(p,alpha=1.7_dp,info=info)
  call assert_true(info==sn_ok,'qsn status')
  call assert_close(q,x,2.0e-7_dp,'qsn inversion')
  call assert_close(dst(0.4_dp,alpha=0.0_dp,nu=7.0_dp),student_t_pdf(0.4_dp,7.0_dp),2.0e-13_dp,'dst t')
  p=0.31_dp
  q=qsc(p,alpha=-1.2_dp,info=info)
  call assert_close(psc(q,alpha=-1.2_dp),p,2.0e-12_dp,'skew cauchy inversion')

  par=sn_uv_params(1.2_dp,2.3_dp,-1.4_dp,0.0_dp)
  call dp_to_cp_sn(par,mean,sd,skew,info)
  call assert_true(info==sn_ok,'dp2cp status')
  call cp_to_dp_sn(mean,sd,skew,back,info)
  call assert_true(info==sn_ok,'cp2dp status')
  call assert_close(back%xi,par%xi,2.0e-9_dp,'cp xi')
  call assert_close(back%omega,par%omega,2.0e-9_dp,'cp omega')
  call assert_close(back%alpha,par%alpha,2.0e-8_dp,'cp alpha')
  op=dp_to_op_uv(par)
  back=op_to_dp_uv(op)
  call assert_close(back%omega,par%omega,2.0e-12_dp,'op omega')
  call assert_close(back%alpha,par%alpha,2.0e-12_dp,'op alpha')
  print '(a)','test_univariate: PASS'
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
end program test_univariate
