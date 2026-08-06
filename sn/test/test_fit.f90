! SPDX-License-Identifier: GPL-2.0-only OR GPL-3.0-only
! Modern Fortran translation of computational code from the R package sn.
program test_fit
  use iso_fortran_env, only : int64
  use sn, only : dp, sn_rng_state, rsn, selm_result, selm_fit, predict_selm, sn_ok
  implicit none
  type(sn_rng_state) :: rng
  type(selm_result) :: fit
  real(dp),allocatable :: x(:,:),y(:),e(:),pred(:),meanpred(:)
  integer :: n,i,info
  n=140
  allocate(x(n,2),y(n),e(n))
  x(:,1)=1.0_dp
  do i=1,n
    x(i,2)=-1.0_dp+2.0_dp*real(i-1,dp)/real(n-1,dp)
  end do
  call rng%seed(20260804_int64)
  call rsn(rng,e,omega=0.65_dp,alpha=2.0_dp,info=info)
  y=1.1_dp+2.4_dp*x(:,2)+e
  call selm_fit(x,y,'SN',fit,penalty='NONE',max_iter=3500,tol=2.0e-7_dp)
  call assert_true(fit%status==sn_ok .and. fit%converged,'selm fit')
  call assert_close(fit%beta(2),2.4_dp,0.12_dp,'slope recovery')
  call assert_true(fit%omega>0.3_dp .and. fit%omega<1.2_dp,'scale range')
  call predict_selm(fit,x,pred,meanpred,info)
  call assert_true(info==sn_ok .and. size(pred)==n .and. size(meanpred)==n,'prediction')
  call assert_close(sum(fit%residuals+pred-y)/real(n,dp),0.0_dp,1.0e-10_dp,'residual identity')
  print '(a)','test_fit: PASS'
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
end program test_fit
