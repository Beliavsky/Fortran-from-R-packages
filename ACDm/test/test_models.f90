! SPDX-License-Identifier: GPL-3.0-or-later
program test_models
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use acdm
  implicit none
  integer :: failures, model, st, i, n, burn
  type(acd_order) :: ord
  type(rng_state) :: rng
  real(dp), allocatable :: par(:),x(:),mu(:),res(:),err(:)
  real(dp),parameter::bp(2)=[0.8_dp,1.2_dp]
  real(dp)::ll
  failures=0;n=180;burn=60
  allocate(x(n),mu(n),res(n),err(n+burn))
  do i=1,n+burn
    err(i)=0.65_dp+0.70_dp*real(mod(17*i,101),dp)/100.0_dp
  end do
  call seed_rng(rng,7781)
  do model=MODEL_ACD,MODEL_TAMACD
    ord=acd_order(1,0,1)
    if(model==MODEL_AMACD.or.model==MODEL_TAMACD)ord%r=1
    par=default_model_parameters(model,ord,1.0_dp,merge(2,0,needs_breaks(model)))
    if(needs_breaks(model))then
      call simulate_acd(n,model,ord,par,DIST_EXPONENTIAL,[real(dp)::],x,st,rng,burn=burn,errors=err,breakpoints=bp)
      call assert_true(st==ACDM_SUCCESS,'simulation model '//model_name(model))
      if(st==ACDM_SUCCESS)then
        call filter_acd(x,model,ord,par,mu,res,st,breakpoints=bp)
        call assert_true(st==ACDM_SUCCESS,'filter model '//model_name(model))
        ll=acd_loglik(x,model,ord,par,DIST_EXPONENTIAL,[real(dp)::],.true.,mu,res,st,breakpoints=bp)
      end if
    else
      call simulate_acd(n,model,ord,par,DIST_EXPONENTIAL,[real(dp)::],x,st,rng,burn=burn,errors=err)
      call assert_true(st==ACDM_SUCCESS,'simulation model '//model_name(model))
      if(st==ACDM_SUCCESS)then
        call filter_acd(x,model,ord,par,mu,res,st)
        call assert_true(st==ACDM_SUCCESS,'filter model '//model_name(model))
        ll=acd_loglik(x,model,ord,par,DIST_EXPONENTIAL,[real(dp)::],.true.,mu,res,st)
      end if
    end if
    if(st==ACDM_SUCCESS)then
      call assert_true(all(ieee_is_finite(x)).and.all(x>0.0_dp),'positive x '//model_name(model))
      call assert_true(all(ieee_is_finite(mu)).and.all(mu>0.0_dp),'positive mu '//model_name(model))
      call assert_true(ieee_is_finite(ll),'finite likelihood '//model_name(model))
    end if
  end do
  call test_acd_reference
  call test_new_days
  call test_exogenous
  if(failures>0) error stop 'test_models failed'
  print '(a)','test_models: PASS'
contains
  logical function needs_breaks(m)
    integer,intent(in)::m
    needs_breaks=m==MODEL_TACD.or.m==MODEL_TAMACD.or.m==MODEL_SNIACD.or.m==MODEL_LSNIACD
  end function
  subroutine test_acd_reference
    real(dp)::xx(6),mm(6),rr(6),pp(3),meanx
    integer::s,j
    type(acd_order)::o
    xx=[1.0_dp,2.0_dp,3.0_dp,4.0_dp,2.0_dp,1.5_dp];pp=[0.2_dp,0.3_dp,0.4_dp]
    o=acd_order(1,0,1);call filter_acd(xx,MODEL_ACD,o,pp,mm,rr,s)
    call assert_true(s==ACDM_SUCCESS,'fixed ACD status')
    meanx=sum(xx)/6.0_dp
    call assert_close(mm(1),meanx,1e-13_dp,'initial mean')
    do j=2,6
      call assert_close(mm(j),0.2_dp+0.3_dp*xx(j-1)+0.4_dp*mm(j-1),1e-13_dp,'ACD recurrence')
    end do
    call assert_close(rr(6),xx(6)/mm(6),1e-13_dp,'ACD residual')
  end subroutine
  subroutine test_new_days
    real(dp)::xx(8),mm(8),rr(8),pp(3),meanx
    integer::s,nd(1)
    type(acd_order)::o
    xx=[1._dp,1.2_dp,0.9_dp,1.1_dp,3._dp,2.8_dp,3.2_dp,3.1_dp]
    pp=[0.1_dp,0.2_dp,0.7_dp];o=acd_order(1,0,1);nd=[5]
    call filter_acd(xx,MODEL_ACD,o,pp,mm,rr,s,new_day=nd)
    meanx=sum(xx)/8.0_dp
    call assert_close(mm(5),meanx,1e-13_dp,'new day reset')
  end subroutine
  subroutine test_exogenous
    real(dp)::xx(6),mm(6),rr(6),pp(4),exo(6,1),base
    integer::s
    type(acd_order)::o
    xx=1.0_dp;exo(:,1)=[0._dp,1._dp,2._dp,3._dp,4._dp,5._dp]
    pp=[0.2_dp,0.1_dp,0.5_dp,0.05_dp];o=acd_order(1,0,1)
    call filter_acd(xx,MODEL_ACD,o,pp,mm,rr,s,exogenous=exo)
    call assert_true(s==ACDM_SUCCESS,'exogenous status')
    base=0.2_dp+0.1_dp*xx(2)+0.5_dp*mm(2)+0.05_dp*exo(3,1)
    call assert_close(mm(3),base,1e-13_dp,'exogenous recurrence')
  end subroutine
  subroutine assert_close(a,b,tol,label)
    real(dp),intent(in)::a,b,tol
    character(*),intent(in)::label
    if(abs(a-b)>tol*max(1._dp,abs(b)))then
      failures=failures+1;print *, 'FAIL ',trim(label),a,b
    end if
  end subroutine
  subroutine assert_true(ok,label)
    logical,intent(in)::ok
    character(*),intent(in)::label
    if(.not.ok)then;failures=failures+1;print *,'FAIL ',trim(label);end if
  end subroutine
end program
