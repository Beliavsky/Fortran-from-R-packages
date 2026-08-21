! SPDX-License-Identifier: GPL-2.0-or-later
module gb2_nonlinear_fit
  use gb2_kinds, only : dp
  use gb2_optimizer, only : optimization_result, bfgs_minimize
  use gb2_indicators, only : main2_gb2
  use gb2_distribution, only : qgb2
  implicit none
  private
  public :: nlsfit_gb2
  real(dp), save :: ctx_ei(4),ctx_wt(4),ctx_lo,ctx_hi
contains
  subroutine nlsfit_gb2(med,ei4,par0,result,cva,ei4w)
    real(dp), intent(in) :: med,ei4(4),par0(4)
    type(optimization_result), intent(out) :: result
    real(dp), intent(in), optional :: cva,ei4w(4)
    real(dp) :: cv,lo,hi,t0(3),pred(6),a,ap,aq,p,q,b,qmed
    type(optimization_result) :: raw
    cv=1.0_dp
    if(present(cva)) cv=cva
    lo=par0(1)*max(0.2_dp,1.0_dp-2.0_dp*cv)
    hi=par0(1)*min(2.0_dp,1.0_dp+2.0_dp*cv)
    ctx_ei=ei4
    if(present(ei4w)) then
    ctx_wt=ei4w
    else
    ctx_wt=1.0_dp/max(abs(ei4),1.0e-8_dp)
    end if
    ctx_lo=lo
    ctx_hi=hi
    t0(1)=logit((par0(1)-lo)/(hi-lo))
    t0(2)=logit((par0(1)*par0(3)-1.0_dp)/99.0_dp)
    t0(3)=logit((par0(1)*par0(4)-2.0_dp)/98.0_dp)
    call bfgs_minimize(nls_obj,nls_grad,t0,raw,1000,1.0e-7_dp)
    call unpack(raw%par,a,ap,aq)
    p=ap/a
    q=aq/a
    call main2_gb2(0.6_dp,a,1.0_dp,ap,aq,pred)
    qmed=qgb2(0.5_dp,a,1.0_dp,p,q)
    b=med/qmed
    b=max(0.01_dp,min(2.0_dp*par0(2),b))
    allocate(result%par(4))
    result%par=[a,b,p,q]
    result%value=raw%value
    result%iterations=raw%iterations
    result%converged=raw%converged
  end subroutine nlsfit_gb2
  pure real(dp) function logistic(t) result(v)
    real(dp), intent(in) :: t
    if(t>=0.0_dp) then
    v=1.0_dp/(1.0_dp+exp(-t))
    else
    v=exp(t)/(1.0_dp+exp(t))
    end if
  end function logistic
  pure real(dp) function logit(p) result(v)
    real(dp), intent(in) :: p
    real(dp) :: pp
    pp=max(1.0e-8_dp,min(1.0_dp-1.0e-8_dp,p))
    v=log(pp/(1.0_dp-pp))
  end function logit
  subroutine unpack(t,a,ap,aq)
    real(dp), intent(in) :: t(:)
    real(dp), intent(out) :: a,ap,aq
    a=ctx_lo+(ctx_hi-ctx_lo)*logistic(t(1))
    ap=1.0_dp+99.0_dp*logistic(t(2))
    aq=2.0_dp+98.0_dp*logistic(t(3))
  end subroutine unpack
  real(dp) function nls_obj(t) result(v)
    real(dp), intent(in) :: t(:)
    real(dp) :: a,ap,aq,pred(6),r(4)
    call unpack(t,a,ap,aq)
    call main2_gb2(0.6_dp,a,1.0_dp,ap,aq,pred)
    r=pred(3:6)-ctx_ei
    v=sum(ctx_wt*r*r)
  end function nls_obj
  subroutine nls_grad(t,g)
    real(dp), intent(in) :: t(:)
    real(dp), intent(out) :: g(:)
    real(dp) :: tp(size(t)),tm(size(t)),h
    integer :: j
    do j=1,size(t)
      h=1.0e-5_dp*max(1.0_dp,abs(t(j)))
      tp=t
      tm=t
      tp(j)=tp(j)+h
      tm(j)=tm(j)-h
      g(j)=(nls_obj(tp)-nls_obj(tm))/(2.0_dp*h)
    end do
  end subroutine nls_grad
end module gb2_nonlinear_fit
