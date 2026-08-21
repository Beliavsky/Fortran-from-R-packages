! SPDX-License-Identifier: GPL-2.0-or-later
module gb2_likelihood
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use gb2_kinds, only : dp, pi
  use gb2_special, only : log_beta, digamma_fn, trigamma_fn, quiet_nan
  use gb2_optimizer, only : optimization_result, bfgs_minimize
  implicit none
  private
  public :: logf_gb2, dlogf_gb2, d2logf_gb2, loglik_gb2, scores_gb2, info_gb2
  public :: fisk_start, prof_gb2, profile_loglik_gb2, profile_scores_gb2
  public :: fit_gb2_full, fit_gb2_profile
  real(dp), allocatable, save :: ctx_x(:),ctx_w(:)
contains
  pure real(dp) function logf_gb2(x,shape1,scale,shape2,shape3) result(v)
    real(dp), intent(in) :: x,shape1,scale,shape2,shape3
    real(dp) :: t,lx,l1p
    if(x<=0.0_dp .or. shape1<=0.0_dp .or. scale<=0.0_dp .or. &
      shape2<=0.0_dp .or. shape3<=0.0_dp) then
      v=-huge(1.0_dp)
      return
    end if
    lx=log(x/scale)
    t=shape1*lx
    if(t>40.0_dp) then
    l1p=t+log(1.0_dp+exp(-t))
    else
    l1p=log(1.0_dp+exp(t))
    end if
    v=log(shape1/scale)-log_beta(shape2,shape3)+(shape1*shape2-1.0_dp)*lx-(shape2+shape3)*l1p
  end function logf_gb2

  pure subroutine dlogf_gb2(x,shape1,scale,shape2,shape3,g)
    real(dp), intent(in) :: x,shape1,scale,shape2,shape3
    real(dp), intent(out) :: g(4)
    real(dp) :: logy,h1,l1p
    logy=shape1*log(x/scale)
    if(logy>=0.0_dp) then
      h1=1.0_dp/(1.0_dp+exp(-logy))
      l1p=logy+log(1.0_dp+exp(-logy))
    else
      h1=exp(logy)/(1.0_dp+exp(logy))
      l1p=log(1.0_dp+exp(logy))
    end if
    g(1)=(1.0_dp+shape2*logy-(shape2+shape3)*h1*logy)/shape1
    g(2)=shape1*(-shape2+(shape2+shape3)*h1)/scale
    g(3)=digamma_fn(shape2+shape3)-digamma_fn(shape2)+logy-l1p
    g(4)=digamma_fn(shape2+shape3)-digamma_fn(shape3)-l1p
  end subroutine dlogf_gb2

  pure subroutine d2logf_gb2(x,shape1,scale,shape2,shape3,d)
    real(dp), intent(in) :: x,shape1,scale,shape2,shape3
    real(dp), intent(out) :: d(4,4)
    real(dp) :: y,logy,h1,h2,h3,h4,trip,triq,trippq
    logy=shape1*log(x/scale)
    if(logy>700.0_dp) then
    y=huge(1.0_dp)
    h1=1.0_dp
    h2=0.0_dp
    else
    y=exp(logy)
    h1=y/(1.0_dp+y)
    h2=logy/(1.0_dp+y)
    end if
    h3=h1*h2
    h4=h3*logy
    trip=trigamma_fn(shape2)
    triq=trigamma_fn(shape3)
    trippq=trigamma_fn(shape2+shape3)
    d=0.0_dp
    d(1,1)=-1.0_dp/shape1**2-(shape2+shape3)*h4/shape1**2
    d(1,2)=-shape2/scale+(shape2+shape3)*(h1+h3)/scale
    d(1,3)=h2/shape1
    d(1,4)=-(logy/shape1)*h1
    d(2,1)=d(1,2)
    if(y>=huge(1.0_dp)/2.0_dp) then
      d(2,2)=shape1*shape2/scale**2-shape1*(shape2+shape3)*h1/scale**2
    else
      d(2,2)=shape1*shape2/scale**2-shape1*(shape2+shape3)*(h1+shape1*h1/(y+1.0_dp))/scale**2
    end if
    d(2,3)=-shape1/(scale*(y+1.0_dp))
    d(2,4)=shape1*h1/scale
    d(3,1)=d(1,3)
    d(3,2)=d(2,3)
    d(3,3)=trippq-trip
    d(3,4)=trippq
    d(4,1)=d(1,4)
    d(4,2)=d(2,4)
    d(4,3)=d(3,4)
    d(4,4)=trippq-triq
  end subroutine d2logf_gb2

  real(dp) function loglik_gb2(x,par,w,hs) result(v)
    real(dp), intent(in) :: x(:),par(4)
    real(dp), intent(in), optional :: w(:),hs(:)
    real(dp) :: wi,sumw
    integer :: i
    if(any(par<=0.0_dp) .or. any(x<=0.0_dp)) then
    v=-huge(1.0_dp)
    return
    end if
    sumw=0.0_dp
    v=0.0_dp
    do i=1,size(x)
      wi=1.0_dp
      if(present(w)) wi=w(i)
      if(present(hs)) wi=wi*hs(i)
      sumw=sumw+wi
      v=v+wi*logf_gb2(x(i),par(1),par(2),par(3),par(4))
    end do
    if(sumw<=0.0_dp) then
    v=quiet_nan()
    else
    v=v/sumw
    end if
  end function loglik_gb2

  subroutine scores_gb2(x,par,g,w,hs)
    real(dp), intent(in) :: x(:),par(4)
    real(dp), intent(out) :: g(4)
    real(dp), intent(in), optional :: w(:),hs(:)
    real(dp) :: gi(4),wi,sumw
    integer :: i
    g=0.0_dp
    sumw=0.0_dp
    do i=1,size(x)
      wi=1.0_dp
      if(present(w)) wi=w(i)
      if(present(hs)) wi=wi*hs(i)
      call dlogf_gb2(x(i),par(1),par(2),par(3),par(4),gi)
      g=g+wi*gi
      sumw=sumw+wi
    end do
    g=g/sumw
  end subroutine scores_gb2

  pure subroutine info_gb2(shape1,scale,shape2,shape3,info)
    real(dp), intent(in) :: shape1,scale,shape2,shape3
    real(dp), intent(out) :: info(4,4)
    real(dp) :: psipq,trip,triq,trippq,tripq
    psipq=digamma_fn(shape2)-digamma_fn(shape3)
    trip=trigamma_fn(shape2)
    triq=trigamma_fn(shape3)
    trippq=trigamma_fn(shape2+shape3)
    tripq=trip+triq
    info=0.0_dp
    info(1,1)=(1.0_dp+(shape2*shape3/(1.0_dp+shape2+shape3))* &
      (tripq+(psipq-(shape2-shape3)/(shape2*shape3))**2 &
      -(shape2**2+shape3**2)/(shape2*shape3)**2))/shape1**2
    info(1,2)=(shape2-shape3-shape2*shape3*psipq)/(scale*(1.0_dp+shape2+shape3))
    info(2,1)=info(1,2)
    info(2,2)=shape1**2*shape2*shape3/(scale**2*(1.0_dp+shape2+shape3))
    info(2,3)=shape1*shape3/(scale*(shape2+shape3))
    info(3,2)=info(2,3)
    info(2,4)=-shape1*shape2/(scale*(shape2+shape3))
    info(4,2)=info(2,4)
    info(1,3)=-(shape3*psipq-1.0_dp)/(shape1*(shape2+shape3))
    info(3,1)=info(1,3)
    info(3,3)=trip-trippq
    info(1,4)=(shape2*psipq+1.0_dp)/(shape1*(shape2+shape3))
    info(4,1)=info(1,4)
    info(3,4)=-trippq
    info(4,3)=info(3,4)
    info(4,4)=triq-trippq
  end subroutine info_gb2

  subroutine fisk_start(z,par,w,hs)
    real(dp), intent(in) :: z(:)
    real(dp), intent(out) :: par(4)
    real(dp), intent(in), optional :: w(:),hs(:)
    real(dp) :: sw,mlz,vlz,wi,lz
    integer :: i
    sw=0.0_dp
    mlz=0.0_dp
    do i=1,size(z)
    wi=1.0_dp
    if(present(w)) wi=w(i)
    if(present(hs)) wi=wi*hs(i)
    sw=sw+wi
    mlz=mlz+wi*log(z(i))
    end do
    mlz=mlz/sw
    vlz=0.0_dp
    do i=1,size(z)
    wi=1.0_dp
    if(present(w)) wi=w(i)
    if(present(hs)) wi=wi*hs(i)
    lz=log(z(i))
    vlz=vlz+wi*(lz-mlz)**2
    end do
    vlz=vlz/sw
    par=[pi/sqrt(3.0_dp*vlz),exp(mlz),1.0_dp,1.0_dp]
  end subroutine fisk_start

  subroutine prof_gb2(x,shape1,scale,pars,w)
    real(dp), intent(in) :: x(:),shape1,scale
    real(dp), intent(out) :: pars(6)
    real(dp), intent(in), optional :: w(:)
    real(dp) :: sw,wi,y,ly,slog,sloga,slogaa,r,s
    integer :: i
    sw=0.0_dp
    slog=0.0_dp
    sloga=0.0_dp
    slogaa=0.0_dp
    r=0.0_dp
    do i=1,size(x)
      wi=1.0_dp
      if(present(w)) wi=w(i)
      y=(x(i)/scale)**shape1
      ly=log(y)
      sw=sw+wi
      slog=slog+wi*ly
      sloga=sloga+wi*log(1.0_dp+y)
      slogaa=slogaa+wi*ly*y/(1.0_dp+y)
      r=r+wi*y/(1.0_dp+y)
    end do
    slog=slog/sw
    sloga=sloga/sw
    slogaa=slogaa/sw
    r=r/sw
    s=1.0_dp/(slogaa-r*slog)
    pars=[r,s,r*s,(1.0_dp-r)*s,slog,sloga]
  end subroutine prof_gb2

  real(dp) function profile_loglik_gb2(x,shape1,scale,w) result(v)
    real(dp), intent(in) :: x(:),shape1,scale
    real(dp), intent(in), optional :: w(:)
    real(dp) :: p(6)
    call prof_gb2(x,shape1,scale,p,w)
    if(p(3)<=0.0_dp .or. p(4)<=0.0_dp) then
    v=-huge(1.0_dp)
    return
    end if
    v=-log_beta(p(3),p(4))+log(shape1/scale)+(p(3)-1.0_dp/shape1)*p(5)-(p(3)+p(4))*p(6)
  end function profile_loglik_gb2

  subroutine profile_scores_gb2(x,shape1,scale,g,w)
    real(dp), intent(in) :: x(:),shape1,scale
    real(dp), intent(out) :: g(2)
    real(dp), intent(in), optional :: w(:)
    real(dp) :: p(6),sw,wi,y,ly,drda,dsda,drdb,dsdb,dldr,dlds
    integer :: i
    call prof_gb2(x,shape1,scale,p,w)
    sw=0.0_dp
    drda=0.0_dp
    dsda=0.0_dp
    drdb=0.0_dp
    dsdb=0.0_dp
    do i=1,size(x)
      wi=1.0_dp
      if(present(w)) wi=w(i)
      y=(x(i)/scale)**shape1
      ly=log(y)
      sw=sw+wi
      drda=drda+wi*y*ly/(1.0_dp+y)**2
      dsda=dsda+wi*(y*ly/(1.0_dp+y)**2+y/(1.0_dp+y))*(ly-p(5))
      drdb=drdb+wi*y/(1.0_dp+y)**2
      dsdb=dsdb+wi*y*(ly-p(5))/(1.0_dp+y)**2
    end do
    drda=drda/(shape1*sw)
    dsda=-p(2)**2*dsda/(shape1*sw)
    drdb=-shape1*drdb/(scale*sw)
    dsdb=p(2)**2*shape1*dsdb/(scale*sw)
    dldr=p(2)*(p(5)-digamma_fn(p(3))+digamma_fn(p(4)))
    dlds=p(1)*(p(5)-digamma_fn(p(3)))+digamma_fn(p(2))-(1.0_dp-p(1))*digamma_fn(p(4))-p(6)
    g=[dldr*drda+dlds*dsda,dldr*drdb+dlds*dsdb]
  end subroutine profile_scores_gb2

  subroutine fit_gb2_full(x,result,w,hs,maxiter,tol)
    real(dp), intent(in) :: x(:)
    type(optimization_result), intent(out) :: result
    real(dp), intent(in), optional :: w(:),hs(:)
    integer, intent(in), optional :: maxiter
    real(dp), intent(in), optional :: tol
    type(optimization_result) :: raw
    real(dp) :: p0(4)
    integer :: nm
    real(dp) :: tt
    call set_context(x,w,hs)
    call fisk_start(x,p0,ctx_w)
    nm=500
    if(present(maxiter)) nm=maxiter
    tt=1.0e-8_dp
    if(present(tol)) tt=tol
    call bfgs_minimize(full_obj,full_grad,log(p0),raw,nm,tt)
    allocate(result%par(4))
    result%par=exp(raw%par)
    result%value=raw%value
    result%iterations=raw%iterations
    result%converged=raw%converged
  end subroutine fit_gb2_full

  subroutine fit_gb2_profile(x,result,w,maxiter,tol)
    real(dp), intent(in) :: x(:)
    type(optimization_result), intent(out) :: result
    real(dp), intent(in), optional :: w(:)
    integer, intent(in), optional :: maxiter
    real(dp), intent(in), optional :: tol
    type(optimization_result) :: raw
    real(dp) :: p0(4),pp(6),theta0(2),tt
    integer :: nm
    call set_context(x,w)
    call fisk_start(x,p0,ctx_w)
    theta0=log(p0(1:2))
    nm=500
    if(present(maxiter)) nm=maxiter
    tt=1.0e-10_dp
    if(present(tol)) tt=tol
    call bfgs_minimize(prof_obj,prof_grad,theta0,raw,nm,tt)
    call prof_gb2(ctx_x,exp(raw%par(1)),exp(raw%par(2)),pp,ctx_w)
    allocate(result%par(4))
    result%par=[exp(raw%par(1)),exp(raw%par(2)),pp(3),pp(4)]
    result%value=raw%value
    result%iterations=raw%iterations
    result%converged=raw%converged
  end subroutine fit_gb2_profile

  subroutine set_context(x,w,hs)
    real(dp), intent(in) :: x(:)
    real(dp), intent(in), optional :: w(:),hs(:)
    integer :: i
    if(allocated(ctx_x)) deallocate(ctx_x)
    if(allocated(ctx_w)) deallocate(ctx_w)
    allocate(ctx_x,source=x)
    allocate(ctx_w(size(x)))
    ctx_w=1.0_dp
    if(present(w)) ctx_w=w
    if(present(hs)) then
    do i=1,size(x)
    ctx_w(i)=ctx_w(i)*hs(i)
    end do
    end if
  end subroutine set_context

  real(dp) function full_obj(theta) result(v)
    real(dp), intent(in) :: theta(:)
    real(dp) :: p(4)
    p=exp(theta)
    v=-loglik_gb2(ctx_x,p,ctx_w)
  end function full_obj
  subroutine full_grad(theta,g)
    real(dp), intent(in) :: theta(:)
    real(dp), intent(out) :: g(:)
    real(dp) :: p(4),s(4)
    p=exp(theta)
    call scores_gb2(ctx_x,p,s,ctx_w)
    g=-s*p
  end subroutine full_grad
  real(dp) function prof_obj(theta) result(v)
    real(dp), intent(in) :: theta(:)
    v=-profile_loglik_gb2(ctx_x,exp(theta(1)),exp(theta(2)),ctx_w)
  end function prof_obj
  subroutine prof_grad(theta,g)
    real(dp), intent(in) :: theta(:)
    real(dp), intent(out) :: g(:)
    real(dp) :: s(2),a,b
    a=exp(theta(1))
    b=exp(theta(2))
    call profile_scores_gb2(ctx_x,a,b,s,ctx_w)
    g=-s*[a,b]
  end subroutine prof_grad
end module gb2_likelihood
