! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2018 Marius Hofert, Erik Hintz and Christiane Lemieux
module nvmix_skewt
  use nvmix_kinds, only : dp,i8,log_two_pi
  use nvmix_types, only : integration_control,sample_result
  use nvmix_special, only : normal_cdf,normal_quantile,inverse_gamma_quantile,student_quantile
  use nvmix_random, only : halton,seed_random,random_normal,random_gamma
  use nvmix_linalg, only : cholesky_lower,forward_solve
  implicit none
  private
  public :: dskewt,pskewt1d,qskewt1d,rskewt,dskewt_copula,rskewt_copula
contains
  real(dp) function dskewt(x,gamma,df,loc,scale,control,log_density) result(v)
    real(dp), intent(in) :: x(:),gamma(:),df,loc(:),scale(:,:)
    type(integration_control), intent(in), optional :: control
    logical, intent(in), optional :: log_density
    type(integration_control) :: ctrl
    real(dp), allocatable :: l(:,:),z(:),delta(:),logs(:)
    real(dp) :: w,logdet,mx,s
    integer :: d,n,i,j
    logical :: ok,lg
    ctrl=integration_control(); if(present(control))ctrl=control
    d=size(x); lg=.false.; if(present(log_density))lg=log_density
    if(size(gamma)/=d .or. size(loc)/=d .or. size(scale,1)/=d .or. size(scale,2)/=d .or. df<=0.0_dp)then
      if(lg)then; v=-huge(1.0_dp); else; v=0.0_dp; end if; return
    end if
    call cholesky_lower(scale,l,ok); if(.not.ok)then; if(lg)then; v=-huge(1.0_dp); else; v=0.0_dp; end if; return; end if
    logdet=0.0_dp; do j=1,d; logdet=logdet+2.0_dp*log(l(j,j)); end do
    n=max(512,ctrl%samples); allocate(z(d),delta(d),logs(n))
    do i=1,n
      w=inverse_gamma_quantile(halton(i,2),df)
      delta=(x-loc-gamma*w)/sqrt(w)
      call forward_solve(l,delta,z,ok)
      logs(i)=-0.5_dp*(real(d,dp)*log_two_pi+logdet+real(d,dp)*log(w)+dot_product(z,z))
    end do
    mx=maxval(logs); s=sum(exp(logs-mx)); v=mx+log(s/real(n,dp))
    if(.not.lg)v=exp(v)
  end function

  real(dp) function pskewt1d(q,gamma,df,loc,scale,control) result(p)
    real(dp), intent(in) :: q,gamma,df,loc,scale
    type(integration_control), intent(in), optional :: control
    type(integration_control) :: ctrl
    real(dp) :: w,s
    integer :: i,n
    ctrl=integration_control(); if(present(control))ctrl=control
    if(scale<=0.0_dp .or. df<=0.0_dp)then; p=0.0_dp; return; end if
    if(abs(gamma)<=1.0e-14_dp)then
      p=0.5_dp+sign(0.5_dp,q-loc)*(1.0_dp-2.0_dp*0.5_dp) ! overwritten below
      p=student_cdf_local((q-loc)/sqrt(scale),df); return
    end if
    n=max(256,ctrl%samples); s=0.0_dp
    do i=1,n
      w=inverse_gamma_quantile(halton(i,2),df)
      s=s+normal_cdf((q-loc-gamma*w)/sqrt(scale*w))
    end do
    p=s/real(n,dp); p=min(1.0_dp,max(0.0_dp,p))
  contains
    real(dp) function student_cdf_local(x,nu) result(value)
      use nvmix_special, only : student_cdf
      real(dp), intent(in) :: x,nu
      value=student_cdf(x,nu)
    end function
  end function

  real(dp) function qskewt1d(prob,gamma,df,loc,scale,control) result(q)
    real(dp), intent(in) :: prob,gamma,df,loc,scale
    type(integration_control), intent(in), optional :: control
    type(integration_control) :: ctrl
    real(dp), allocatable :: w(:)
    real(dp) :: lo,hi,mid,sd
    integer :: i,n
    ctrl=integration_control(); if(present(control))ctrl=control
    if(prob<=0.0_dp)then; q=-huge(1.0_dp); return; end if
    if(prob>=1.0_dp)then; q=huge(1.0_dp); return; end if
    if(abs(gamma)<=1.0e-14_dp)then
      q=loc+sqrt(scale)*student_quantile(prob,df); return
    end if
    n=max(256,ctrl%samples); allocate(w(n))
    do i=1,n; w(i)=inverse_gamma_quantile(halton(i,2),df); end do
    sd=sqrt(scale); lo=loc+gamma*max(1.0_dp,df/max(df-2.0_dp,0.1_dp))-sd
    hi=lo+2.0_dp*sd
    do while(cdf_cached(lo)>prob); lo=loc+2.0_dp*(lo-loc); end do
    do while(cdf_cached(hi)<prob); hi=loc+2.0_dp*(hi-loc); end do
    do i=1,80
      mid=0.5_dp*(lo+hi)
      if(cdf_cached(mid)<prob)then; lo=mid; else; hi=mid; end if
    end do
    q=0.5_dp*(lo+hi)
  contains
    real(dp) function cdf_cached(x) result(value)
      real(dp), intent(in) :: x
      integer :: k
      value=0.0_dp
      do k=1,n
        value=value+normal_cdf((x-loc-gamma*w(k))/sqrt(scale*w(k)))
      end do
      value=value/real(n,dp)
    end function
  end function

  function rskewt(n,gamma,df,loc,scale,seed) result(result)
    integer, intent(in) :: n
    real(dp), intent(in) :: gamma(:),df,loc(:),scale(:,:)
    integer(i8), intent(in), optional :: seed
    type(sample_result) :: result
    real(dp), allocatable :: l(:,:),z(:),y(:)
    real(dp) :: w
    integer :: i,j,d
    logical :: ok
    d=size(loc); allocate(result%x(max(0,n),d))
    if(size(gamma)/=d .or. size(scale,1)/=d .or. size(scale,2)/=d .or. n<1)then
      result%ok=.false.; result%message='invalid skew-t dimensions'; return
    end if
    call cholesky_lower(scale,l,ok)
    if(.not.ok)then
      result%ok=.false.
      result%message='scale is not positive definite'
      return
    end if
    if(present(seed))call seed_random(seed)
    allocate(z(d),y(d))
    do i=1,n
      w=1.0_dp/random_gamma(0.5_dp*df,2.0_dp/df)
      do j=1,d; z(j)=random_normal(); end do
      y=matmul(l,z)
      result%x(i,:)=loc+gamma*w+sqrt(w)*y
    end do
  end function

  real(dp) function dskewt_copula(u,gamma,df,scale,control,log_density) result(v)
    real(dp), intent(in) :: u(:),gamma(:),df,scale(:,:)
    type(integration_control), intent(in), optional :: control
    logical, intent(in), optional :: log_density
    type(integration_control) :: ctrl
    real(dp), allocatable :: x(:),loc(:),one_scale(:,:),one_gamma(:)
    real(dp) :: ld,marg
    integer :: d,j
    logical :: lg
    ctrl=integration_control(); if(present(control))ctrl=control
    d=size(u); allocate(x(d),loc(d)); loc=0.0_dp; lg=.false.; if(present(log_density))lg=log_density
    if(any(u<=0.0_dp) .or. any(u>=1.0_dp))then; if(lg)then; v=-huge(1.0_dp); else; v=0.0_dp; end if; return; end if
    do j=1,d; x(j)=qskewt1d(u(j),gamma(j),df,0.0_dp,scale(j,j),ctrl); end do
    ld=dskewt(x,gamma,df,loc,scale,ctrl,.true.)
    do j=1,d
      allocate(one_scale(1,1),one_gamma(1)); one_scale(1,1)=scale(j,j); one_gamma(1)=gamma(j)
      marg=dskewt([x(j)],one_gamma,df,[0.0_dp],one_scale,ctrl,.false.)
      ld=ld-log(max(marg,tiny(1.0_dp))); deallocate(one_scale,one_gamma)
    end do
    if(lg)then; v=ld; else; v=exp(ld); end if
  end function

  function rskewt_copula(n,gamma,df,scale,seed,control) result(result)
    integer, intent(in) :: n
    real(dp), intent(in) :: gamma(:),df,scale(:,:)
    integer(i8), intent(in), optional :: seed
    type(integration_control), intent(in), optional :: control
    type(sample_result) :: result
    type(sample_result) :: raw
    type(integration_control) :: ctrl
    real(dp), allocatable :: loc(:)
    integer :: i,j,d
    ctrl=integration_control(); if(present(control))ctrl=control
    d=size(gamma); allocate(loc(d)); loc=0.0_dp
    if(present(seed))then; raw=rskewt(n,gamma,df,loc,scale,seed); else; raw=rskewt(n,gamma,df,loc,scale); end if
    allocate(result%x(n,d)); result%ok=raw%ok; result%message=raw%message; if(.not.raw%ok)return
    do i=1,n; do j=1,d; result%x(i,j)=pskewt1d(raw%x(i,j),gamma(j),df,0.0_dp,scale(j,j),ctrl); end do; end do
  end function
end module nvmix_skewt
