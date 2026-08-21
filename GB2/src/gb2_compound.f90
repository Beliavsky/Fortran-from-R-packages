! SPDX-License-Identifier: GPL-2.0-or-later
module gb2_compound
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use gb2_kinds, only : dp
  use gb2_special, only : gamma_cdf, gamma_quantile, adaptive_integral, quiet_nan
  use gb2_distribution, only : dgb2
  use gb2_moments, only : moment_gb2
  use gb2_optimizer, only : optimization_result, bfgs_minimize
  implicit none
  private
  public :: fg_cgb2, dl_cgb2, pl_cgb2, dcgb2, pcgb2, prcgb2
  public :: component_moments_cgb2, moment_cgb2, incomplete_moment_cgb2
  public :: qcgb2, rcgb2, vofp_cgb2, pofv_cgb2, loglik_cgb2, scores_cgb2, scoreu_cgb2
  public :: hess_cgb2, fit_cgb2
  real(dp), allocatable, save :: ctx_fac(:,:),ctx_w(:)
  real(dp), allocatable, save :: int_pl0(:),int_pl(:)
  real(dp), save :: int_shape1=1.0_dp,int_scale=1.0_dp,int_shape2=1.0_dp,int_shape3=1.0_dp
  integer, save :: int_component=1
  character(len=1), save :: int_decomp='r'
contains
  subroutine validate_mix(pl0,pl)
    real(dp), intent(in) :: pl0(:),pl(:)
    if(size(pl0)/=size(pl) .or. size(pl0)<2) error stop 'compound GB2: incompatible mixture vectors'
    if(any(pl0<=0.0_dp) .or. any(pl<0.0_dp)) error stop 'compound GB2: probabilities must be nonnegative (pl0 strictly positive)'
    if(abs(sum(pl0)-1.0_dp)>1.0e-8_dp .or. abs(sum(pl)-1.0_dp)>1.0e-8_dp) error stop 'compound GB2: probabilities must sum to one'
  end subroutine validate_mix

  subroutine fg_cgb2(x,shape1,scale,shape2,shape3,pl0,fac,decomp)
    real(dp), intent(in) :: x(:),shape1,scale,shape2,shape3,pl0(:)
    real(dp), intent(out) :: fac(:,:)
    character(len=*), intent(in), optional :: decomp
    character(len=1) :: dec
    real(dp) :: sh,a0,pq,t,u1,u2,cum,den,num
    integer :: i,l,nc
    nc=size(pl0)
    if(any(shape(fac)/=[size(x),nc])) error stop 'fg_cgb2: fac shape mismatch'
    dec='r'
    if(present(decomp)) dec=decomp(1:1)
    if(dec=='r' .or. dec=='R') then
      sh=shape3
      a0=shape1
    else if(dec=='l' .or. dec=='L') then
      sh=shape2
      a0=-shape1
    else
      error stop 'fg_cgb2: decomp must be r or l'
    end if
    pq=shape2+shape3
    cum=0.0_dp
    u1=0.0_dp
    do l=1,nc
      cum=cum+pl0(l)
      if(l==nc) cum=1.0_dp
      u2=gamma_quantile(cum,sh)
      den=gamma_interval(u1,u2,sh)
      do i=1,size(x)
        if(x(i)<=0.0_dp) then
          if(a0>0.0_dp) then
          t=1.0_dp
          else
          t=huge(1.0_dp)
          end if
        else
          t=1.0_dp+exp(clamp_exp_arg(a0*log(x(i)/scale)))
        end if
        num=gamma_interval_scaled(u1,u2,t,pq)
        fac(i,l)=num/den
      end do
      u1=u2
    end do
  end subroutine fg_cgb2

  pure real(dp) function clamp_exp_arg(z) result(v)
    real(dp), intent(in) :: z
    v=max(log(tiny(1.0_dp)),min(log(huge(1.0_dp))-2.0_dp,z))
  end function clamp_exp_arg

  pure real(dp) function gamma_interval(lo,hi,a) result(v)
    real(dp), intent(in) :: lo,hi,a
    if(hi>=huge(1.0_dp)/2.0_dp) then
    v=1.0_dp-gamma_cdf(lo,a)
    else
    v=gamma_cdf(hi,a)-gamma_cdf(lo,a)
    end if
  end function gamma_interval

  pure real(dp) function gamma_interval_scaled(lo,hi,t,a) result(v)
    real(dp), intent(in) :: lo,hi,t,a
    real(dp) :: l2,h2
    if(lo<=0.0_dp) then
    l2=0.0_dp
    else if(t>huge(1.0_dp)/lo) then
    l2=huge(1.0_dp)
    else
    l2=t*lo
    end if
    if(hi>=huge(1.0_dp)/2.0_dp .or. t>huge(1.0_dp)/max(hi,tiny(1.0_dp))) then
    h2=huge(1.0_dp)
    else
    h2=t*hi
    end if
    v=gamma_interval(l2,h2,a)
  end function gamma_interval_scaled

  subroutine dl_cgb2(x,shape1,scale,shape2,shape3,pl0,dens,decomp)
    real(dp), intent(in) :: x(:),shape1,scale,shape2,shape3,pl0(:)
    real(dp), intent(out) :: dens(:,:)
    character(len=*), intent(in), optional :: decomp
    real(dp), allocatable :: fac(:,:)
    integer :: i,l
    allocate(fac(size(x),size(pl0)))
    call fg_cgb2(x,shape1,scale,shape2,shape3,pl0,fac,decomp)
    do l=1,size(pl0)
    do i=1,size(x)
    dens(i,l)=dgb2(x(i),shape1,scale,shape2,shape3)*fac(i,l)
    end do
    end do
  end subroutine dl_cgb2

  subroutine pl_cgb2(y,shape1,scale,shape2,shape3,pl0,fl,decomp,tol)
    real(dp), intent(in) :: y(:),shape1,scale,shape2,shape3,pl0(:)
    real(dp), intent(out) :: fl(:,:)
    character(len=*), intent(in), optional :: decomp
    real(dp), intent(in), optional :: tol
    character(len=1) :: dec
    integer :: i,l
    if(any(shape(fl)/=[size(y),size(pl0)])) error stop 'pl_cgb2: shape mismatch'
    dec='r'
    if(present(decomp)) dec=decomp(1:1)
    do i=1,size(y)
      if(.not.ieee_is_finite(y(i))) then
        fl(i,:)=1.0_dp
      else
        do l=1,size(pl0)
        fl(i,l)=component_cdf(y(i),l,shape1,scale,shape2,shape3,pl0,dec,tol)
        end do
      end if
    end do
  end subroutine pl_cgb2

  real(dp) function dcgb2(x,shape1,scale,shape2,shape3,pl0,pl,decomp) result(v)
    real(dp), intent(in) :: x,shape1,scale,shape2,shape3,pl0(:),pl(:)
    character(len=*), intent(in), optional :: decomp
    real(dp), allocatable :: f(:,:)
    call validate_mix(pl0,pl)
    if(.not.ieee_is_finite(x)) then
    v=0.0_dp
    return
    end if
    allocate(f(1,size(pl0)))
    call dl_cgb2([x],shape1,scale,shape2,shape3,pl0,f,decomp)
    v=dot_product(f(1,:),pl)
  end function dcgb2

  subroutine set_integral_context(shape1,scale,shape2,shape3,pl0,decomp,pl)
    real(dp), intent(in) :: shape1,scale,shape2,shape3,pl0(:)
    character(len=*), intent(in) :: decomp
    real(dp), intent(in), optional :: pl(:)
    int_shape1=shape1
    int_scale=scale
    int_shape2=shape2
    int_shape3=shape3
    int_decomp=decomp(1:1)
    if(allocated(int_pl0)) deallocate(int_pl0)
    allocate(int_pl0(size(pl0)))
    int_pl0=pl0
    if(allocated(int_pl)) deallocate(int_pl)
    if(present(pl)) then
      allocate(int_pl(size(pl)))
      int_pl=pl
    end if
  end subroutine set_integral_context

  real(dp) function mixture_integrand(x) result(f)
    real(dp), intent(in) :: x
    f=dcgb2(x,int_shape1,int_scale,int_shape2,int_shape3,int_pl0,int_pl,int_decomp)
  end function mixture_integrand

  real(dp) function component_integrand(x) result(f)
    real(dp), intent(in) :: x
    real(dp) :: fa(1,size(int_pl0))
    call dl_cgb2([x],int_shape1,int_scale,int_shape2,int_shape3,int_pl0,fa,int_decomp)
    f=fa(1,int_component)
  end function component_integrand

  real(dp) function pcgb2(y,shape1,scale,shape2,shape3,pl0,pl,decomp,tol) result(v)
    real(dp), intent(in) :: y,shape1,scale,shape2,shape3,pl0(:),pl(:)
    character(len=*), intent(in), optional :: decomp
    real(dp), intent(in), optional :: tol
    real(dp) :: tt
    character(len=1) :: dec
    call validate_mix(pl0,pl)
    if(y<=0.0_dp) then
    v=0.0_dp
    return
    end if
    if(.not.ieee_is_finite(y)) then
    v=1.0_dp
    return
    end if
    tt=1.0e-8_dp
    if(present(tol)) tt=tol
    dec='r'
    if(present(decomp)) dec=decomp(1:1)
    call set_integral_context(shape1,scale,shape2,shape3,pl0,dec,pl)
    v=adaptive_integral(mixture_integrand,0.0_dp,y,tt)
    v=max(0.0_dp,min(1.0_dp,v))
  end function pcgb2

  real(dp) function prcgb2(y1,y2,shape1,scale,shape2,shape3,pl0,pl,decomp,tol) result(v)
    real(dp), intent(in) :: y1,y2,shape1,scale,shape2,shape3,pl0(:),pl(:)
    character(len=*), intent(in), optional :: decomp
    real(dp), intent(in), optional :: tol
    real(dp) :: lo,hi
    lo=min(y1,y2)
    hi=max(y1,y2)
    if(abs(hi-lo)<=tiny(1.0_dp)) then
    v=0.0_dp
    return
    end if
    if(.not.ieee_is_finite(hi)) then
      v=1.0_dp-pcgb2(lo,shape1,scale,shape2,shape3,pl0,pl,decomp,tol)
    else
      v=pcgb2(hi,shape1,scale,shape2,shape3,pl0,pl,decomp,tol) &
        -pcgb2(lo,shape1,scale,shape2,shape3,pl0,pl,decomp,tol)
    end if
    v=max(0.0_dp,min(1.0_dp,v))
  end function prcgb2

  subroutine component_moments_cgb2(k,shape1,scale,shape2,shape3,pl0,mom,decomp)
    real(dp), intent(in) :: k,shape1,scale,shape2,shape3,pl0(:)
    real(dp), intent(out) :: mom(:)
    character(len=*), intent(in), optional :: decomp
    character(len=1) :: dec
    real(dp) :: sh,a0,shk,u1,u2,cum,eg,den,num
    integer :: l,nc
    nc=size(pl0)
    if(size(mom)/=nc) error stop 'component_moments_cgb2: shape mismatch'
    dec='r'
    if(present(decomp)) dec=decomp(1:1)
    if(dec=='r' .or. dec=='R') then
    sh=shape3
    a0=shape1
    else
    sh=shape2
    a0=-shape1
    end if
    shk=sh-k/a0
    if(shk<=0.0_dp) then
    mom=quiet_nan()
    return
    end if
    eg=moment_gb2(k,shape1,scale,shape2,shape3)
    cum=0.0_dp
    u1=0.0_dp
    do l=1,nc
      cum=cum+pl0(l)
      if(l==nc) cum=1.0_dp
      u2=gamma_quantile(cum,sh)
      den=gamma_interval(u1,u2,sh)
      num=gamma_interval(u1,u2,shk)
      mom(l)=eg*num/den
      u1=u2
    end do
  end subroutine component_moments_cgb2

  real(dp) function moment_cgb2(k,shape1,scale,shape2,shape3,pl0,pl,decomp) result(v)
    real(dp), intent(in) :: k,shape1,scale,shape2,shape3,pl0(:),pl(:)
    character(len=*), intent(in), optional :: decomp
    real(dp), allocatable :: m(:)
    call validate_mix(pl0,pl)
    allocate(m(size(pl0)))
    call component_moments_cgb2(k,shape1,scale,shape2,shape3,pl0,m,decomp)
    v=dot_product(pl,m)
  end function moment_cgb2

  real(dp) function incomplete_moment_cgb2(x,k,shape1,scale,shape2,shape3,pl0,pl,decomp,tol) result(v)
    real(dp), intent(in) :: x,k,shape1,scale,shape2,shape3,pl0(:),pl(:)
    character(len=*), intent(in), optional :: decomp
    real(dp), intent(in), optional :: tol
    real(dp), allocatable :: ek(:),ppl0(:),fk(:)
    real(dp) :: sh,a0,shk,u1,u2,cum,denom,num
    character(len=1) :: dec
    integer :: l,nc
    call validate_mix(pl0,pl)
    nc=size(pl0)
    allocate(ek(nc),ppl0(nc),fk(nc))
    dec='r'
    if(present(decomp)) dec=decomp(1:1)
    if(dec=='r' .or. dec=='R') then
    sh=shape3
    a0=shape1
    else
    sh=shape2
    a0=-shape1
    end if
    shk=sh-k/a0
    if(shk<=0.0_dp) then
    v=quiet_nan()
    return
    end if
    call component_moments_cgb2(k,shape1,scale,shape2,shape3,pl0,ek,dec)
    cum=0.0_dp
    u1=0.0_dp
    do l=1,nc
      cum=cum+pl0(l)
      if(l==nc) cum=1.0_dp
      u2=gamma_quantile(cum,sh)
      ppl0(l)=gamma_interval(u1,u2,shk)
      u1=u2
    end do
    ! ppl0 are the interval probabilities under the shifted gamma law and sum to one.
    do l=1,nc
      fk(l)=component_cdf(x,l,shape1,scale,shape2+k/shape1,shape3-k/shape1,ppl0,dec,tol)
    end do
    num=dot_product(pl,ek*fk)
    denom=dot_product(pl,ek)
    v=num/denom
  end function incomplete_moment_cgb2

  real(dp) function component_cdf(x,l,shape1,scale,shape2,shape3,pl0,decomp,tol) result(v)
    real(dp), intent(in) :: x,shape1,scale,shape2,shape3,pl0(:)
    integer, intent(in) :: l
    character(len=*), intent(in) :: decomp
    real(dp), intent(in), optional :: tol
    real(dp) :: tt
    if(x<=0.0_dp) then
    v=0.0_dp
    return
    end if
    tt=1.0e-8_dp
    if(present(tol)) tt=tol
    call set_integral_context(shape1,scale,shape2,shape3,pl0,decomp)
    int_component=l
    v=adaptive_integral(component_integrand,0.0_dp,x,tt)
  end function component_cdf

  real(dp) function qcgb2(prob,shape1,scale,shape2,shape3,pl0,pl,decomp,tol,maxiter) result(q)
    real(dp), intent(in) :: prob,shape1,scale,shape2,shape3,pl0(:),pl(:)
    character(len=*), intent(in), optional :: decomp
    real(dp), intent(in), optional :: tol
    integer, intent(in), optional :: maxiter
    real(dp) :: lo,hi,f,pcur,dens,trial,tt
    integer :: it,mi
    character(len=1) :: dec
    if(prob<0.0_dp .or. prob>1.0_dp) then
    q=quiet_nan()
    return
    end if
    if(prob<=0.0_dp) then
    q=0.0_dp
    return
    end if
    if(prob>=1.0_dp) then
    q=huge(1.0_dp)
    return
    end if
    tt=1.0e-8_dp
    if(present(tol)) tt=tol
    mi=80
    if(present(maxiter)) mi=maxiter
    dec='r'
    if(present(decomp)) dec=decomp(1:1)
    lo=0.0_dp
    hi=max(scale,moment_cgb2(1.0_dp,shape1,scale,shape2,shape3,pl0,pl,dec))
    if(.not.ieee_is_finite(hi) .or. hi<=0.0_dp) hi=scale
    do while(pcgb2(hi,shape1,scale,shape2,shape3,pl0,pl,dec,tt)<prob)
      hi=2.0_dp*hi
      if(hi>huge(1.0_dp)/4.0_dp) exit
    end do
    q=0.5_dp*(lo+hi)
    do it=1,mi
      pcur=pcgb2(q,shape1,scale,shape2,shape3,pl0,pl,dec,tt)
      f=pcur-prob
      if(f>0.0_dp) then
      hi=q
      else
      lo=q
      end if
      if(abs(f)<max(tt,1.0e-10_dp) .or. hi-lo<max(tt,1.0e-10_dp)*max(1.0_dp,q)) exit
      dens=dcgb2(q,shape1,scale,shape2,shape3,pl0,pl,dec)
      if(dens>tiny(1.0_dp)) then
      trial=q-f/dens
      else
      trial=0.5_dp*(lo+hi)
      end if
      if(trial<=lo .or. trial>=hi .or. .not.ieee_is_finite(trial)) trial=0.5_dp*(lo+hi)
      q=trial
    end do
  end function qcgb2

  subroutine rcgb2(n,shape1,scale,shape2,shape3,pl0,pl,x,decomp,tol,maxiter)
    integer, intent(in) :: n
    real(dp), intent(in) :: shape1,scale,shape2,shape3,pl0(:),pl(:)
    real(dp), intent(out) :: x(:)
    character(len=*), intent(in), optional :: decomp
    real(dp), intent(in), optional :: tol
    integer, intent(in), optional :: maxiter
    real(dp) :: u
    integer :: i
    if(size(x)/=n) error stop 'rcgb2: size mismatch'
    do i=1,n
    call random_number(u)
    x(i)=qcgb2(u,shape1,scale,shape2,shape3,pl0,pl,decomp,tol,maxiter)
    end do
  end subroutine rcgb2

  subroutine vofp_cgb2(pl,vl)
    real(dp), intent(in) :: pl(:)
    real(dp), intent(out) :: vl(:)
    if(size(vl)/=size(pl)-1) error stop 'vofp_cgb2: size mismatch'
    vl=log(pl(1:size(vl))/pl(size(pl)))
  end subroutine vofp_cgb2

  subroutine pofv_cgb2(vl,pl)
    real(dp), intent(in) :: vl(:)
    real(dp), intent(out) :: pl(:)
    real(dp) :: m,s
    if(size(pl)/=size(vl)+1) error stop 'pofv_cgb2: size mismatch'
    m=max(0.0_dp,maxval(vl))
    s=exp(-m)+sum(exp(vl-m))
    pl(1:size(vl))=exp(vl-m)/s
    pl(size(pl))=exp(-m)/s
  end subroutine pofv_cgb2

  real(dp) function loglik_cgb2(fac,pl,w) result(v)
    real(dp), intent(in) :: fac(:,:),pl(:)
    real(dp), intent(in), optional :: w(:)
    real(dp) :: sw,wi,mix
    integer :: i
    sw=0.0_dp
    v=0.0_dp
    do i=1,size(fac,1)
      wi=1.0_dp
      if(present(w)) wi=w(i)
      mix=dot_product(fac(i,:),pl)
      if(mix<=0.0_dp) then
        v=-huge(1.0_dp)
        return
      end if
      sw=sw+wi
      v=v+wi*log(mix)
    end do
    v=v/sw
  end function loglik_cgb2

  subroutine scoreu_cgb2(fac,pl,u)
    real(dp), intent(in) :: fac(:,:),pl(:)
    real(dp), intent(out) :: u(:,:)
    real(dp) :: denom
    integer :: i,j,ncomp
    ncomp=size(pl)
    if(any(shape(u)/=[size(fac,1),ncomp-1])) error stop 'scoreu_cgb2: shape mismatch'
    do i=1,size(fac,1)
    denom=dot_product(fac(i,:),pl)
    do j=1,ncomp-1
    u(i,j)=pl(j)*(fac(i,j)/denom-1.0_dp)
    end do
    end do
  end subroutine scoreu_cgb2

  subroutine scores_cgb2(fac,pl,g,w)
    real(dp), intent(in) :: fac(:,:),pl(:)
    real(dp), intent(out) :: g(:)
    real(dp), intent(in), optional :: w(:)
    real(dp), allocatable :: u(:,:)
    real(dp) :: sw,wi
    integer :: i
    allocate(u(size(fac,1),size(pl)-1))
    call scoreu_cgb2(fac,pl,u)
    g=0.0_dp
    sw=0.0_dp
    do i=1,size(fac,1)
    wi=1.0_dp
    if(present(w)) wi=w(i)
    g=g+wi*u(i,:)
    sw=sw+wi
    end do
    g=g/sw
  end subroutine scores_cgb2

  subroutine hess_cgb2(u,pl,h,w)
    real(dp), intent(in) :: u(:,:),pl(:)
    real(dp), intent(out) :: h(:,:)
    real(dp), intent(in), optional :: w(:)
    real(dp), allocatable :: sumsc(:)
    real(dp) :: wi
    integer :: i,j,k,l1
    l1=size(pl)-1
    allocate(sumsc(l1))
    sumsc=0.0_dp
    h=0.0_dp
    do k=1,size(u,1)
      wi=1.0_dp
      if(present(w)) wi=w(k)
      sumsc=sumsc+wi*u(k,:)
      do j=1,l1
        do i=1,l1
          h(i,j)=h(i,j)-wi*u(k,i)*u(k,j)
        end do
      end do
    end do
    do j=1,l1
    do i=1,l1
    h(i,j)=h(i,j)-pl(i)*sumsc(j)-sumsc(i)*pl(j)
    end do
    end do
    do i=1,l1
    h(i,i)=h(i,i)+sumsc(i)
    end do
  end subroutine hess_cgb2

  subroutine fit_cgb2(fac,pl0,result,w,maxiter,tol)
    real(dp), intent(in) :: fac(:,:),pl0(:)
    type(optimization_result), intent(out) :: result
    real(dp), intent(in), optional :: w(:)
    integer, intent(in), optional :: maxiter
    real(dp), intent(in), optional :: tol
    real(dp), allocatable :: v0(:),pl(:)
    type(optimization_result) :: raw
    integer :: nm
    real(dp) :: tt
    call set_context(fac,w)
    allocate(v0(size(pl0)-1))
    call vofp_cgb2(pl0,v0)
    nm=300
    if(present(maxiter)) nm=maxiter
    tt=1.0e-9_dp
    if(present(tol)) tt=tol
    call bfgs_minimize(mix_obj,mix_grad,v0,raw,nm,tt)
    allocate(pl(size(pl0)))
    call pofv_cgb2(raw%par,pl)
    allocate(result%par(size(pl)))
    result%par=pl
    result%value=raw%value
    result%iterations=raw%iterations
    result%converged=raw%converged
  end subroutine fit_cgb2

  subroutine set_context(fac,w)
    real(dp), intent(in) :: fac(:,:)
    real(dp), intent(in), optional :: w(:)
    if(allocated(ctx_fac)) deallocate(ctx_fac)
    if(allocated(ctx_w)) deallocate(ctx_w)
    allocate(ctx_fac,source=fac)
    allocate(ctx_w(size(fac,1)))
    ctx_w=1.0_dp
    if(present(w)) ctx_w=w
  end subroutine set_context
  real(dp) function mix_obj(vl) result(v)
    real(dp), intent(in) :: vl(:)
    real(dp) :: pl(size(vl)+1)
    call pofv_cgb2(vl,pl)
    v=-loglik_cgb2(ctx_fac,pl,ctx_w)
  end function mix_obj
  subroutine mix_grad(vl,g)
    real(dp), intent(in) :: vl(:)
    real(dp), intent(out) :: g(:)
    real(dp) :: pl(size(vl)+1)
    call pofv_cgb2(vl,pl)
    call scores_cgb2(ctx_fac,pl,g,ctx_w)
    g=-g
  end subroutine mix_grad
end module gb2_compound
