! SPDX-License-Identifier: GPL-2.0-or-later
module tmvtnorm_marginals
  use mvtnorm_kinds, only : dp
  use mvtnorm_types, only : probability_control, probability_result
  use mvtnorm_probabilities, only : mvn_prob => pmvnorm, mvt_prob => pmvt
  use mvtnorm_distributions, only : dmvnorm_one
  use mvtnorm_linalg, only : inverse_spd
  use mvtnorm_special, only : normal_pdf, normal_cdf
  implicit none
  private
  public :: dtmvnorm_marginal, dtmvnorm_marginal2
  public :: ptmvnorm_marginal, qtmvnorm_marginal
  public :: ptmvt_marginal
  public :: tail_lower, tail_upper, tail_both
  integer,parameter :: tail_lower=1, tail_upper=2, tail_both=3

contains

  real(dp) function dtmvnorm_marginal(xn,nidx,mean,sigma,lower,upper,log_density,control) result(v)
    real(dp),intent(in)::xn,mean(:),sigma(:,:),lower(:),upper(:)
    integer,intent(in)::nidx
    logical,intent(in),optional::log_density
    type(probability_control),intent(in),optional::control
    integer::d,m,i,k
    integer,allocatable::idx(:)
    real(dp),allocatable::muo(:),lo(:),up(:),cm(:),cc(:,:),cross(:)
    real(dp)::den,condp,ld,varn
    logical::ll
    type(probability_result)::pr
    ll=.false.
    if(present(log_density)) ll=log_density
    d=size(mean)
    if(nidx<1 .or. nidx>d .or. xn<lower(nidx) .or. xn>upper(nidx) .or. abs(xn)>=huge(1.0_dp)/2.0_dp) then
      v=merge(-huge(1.0_dp),0.0_dp,ll)
      return
    end if
    if(present(control)) then
    pr=mvn_prob(lower,upper,mean,sigma,control)
    else
    pr=mvn_prob(lower,upper,mean,sigma)
    end if
    den=pr%value
    if(den<=0.0_dp) then
    v=merge(-huge(1.0_dp),0.0_dp,ll)
    return
    end if
    varn=sigma(nidx,nidx)
    if(d==1) then
      ld=log(normal_pdf((xn-mean(1))/sqrt(varn))/sqrt(varn))-log(den)
      v=merge(ld,exp(ld),ll)
      return
    end if
    allocate(idx(d-1))
    k=0
    do i=1,d
    if(i/=nidx) then
    k=k+1
    idx(k)=i
    end if
    end do
    m=d-1
    allocate(muo(m),lo(m),up(m),cm(m),cc(m,m),cross(m))
    muo=mean(idx)
    lo=lower(idx)
    up=upper(idx)
    cross=sigma(idx,nidx)
    cm=muo+cross*(xn-mean(nidx))/varn
    cc=sigma(idx,idx)-outer(cross,cross)/varn
    if(present(control)) then
    pr=mvn_prob(lo,up,cm,cc,control)
    else
    pr=mvn_prob(lo,up,cm,cc)
    end if
    condp=max(0.0_dp,pr%value)
    if(condp<=0.0_dp) then
    v=merge(-huge(1.0_dp),0.0_dp,ll)
    return
    end if
    ld=log(normal_pdf((xn-mean(nidx))/sqrt(varn))/sqrt(varn))+log(condp)-log(den)
    v=merge(ld,exp(ld),ll)
  end function dtmvnorm_marginal

  real(dp) function dtmvnorm_marginal2(xq,xr,q,r,mean,sigma,lower,upper,log_density,control) result(v)
    real(dp),intent(in)::xq,xr,mean(:),sigma(:,:),lower(:),upper(:)
    integer,intent(in)::q,r
    logical,intent(in),optional::log_density
    type(probability_control),intent(in),optional::control
    integer::d,m,i,k
    integer,allocatable::idx(:)
    real(dp)::den,condp,ld
    real(dp)::xs(2),mus(2),ss(2,2)
    real(dp),allocatable::so(:,:),os(:,:),oo(:,:),invss(:,:),cm(:),cc(:,:),lo(:),up(:)
    logical::ll,ok
    character(len=256)::msg
    type(probability_result)::pr
    ll=.false.
    if(present(log_density)) ll=log_density
    d=size(mean)
    if(q<1 .or. r<1 .or. q>d .or. r>d .or. q==r) then
    v=merge(-huge(1.0_dp),0.0_dp,ll)
    return
    end if
    if(xq<lower(q) .or. xq>upper(q) .or. xr<lower(r) .or. xr>upper(r) .or. &
       abs(xq)>=huge(1.0_dp)/2.0_dp .or. abs(xr)>=huge(1.0_dp)/2.0_dp) then
      v=merge(-huge(1.0_dp),0.0_dp,ll)
      return
    end if
    if(present(control)) then
    pr=mvn_prob(lower,upper,mean,sigma,control)
    else
    pr=mvn_prob(lower,upper,mean,sigma)
    end if
    den=pr%value
    if(den<=0.0_dp) then
    v=merge(-huge(1.0_dp),0.0_dp,ll)
    return
    end if
    xs=[xq,xr]
    mus=[mean(q),mean(r)]
    ss=reshape([sigma(q,q),sigma(r,q),sigma(q,r),sigma(r,r)],[2,2])
    ld=dmvnorm_one(xs,mus,ss,.true.)
    if(d==2) then
    v=merge(ld-log(den),exp(ld)/den,ll)
    return
    end if
    allocate(idx(d-2))
    k=0
    do i=1,d
    if(i/=q .and. i/=r) then
    k=k+1
    idx(k)=i
    end if
    end do
    m=d-2
    allocate(so(m,2),os(2,m),oo(m,m),cm(m),cc(m,m),lo(m),up(m))
    so(:,1)=sigma(idx,q)
    so(:,2)=sigma(idx,r)
    os=transpose(so)
    oo=sigma(idx,idx)
    call inverse_spd(ss,invss,ok,msg)
    if(.not.ok) then
    v=merge(-huge(1.0_dp),0.0_dp,ll)
    return
    end if
    cm=mean(idx)+matmul(so,matmul(invss,xs-mus))
    cc=oo-matmul(so,matmul(invss,os))
    lo=lower(idx)
    up=upper(idx)
    if(present(control)) then
    pr=mvn_prob(lo,up,cm,cc,control)
    else
    pr=mvn_prob(lo,up,cm,cc)
    end if
    condp=max(0.0_dp,pr%value)
    if(condp<=0.0_dp) then
    v=merge(-huge(1.0_dp),0.0_dp,ll)
    return
    end if
    ld=ld+log(condp)-log(den)
    v=merge(ld,exp(ld),ll)
  end function dtmvnorm_marginal2

  real(dp) function ptmvnorm_marginal(xn,nidx,mean,sigma,lower,upper,control) result(v)
    real(dp),intent(in)::xn,mean(:),sigma(:,:),lower(:),upper(:)
    integer,intent(in)::nidx
    type(probability_control),intent(in),optional::control
    real(dp),allocatable::up2(:)
    type(probability_result)::a,b
    if(xn<=lower(nidx)) then
    v=0.0_dp
    return
    end if
    if(xn>=upper(nidx)) then
    v=1.0_dp
    return
    end if
    up2=upper
    up2(nidx)=xn
    if(present(control)) then
    a=mvn_prob(lower,up2,mean,sigma,control)
    b=mvn_prob(lower,upper,mean,sigma,control)
    else
    a=mvn_prob(lower,up2,mean,sigma)
    b=mvn_prob(lower,upper,mean,sigma)
    end if
    if(b%value<=0.0_dp) then
    v=0.0_dp
    else
    v=max(0.0_dp,min(1.0_dp,a%value/b%value))
    end if
  end function ptmvnorm_marginal

  real(dp) function ptmvt_marginal(xn,nidx,mean,sigma,df,lower,upper,control) result(v)
    real(dp),intent(in)::xn,mean(:),sigma(:,:),df,lower(:),upper(:)
    integer,intent(in)::nidx
    type(probability_control),intent(in),optional::control
    real(dp),allocatable::up2(:)
    type(probability_result)::a,b
    if(xn<=lower(nidx)) then
    v=0.0_dp
    return
    end if
    if(xn>=upper(nidx)) then
    v=1.0_dp
    return
    end if
    up2=upper
    up2(nidx)=xn
    if(present(control)) then
    a=mvt_prob(lower,up2,mean,sigma,df,control)
    b=mvt_prob(lower,upper,mean,sigma,df,control)
    else
    a=mvt_prob(lower,up2,mean,sigma,df)
    b=mvt_prob(lower,upper,mean,sigma,df)
    end if
    if(b%value<=0.0_dp) then
    v=0.0_dp
    else
    v=max(0.0_dp,min(1.0_dp,a%value/b%value))
    end if
  end function ptmvt_marginal

  real(dp) function qtmvnorm_marginal(p,nidx,mean,sigma,lower,upper,tail,interval,tol,maxit,control) result(q)
    real(dp),intent(in)::p,mean(:),sigma(:,:),lower(:),upper(:)
    integer,intent(in)::nidx
    integer,intent(in),optional::tail,maxit
    real(dp),intent(in),optional::interval(2),tol
    type(probability_control),intent(in),optional::control
    integer::tt,mi,it
    real(dp)::a,b,m,fa,fm,eps
    tt=tail_lower
    if(present(tail)) tt=tail
    mi=100
    if(present(maxit)) mi=maxit
    eps=1.0e-8_dp
    if(present(tol)) eps=tol
    if(present(interval)) then
    a=interval(1)
    b=interval(2)
    else
    a=merge(lower(nidx),mean(nidx)-10.0_dp*sqrt(sigma(nidx,nidx)),lower(nidx)>-huge(1.0_dp)/4.0_dp)
          b=merge(upper(nidx),mean(nidx)+10.0_dp*sqrt(sigma(nidx,nidx)),upper(nidx)< huge(1.0_dp)/4.0_dp)
          end if
    fa=qobj(a)
    do it=1,mi
      m=0.5_dp*(a+b)
      fm=qobj(m)
      if(abs(fm)<=eps .or. abs(b-a)<=eps*max(1.0_dp,abs(m))) exit
      if(fa*fm<=0.0_dp) then
      b=m
      else
      a=m
      fa=fm
      end if
    end do
    q=0.5_dp*(a+b)
  contains
    real(dp) function qobj(x) result(f)
      real(dp),intent(in)::x
      real(dp),allocatable::lo(:),up(:)
      type(probability_result)::aa,bb
      lo=lower
      up=upper
      select case(tt)
      case(tail_lower); up(nidx)=x
      case(tail_upper); lo(nidx)=x
      case default
      lo(nidx)=-abs(x)
      up(nidx)=abs(x)
      end select
      if(any(lo>=up)) then
      f=-p
      return
      end if
      if(present(control)) then
        aa=mvn_prob(max(lo,lower),min(up,upper),mean,sigma,control)
        bb=mvn_prob(lower,upper,mean,sigma,control)
      else
        aa=mvn_prob(max(lo,lower),min(up,upper),mean,sigma)
        bb=mvn_prob(lower,upper,mean,sigma)
      end if
      if(bb%value<=0.0_dp) then
      f=-p
      else
      f=aa%value/bb%value-p
      end if
    end function qobj
  end function qtmvnorm_marginal

  pure function outer(a,b) result(c)
    real(dp),intent(in)::a(:),b(:)
    real(dp)::c(size(a),size(b))
    c=spread(a,2,size(b))*spread(b,1,size(a))
  end function outer

end module tmvtnorm_marginals
