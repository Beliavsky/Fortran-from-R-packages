! SPDX-License-Identifier: GPL-2.0-or-later
module flexsurv_standardize_advanced
  use flexsurv_kinds, only : dp
  use flexsurv_fit, only : flexsurv_spec, predict_survival, predict_hazard
  use flexsurv_standardize, only : standsurv_survival
  use flexsurv_math, only : integrate_gauss_legendre, normal_quantile
  use relsurv_ratetable, only : ratetable_type, expected_survival, pystep2
  use numderiv, only : grad
  implicit none
  private

  integer,parameter,public::stand_transform_none=0,stand_transform_log=1
  integer,parameter,public::stand_transform_loglog=2,stand_transform_logit=3
  integer,parameter,public::stand_contrast_difference=1,stand_contrast_ratio=2

  type,public::standsurv_ci
    real(dp)::estimate=0.0_dp
    real(dp)::se=0.0_dp
    real(dp)::lower=0.0_dp
    real(dp)::upper=0.0_dp
  end type standsurv_ci

  public::standsurv_acsurvival,standsurv_achazard,standsurv_acrmst,standsurv_acquantile
  public::standsurv_delta_survival,standsurv_delta_acsurvival,standsurv_delta_contrast
  public::standsurv_bootstrap_survival,standsurv_bootstrap_acsurvival
  public::population_hazard_at,stand_transform,stand_inv_transform

contains

  real(dp) function standsurv_acsurvival(spec,theta,t,tab,rate_x,weights,scale_ratetable) result(s)
    type(flexsurv_spec),intent(in)::spec
    real(dp),intent(in)::theta(:),t
    type(ratetable_type),intent(in)::tab
    real(dp),intent(in)::rate_x(:,:)
    real(dp),intent(in),optional::weights(:),scale_ratetable
    real(dp),allocatable::tv(:),es(:)
    real(dp)::w,sw,sc
    integer::i,n
    n=size(rate_x,1);sc=1.0_dp;if(present(scale_ratetable))sc=scale_ratetable
    allocate(tv(n),es(n));tv=max(0.0_dp,t*sc);es=expected_survival(tab,rate_x,tv)
    s=0.0_dp;sw=0.0_dp
    do i=1,n
      w=1.0_dp;if(present(weights))w=weights(i)
      s=s+w*predict_survival(spec,theta,i,t)*es(i);sw=sw+w
    end do
    if(sw>0.0_dp)s=s/sw
  end function standsurv_acsurvival

  real(dp) function standsurv_achazard(spec,theta,t,tab,rate_x,weights,scale_ratetable, &
      reproduce_upstream_weighted_bug) result(h)
    type(flexsurv_spec),intent(in)::spec
    real(dp),intent(in)::theta(:),t
    type(ratetable_type),intent(in)::tab
    real(dp),intent(in)::rate_x(:,:)
    real(dp),intent(in),optional::weights(:),scale_ratetable
    logical,intent(in),optional::reproduce_upstream_weighted_bug
    real(dp),allocatable::tv(:),es(:),eh(:)
    real(dp)::w,den,rs,exh,sc,totalh
    logical::bug
    integer::i,n
    n=size(rate_x,1);sc=1.0_dp;if(present(scale_ratetable))sc=scale_ratetable
    allocate(tv(n),es(n),eh(n));tv=max(0.0_dp,t*sc);es=expected_survival(tab,rate_x,tv)
    call population_hazard_at(tab,rate_x,t*sc,eh);eh=eh*sc
    bug=.true.;if(present(reproduce_upstream_weighted_bug))bug=reproduce_upstream_weighted_bug
    h=0.0_dp;den=0.0_dp
    do i=1,n
      w=1.0_dp;if(present(weights))w=weights(i)
      rs=predict_survival(spec,theta,i,t);exh=predict_hazard(spec,theta,i,t)
      totalh=exh+eh(i)
      if(present(weights).and.bug)totalh=exh*eh(i)
      h=h+w*rs*es(i)*totalh;den=den+w*rs*es(i)
    end do
    if(den>0.0_dp)h=h/den
  end function standsurv_achazard

  real(dp) function standsurv_acrmst(spec,theta,t,tab,rate_x,weights,scale_ratetable,nquad) result(r)
    type(flexsurv_spec),intent(in)::spec
    real(dp),intent(in)::theta(:),t
    type(ratetable_type),intent(in)::tab
    real(dp),intent(in)::rate_x(:,:)
    real(dp),intent(in),optional::weights(:),scale_ratetable
    integer,intent(in),optional::nquad
    integer::nq
    if(t<=0.0_dp)then;r=0.0_dp;return;end if
    nq=100;if(present(nquad))nq=nquad
    r=integrate_gauss_legendre(fn,0.0_dp,t,nq)
  contains
    real(dp) function fn(x) result(v)
      real(dp),intent(in)::x
      v=standsurv_acsurvival(spec,theta,x,tab,rate_x,weights,scale_ratetable)
    end function fn
  end function standsurv_acrmst

  real(dp) function standsurv_acquantile(spec,theta,p,tab,rate_x,weights,scale_ratetable,interval) result(q)
    type(flexsurv_spec),intent(in)::spec
    real(dp),intent(in)::theta(:),p
    type(ratetable_type),intent(in)::tab
    real(dp),intent(in)::rate_x(:,:)
    real(dp),intent(in),optional::weights(:),scale_ratetable,interval(:)
    real(dp)::lo,hi,mid,target
    integer::it
    if(p<=0.0_dp)then;q=0.0_dp;return;end if
    target=1.0_dp-p;lo=0.0_dp;hi=1.0_dp
    if(present(interval))then
      if(size(interval)>=2)then;lo=max(0.0_dp,interval(1));hi=max(lo,interval(2));end if
    else
      do while(standsurv_acsurvival(spec,theta,hi,tab,rate_x,weights,scale_ratetable)>target.and.hi<1.0e12_dp)
        hi=2.0_dp*hi
      end do
    end if
    do it=1,120
      mid=0.5_dp*(lo+hi)
      if(standsurv_acsurvival(spec,theta,mid,tab,rate_x,weights,scale_ratetable)>target)then
        lo=mid
      else
        hi=mid
      end if
    end do
    q=0.5_dp*(lo+hi)
  end function standsurv_acquantile

  subroutine population_hazard_at(tab,x,time,haz)
    type(ratetable_type),intent(in)::tab
    real(dp),intent(in)::x(:,:),time
    real(dp),intent(out)::haz(size(x,1))
    real(dp)::d(tab%ndim),wt
    integer::i,j,idx,idx2
    do i=1,size(x,1)
      d=x(i,:)
      do j=1,tab%ndim
        if(tab%factor(j)/=1)d(j)=d(j)+time
      end do
      call pystep2(tab,d,idx,idx2,wt)
      haz(i)=wt*tab%rate(idx)+(1.0_dp-wt)*tab%rate(idx2)
    end do
  end subroutine population_hazard_at

  function standsurv_delta_survival(spec,theta,cov,t,weights,transform,cl) result(out)
    type(flexsurv_spec),intent(in)::spec
    real(dp),intent(in)::theta(:),cov(:,:),t
    real(dp),intent(in),optional::weights(:),cl
    integer,intent(in),optional::transform
    type(standsurv_ci)::out
    real(dp),allocatable::g(:)
    real(dp)::v,z,lev,estt
    integer::trn
    trn=stand_transform_none;if(present(transform))trn=transform
    lev=0.95_dp;if(present(cl))lev=cl
    allocate(g(size(theta)));call grad(fn,theta,g)
    v=max(0.0_dp,dot_product(g,matmul(cov,g)));out%estimate=standsurv_survival(spec,theta,t,weights)
    out%se=sqrt(v);z=normal_quantile(0.5_dp+0.5_dp*lev);estt=stand_transform(out%estimate,trn)
    out%lower=stand_inv_transform(estt-z*out%se,trn);out%upper=stand_inv_transform(estt+z*out%se,trn)
  contains
    real(dp) function fn(x) result(vv)
      real(dp),intent(in)::x(:)
      vv=stand_transform(standsurv_survival(spec,x,t,weights),trn)
    end function fn
  end function standsurv_delta_survival

  function standsurv_delta_acsurvival(spec,theta,cov,t,tab,rate_x,weights,scale_ratetable,transform,cl) result(out)
    type(flexsurv_spec),intent(in)::spec
    real(dp),intent(in)::theta(:),cov(:,:),t
    type(ratetable_type),intent(in)::tab
    real(dp),intent(in)::rate_x(:,:)
    real(dp),intent(in),optional::weights(:),scale_ratetable,cl
    integer,intent(in),optional::transform
    type(standsurv_ci)::out
    real(dp),allocatable::g(:)
    real(dp)::v,z,lev,estt
    integer::trn
    trn=stand_transform_none;if(present(transform))trn=transform
    lev=0.95_dp;if(present(cl))lev=cl
    allocate(g(size(theta)));call grad(fn,theta,g)
    v=max(0.0_dp,dot_product(g,matmul(cov,g)))
    out%estimate=standsurv_acsurvival(spec,theta,t,tab,rate_x,weights,scale_ratetable)
    out%se=sqrt(v);z=normal_quantile(0.5_dp+0.5_dp*lev);estt=stand_transform(out%estimate,trn)
    out%lower=stand_inv_transform(estt-z*out%se,trn);out%upper=stand_inv_transform(estt+z*out%se,trn)
  contains
    real(dp) function fn(x) result(vv)
      real(dp),intent(in)::x(:)
      vv=stand_transform(standsurv_acsurvival(spec,x,t,tab,rate_x,weights,scale_ratetable),trn)
    end function fn
  end function standsurv_delta_acsurvival

  function standsurv_delta_contrast(spec1,spec0,theta,cov,t,kind,weights1,weights0,transform,cl) result(out)
    type(flexsurv_spec),intent(in)::spec1,spec0
    real(dp),intent(in)::theta(:),cov(:,:),t
    integer,intent(in)::kind
    real(dp),intent(in),optional::weights1(:),weights0(:),cl
    integer,intent(in),optional::transform
    type(standsurv_ci)::out
    real(dp),allocatable::g(:)
    real(dp)::a,b,v,z,lev,estt
    integer::trn
    trn=stand_transform_none;if(present(transform))trn=transform
    lev=0.95_dp;if(present(cl))lev=cl
    a=standsurv_survival(spec1,theta,t,weights1);b=standsurv_survival(spec0,theta,t,weights0)
    out%estimate=contrast_value(a,b,kind);allocate(g(size(theta)));call grad(fn,theta,g)
    v=max(0.0_dp,dot_product(g,matmul(cov,g)));out%se=sqrt(v);z=normal_quantile(0.5_dp+0.5_dp*lev)
    estt=stand_transform(out%estimate,trn);out%lower=stand_inv_transform(estt-z*out%se,trn)
    out%upper=stand_inv_transform(estt+z*out%se,trn)
  contains
    real(dp) function fn(x) result(vv)
      real(dp),intent(in)::x(:)
      real(dp)::aa,bb
      aa=standsurv_survival(spec1,x,t,weights1);bb=standsurv_survival(spec0,x,t,weights0)
      vv=stand_transform(contrast_value(aa,bb,kind),trn)
    end function fn
  end function standsurv_delta_contrast

  subroutine standsurv_bootstrap_survival(spec,draws,t,est,se,lower,upper,weights,cl)
    type(flexsurv_spec),intent(in)::spec
    real(dp),intent(in)::draws(:,:),t(:)
    real(dp),intent(out)::est(size(t)),se(size(t)),lower(size(t)),upper(size(t))
    real(dp),intent(in),optional::weights(:),cl
    real(dp),allocatable::sim(:,:)
    real(dp)::lev
    integer::i,b
    lev=0.95_dp;if(present(cl))lev=cl;allocate(sim(size(draws,2),size(t)))
    do b=1,size(draws,2);do i=1,size(t);sim(b,i)=standsurv_survival(spec,draws(:,b),t(i),weights);end do;end do
    call summarize_boot(sim,lev,est,se,lower,upper)
  end subroutine standsurv_bootstrap_survival

  subroutine standsurv_bootstrap_acsurvival(spec,draws,t,tab,rate_x,est,se,lower,upper,weights,scale_ratetable,cl)
    type(flexsurv_spec),intent(in)::spec
    real(dp),intent(in)::draws(:,:),t(:)
    type(ratetable_type),intent(in)::tab
    real(dp),intent(in)::rate_x(:,:)
    real(dp),intent(out)::est(size(t)),se(size(t)),lower(size(t)),upper(size(t))
    real(dp),intent(in),optional::weights(:),scale_ratetable,cl
    real(dp),allocatable::sim(:,:)
    real(dp)::lev
    integer::i,b
    lev=0.95_dp;if(present(cl))lev=cl;allocate(sim(size(draws,2),size(t)))
    do b=1,size(draws,2);do i=1,size(t)
      sim(b,i)=standsurv_acsurvival(spec,draws(:,b),t(i),tab,rate_x,weights,scale_ratetable)
    end do;end do
    call summarize_boot(sim,lev,est,se,lower,upper)
  end subroutine standsurv_bootstrap_acsurvival

  pure real(dp) function stand_transform(x,kind) result(y)
    real(dp),intent(in)::x
    integer,intent(in)::kind
    select case(kind)
    case(stand_transform_log);y=log(max(x,tiny(1.0_dp)))
    case(stand_transform_loglog);y=log(max(-log(max(1.0_dp-x,tiny(1.0_dp))),tiny(1.0_dp)))
    case(stand_transform_logit);y=log(max(x,tiny(1.0_dp))/max(1.0_dp-x,tiny(1.0_dp)))
    case default;y=x
    end select
  end function stand_transform

  pure real(dp) function stand_inv_transform(x,kind) result(y)
    real(dp),intent(in)::x
    integer,intent(in)::kind
    select case(kind)
    case(stand_transform_log);y=exp(x)
    case(stand_transform_loglog);y=1.0_dp-exp(-exp(x))
    case(stand_transform_logit);if(x>=0.0_dp)then;y=1.0_dp/(1.0_dp+exp(-x));else;y=exp(x)/(1.0_dp+exp(x));end if
    case default;y=x
    end select
  end function stand_inv_transform

  pure real(dp) function contrast_value(a,b,kind) result(v)
    real(dp),intent(in)::a,b
    integer,intent(in)::kind
    if(kind==stand_contrast_ratio)then
      if(abs(b)>tiny(1.0_dp))then;v=a/b;else;v=huge(1.0_dp);end if
    else
      v=a-b
    end if
  end function contrast_value

  subroutine summarize_boot(sim,cl,est,se,lower,upper)
    real(dp),intent(in)::sim(:,:),cl
    real(dp),intent(out)::est(size(sim,2)),se(size(sim,2)),lower(size(sim,2)),upper(size(sim,2))
    real(dp),allocatable::x(:)
    real(dp)::mu
    integer::j
    allocate(x(size(sim,1)))
    do j=1,size(sim,2)
      x=sim(:,j);mu=sum(x)/real(size(x),dp);est(j)=mu
      if(size(x)>1)then;se(j)=sqrt(sum((x-mu)**2)/real(size(x)-1,dp));else;se(j)=0.0_dp;end if
      call quant_pair(x,(1.0_dp-cl)/2.0_dp,lower(j),upper(j))
    end do
  end subroutine summarize_boot

  subroutine quant_pair(x,p,lo,hi)
    real(dp),intent(in)::x(:),p
    real(dp),intent(out)::lo,hi
    real(dp),allocatable::y(:)
    integer::n,il,iu
    allocate(y(size(x)));y=x;call sort_real(y);n=size(y)
    il=max(1,min(n,nint(1.0_dp+p*real(n-1,dp))))
    iu=max(1,min(n,nint(1.0_dp+(1.0_dp-p)*real(n-1,dp))))
    lo=y(il);hi=y(iu)
  end subroutine quant_pair

  subroutine sort_real(x)
    real(dp),intent(inout)::x(:)
    integer::i,j
    real(dp)::key
    do i=2,size(x)
      key=x(i);j=i-1
      do while(j>=1)
        if(x(j)<=key)exit
        x(j+1)=x(j);j=j-1
      end do
      x(j+1)=key
    end do
  end subroutine sort_real

end module flexsurv_standardize_advanced
