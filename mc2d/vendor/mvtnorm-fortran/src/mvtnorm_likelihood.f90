! SPDX-License-Identifier: GPL-2.0-only
module mvtnorm_likelihood
  use mvtnorm_kinds, only : dp
  use mvtnorm_types, only : likelihood_result, probability_control, conditional_result, probability_result
  use mvtnorm_distributions, only : dmvnorm_one
  use mvtnorm_probabilities, only : pmvnorm
  use mvtnorm_conditioning, only : marginal_mvnormal, conditional_mvnormal
  use mvtnorm_triangular, only : cov2chol, chol2cov, pack_lower, unpack_lower
  use mvtnorm_special, only : normal_cdf, normal_pdf
  implicit none
  private
  public :: ldmvnorm, lpmvnorm, ldpmvnorm
  public :: sldmvnorm, slpmvnorm, sldpmvnorm
  public :: lpRR, slpRR, standardize_cholesky

contains

  function ldmvnorm(obs,mean,covariance,individual) result(res)
    real(dp),intent(in)::obs(:,:),mean(:),covariance(:,:)
    logical,intent(in),optional::individual
    type(likelihood_result)::res
    integer::i,n
    logical::ind
    n=size(obs,1); ind=.true.; if(present(individual)) ind=individual
    allocate(res%loglik(n))
    if(size(obs,2)/=size(mean)) then; res%message='non-conforming dimensions'; return; end if
    do i=1,n
      res%loglik(i)=dmvnorm_one(obs(i,:),mean,covariance,.true.)
    end do
    res%total=sum(res%loglik); res%ok=.true.
    if(.not.ind) then
      res%loglik= [res%total]
    end if
  end function ldmvnorm

  function lpmvnorm(lower,upper,mean,covariance,control,individual) result(res)
    real(dp),intent(in)::lower(:,:),upper(:,:),mean(:),covariance(:,:)
    type(probability_control),intent(in),optional::control
    logical,intent(in),optional::individual
    type(likelihood_result)::res
    integer::i,n
    logical::ind
    real(dp)::p
    type(probability_result) :: pr
    n=size(lower,1); ind=.true.; if(present(individual)) ind=individual
    if(size(upper,1)/=n .or. size(lower,2)/=size(mean) .or. size(upper,2)/=size(mean)) then
      res%message='non-conforming dimensions'; return
    end if
    allocate(res%loglik(n))
    do i=1,n
      pr=pmvnorm(lower(i,:),upper(i,:),mean,covariance,control)
      p=pr%value
      res%loglik(i)=log(max(p,tiny(1.0_dp)))
    end do
    res%total=sum(res%loglik); res%ok=.true.
    if(.not.ind) res%loglik=[res%total]
  end function lpmvnorm

  function ldpmvnorm(obs,lower,upper,mean,covariance,control,individual) result(res)
    real(dp),intent(in)::obs(:,:),lower(:,:),upper(:,:),mean(:),covariance(:,:)
    type(probability_control),intent(in),optional::control
    logical,intent(in),optional::individual
    type(likelihood_result)::res
    integer,allocatable::exact_idx(:)
    type(conditional_result)::marg,cond
    integer::n,c,d,i,j
    real(dp)::p
    type(probability_result) :: pr
    logical::ind
    n=size(obs,1); c=size(obs,2); d=size(lower,2); ind=.true.; if(present(individual)) ind=individual
    if(size(lower,1)/=n .or. size(upper,1)/=n .or. size(upper,2)/=d .or. size(mean)/=c+d) then
      res%message='non-conforming mixed-data dimensions'; return
    end if
    allocate(exact_idx(c)); do j=1,c; exact_idx(j)=j; end do
    marg=marginal_mvnormal(mean,covariance,exact_idx)
    if(.not.marg%ok) then; res%message=marg%message; return; end if
    allocate(res%loglik(n))
    do i=1,n
      cond=conditional_mvnormal(mean,covariance,exact_idx,obs(i,:))
      if(.not.cond%ok) then; res%message=cond%message; return; end if
      pr=pmvnorm(lower(i,:),upper(i,:),cond%mean,cond%covariance,control)
      p=pr%value
      res%loglik(i)=dmvnorm_one(obs(i,:),marg%mean,marg%covariance,.true.)+log(max(p,tiny(1.0_dp)))
    end do
    res%total=sum(res%loglik); res%ok=.true.
    if(.not.ind) res%loglik=[res%total]
  end function ldpmvnorm

  function sldmvnorm(obs,mean,covariance,step) result(res)
    real(dp),intent(in)::obs(:,:),mean(:),covariance(:,:)
    real(dp),intent(in),optional::step
    type(likelihood_result)::res
    res=finite_score_density(obs,mean,covariance,step)
  end function sldmvnorm

  function slpmvnorm(lower,upper,mean,covariance,control,step) result(res)
    real(dp),intent(in)::lower(:,:),upper(:,:),mean(:),covariance(:,:)
    type(probability_control),intent(in),optional::control
    real(dp),intent(in),optional::step
    type(likelihood_result)::res
    res=finite_score_probability(lower,upper,mean,covariance,control,step)
  end function slpmvnorm

  function sldpmvnorm(obs,lower,upper,mean,covariance,control,step) result(res)
    real(dp),intent(in)::obs(:,:),lower(:,:),upper(:,:),mean(:),covariance(:,:)
    type(probability_control),intent(in),optional::control
    real(dp),intent(in),optional::step
    type(likelihood_result)::res
    type(likelihood_result)::base,rp,rm
    real(dp),allocatable::l(:,:),v(:),vp(:),vm(:),mp(:),mm(:),cp(:,:),cm(:,:)
    logical::ok
    character(len=256)::message
    real(dp)::h
    integer::p,k,j,npar,n
    p=size(mean); n=size(obs,1); h=sqrt(epsilon(1.0_dp)); if(present(step)) h=step
    base=ldpmvnorm(obs,lower,upper,mean,covariance,control,.true.)
    if(.not.base%ok) then; res=base; return; end if
    l=cov2chol(covariance,ok,message); v=pack_lower(l,.true.,.false.)
    npar=p+size(v); allocate(res%score(npar,n)); res%score=0.0_dp
    allocate(mp(p),mm(p)); mp=mean; mm=mean
    do k=1,p
      mp=mean; mm=mean; mp(k)=mp(k)+h; mm(k)=mm(k)-h
      rp=ldpmvnorm(obs,lower,upper,mp,covariance,control,.true.)
      rm=ldpmvnorm(obs,lower,upper,mm,covariance,control,.true.)
      res%score(k,:)=(rp%loglik-rm%loglik)/(2.0_dp*h)
    end do
    do j=1,size(v)
      vp=v; vm=v; h=sqrt(epsilon(1.0_dp))*max(1.0_dp,abs(v(j)))
      vp(j)=vp(j)+h; vm(j)=vm(j)-h
      cp=chol2cov(unpack_lower(vp,p,.true.,.false.))
      cm=chol2cov(unpack_lower(vm,p,.true.,.false.))
      rp=ldpmvnorm(obs,lower,upper,mean,cp,control,.true.)
      rm=ldpmvnorm(obs,lower,upper,mean,cm,control,.true.)
      res%score(p+j,:)=(rp%loglik-rm%loglik)/(2.0_dp*h)
    end do
    res%loglik=base%loglik; res%total=base%total; res%ok=.true.
  end function sldpmvnorm

  function finite_score_density(obs,mean,covariance,step) result(res)
    real(dp),intent(in)::obs(:,:),mean(:),covariance(:,:)
    real(dp),intent(in),optional::step
    type(likelihood_result)::res
    type(likelihood_result)::base,rp,rm
    real(dp),allocatable::l(:,:),v(:),vp(:),vm(:),mp(:),mm(:),cp(:,:),cm(:,:)
    logical::ok; character(len=256)::message
    real(dp)::h
    integer::p,k,j,npar,n
    p=size(mean); n=size(obs,1); h=sqrt(epsilon(1.0_dp)); if(present(step)) h=step
    base=ldmvnorm(obs,mean,covariance,.true.); if(.not.base%ok) then; res=base; return; end if
    l=cov2chol(covariance,ok,message); v=pack_lower(l,.true.,.false.); npar=p+size(v)
    allocate(res%score(npar,n),mp(p),mm(p)); res%score=0.0_dp
    do k=1,p
      mp=mean; mm=mean; mp(k)=mp(k)+h; mm(k)=mm(k)-h
      rp=ldmvnorm(obs,mp,covariance,.true.); rm=ldmvnorm(obs,mm,covariance,.true.)
      res%score(k,:)=(rp%loglik-rm%loglik)/(2.0_dp*h)
    end do
    do j=1,size(v)
      vp=v; vm=v; h=sqrt(epsilon(1.0_dp))*max(1.0_dp,abs(v(j)))
      vp(j)=vp(j)+h; vm(j)=vm(j)-h
      cp=chol2cov(unpack_lower(vp,p,.true.,.false.)); cm=chol2cov(unpack_lower(vm,p,.true.,.false.))
      rp=ldmvnorm(obs,mean,cp,.true.); rm=ldmvnorm(obs,mean,cm,.true.)
      res%score(p+j,:)=(rp%loglik-rm%loglik)/(2.0_dp*h)
    end do
    res%loglik=base%loglik; res%total=base%total; res%ok=.true.
  end function finite_score_density

  function finite_score_probability(lower,upper,mean,covariance,control,step) result(res)
    real(dp),intent(in)::lower(:,:),upper(:,:),mean(:),covariance(:,:)
    type(probability_control),intent(in),optional::control
    real(dp),intent(in),optional::step
    type(likelihood_result)::res
    type(likelihood_result)::base,rp,rm
    real(dp),allocatable::l(:,:),v(:),vp(:),vm(:),mp(:),mm(:),cp(:,:),cm(:,:)
    logical::ok; character(len=256)::message
    real(dp)::h
    integer::p,k,j,npar,n
    p=size(mean); n=size(lower,1); h=1.0e-5_dp; if(present(step)) h=step
    base=lpmvnorm(lower,upper,mean,covariance,control,.true.); if(.not.base%ok) then; res=base; return; end if
    l=cov2chol(covariance,ok,message); v=pack_lower(l,.true.,.false.); npar=p+size(v)
    allocate(res%score(npar,n),mp(p),mm(p)); res%score=0.0_dp
    do k=1,p
      mp=mean; mm=mean; mp(k)=mp(k)+h; mm(k)=mm(k)-h
      rp=lpmvnorm(lower,upper,mp,covariance,control,.true.); rm=lpmvnorm(lower,upper,mm,covariance,control,.true.)
      res%score(k,:)=(rp%loglik-rm%loglik)/(2.0_dp*h)
    end do
    do j=1,size(v)
      vp=v; vm=v; h=1.0e-5_dp*max(1.0_dp,abs(v(j)))
      vp(j)=vp(j)+h; vm(j)=vm(j)-h
      cp=chol2cov(unpack_lower(vp,p,.true.,.false.)); cm=chol2cov(unpack_lower(vm,p,.true.,.false.))
      rp=lpmvnorm(lower,upper,mean,cp,control,.true.); rm=lpmvnorm(lower,upper,mean,cm,control,.true.)
      res%score(p+j,:)=(rp%loglik-rm%loglik)/(2.0_dp*h)
    end do
    res%loglik=base%loglik; res%total=base%total; res%ok=.true.
  end function finite_score_probability

  real(dp) function lpRR(lower,upper,mean,b,d,z,weights,log_probability) result(v)
    real(dp),intent(in)::lower(:),upper(:),mean(:),b(:,:),d(:),z(:,:),weights(:)
    logical,intent(in),optional::log_probability
    logical::logp
    real(dp)::s,prodv,lo,up
    integer::m,j
    logp=.true.; if(present(log_probability)) logp=log_probability
    s=0.0_dp
    do m=1,size(z,2)
      prodv=1.0_dp
      do j=1,size(b,1)
        lo=(lower(j)-mean(j)-dot_product(b(j,:),z(:,m)))/sqrt(d(j))
        up=(upper(j)-mean(j)-dot_product(b(j,:),z(:,m)))/sqrt(d(j))
        prodv=prodv*max(0.0_dp,normal_cdf(up)-normal_cdf(lo))
      end do
      s=s+weights(min(m,size(weights)))*prodv
    end do
    if(logp) then; v=log(max(s,tiny(1.0_dp))); else; v=s; end if
  end function lpRR

  function slpRR(lower,upper,mean,b,d,z,weights,step) result(score)
    real(dp),intent(in)::lower(:),upper(:),mean(:),b(:,:),d(:),z(:,:),weights(:)
    real(dp),intent(in),optional::step
    real(dp),allocatable::score(:)
    real(dp),allocatable::mp(:),mm(:),bp(:,:),bm(:,:),dpv(:),dmv(:)
    real(dp)::h
    integer::j,k,idx,npar
    h=1.0e-5_dp; if(present(step)) h=step
    npar=size(mean)+size(b)+size(d); allocate(score(npar)); idx=0
    do j=1,size(mean)
      mp=mean; mm=mean; mp(j)=mp(j)+h; mm(j)=mm(j)-h; idx=idx+1
      score(idx)=(lpRR(lower,upper,mp,b,d,z,weights)-lpRR(lower,upper,mm,b,d,z,weights))/(2.0_dp*h)
    end do
    do k=1,size(b,2); do j=1,size(b,1)
      bp=b; bm=b; bp(j,k)=bp(j,k)+h; bm(j,k)=bm(j,k)-h; idx=idx+1
      score(idx)=(lpRR(lower,upper,mean,bp,d,z,weights)-lpRR(lower,upper,mean,bm,d,z,weights))/(2.0_dp*h)
    end do; end do
    do j=1,size(d)
      dpv=d; dmv=d; dpv(j)=dpv(j)+h; dmv(j)=max(tiny(1.0_dp),dmv(j)-h); idx=idx+1
      score(idx)=(lpRR(lower,upper,mean,b,dpv,z,weights)-lpRR(lower,upper,mean,b,dmv,z,weights))/(dpv(j)-dmv(j))
    end do
  end function slpRR

  function standardize_cholesky(chol) result(out)
    real(dp),intent(in)::chol(:,:)
    real(dp),allocatable::out(:,:)
    integer::i
    out=chol
    do i=1,size(chol,1)
      out(i,:)=out(i,:)/sqrt(dot_product(chol(i,:),chol(i,:)))
    end do
  end function standardize_cholesky

end module mvtnorm_likelihood
