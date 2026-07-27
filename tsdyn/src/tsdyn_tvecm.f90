! SPDX-License-Identifier: GPL-2.0-or-later
! Numerical translation derived from tsDyn, copyright its original authors.
! This file may be redistributed and/or modified under GPL version 2 or later.
module tsdyn_tvecm
  use tsdyn_kinds, only: dp, n_deterministic, include_none, include_const, include_trend, include_both
  use tsdyn_linalg, only: ols_fit, covariance_matrix, inverse_matrix
  use tsdyn_utils, only: build_deterministic, percentile_grid, random_normal_vector
  use tsdyn_vecm, only: vecm_model, fit_vecm
  implicit none
  private
  public :: tvecm_model, fit_tvecm, forecast_tvecm, simulate_tvecm, tvecm_regimes, tvecm_lr_statistic

  type :: tvecm_model
    integer :: nvar=0,lag_diff=0,nthresh=0,nregime=0,include=include_const,nobs=0,rank=0
    character(len=8) :: common='All'
    real(dp),allocatable :: beta(:),thresholds(:),ect_coefficients(:,:),short_coefficients(:,:,:)
    real(dp),allocatable :: fitted(:,:),residuals(:,:),sigma(:,:),ect(:)
    integer,allocatable :: regimes(:),regime_counts(:)
    real(dp)::ssr=0.0_dp,loglik=0.0_dp,aic=0.0_dp,bic=0.0_dp
  end type tvecm_model
contains
  pure integer function classify(z,thr) result(r)
    real(dp),intent(in)::z,thr(:)
    integer::j
    r=size(thr)+1
    do j=1,size(thr);if(z<=thr(j))then;r=j;return;end if;end do
  end function classify

  subroutine build_base(y,p,include,beta,x1,short,target,ect,info)
    real(dp),intent(in)::y(:,:),beta(:)
    integer,intent(in)::p,include
    real(dp),allocatable,intent(out)::x1(:,:),short(:,:),target(:,:),ect(:)
    integer,intent(out)::info
    real(dp),allocatable::det(:,:)
    integer::n,k,nobs,nd,i,j,t
    n=size(y,1);k=size(y,2);nobs=n-p-1;nd=n_deterministic(include)
    if(size(beta)/=k.or.nobs<2)then;info=-1;allocate(x1(0,0),short(0,0),target(0,0),ect(0));return;end if
    allocate(x1(nobs,k),short(nobs,nd+k*p),target(nobs,k),ect(nobs));call build_deterministic(nobs,include,det,trend_start=real(p+2,dp))
    if(nd>0)short(:,1:nd)=det
    do i=1,nobs
      t=p+1+i;x1(i,:)=y(t-1,:);target(i,:)=y(t,:)-y(t-1,:);ect(i)=dot_product(x1(i,:),beta)
      do j=1,p;short(i,nd+(j-1)*k+1:nd+j*k)=y(t-j,:)-y(t-j-1,:);end do
    end do
    info=0
  end subroutine build_base

  subroutine build_design(y,p,include,beta,thresholds,common,design,target,ect,reg,counts,info)
    real(dp),intent(in)::y(:,:),beta(:),thresholds(:)
    integer,intent(in)::p,include
    character(len=*),intent(in)::common
    real(dp),allocatable,intent(out)::design(:,:),target(:,:),ect(:)
    integer,allocatable,intent(out)::reg(:),counts(:)
    integer,intent(out)::info
    real(dp),allocatable::x1(:,:),short(:,:)
    integer::nobs,nreg,nshort,i,r,off
    call build_base(y,p,include,beta,x1,short,target,ect,info);if(info/=0)then;allocate(design(0,0),reg(0),counts(0));return;end if
    nobs=size(target,1);nreg=size(thresholds)+1;nshort=size(short,2);allocate(reg(nobs),counts(nreg));counts=0
    if(trim(common)=='only_ECT'.or.trim(common)=='only_ect')then
      allocate(design(nobs,nreg+nshort));design=0.0_dp
      do i=1,nobs;r=classify(ect(i),thresholds);reg(i)=r;counts(r)=counts(r)+1;design(i,r)=ect(i);if(nshort>0)design(i,nreg+1:)=short(i,:);end do
    else
      allocate(design(nobs,nreg*(1+nshort)));design=0.0_dp
      do i=1,nobs;r=classify(ect(i),thresholds);reg(i)=r;counts(r)=counts(r)+1;off=(r-1)*(1+nshort);design(i,off+1)=ect(i);if(nshort>0)design(i,off+2:off+1+nshort)=short(i,:);end do
    end if
    info=0
  end subroutine build_design

  subroutine evaluate(y,p,include,beta,thresholds,common,trimv,ssr,bcoef,fit,res,ect,reg,counts,rank,info)
    real(dp),intent(in)::y(:,:),beta(:),thresholds(:),trimv
    integer,intent(in)::p,include
    character(len=*),intent(in)::common
    real(dp),intent(out)::ssr
    real(dp),allocatable,intent(out)::bcoef(:,:),fit(:,:),res(:,:),ect(:)
    integer,allocatable,intent(out)::reg(:),counts(:)
    integer,intent(out)::rank,info
    real(dp),allocatable::design(:,:),target(:,:)
    integer::r,minpar,nobs
    call build_design(y,p,include,beta,thresholds,common,design,target,ect,reg,counts,info)
    if(info/=0)then;ssr=huge(1.0_dp);allocate(bcoef(0,0),fit(0,0),res(0,0));rank=0;return;end if
    nobs=size(target,1);minpar=1+merge(0,size(target,2)*p+n_deterministic(include),trim(common)=='only_ECT'.or.trim(common)=='only_ect')
    do r=1,size(counts)
      if(real(counts(r),dp)<trimv*real(nobs,dp).or.counts(r)<=minpar)then;info=1;ssr=huge(1.0_dp);allocate(bcoef(0,0),fit(0,0),res(0,0));rank=0;return;end if
    end do
    call ols_fit(design,target,bcoef,fit,res,rank,ssr,info)
  end subroutine evaluate

  subroutine fit_tvecm(y,p,nthresh,include,model,info,beta_fixed,thresholds,common,trim_fraction,ngrid_beta,ngrid_threshold)
    real(dp),intent(in)::y(:,:)
    integer,intent(in)::p,nthresh,include
    type(tvecm_model),intent(out)::model
    integer,intent(out)::info
    real(dp),intent(in),optional::beta_fixed(:),thresholds(:),trim_fraction
    character(len=*),intent(in),optional::common
    integer,intent(in),optional::ngrid_beta,ngrid_threshold
    integer::k,nreg,nb,nt,ib,i,j,istat,rank,r,off,nshort
    real(dp)::tr,best,ss,det_sigma
    character(len=8)::com
    real(dp),allocatable::betas(:,:),beta(:),thr(:),grid(:),bestbeta(:),bestthr(:),bcoef(:,:),fit(:,:),res(:,:),ect(:),bb(:,:),ff(:,:),rr(:,:),ee(:)
    integer,allocatable::reg(:),counts(:),breg(:),bcounts(:)
    k=size(y,2);nreg=nthresh+1;tr=0.05_dp;if(present(trim_fraction))tr=trim_fraction;nb=21;if(present(ngrid_beta))nb=max(3,ngrid_beta);nt=30;if(present(ngrid_threshold))nt=max(3,ngrid_threshold)
    allocate(betas(k,0))
    com='All';if(present(common))com=adjustl(common)
    if(nthresh<1.or.nthresh>2)then;info=-1;return;end if
    if(present(beta_fixed))then
      if(size(beta_fixed)/=k)then;info=-2;return;end if
      deallocate(betas);allocate(betas(k,1));betas(:,1)=beta_fixed
    else
      if(k/=2)then;info=-3;return;end if
      deallocate(betas);call beta_grid_2d(y,nb,betas,istat);if(istat/=0)then;info=istat;return;end if
    end if
    allocate(bestbeta(k),bestthr(nthresh),beta(k),thr(nthresh));best=huge(1.0_dp);bestbeta=0.0_dp;bestthr=0.0_dp
    do ib=1,size(betas,2)
      beta=betas(:,ib)
      if(present(thresholds))then
        if(size(thresholds)/=nthresh)then;info=-4;return;end if
        call evaluate(y,p,include,beta,thresholds,com,tr,ss,bcoef,fit,res,ect,reg,counts,rank,istat)
        if(istat==0.and.ss<best)then;best=ss;bestbeta=beta;bestthr=thresholds;call move_all(bcoef,fit,res,ect,reg,counts,bb,ff,rr,ee,breg,bcounts);end if
      else
        call percentile_grid_for_beta(y,p,beta,tr,nt,grid,istat);if(istat/=0)cycle
        if(nthresh==1)then
          do i=1,size(grid);thr(1)=grid(i);call try_candidate(beta,thr);end do
        else
          do i=1,size(grid)-1;do j=i+1,size(grid);thr=[grid(i),grid(j)];call try_candidate(beta,thr);end do;end do
        end if
      end if
    end do
    if(best>=0.5_dp*huge(1.0_dp))then;info=2;return;end if
    nshort=n_deterministic(include)+k*p
    model%nvar=k;model%lag_diff=p;model%nthresh=nthresh;model%nregime=nreg;model%include=include;model%common=com;model%nobs=size(ff,1);model%rank=rank;model%ssr=best
    allocate(model%beta(k),model%thresholds(nthresh),model%ect_coefficients(nreg,k),model%short_coefficients(nshort,nreg,k));model%beta=bestbeta;model%thresholds=bestthr;model%ect_coefficients=0.0_dp;model%short_coefficients=0.0_dp
    if(trim(com)=='only_ECT'.or.trim(com)=='only_ect')then
      do r=1,nreg;model%ect_coefficients(r,:)=bb(r,:);if(nshort>0)model%short_coefficients(:,r,:)=bb(nreg+1:,:);end do
    else
      do r=1,nreg;off=(r-1)*(1+nshort);model%ect_coefficients(r,:)=bb(off+1,:);if(nshort>0)model%short_coefficients(:,r,:)=bb(off+2:off+1+nshort,:);end do
    end if
    allocate(model%fitted(size(ff,1),k),model%residuals(size(rr,1),k),model%ect(size(ee)),model%regimes(size(breg)),model%regime_counts(size(bcounts)))
    model%fitted=ff;model%residuals=rr;model%ect=ee;model%regimes=breg;model%regime_counts=bcounts;call covariance_matrix(rr,model%sigma,center=.false.)
    det_sigma=determinant_spd(model%sigma);model%loglik=-0.5_dp*real(model%nobs,dp)*(real(k,dp)*(log(2.0_dp*acos(-1.0_dp))+1.0_dp)+log(max(det_sigma,tiny(1.0_dp))))
    model%aic=-2.0_dp*model%loglik+2.0_dp*real(rank+nthresh+k,dp);model%bic=-2.0_dp*model%loglik+log(real(model%nobs,dp))*real(rank+nthresh+k,dp);info=0
  contains
    subroutine try_candidate(beta,thr)
      real(dp),intent(in)::beta(:),thr(:)
      call evaluate(y,p,include,beta,thr,com,tr,ss,bcoef,fit,res,ect,reg,counts,rank,istat)
      if(istat==0.and.ss<best)then;best=ss;bestbeta=beta;bestthr=thr;call move_all(bcoef,fit,res,ect,reg,counts,bb,ff,rr,ee,breg,bcounts);end if
    end subroutine try_candidate
  end subroutine fit_tvecm

  subroutine beta_grid_2d(y,nb,betas,info)
    real(dp),intent(in)::y(:,:)
    integer,intent(in)::nb
    real(dp),allocatable,intent(out)::betas(:,:)
    integer,intent(out)::info
    real(dp),allocatable::x(:,:),ym(:,:),b(:,:),fit(:,:),res(:,:),xtx(:,:),inv(:,:)
    real(dp)::ssr,se,center
    integer::rank,i,istat
    allocate(x(size(y,1),1),ym(size(y,1),1));x(:,1)=y(:,2);ym(:,1)=y(:,1)
    call ols_fit(x,ym,b,fit,res,rank,ssr,istat);if(istat/=0)then;info=istat;return;end if
    center=b(1,1);allocate(xtx(1,1));xtx=matmul(transpose(x),x);call inverse_matrix(xtx,inv,istat);if(istat/=0)then;info=istat;return;end if
    se=sqrt(max(ssr/real(max(1,size(y,1)-1),dp)*inv(1,1),tiny(1.0_dp)))
    allocate(betas(2,nb));do i=1,nb;betas(1,i)=1.0_dp;betas(2,i)=-(center-2.0_dp*se+4.0_dp*se*real(i-1,dp)/real(max(1,nb-1),dp));end do;info=0
  end subroutine beta_grid_2d

  subroutine percentile_grid_for_beta(y,p,beta,trimv,ngrid,grid,info)
    real(dp),intent(in)::y(:,:),beta(:),trimv
    integer,intent(in)::p,ngrid
    real(dp),allocatable,intent(out)::grid(:)
    integer,intent(out)::info
    real(dp),allocatable::ect(:)
    integer::i,t,nobs
    nobs=size(y,1)-p-1;if(nobs<1)then;info=-1;allocate(grid(0));return;end if
    allocate(ect(nobs));do i=1,nobs;t=p+1+i;ect(i)=dot_product(y(t-1,:),beta);end do
    call percentile_grid(ect,trimv,ngrid,grid);info=0
  end subroutine percentile_grid_for_beta

  subroutine move_all(bcoef,fit,res,ect,reg,counts,bb,ff,rr,ee,breg,bcounts)
    real(dp),allocatable,intent(inout)::bcoef(:,:),fit(:,:),res(:,:),ect(:),bb(:,:),ff(:,:),rr(:,:),ee(:)
    integer,allocatable,intent(inout)::reg(:),counts(:),breg(:),bcounts(:)
    if(allocated(bb))deallocate(bb);if(allocated(ff))deallocate(ff);if(allocated(rr))deallocate(rr);if(allocated(ee))deallocate(ee);if(allocated(breg))deallocate(breg);if(allocated(bcounts))deallocate(bcounts)
    call move_alloc(bcoef,bb);call move_alloc(fit,ff);call move_alloc(res,rr);call move_alloc(ect,ee);call move_alloc(reg,breg);call move_alloc(counts,bcounts)
  end subroutine move_all

  real(dp) function determinant_spd(a) result(det)
    use tsdyn_linalg, only: cholesky_lower
    real(dp),intent(in)::a(:,:)
    real(dp),allocatable::l(:,:)
    integer::i,istat
    call cholesky_lower(a,l,istat);if(istat/=0)then;det=tiny(1.0_dp);return;end if;det=1.0_dp;do i=1,size(l,1);det=det*l(i,i)**2;end do
  end function determinant_spd

  pure function tvecm_delta(model,work,t,r) result(dy)
    type(tvecm_model),intent(in)::model
    real(dp),intent(in)::work(:,:)
    integer,intent(in)::t,r
    real(dp)::dy(model%nvar),ect
    real(dp),allocatable::short(:)
    integer::nd,j,k
    k=model%nvar;nd=n_deterministic(model%include);allocate(short(nd+k*model%lag_diff));short=0.0_dp
    ect=dot_product(work(t-1,:),model%beta)
    select case(model%include)
    case(include_const);short(1)=1.0_dp
    case(include_trend);short(1)=real(t,dp)
    case(include_both);short(1)=1.0_dp;short(2)=real(t,dp)
    end select
    do j=1,model%lag_diff;short(nd+(j-1)*k+1:nd+j*k)=work(t-j,:)-work(t-j-1,:);end do
    dy=ect*model%ect_coefficients(r,:)
    if(size(short)>0)dy=dy+matmul(short,model%short_coefficients(:,r,:))
  end function tvecm_delta

  subroutine forecast_tvecm(model,history,h,forecast,regimes,info)
    type(tvecm_model),intent(in)::model
    real(dp),intent(in)::history(:,:)
    integer,intent(in)::h
    real(dp),allocatable,intent(out)::forecast(:,:)
    integer,allocatable,intent(out),optional::regimes(:)
    integer,intent(out)::info
    real(dp),allocatable::work(:,:)
    real(dp)::ect
    integer::t,idx,r
    if(size(history,1)<model%lag_diff+1.or.size(history,2)/=model%nvar.or.h<1)then;info=-1;allocate(forecast(0,0));if(present(regimes))allocate(regimes(0));return;end if
    allocate(work(size(history,1)+h,model%nvar),forecast(h,model%nvar));work(1:size(history,1),:)=history;if(present(regimes))allocate(regimes(h))
    do t=1,h;idx=size(history,1)+t;ect=dot_product(work(idx-1,:),model%beta);r=classify(ect,model%thresholds);work(idx,:)=work(idx-1,:)+tvecm_delta(model,work,idx,r);forecast(t,:)=work(idx,:);if(present(regimes))regimes(t)=r;end do
    info=0
  end subroutine forecast_tvecm

  subroutine simulate_tvecm(model,n,y,info,innov,start)
    type(tvecm_model),intent(in)::model
    integer,intent(in)::n
    real(dp),allocatable,intent(out)::y(:,:)
    integer,intent(out)::info
    real(dp),intent(in),optional::innov(:,:),start(:,:)
    real(dp),allocatable::work(:,:),e(:)
    real(dp)::ect
    integer::p,t,r
    p=model%lag_diff;if(n<1)then;info=-1;allocate(y(0,0));return;end if
    allocate(work(n+p+1,model%nvar),e(model%nvar));work=0.0_dp
    if(present(start))then;if(size(start,1)<p+1.or.size(start,2)/=model%nvar)then;info=-2;allocate(y(0,0));return;end if;work(1:p+1,:)=start(size(start,1)-p:,:);end if
    do t=p+2,n+p+1;ect=dot_product(work(t-1,:),model%beta);r=classify(ect,model%thresholds);if(present(innov))then;e=innov(t-p-1,:);else;call random_normal_vector(e);end if;work(t,:)=work(t-1,:)+tvecm_delta(model,work,t,r)+e;end do
    allocate(y(n,model%nvar));y=work(p+2:,:);info=0
  end subroutine simulate_tvecm

  subroutine tvecm_regimes(model,y,regimes,info)
    type(tvecm_model),intent(in)::model
    real(dp),intent(in)::y(:,:)
    integer,allocatable,intent(out)::regimes(:)
    integer,intent(out)::info
    integer::i,t
    if(size(y,1)<=model%lag_diff+1)then;info=-1;allocate(regimes(0));return;end if
    allocate(regimes(size(y,1)-model%lag_diff-1));do i=1,size(regimes);t=model%lag_diff+1+i;regimes(i)=classify(dot_product(y(t-1,:),model%beta),model%thresholds);end do;info=0
  end subroutine tvecm_regimes

  subroutine tvecm_lr_statistic(y,p,include,beta,statistic,threshold,info,ngrid)
    real(dp),intent(in)::y(:,:),beta(:)
    integer,intent(in)::p,include
    real(dp),intent(out)::statistic,threshold
    integer,intent(out)::info
    integer,intent(in),optional::ngrid
    type(vecm_model)::lin
    type(tvecm_model)::tv
    real(dp),allocatable::bm(:,:)
    integer::istat,ng
    ng=30;if(present(ngrid))ng=ngrid;allocate(bm(size(beta),1));bm(:,1)=beta
    call fit_vecm(y,p,1,include,'fixed',lin,istat,beta_fixed=bm);if(istat/=0)then;info=istat;return;end if
    call fit_tvecm(y,p,1,include,tv,istat,beta_fixed=beta,ngrid_threshold=ng);if(istat/=0)then;info=istat;return;end if
    statistic=real(tv%nobs,dp)*log(lin%ssr/tv%ssr);threshold=tv%thresholds(1);info=0
  end subroutine tvecm_lr_statistic
end module tsdyn_tvecm
