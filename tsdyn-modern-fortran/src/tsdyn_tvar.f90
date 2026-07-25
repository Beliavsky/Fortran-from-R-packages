! SPDX-License-Identifier: GPL-2.0-or-later
! Numerical translation derived from tsDyn, copyright its original authors.
! This file may be redistributed and/or modified under GPL version 2 or later.
module tsdyn_tvar
  use tsdyn_kinds, only: dp, n_deterministic, include_none, include_const, include_trend, include_both
  use tsdyn_linalg, only: ols_fit, covariance_matrix
  use tsdyn_utils, only: build_deterministic, percentile_grid, random_normal_vector
  use tsdyn_var, only: var_model, fit_var
  implicit none
  private
  public :: tvar_model, fit_tvar, forecast_tvar, simulate_tvar, tvar_regimes, tvar_lr_statistic

  type :: tvar_model
    integer :: nvar=0,order=0,nthresh=0,nregime=0,include=include_const,th_delay=1,nobs=0,rank=0
    logical :: common_deterministic=.false.
    character(len=4) :: transition='TAR'
    real(dp),allocatable :: thresholds(:),transition_weights(:),coefficients(:,:,:)
    real(dp),allocatable :: fitted(:,:),residuals(:,:),sigma(:,:),transition_values(:)
    integer,allocatable :: regimes(:),regime_counts(:)
    real(dp)::ssr=0.0_dp,loglik=0.0_dp,aic=0.0_dp,bic=0.0_dp
  end type tvar_model
contains
  pure integer function classify(z,thr) result(r)
    real(dp),intent(in)::z,thr(:)
    integer::j
    r=size(thr)+1
    do j=1,size(thr);if(z<=thr(j))then;r=j;return;end if;end do
  end function classify

  subroutine build_tvar_design(y,p,include,delay,transition,weights,thresholds,common_det,design,target,z,reg,counts,info)
    real(dp),intent(in)::y(:,:),weights(:),thresholds(:)
    integer,intent(in)::p,include,delay
    character(len=*),intent(in)::transition
    logical,intent(in)::common_det
    real(dp),allocatable,intent(out)::design(:,:),target(:,:),z(:)
    integer,allocatable,intent(out)::reg(:),counts(:)
    integer,intent(out)::info
    integer::n,k,nobs,nd,nreg,npar,i,j,t,r,off
    real(dp),allocatable::det(:,:)
    n=size(y,1);k=size(y,2);nobs=n-p;nd=n_deterministic(include);nreg=size(thresholds)+1
    if(p<1.or.delay<1.or.delay>p.or.size(weights)/=k.or.nobs<1)then;info=-1;allocate(design(0,0),target(0,0),z(0),reg(0),counts(0));return;end if
    if(common_det)then;npar=nd+nreg*k*p;else;npar=nreg*(nd+k*p);end if
    allocate(design(nobs,npar),target(nobs,k),z(nobs),reg(nobs),counts(nreg));design=0.0_dp;counts=0
    call build_deterministic(nobs,include,det,trend_start=real(p+1,dp))
    do i=1,nobs
      t=p+i;target(i,:)=y(t,:);z(i)=dot_product(y(t-delay,:),weights)
      if(trim(transition)=='MTAR'.or.trim(transition)=='mtar')z(i)=dot_product(y(t-delay,:)-y(t-delay-1,:),weights)
      r=classify(z(i),thresholds);reg(i)=r;counts(r)=counts(r)+1
      if(common_det)then
        if(nd>0)design(i,1:nd)=det(i,:);off=nd+(r-1)*k*p
      else
        off=(r-1)*(nd+k*p);if(nd>0)design(i,off+1:off+nd)=det(i,:);off=off+nd
      end if
      do j=1,p;design(i,off+(j-1)*k+1:off+j*k)=y(t-j,:);end do
    end do
    info=0
  end subroutine build_tvar_design

  subroutine eval_tvar(y,p,include,delay,transition,weights,thresholds,common_det,trim,ssr,beta,fit,res,z,reg,counts,rank,info)
    real(dp),intent(in)::y(:,:),weights(:),thresholds(:),trim
    integer,intent(in)::p,include,delay
    character(len=*),intent(in)::transition
    logical,intent(in)::common_det
    real(dp),intent(out)::ssr
    real(dp),allocatable,intent(out)::beta(:,:),fit(:,:),res(:,:),z(:)
    integer,allocatable,intent(out)::reg(:),counts(:)
    integer,intent(out)::rank,info
    real(dp),allocatable::design(:,:),target(:,:)
    integer::r,istat,minpar,nobs
    call build_tvar_design(y,p,include,delay,transition,weights,thresholds,common_det,design,target,z,reg,counts,istat)
    if(istat/=0)then;info=istat;ssr=huge(1.0_dp);allocate(beta(0,0),fit(0,0),res(0,0));rank=0;return;end if
    nobs=size(target,1);minpar=size(y,2)*p+merge(0,n_deterministic(include),common_det)
    do r=1,size(counts)
      if(real(counts(r),dp)<trim*real(nobs,dp).or.counts(r)<=minpar)then;info=1;ssr=huge(1.0_dp);allocate(beta(0,0),fit(0,0),res(0,0));rank=0;return;end if
    end do
    call ols_fit(design,target,beta,fit,res,rank,ssr,info)
  end subroutine eval_tvar

  subroutine fit_tvar(y,p,include,nthresh,model,info,thresholds,th_delay,transition,weights,trim,ngrid,common_deterministic)
    real(dp),intent(in)::y(:,:)
    integer,intent(in)::p,include,nthresh
    type(tvar_model),intent(out)::model
    integer,intent(out)::info
    real(dp),intent(in),optional::thresholds(:),weights(:),trim
    integer,intent(in),optional::th_delay,ngrid
    character(len=*),intent(in),optional::transition
    logical,intent(in),optional::common_deterministic
    integer::k,nreg,delay,ng,i,j,istat,rank,nd,r,off
    real(dp)::tr,best,ss,det_sigma
    logical::common_det
    character(len=4)::trans
    real(dp),allocatable::w(:),grid(:),thr(:),bestthr(:),beta(:,:),fit(:,:),res(:,:),z(:),bb(:,:),ff(:,:),rr(:,:),zz(:)
    integer,allocatable::reg(:),counts(:),breg(:),bcounts(:)
    k=size(y,2);nreg=nthresh+1;delay=1;if(present(th_delay))delay=th_delay;ng=40;if(present(ngrid))ng=max(3,ngrid)
    tr=0.1_dp;if(present(trim))tr=trim;trans='TAR';if(present(transition))trans=adjustl(transition);common_det=.false.;if(present(common_deterministic))common_det=common_deterministic
    if(nthresh<1.or.nthresh>2)then;info=-1;return;end if
    allocate(w(k));w=0.0_dp;w(1)=1.0_dp;if(present(weights))then;if(size(weights)/=k)then;info=-2;return;end if;w=weights;end if
    best=huge(1.0_dp);allocate(bestthr(nthresh),thr(nthresh));bestthr=0.0_dp
    if(present(thresholds))then
      if(size(thresholds)/=nthresh)then;info=-3;return;end if
      call eval_tvar(y,p,include,delay,trans,w,thresholds,common_det,tr,best,bb,ff,rr,zz,breg,bcounts,rank,istat);bestthr=thresholds
      if(istat/=0)then;info=istat;return;end if
    else
      call tvar_grid(y,p,delay,trans,w,tr,ng,grid,istat);if(istat/=0)then;info=istat;return;end if
      if(nthresh==1)then
        do i=1,size(grid);thr(1)=grid(i);call try_threshold(thr);end do
      else
        do i=1,size(grid)-1;do j=i+1,size(grid);thr=[grid(i),grid(j)];call try_threshold(thr);end do;end do
      end if
      if(best>=0.5_dp*huge(1.0_dp))then;info=2;return;end if
    end if
    nd=n_deterministic(include);model%nvar=k;model%order=p;model%nthresh=nthresh;model%nregime=nreg;model%include=include;model%th_delay=delay;model%common_deterministic=common_det;model%transition=trans
    allocate(model%thresholds(nthresh),model%transition_weights(k),model%coefficients(nd+k*p,nreg,k));model%thresholds=bestthr;model%transition_weights=w;model%coefficients=0.0_dp
    if(common_det)then
      do r=1,nreg
        if(nd>0)model%coefficients(1:nd,r,:)=bb(1:nd,:)
        off=nd+(r-1)*k*p;model%coefficients(nd+1:nd+k*p,r,:)=bb(off+1:off+k*p,:)
      end do
    else
      do r=1,nreg;off=(r-1)*(nd+k*p);model%coefficients(:,r,:)=bb(off+1:off+nd+k*p,:);end do
    end if
    model%nobs=size(ff,1);model%rank=rank;model%ssr=best
    allocate(model%fitted(size(ff,1),k),model%residuals(size(rr,1),k),model%transition_values(size(zz)),model%regimes(size(breg)),model%regime_counts(size(bcounts)))
    model%fitted=ff;model%residuals=rr;model%transition_values=zz;model%regimes=breg;model%regime_counts=bcounts;call covariance_matrix(rr,model%sigma,center=.false.)
    det_sigma=determinant_spd(model%sigma);model%loglik=-0.5_dp*real(model%nobs,dp)*(real(k,dp)*(log(2.0_dp*acos(-1.0_dp))+1.0_dp)+log(max(det_sigma,tiny(1.0_dp))))
    model%aic=-2.0_dp*model%loglik+2.0_dp*real(rank+nthresh+k*(k+1)/2,dp);model%bic=-2.0_dp*model%loglik+log(real(model%nobs,dp))*real(rank+nthresh+k*(k+1)/2,dp);info=0
  contains
    subroutine try_threshold(thr)
      real(dp),intent(in)::thr(:)
      call eval_tvar(y,p,include,delay,trans,w,thr,common_det,tr,ss,beta,fit,res,z,reg,counts,rank,istat)
      if(istat==0.and.ss<best)then
        best=ss;bestthr=thr;call move_all(beta,fit,res,z,reg,counts,bb,ff,rr,zz,breg,bcounts)
      end if
    end subroutine try_threshold
  end subroutine fit_tvar

  subroutine move_all(beta,fit,res,z,reg,counts,bb,ff,rr,zz,breg,bcounts)
    real(dp),allocatable,intent(inout)::beta(:,:),fit(:,:),res(:,:),z(:),bb(:,:),ff(:,:),rr(:,:),zz(:)
    integer,allocatable,intent(inout)::reg(:),counts(:),breg(:),bcounts(:)
    if(allocated(bb))deallocate(bb);if(allocated(ff))deallocate(ff);if(allocated(rr))deallocate(rr);if(allocated(zz))deallocate(zz);if(allocated(breg))deallocate(breg);if(allocated(bcounts))deallocate(bcounts)
    call move_alloc(beta,bb);call move_alloc(fit,ff);call move_alloc(res,rr);call move_alloc(z,zz);call move_alloc(reg,breg);call move_alloc(counts,bcounts)
  end subroutine move_all

  subroutine tvar_grid(y,p,delay,transition,w,trimv,ngrid,grid,info)
    real(dp),intent(in)::y(:,:),w(:),trimv
    integer,intent(in)::p,delay,ngrid
    character(len=*),intent(in)::transition
    real(dp),allocatable,intent(out)::grid(:)
    integer,intent(out)::info
    real(dp),allocatable::z(:)
    integer::i,t
    allocate(z(size(y,1)-p))
    do i=1,size(z);t=p+i;z(i)=dot_product(y(t-delay,:),w);if(trim(transition)=='MTAR'.or.trim(transition)=='mtar')z(i)=dot_product(y(t-delay,:)-y(t-delay-1,:),w);end do
    call percentile_grid(z,trimv,ngrid,grid);info=0
  end subroutine tvar_grid

  real(dp) function determinant_spd(a) result(det)
    use tsdyn_linalg, only: cholesky_lower
    real(dp),intent(in)::a(:,:)
    real(dp),allocatable::l(:,:)
    integer::i,istat
    call cholesky_lower(a,l,istat);if(istat/=0)then;det=tiny(1.0_dp);return;end if;det=1.0_dp;do i=1,size(l,1);det=det*l(i,i)**2;end do
  end function determinant_spd

  pure function tvar_value(model,work,t,r) result(v)
    type(tvar_model),intent(in)::model
    real(dp),intent(in)::work(:,:)
    integer,intent(in)::t,r
    real(dp)::v(model%nvar)
    integer::nd,j,k
    nd=n_deterministic(model%include);k=model%nvar;v=0.0_dp
    select case(model%include)
    case(include_const);v=model%coefficients(1,r,:)
    case(include_trend);v=real(t,dp)*model%coefficients(1,r,:)
    case(include_both);v=model%coefficients(1,r,:)+real(t,dp)*model%coefficients(2,r,:)
    end select
    do j=1,model%order;v=v+matmul(work(t-j,:),model%coefficients(nd+(j-1)*k+1:nd+j*k,r,:));end do
  end function tvar_value

  subroutine forecast_tvar(model,history,h,forecast,regimes,info)
    type(tvar_model),intent(in)::model
    real(dp),intent(in)::history(:,:)
    integer,intent(in)::h
    real(dp),allocatable,intent(out)::forecast(:,:)
    integer,allocatable,intent(out),optional::regimes(:)
    integer,intent(out)::info
    real(dp),allocatable::work(:,:)
    real(dp)::z
    integer::t,idx,r
    if(size(history,1)<model%order+1.or.size(history,2)/=model%nvar.or.h<1)then;info=-1;allocate(forecast(0,0));if(present(regimes))allocate(regimes(0));return;end if
    allocate(work(size(history,1)+h,model%nvar),forecast(h,model%nvar));work(1:size(history,1),:)=history;if(present(regimes))allocate(regimes(h))
    do t=1,h;idx=size(history,1)+t;z=dot_product(work(idx-model%th_delay,:),model%transition_weights);if(trim(model%transition)=='MTAR')z=dot_product(work(idx-model%th_delay,:)-work(idx-model%th_delay-1,:),model%transition_weights);r=classify(z,model%thresholds);work(idx,:)=tvar_value(model,work,idx,r);forecast(t,:)=work(idx,:);if(present(regimes))regimes(t)=r;end do
    info=0
  end subroutine forecast_tvar

  subroutine simulate_tvar(model,n,y,info,innov,start)
    type(tvar_model),intent(in)::model
    integer,intent(in)::n
    real(dp),allocatable,intent(out)::y(:,:)
    integer,intent(out)::info
    real(dp),intent(in),optional::innov(:,:),start(:,:)
    real(dp),allocatable::work(:,:),e(:)
    real(dp)::z
    integer::p,t,r
    p=model%order;if(n<1)then;info=-1;allocate(y(0,0));return;end if
    allocate(work(n+p,model%nvar),e(model%nvar));work=0.0_dp
    if(present(start))then;if(size(start,1)<p.or.size(start,2)/=model%nvar)then;info=-2;allocate(y(0,0));return;end if;work(1:p,:)=start(size(start,1)-p+1:,:);end if
    do t=p+1,n+p;z=dot_product(work(t-model%th_delay,:),model%transition_weights);if(trim(model%transition)=='MTAR')z=dot_product(work(t-model%th_delay,:)-work(t-model%th_delay-1,:),model%transition_weights);r=classify(z,model%thresholds);if(present(innov))then;e=innov(t-p,:);else;call random_normal_vector(e);end if;work(t,:)=tvar_value(model,work,t,r)+e;end do
    allocate(y(n,model%nvar));y=work(p+1:,:);info=0
  end subroutine simulate_tvar

  subroutine tvar_regimes(model,y,regimes,info)
    type(tvar_model),intent(in)::model
    real(dp),intent(in)::y(:,:)
    integer,allocatable,intent(out)::regimes(:)
    integer,intent(out)::info
    real(dp)::z
    integer::i,t
    if(size(y,1)<=model%order)then;info=-1;allocate(regimes(0));return;end if
    allocate(regimes(size(y,1)-model%order))
    do i=1,size(regimes);t=model%order+i;z=dot_product(y(t-model%th_delay,:),model%transition_weights);if(trim(model%transition)=='MTAR')z=dot_product(y(t-model%th_delay,:)-y(t-model%th_delay-1,:),model%transition_weights);regimes(i)=classify(z,model%thresholds);end do
    info=0
  end subroutine tvar_regimes

  subroutine tvar_lr_statistic(y,p,include,statistic,threshold,info,ngrid)
    real(dp),intent(in)::y(:,:)
    integer,intent(in)::p,include
    real(dp),intent(out)::statistic,threshold
    integer,intent(out)::info
    integer,intent(in),optional::ngrid
    type(var_model)::lin
    type(tvar_model)::tv
    integer::istat,ng
    ng=40;if(present(ngrid))ng=ngrid;call fit_var(y,p,include,lin,istat);if(istat/=0)then;info=istat;return;end if
    call fit_tvar(y,p,include,1,tv,istat,ngrid=ng);if(istat/=0)then;info=istat;return;end if
    statistic=real(tv%nobs,dp)*log(lin%ssr/tv%ssr);threshold=tv%thresholds(1);info=0
  end subroutine tvar_lr_statistic
end module tsdyn_tvar
