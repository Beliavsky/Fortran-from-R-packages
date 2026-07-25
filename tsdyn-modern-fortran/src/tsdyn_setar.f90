! SPDX-License-Identifier: GPL-2.0-or-later
! Numerical translation derived from tsDyn, copyright its original authors.
! This file may be redistributed and/or modified under GPL version 2 or later.
module tsdyn_setar
  use tsdyn_kinds, only: dp, n_deterministic, include_none, include_const, include_trend, include_both
  use tsdyn_linalg, only: ols_fit
  use tsdyn_utils, only: build_deterministic, percentile_grid, random_normal
  implicit none
  private
  public :: setar_model, fit_setar, simulate_setar, forecast_setar
  public :: setar_regimes, select_setar_orders, setar_lr_statistic

  type :: setar_model
    integer :: nregime=0, nthresh=0, pmax=0, include=include_const, th_delay=0
    character(len=4) :: transition='TAR'
    integer, allocatable :: orders(:)
    real(dp), allocatable :: thresholds(:), coefficients(:,:)
    real(dp), allocatable :: fitted(:), residuals(:), transition_values(:)
    integer, allocatable :: regimes(:), regime_counts(:)
    real(dp) :: ssr=0.0_dp, sigma2=0.0_dp, aic=0.0_dp, bic=0.0_dp
    integer :: rank=0,nobs=0
  end type setar_model
contains
  pure integer function classify_regime(z,thresholds) result(r)
    real(dp),intent(in)::z,thresholds(:)
    integer::j
    r=size(thresholds)+1
    do j=1,size(thresholds)
      if(z<=thresholds(j))then;r=j;return;end if
    end do
  end function classify_regime

  subroutine build_setar_design(x,orders,include,th_delay,transition,thresholds,design,yy,z,reg,counts,info)
    real(dp),intent(in)::x(:),thresholds(:)
    integer,intent(in)::orders(:),include,th_delay
    character(len=*),intent(in)::transition
    real(dp),allocatable,intent(out)::design(:,:),yy(:),z(:)
    integer,allocatable,intent(out)::reg(:),counts(:)
    integer,intent(out)::info
    integer::pmax,nreg,nd,nobs,i,j,t,r,offset,npar
    real(dp),allocatable::det(:,:)
    pmax=maxval(orders);nreg=size(orders);nd=n_deterministic(include)
    nobs=size(x)-pmax
    if(nobs<1.or.size(thresholds)/=nreg-1.or.th_delay<0.or.th_delay>=pmax)then
      info=-1;allocate(design(0,0),yy(0),z(0),reg(0),counts(0));return
    end if
    npar=sum(orders+nd)
    allocate(design(nobs,npar),yy(nobs),z(nobs),reg(nobs),counts(nreg));design=0.0_dp;counts=0
    call build_deterministic(nobs,include,det,trend_start=real(pmax+1,dp))
    do i=1,nobs
      t=pmax+i;yy(i)=x(t)
      z(i)=x(t-1-th_delay)
      if(trim(transition)=='MTAR'.or.trim(transition)=='mtar')z(i)=x(t-1-th_delay)-x(t-2-th_delay)
      r=classify_regime(z(i),thresholds);reg(i)=r;counts(r)=counts(r)+1
      offset=0
      do j=1,r-1;offset=offset+nd+orders(j);end do
      if(nd>0)design(i,offset+1:offset+nd)=det(i,:)
      do j=1,orders(r);design(i,offset+nd+j)=x(t-j);end do
    end do
    info=0
  end subroutine build_setar_design

  subroutine evaluate_thresholds(x,orders,include,th_delay,transition,thresholds,trim,ssr,beta,fit,res,reg,z,counts,rank,info)
    real(dp),intent(in)::x(:),thresholds(:),trim
    integer,intent(in)::orders(:),include,th_delay
    character(len=*),intent(in)::transition
    real(dp),intent(out)::ssr
    real(dp),allocatable,intent(out)::beta(:,:),fit(:,:),res(:,:),z(:)
    integer,allocatable,intent(out)::reg(:),counts(:)
    integer,intent(out)::rank,info
    real(dp),allocatable::design(:,:),yy(:),ym(:,:)
    integer::nobs,r,istat
    call build_setar_design(x,orders,include,th_delay,transition,thresholds,design,yy,z,reg,counts,istat)
    if(istat/=0)then;info=istat;ssr=huge(1.0_dp);allocate(beta(0,0),fit(0,0),res(0,0));rank=0;return;end if
    nobs=size(yy)
    do r=1,size(counts)
      if(real(counts(r),dp)<trim*real(nobs,dp).or.counts(r)<=orders(r)+n_deterministic(include))then
        info=1;ssr=huge(1.0_dp);allocate(beta(0,0),fit(0,0),res(0,0));rank=0;return
      end if
    end do
    allocate(ym(nobs,1));ym(:,1)=yy
    call ols_fit(design,ym,beta,fit,res,rank,ssr,info)
  end subroutine evaluate_thresholds

  subroutine fit_setar(x,orders,include,nthresh,model,info,thresholds,th_delay,transition,trim,ngrid)
    real(dp),intent(in)::x(:)
    integer,intent(in)::orders(:),include,nthresh
    type(setar_model),intent(out)::model
    integer,intent(out)::info
    real(dp),intent(in),optional::thresholds(:),trim
    integer,intent(in),optional::th_delay,ngrid
    character(len=*),intent(in),optional::transition
    integer::delay,ng,i,j,istat,rank,nd,nreg,pmax,offset,r
    real(dp)::tr,best,ss
    character(len=4)::trans
    real(dp),allocatable::grid(:),thr(:),bestthr(:),beta(:,:),fit(:,:),res(:,:),z(:)
    real(dp),allocatable::bb(:,:),ff(:,:),rr(:,:),zz(:)
    integer,allocatable::reg(:),counts(:),breg(:),bcounts(:)

    nreg=nthresh+1;pmax=maxval(orders);nd=n_deterministic(include)
    if(size(orders)/=nreg.or.nthresh<1.or.nthresh>2)then;info=-1;return;end if
    delay=0;if(present(th_delay))delay=th_delay
    ng=40;if(present(ngrid))ng=max(3,ngrid)
    tr=0.15_dp;if(present(trim))tr=trim
    trans='TAR';if(present(transition))trans=adjustl(transition)
    best=huge(1.0_dp);allocate(bestthr(nthresh));bestthr=0.0_dp
    if(present(thresholds))then
      if(size(thresholds)/=nthresh)then;info=-2;return;end if
      call evaluate_thresholds(x,orders,include,delay,trans,thresholds,tr,best,bb,ff,rr,breg,zz,bcounts,rank,istat)
      if(istat/=0)then;info=istat;return;end if
      bestthr=thresholds
    else
      ! Candidate transition values from the aligned sample.
      allocate(thr(nthresh));thr=0.0_dp
      call build_candidate_grid(x,pmax,delay,trans,tr,ng,grid,istat)
      if(istat/=0)then;info=istat;return;end if
      if(nthresh==1)then
        do i=1,size(grid)
          thr(1)=grid(i)
          call evaluate_thresholds(x,orders,include,delay,trans,thr,tr,ss,beta,fit,res,reg,z,counts,rank,istat)
          if(istat==0.and.ss<best)then
            best=ss;bestthr=thr;call move_fit(beta,fit,res,reg,z,counts,bb,ff,rr,breg,zz,bcounts)
          end if
        end do
      else
        do i=1,size(grid)-1
          do j=i+1,size(grid)
            thr=[grid(i),grid(j)]
            call evaluate_thresholds(x,orders,include,delay,trans,thr,tr,ss,beta,fit,res,reg,z,counts,rank,istat)
            if(istat==0.and.ss<best)then
              best=ss;bestthr=thr;call move_fit(beta,fit,res,reg,z,counts,bb,ff,rr,breg,zz,bcounts)
            end if
          end do
        end do
      end if
      if(best>=0.5_dp*huge(1.0_dp))then;info=2;return;end if
    end if
    model%nregime=nreg;model%nthresh=nthresh;model%pmax=pmax;model%include=include;model%th_delay=delay;model%transition=trans
    allocate(model%orders(nreg),model%thresholds(nthresh),model%coefficients(nd+pmax,nreg))
    model%orders=orders;model%thresholds=bestthr;model%coefficients=0.0_dp
    offset=0
    do r=1,nreg
      model%coefficients(1:nd+orders(r),r)=bb(offset+1:offset+nd+orders(r),1)
      offset=offset+nd+orders(r)
    end do
    allocate(model%fitted(size(ff,1)),model%residuals(size(rr,1)),model%transition_values(size(zz)))
    allocate(model%regimes(size(breg)),model%regime_counts(size(bcounts)))
    model%fitted=ff(:,1);model%residuals=rr(:,1);model%transition_values=zz;model%regimes=breg;model%regime_counts=bcounts
    model%ssr=best;model%nobs=size(ff,1);model%rank=rank;model%sigma2=best/real(max(1,model%nobs-rank),dp)
    model%aic=real(model%nobs,dp)*log(max(best/real(model%nobs,dp),tiny(1.0_dp)))+2.0_dp*real(rank+nthresh+1,dp)
    model%bic=real(model%nobs,dp)*log(max(best/real(model%nobs,dp),tiny(1.0_dp)))+log(real(model%nobs,dp))*real(rank+nthresh+1,dp)
    info=0
  contains
    subroutine move_fit(beta,fit,res,reg,z,counts,bb,ff,rr,breg,zz,bcounts)
      real(dp),allocatable,intent(inout)::beta(:,:),fit(:,:),res(:,:),z(:),bb(:,:),ff(:,:),rr(:,:),zz(:)
      integer,allocatable,intent(inout)::reg(:),counts(:),breg(:),bcounts(:)
      if(allocated(bb))deallocate(bb);if(allocated(ff))deallocate(ff);if(allocated(rr))deallocate(rr);if(allocated(zz))deallocate(zz)
      if(allocated(breg))deallocate(breg);if(allocated(bcounts))deallocate(bcounts)
      call move_alloc(beta,bb);call move_alloc(fit,ff);call move_alloc(res,rr);call move_alloc(z,zz);call move_alloc(reg,breg);call move_alloc(counts,bcounts)
    end subroutine move_fit
  end subroutine fit_setar

  subroutine build_candidate_grid(x,pmax,delay,transition,trimv,ngrid,grid,info)
    real(dp),intent(in)::x(:),trimv
    integer,intent(in)::pmax,delay,ngrid
    character(len=*),intent(in)::transition
    real(dp),allocatable,intent(out)::grid(:)
    integer,intent(out)::info
    real(dp),allocatable::z(:)
    integer::i,t,nobs
    nobs=size(x)-pmax
    if(nobs<1.or.delay>=pmax)then;info=-1;allocate(grid(0));return;end if
    allocate(z(nobs))
    do i=1,nobs;t=pmax+i;z(i)=x(t-1-delay);if(trim(transition)=='MTAR'.or.trim(transition)=='mtar')z(i)=x(t-1-delay)-x(t-2-delay);end do
    call percentile_grid(z,trimv,ngrid,grid);info=0
  end subroutine build_candidate_grid

  subroutine setar_regimes(model,x,regimes,info)
    type(setar_model),intent(in)::model
    real(dp),intent(in)::x(:)
    integer,allocatable,intent(out)::regimes(:)
    integer,intent(out)::info
    integer::i,t
    real(dp)::z
    if(size(x)<=model%pmax)then;info=-1;allocate(regimes(0));return;end if
    allocate(regimes(size(x)-model%pmax))
    do i=1,size(regimes)
      t=model%pmax+i;z=x(t-1-model%th_delay)
      if(trim(model%transition)=='MTAR')z=x(t-1-model%th_delay)-x(t-2-model%th_delay)
      regimes(i)=classify_regime(z,model%thresholds)
    end do
    info=0
  end subroutine setar_regimes

  pure real(dp) function regime_value(model,work,t,r) result(v)
    type(setar_model),intent(in)::model
    real(dp),intent(in)::work(:)
    integer,intent(in)::t,r
    integer::nd,j
    nd=n_deterministic(model%include);v=0.0_dp
    select case(model%include)
    case(include_const);v=model%coefficients(1,r)
    case(include_trend);v=model%coefficients(1,r)*real(t,dp)
    case(include_both);v=model%coefficients(1,r)+model%coefficients(2,r)*real(t,dp)
    end select
    do j=1,model%orders(r);v=v+model%coefficients(nd+j,r)*work(t-j);end do
  end function regime_value

  subroutine forecast_setar(model,history,h,forecast,regimes,info)
    type(setar_model),intent(in)::model
    real(dp),intent(in)::history(:)
    integer,intent(in)::h
    real(dp),allocatable,intent(out)::forecast(:)
    integer,allocatable,intent(out),optional::regimes(:)
    integer,intent(out)::info
    real(dp),allocatable::work(:)
    integer::t,r,idx
    real(dp)::z
    if(size(history)<model%pmax+1.or.h<1)then;info=-1;allocate(forecast(0));if(present(regimes))allocate(regimes(0));return;end if
    allocate(work(size(history)+h),forecast(h));work(1:size(history))=history
    if(present(regimes))allocate(regimes(h))
    do t=1,h
      idx=size(history)+t;z=work(idx-1-model%th_delay)
      if(trim(model%transition)=='MTAR')z=work(idx-1-model%th_delay)-work(idx-2-model%th_delay)
      r=classify_regime(z,model%thresholds);work(idx)=regime_value(model,work,idx,r);forecast(t)=work(idx)
      if(present(regimes))regimes(t)=r
    end do
    info=0
  end subroutine forecast_setar

  subroutine simulate_setar(model,n,x,innov,start,info)
    type(setar_model),intent(in)::model
    integer,intent(in)::n
    real(dp),allocatable,intent(out)::x(:)
    real(dp),intent(in),optional::innov(:),start(:)
    integer,intent(out)::info
    real(dp),allocatable::work(:)
    integer::p,t,r
    real(dp)::z,e
    p=model%pmax
    if(n<1)then;info=-1;allocate(x(0));return;end if
    allocate(work(n+p));work=0.0_dp
    if(present(start))then;if(size(start)<p)then;info=-2;allocate(x(0));return;end if;work(1:p)=start(size(start)-p+1:);end if
    do t=p+1,n+p
      z=work(t-1-model%th_delay);if(trim(model%transition)=='MTAR')z=work(t-1-model%th_delay)-work(t-2-model%th_delay)
      r=classify_regime(z,model%thresholds);e=sqrt(max(model%sigma2,0.0_dp))*random_normal();if(present(innov))e=innov(t-p)
      work(t)=regime_value(model,work,t,r)+e
    end do
    allocate(x(n));x=work(p+1:);info=0
  end subroutine simulate_setar

  subroutine select_setar_orders(x,pmax,include,criterion,best_orders,best_model,info,ngrid)
    real(dp),intent(in)::x(:)
    integer,intent(in)::pmax,include
    character(len=*),intent(in)::criterion
    integer,allocatable,intent(out)::best_orders(:)
    type(setar_model),intent(out)::best_model
    integer,intent(out)::info
    integer,intent(in),optional::ngrid
    integer::pl,ph,istat,ng
    type(setar_model)::m
    real(dp)::score,best
    ng=30;if(present(ngrid))ng=ngrid;best=huge(1.0_dp);allocate(best_orders(2));best_orders=1
    do pl=1,pmax;do ph=1,pmax
      call fit_setar(x,[pl,ph],include,1,m,istat,trim=0.15_dp,ngrid=ng)
      if(istat==0)then
        if(trim(criterion)=='AIC'.or.trim(criterion)=='aic')then;score=m%aic;else;score=m%bic;end if
        if(score<best)then;best=score;best_orders=[pl,ph];end if
      end if
    end do;end do
    if(best>=0.5_dp*huge(1.0_dp))then
      info=1
    else
      call fit_setar(x,best_orders,include,1,best_model,istat,trim=0.15_dp,ngrid=ng)
      info=istat
    end if
  end subroutine select_setar_orders

  subroutine setar_lr_statistic(x,p,include,statistic,threshold,info,ngrid)
    use tsdyn_ar, only: ar_model,fit_ar
    real(dp),intent(in)::x(:)
    integer,intent(in)::p,include
    real(dp),intent(out)::statistic,threshold
    integer,intent(out)::info
    integer,intent(in),optional::ngrid
    type(ar_model)::lin
    type(setar_model)::tar
    integer::istat,ng
    ng=40;if(present(ngrid))ng=ngrid
    call fit_ar(x,p,include,'level',lin,istat);if(istat/=0)then;info=istat;return;end if
    call fit_setar(x,[p,p],include,1,tar,istat,ngrid=ng);if(istat/=0)then;info=istat;return;end if
    statistic=real(tar%nobs,dp)*log(lin%ssr/tar%ssr);threshold=tar%thresholds(1);info=0
  end subroutine setar_lr_statistic
end module tsdyn_setar
