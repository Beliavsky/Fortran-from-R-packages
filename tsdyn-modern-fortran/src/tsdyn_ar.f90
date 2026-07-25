! SPDX-License-Identifier: GPL-2.0-or-later
! Numerical translation derived from tsDyn, copyright its original authors.
! This file may be redistributed and/or modified under GPL version 2 or later.
module tsdyn_ar
  use tsdyn_kinds, only: dp, n_deterministic, include_none, include_const, include_trend, include_both
  use tsdyn_linalg, only: ols_fit, inverse_matrix, general_eigenvalues
  use tsdyn_utils, only: build_deterministic, lag_embed_univariate, random_normal, quantile_linear
  implicit none
  private
  public :: ar_model, fit_ar, simulate_ar, forecast_ar, ar_roots
  public :: select_ar_order, accuracy_metrics, compute_accuracy, residual_bootstrap_ar

  type :: ar_model
    integer :: order = 0
    integer :: include = include_const
    character(len=8) :: specification = 'level'
    real(dp), allocatable :: coefficients(:)
    real(dp), allocatable :: fitted(:), residuals(:)
    real(dp) :: sigma2 = 0.0_dp
    real(dp) :: ssr = 0.0_dp
    real(dp) :: aic = 0.0_dp
    real(dp) :: bic = 0.0_dp
    integer :: rank = 0
    integer :: nobs = 0
    logical :: stable = .false.
  end type ar_model

  type :: accuracy_metrics
    real(dp) :: me=0.0_dp, rmse=0.0_dp, mae=0.0_dp, mpe=0.0_dp, mape=0.0_dp
  end type accuracy_metrics
contains
  subroutine fit_ar(x,p,include,specification,model,info)
    real(dp),intent(in)::x(:)
    integer,intent(in)::p,include
    character(len=*),intent(in),optional::specification
    type(ar_model),intent(out)::model
    integer,intent(out)::info
    character(len=8)::spec
    real(dp),allocatable::lags(:,:),yy(:),dmat(:,:),design(:,:),ymat(:,:)
    real(dp),allocatable::beta(:,:),fit(:,:),res(:,:),dx(:)
    real(dp),allocatable::wr(:),wi(:)
    integer::istat,nobs,nd,k,j
    real(dp)::ssr

    spec='level'; if(present(specification))spec=adjustl(specification)
    model%order=p; model%include=include; model%specification=spec
    nd=n_deterministic(include)
    if(p<0.or.nd<0)then;info=-1;return;end if
    select case(trim(spec))
    case('level')
      if(p==0)then
        nobs=size(x); allocate(lags(nobs,0),yy(nobs)); yy=x
      else
        call lag_embed_univariate(x,p,1,1,lags,yy,istat)
        if(istat/=0)then;info=istat;return;end if
      end if
      call build_deterministic(size(yy),include,dmat)
      allocate(design(size(yy),nd+p))
      if(nd>0)design(:,1:nd)=dmat
      if(p>0)design(:,nd+1:)=lags
    case('diff')
      if(size(x)<p+2)then;info=-2;return;end if
      allocate(dx(size(x)-1)); dx=x(2:)-x(:size(x)-1)
      if(p==0)then
        allocate(lags(size(dx),0),yy(size(dx))); yy=dx
      else
        call lag_embed_univariate(dx,p,1,1,lags,yy,istat)
      end if
      call build_deterministic(size(yy),include,dmat)
      allocate(design(size(yy),nd+p))
      if(nd>0)design(:,1:nd)=dmat
      if(p>0)design(:,nd+1:)=lags
    case('ADF','adf')
      if(size(x)<p+3)then;info=-3;return;end if
      allocate(dx(size(x)-1)); dx=x(2:)-x(:size(x)-1)
      nobs=size(x)-p-1
      allocate(yy(nobs),design(nobs,nd+1+p))
      yy=dx(p+1:)
      call build_deterministic(nobs,include,dmat)
      if(nd>0)design(:,1:nd)=dmat
      design(:,nd+1)=x(p+1:size(x)-1)
      do j=1,p
        design(:,nd+1+j)=dx(p+1-j:size(dx)-j)
      end do
    case default
      info=-4;return
    end select
    if(size(design,2)==0)then;info=-5;return;end if
    allocate(ymat(size(yy),1)); ymat(:,1)=yy
    call ols_fit(design,ymat,beta,fit,res,model%rank,ssr,istat)
    if(istat/=0)then;info=istat;return;end if
    allocate(model%coefficients(size(beta,1)),model%fitted(size(yy)),model%residuals(size(yy)))
    model%coefficients=beta(:,1);model%fitted=fit(:,1);model%residuals=res(:,1)
    model%ssr=ssr; model%nobs=size(yy); k=size(beta,1)
    model%sigma2=ssr/real(max(1,model%nobs-model%rank),dp)
    model%aic=real(model%nobs,dp)*log(max(ssr/real(model%nobs,dp),tiny(1.0_dp)))+2.0_dp*real(k+1,dp)
    model%bic=real(model%nobs,dp)*log(max(ssr/real(model%nobs,dp),tiny(1.0_dp)))+log(real(model%nobs,dp))*real(k+1,dp)
    model%stable=.false.
    if(trim(spec)=='level'.and.p>0)then
      call ar_roots(model%coefficients(nd+1:nd+p),wr,wi,istat)
      if(istat==0)model%stable=all(sqrt(wr*wr+wi*wi)>1.0_dp)
    else if(p==0)then
      model%stable=.true.
    end if
    info=0
  end subroutine fit_ar

  subroutine ar_roots(phi,wr,wi,info)
    real(dp),intent(in)::phi(:)
    real(dp),allocatable,intent(out)::wr(:),wi(:)
    integer,intent(out)::info
    real(dp),allocatable::comp(:,:)
    integer::p,i
    p=size(phi)
    if(p==0)then;allocate(wr(0),wi(0));info=0;return;end if
    allocate(comp(p,p));comp=0.0_dp
    comp(1,:)=phi
    do i=2,p;comp(i,i-1)=1.0_dp;end do
    call general_eigenvalues(comp,wr,wi,info)
    if(info==0)then
      block
        real(dp), allocatable :: den(:), wr0(:), wi0(:)
        allocate(den(size(wr)),wr0(size(wr)),wi0(size(wi)))
        wr0=wr;wi0=wi;den=wr0*wr0+wi0*wi0
        where(den>tiny(1.0_dp))
          wr=wr0/den
          wi=-wi0/den
        elsewhere
          wr=huge(1.0_dp);wi=0.0_dp
        end where
      end block
    end if
  end subroutine ar_roots

  subroutine simulate_ar(phi,n,x,intercept,trend,innov,start,info)
    real(dp),intent(in)::phi(:)
    integer,intent(in)::n
    real(dp),allocatable,intent(out)::x(:)
    real(dp),intent(in),optional::intercept,trend,innov(:),start(:)
    integer,intent(out)::info
    integer::p,t,j
    real(dp)::c,b,e
    p=size(phi); c=0.0_dp;b=0.0_dp
    if(present(intercept))c=intercept;if(present(trend))b=trend
    if(n<1)then;info=-1;allocate(x(0));return;end if
    allocate(x(n+p));x=0.0_dp
    if(present(start))then
      if(size(start)<p)then;info=-2;return;end if
      x(1:p)=start(size(start)-p+1:)
    end if
    do t=p+1,n+p
      e=random_normal();if(present(innov))e=innov(t-p)
      x(t)=c+b*real(t-p,dp)+e
      do j=1,p;x(t)=x(t)+phi(j)*x(t-j);end do
    end do
    x=x(p+1:);info=0
  end subroutine simulate_ar

  subroutine forecast_ar(model,history,h,forecast,info)
    type(ar_model),intent(in)::model
    real(dp),intent(in)::history(:)
    integer,intent(in)::h
    real(dp),allocatable,intent(out)::forecast(:)
    integer,intent(out)::info
    real(dp),allocatable::work(:)
    integer::p,nd,t,j
    real(dp)::v
    p=model%order;nd=n_deterministic(model%include)
    if(trim(model%specification)/='level'.or.size(history)<p.or.h<1)then;info=-1;allocate(forecast(0));return;end if
    allocate(work(size(history)+h));work(1:size(history))=history;allocate(forecast(h))
    do t=1,h
      v=0.0_dp
      select case(model%include)
      case(include_const);v=model%coefficients(1)
      case(include_trend);v=model%coefficients(1)*real(model%nobs+t,dp)
      case(include_both);v=model%coefficients(1)+model%coefficients(2)*real(model%nobs+t,dp)
      end select
      do j=1,p;v=v+model%coefficients(nd+j)*work(size(history)+t-j);end do
      work(size(history)+t)=v;forecast(t)=v
    end do
    info=0
  end subroutine forecast_ar

  subroutine select_ar_order(x,pmax,include,criterion,best_order,scores,info)
    real(dp),intent(in)::x(:)
    integer,intent(in)::pmax,include
    character(len=*),intent(in)::criterion
    integer,intent(out)::best_order,info
    real(dp),allocatable,intent(out)::scores(:)
    type(ar_model)::m
    integer::p,istat
    allocate(scores(0:pmax));scores=huge(1.0_dp);best_order=0
    do p=0,pmax
      call fit_ar(x,p,include,'level',m,istat)
      if(istat==0)then
        if(trim(criterion)=='AIC'.or.trim(criterion)=='aic')then;scores(p)=m%aic;else;scores(p)=m%bic;end if
      end if
    end do
    best_order=minloc(scores,dim=1)-1;info=0
  end subroutine select_ar_order

  pure function compute_accuracy(fit,true_values) result(a)
    real(dp),intent(in)::fit(:),true_values(:)
    type(accuracy_metrics)::a
    real(dp),allocatable::e(:)
    integer::n
    n=min(size(fit),size(true_values));if(n==0)return
    allocate(e(n));e=fit(1:n)-true_values(1:n)
    a%me=sum(e)/real(n,dp);a%rmse=sqrt(sum(e*e)/real(n,dp));a%mae=sum(abs(e))/real(n,dp)
    where(abs(true_values(1:n))>tiny(1.0_dp))
      e=100.0_dp*e/true_values(1:n)
    elsewhere
      e=0.0_dp
    end where
    a%mpe=sum(e)/real(n,dp);a%mape=sum(abs(e))/real(n,dp)
  end function compute_accuracy

  subroutine residual_bootstrap_ar(model,history,h,nsim,paths,info)
    type(ar_model),intent(in)::model
    real(dp),intent(in)::history(:)
    integer,intent(in)::h,nsim
    real(dp),allocatable,intent(out)::paths(:,:)
    integer,intent(out)::info
    real(dp),allocatable::work(:)
    integer::s,t,j,p,nd,idx
    real(dp)::u,v
    p=model%order;nd=n_deterministic(model%include)
    if(size(history)<p.or.h<1.or.nsim<1)then;info=-1;allocate(paths(0,0));return;end if
    allocate(paths(h,nsim),work(size(history)+h))
    do s=1,nsim
      work(1:size(history))=history
      do t=1,h
        call random_number(u);idx=min(size(model%residuals),1+int(u*real(size(model%residuals),dp)))
        v=model%residuals(idx)
        select case(model%include)
        case(include_const);v=v+model%coefficients(1)
        case(include_trend);v=v+model%coefficients(1)*real(model%nobs+t,dp)
        case(include_both);v=v+model%coefficients(1)+model%coefficients(2)*real(model%nobs+t,dp)
        end select
        do j=1,p;v=v+model%coefficients(nd+j)*work(size(history)+t-j);end do
        work(size(history)+t)=v;paths(t,s)=v
      end do
    end do
    info=0
  end subroutine residual_bootstrap_ar
end module tsdyn_ar
