! SPDX-License-Identifier: GPL-2.0-or-later
! Numerical translation derived from tsDyn, copyright its original authors.
! This file may be redistributed and/or modified under GPL version 2 or later.
module tsdyn_var
  use tsdyn_kinds, only: dp, n_deterministic, include_none, include_const, include_trend, include_both
  use tsdyn_linalg, only: ols_fit, covariance_matrix, cholesky_lower
  use tsdyn_utils, only: build_deterministic, lag_embed_multivariate, differences, random_normal_vector
  implicit none
  private
  public :: var_model, fit_var, simulate_var, forecast_var, var_companion
  public :: impulse_response_var, fevd_var, select_var_lag, bootstrap_var

  type :: var_model
    integer :: order=0,nvar=0,include=include_const,nobs=0,rank=0
    character(len=8) :: specification='level'
    real(dp),allocatable :: coefficients(:,:),fitted(:,:),residuals(:,:),sigma(:,:)
    real(dp)::ssr=0.0_dp,loglik=0.0_dp,aic=0.0_dp,bic=0.0_dp
  end type var_model
contains
  subroutine fit_var(y,p,include,model,info,specification)
    real(dp),intent(in)::y(:,:)
    integer,intent(in)::p,include
    type(var_model),intent(out)::model
    integer,intent(out)::info
    character(len=*),intent(in),optional::specification
    character(len=8)::spec
    real(dp),allocatable::xlag(:,:),target(:,:),dmat(:,:),design(:,:),beta(:,:),fit(:,:),res(:,:),dy(:,:)
    integer::nd,istat,k,npar
    real(dp)::ssr,det
    spec='level';if(present(specification))spec=adjustl(specification)
    nd=n_deterministic(include);k=size(y,2)
    if(p<1.or.size(y,1)<=p.or.nd<0)then;info=-1;return;end if
    select case(trim(spec))
    case('level')
      call lag_embed_multivariate(y,p,xlag,target,istat)
    case('diff')
      call differences(y,dy);call lag_embed_multivariate(dy,p,xlag,target,istat)
    case('ADF','adf')
      call build_adf_data(y,p,xlag,target,istat)
    case default
      info=-2;return
    end select
    if(istat/=0)then;info=istat;return;end if
    call build_deterministic(size(target,1),include,dmat)
    allocate(design(size(target,1),nd+size(xlag,2)))
    if(nd>0)design(:,1:nd)=dmat
    design(:,nd+1:)=xlag
    call ols_fit(design,target,beta,fit,res,model%rank,ssr,istat)
    if(istat/=0)then;info=istat;return;end if
    call covariance_matrix(res,model%sigma,center=.false.)
    model%order=p;model%nvar=k;model%include=include;model%specification=spec;model%nobs=size(target,1)
    allocate(model%coefficients(size(beta,1),k),model%fitted(size(fit,1),k),model%residuals(size(res,1),k))
    model%coefficients=beta;model%fitted=fit;model%residuals=res;model%ssr=ssr
    det=determinant_spd(model%sigma);model%loglik=-0.5_dp*real(model%nobs,dp)*(real(k,dp)*(log(2.0_dp*acos(-1.0_dp))+1.0_dp)+log(max(det,tiny(1.0_dp))))
    npar=size(beta)+k*(k+1)/2
    model%aic=-2.0_dp*model%loglik+2.0_dp*real(npar,dp)
    model%bic=-2.0_dp*model%loglik+log(real(model%nobs,dp))*real(npar,dp)
    info=0
  end subroutine fit_var

  subroutine build_adf_data(y,p,xlag,target,info)
    real(dp),intent(in)::y(:,:)
    integer,intent(in)::p
    real(dp),allocatable,intent(out)::xlag(:,:),target(:,:)
    integer,intent(out)::info
    real(dp),allocatable::dy(:,:)
    integer::n,k,nobs,i,j,t
    call differences(y,dy);n=size(y,1);k=size(y,2);nobs=n-p-1
    if(nobs<1)then;info=-1;allocate(xlag(0,0),target(0,0));return;end if
    allocate(target(nobs,k),xlag(nobs,k*(p+1)))
    do i=1,nobs
      t=p+1+i;target(i,:)=y(t,:)-y(t-1,:);xlag(i,1:k)=y(t-1,:)
      do j=1,p;xlag(i,j*k+1:(j+1)*k)=dy(t-1-j,:);end do
    end do
    info=0
  end subroutine build_adf_data

  real(dp) function determinant_spd(a) result(det)
    real(dp),intent(in)::a(:,:)
    real(dp),allocatable::l(:,:)
    integer::info,i
    call cholesky_lower(a,l,info)
    if(info/=0)then;det=tiny(1.0_dp);else;det=1.0_dp;do i=1,size(l,1);det=det*l(i,i)**2;end do;end if
  end function determinant_spd

  subroutine simulate_var(coef,p,include,n,y,info,innov,start,trend_start)
    real(dp),intent(in)::coef(:,:)
    integer,intent(in)::p,include,n
    real(dp),allocatable,intent(out)::y(:,:)
    integer,intent(out)::info
    real(dp),intent(in),optional::innov(:,:),start(:,:),trend_start
    integer::k,nd,t,j
    real(dp)::tr
    real(dp),allocatable::work(:,:),e(:)
    k=size(coef,2);nd=n_deterministic(include);tr=1.0_dp;if(present(trend_start))tr=trend_start
    if(size(coef,1)/=nd+k*p.or.n<1)then;info=-1;allocate(y(0,0));return;end if
    allocate(work(n+p,k));work=0.0_dp
    if(present(start))then;if(size(start,1)<p.or.size(start,2)/=k)then;info=-2;allocate(y(0,0));return;end if;work(1:p,:)=start(size(start,1)-p+1:,:);end if
    allocate(e(k))
    do t=p+1,n+p
      e=0.0_dp;if(present(innov))then;e=innov(t-p,:);else;call random_normal_vector(e);end if
      work(t,:)=e
      select case(include)
      case(include_const);work(t,:)=work(t,:)+coef(1,:)
      case(include_trend);work(t,:)=work(t,:)+(tr+real(t-p-1,dp))*coef(1,:)
      case(include_both);work(t,:)=work(t,:)+coef(1,:)+(tr+real(t-p-1,dp))*coef(2,:)
      end select
      do j=1,p;work(t,:)=work(t,:)+matmul(work(t-j,:),coef(nd+(j-1)*k+1:nd+j*k,:));end do
    end do
    allocate(y(n,k));y=work(p+1:,:);info=0
  end subroutine simulate_var

  subroutine forecast_var(model,history,h,forecast,info,exogen_future)
    type(var_model),intent(in)::model
    real(dp),intent(in)::history(:,:)
    integer,intent(in)::h
    real(dp),allocatable,intent(out)::forecast(:,:)
    integer,intent(out)::info
    real(dp),intent(in),optional::exogen_future(:,:)
    real(dp),allocatable::work(:,:)
    integer::k,p,nd,t,j,n_hist,lag_index
    if(trim(model%specification)/='level')then;info=-1;allocate(forecast(0,0));return;end if
    k=model%nvar;p=model%order;nd=n_deterministic(model%include);n_hist=size(history,1)
    if(size(history,1)<p.or.size(history,2)/=k.or.h<1)then;info=-2;allocate(forecast(0,0));return;end if
    allocate(work(n_hist+h,k),forecast(h,k));work(1:n_hist,:)=history
    do t=1,h
      work(n_hist+t,:)=0.0_dp
      select case(model%include)
      case(include_const);work(n_hist+t,:)=model%coefficients(1,:)
      case(include_trend);work(n_hist+t,:)=real(model%nobs+t,dp)*model%coefficients(1,:)
      case(include_both);work(n_hist+t,:)=model%coefficients(1,:)+real(model%nobs+t,dp)*model%coefficients(2,:)
      end select
      do j=1,p
        lag_index=n_hist+t-j
        if(lag_index<1)then;info=-3;return;end if
        work(n_hist+t,:)=work(n_hist+t,:)+matmul(work(lag_index,:),model%coefficients(nd+(j-1)*k+1:nd+j*k,:))
      end do
      forecast(t,:)=work(n_hist+t,:)
    end do
    if(present(exogen_future))continue
    info=0
  end subroutine forecast_var

  subroutine var_companion(model,companion,info)
    type(var_model),intent(in)::model
    real(dp),allocatable,intent(out)::companion(:,:)
    integer,intent(out)::info
    integer::k,p,nd,j
    k=model%nvar;p=model%order;nd=n_deterministic(model%include)
    if(trim(model%specification)/='level')then;info=-1;allocate(companion(0,0));return;end if
    allocate(companion(k*p,k*p));companion=0.0_dp
    do j=1,p;companion(1:k,(j-1)*k+1:j*k)=transpose(model%coefficients(nd+(j-1)*k+1:nd+j*k,:));end do
    do j=2,p;companion((j-1)*k+1:j*k,(j-2)*k+1:(j-1)*k)=identity(k);end do
    info=0
  end subroutine var_companion

  pure function identity(n) result(a)
    integer,intent(in)::n
    real(dp)::a(n,n)
    integer::i
    a=0.0_dp;do i=1,n;a(i,i)=1.0_dp;end do
  end function identity

  subroutine impulse_response_var(model,h,irf,info,orthogonal)
    type(var_model),intent(in)::model
    integer,intent(in)::h
    real(dp),allocatable,intent(out)::irf(:,:,:)
    integer,intent(out)::info
    logical,intent(in),optional::orthogonal
    real(dp),allocatable::comp(:,:),power(:,:),impact(:,:),l(:,:)
    logical::ortho
    integer::k,p,s,istat
    k=model%nvar;p=model%order;ortho=.false.;if(present(orthogonal))ortho=orthogonal
    call var_companion(model,comp,istat);if(istat/=0)then;info=istat;return;end if
    allocate(impact(k,k));impact=identity(k)
    if(ortho)then;call cholesky_lower(model%sigma,l,istat);if(istat/=0)then;info=istat;return;end if;impact=l;end if
    allocate(irf(0:h,k,k),power(k*p,k*p));power=identity(k*p)
    do s=0,h;irf(s,:,:)=matmul(power(1:k,1:k),impact);power=matmul(power,comp);end do
    info=0
  end subroutine impulse_response_var

  subroutine fevd_var(model,h,fevd,info)
    type(var_model),intent(in)::model
    integer,intent(in)::h
    real(dp),allocatable,intent(out)::fevd(:,:,:)
    integer,intent(out)::info
    real(dp),allocatable::irf(:,:,:)
    integer::i,j,s,k,istat
    real(dp)::den
    call impulse_response_var(model,h-1,irf,istat,orthogonal=.true.);if(istat/=0)then;info=istat;return;end if
    k=model%nvar;allocate(fevd(h,k,k));fevd=0.0_dp
    do s=1,h
      do i=1,k
        den=sum(irf(0:s-1,i,:)**2)
        do j=1,k;fevd(s,i,j)=sum(irf(0:s-1,i,j)**2)/max(den,tiny(1.0_dp));end do
      end do
    end do
    info=0
  end subroutine fevd_var

  subroutine select_var_lag(y,pmax,include,criterion,best_p,scores,info)
    real(dp),intent(in)::y(:,:)
    integer,intent(in)::pmax,include
    character(len=*),intent(in)::criterion
    integer,intent(out)::best_p,info
    real(dp),allocatable,intent(out)::scores(:)
    type(var_model)::m
    integer::p,istat
    allocate(scores(pmax));scores=huge(1.0_dp)
    do p=1,pmax
      call fit_var(y,p,include,m,istat)
      if(istat==0)then;if(trim(criterion)=='AIC'.or.trim(criterion)=='aic')then;scores(p)=m%aic;else;scores(p)=m%bic;end if;end if
    end do
    best_p=minloc(scores,dim=1);info=0
  end subroutine select_var_lag

  subroutine bootstrap_var(model,history,h,nsim,paths,info,wild)
    type(var_model),intent(in)::model
    real(dp),intent(in)::history(:,:)
    integer,intent(in)::h,nsim
    real(dp),allocatable,intent(out)::paths(:,:,:)
    integer,intent(out)::info
    logical,intent(in),optional::wild
    real(dp),allocatable::work(:,:),e(:)
    real(dp)::u,signv
    integer::s,t,j,idx,k,p,nd,n_hist,lag_index
    logical::do_wild
    k=model%nvar;p=model%order;nd=n_deterministic(model%include);n_hist=size(history,1);do_wild=.false.;if(present(wild))do_wild=wild
    if(size(history,1)<p.or.size(history,2)/=k.or.h<1.or.nsim<1)then;info=-1;allocate(paths(0,0,0));return;end if
    allocate(paths(h,k,nsim),work(n_hist+h,k),e(k))
    do s=1,nsim
      work(1:n_hist,:)=history
      do t=1,h
        call random_number(u);idx=min(model%nobs,1+int(u*real(model%nobs,dp)));e=model%residuals(idx,:)
        if(do_wild)then;call random_number(u);signv=merge(1.0_dp,-1.0_dp,u>=0.5_dp);e=signv*e;end if
        work(n_hist+t,:)=e
        select case(model%include)
        case(include_const);work(n_hist+t,:)=work(n_hist+t,:)+model%coefficients(1,:)
        case(include_trend);work(n_hist+t,:)=work(n_hist+t,:)+real(model%nobs+t,dp)*model%coefficients(1,:)
        case(include_both);work(n_hist+t,:)=work(n_hist+t,:)+model%coefficients(1,:)+real(model%nobs+t,dp)*model%coefficients(2,:)
        end select
        do j=1,p
          lag_index=n_hist+t-j
          if(lag_index<1)then;info=-2;return;end if
          work(n_hist+t,:)=work(n_hist+t,:)+matmul(work(lag_index,:),model%coefficients(nd+(j-1)*k+1:nd+j*k,:))
        end do
        paths(t,:,s)=work(n_hist+t,:)
      end do
    end do
    info=0
  end subroutine bootstrap_var
end module tsdyn_var
