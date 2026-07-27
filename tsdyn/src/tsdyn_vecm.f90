! SPDX-License-Identifier: GPL-2.0-or-later
! Numerical translation derived from tsDyn, copyright its original authors.
! This file may be redistributed and/or modified under GPL version 2 or later.
module tsdyn_vecm
  use tsdyn_kinds, only: dp, n_deterministic, include_none, include_const, include_trend, include_both
  use tsdyn_linalg, only: ols_fit, inverse_matrix, covariance_matrix, general_eigen
  use tsdyn_utils, only: build_deterministic, differences, random_normal_vector
  use tsdyn_var, only: var_model, impulse_response_var
  implicit none
  private
  public :: vecm_model, fit_vecm, simulate_vecm, forecast_vecm
  public :: vecm_var_coefficients, impulse_response_vecm, select_vecm_rank

  type :: vecm_model
    integer :: nvar=0,rank_coint=0,lag_diff=0,include=include_const,nobs=0,rank_reg=0
    character(len=10) :: method='2OLS'
    real(dp),allocatable :: beta(:,:),alpha(:,:),coefficients(:,:),fitted(:,:),residuals(:,:),sigma(:,:)
    real(dp),allocatable :: eigenvalues(:)
    real(dp)::ssr=0.0_dp,loglik=0.0_dp,aic=0.0_dp,bic=0.0_dp
  end type vecm_model
contains
  subroutine build_vecm_data(y,p,include,xminus1,delta_lags,dy_target,det,design_short,info)
    real(dp),intent(in)::y(:,:)
    integer,intent(in)::p,include
    real(dp),allocatable,intent(out)::xminus1(:,:),delta_lags(:,:),dy_target(:,:),det(:,:),design_short(:,:)
    integer,intent(out)::info
    real(dp),allocatable::dy(:,:)
    integer::n,k,nobs,nd,i,j,t
    n=size(y,1);k=size(y,2);nd=n_deterministic(include);nobs=n-p-1
    if(p<0.or.nobs<2)then;info=-1;allocate(xminus1(0,0),delta_lags(0,0),dy_target(0,0),det(0,0),design_short(0,0));return;end if
    call differences(y,dy)
    allocate(xminus1(nobs,k),dy_target(nobs,k),delta_lags(nobs,k*p))
    do i=1,nobs
      t=p+1+i;xminus1(i,:)=y(t-1,:);dy_target(i,:)=y(t,:)-y(t-1,:)
      do j=1,p;delta_lags(i,(j-1)*k+1:j*k)=dy(t-1-j,:);end do
    end do
    call build_deterministic(nobs,include,det,trend_start=real(p+2,dp))
    allocate(design_short(nobs,nd+k*p))
    if(nd>0)design_short(:,1:nd)=det
    if(p>0)design_short(:,nd+1:)=delta_lags
    info=0
  end subroutine build_vecm_data

  subroutine fit_vecm(y,p,rank_coint,include,method,model,info,beta_fixed)
    real(dp),intent(in)::y(:,:)
    integer,intent(in)::p,rank_coint,include
    character(len=*),intent(in)::method
    type(vecm_model),intent(out)::model
    integer,intent(out)::info
    real(dp),intent(in),optional::beta_fixed(:,:)
    real(dp),allocatable::x1(:,:),dl(:,:),dyt(:,:),det(:,:),z(:,:),beta(:,:),ect(:,:),design(:,:)
    real(dp),allocatable::bcoef(:,:),fit(:,:),res(:,:),u(:,:),v(:,:),bu(:,:),bv(:,:),fu(:,:),fv(:,:),ru(:,:),rv(:,:)
    real(dp),allocatable::s00(:,:),s11(:,:),s01(:,:),s00i(:,:),s11i(:,:),mat(:,:),wr(:),wi(:),vr(:,:)
    integer,allocatable::ord(:)
    integer::istat,k,nobs,j,rrank,ranku,rankv,npar
    real(dp)::ssr,ssu,ssv,det_sigma
    character(len=10)::meth
    k=size(y,2);meth=adjustl(method)
    if(rank_coint<1.or.rank_coint>=k)then;info=-1;return;end if
    call build_vecm_data(y,p,include,x1,dl,dyt,det,z,istat);if(istat/=0)then;info=istat;return;end if
    nobs=size(dyt,1)
    if(present(beta_fixed))then
      if(size(beta_fixed,1)/=k.or.size(beta_fixed,2)/=rank_coint)then;info=-2;return;end if
      allocate(beta(k,rank_coint));beta=beta_fixed;meth='fixed'
    else if(trim(meth)=='2OLS'.or.trim(meth)=='2ols')then
      if(rank_coint/=1)then;info=-3;return;end if
      call engle_granger_beta(y,beta,istat);if(istat/=0)then;info=istat;return;end if
    else if(trim(meth)=='ML'.or.trim(meth)=='ml'.or.trim(meth)=='Johansen'.or.trim(meth)=='johansen')then
      if(size(z,2)>0)then
        call ols_fit(z,dyt,bu,fu,ru,ranku,ssu,istat);if(istat/=0)then;info=istat;return;end if
        call ols_fit(z,x1,bv,fv,rv,rankv,ssv,istat);if(istat/=0)then;info=istat;return;end if
        u=ru;v=rv
      else
        allocate(u(nobs,k),v(nobs,k));u=dyt;v=x1
      end if
      allocate(s00(k,k),s11(k,k),s01(k,k))
      s00=matmul(transpose(u),u)/real(nobs,dp);s11=matmul(transpose(v),v)/real(nobs,dp);s01=matmul(transpose(u),v)/real(nobs,dp)
      call inverse_matrix(s00,s00i,istat);if(istat/=0)then;info=10+istat;return;end if
      call inverse_matrix(s11,s11i,istat);if(istat/=0)then;info=20+istat;return;end if
      mat=matmul(s11i,matmul(transpose(s01),matmul(s00i,s01)))
      call general_eigen(mat,wr,wi,vr,istat);if(istat/=0)then;info=30+istat;return;end if
      call order_descending(wr,ord)
      allocate(beta(k,rank_coint),model%eigenvalues(k))
      model%eigenvalues=wr(ord)
      do j=1,rank_coint
        beta(:,j)=vr(:,ord(j))/max(abs(vr(1,ord(j))),tiny(1.0_dp))
        if(beta(1,j)<0.0_dp)beta(:,j)=-beta(:,j)
      end do
    else
      info=-4;return
    end if
    ect=matmul(x1,beta)
    allocate(design(nobs,rank_coint+size(z,2)));design(:,1:rank_coint)=ect
    if(size(z,2)>0)design(:,rank_coint+1:)=z
    call ols_fit(design,dyt,bcoef,fit,res,rrank,ssr,istat);if(istat/=0)then;info=istat;return;end if
    model%nvar=k;model%rank_coint=rank_coint;model%lag_diff=p;model%include=include;model%method=meth;model%nobs=nobs;model%rank_reg=rrank
    allocate(model%beta(k,rank_coint),model%alpha(k,rank_coint),model%coefficients(size(bcoef,1),k),model%fitted(nobs,k),model%residuals(nobs,k))
    model%beta=beta;model%alpha=transpose(bcoef(1:rank_coint,:));model%coefficients=bcoef;model%fitted=fit;model%residuals=res;model%ssr=ssr
    call covariance_matrix(res,model%sigma,center=.false.)
    det_sigma=determinant_general(model%sigma)
    model%loglik=-0.5_dp*real(nobs,dp)*(real(k,dp)*(log(2.0_dp*acos(-1.0_dp))+1.0_dp)+log(max(det_sigma,tiny(1.0_dp))))
    npar=size(bcoef)+k*rank_coint+k*(k+1)/2
    model%aic=-2.0_dp*model%loglik+2.0_dp*real(npar,dp);model%bic=-2.0_dp*model%loglik+log(real(nobs,dp))*real(npar,dp)
    info=0
  end subroutine fit_vecm

  subroutine engle_granger_beta(y,beta,info)
    real(dp),intent(in)::y(:,:)
    real(dp),allocatable,intent(out)::beta(:,:)
    integer,intent(out)::info
    real(dp),allocatable::b(:,:),fit(:,:),res(:,:),ym(:,:)
    integer::rank,k
    real(dp)::ssr
    k=size(y,2);allocate(ym(size(y,1),1));ym(:,1)=y(:,1)
    call ols_fit(y(:,2:),ym,b,fit,res,rank,ssr,info)
    if(info/=0)return
    allocate(beta(k,1));beta(1,1)=1.0_dp;beta(2:,1)=-b(:,1)
  end subroutine engle_granger_beta

  subroutine order_descending(x,ord)
    real(dp),intent(in)::x(:)
    integer,allocatable,intent(out)::ord(:)
    integer::i,j,key
    allocate(ord(size(x)));ord=[(i,i=1,size(x))]
    do i=2,size(x);key=ord(i);j=i-1;do while(j>=1.and.x(ord(j))<x(key));ord(j+1)=ord(j);j=j-1;end do;ord(j+1)=key;end do
  end subroutine order_descending

  real(dp) function determinant_general(a) result(det)
    use tsdyn_linalg, only: cholesky_lower
    real(dp),intent(in)::a(:,:)
    real(dp),allocatable::l(:,:)
    integer::i,istat
    call cholesky_lower(a,l,istat);if(istat/=0)then;det=tiny(1.0_dp);return;end if
    det=1.0_dp;do i=1,size(l,1);det=det*l(i,i)**2;end do
  end function determinant_general

  subroutine forecast_vecm(model,history,h,forecast,info)
    type(vecm_model),intent(in)::model
    real(dp),intent(in)::history(:,:)
    integer,intent(in)::h
    real(dp),allocatable,intent(out)::forecast(:,:)
    integer,intent(out)::info
    real(dp),allocatable::work(:,:),reg(:)
    integer::k,p,nd,t,j,off
    k=model%nvar;p=model%lag_diff;nd=n_deterministic(model%include)
    if(size(history,1)<p+1.or.size(history,2)/=k.or.h<1)then;info=-1;allocate(forecast(0,0));return;end if
    allocate(work(size(history,1)+h,k),forecast(h,k),reg(size(model%coefficients,1)));work(1:size(history,1),:)=history
    do t=1,h
      reg=0.0_dp;do j=1,model%rank_coint; reg(j)=dot_product(work(size(history,1)+t-1,:),model%beta(:,j)); end do
      off=model%rank_coint
      select case(model%include)
      case(include_const);reg(off+1)=1.0_dp
      case(include_trend);reg(off+1)=real(model%nobs+t,dp)
      case(include_both);reg(off+1)=1.0_dp;reg(off+2)=real(model%nobs+t,dp)
      end select
      off=off+nd
      do j=1,p
        reg(off+(j-1)*k+1:off+j*k)=work(size(history,1)+t-j,:)-work(size(history,1)+t-j-1,:)
      end do
      work(size(history,1)+t,:)=work(size(history,1)+t-1,:)+matmul(reg,model%coefficients)
      forecast(t,:)=work(size(history,1)+t,:)
    end do
    info=0
  end subroutine forecast_vecm

  subroutine simulate_vecm(model,n,y,info,innov,start)
    type(vecm_model),intent(in)::model
    integer,intent(in)::n
    real(dp),allocatable,intent(out)::y(:,:)
    integer,intent(out)::info
    real(dp),intent(in),optional::innov(:,:),start(:,:)
    real(dp),allocatable::work(:,:),reg(:),e(:)
    integer::k,p,nd,t,j,off
    k=model%nvar;p=model%lag_diff;nd=n_deterministic(model%include)
    if(n<1)then;info=-1;allocate(y(0,0));return;end if
    allocate(work(n+p+1,k),reg(size(model%coefficients,1)),e(k));work=0.0_dp
    if(present(start))then;if(size(start,1)<p+1.or.size(start,2)/=k)then;info=-2;allocate(y(0,0));return;end if;work(1:p+1,:)=start(size(start,1)-p:,:);end if
    do t=p+2,n+p+1
      reg=0.0_dp;do j=1,model%rank_coint; reg(j)=dot_product(work(t-1,:),model%beta(:,j)); end do;off=model%rank_coint
      select case(model%include)
      case(include_const);reg(off+1)=1.0_dp
      case(include_trend);reg(off+1)=real(t-p-1,dp)
      case(include_both);reg(off+1)=1.0_dp;reg(off+2)=real(t-p-1,dp)
      end select
      off=off+nd
      do j=1,p;reg(off+(j-1)*k+1:off+j*k)=work(t-j,:)-work(t-j-1,:);end do
      if(present(innov))then;e=innov(t-p-1,:);else;call random_normal_vector(e);end if
      work(t,:)=work(t-1,:)+matmul(reg,model%coefficients)+e
    end do
    allocate(y(n,k));y=work(p+2:,:);info=0
  end subroutine simulate_vecm

  subroutine vecm_var_coefficients(model,a,info)
    type(vecm_model),intent(in)::model
    real(dp),allocatable,intent(out)::a(:,:,:)
    integer,intent(out)::info
    integer::k,p,nd,j,off
    real(dp),allocatable::pi(:,:),gamma(:,:,:)
    k=model%nvar;p=model%lag_diff;nd=n_deterministic(model%include);off=model%rank_coint+nd
    allocate(pi(k,k));pi=matmul(model%alpha,transpose(model%beta))
    allocate(gamma(k,k,max(1,p)));gamma=0.0_dp
    do j=1,p;gamma(:,:,j)=transpose(model%coefficients(off+(j-1)*k+1:off+j*k,:));end do
    allocate(a(k,k,p+1));a=0.0_dp
    if(p==0)then;a(:,:,1)=identity(k)+pi
    else
      a(:,:,1)=identity(k)+pi+gamma(:,:,1)
      do j=2,p;a(:,:,j)=gamma(:,:,j)-gamma(:,:,j-1);end do
      a(:,:,p+1)=-gamma(:,:,p)
    end if
    info=0
  end subroutine vecm_var_coefficients

  pure function identity(n) result(a)
    integer,intent(in)::n
    real(dp)::a(n,n)
    integer::i
    a=0.0_dp;do i=1,n;a(i,i)=1.0_dp;end do
  end function identity

  subroutine impulse_response_vecm(model,h,irf,info,orthogonal)
    type(vecm_model),intent(in)::model
    integer,intent(in)::h
    real(dp),allocatable,intent(out)::irf(:,:,:)
    integer,intent(out)::info
    logical,intent(in),optional::orthogonal
    real(dp),allocatable::a(:,:,:)
    type(var_model)::vm
    integer::istat,j,k,p
    call vecm_var_coefficients(model,a,istat);if(istat/=0)then;info=istat;return;end if
    k=model%nvar;p=size(a,3);vm%nvar=k;vm%order=p;vm%include=include_none;vm%specification='level';allocate(vm%coefficients(k*p,k),vm%sigma(k,k));vm%sigma=model%sigma
    do j=1,p;vm%coefficients((j-1)*k+1:j*k,:)=transpose(a(:,:,j));end do
    call impulse_response_var(vm,h,irf,info,orthogonal)
  end subroutine impulse_response_vecm

  subroutine select_vecm_rank(y,p,include,rmax,criterion,best_rank,scores,info)
    real(dp),intent(in)::y(:,:)
    integer,intent(in)::p,include,rmax
    character(len=*),intent(in)::criterion
    integer,intent(out)::best_rank,info
    real(dp),allocatable,intent(out)::scores(:)
    type(vecm_model)::m
    integer::r,istat
    allocate(scores(rmax));scores=huge(1.0_dp)
    do r=1,rmax
      call fit_vecm(y,p,r,include,'ML',m,istat)
      if(istat==0)then;if(trim(criterion)=='AIC'.or.trim(criterion)=='aic')then;scores(r)=m%aic;else;scores(r)=m%bic;end if;end if
    end do
    best_rank=minloc(scores,dim=1);info=0
  end subroutine select_vecm_rank
end module tsdyn_vecm
