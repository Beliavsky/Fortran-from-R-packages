! SPDX-License-Identifier: GPL-2.0-or-later
! Numerical translation derived from tsDyn, copyright its original authors.
! This file may be redistributed and/or modified under GPL version 2 or later.
module tsdyn_lstar
  use tsdyn_kinds, only: dp, n_deterministic, include_none, include_const, include_trend, include_both
  use tsdyn_linalg, only: ols_fit
  use tsdyn_utils, only: sigmoid, percentile_grid, build_deterministic, random_normal
  implicit none
  private
  public :: lstar_model, fit_lstar, forecast_lstar, simulate_lstar, lstar_transition

  type :: lstar_model
    integer :: p_low=0,p_high=0,pmax=0,include=include_const,th_delay=0
    real(dp) :: gamma=1.0_dp,threshold=0.0_dp,ssr=0.0_dp,sigma2=0.0_dp,aic=0.0_dp,bic=0.0_dp
    real(dp),allocatable :: coef_low(:),coef_high(:),fitted(:),residuals(:),transition_values(:)
    integer :: nobs=0,rank=0
  end type lstar_model
contains
  pure elemental real(dp) function lstar_transition(z,gamma,threshold) result(g)
    real(dp),intent(in)::z,gamma,threshold
    g=sigmoid(gamma*(z-threshold))
  end function lstar_transition

  subroutine build_lstar_data(x,p_low,p_high,include,delay,gamma,threshold,design,yy,z,g,info)
    real(dp),intent(in)::x(:),gamma,threshold
    integer,intent(in)::p_low,p_high,include,delay
    real(dp),allocatable,intent(out)::design(:,:),yy(:),z(:),g(:)
    integer,intent(out)::info
    integer::pmax,nobs,nd,i,j,t,nl,nh
    real(dp),allocatable::det(:,:)
    pmax=max(p_low,p_high);nobs=size(x)-pmax;nd=n_deterministic(include);nl=nd+p_low;nh=nd+p_high
    if(nobs<1.or.delay<0.or.delay>=pmax)then;info=-1;allocate(design(0,0),yy(0),z(0),g(0));return;end if
    allocate(design(nobs,nl+nh),yy(nobs),z(nobs),g(nobs));design=0.0_dp
    call build_deterministic(nobs,include,det,trend_start=real(pmax+1,dp))
    do i=1,nobs
      t=pmax+i;yy(i)=x(t);z(i)=x(t-1-delay);g(i)=lstar_transition(z(i),gamma,threshold)
      if(nd>0)then
        design(i,1:nd)=det(i,:)
        design(i,nl+1:nl+nd)=det(i,:)*g(i)
      end if
      do j=1,p_low;design(i,nd+j)=x(t-j);end do
      do j=1,p_high;design(i,nl+nd+j)=x(t-j)*g(i);end do
    end do
    info=0
  end subroutine build_lstar_data

  subroutine evaluate_lstar(x,p_low,p_high,include,delay,gamma,threshold,ssr,beta,fit,res,z,g,rank,info)
    real(dp),intent(in)::x(:),gamma,threshold
    integer,intent(in)::p_low,p_high,include,delay
    real(dp),intent(out)::ssr
    real(dp),allocatable,intent(out)::beta(:,:),fit(:,:),res(:,:),z(:),g(:)
    integer,intent(out)::rank,info
    real(dp),allocatable::d(:,:),yy(:),ym(:,:)
    call build_lstar_data(x,p_low,p_high,include,delay,gamma,threshold,d,yy,z,g,info)
    if(info/=0)then;ssr=huge(1.0_dp);allocate(beta(0,0),fit(0,0),res(0,0));rank=0;return;end if
    allocate(ym(size(yy),1));ym(:,1)=yy
    call ols_fit(d,ym,beta,fit,res,rank,ssr,info)
    if(min(sum(g)/real(size(g),dp),1.0_dp-sum(g)/real(size(g),dp))<0.02_dp)ssr=ssr+1.0e6_dp
  end subroutine evaluate_lstar

  subroutine fit_lstar(x,p_low,p_high,include,model,info,th_delay,gamma,threshold,ngamma,nthreshold,trim)
    real(dp),intent(in)::x(:)
    integer,intent(in)::p_low,p_high,include
    type(lstar_model),intent(out)::model
    integer,intent(out)::info
    integer,intent(in),optional::th_delay,ngamma,nthreshold
    real(dp),intent(in),optional::gamma,threshold,trim
    integer::delay,ng,nt,i,j,istat,rank,nd,nl
    real(dp)::tr,gam,th,best,ss,bg,bt,stepg,stept
    real(dp),allocatable::zraw(:),ths(:),gams(:),beta(:,:),fit(:,:),res(:,:),z(:),gv(:)
    real(dp),allocatable::bb(:,:),ff(:,:),rr(:,:),zz(:),gg(:)
    integer::pmax,t
    delay=0;if(present(th_delay))delay=th_delay
    ng=30;if(present(ngamma))ng=max(2,ngamma);nt=60;if(present(nthreshold))nt=max(2,nthreshold)
    tr=0.1_dp;if(present(trim))tr=trim
    pmax=max(p_low,p_high)
    if(size(x)<=pmax+4.or.delay>=pmax)then;info=-1;return;end if
    allocate(zraw(size(x)-pmax))
    do i=1,size(zraw);t=pmax+i;zraw(i)=x(t-1-delay);end do
    call percentile_grid(zraw,tr,nt,ths)
    allocate(gams(ng))
    if(ng==1)then;gams(1)=1.0_dp;else
      do i=1,ng;gams(i)=exp(log(0.5_dp)+real(i-1,dp)/real(ng-1,dp)*(log(100.0_dp)-log(0.5_dp)));end do
    end if
    best=huge(1.0_dp);bg=1.0_dp;bt=ths(max(1,size(ths)/2))
    if(present(gamma).and.present(threshold))then
      call evaluate_lstar(x,p_low,p_high,include,delay,gamma,threshold,best,bb,ff,rr,zz,gg,rank,istat)
      bg=gamma;bt=threshold
    else
      do i=1,size(gams)
        do j=1,size(ths)
          call evaluate_lstar(x,p_low,p_high,include,delay,gams(i),ths(j),ss,beta,fit,res,z,gv,rank,istat)
          if(istat==0.and.ss<best)then
            best=ss;bg=gams(i);bt=ths(j)
            call move_arrays(beta,fit,res,z,gv,bb,ff,rr,zz,gg)
          end if
        end do
      end do
      ! Deterministic pattern refinement in log-gamma and threshold.
      stepg=0.5_dp;stept=max(1.0e-8_dp,0.1_dp*(maxval(zraw)-minval(zraw)))
      do j=1,30
        do i=-1,1
          gam=max(1.0e-4_dp,bg*exp(stepg*real(i,dp)))
          call try_point(gam,bt)
          th=bt+stept*real(i,dp)
          call try_point(bg,th)
        end do
        stepg=stepg*0.7_dp;stept=stept*0.7_dp
      end do
    end if
    if(best>=0.5_dp*huge(1.0_dp))then;info=2;return;end if
    nd=n_deterministic(include);nl=nd+p_low
    model%p_low=p_low;model%p_high=p_high;model%pmax=pmax;model%include=include;model%th_delay=delay
    model%gamma=bg;model%threshold=bt;model%ssr=best;model%nobs=size(ff,1);model%rank=rank
    model%sigma2=best/real(max(1,model%nobs-rank),dp)
    allocate(model%coef_low(nl),model%coef_high(nd+p_high),model%fitted(model%nobs),model%residuals(model%nobs),model%transition_values(model%nobs))
    model%coef_low=bb(1:nl,1);model%coef_high=bb(nl+1:,1);model%fitted=ff(:,1);model%residuals=rr(:,1);model%transition_values=gg
    model%aic=real(model%nobs,dp)*log(max(best/real(model%nobs,dp),tiny(1.0_dp)))+2.0_dp*real(rank+3,dp)
    model%bic=real(model%nobs,dp)*log(max(best/real(model%nobs,dp),tiny(1.0_dp)))+log(real(model%nobs,dp))*real(rank+3,dp)
    info=0
  contains
    subroutine try_point(gam,th)
      real(dp),intent(in)::gam,th
      call evaluate_lstar(x,p_low,p_high,include,delay,gam,th,ss,beta,fit,res,z,gv,rank,istat)
      if(istat==0.and.ss<best)then
        best=ss;bg=gam;bt=th;call move_arrays(beta,fit,res,z,gv,bb,ff,rr,zz,gg)
      end if
    end subroutine try_point
    subroutine move_arrays(beta,fit,res,z,g,bb,ff,rr,zz,gg)
      real(dp),allocatable,intent(inout)::beta(:,:),fit(:,:),res(:,:),z(:),g(:),bb(:,:),ff(:,:),rr(:,:),zz(:),gg(:)
      if(allocated(bb))deallocate(bb);if(allocated(ff))deallocate(ff);if(allocated(rr))deallocate(rr);if(allocated(zz))deallocate(zz);if(allocated(gg))deallocate(gg)
      call move_alloc(beta,bb);call move_alloc(fit,ff);call move_alloc(res,rr);call move_alloc(z,zz);call move_alloc(g,gg)
    end subroutine move_arrays
  end subroutine fit_lstar

  pure real(dp) function lstar_value(model,work,t) result(v)
    type(lstar_model),intent(in)::model
    real(dp),intent(in)::work(:)
    integer,intent(in)::t
    integer::j,nd
    real(dp)::low,high,g,z
    nd=n_deterministic(model%include);low=0.0_dp;high=0.0_dp
    select case(model%include)
    case(include_const);low=model%coef_low(1);high=model%coef_high(1)
    case(include_trend);low=model%coef_low(1)*real(t,dp);high=model%coef_high(1)*real(t,dp)
    case(include_both);low=model%coef_low(1)+model%coef_low(2)*real(t,dp);high=model%coef_high(1)+model%coef_high(2)*real(t,dp)
    end select
    do j=1,model%p_low;low=low+model%coef_low(nd+j)*work(t-j);end do
    do j=1,model%p_high;high=high+model%coef_high(nd+j)*work(t-j);end do
    z=work(t-1-model%th_delay);g=lstar_transition(z,model%gamma,model%threshold)
    v=low+high*g
  end function lstar_value

  subroutine forecast_lstar(model,history,h,forecast,info)
    type(lstar_model),intent(in)::model
    real(dp),intent(in)::history(:)
    integer,intent(in)::h
    real(dp),allocatable,intent(out)::forecast(:)
    integer,intent(out)::info
    real(dp),allocatable::work(:)
    integer::i,t
    if(size(history)<model%pmax.or.h<1)then;info=-1;allocate(forecast(0));return;end if
    allocate(work(size(history)+h),forecast(h));work(1:size(history))=history
    do i=1,h;t=size(history)+i;work(t)=lstar_value(model,work,t);forecast(i)=work(t);end do
    info=0
  end subroutine forecast_lstar

  subroutine simulate_lstar(model,n,x,innov,start,info)
    type(lstar_model),intent(in)::model
    integer,intent(in)::n
    real(dp),allocatable,intent(out)::x(:)
    real(dp),intent(in),optional::innov(:),start(:)
    integer,intent(out)::info
    real(dp),allocatable::work(:)
    integer::p,t
    real(dp)::e
    p=model%pmax;if(n<1)then;info=-1;allocate(x(0));return;end if
    allocate(work(n+p));work=0.0_dp
    if(present(start))then;if(size(start)<p)then;info=-2;allocate(x(0));return;end if;work(1:p)=start(size(start)-p+1:);end if
    do t=p+1,n+p;e=sqrt(max(model%sigma2,0.0_dp))*random_normal();if(present(innov))e=innov(t-p);work(t)=lstar_value(model,work,t)+e;end do
    allocate(x(n));x=work(p+1:);info=0
  end subroutine simulate_lstar
end module tsdyn_lstar
