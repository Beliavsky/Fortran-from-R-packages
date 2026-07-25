! SPDX-License-Identifier: GPL-2.0-or-later
! Numerical translation derived from tsDyn, copyright its original authors.
! This file may be redistributed and/or modified under GPL version 2 or later.
module tsdyn_tests
  use tsdyn_kinds, only: dp, include_none, include_const, include_trend, include_both, n_deterministic
  use tsdyn_linalg, only: covariance_matrix, symmetric_eigen, ols_fit, inverse_matrix
  use tsdyn_utils, only: lag_embed_univariate, seed_random, random_normal, build_deterministic, quantile_linear
  use tsdyn_ar, only: ar_model, fit_ar, select_ar_order, simulate_ar
  use tsdyn_vecm, only: vecm_model
  implicit none
  private
  public :: correlation_integral, delta_statistic, delta_linear_statistic
  public :: delta_shuffle_test, delta_linearity_test
  public :: bbc_test_result, bbc_unit_root_test
  public :: kapshin_result, kapshin_test
  public :: johansen_statistics

  type :: bbc_test_result
    real(dp) :: lr=0.0_dp,wald=0.0_dp,lm=0.0_dp
    real(dp) :: threshold_lr=0.0_dp,threshold_wald=0.0_dp,threshold_lm=0.0_dp
  end type bbc_test_result
  type :: kapshin_result
    real(dp) :: sup=0.0_dp,average=0.0_dp,exp_average=0.0_dp
    integer :: valid_pairs=0
  end type kapshin_result
contains
  real(dp) function correlation_integral(x,m,d,eps,theiler) result(c)
    real(dp),intent(in)::x(:),eps
    integer,intent(in)::m,d,theiler
    real(dp),allocatable::xx(:,:),yy(:)
    integer::info,i,j,n,cnt,den
    call lag_embed_univariate(x,m,d,1,xx,yy,info)
    if(info/=0)then;c=0.0_dp;return;end if
    n=size(xx,1);cnt=0;den=0
    do i=1,n-1
      do j=i+1,n
        if(abs(i-j)<=theiler)cycle
        den=den+1
        if(sqrt(sum((xx(i,:)-xx(j,:))**2))<=eps)cnt=cnt+1
      end do
    end do
    if(den>0)then;c=real(cnt,dp)/real(den,dp);else;c=0.0_dp;end if
  end function correlation_integral

  real(dp) function delta_statistic(x,m,d,eps) result(v)
    real(dp),intent(in)::x(:),eps
    integer,intent(in)::m,d
    real(dp)::cm1,cm,cp1
    if(m<2)then;v=0.0_dp;return;end if
    cm1=correlation_integral(x,m-1,d,eps,1);cm=correlation_integral(x,m,d,eps,1);cp1=correlation_integral(x,m+1,d,eps,1)
    if(cm1*cp1>tiny(1.0_dp))then;v=1.0_dp-cm*cm/(cm1*cp1);else;v=0.0_dp;end if
  end function delta_statistic

  real(dp) function delta_linear_statistic(x,m,d) result(v)
    real(dp),intent(in)::x(:)
    integer,intent(in)::m,d
    real(dp),allocatable::x1(:,:),x2(:,:),y(:),c1(:,:),c2(:,:),e1(:),e2(:),q1(:,:),q2(:,:)
    integer::info
    call lag_embed_univariate(x,m+1,d,1,x1,y,info);if(info/=0)then;v=0.0_dp;return;end if
    call lag_embed_univariate(x,m,d,1,x2,y,info);if(info/=0)then;v=0.0_dp;return;end if
    call covariance_matrix(x1,c1);call covariance_matrix(x2,c2)
    call symmetric_eigen(c1,e1,q1,info);if(info/=0)then;v=0.0_dp;return;end if
    call symmetric_eigen(c2,e2,q2,info);if(info/=0)then;v=0.0_dp;return;end if
    v=1.0_dp-maxval(e1)/max(maxval(e2),tiny(1.0_dp))
  end function delta_linear_statistic

  subroutine delta_shuffle_test(x,m,d,eps,b,pvalue,observed,info,seed)
    real(dp),intent(in)::x(:),eps
    integer,intent(in)::m,d,b
    real(dp),intent(out)::pvalue,observed
    integer,intent(out)::info
    integer,intent(in),optional::seed
    real(dp),allocatable::xb(:)
    real(dp)::u,s
    integer::r,i,j,count
    if(b<1)then;info=-1;return;end if
    if(present(seed))call seed_random(seed)
    observed=delta_statistic(x,m,d,eps);allocate(xb(size(x)));count=0
    do r=1,b
      xb=x
      do i=size(xb),2,-1;call random_number(u);j=1+int(u*real(i,dp));s=xb(i);xb(i)=xb(j);xb(j)=s;end do
      if(delta_statistic(xb,m,d,eps)>=observed)count=count+1
    end do
    pvalue=real(1+count,dp)/real(1+b,dp);info=0
  end subroutine delta_shuffle_test

  subroutine delta_linearity_test(x,m,d,eps,b,pvalue,observed,info,seed,pmax)
    real(dp),intent(in)::x(:),eps
    integer,intent(in)::m,d,b
    real(dp),intent(out)::pvalue,observed
    integer,intent(out)::info
    integer,intent(in),optional::seed,pmax
    integer::pm,best,r,istat,count,nd
    real(dp),allocatable::scores(:),xb(:),innov(:)
    type(ar_model)::ar
    real(dp)::statb
    pm=min(10,max(1,size(x)/10));if(present(pmax))pm=pmax;if(present(seed))call seed_random(seed)
    call select_ar_order(x,pm,include_const,'AIC',best,scores,istat);if(istat/=0)then;info=istat;return;end if
    call fit_ar(x,best,include_const,'level',ar,istat);if(istat/=0)then;info=istat;return;end if
    observed=delta_statistic(x,m,d,eps)-delta_linear_statistic(x,m,d);allocate(innov(size(x)));count=0;nd=n_deterministic(ar%include)
    do r=1,b
      call fill_normal(innov,sqrt(max(ar%sigma2,tiny(1.0_dp))))
      call simulate_ar(ar%coefficients(nd+1:),size(x),xb,intercept=ar%coefficients(1),innov=innov,info=istat)
      if(istat/=0)cycle
      statb=delta_statistic(xb,m,d,eps)-delta_linear_statistic(xb,m,d)
      if(statb>=observed)count=count+1
    end do
    pvalue=real(1+count,dp)/real(1+b,dp);info=0
  end subroutine delta_linearity_test

  subroutine fill_normal(x,sd)
    real(dp),intent(out)::x(:)
    real(dp),intent(in)::sd
    integer::i
    do i=1,size(x);x(i)=sd*random_normal();end do
  end subroutine fill_normal

  subroutine bbc_unit_root_test(x,m,trim_fraction,result,info,ngrid)
    real(dp),intent(in)::x(:),trim_fraction
    integer,intent(in)::m
    type(bbc_test_result),intent(out)::result
    integer,intent(out)::info
    integer,intent(in),optional::ngrid
    real(dp),allocatable::xd(:),dy(:),xlag(:,:),grid(:)
    real(dp)::gam,lr,wald,lm
    integer::nobs,i,ng,istat
    ng=40;if(present(ngrid))ng=max(5,ngrid);if(size(x)<m+10)then;info=-1;return;end if
    allocate(xd(size(x)));xd=x-sum(x)/real(size(x),dp);nobs=size(x)-m-1;allocate(dy(nobs),xlag(nobs,m+1))
    do i=1,nobs
      dy(i)=xd(m+1+i)-xd(m+i);xlag(i,1)=xd(m+i)
      call fill_diff_lags(xd,m+i,m,xlag(i,2:))
    end do
    allocate(grid(ng));do i=1,ng;grid(i)=abs(quantile_linear(xlag(:,1),trim_fraction+(0.5_dp-trim_fraction)*real(i-1,dp)/real(max(1,ng-1),dp)));end do
    result%lr=-huge(1.0_dp);result%wald=-huge(1.0_dp);result%lm=-huge(1.0_dp)
    do i=1,ng
      gam=max(grid(i),tiny(1.0_dp));lr=0.0_dp;wald=0.0_dp;lm=0.0_dp
      call bbc_at_threshold(dy,xlag,gam,trim_fraction,lr,wald,lm,istat);if(istat/=0)cycle
      if(lr>result%lr)then;result%lr=lr;result%threshold_lr=gam;end if
      if(wald>result%wald)then;result%wald=wald;result%threshold_wald=gam;end if
      if(lm>result%lm)then;result%lm=lm;result%threshold_lm=gam;end if
    end do
    if(result%lr<0.0_dp)then;info=2;else;info=0;end if
  end subroutine bbc_unit_root_test

  subroutine fill_diff_lags(x,t,m,out)
    real(dp),intent(in)::x(:)
    integer,intent(in)::t,m
    real(dp),intent(out)::out(:)
    integer::j
    do j=1,m;out(j)=x(t-j+1)-x(t-j);end do
  end subroutine fill_diff_lags

  subroutine bbc_at_threshold(dy,xlag,gam,trimv,lr,wald,lm,info)
    real(dp),intent(in)::dy(:),xlag(:,:),gam,trimv
    real(dp),intent(out)::lr,wald,lm
    integer,intent(out)::info
    real(dp),allocatable::xu(:,:),xr(:,:),ym(:,:),bu(:,:),fu(:,:),ru(:,:),br(:,:),fr(:,:),rr(:,:),xtxi(:,:),rmat(:,:),middle(:,:),rinv(:,:)
    integer::n,m,i,r,ranku,rankr,istat,off
    real(dp)::ssu,ssr,s2
    lr=0.0_dp;wald=0.0_dp;lm=0.0_dp
    n=size(dy);m=size(xlag,2)-1
    allocate(xu(n,3*(m+2)),xr(n,3*(m+1)),ym(n,1));xu=0.0_dp;xr=0.0_dp;ym(:,1)=dy
    do i=1,n
      if(xlag(i,1)<=-gam)then;r=1;else if(xlag(i,1)>gam)then;r=3;else;r=2;end if
      off=(r-1)*(m+2);xu(i,off+1)=1.0_dp;xu(i,off+2:off+m+2)=xlag(i,:)
      off=(r-1)*(m+1);xr(i,off+1)=1.0_dp;xr(i,off+2:off+m+1)=xlag(i,2:)
    end do
    do r=1,3
      if(real(count(any(abs(xu(:,(r-1)*(m+2)+1:r*(m+2)))>0.0_dp,dim=2)),dp)<trimv*real(n,dp))then;info=1;return;end if
    end do
    call ols_fit(xu,ym,bu,fu,ru,ranku,ssu,istat);if(istat/=0)then;info=istat;return;end if
    call ols_fit(xr,ym,br,fr,rr,rankr,ssr,istat);if(istat/=0)then;info=istat;return;end if
    lr=real(n,dp)*log(ssr/ssu);s2=ssu/real(max(1,n-ranku),dp)
    call inverse_matrix(matmul(transpose(xu),xu),xtxi,istat);if(istat/=0)then;info=istat;return;end if
    allocate(rmat(3,size(xu,2)));rmat=0.0_dp
    do r=1,3;rmat(r,(r-1)*(m+2)+2)=1.0_dp;end do
    middle=matmul(rmat,matmul(xtxi,transpose(rmat)));call inverse_matrix(middle,rinv,istat);if(istat/=0)then;info=istat;return;end if
    wald=dot_product(matmul(rmat,bu(:,1)),matmul(rinv,matmul(rmat,bu(:,1))))/max(s2,tiny(1.0_dp))
    lm=dot_product(rr(:,1),matmul(xu,matmul(xtxi,matmul(transpose(xu),rr(:,1)))))/max(ssr/real(n,dp),tiny(1.0_dp));info=0
  end subroutine bbc_at_threshold

  subroutine kapshin_test(x,m,include,cpar,delta_power,min_obs_mid,result,info,npoints)
    real(dp),intent(in)::x(:),cpar,delta_power
    integer,intent(in)::m,include,min_obs_mid
    type(kapshin_result),intent(out)::result
    integer,intent(out)::info
    integer,intent(in),optional::npoints
    real(dp),allocatable::xd(:),det(:,:),ym(:,:),b(:,:),fit(:,:),res(:,:),dy(:),level(:),dlags(:,:),ths(:)
    integer::n,nd,np,i,j,istat,rank,countv
    real(dp)::ssr,p1,p2,w
    n=size(x);nd=n_deterministic(include);allocate(xd(n));xd=x
    if(include/=include_none)then
      call build_deterministic(n,include,det);allocate(ym(n,1));ym(:,1)=x;call ols_fit(det,ym,b,fit,res,rank,ssr,istat);if(istat/=0)then;info=istat;return;end if;xd=res(:,1)
    end if
    call build_adf_arrays(xd,m,dy,level,dlags,istat);if(istat/=0)then;info=istat;return;end if
    p1=0.5_dp-cpar/real(size(level),dp)**delta_power;p2=0.5_dp+cpar/real(size(level),dp)**delta_power;np=max(5,int((p2-p1)*real(size(level),dp)));if(present(npoints))np=max(3,npoints)
    allocate(ths(np));do i=1,np;ths(i)=quantile_linear(level,p1+(p2-p1)*real(i-1,dp)/real(max(1,np-1),dp));end do
    result%sup=-huge(1.0_dp);result%average=0.0_dp;result%exp_average=0.0_dp;countv=0
    do i=1,(np+1)/2
      do j=(np+1)/2,np
        w=0.0_dp;call kapshin_wald(dy,level,dlags,ths(i),ths(j),min_obs_mid,w,istat);if(istat/=0)cycle
        result%sup=max(result%sup,w);result%average=result%average+w;result%exp_average=result%exp_average+exp(min(350.0_dp,w/2.0_dp));countv=countv+1
      end do
    end do
    result%valid_pairs=countv;if(countv>0)then;result%average=result%average/real(countv,dp);result%exp_average=result%exp_average/real(countv,dp);info=0;else;info=2;end if
  end subroutine kapshin_test

  subroutine build_adf_arrays(x,m,dy,level,dlags,info)
    real(dp),intent(in)::x(:)
    integer,intent(in)::m
    real(dp),allocatable,intent(out)::dy(:),level(:),dlags(:,:)
    integer,intent(out)::info
    integer::nobs,i,t
    nobs=size(x)-m-1;if(nobs<2)then;info=-1;allocate(dy(0),level(0),dlags(0,0));return;end if
    allocate(dy(nobs),level(nobs),dlags(nobs,m));do i=1,nobs;t=m+1+i;dy(i)=x(t)-x(t-1);level(i)=x(t-1);call fill_diff_lags(x,t-1,m,dlags(i,:));end do;info=0
  end subroutine build_adf_arrays

  subroutine kapshin_wald(dy,level,dlags,g1,g2,minmid,w,info)
    real(dp),intent(in)::dy(:),level(:),dlags(:,:),g1,g2
    integer,intent(in)::minmid
    real(dp),intent(out)::w
    integer,intent(out)::info
    real(dp),allocatable::x(:,:),ym(:,:),b(:,:),fit(:,:),res(:,:),xtxi(:,:),rmat(:,:),mid(:,:),mi(:,:)
    integer::i,n,m,rank,istat
    real(dp)::ssr,s2
    w=0.0_dp
    n=size(dy);m=size(dlags,2);if(count(level>g1.and.level<=g2)<minmid)then;info=1;return;end if
    allocate(x(n,m+2),ym(n,1));ym(:,1)=dy
    do i=1,n;x(i,1)=merge(level(i),0.0_dp,level(i)<=g1);x(i,2)=merge(level(i),0.0_dp,level(i)>g2);x(i,3:)=dlags(i,:);end do
    call ols_fit(x,ym,b,fit,res,rank,ssr,istat);if(istat/=0)then;info=istat;return;end if
    call inverse_matrix(matmul(transpose(x),x),xtxi,istat);if(istat/=0)then;info=istat;return;end if
    allocate(rmat(2,m+2));rmat=0.0_dp;rmat(1,1)=1.0_dp;rmat(2,2)=1.0_dp;mid=matmul(rmat,matmul(xtxi,transpose(rmat)));call inverse_matrix(mid,mi,istat);if(istat/=0)then;info=istat;return;end if
    s2=ssr/real(max(1,n-rank),dp);w=dot_product(matmul(rmat,b(:,1)),matmul(mi,matmul(rmat,b(:,1))))/max(s2,tiny(1.0_dp));info=0
  end subroutine kapshin_wald

  subroutine johansen_statistics(model,trace_stats,max_eigen_stats,info)
    type(vecm_model),intent(in)::model
    real(dp),allocatable,intent(out)::trace_stats(:),max_eigen_stats(:)
    integer,intent(out)::info
    integer::k,r
    real(dp),allocatable::lam(:)
    if(.not.allocated(model%eigenvalues))then;info=-1;allocate(trace_stats(0),max_eigen_stats(0));return;end if
    k=size(model%eigenvalues);allocate(lam(k));lam=max(0.0_dp,min(1.0_dp-1.0e-12_dp,model%eigenvalues));allocate(trace_stats(0:k-1),max_eigen_stats(0:k-1))
    do r=0,k-1;trace_stats(r)=-real(model%nobs,dp)*sum(log(1.0_dp-lam(r+1:k)));max_eigen_stats(r)=-real(model%nobs,dp)*log(1.0_dp-lam(r+1));end do;info=0
  end subroutine johansen_statistics
end module tsdyn_tests
