! SPDX-License-Identifier: GPL-2.0-or-later
! This file is part of a modern Fortran translation of PerformanceAnalytics.
! Copyright (C) original PerformanceAnalytics authors and translation contributors.
! Distributed under GNU GPL version 2 or, at your option, any later version.
module advanced_moments_mod
  use kinds_mod, only: dp
  use statistics_mod, only: mean_value, variance_value, solve_linear_system
  use comoments_mod, only: covariance_matrix, coskewness_matrix, cokurtosis_matrix
  implicit none
  private
  public :: mca_result, nce_result
  public :: structured_covariance, structured_coskewness, structured_cokurtosis
  public :: shrink_covariance, shrink_coskewness, shrink_cokurtosis
  public :: multi_target_shrink_covariance, multi_target_shrink_coskewness, multi_target_shrink_cokurtosis
  public :: ewma_coskewness, ewma_cokurtosis
  public :: m3_mca, m4_mca, nearest_comoment_estimator

  type :: mca_result
    integer :: components = 0
    integer :: iterations = 0
    logical :: converged = .false.
    real(dp), allocatable :: directions(:,:)
    real(dp), allocatable :: m3(:,:)
    real(dp), allocatable :: m4(:,:)
  end type mca_result

  type :: nce_result
    integer :: factors = 0
    integer :: iterations = 0
    logical :: converged = .false.
    real(dp) :: objective = huge(1.0_dp)
    real(dp), allocatable :: loadings(:,:)
    real(dp), allocatable :: residual_variance(:)
    real(dp), allocatable :: factor_skewness(:)
    real(dp), allocatable :: residual_third(:)
    real(dp), allocatable :: factor_kurtosis(:)
    real(dp), allocatable :: residual_fourth(:)
    real(dp), allocatable :: covariance(:,:)
    real(dp), allocatable :: coskewness(:,:)
    real(dp), allocatable :: cokurtosis(:,:)
  end type nce_result
contains
  subroutine centered_data(r, x)
    real(dp), intent(in) :: r(:,:)
    real(dp), allocatable, intent(out) :: x(:,:)
    integer :: j
    allocate(x(size(r,1),size(r,2)))
    do j = 1, size(r,2)
      x(:,j) = r(:,j) - mean_value(r(:,j))
    end do
  end subroutine centered_data

  pure function lower_string(text) result(out)
    character(len=*), intent(in) :: text
    character(len=len(text)) :: out
    integer :: i, c
    out = text
    do i = 1, len(text)
      c = iachar(out(i:i))
      if (c >= iachar('A') .and. c <= iachar('Z')) out(i:i) = achar(c + 32)
    end do
  end function lower_string

  subroutine symmetric_eigen_jacobi(a, values, vectors, ok)
    real(dp), intent(in) :: a(:,:)
    real(dp), allocatable, intent(out) :: values(:), vectors(:,:)
    logical, intent(out) :: ok
    real(dp), allocatable :: b(:,:)
    real(dp) :: app, aqq, apq, tau, t, c, s, bip, biq, vip, viq, off
    integer :: n, i, j, p, q, iter, maxit, imax
    n = size(a,1)
    allocate(values(n), vectors(n,n), b(n,n))
    ok = size(a,2) == n
    if (.not. ok) return
    b = 0.5_dp * (a + transpose(a))
    vectors = 0.0_dp
    do i = 1, n
      vectors(i,i) = 1.0_dp
    end do
    maxit = max(50, 40*n*n)
    do iter = 1, maxit
      off = 0.0_dp; p = 1; q = min(2,n)
      do i = 1, n-1
        do j = i+1, n
          if (abs(b(i,j)) > off) then
            off = abs(b(i,j)); p = i; q = j
          end if
        end do
      end do
      if (off <= 1.0e-12_dp * max(1.0_dp,maxval(abs(b)))) exit
      app = b(p,p); aqq = b(q,q); apq = b(p,q)
      tau = (aqq-app)/(2.0_dp*apq)
      if (tau >= 0.0_dp) then
        t = 1.0_dp/(tau+sqrt(1.0_dp+tau*tau))
      else
        t = -1.0_dp/(-tau+sqrt(1.0_dp+tau*tau))
      end if
      c = 1.0_dp/sqrt(1.0_dp+t*t); s = t*c
      do i = 1, n
        if (i /= p .and. i /= q) then
          bip = b(i,p); biq = b(i,q)
          b(i,p) = c*bip-s*biq; b(p,i)=b(i,p)
          b(i,q) = s*bip+c*biq; b(q,i)=b(i,q)
        end if
      end do
      b(p,p)=c*c*app-2.0_dp*s*c*apq+s*s*aqq
      b(q,q)=s*s*app+2.0_dp*s*c*apq+c*c*aqq
      b(p,q)=0.0_dp; b(q,p)=0.0_dp
      do i=1,n
        vip=vectors(i,p);viq=vectors(i,q)
        vectors(i,p)=c*vip-s*viq
        vectors(i,q)=s*vip+c*viq
      end do
    end do
    values = [(b(i,i),i=1,n)]
    do i = 1, n-1
      imax = i
      do j = i+1, n
        if (values(j) > values(imax)) imax = j
      end do
      if (imax /= i) then
        app=values(i);values(i)=values(imax);values(imax)=app
        b(:,1)=vectors(:,i);vectors(:,i)=vectors(:,imax);vectors(:,imax)=b(:,1)
      end if
    end do
    ok = .true.
  end subroutine symmetric_eigen_jacobi

  subroutine factor_loadings(r, factors, beta, residual)
    real(dp), intent(in) :: r(:,:), factors(:,:)
    real(dp), allocatable, intent(out) :: beta(:,:), residual(:,:)
    real(dp), allocatable :: x(:,:), f(:,:), gram(:,:), rhs(:), coef(:)
    logical :: ok
    integer :: p, k, j
    call centered_data(r,x); call centered_data(factors,f)
    p=size(r,2);k=size(factors,2)
    allocate(beta(p,k),residual(size(r,1),p),gram(k,k),rhs(k),coef(k))
    gram=matmul(transpose(f),f)
    do j=1,p
      rhs=matmul(transpose(f),x(:,j))
      call solve_linear_system(gram,rhs,coef,ok)
      if(ok)then;beta(j,:)=coef;else;beta(j,:)=0.0_dp;end if
    end do
    residual=x-matmul(f,transpose(beta))
  end subroutine factor_loadings

  subroutine pca_factors(r, k, beta, factors, residual)
    real(dp), intent(in) :: r(:,:)
    integer, intent(in) :: k
    real(dp), allocatable, intent(out) :: beta(:,:), factors(:,:), residual(:,:)
    real(dp), allocatable :: x(:,:), cov(:,:), eig(:), u(:,:)
    logical :: ok
    integer :: p, kk, j
    call centered_data(r,x);p=size(r,2);kk=max(0,min(k,p))
    allocate(cov(p,p));call covariance_matrix(x,cov,.false.)
    call symmetric_eigen_jacobi(cov,eig,u,ok)
    allocate(beta(p,kk),factors(size(r,1),kk),residual(size(r,1),p))
    if(.not.ok .or. kk==0)then
      beta=0.0_dp;factors=0.0_dp;residual=x;return
    end if
    do j=1,kk
      beta(:,j)=u(:,j)*sqrt(max(eig(j),0.0_dp))
      if(eig(j)>1.0e-14_dp)then
        factors(:,j)=matmul(x,u(:,j))/sqrt(eig(j))
      else
        factors(:,j)=0.0_dp
      end if
    end do
    residual=x-matmul(factors,transpose(beta))
  end subroutine pca_factors

  subroutine structured_covariance(r, kind, out, factor)
    real(dp), intent(in) :: r(:,:)
    character(len=*), intent(in) :: kind
    real(dp), intent(out) :: out(:,:)
    real(dp), intent(in), optional :: factor(:,:)
    real(dp), allocatable :: x(:,:), sample(:,:), beta(:,:), resid(:,:), f(:,:), fcov(:,:), sd(:)
    real(dp) :: avgvar, rho
    integer :: p, i, j
    character(len=:), allocatable :: key
    call centered_data(r,x);p=size(r,2);allocate(sample(p,p));call covariance_matrix(x,sample,.false.)
    key=trim(lower_string(kind));out=0.0_dp
    select case(key)
    case('indep','independent')
      do i=1,p;out(i,i)=sample(i,i);end do
    case('indepid','independent_equal')
      avgvar=sum([(sample(i,i),i=1,p)])/real(p,dp)
      do i=1,p;out(i,i)=avgvar;end do
    case('cc','constant_correlation')
      allocate(sd(p));do i=1,p;sd(i)=sqrt(max(sample(i,i),0.0_dp));end do
      rho=0.0_dp
      if(p>1)then
        do i=1,p-1;do j=i+1,p
          if(sd(i)*sd(j)>0.0_dp)rho=rho+sample(i,j)/(sd(i)*sd(j))
        end do;end do
        rho=rho/real(p*(p-1)/2,dp)
      end if
      do i=1,p;do j=1,p
        if(i==j)then;out(i,j)=sd(i)*sd(i);else;out(i,j)=rho*sd(i)*sd(j);end if
      end do;end do
    case('observedfactor','factor')
      if(present(factor))then
        call factor_loadings(r,factor,beta,resid);call centered_data(factor,f)
        allocate(fcov(size(f,2),size(f,2)));call covariance_matrix(f,fcov,.false.)
        out=matmul(matmul(beta,fcov),transpose(beta))
        do i=1,p;out(i,i)=out(i,i)+sum(resid(:,i)**2)/real(size(r,1),dp);end do
      else
        out=sample
      end if
    case('latent1factor','pca')
      call pca_factors(r,1,beta,f,resid)
      out=matmul(beta,transpose(beta))
      do i=1,p;out(i,i)=out(i,i)+sum(resid(:,i)**2)/real(size(r,1),dp);end do
    case default
      out=sample
    end select
  end subroutine structured_covariance

  subroutine factor_model_moments(r, factors, m3, m4)
    real(dp), intent(in) :: r(:,:), factors(:,:)
    real(dp), intent(out) :: m3(:,:),m4(:,:)
    real(dp),allocatable::beta(:,:),resid(:,:),f(:,:),fm3(:,:),fm4(:,:),fcov(:,:),ecov(:,:)
    real(dp)::v
    integer::p,k,i,j,l,h,a,b,c,d,col
    call factor_loadings(r,factors,beta,resid);call centered_data(factors,f)
    p=size(r,2);k=size(f,2)
    allocate(fm3(k,k*k),fm4(k,k*k*k),fcov(k,k),ecov(p,p))
    call coskewness_matrix(f,fm3);call cokurtosis_matrix(f,fm4)
    call covariance_matrix(f,fcov,.false.);call covariance_matrix(resid,ecov,.false.)
    m3=0.0_dp;m4=0.0_dp
    do i=1,p;do j=1,p;do l=1,p
      col=(j-1)*p+l;v=0.0_dp
      do a=1,k;do b=1,k;do c=1,k
        v=v+beta(i,a)*beta(j,b)*beta(l,c)*fm3(a,(b-1)*k+c)
      end do;end do;end do
      if(i==j .and. j==l)v=v+sum(resid(:,i)**3)/real(size(r,1),dp)
      m3(i,col)=v
    end do;end do;end do
    do i=1,p;do j=1,p;do l=1,p;do h=1,p
      col=((j-1)*p+l-1)*p+h;v=0.0_dp
      do a=1,k;do b=1,k;do c=1,k;do d=1,k
        v=v+beta(i,a)*beta(j,b)*beta(l,c)*beta(h,d)*fm4(a,((b-1)*k+c-1)*k+d)
      end do;end do;end do;end do
      v=v+ecov(i,j)*ecov(l,h)+ecov(i,l)*ecov(j,h)+ecov(i,h)*ecov(j,l)
      if(i==j .and. j==l .and. l==h) v=v+sum(resid(:,i)**4)/real(size(r,1),dp)-3.0_dp*ecov(i,i)**2
      m4(i,col)=v
    end do;end do;end do;end do
  end subroutine factor_model_moments

  subroutine structured_coskewness(r, kind, out, factor, unbiased_marg)
    real(dp),intent(in)::r(:,:)
    character(len=*),intent(in)::kind
    real(dp),intent(out)::out(:,:)
    real(dp),intent(in),optional::factor(:,:)
    logical,intent(in),optional::unbiased_marg
    real(dp),allocatable::x(:,:),sample(:,:),m4tmp(:,:),beta(:,:),f(:,:),resid(:,:),sd(:)
    real(dp)::avgskew,sk,base,skew_scale
    integer::p,i,j,k,col,n
    logical::unbiased
    character(len=:),allocatable::key
    call centered_data(r,x);p=size(r,2);n=size(r,1);allocate(sample(p,p*p));call coskewness_matrix(x,sample)
    unbiased=.false.;if(present(unbiased_marg))unbiased=unbiased_marg
    skew_scale=1.0_dp;if(unbiased .and. n>2)skew_scale=real(n*n,dp)/real((n-1)*(n-2),dp)
    key=trim(lower_string(kind));out=0.0_dp
    select case(key)
    case('indep','independent')
      do i=1,p;out(i,(i-1)*p+i)=skew_scale*sum(x(:,i)**3)/real(n,dp);end do
    case('indepid','independent_equal')
      avgskew=skew_scale*sum([(sum(x(:,i)**3)/real(n,dp),i=1,p)])/real(p,dp)
      do i=1,p;out(i,(i-1)*p+i)=avgskew;end do
    case('centralsymmetric','cs','zero')
      out=0.0_dp
    case('observedfactor','factor')
      if(present(factor))then
        allocate(m4tmp(p,p*p*p));call factor_model_moments(r,factor,out,m4tmp)
      else;out=sample;end if
    case('latent1factor','pca')
      call pca_factors(r,1,beta,f,resid);allocate(m4tmp(p,p*p*p));call factor_model_moments(r,f,out,m4tmp)
    case('simaan')
      allocate(sd(p));do i=1,p
        sk=sum(x(:,i)**3)/real(size(x,1),dp);sd(i)=sign(abs(sk)**(1.0_dp/3.0_dp),sk)
      end do
      do i=1,p;do j=1,p;do k=1,p;out(i,(j-1)*p+k)=sd(i)*sd(j)*sd(k);end do;end do;end do
    case('cc','constant_correlation')
      allocate(sd(p));do i=1,p;sd(i)=sqrt(max(sum(x(:,i)**2)/real(size(x,1),dp),tiny(1.0_dp)));end do
      base=0.0_dp
      do i=1,p;do j=1,p;do k=1,p
        if(.not.(i==j .and. j==k))base=base+sample(i,(j-1)*p+k)/(sd(i)*sd(j)*sd(k))
      end do;end do;end do
      if(p**3-p>0)base=base/real(p**3-p,dp)
      do i=1,p;do j=1,p;do k=1,p
        col=(j-1)*p+k
        if(i==j .and. j==k)then;out(i,col)=sample(i,col);else;out(i,col)=base*sd(i)*sd(j)*sd(k);end if
      end do;end do;end do
    case default
      out=sample
    end select
  end subroutine structured_coskewness

  subroutine structured_cokurtosis(r, kind, out, factor)
    real(dp),intent(in)::r(:,:)
    character(len=*),intent(in)::kind
    real(dp),intent(out)::out(:,:)
    real(dp),intent(in),optional::factor(:,:)
    real(dp),allocatable::x(:,:),sample(:,:),m3tmp(:,:),cov(:,:),beta(:,:),f(:,:),resid(:,:)
    real(dp)::avgvar,avgfour,v
    integer::p,i,j,k,l,col
    character(len=:),allocatable::key
    call centered_data(r,x);p=size(r,2);allocate(sample(p,p*p*p),cov(p,p));call cokurtosis_matrix(x,sample)
    call covariance_matrix(x,cov,.false.);key=trim(lower_string(kind));out=0.0_dp
    select case(key)
    case('indep','independent')
      do i=1,p;do j=1,p;do k=1,p;do l=1,p
        col=((j-1)*p+k-1)*p+l
        v=cov(i,j)*cov(k,l)+cov(i,k)*cov(j,l)+cov(i,l)*cov(j,k)
        if(i==j .and. j==k .and. k==l)v=sum(x(:,i)**4)/real(size(x,1),dp)
        out(i,col)=v
      end do;end do;end do;end do
    case('indepid','independent_equal')
      avgvar=sum([(cov(i,i),i=1,p)])/real(p,dp)
      avgfour=sum([(sum(x(:,i)**4)/real(size(x,1),dp),i=1,p)])/real(p,dp)
      do i=1,p;do j=1,p;do k=1,p;do l=1,p
        col=((j-1)*p+k-1)*p+l;v=0.0_dp
        if(i==j .and. k==l)v=v+avgvar*avgvar
        if(i==k .and. j==l)v=v+avgvar*avgvar
        if(i==l .and. j==k)v=v+avgvar*avgvar
        if(i==j .and. j==k .and. k==l)v=avgfour
        out(i,col)=v
      end do;end do;end do;end do
    case('observedfactor','factor')
      if(present(factor))then
        allocate(m3tmp(p,p*p));call factor_model_moments(r,factor,m3tmp,out)
      else;out=sample;end if
    case('latent1factor','pca')
      call pca_factors(r,1,beta,f,resid);allocate(m3tmp(p,p*p));call factor_model_moments(r,f,m3tmp,out)
    case('cc','constant_correlation')
      call constant_correlation_m4(sample,cov,out)
    case default
      out=sample
    end select
  end subroutine structured_cokurtosis


  subroutine constant_correlation_m4(sample,cov,out)
    real(dp),intent(in)::sample(:,:),cov(:,:)
    real(dp),intent(out)::out(:,:)
    real(dp),allocatable::sd(:),sumcat(:)
    integer,allocatable::ncat(:)
    integer::p,i,j,k,l,col,cat
    real(dp)::den,coef
    p=size(cov,1);allocate(sd(p),sumcat(5),ncat(5));sumcat=0.0_dp;ncat=0
    do i=1,p;sd(i)=sqrt(max(cov(i,i),tiny(1.0_dp)));end do
    do i=1,p;do j=i,p;do k=j,p;do l=k,p
      if(i==l)cycle
      cat=multiplicity_category4(i,j,k,l);den=sd(i)*sd(j)*sd(k)*sd(l)
      if(den>0.0_dp)then;sumcat(cat)=sumcat(cat)+sample(i,((j-1)*p+k-1)*p+l)/den;ncat(cat)=ncat(cat)+1;end if
    end do;end do;end do;end do
    out=0.0_dp
    do i=1,p;do j=1,p;do k=1,p;do l=1,p
      col=((j-1)*p+k-1)*p+l
      if(i==j .and. j==k .and. k==l)then
        out(i,col)=sample(i,col)
      else
        cat=multiplicity_category4(i,j,k,l)
        if(ncat(cat)>0)then;coef=sumcat(cat)/real(ncat(cat),dp);else;coef=0.0_dp;end if
        out(i,col)=coef*sd(i)*sd(j)*sd(k)*sd(l)
      end if
    end do;end do;end do;end do
  end subroutine constant_correlation_m4

  pure integer function multiplicity_category4(i,j,k,l) result(cat)
    integer,intent(in)::i,j,k,l
    integer::x(4),a,b,counts(4),nuniq,maxcount
    x=[i,j,k,l];counts=0;nuniq=0
    do a=1,4
      do b=1,a-1
        if(x(a)==x(b))exit
      end do
      if(b==a)then
        nuniq=nuniq+1;counts(nuniq)=count(x==x(a))
      end if
    end do
    maxcount=maxval(counts)
    if(maxcount==3)then;cat=1
    else if(maxcount==2 .and. nuniq==2)then;cat=2
    else if(maxcount==2)then;cat=3
    else if(nuniq==4)then;cat=4
    else;cat=5
    end if
  end function multiplicity_category4

  subroutine shrink_covariance(r, target_kind, out, lambda, factor)
    real(dp),intent(in)::r(:,:)
    character(len=*),intent(in)::target_kind
    real(dp),intent(out)::out(:,:)
    real(dp),intent(out)::lambda
    real(dp),intent(in),optional::factor(:,:)
    real(dp),allocatable::x(:,:),s(:,:),targ(:,:),obs(:,:)
    real(dp)::num,den
    integer::n,p,i
    call centered_data(r,x);n=size(x,1);p=size(x,2);allocate(s(p,p),targ(p,p),obs(p,p))
    call covariance_matrix(x,s,.false.)
    if(present(factor))then;call structured_covariance(r,target_kind,targ,factor);else;call structured_covariance(r,target_kind,targ);end if
    num=0.0_dp
    do i=1,n
      obs=spread(x(i,:),2,p)*spread(x(i,:),1,p)
      num=num+sum((obs-s)**2)
    end do
    num=num/real(max(n*n,1),dp);den=sum((s-targ)**2)
    if(den<=tiny(1.0_dp))then;lambda=0.0_dp;else;lambda=max(0.0_dp,min(1.0_dp,num/den));end if
    out=(1.0_dp-lambda)*s+lambda*targ
  end subroutine shrink_covariance

  subroutine shrink_coskewness(r, target_kind, out, lambda, factor)
    real(dp),intent(in)::r(:,:)
    character(len=*),intent(in)::target_kind
    real(dp),intent(out)::out(:,:)
    real(dp),intent(out)::lambda
    real(dp),intent(in),optional::factor(:,:)
    real(dp),allocatable::x(:,:),s(:,:),targ(:,:),obs(:,:)
    real(dp)::num,den
    integer::n,p,i,j,k,l,col
    call centered_data(r,x);n=size(x,1);p=size(x,2);allocate(s(p,p*p),targ(p,p*p),obs(p,p*p))
    call coskewness_matrix(x,s)
    if(present(factor))then;call structured_coskewness(r,target_kind,targ,factor);else;call structured_coskewness(r,target_kind,targ);end if
    num=0.0_dp
    do l=1,n
      do i=1,p;do j=1,p;do k=1,p;col=(j-1)*p+k;obs(i,col)=x(l,i)*x(l,j)*x(l,k);end do;end do;end do
      num=num+sum((obs-s)**2)
    end do
    num=num/real(max(n*n,1),dp);den=sum((s-targ)**2)
    if(den<=tiny(1.0_dp))then;lambda=0.0_dp;else;lambda=max(0.0_dp,min(1.0_dp,num/den));end if
    out=(1.0_dp-lambda)*s+lambda*targ
  end subroutine shrink_coskewness

  subroutine shrink_cokurtosis(r, target_kind, out, lambda, factor)
    real(dp),intent(in)::r(:,:)
    character(len=*),intent(in)::target_kind
    real(dp),intent(out)::out(:,:)
    real(dp),intent(out)::lambda
    real(dp),intent(in),optional::factor(:,:)
    real(dp),allocatable::x(:,:),s(:,:),targ(:,:),obs(:,:)
    real(dp)::num,den
    integer::n,p,t,i,j,k,l,col
    call centered_data(r,x);n=size(x,1);p=size(x,2);allocate(s(p,p*p*p),targ(p,p*p*p),obs(p,p*p*p))
    call cokurtosis_matrix(x,s)
    if(present(factor))then;call structured_cokurtosis(r,target_kind,targ,factor);else;call structured_cokurtosis(r,target_kind,targ);end if
    num=0.0_dp
    do t=1,n
      do i=1,p;do j=1,p;do k=1,p;do l=1,p
        col=((j-1)*p+k-1)*p+l;obs(i,col)=x(t,i)*x(t,j)*x(t,k)*x(t,l)
      end do;end do;end do;end do
      num=num+sum((obs-s)**2)
    end do
    num=num/real(max(n*n,1),dp);den=sum((s-targ)**2)
    if(den<=tiny(1.0_dp))then;lambda=0.0_dp;else;lambda=max(0.0_dp,min(1.0_dp,num/den));end if
    out=(1.0_dp-lambda)*s+lambda*targ
  end subroutine shrink_cokurtosis


  subroutine project_simplex_leq_one(x)
    real(dp),intent(inout)::x(:)
    real(dp)::lo,hi,mid,s
    integer::it
    x=max(x,0.0_dp)
    if(sum(x)<=1.0_dp)return
    lo=minval(x)-1.0_dp;hi=maxval(x)
    do it=1,80
      mid=0.5_dp*(lo+hi);s=sum(max(x-mid,0.0_dp))
      if(s>1.0_dp)then;lo=mid;else;hi=mid;end if
    end do
    x=max(x-hi,0.0_dp)
  end subroutine project_simplex_leq_one

  subroutine solve_shrinkage_weights(dmat,variance_term,lambda)
    real(dp),intent(in)::dmat(:,:),variance_term
    real(dp),intent(out)::lambda(:)
    real(dp),allocatable::a(:,:),grad(:),old(:)
    real(dp)::step
    integer::m,it
    m=size(dmat,2);allocate(a(m,m),grad(m),old(m));a=matmul(transpose(dmat),dmat)
    step=1.0_dp/max(maxval(sum(abs(a),dim=2)),1.0e-12_dp);lambda=0.0_dp
    do it=1,1000
      old=lambda;grad=matmul(a,lambda)-variance_term
      lambda=lambda-step*grad;call project_simplex_leq_one(lambda)
      if(maxval(abs(lambda-old))<1.0e-10_dp)exit
    end do
  end subroutine solve_shrinkage_weights

  subroutine multi_target_shrink_covariance(r,target_kinds,out,lambdas,factor)
    real(dp),intent(in)::r(:,:)
    character(len=*),intent(in)::target_kinds(:)
    real(dp),intent(out)::out(:,:),lambdas(:)
    real(dp),intent(in),optional::factor(:,:)
    real(dp),allocatable::x(:,:),sample(:,:),target(:,:),dmat(:,:),obs(:,:)
    real(dp)::varterm
    integer::n,p,m,t,j
    call centered_data(r,x);n=size(x,1);p=size(x,2);m=size(target_kinds)
    allocate(sample(p,p),target(p,p),dmat(p*p,m),obs(p,p));call covariance_matrix(x,sample,.false.);varterm=0.0_dp
    do t=1,n;obs=spread(x(t,:),2,p)*spread(x(t,:),1,p);varterm=varterm+sum((obs-sample)**2);end do
    varterm=varterm/real(max(n*n,1),dp)
    do j=1,m
      if(present(factor))then;call structured_covariance(r,target_kinds(j),target,factor);else;call structured_covariance(r,target_kinds(j),target);end if
      dmat(:,j)=reshape(target-sample,[p*p])
    end do
    call solve_shrinkage_weights(dmat,varterm,lambdas(:m));out=sample
    do j=1,m;out=out+lambdas(j)*reshape(dmat(:,j),[p,p]);end do
  end subroutine multi_target_shrink_covariance

  subroutine multi_target_shrink_coskewness(r,target_kinds,out,lambdas,factor)
    real(dp),intent(in)::r(:,:)
    character(len=*),intent(in)::target_kinds(:)
    real(dp),intent(out)::out(:,:),lambdas(:)
    real(dp),intent(in),optional::factor(:,:)
    real(dp),allocatable::x(:,:),sample(:,:),target(:,:),dmat(:,:),obs(:,:)
    real(dp)::varterm
    integer::n,p,m,t,j,i,a,b,col
    call centered_data(r,x);n=size(x,1);p=size(x,2);m=size(target_kinds)
    allocate(sample(p,p*p),target(p,p*p),dmat(p*p*p,m),obs(p,p*p));call coskewness_matrix(x,sample);varterm=0.0_dp
    do t=1,n
      do i=1,p;do a=1,p;do b=1,p;col=(a-1)*p+b;obs(i,col)=x(t,i)*x(t,a)*x(t,b);end do;end do;end do
      varterm=varterm+sum((obs-sample)**2)
    end do
    varterm=varterm/real(max(n*n,1),dp)
    do j=1,m
      if(present(factor))then;call structured_coskewness(r,target_kinds(j),target,factor);else;call structured_coskewness(r,target_kinds(j),target);end if
      dmat(:,j)=reshape(target-sample,[p*p*p])
    end do
    call solve_shrinkage_weights(dmat,varterm,lambdas(:m));out=sample
    do j=1,m;out=out+lambdas(j)*reshape(dmat(:,j),[p,p*p]);end do
  end subroutine multi_target_shrink_coskewness

  subroutine multi_target_shrink_cokurtosis(r,target_kinds,out,lambdas,factor)
    real(dp),intent(in)::r(:,:)
    character(len=*),intent(in)::target_kinds(:)
    real(dp),intent(out)::out(:,:),lambdas(:)
    real(dp),intent(in),optional::factor(:,:)
    real(dp),allocatable::x(:,:),sample(:,:),target(:,:),dmat(:,:),obs(:,:)
    real(dp)::varterm
    integer::n,p,m,t,j,i,a,b,c,col
    call centered_data(r,x);n=size(x,1);p=size(x,2);m=size(target_kinds)
    allocate(sample(p,p*p*p),target(p,p*p*p),dmat(p*p*p*p,m),obs(p,p*p*p));call cokurtosis_matrix(x,sample);varterm=0.0_dp
    do t=1,n
      do i=1,p;do a=1,p;do b=1,p;do c=1,p;col=((a-1)*p+b-1)*p+c;obs(i,col)=x(t,i)*x(t,a)*x(t,b)*x(t,c);end do;end do;end do;end do
      varterm=varterm+sum((obs-sample)**2)
    end do
    varterm=varterm/real(max(n*n,1),dp)
    do j=1,m
      if(present(factor))then;call structured_cokurtosis(r,target_kinds(j),target,factor);else;call structured_cokurtosis(r,target_kinds(j),target);end if
      dmat(:,j)=reshape(target-sample,[p*p*p*p])
    end do
    call solve_shrinkage_weights(dmat,varterm,lambdas(:m));out=sample
    do j=1,m;out=out+lambdas(j)*reshape(dmat(:,j),[p,p*p*p]);end do
  end subroutine multi_target_shrink_cokurtosis

  subroutine ewma_coskewness(r,lambda,out,last)
    real(dp),intent(in)::r(:,:),lambda
    real(dp),intent(out)::out(:,:)
    real(dp),intent(in),optional::last(:,:)
    real(dp),allocatable::mu(:),d(:)
    integer::p,t,i,j,k,col
    p=size(r,2);allocate(mu(p),d(p));mu=0.0_dp
    if(present(last))then;out=last;else;out=0.0_dp;end if
    do t=1,size(r,1)
      d=r(t,:)-mu
      do i=1,p;do j=1,p;do k=1,p;col=(j-1)*p+k
        out(i,col)=lambda*out(i,col)+(1.0_dp-lambda)*d(i)*d(j)*d(k)
      end do;end do;end do
      mu=lambda*mu+(1.0_dp-lambda)*r(t,:)
    end do
  end subroutine ewma_coskewness

  subroutine ewma_cokurtosis(r,lambda,out,last)
    real(dp),intent(in)::r(:,:),lambda
    real(dp),intent(out)::out(:,:)
    real(dp),intent(in),optional::last(:,:)
    real(dp),allocatable::mu(:),d(:)
    integer::p,t,i,j,k,l,col
    p=size(r,2);allocate(mu(p),d(p));mu=0.0_dp
    if(present(last))then;out=last;else;out=0.0_dp;end if
    do t=1,size(r,1)
      d=r(t,:)-mu
      do i=1,p;do j=1,p;do k=1,p;do l=1,p;col=((j-1)*p+k-1)*p+l
        out(i,col)=lambda*out(i,col)+(1.0_dp-lambda)*d(i)*d(j)*d(k)*d(l)
      end do;end do;end do;end do
      mu=lambda*mu+(1.0_dp-lambda)*r(t,:)
    end do
  end subroutine ewma_cokurtosis

  subroutine leading_eigenvectors(g,k,u)
    real(dp),intent(in)::g(:,:)
    integer,intent(in)::k
    real(dp),intent(out)::u(:,:)
    real(dp),allocatable::eval(:),evec(:,:)
    logical::ok
    call symmetric_eigen_jacobi(g,eval,evec,ok)
    if(ok)then;u=evec(:,1:k);else;u=0.0_dp;end if
  end subroutine leading_eigenvectors

  subroutine m3_gram(m3,u,g)
    real(dp),intent(in)::m3(:,:),u(:,:)
    real(dp),intent(out)::g(:,:)
    real(dp),allocatable::z(:,:)
    integer::p,k,i,j,l,a,b,col
    p=size(u,1);k=size(u,2);allocate(z(p,k*k));z=0.0_dp
    do i=1,p;do a=1,k;do b=1,k
      col=(a-1)*k+b
      do j=1,p;do l=1,p;z(i,col)=z(i,col)+m3(i,(j-1)*p+l)*u(j,a)*u(l,b);end do;end do
    end do;end do;end do
    g=matmul(z,transpose(z))
  end subroutine m3_gram

  subroutine m4_gram(m4,u,g)
    real(dp),intent(in)::m4(:,:),u(:,:)
    real(dp),intent(out)::g(:,:)
    real(dp),allocatable::z(:,:)
    integer::p,k,i,j,l,h,a,b,c,col,q
    p=size(u,1);k=size(u,2);allocate(z(p,k*k*k));z=0.0_dp
    do i=1,p;do a=1,k;do b=1,k;do c=1,k
      q=((a-1)*k+b-1)*k+c
      do j=1,p;do l=1,p;do h=1,p
        col=((j-1)*p+l-1)*p+h;z(i,q)=z(i,q)+m4(i,col)*u(j,a)*u(l,b)*u(h,c)
      end do;end do;end do
    end do;end do;end do;end do
    g=matmul(z,transpose(z))
  end subroutine m4_gram

  subroutine reconstruct_m3(m3,u,out)
    real(dp),intent(in)::m3(:,:),u(:,:)
    real(dp),intent(out)::out(:,:)
    real(dp),allocatable::core(:,:,:)
    real(dp)::v
    integer::p,k,i,j,l,a,b,c,ii,jj,ll,col
    p=size(u,1);k=size(u,2);allocate(core(k,k,k));core=0.0_dp
    do a=1,k;do b=1,k;do c=1,k
      do ii=1,p;do jj=1,p;do ll=1,p
        core(a,b,c)=core(a,b,c)+u(ii,a)*u(jj,b)*u(ll,c)*m3(ii,(jj-1)*p+ll)
      end do;end do;end do
    end do;end do;end do
    do i=1,p;do j=1,p;do l=1,p;v=0.0_dp
      do a=1,k;do b=1,k;do c=1,k;v=v+u(i,a)*u(j,b)*u(l,c)*core(a,b,c);end do;end do;end do
      col=(j-1)*p+l;out(i,col)=v
    end do;end do;end do
  end subroutine reconstruct_m3

  subroutine reconstruct_m4(m4,u,out)
    real(dp),intent(in)::m4(:,:),u(:,:)
    real(dp),intent(out)::out(:,:)
    real(dp),allocatable::core(:,:,:,:)
    real(dp)::v
    integer::p,k,i,j,l,h,a,b,c,d,ii,jj,ll,hh,col
    p=size(u,1);k=size(u,2);allocate(core(k,k,k,k));core=0.0_dp
    do a=1,k;do b=1,k;do c=1,k;do d=1,k
      do ii=1,p;do jj=1,p;do ll=1,p;do hh=1,p
        core(a,b,c,d)=core(a,b,c,d)+u(ii,a)*u(jj,b)*u(ll,c)*u(hh,d)*m4(ii,((jj-1)*p+ll-1)*p+hh)
      end do;end do;end do;end do
    end do;end do;end do;end do
    do i=1,p;do j=1,p;do l=1,p;do h=1,p;v=0.0_dp
      do a=1,k;do b=1,k;do c=1,k;do d=1,k
        v=v+u(i,a)*u(j,b)*u(l,c)*u(h,d)*core(a,b,c,d)
      end do;end do;end do;end do
      col=((j-1)*p+l-1)*p+h;out(i,col)=v
    end do;end do;end do;end do
  end subroutine reconstruct_m4

  subroutine m3_mca(r,k,result,max_iterations,tolerance,m3_input)
    real(dp),intent(in)::r(:,:)
    integer,intent(in)::k
    type(mca_result),intent(out)::result
    integer,intent(in),optional::max_iterations
    real(dp),intent(in),optional::tolerance,m3_input(:,:)
    real(dp),allocatable::m3(:,:),cov(:,:),u(:,:),unew(:,:),g(:,:)
    real(dp)::tol,diff
    integer::p,kk,it,maxit
    p=size(r,2);kk=max(1,min(k,p));maxit=200;if(present(max_iterations))maxit=max_iterations
    tol=1.0e-6_dp;if(present(tolerance))tol=tolerance
    allocate(m3(p,p*p),cov(p,p),u(p,kk),unew(p,kk),g(p,p))
    if(present(m3_input))then;m3=m3_input;else;call coskewness_matrix(r,m3);end if
    call covariance_matrix(r,cov,.false.);call leading_eigenvectors(cov,kk,u)
    result%converged=.false.
    do it=1,maxit
      call m3_gram(m3,u,g);call leading_eigenvectors(g,kk,unew)
      diff=sqrt(sum((abs(unew)-abs(u))**2));u=unew
      if(diff<tol)then;result%converged=.true.;exit;end if
    end do
    result%components=kk;result%iterations=min(it,maxit)
    allocate(result%directions(p,kk),result%m3(p,p*p),result%m4(0,0));result%directions=u
    call reconstruct_m3(m3,u,result%m3)
  end subroutine m3_mca

  subroutine m4_mca(r,k,result,max_iterations,tolerance,m4_input)
    real(dp),intent(in)::r(:,:)
    integer,intent(in)::k
    type(mca_result),intent(out)::result
    integer,intent(in),optional::max_iterations
    real(dp),intent(in),optional::tolerance,m4_input(:,:)
    real(dp),allocatable::m4(:,:),cov(:,:),u(:,:),unew(:,:),g(:,:)
    real(dp)::tol,diff
    integer::p,kk,it,maxit
    p=size(r,2);kk=max(1,min(k,p));maxit=200;if(present(max_iterations))maxit=max_iterations
    tol=1.0e-6_dp;if(present(tolerance))tol=tolerance
    allocate(m4(p,p*p*p),cov(p,p),u(p,kk),unew(p,kk),g(p,p))
    if(present(m4_input))then;m4=m4_input;else;call cokurtosis_matrix(r,m4);end if
    call covariance_matrix(r,cov,.false.);call leading_eigenvectors(cov,kk,u)
    result%converged=.false.
    do it=1,maxit
      call m4_gram(m4,u,g);call leading_eigenvectors(g,kk,unew)
      diff=sqrt(sum((abs(unew)-abs(u))**2));u=unew
      if(diff<tol)then;result%converged=.true.;exit;end if
    end do
    result%components=kk;result%iterations=min(it,maxit)
    allocate(result%directions(p,kk),result%m4(p,p*p*p),result%m3(0,0));result%directions=u
    call reconstruct_m4(m4,u,result%m4)
  end subroutine m4_mca

  subroutine nce_model_moments(b,d,fs,es,fk,ek,m2,m3,m4)
    real(dp),intent(in)::b(:,:),d(:),fs(:),es(:),fk(:),ek(:)
    real(dp),intent(out)::m2(:,:),m3(:,:),m4(:,:)
    real(dp),allocatable::s(:,:)
    real(dp)::v
    integer::p,k,i,j,l,h,a,col
    p=size(d);k=size(b,2);allocate(s(p,p));s=matmul(b,transpose(b));m2=s
    do i=1,p;m2(i,i)=m2(i,i)+d(i);end do
    m3=0.0_dp
    do i=1,p;do j=1,p;do l=1,p;v=0.0_dp
      do a=1,k;v=v+fs(a)*b(i,a)*b(j,a)*b(l,a);end do
      if(i==j .and. j==l)v=v+es(i)
      m3(i,(j-1)*p+l)=v
    end do;end do;end do
    do i=1,p;do j=1,p;do l=1,p;do h=1,p
      col=((j-1)*p+l-1)*p+h
      v=m2(i,j)*m2(l,h)+m2(i,l)*m2(j,h)+m2(i,h)*m2(j,l)
      do a=1,k;v=v+(fk(a)-3.0_dp)*b(i,a)*b(j,a)*b(l,a)*b(h,a);end do
      if(i==j .and. j==l .and. l==h)v=v+ek(i)-3.0_dp*d(i)*d(i)
      m4(i,col)=v
    end do;end do;end do;end do
  end subroutine nce_model_moments

  real(dp) function nce_objective(b,d,fs,es,fk,ek,s2,s3,s4,include,w2,w3,w4) result(obj)
    real(dp),intent(in)::b(:,:),d(:),fs(:),es(:),fk(:),ek(:),s2(:,:),s3(:,:),s4(:,:)
    logical,intent(in)::include(3)
    real(dp),intent(in)::w2,w3,w4
    real(dp),allocatable::m2(:,:),m3(:,:),m4(:,:)
    real(dp)::den
    integer::p
    p=size(d);allocate(m2(p,p),m3(p,p*p),m4(p,p*p*p));call nce_model_moments(b,d,fs,es,fk,ek,m2,m3,m4)
    obj=0.0_dp
    if(include(1))then;den=max(sum(s2*s2),1.0e-20_dp);obj=obj+w2*sum((s2-m2)**2)/den;end if
    if(include(2))then;den=max(sum(s3*s3),1.0e-20_dp);obj=obj+w3*sum((s3-m3)**2)/den;end if
    if(include(3))then;den=max(sum(s4*s4),1.0e-20_dp);obj=obj+w4*sum((s4-m4)**2)/den;end if
  end function nce_objective

  subroutine project_nce(d,fs,es,fk,ek)
    real(dp),intent(inout)::d(:),fs(:),es(:),fk(:),ek(:)
    integer::i
    d=max(d,1.0e-10_dp)
    do i=1,size(fk);fk(i)=max(fk(i),fs(i)*fs(i)+1.0_dp+1.0e-8_dp);end do
    do i=1,size(ek);ek(i)=max(ek(i),d(i)*d(i)+es(i)*es(i)/d(i)+1.0e-12_dp);end do
  end subroutine project_nce


  subroutine nce_order_variances(x,s2,s3,s4,v)
    real(dp),intent(in)::x(:,:),s2(:,:),s3(:,:),s4(:,:)
    real(dp),intent(out)::v(3)
    real(dp),allocatable::o2(:,:),o3(:,:),o4(:,:)
    integer::n,p,t,i,j,k,l,col
    n=size(x,1);p=size(x,2);allocate(o2(p,p),o3(p,p*p),o4(p,p*p*p));v=0.0_dp
    do t=1,n
      o2=spread(x(t,:),2,p)*spread(x(t,:),1,p)
      do i=1,p;do j=1,p;do k=1,p;col=(j-1)*p+k;o3(i,col)=x(t,i)*x(t,j)*x(t,k);end do;end do;end do
      do i=1,p;do j=1,p;do k=1,p;do l=1,p;col=((j-1)*p+k-1)*p+l;o4(i,col)=x(t,i)*x(t,j)*x(t,k)*x(t,l);end do;end do;end do;end do
      v(1)=v(1)+sum((o2-s2)**2)/real(size(s2),dp)
      v(2)=v(2)+sum((o3-s3)**2)/real(size(s3),dp)
      v(3)=v(3)+sum((o4-s4)**2)/real(size(s4),dp)
    end do
    v=v/real(max(n-1,1),dp)
  end subroutine nce_order_variances

  subroutine nearest_comoment_estimator(r,k,result,include_orders,max_iterations,tolerance,weights,weight_mode,ridge_alpha)
    real(dp),intent(in)::r(:,:)
    integer,intent(in)::k
    type(nce_result),intent(out)::result
    logical,intent(in),optional::include_orders(3)
    integer,intent(in),optional::max_iterations
    real(dp),intent(in),optional::tolerance,weights(3),ridge_alpha
    character(len=*),intent(in),optional::weight_mode
    real(dp),allocatable::x(:,:),s2(:,:),s3(:,:),s4(:,:),b(:,:),f(:,:),resid(:,:)
    real(dp),allocatable::d(:),fs(:),es(:),fk(:),ek(:),steps(:)
    real(dp)::obj,best,old,trial,tol,w(3),scale,ov(3),alpha_w,avg_v
    character(len=24)::wm
    logical::inc(3),improved
    integer::p,kk,n,it,maxit,i,j
    call centered_data(r,x);n=size(r,1);p=size(r,2);kk=max(0,min(k,p))
    inc=.true.;if(present(include_orders))inc=include_orders
    maxit=120;if(present(max_iterations))maxit=max_iterations
    tol=1.0e-7_dp;if(present(tolerance))tol=tolerance
    w=1.0_dp
    allocate(s2(p,p),s3(p,p*p),s4(p,p*p*p));call covariance_matrix(x,s2,.false.)
    call coskewness_matrix(x,s3);call cokurtosis_matrix(x,s4)
    if(present(weights))then
      w=weights
    else if(present(weight_mode))then
      wm=trim(lower_string(weight_mode));call nce_order_variances(x,s2,s3,s4,ov)
      select case(trim(wm))
      case('diagonal','d')
        w=1.0_dp/max(ov,1.0e-16_dp)
      case('ridge_diagonal','ridged')
        alpha_w=0.1_dp;if(present(ridge_alpha))alpha_w=max(0.0_dp,min(1.0_dp,ridge_alpha))
        avg_v=sum(ov)/3.0_dp;w=1.0_dp/max((1.0_dp-alpha_w)*ov+alpha_w*avg_v,1.0e-16_dp)
      case('ridge_identity','ridgei')
        alpha_w=0.1_dp;if(present(ridge_alpha))alpha_w=max(0.0_dp,min(1.0_dp,ridge_alpha))
        avg_v=sum(ov)/3.0_dp;w=1.0_dp/max((1.0_dp-alpha_w)*ov+alpha_w*avg_v,1.0e-16_dp)
      case default
        w=1.0_dp
      end select
      if(sum(w)>0.0_dp)w=3.0_dp*w/sum(w)
    end if
    call pca_factors(x,kk,b,f,resid)
    allocate(d(p),fs(kk),es(p),fk(kk),ek(p),steps(6))
    do i=1,p
      d(i)=max(sum(resid(:,i)**2)/real(n,dp),1.0e-8_dp)
      es(i)=sum(resid(:,i)**3)/real(n,dp)
      ek(i)=sum(resid(:,i)**4)/real(n,dp)
    end do
    do j=1,kk
      fs(j)=sum((f(:,j)-mean_value(f(:,j)))**3)/real(n,dp)
      fk(j)=sum((f(:,j)-mean_value(f(:,j)))**4)/real(n,dp)
    end do
    call project_nce(d,fs,es,fk,ek)
    scale=sqrt(max(sum([(s2(i,i),i=1,p)])/real(max(p,1),dp),1.0e-8_dp))
    steps=[0.08_dp*scale,0.08_dp*scale*scale,0.15_dp,0.08_dp*scale**3,0.25_dp,0.08_dp*scale**4]
    obj=nce_objective(b,d,fs,es,fk,ek,s2,s3,s4,inc,w(1),w(2),w(3))
    result%converged=.false.
    do it=1,maxit
      improved=.false.
      do i=1,p;do j=1,kk
        old=b(i,j);best=obj
        b(i,j)=old+steps(1);trial=nce_objective(b,d,fs,es,fk,ek,s2,s3,s4,inc,w(1),w(2),w(3))
        if(trial<best)then;best=trial;else;b(i,j)=old-steps(1);trial=nce_objective(b,d,fs,es,fk,ek,s2,s3,s4,inc,w(1),w(2),w(3));if(trial<best)then;best=trial;else;b(i,j)=old;end if;end if
        if(best<obj)then;obj=best;improved=.true.;end if
      end do;end do
      do i=1,p
        old=d(i);best=obj;d(i)=max(1.0e-10_dp,old+steps(2));call project_nce(d,fs,es,fk,ek);trial=nce_objective(b,d,fs,es,fk,ek,s2,s3,s4,inc,w(1),w(2),w(3))
        if(trial<best)then;best=trial;else;d(i)=max(1.0e-10_dp,old-steps(2));call project_nce(d,fs,es,fk,ek);trial=nce_objective(b,d,fs,es,fk,ek,s2,s3,s4,inc,w(1),w(2),w(3));if(trial<best)then;best=trial;else;d(i)=old;call project_nce(d,fs,es,fk,ek);end if;end if
        if(best<obj)then;obj=best;improved=.true.;end if
        old=es(i);best=obj;es(i)=old+steps(4);call project_nce(d,fs,es,fk,ek);trial=nce_objective(b,d,fs,es,fk,ek,s2,s3,s4,inc,w(1),w(2),w(3))
        if(trial<best)then;best=trial;else;es(i)=old-steps(4);call project_nce(d,fs,es,fk,ek);trial=nce_objective(b,d,fs,es,fk,ek,s2,s3,s4,inc,w(1),w(2),w(3));if(trial<best)then;best=trial;else;es(i)=old;call project_nce(d,fs,es,fk,ek);end if;end if
        if(best<obj)then;obj=best;improved=.true.;end if
        old=ek(i);best=obj;ek(i)=old+steps(6);call project_nce(d,fs,es,fk,ek);trial=nce_objective(b,d,fs,es,fk,ek,s2,s3,s4,inc,w(1),w(2),w(3))
        if(trial<best)then;best=trial;else;ek(i)=max(d(i)*d(i)+es(i)*es(i)/d(i)+1.0e-12_dp,old-steps(6));trial=nce_objective(b,d,fs,es,fk,ek,s2,s3,s4,inc,w(1),w(2),w(3));if(trial<best)then;best=trial;else;ek(i)=old;end if;end if
        if(best<obj)then;obj=best;improved=.true.;end if
      end do
      do j=1,kk
        old=fs(j);best=obj;fs(j)=old+steps(3);call project_nce(d,fs,es,fk,ek);trial=nce_objective(b,d,fs,es,fk,ek,s2,s3,s4,inc,w(1),w(2),w(3))
        if(trial<best)then;best=trial;else;fs(j)=old-steps(3);call project_nce(d,fs,es,fk,ek);trial=nce_objective(b,d,fs,es,fk,ek,s2,s3,s4,inc,w(1),w(2),w(3));if(trial<best)then;best=trial;else;fs(j)=old;call project_nce(d,fs,es,fk,ek);end if;end if
        if(best<obj)then;obj=best;improved=.true.;end if
        old=fk(j);best=obj;fk(j)=old+steps(5);call project_nce(d,fs,es,fk,ek);trial=nce_objective(b,d,fs,es,fk,ek,s2,s3,s4,inc,w(1),w(2),w(3))
        if(trial<best)then;best=trial;else;fk(j)=max(fs(j)*fs(j)+1.0_dp+1.0e-8_dp,old-steps(5));trial=nce_objective(b,d,fs,es,fk,ek,s2,s3,s4,inc,w(1),w(2),w(3));if(trial<best)then;best=trial;else;fk(j)=old;end if;end if
        if(best<obj)then;obj=best;improved=.true.;end if
      end do
      if(.not.improved)steps=0.5_dp*steps
      if(maxval(steps)<tol)then;result%converged=.true.;exit;end if
    end do
    result%factors=kk;result%iterations=min(it,maxit);result%objective=obj
    allocate(result%loadings(p,kk),result%residual_variance(p),result%factor_skewness(kk),result%residual_third(p), &
      result%factor_kurtosis(kk),result%residual_fourth(p),result%covariance(p,p),result%coskewness(p,p*p),result%cokurtosis(p,p*p*p))
    result%loadings=b;result%residual_variance=d;result%factor_skewness=fs;result%residual_third=es
    result%factor_kurtosis=fk;result%residual_fourth=ek
    call nce_model_moments(b,d,fs,es,fk,ek,result%covariance,result%coskewness,result%cokurtosis)
  end subroutine nearest_comoment_estimator
end module advanced_moments_mod
