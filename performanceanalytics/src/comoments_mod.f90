! SPDX-License-Identifier: GPL-2.0-or-later
! This file is part of a modern Fortran translation of PerformanceAnalytics.
! Copyright (C) original PerformanceAnalytics authors and translation contributors.
! Distributed under GNU GPL version 2 or, at your option, any later version.
module comoments_mod
  use kinds_mod, only: dp
  use statistics_mod, only: mean_value
  implicit none
  private
  public :: covariance_matrix, coskewness_matrix, cokurtosis_matrix
  public :: coskewness_unique, cokurtosis_unique, m3_vec_to_mat, m4_vec_to_mat
  public :: m3_mat_to_vec, m4_mat_to_vec, m3_inner_product, m4_inner_product
  public :: portfolio_mean, portfolio_variance, portfolio_third_moment, portfolio_fourth_moment
  public :: portfolio_skewness, portfolio_kurtosis, ewma_covariance
  public :: covariation, coskewness_value, cokurtosis_value, beta_coskewness, beta_cokurtosis
contains
  subroutine covariance_matrix(r,cov,unbiased)
    real(dp),intent(in)::r(:,:)
    real(dp),intent(out)::cov(:,:)
    logical,intent(in),optional::unbiased
    real(dp),allocatable::xc(:,:)
    real(dp)::den
    integer::n,p,j
    logical::u
    n=size(r,1);p=size(r,2);u=.true.;if(present(unbiased))u=unbiased
    allocate(xc(n,p))
    do j=1,p;xc(:,j)=r(:,j)-mean_value(r(:,j));end do
    den=real(n-merge(1,0,u),dp)
    if(den<=0.0_dp)then;cov=0.0_dp;else;cov=matmul(transpose(xc),xc)/den;end if
  end subroutine covariance_matrix

  pure real(dp) function covariation(x,y) result(v)
    real(dp),intent(in)::x(:),y(:)
    integer::n
    n=min(size(x),size(y))
    if(n==0)then;v=0.0_dp;else;v=sum((x(:n)-mean_value(x(:n)))*(y(:n)-mean_value(y(:n))))/real(n,dp);end if
  end function covariation

  pure real(dp) function coskewness_value(x,y,p1,p2,normalize) result(v)
    real(dp),intent(in)::x(:),y(:)
    integer,intent(in)::p1,p2
    logical,intent(in),optional::normalize
    integer::n
    real(dp)::den
    logical::norm
    n=min(size(x),size(y));norm=.false.;if(present(normalize))norm=normalize
    if(n==0)then;v=0.0_dp;return;end if
    v=sum((x(:n)-mean_value(x(:n)))**p1*(y(:n)-mean_value(y(:n)))**p2)/real(n,dp)
    if(norm)then
      den=sum((y(:n)-mean_value(y(:n)))**(p1+p2))/real(n,dp)
      if(abs(den)>tiny(1.0_dp))v=v/den
    end if
  end function coskewness_value

  pure real(dp) function cokurtosis_value(x,y,p1,p2,normalize) result(v)
    real(dp),intent(in)::x(:),y(:)
    integer,intent(in)::p1,p2
    logical,intent(in),optional::normalize
    if(present(normalize))then;v=coskewness_value(x,y,p1,p2,normalize);else;v=coskewness_value(x,y,p1,p2);end if
  end function cokurtosis_value

  subroutine coskewness_matrix(r,m3,unbiased)
    real(dp),intent(in)::r(:,:)
    real(dp),intent(out)::m3(:,:)
    logical,intent(in),optional::unbiased
    real(dp),allocatable::xc(:,:)
    real(dp)::c
    integer::n,p,i,j,k,t,col
    logical::u
    n=size(r,1);p=size(r,2);u=.false.;if(present(unbiased))u=unbiased
    allocate(xc(n,p));do j=1,p;xc(:,j)=r(:,j)-mean_value(r(:,j));end do
    if(u .and. n>2)then;c=real(n,dp)/real((n-1)*(n-2),dp);else;c=1.0_dp/real(max(n,1),dp);end if
    m3=0.0_dp
    do i=1,p;do j=1,p;do k=1,p
      col=(j-1)*p+k
      do t=1,n;m3(i,col)=m3(i,col)+xc(t,i)*xc(t,j)*xc(t,k);end do
      m3(i,col)=c*m3(i,col)
    end do;end do;end do
  end subroutine coskewness_matrix

  subroutine cokurtosis_matrix(r,m4)
    real(dp),intent(in)::r(:,:)
    real(dp),intent(out)::m4(:,:)
    real(dp),allocatable::xc(:,:)
    integer::n,p,i,j,k,l,t,col
    n=size(r,1);p=size(r,2);allocate(xc(n,p));do j=1,p;xc(:,j)=r(:,j)-mean_value(r(:,j));end do
    m4=0.0_dp
    do i=1,p;do j=1,p;do k=1,p;do l=1,p
      col=((j-1)*p+k-1)*p+l
      do t=1,n;m4(i,col)=m4(i,col)+xc(t,i)*xc(t,j)*xc(t,k)*xc(t,l);end do
      m4(i,col)=m4(i,col)/real(max(n,1),dp)
    end do;end do;end do;end do
  end subroutine cokurtosis_matrix

  subroutine coskewness_unique(r,v,unbiased)
    real(dp),intent(in)::r(:,:)
    real(dp),intent(out)::v(:)
    logical,intent(in),optional::unbiased
    real(dp),allocatable::xc(:,:)
    real(dp)::c,s
    integer::n,p,i,j,k,t,q
    logical::u
    n=size(r,1);p=size(r,2);u=.false.;if(present(unbiased))u=unbiased
    allocate(xc(n,p));do j=1,p;xc(:,j)=r(:,j)-mean_value(r(:,j));end do
    if(u .and. n>2)then;c=real(n,dp)/real((n-1)*(n-2),dp);else;c=1.0_dp/real(max(n,1),dp);end if
    q=0
    do i=1,p;do j=i,p;do k=j,p
      q=q+1;s=0.0_dp;do t=1,n;s=s+xc(t,i)*xc(t,j)*xc(t,k);end do
      if(q<=size(v))v(q)=c*s
    end do;end do;end do
  end subroutine coskewness_unique

  subroutine cokurtosis_unique(r,v)
    real(dp),intent(in)::r(:,:)
    real(dp),intent(out)::v(:)
    real(dp),allocatable::xc(:,:)
    real(dp)::s
    integer::n,p,i,j,k,l,t,q
    n=size(r,1);p=size(r,2);allocate(xc(n,p));do j=1,p;xc(:,j)=r(:,j)-mean_value(r(:,j));end do
    q=0
    do i=1,p;do j=i,p;do k=j,p;do l=k,p
      q=q+1;s=0.0_dp;do t=1,n;s=s+xc(t,i)*xc(t,j)*xc(t,k)*xc(t,l);end do
      if(q<=size(v))v(q)=s/real(max(n,1),dp)
    end do;end do;end do;end do
  end subroutine cokurtosis_unique

  subroutine m3_vec_to_mat(v,p,m3)
    real(dp),intent(in)::v(:)
    integer,intent(in)::p
    real(dp),intent(out)::m3(:,:)
    integer::i,j,k,q
    q=0;m3=0.0_dp
    do i=1,p;do j=i,p;do k=j,p
      q=q+1
      call set_m3_permutations(m3,p,i,j,k,v(q))
    end do;end do;end do
  end subroutine m3_vec_to_mat

  subroutine set_m3_permutations(m3,p,i,j,k,value)
    real(dp),intent(inout)::m3(:,:)
    integer,intent(in)::p,i,j,k
    real(dp),intent(in)::value
    integer::a(6),b(6),c(6),q,col
    a=[i,i,j,j,k,k];b=[j,k,i,k,i,j];c=[k,j,k,i,j,i]
    do q=1,6;col=(b(q)-1)*p+c(q);m3(a(q),col)=value;end do
  end subroutine set_m3_permutations

  subroutine m4_vec_to_mat(v,p,m4)
    real(dp),intent(in)::v(:)
    integer,intent(in)::p
    real(dp),intent(out)::m4(:,:)
    integer::i,j,k,l,q,a,b,c,d,col
    q=0;m4=0.0_dp
    do i=1,p;do j=i,p;do k=j,p;do l=k,p
      q=q+1
      do a=1,p;do b=1,p;do c=1,p;do d=1,p
        if(same_multiset4([a,b,c,d],[i,j,k,l]))then
          col=((b-1)*p+c-1)*p+d;m4(a,col)=v(q)
        end if
      end do;end do;end do;end do
    end do;end do;end do;end do
  end subroutine m4_vec_to_mat

  pure logical function same_multiset4(x,y) result(ok)
    integer,intent(in)::x(4),y(4)
    integer::a(4),b(4),i,j,tmp
    a=x;b=y
    do i=2,4
      tmp=a(i);j=i-1
      do while(j>=1)
        if(a(j)<=tmp) exit
        a(j+1)=a(j);j=j-1
      end do
      a(j+1)=tmp
    end do
    do i=2,4
      tmp=b(i);j=i-1
      do while(j>=1)
        if(b(j)<=tmp) exit
        b(j+1)=b(j);j=j-1
      end do
      b(j+1)=tmp
    end do
    ok=all(a==b)
  end function same_multiset4

  subroutine m3_mat_to_vec(m3,p,v)
    real(dp),intent(in)::m3(:,:)
    integer,intent(in)::p
    real(dp),intent(out)::v(:)
    integer::i,j,k,q
    q=0
    do i=1,p;do j=i,p;do k=j,p;q=q+1;v(q)=m3(i,(j-1)*p+k);end do;end do;end do
  end subroutine m3_mat_to_vec

  subroutine m4_mat_to_vec(m4,p,v)
    real(dp),intent(in)::m4(:,:)
    integer,intent(in)::p
    real(dp),intent(out)::v(:)
    integer::i,j,k,l,q
    q=0
    do i=1,p;do j=i,p;do k=j,p;do l=k,p;q=q+1;v(q)=m4(i,((j-1)*p+k-1)*p+l);end do;end do;end do;end do
  end subroutine m4_mat_to_vec

  pure real(dp) function m3_inner_product(v1,v2,p) result(s)
    real(dp),intent(in)::v1(:),v2(:)
    integer,intent(in)::p
    integer::i,j,k,q,mult
    s=0.0_dp;q=0
    do i=1,p;do j=i,p;do k=j,p
      q=q+1
      if(i==k)then;mult=1;else if(i==j .or. j==k)then;mult=3;else;mult=6;end if
      s=s+real(mult,dp)*v1(q)*v2(q)
    end do;end do;end do
  end function m3_inner_product

  pure real(dp) function m4_inner_product(v1,v2,p) result(s)
    real(dp),intent(in)::v1(:),v2(:)
    integer,intent(in)::p
    integer::i,j,k,l,q,mult
    s=0.0_dp;q=0
    do i=1,p;do j=i,p;do k=j,p;do l=k,p
      q=q+1
      if(i==l)then
        mult=1
      else if((i==k .and. k<l) .or. (i<j .and. j==l))then
        mult=4
      else if(i==j .and. k==l .and. j<k)then
        mult=6
      else if(i==j .or. j==k .or. k==l)then
        mult=12
      else
        mult=24
      end if
      s=s+real(mult,dp)*v1(q)*v2(q)
    end do;end do;end do;end do
  end function m4_inner_product

  pure real(dp) function portfolio_mean(w,mu) result(v)
    real(dp),intent(in)::w(:),mu(:)
    v=dot_product(w,mu)
  end function portfolio_mean

  pure real(dp) function portfolio_variance(w,cov) result(v)
    real(dp),intent(in)::w(:),cov(:,:)
    v=dot_product(w,matmul(cov,w))
  end function portfolio_variance

  pure real(dp) function portfolio_third_moment(w,m3) result(v)
    real(dp),intent(in)::w(:),m3(:,:)
    integer::p,i,j,k
    p=size(w);v=0.0_dp
    do i=1,p;do j=1,p;do k=1,p;v=v+w(i)*w(j)*w(k)*m3(i,(j-1)*p+k);end do;end do;end do
  end function portfolio_third_moment

  pure real(dp) function portfolio_fourth_moment(w,m4) result(v)
    real(dp),intent(in)::w(:),m4(:,:)
    integer::p,i,j,k,l
    p=size(w);v=0.0_dp
    do i=1,p;do j=1,p;do k=1,p;do l=1,p
      v=v+w(i)*w(j)*w(k)*w(l)*m4(i,((j-1)*p+k-1)*p+l)
    end do;end do;end do;end do
  end function portfolio_fourth_moment

  pure real(dp) function portfolio_skewness(w,cov,m3) result(v)
    real(dp),intent(in)::w(:),cov(:,:),m3(:,:)
    real(dp)::s2
    s2=portfolio_variance(w,cov)
    if(s2<=tiny(1.0_dp))then;v=0.0_dp;else;v=portfolio_third_moment(w,m3)/s2**1.5_dp;end if
  end function portfolio_skewness

  pure real(dp) function portfolio_kurtosis(w,cov,m4) result(v)
    real(dp),intent(in)::w(:),cov(:,:),m4(:,:)
    real(dp)::s2
    s2=portfolio_variance(w,cov)
    if(s2<=tiny(1.0_dp))then;v=0.0_dp;else;v=portfolio_fourth_moment(w,m4)/(s2*s2);end if
  end function portfolio_kurtosis

  subroutine ewma_covariance(r,lambda,cov)
    real(dp),intent(in)::r(:,:),lambda
    real(dp),intent(out)::cov(:,:)
    real(dp),allocatable::mu(:),d(:)
    real(dp)::weight,total
    integer::n,p,t
    n=size(r,1);p=size(r,2);allocate(mu(p),d(p));mu=0.0_dp;cov=0.0_dp;total=0.0_dp
    do t=1,n
      weight=(1.0_dp-lambda)*lambda**real(n-t,dp);mu=mu+weight*r(t,:);total=total+weight
    end do
    if(total>0.0_dp)mu=mu/total
    do t=1,n
      weight=(1.0_dp-lambda)*lambda**real(n-t,dp);d=r(t,:)-mu
      cov=cov+weight*spread(d,2,p)*spread(d,1,p)
    end do
    if(total>0.0_dp)cov=cov/total
  end subroutine ewma_covariance

  pure real(dp) function beta_coskewness(asset,market) result(v)
    real(dp),intent(in)::asset(:),market(:)
    real(dp)::den
    den=sum((market-mean_value(market))**3)/real(max(size(market),1),dp)
    if(abs(den)<=tiny(1.0_dp))then;v=0.0_dp;else;v=coskewness_value(asset,market,1,2)/den;end if
  end function beta_coskewness

  pure real(dp) function beta_cokurtosis(asset,market) result(v)
    real(dp),intent(in)::asset(:),market(:)
    real(dp)::den
    den=sum((market-mean_value(market))**4)/real(max(size(market),1),dp)
    if(abs(den)<=tiny(1.0_dp))then;v=0.0_dp;else;v=cokurtosis_value(asset,market,1,3)/den;end if
  end function beta_cokurtosis
end module comoments_mod
