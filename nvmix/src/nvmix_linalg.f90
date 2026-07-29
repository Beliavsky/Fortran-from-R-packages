! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2018 Marius Hofert, Erik Hintz and Christiane Lemieux
module nvmix_linalg
  use nvmix_kinds, only : dp
  implicit none
  private
  public :: cholesky_lower, forward_solve, backward_solve_transpose, solve_spd
  public :: quadratic_form_spd, sample_mean_covariance, covariance_to_correlation
  public :: inverse_matrix, symmetric_eigen_jacobi, nearest_psd
contains
  subroutine cholesky_lower(a,l,ok)
    real(dp), intent(in) :: a(:,:)
    real(dp), allocatable, intent(out) :: l(:,:)
    logical, intent(out) :: ok
    integer :: n,i,j,k
    real(dp) :: s
    n=size(a,1); allocate(l(n,n)); l=0.0_dp; ok=.false.
    if(size(a,2)/=n)return
    do i=1,n
      do j=1,i
        s=a(i,j)
        do k=1,j-1; s=s-l(i,k)*l(j,k); end do
        if(i==j) then
          if(s<=100.0_dp*epsilon(1.0_dp)*max(1.0_dp,maxval(abs(a))))return
          l(i,j)=sqrt(s)
        else
          l(i,j)=s/l(j,j)
        end if
      end do
    end do
    ok=.true.
  end subroutine
  subroutine forward_solve(l,b,x,ok)
    real(dp), intent(in) :: l(:,:),b(:)
    real(dp), intent(out) :: x(:)
    logical, intent(out) :: ok
    integer :: i,j,n
    n=size(b); ok=.false.; if(size(l,1)/=n .or. size(l,2)/=n .or. size(x)/=n)return
    do i=1,n
      if(abs(l(i,i))<=tiny(1.0_dp))return
      x(i)=b(i)
      do j=1,i-1; x(i)=x(i)-l(i,j)*x(j); end do
      x(i)=x(i)/l(i,i)
    end do
    ok=.true.
  end subroutine
  subroutine backward_solve_transpose(l,b,x,ok)
    real(dp), intent(in) :: l(:,:),b(:)
    real(dp), intent(out) :: x(:)
    logical, intent(out) :: ok
    integer :: i,j,n
    n=size(b); ok=.false.; if(size(l,1)/=n .or. size(l,2)/=n .or. size(x)/=n)return
    do i=n,1,-1
      if(abs(l(i,i))<=tiny(1.0_dp))return
      x(i)=b(i)
      do j=i+1,n; x(i)=x(i)-l(j,i)*x(j); end do
      x(i)=x(i)/l(i,i)
    end do
    ok=.true.
  end subroutine
  subroutine solve_spd(a,b,x,ok)
    real(dp), intent(in) :: a(:,:),b(:)
    real(dp), intent(out) :: x(:)
    logical, intent(out) :: ok
    real(dp), allocatable :: l(:,:),y(:)
    logical :: ok1
    allocate(y(size(b))); call cholesky_lower(a,l,ok); if(.not.ok)return
    call forward_solve(l,b,y,ok); if(.not.ok)return
    call backward_solve_transpose(l,y,x,ok1); ok=ok1
  end subroutine
  real(dp) function quadratic_form_spd(a,x,logdet,ok) result(q)
    real(dp), intent(in) :: a(:,:),x(:)
    real(dp), intent(out), optional :: logdet
    logical, intent(out), optional :: ok
    real(dp), allocatable :: l(:,:),y(:)
    logical :: good
    integer :: i
    call cholesky_lower(a,l,good)
    if(.not.good) then; q=huge(1.0_dp); if(present(ok))ok=.false.; if(present(logdet))logdet=huge(1.0_dp); return; end if
    allocate(y(size(x))); call forward_solve(l,x,y,good); q=dot_product(y,y)
    if(present(logdet))then; logdet=0.0_dp; do i=1,size(x); logdet=logdet+2.0_dp*log(l(i,i)); end do; end if
    if(present(ok))ok=good
  end function
  subroutine sample_mean_covariance(x,mean,covariance,ok)
    real(dp), intent(in) :: x(:,:)
    real(dp), allocatable, intent(out) :: mean(:),covariance(:,:)
    logical, intent(out) :: ok
    integer :: n,d,i,j,k
    real(dp), allocatable :: delta(:)
    n=size(x,1); d=size(x,2); allocate(mean(d),covariance(d,d),delta(d)); mean=0.0_dp; covariance=0.0_dp
    ok=n>=2 .and. d>=1; if(.not.ok)return
    do i=1,n; mean=mean+x(i,:); end do; mean=mean/real(n,dp)
    do i=1,n; delta=x(i,:)-mean; do j=1,d; do k=1,j; covariance(j,k)=covariance(j,k)+delta(j)*delta(k); end do; end do; end do
    covariance=covariance/real(n-1,dp)
    do j=1,d; do k=1,j-1; covariance(k,j)=covariance(j,k); end do; end do
  end subroutine
  subroutine covariance_to_correlation(covariance,correlation,sd,ok)
    real(dp), intent(in) :: covariance(:,:)
    real(dp), allocatable, intent(out) :: correlation(:,:),sd(:)
    logical, intent(out) :: ok
    integer :: d,i,j
    d=size(covariance,1); allocate(correlation(d,d),sd(d)); ok=size(covariance,2)==d
    if(.not.ok)return
    do i=1,d
      if(covariance(i,i)<=0.0_dp)then; ok=.false.; return; end if
      sd(i)=sqrt(covariance(i,i))
    end do
    do i=1,d; do j=1,d; correlation(i,j)=covariance(i,j)/(sd(i)*sd(j)); end do; end do
    do i=1,d; correlation(i,i)=1.0_dp; end do
  end subroutine
  subroutine inverse_matrix(a,ainv,ok)
    real(dp), intent(in) :: a(:,:)
    real(dp), allocatable, intent(out) :: ainv(:,:)
    logical, intent(out) :: ok
    real(dp), allocatable :: aug(:,:),row(:)
    real(dp) :: pivot,factor
    integer :: n,i,j,k,p
    n=size(a,1); allocate(ainv(n,n),aug(n,2*n),row(2*n)); ok=.false.; if(size(a,2)/=n)return
    aug(:,1:n)=a; aug(:,n+1:2*n)=0.0_dp; do i=1,n; aug(i,n+i)=1.0_dp; end do
    do i=1,n
      p=i; do k=i+1,n; if(abs(aug(k,i))>abs(aug(p,i)))p=k; end do
      if(abs(aug(p,i))<=100.0_dp*epsilon(1.0_dp))return
      if(p/=i)then; row=aug(i,:); aug(i,:)=aug(p,:); aug(p,:)=row; end if
      pivot=aug(i,i); aug(i,:)=aug(i,:)/pivot
      do j=1,n
        if(j==i)cycle; factor=aug(j,i); aug(j,:)=aug(j,:)-factor*aug(i,:)
      end do
    end do
    ainv=aug(:,n+1:2*n); ok=.true.
  end subroutine
  subroutine symmetric_eigen_jacobi(a,values,vectors,ok)
    real(dp), intent(in) :: a(:,:)
    real(dp), allocatable, intent(out) :: values(:),vectors(:,:)
    logical, intent(out) :: ok
    real(dp), allocatable :: b(:,:)
    real(dp) :: app,aqq,apq,tau,t,c,s,bkp,bkq,vkp,vkq
    integer :: n,p,q,k,iter
    n=size(a,1); allocate(b(n,n),values(n),vectors(n,n)); b=0.5_dp*(a+transpose(a)); vectors=0.0_dp
    do k=1,n; vectors(k,k)=1.0_dp; end do; ok=.false.
    do iter=1,100*n*n
      p=1; q=min(2,n); apq=0.0_dp
      do k=1,n; if(k<n)then
        do q=k+1,n; if(abs(b(k,q))>abs(apq))then; apq=b(k,q); p=k; end if; end do
      end if; end do
      q=p+1; apq=0.0_dp
      do k=p+1,n; if(abs(b(p,k))>abs(apq))then; apq=b(p,k); q=k; end if; end do
      if(abs(apq)<1.0e-12_dp*max(1.0_dp,maxval(abs(b))))then; ok=.true.; exit; end if
      app=b(p,p); aqq=b(q,q); tau=(aqq-app)/(2.0_dp*apq)
      if(tau>=0.0_dp)then; t=1.0_dp/(tau+sqrt(1.0_dp+tau*tau)); else; t=-1.0_dp/(-tau+sqrt(1.0_dp+tau*tau)); end if
      c=1.0_dp/sqrt(1.0_dp+t*t); s=t*c
      do k=1,n
        if(k/=p .and. k/=q)then
          bkp=b(k,p); bkq=b(k,q); b(k,p)=c*bkp-s*bkq; b(p,k)=b(k,p); b(k,q)=s*bkp+c*bkq; b(q,k)=b(k,q)
        end if
        vkp=vectors(k,p); vkq=vectors(k,q); vectors(k,p)=c*vkp-s*vkq; vectors(k,q)=s*vkp+c*vkq
      end do
      b(p,p)=app-t*apq; b(q,q)=aqq+t*apq; b(p,q)=0.0_dp; b(q,p)=0.0_dp
    end do
    do k=1,n; values(k)=b(k,k); end do
  end subroutine
  subroutine nearest_psd(a,out,floor_value,ok)
    real(dp), intent(in) :: a(:,:)
    real(dp), allocatable, intent(out) :: out(:,:)
    real(dp), intent(in), optional :: floor_value
    logical, intent(out) :: ok
    real(dp), allocatable :: vals(:),vecs(:,:),diagv(:,:)
    real(dp) :: fl
    integer :: n,i
    fl=1.0e-10_dp; if(present(floor_value))fl=floor_value
    call symmetric_eigen_jacobi(a,vals,vecs,ok); if(.not.ok)return
    n=size(vals); allocate(diagv(n,n),out(n,n)); diagv=0.0_dp
    do i=1,n; diagv(i,i)=max(vals(i),fl); end do
    out=matmul(vecs,matmul(diagv,transpose(vecs))); out=0.5_dp*(out+transpose(out))
  end subroutine
end module nvmix_linalg
