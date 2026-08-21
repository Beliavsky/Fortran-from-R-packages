! SPDX-License-Identifier: GPL-2.0-or-later
module mnb_math
  use mnb_kinds, only : dp
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  implicit none
  private
  public :: normal_cdf, normal_quantile, gamma_rng, poisson_rng, set_mnb_seed
  public :: invert_matrix, symmetric_eigen_jacobi, covariance_from_hessian
  public :: numerical_hessian, numerical_gradient, sample_sd_unique

  abstract interface
    function scalar_function(x) result(f)
      import dp
      real(dp), intent(in) :: x(:)
      real(dp) :: f
    end function scalar_function
  end interface

contains

  pure real(dp) function normal_cdf(x) result(p)
    real(dp), intent(in) :: x
    p = 0.5_dp * erfc(-x/sqrt(2.0_dp))
  end function normal_cdf

  pure real(dp) function normal_quantile(p) result(x)
    real(dp), intent(in) :: p
    real(dp), parameter :: a1=-3.969683028665376d1, a2=2.209460984245205d2
    real(dp), parameter :: a3=-2.759285104469687d2, a4=1.383577518672690d2
    real(dp), parameter :: a5=-3.066479806614716d1, a6=2.506628277459239d0
    real(dp), parameter :: b1=-5.447609879822406d1, b2=1.615858368580409d2
    real(dp), parameter :: b3=-1.556989798598866d2, b4=6.680131188771972d1
    real(dp), parameter :: b5=-1.328068155288572d1
    real(dp), parameter :: c1=-7.784894002430293d-3, c2=-3.223964580411365d-1
    real(dp), parameter :: c3=-2.400758277161838d0, c4=-2.549732539343734d0
    real(dp), parameter :: c5=4.374664141464968d0, c6=2.938163982698783d0
    real(dp), parameter :: d1=7.784695709041462d-3, d2=3.224671290700398d-1
    real(dp), parameter :: d3=2.445134137142996d0, d4=3.754408661907416d0
    real(dp), parameter :: plow=0.02425_dp, phigh=1.0_dp-plow
    real(dp) :: q, r
    if (p <= 0.0_dp) then
      x = -huge(1.0_dp)
    else if (p >= 1.0_dp) then
      x = huge(1.0_dp)
    else if (p < plow) then
      q = sqrt(-2.0_dp*log(p))
      x = (((((c1*q+c2)*q+c3)*q+c4)*q+c5)*q+c6) / ((((d1*q+d2)*q+d3)*q+d4)*q+1.0_dp)
    else if (p <= phigh) then
      q = p-0.5_dp; r=q*q
      x = (((((a1*r+a2)*r+a3)*r+a4)*r+a5)*r+a6)*q / (((((b1*r+b2)*r+b3)*r+b4)*r+b5)*r+1.0_dp)
    else
      q = sqrt(-2.0_dp*log(1.0_dp-p))
      x = -(((((c1*q+c2)*q+c3)*q+c4)*q+c5)*q+c6) / ((((d1*q+d2)*q+d3)*q+d4)*q+1.0_dp)
    end if
  end function normal_quantile

  subroutine set_mnb_seed(seed)
    integer, intent(in) :: seed
    integer :: n, i
    integer, allocatable :: s(:)
    call random_seed(size=n); allocate(s(n))
    do i=1,n
      s(i) = modulo(seed + 104729*i, huge(1)-1)
      if (s(i) == 0) s(i)=i
    end do
    call random_seed(put=s)
  end subroutine set_mnb_seed

  real(dp) function randn() result(z)
    real(dp) :: u1,u2
    call random_number(u1); call random_number(u2)
    u1=max(u1,tiny(1.0_dp)); z=sqrt(-2.0_dp*log(u1))*cos(2.0_dp*acos(-1.0_dp)*u2)
  end function randn

  recursive real(dp) function gamma_rng(shape) result(x)
    real(dp), intent(in) :: shape
    real(dp) :: d,c,z,u,v
    if (shape <= 0.0_dp) then
      x=0.0_dp; return
    end if
    if (shape < 1.0_dp) then
      call random_number(u)
      x=gamma_rng(shape+1.0_dp)*u**(1.0_dp/shape)
      return
    end if
    d=shape-1.0_dp/3.0_dp; c=1.0_dp/sqrt(9.0_dp*d)
    do
      z=randn(); v=1.0_dp+c*z
      if (v <= 0.0_dp) cycle
      v=v*v*v; call random_number(u)
      if (u < 1.0_dp-0.0331_dp*z**4) exit
      if (log(max(u,tiny(1.0_dp))) < 0.5_dp*z*z+d*(1.0_dp-v+log(v))) exit
    end do
    x=d*v
  end function gamma_rng

  integer function poisson_rng(lambda) result(k)
    real(dp), intent(in) :: lambda
    real(dp) :: l,p,u,z
    integer :: j
    if (lambda <= 0.0_dp) then; k=0; return; end if
    if (lambda < 30.0_dp) then
      l=exp(-lambda); p=1.0_dp; j=0
      do
        j=j+1; call random_number(u); p=p*u
        if (p <= l) exit
      end do
      k=j-1
    else
      do
        z=randn(); k=nint(lambda+sqrt(lambda)*z)
        if (k >= 0) exit
      end do
    end if
  end function poisson_rng

  subroutine numerical_gradient(fn,x,g)
    procedure(scalar_function) :: fn
    real(dp), intent(in) :: x(:)
    real(dp), intent(out) :: g(:)
    real(dp) :: h,fp,fm
    real(dp), allocatable :: xp(:),xm(:)
    integer :: i
    allocate(xp(size(x)),xm(size(x)))
    do i=1,size(x)
      h=1.0e-5_dp*max(1.0_dp,abs(x(i)))
      xp=x; xm=x; xp(i)=x(i)+h; xm(i)=x(i)-h
      fp=fn(xp); fm=fn(xm); g(i)=(fp-fm)/(2.0_dp*h)
    end do
  end subroutine numerical_gradient

  subroutine numerical_hessian(fn,x,hess)
    procedure(scalar_function) :: fn
    real(dp), intent(in) :: x(:)
    real(dp), intent(out) :: hess(:,:)
    real(dp), allocatable :: xpp(:),xpm(:),xmp(:),xmm(:),xp(:),xm(:)
    real(dp) :: hi,hj,f0,fp,fm
    integer :: i,j,n
    n=size(x); allocate(xpp(n),xpm(n),xmp(n),xmm(n),xp(n),xm(n)); f0=fn(x)
    do i=1,n
      hi=1.0e-4_dp*max(1.0_dp,abs(x(i))); xp=x; xm=x
      xp(i)=x(i)+hi; xm(i)=x(i)-hi; fp=fn(xp); fm=fn(xm)
      hess(i,i)=(fp-2.0_dp*f0+fm)/(hi*hi)
      do j=1,i-1
        hj=1.0e-4_dp*max(1.0_dp,abs(x(j)))
        xpp=x; xpm=x; xmp=x; xmm=x
        xpp(i)=x(i)+hi; xpp(j)=x(j)+hj
        xpm(i)=x(i)+hi; xpm(j)=x(j)-hj
        xmp(i)=x(i)-hi; xmp(j)=x(j)+hj
        xmm(i)=x(i)-hi; xmm(j)=x(j)-hj
        hess(i,j)=(fn(xpp)-fn(xpm)-fn(xmp)+fn(xmm))/(4.0_dp*hi*hj)
        hess(j,i)=hess(i,j)
      end do
    end do
  end subroutine numerical_hessian

  subroutine invert_matrix(a,ainv,ok)
    real(dp), intent(in) :: a(:,:)
    real(dp), intent(out) :: ainv(:,:)
    logical, intent(out) :: ok
    real(dp), allocatable :: aug(:,:),tmp(:)
    real(dp) :: piv,fac
    integer :: n,i,k,ip
    n=size(a,1); allocate(aug(n,2*n),tmp(2*n)); aug=0.0_dp
    aug(:,1:n)=a
    do i=1,n; aug(i,n+i)=1.0_dp; end do
    ok=.true.
    do k=1,n
      ip=k
      do i=k+1,n
        if(abs(aug(i,k))>abs(aug(ip,k))) ip=i
      end do
      if(abs(aug(ip,k)) <= 100.0_dp*epsilon(1.0_dp)) then; ok=.false.; ainv=0.0_dp; return; end if
      if(ip/=k) then; tmp=aug(k,:); aug(k,:)=aug(ip,:); aug(ip,:)=tmp; end if
      piv=aug(k,k); aug(k,:)=aug(k,:)/piv
      do i=1,n
        if(i==k) cycle
        fac=aug(i,k); aug(i,:)=aug(i,:)-fac*aug(k,:)
      end do
    end do
    ainv=aug(:,n+1:2*n)
  end subroutine invert_matrix

  subroutine covariance_from_hessian(hess,cov,ok)
    real(dp), intent(in) :: hess(:,:)
    real(dp), intent(out) :: cov(:,:)
    logical, intent(out) :: ok
    call invert_matrix(-hess,cov,ok)
  end subroutine covariance_from_hessian

  subroutine symmetric_eigen_jacobi(a,values,vectors)
    real(dp), intent(in) :: a(:,:)
    real(dp), intent(out) :: values(:),vectors(:,:)
    real(dp), allocatable :: b(:,:)
    real(dp) :: app,aqq,apq,tau,t,c,s,bip,biq,vip,viq
    integer :: n,p,q,i,iter,maxiter,k
    n=size(a,1); allocate(b(n,n)); b=a; vectors=0.0_dp
    do i=1,n; vectors(i,i)=1.0_dp; end do
    maxiter=max(100,50*n*n)
    do iter=1,maxiter
      p=1;q=min(2,n)
      if(n==1) exit
      do i=1,n
        do k=i+1,n
          if(abs(b(i,k))>abs(b(p,q))) then;p=i;q=k;end if
        end do
      end do
      if(abs(b(p,q)) < 1.0e-12_dp) exit
      app=b(p,p); aqq=b(q,q); apq=b(p,q); tau=(aqq-app)/(2.0_dp*apq)
      t=sign(1.0_dp,tau)/(abs(tau)+sqrt(1.0_dp+tau*tau)); c=1.0_dp/sqrt(1.0_dp+t*t); s=t*c
      do i=1,n
        if(i/=p .and. i/=q) then
          bip=b(i,p); biq=b(i,q); b(i,p)=c*bip-s*biq; b(p,i)=b(i,p)
          b(i,q)=s*bip+c*biq; b(q,i)=b(i,q)
        end if
        vip=vectors(i,p); viq=vectors(i,q); vectors(i,p)=c*vip-s*viq; vectors(i,q)=s*vip+c*viq
      end do
      b(p,p)=c*c*app-2.0_dp*s*c*apq+s*s*aqq
      b(q,q)=s*s*app+2.0_dp*s*c*apq+c*c*aqq; b(p,q)=0.0_dp;b(q,p)=0.0_dp
    end do
    do i=1,n; values(i)=b(i,i); end do
    call sort_eigen_desc(values,vectors)
  contains
    subroutine sort_eigen_desc(v,z)
      real(dp),intent(inout)::v(:),z(:,:)
      real(dp)::tv;real(dp),allocatable::col(:);integer::ii,jj,imax
      allocate(col(size(z,1)))
      do ii=1,size(v)-1
        imax=ii
        do jj=ii+1,size(v);if(v(jj)>v(imax))imax=jj;end do
        if(imax/=ii)then;tv=v(ii);v(ii)=v(imax);v(imax)=tv;col=z(:,ii);z(:,ii)=z(:,imax);z(:,imax)=col;end if
      end do
    end subroutine
  end subroutine symmetric_eigen_jacobi

  real(dp) function sample_sd_unique(x) result(s)
    real(dp), intent(in) :: x(:)
    real(dp), allocatable :: u(:)
    integer :: i,j,n
    logical :: seen
    allocate(u(size(x))); n=0
    do i=1,size(x)
      seen=.false.
      do j=1,n
        if(abs(x(i)-u(j)) <= 0.0_dp) then; seen=.true.; exit; end if
      end do
      if(.not.seen) then;n=n+1;u(n)=x(i);end if
    end do
    if(n<=1)then;s=0.0_dp;else;s=sqrt(sum((u(1:n)-sum(u(1:n))/real(n,dp))**2)/real(n-1,dp));end if
  end function sample_sd_unique
end module mnb_math
