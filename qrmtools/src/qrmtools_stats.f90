! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2016 Marius Hofert, Kurt Hornik and Alexander J. McNeil
module qrmtools_stats
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_value, ieee_quiet_nan
  use qrmtools_kinds, only : dp, pi
  implicit none
  private
  public :: normal_pdf, normal_cdf, normal_quantile
  public :: student_t_pdf, student_t_cdf, student_t_quantile
  public :: gamma_cdf, beta_cdf, chi_square_cdf
  public :: seed_random, random_normal, random_uniform
  public :: mean_value, variance_value, covariance_matrix
  public :: quantile_type1, sort_increasing, stable_order
  public :: invert_matrix, solve_linear, cholesky_factor
  public :: regularized_beta, numerical_gradient, numerical_hessian
  public :: trapz_integral, bisection_root

  abstract interface
    pure function scalar_function(x) result(value)
      import dp
      real(dp), intent(in) :: x
      real(dp) :: value
    end function scalar_function
    function vector_objective(x) result(value)
      import dp
      real(dp), intent(in) :: x(:)
      real(dp) :: value
    end function vector_objective
  end interface

contains
  pure real(dp) function normal_pdf(x) result(value)
    real(dp), intent(in) :: x
    value = exp(-0.5_dp*x*x)/sqrt(2.0_dp*pi)
  end function normal_pdf

  pure real(dp) function normal_cdf(x) result(value)
    real(dp), intent(in) :: x
    value = 0.5_dp*erfc(-x/sqrt(2.0_dp))
  end function normal_cdf

  pure real(dp) function normal_quantile(p) result(x)
    real(dp), intent(in) :: p
    real(dp), parameter :: a(6)=[-3.969683028665376e1_dp,2.209460984245205e2_dp,&
      -2.759285104469687e2_dp,1.383577518672690e2_dp,-3.066479806614716e1_dp,2.506628277459239_dp]
    real(dp), parameter :: b(5)=[-5.447609879822406e1_dp,1.615858368580409e2_dp,&
      -1.556989798598866e2_dp,6.680131188771972e1_dp,-1.328068155288572e1_dp]
    real(dp), parameter :: c(6)=[-7.784894002430293e-3_dp,-3.223964580411365e-1_dp,&
      -2.400758277161838_dp,-2.549732539343734_dp,4.374664141464968_dp,2.938163982698783_dp]
    real(dp), parameter :: d(4)=[7.784695709041462e-3_dp,3.224671290700398e-1_dp,&
      2.445134137142996_dp,3.754408661907416_dp]
    real(dp), parameter :: plow=0.02425_dp, phigh=1.0_dp-plow
    real(dp) :: q,r
    if (p <= 0.0_dp) then
      x=-huge(1.0_dp)
    else if (p >= 1.0_dp) then
      x=huge(1.0_dp)
    else if (p < plow) then
      q=sqrt(-2.0_dp*log(p))
      x=(((((c(1)*q+c(2))*q+c(3))*q+c(4))*q+c(5))*q+c(6))/&
        ((((d(1)*q+d(2))*q+d(3))*q+d(4))*q+1.0_dp)
    else if (p <= phigh) then
      q=p-0.5_dp; r=q*q
      x=(((((a(1)*r+a(2))*r+a(3))*r+a(4))*r+a(5))*r+a(6))*q/&
        (((((b(1)*r+b(2))*r+b(3))*r+b(4))*r+b(5))*r+1.0_dp)
    else
      q=sqrt(-2.0_dp*log(1.0_dp-p))
      x=-(((((c(1)*q+c(2))*q+c(3))*q+c(4))*q+c(5))*q+c(6))/&
        ((((d(1)*q+d(2))*q+d(3))*q+d(4))*q+1.0_dp)
    end if
  end function normal_quantile

  pure real(dp) function beta_continued_fraction(a,b,x) result(cf)
    real(dp), intent(in) :: a,b,x
    integer, parameter :: max_iter=300
    real(dp), parameter :: eps=3.0e-14_dp, fpmin=1.0e-300_dp
    integer :: m,m2
    real(dp) :: aa,c,d,h,del,qab,qam,qap
    qab=a+b; qap=a+1.0_dp; qam=a-1.0_dp
    c=1.0_dp; d=1.0_dp-qab*x/qap
    if (abs(d)<fpmin) d=fpmin
    d=1.0_dp/d; h=d
    do m=1,max_iter
      m2=2*m
      aa=real(m,dp)*(b-real(m,dp))*x/((qam+real(m2,dp))*(a+real(m2,dp)))
      d=1.0_dp+aa*d; if(abs(d)<fpmin)d=fpmin
      c=1.0_dp+aa/c; if(abs(c)<fpmin)c=fpmin
      d=1.0_dp/d; h=h*d*c
      aa=-(a+real(m,dp))*(qab+real(m,dp))*x/((a+real(m2,dp))*(qap+real(m2,dp)))
      d=1.0_dp+aa*d; if(abs(d)<fpmin)d=fpmin
      c=1.0_dp+aa/c; if(abs(c)<fpmin)c=fpmin
      d=1.0_dp/d; del=d*c; h=h*del
      if(abs(del-1.0_dp)<=eps) exit
    end do
    cf=h
  end function beta_continued_fraction

  pure real(dp) function regularized_beta(x,a,b) result(value)
    real(dp), intent(in) :: x,a,b
    real(dp) :: bt
    if (a<=0.0_dp .or. b<=0.0_dp) then
      value=ieee_value(1.0_dp,ieee_quiet_nan)
    else if (x<=0.0_dp) then
      value=0.0_dp
    else if (x>=1.0_dp) then
      value=1.0_dp
    else
      bt=exp(log_gamma(a+b)-log_gamma(a)-log_gamma(b)+a*log(x)+b*log(1.0_dp-x))
      if (x < (a+1.0_dp)/(a+b+2.0_dp)) then
        value=bt*beta_continued_fraction(a,b,x)/a
      else
        value=1.0_dp-bt*beta_continued_fraction(b,a,1.0_dp-x)/b
      end if
      value=min(max(value,0.0_dp),1.0_dp)
    end if
  end function regularized_beta

  pure real(dp) function beta_cdf(x,a,b) result(value)
    real(dp), intent(in) :: x,a,b
    value=regularized_beta(x,a,b)
  end function beta_cdf

  pure real(dp) function gamma_cdf(x,shape,scale) result(value)
    real(dp), intent(in) :: x,shape,scale
    integer, parameter :: max_iter=300
    real(dp), parameter :: eps=3.0e-14_dp, fpmin=1.0e-300_dp
    integer :: n
    real(dp) :: ap,del,s,b,c,d,h,an,z,gln
    if (x<=0.0_dp) then
      value=0.0_dp; return
    end if
    if(shape<=0.0_dp .or. scale<=0.0_dp) then
      value=ieee_value(1.0_dp,ieee_quiet_nan); return
    end if
    z=x/scale; gln=log_gamma(shape)
    if(z<shape+1.0_dp) then
      ap=shape; del=1.0_dp/shape; s=del
      do n=1,max_iter
        ap=ap+1.0_dp; del=del*z/ap; s=s+del
        if(abs(del)<=abs(s)*eps) exit
      end do
      value=s*exp(-z+shape*log(z)-gln)
    else
      b=z+1.0_dp-shape; c=1.0_dp/fpmin; d=1.0_dp/b; h=d
      do n=1,max_iter
        an=-real(n,dp)*(real(n,dp)-shape); b=b+2.0_dp
        d=an*d+b; if(abs(d)<fpmin)d=fpmin
        c=b+an/c; if(abs(c)<fpmin)c=fpmin
        d=1.0_dp/d; del=d*c; h=h*del
        if(abs(del-1.0_dp)<=eps) exit
      end do
      value=1.0_dp-exp(-z+shape*log(z)-gln)*h
    end if
    value=min(max(value,0.0_dp),1.0_dp)
  end function gamma_cdf

  pure real(dp) function chi_square_cdf(x,df) result(value)
    real(dp), intent(in) :: x,df
    value=gamma_cdf(x,0.5_dp*df,2.0_dp)
  end function chi_square_cdf

  pure real(dp) function student_t_pdf(x,df) result(value)
    real(dp), intent(in) :: x,df
    if(df<=0.0_dp) then
      value=ieee_value(1.0_dp,ieee_quiet_nan)
    else
      value=exp(log_gamma((df+1.0_dp)/2.0_dp)-log_gamma(df/2.0_dp)-&
        0.5_dp*log(df*pi)-0.5_dp*(df+1.0_dp)*log(1.0_dp+x*x/df))
    end if
  end function student_t_pdf

  pure real(dp) function student_t_cdf(x,df) result(value)
    real(dp), intent(in) :: x,df
    real(dp) :: z,b
    if(df<=0.0_dp) then
      value=ieee_value(1.0_dp,ieee_quiet_nan); return
    end if
    if(abs(x)<=tiny(1.0_dp)) then
      value=0.5_dp; return
    end if
    z=df/(df+x*x); b=regularized_beta(z,df/2.0_dp,0.5_dp)
    if(x>0.0_dp) then
      value=1.0_dp-0.5_dp*b
    else
      value=0.5_dp*b
    end if
  end function student_t_cdf

  real(dp) function student_t_quantile(p,df) result(value)
    real(dp), intent(in) :: p,df
    real(dp) :: lo,hi,mid
    integer :: i
    if(p<=0.0_dp) then; value=-huge(1.0_dp); return; end if
    if(p>=1.0_dp) then; value=huge(1.0_dp); return; end if
    lo=-1.0_dp; hi=1.0_dp
    do while(student_t_cdf(lo,df)>p); lo=2.0_dp*lo; end do
    do while(student_t_cdf(hi,df)<p); hi=2.0_dp*hi; end do
    do i=1,120
      mid=0.5_dp*(lo+hi)
      if(student_t_cdf(mid,df)<p) then; lo=mid; else; hi=mid; end if
    end do
    value=0.5_dp*(lo+hi)
  end function student_t_quantile

  subroutine seed_random(seed)
    integer, intent(in) :: seed
    integer :: n,i
    integer, allocatable :: state(:)
    call random_seed(size=n); allocate(state(n))
    do i=1,n
      state(i)=modulo(seed+104729*i+37*i*i,huge(1)-1)
      if(state(i)==0) state(i)=i
    end do
    call random_seed(put=state)
  end subroutine seed_random

  real(dp) function random_uniform() result(value)
    call random_number(value)
    value=max(min(value,1.0_dp-epsilon(1.0_dp)),tiny(1.0_dp))
  end function random_uniform

  real(dp) function random_normal() result(value)
    real(dp) :: u1,u2
    u1=random_uniform(); u2=random_uniform()
    value=sqrt(-2.0_dp*log(u1))*cos(2.0_dp*pi*u2)
  end function random_normal

  pure real(dp) function mean_value(x) result(value)
    real(dp), intent(in) :: x(:)
    if(size(x)==0) then; value=0.0_dp; else; value=sum(x)/real(size(x),dp); end if
  end function mean_value

  pure real(dp) function variance_value(x) result(value)
    real(dp), intent(in) :: x(:)
    real(dp) :: m
    if(size(x)<=1) then; value=0.0_dp; else
      m=mean_value(x); value=sum((x-m)**2)/real(size(x)-1,dp)
    end if
  end function variance_value

  function covariance_matrix(x) result(covariance)
    real(dp), intent(in) :: x(:,:)
    real(dp), allocatable :: covariance(:,:)
    real(dp), allocatable :: means(:)
    integer :: i,j,n,d
    n=size(x,1); d=size(x,2); allocate(covariance(d,d),means(d))
    means=sum(x,dim=1)/real(n,dp)
    do j=1,d; do i=1,d
      covariance(i,j)=sum((x(:,i)-means(i))*(x(:,j)-means(j)))/real(max(n-1,1),dp)
    end do; end do
    covariance=0.5_dp*(covariance+transpose(covariance))
  end function covariance_matrix

  subroutine sort_increasing(x)
    real(dp), intent(inout) :: x(:)
    integer :: i,j
    real(dp) :: key
    do i=2,size(x)
      key=x(i); j=i-1
      do while(j>=1)
        if(x(j)<=key) exit
        x(j+1)=x(j); j=j-1
      end do
      x(j+1)=key
    end do
  end subroutine sort_increasing

  function stable_order(x,descending) result(order)
    real(dp), intent(in) :: x(:)
    logical, intent(in), optional :: descending
    integer, allocatable :: order(:)
    integer :: i,j,key
    logical :: desc
    desc=.false.; if(present(descending))desc=descending
    allocate(order(size(x))); order=[(i,i=1,size(x))]
    do i=2,size(x)
      key=order(i); j=i-1
      do while(j>=1)
        if((.not.desc .and. x(order(j))<=x(key)) .or. (desc .and. x(order(j))>=x(key))) exit
        order(j+1)=order(j); j=j-1
      end do
      order(j+1)=key
    end do
  end function stable_order

  real(dp) function quantile_type1(x,p) result(value)
    real(dp), intent(in) :: x(:),p
    real(dp), allocatable :: y(:)
    integer :: k
    y=x; call sort_increasing(y)
    if(p<=0.0_dp) then; value=y(1)
    else if(p>=1.0_dp) then; value=y(size(y))
    else; k=max(1,ceiling(p*real(size(y),dp))); value=y(k); end if
  end function quantile_type1

  subroutine solve_linear(a,b,x,ok)
    real(dp), intent(in) :: a(:,:),b(:)
    real(dp), allocatable, intent(out) :: x(:)
    logical, intent(out) :: ok
    real(dp), allocatable :: aug(:,:)
    real(dp) :: pivot,factor
    integer :: n,i,k,p
    n=size(b); allocate(aug(n,n+1),x(n)); aug(:,1:n)=a; aug(:,n+1)=b; ok=.true.
    do k=1,n
      p=k
      do i=k+1,n; if(abs(aug(i,k))>abs(aug(p,k)))p=i; end do
      if(abs(aug(p,k))<=1.0e-14_dp*max(1.0_dp,maxval(abs(a)))) then; ok=.false.; x=0.0_dp; return; end if
      if(p/=k) aug([k,p],:)=aug([p,k],:)
      pivot=aug(k,k); aug(k,k:n+1)=aug(k,k:n+1)/pivot
      do i=1,n
        if(i==k) cycle
        factor=aug(i,k); aug(i,k:n+1)=aug(i,k:n+1)-factor*aug(k,k:n+1)
      end do
    end do
    x=aug(:,n+1)
  end subroutine solve_linear

  subroutine invert_matrix(a,inverse,ok)
    real(dp), intent(in) :: a(:,:)
    real(dp), allocatable, intent(out) :: inverse(:,:)
    logical, intent(out) :: ok
    real(dp), allocatable :: e(:),col(:)
    integer :: n,j
    n=size(a,1); allocate(inverse(n,n),e(n)); ok=.true.
    do j=1,n
      e=0.0_dp; e(j)=1.0_dp; call solve_linear(a,e,col,ok)
      if(.not.ok) then; inverse=0.0_dp; return; end if
      inverse(:,j)=col
    end do
  end subroutine invert_matrix

  subroutine cholesky_factor(a,l,ok)
    real(dp), intent(in) :: a(:,:)
    real(dp), allocatable, intent(out) :: l(:,:)
    logical, intent(out) :: ok
    integer :: n,i,j,k
    real(dp) :: s
    n=size(a,1); allocate(l(n,n),source=0.0_dp); ok=.true.
    do i=1,n
      do j=1,i
        s=a(i,j)
        do k=1,j-1; s=s-l(i,k)*l(j,k); end do
        if(i==j) then
          if(s<=0.0_dp) then; ok=.false.; return; end if
          l(i,j)=sqrt(s)
        else
          l(i,j)=s/l(j,j)
        end if
      end do
    end do
  end subroutine cholesky_factor

  function numerical_gradient(f,x) result(g)
    procedure(vector_objective) :: f
    real(dp), intent(in) :: x(:)
    real(dp), allocatable :: g(:),xp(:),xm(:)
    real(dp) :: h
    integer :: i,n
    n=size(x); allocate(g(n),xp(n),xm(n))
    do i=1,n
      h=max(1.0e-6_dp,1.0e-5_dp*max(1.0_dp,abs(x(i))))
      xp=x; xm=x; xp(i)=xp(i)+h; xm(i)=xm(i)-h
      g(i)=(f(xp)-f(xm))/(2.0_dp*h)
    end do
  end function numerical_gradient

  function numerical_hessian(f,x) result(hess)
    procedure(vector_objective) :: f
    real(dp), intent(in) :: x(:)
    real(dp), allocatable :: hess(:,:),xpp(:),xpm(:),xmp(:),xmm(:),xp(:),xm(:)
    real(dp), allocatable :: step(:)
    real(dp) :: f0
    integer :: i,j,n
    n=size(x); allocate(hess(n,n),xpp(n),xpm(n),xmp(n),xmm(n),xp(n),xm(n),step(n))
    step=max(1.0e-5_dp,1.0e-4_dp*max(abs(x),1.0_dp)); f0=f(x)
    do i=1,n
      xp=x; xm=x; xp(i)=xp(i)+step(i); xm(i)=xm(i)-step(i)
      hess(i,i)=(f(xp)-2.0_dp*f0+f(xm))/step(i)**2
      do j=i+1,n
        xpp=x; xpm=x; xmp=x; xmm=x
        xpp(i)=xpp(i)+step(i); xpp(j)=xpp(j)+step(j)
        xpm(i)=xpm(i)+step(i); xpm(j)=xpm(j)-step(j)
        xmp(i)=xmp(i)-step(i); xmp(j)=xmp(j)+step(j)
        xmm(i)=xmm(i)-step(i); xmm(j)=xmm(j)-step(j)
        hess(i,j)=(f(xpp)-f(xpm)-f(xmp)+f(xmm))/(4.0_dp*step(i)*step(j)); hess(j,i)=hess(i,j)
      end do
    end do
  end function numerical_hessian

  real(dp) function trapz_integral(f,a,b,n) result(value)
    procedure(scalar_function) :: f
    real(dp), intent(in) :: a,b
    integer, intent(in), optional :: n
    integer :: m,i
    real(dp) :: h,s,x
    m=2048; if(present(n))m=max(n,2); h=(b-a)/real(m,dp); s=0.5_dp*(f(a)+f(b))
    do i=1,m-1; x=a+h*real(i,dp); s=s+f(x); end do
    value=h*s
  end function trapz_integral

  real(dp) function bisection_root(f,a,b,tolerance,ok) result(root)
    procedure(scalar_function) :: f
    real(dp), intent(in) :: a,b
    real(dp), intent(in), optional :: tolerance
    logical, intent(out), optional :: ok
    real(dp) :: lo,hi,mid,flo,fmid,tol
    integer :: i
    lo=a; hi=b; tol=1.0e-10_dp; if(present(tolerance))tol=tolerance
    flo=f(lo)
    if(flo*f(hi)>0.0_dp) then
      root=ieee_value(1.0_dp,ieee_quiet_nan); if(present(ok))ok=.false.; return
    end if
    do i=1,200
      mid=0.5_dp*(lo+hi); fmid=f(mid)
      if(abs(fmid)<tol .or. abs(hi-lo)<tol) exit
      if(flo*fmid<=0.0_dp) then; hi=mid; else; lo=mid; flo=fmid; end if
    end do
    root=0.5_dp*(lo+hi); if(present(ok))ok=.true.
  end function bisection_root
end module qrmtools_stats
