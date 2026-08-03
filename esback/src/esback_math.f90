! SPDX-License-Identifier: GPL-3.0-only
module esback_math
  use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
  use esback_kinds, only: dp, pi, sqrt_two
  use esback_types, only: rng_state, esback_ok, esback_singular
  implicit none
  private
  public :: normal_pdf, normal_cdf, normal_quantile, chi_square_survival
  public :: sample_mean, sample_variance, sample_sd, empirical_quantile_type7
  public :: rng_seed, rng_uniform, rng_normal, rng_index
  public :: solve_linear, invert_matrix, least_squares, outer_product
  public :: sort_real, median_real, interpolate_linear, empirical_cdf_values

  interface
    subroutine dgesv(n, nrhs, a, lda, ipiv, b, ldb, info)
      import :: dp
      integer, intent(in) :: n, nrhs, lda, ldb
      integer, intent(out) :: ipiv(*)
      integer, intent(out) :: info
      real(dp), intent(inout) :: a(lda,*)
      real(dp), intent(inout) :: b(ldb,*)
    end subroutine dgesv
    subroutine dgetrf(m, n, a, lda, ipiv, info)
      import :: dp
      integer, intent(in) :: m, n, lda
      integer, intent(out) :: ipiv(*)
      integer, intent(out) :: info
      real(dp), intent(inout) :: a(lda,*)
    end subroutine dgetrf
    subroutine dgetri(n, a, lda, ipiv, work, lwork, info)
      import :: dp
      integer, intent(in) :: n, lda, lwork
      integer, intent(in) :: ipiv(*)
      integer, intent(out) :: info
      real(dp), intent(inout) :: a(lda,*)
      real(dp), intent(out) :: work(*)
    end subroutine dgetri
  end interface
contains
  pure elemental real(dp) function normal_pdf(x) result(y)
    real(dp), intent(in) :: x
    y = exp(-0.5_dp*x*x)/sqrt(2.0_dp*pi)
  end function normal_pdf

  pure elemental real(dp) function normal_cdf(x) result(y)
    real(dp), intent(in) :: x
    y = 0.5_dp*erfc(-x/sqrt_two)
  end function normal_cdf

  pure elemental real(dp) function normal_quantile(p) result(x)
    real(dp), intent(in) :: p
    real(dp), parameter :: a(6) = [ &
      -3.969683028665376e1_dp, 2.209460984245205e2_dp, -2.759285104469687e2_dp, &
       1.383577518672690e2_dp, -3.066479806614716e1_dp, 2.506628277459239_dp]
    real(dp), parameter :: b(5) = [ &
      -5.447609879822406e1_dp, 1.615858368580409e2_dp, -1.556989798598866e2_dp, &
       6.680131188771972e1_dp, -1.328068155288572e1_dp]
    real(dp), parameter :: c(6) = [ &
      -7.784894002430293e-3_dp, -3.223964580411365e-1_dp, -2.400758277161838_dp, &
      -2.549732539343734_dp, 4.374664141464968_dp, 2.938163982698783_dp]
    real(dp), parameter :: d(4) = [ &
       7.784695709041462e-3_dp, 3.224671290700398e-1_dp, 2.445134137142996_dp, &
       3.754408661907416_dp]
    real(dp), parameter :: plow=0.02425_dp, phigh=1.0_dp-plow
    real(dp) :: q, r
    if (p <= 0.0_dp) then
      x = -huge(1.0_dp)
    else if (p >= 1.0_dp) then
      x = huge(1.0_dp)
    else if (p < plow) then
      q = sqrt(-2.0_dp*log(p))
      x = (((((c(1)*q+c(2))*q+c(3))*q+c(4))*q+c(5))*q+c(6))/ &
          ((((d(1)*q+d(2))*q+d(3))*q+d(4))*q+1.0_dp)
    else if (p <= phigh) then
      q = p-0.5_dp
      r = q*q
      x = (((((a(1)*r+a(2))*r+a(3))*r+a(4))*r+a(5))*r+a(6))*q/ &
          (((((b(1)*r+b(2))*r+b(3))*r+b(4))*r+b(5))*r+1.0_dp)
    else
      q = sqrt(-2.0_dp*log(1.0_dp-p))
      x = -(((((c(1)*q+c(2))*q+c(3))*q+c(4))*q+c(5))*q+c(6))/ &
           ((((d(1)*q+d(2))*q+d(3))*q+d(4))*q+1.0_dp)
    end if
    if (p > 0.0_dp .and. p < 1.0_dp) then
      x = x - (normal_cdf(x)-p)/max(normal_pdf(x), tiny(1.0_dp))
    end if
  end function normal_quantile

  pure real(dp) function chi_square_survival(x, df) result(p)
    real(dp), intent(in) :: x
    integer, intent(in) :: df
    if (x <= 0.0_dp) then
      p = 1.0_dp
    else if (df <= 0) then
      p = ieee_value(1.0_dp, ieee_quiet_nan)
    else
      p = regularized_gamma_q(0.5_dp*real(df,dp), 0.5_dp*x)
    end if
  end function chi_square_survival

  pure real(dp) function regularized_gamma_q(a, x) result(q)
    real(dp), intent(in) :: a, x
    integer, parameter :: maxit=10000
    real(dp), parameter :: eps=2.0e-15_dp, fpmin=1.0e-300_dp
    integer :: i
    real(dp) :: ap, del, sumv, b, c, d, h, an
    if (x < 0.0_dp .or. a <= 0.0_dp) then
      q = ieee_value(1.0_dp, ieee_quiet_nan)
      return
    end if
    if (x <= tiny(1.0_dp)) then
      q = 1.0_dp
      return
    end if
    if (x < a + 1.0_dp) then
      ap = a
      sumv = 1.0_dp/a
      del = sumv
      do i=1,maxit
        ap = ap + 1.0_dp
        del = del*x/ap
        sumv = sumv + del
        if (abs(del) <= abs(sumv)*eps) exit
      end do
      q = 1.0_dp - sumv*exp(-x+a*log(x)-log_gamma(a))
    else
      b = x + 1.0_dp - a
      c = 1.0_dp/fpmin
      d = 1.0_dp/b
      h = d
      do i=1,maxit
        an = -real(i,dp)*(real(i,dp)-a)
        b = b + 2.0_dp
        d = an*d + b
        if (abs(d) < fpmin) d = fpmin
        c = b + an/c
        if (abs(c) < fpmin) c = fpmin
        d = 1.0_dp/d
        del = d*c
        h = h*del
        if (abs(del-1.0_dp) <= eps) exit
      end do
      q = exp(-x+a*log(x)-log_gamma(a))*h
    end if
    q = min(1.0_dp,max(0.0_dp,q))
  end function regularized_gamma_q

  pure real(dp) function sample_mean(x) result(m)
    real(dp), intent(in) :: x(:)
    if (size(x)==0) then
      m = ieee_value(1.0_dp, ieee_quiet_nan)
    else
      m = sum(x)/real(size(x),dp)
    end if
  end function sample_mean

  pure real(dp) function sample_variance(x) result(v)
    real(dp), intent(in) :: x(:)
    real(dp) :: m
    if (size(x)<2) then
      v = ieee_value(1.0_dp, ieee_quiet_nan)
    else
      m = sample_mean(x)
      v = sum((x-m)**2)/real(size(x)-1,dp)
    end if
  end function sample_variance

  pure real(dp) function sample_sd(x) result(s)
    real(dp), intent(in) :: x(:)
    s = sqrt(max(0.0_dp,sample_variance(x)))
  end function sample_sd

  function empirical_quantile_type7(x,p) result(q)
    real(dp), intent(in) :: x(:), p
    real(dp) :: q
    real(dp), allocatable :: y(:)
    real(dp) :: h, g
    integer :: j, n
    n=size(x)
    if(n==0 .or. p<0.0_dp .or. p>1.0_dp) then
      q=ieee_value(1.0_dp,ieee_quiet_nan); return
    end if
    allocate(y(n)); y=x; call sort_real(y)
    if(n==1) then; q=y(1); return; end if
    h=1.0_dp+real(n-1,dp)*p
    j=floor(h); g=h-real(j,dp)
    if(j>=n) then; q=y(n)
    else; q=(1.0_dp-g)*y(max(1,j))+g*y(j+1)
    end if
  end function empirical_quantile_type7

  subroutine rng_seed(rng, seed)
    type(rng_state), intent(out) :: rng
    integer(kind=8), intent(in) :: seed
    if(seed==0_8) then; rng%state=88172645463393265_8; else; rng%state=abs(seed); end if
  end subroutine rng_seed

  real(dp) function rng_uniform(rng) result(u)
    type(rng_state), intent(inout) :: rng
    integer(kind=8) :: x
    x=rng%state
    x=ieor(x,shiftl(x,13)); x=ieor(x,shiftr(x,7)); x=ieor(x,shiftl(x,17))
    rng%state=x
    u=real(iand(x,int(z'7FFFFFFFFFFFFFFF',kind=8)),dp)/real(huge(1_8),dp)
    if(u<=0.0_dp) u=epsilon(1.0_dp)
    if(u>=1.0_dp) u=1.0_dp-epsilon(1.0_dp)
  end function rng_uniform

  real(dp) function rng_normal(rng) result(z)
    type(rng_state), intent(inout) :: rng
    real(dp) :: u1,u2
    u1=rng_uniform(rng); u2=rng_uniform(rng)
    z=sqrt(-2.0_dp*log(u1))*cos(2.0_dp*pi*u2)
  end function rng_normal

  integer function rng_index(rng,n) result(i)
    type(rng_state), intent(inout) :: rng
    integer,intent(in)::n
    i=1+int(rng_uniform(rng)*real(n,dp)); if(i>n)i=n
  end function rng_index

  subroutine solve_linear(a,b,x,status)
    real(dp),intent(in)::a(:,:),b(:)
    real(dp),allocatable,intent(out)::x(:)
    integer,intent(out)::status
    real(dp),allocatable::aa(:,:),bb(:,:)
    integer,allocatable::ipiv(:)
    integer::n,info
    n=size(b); allocate(aa(n,n),bb(n,1),ipiv(n),x(n)); aa=a; bb(:,1)=b
    call dgesv(n,1,aa,n,ipiv,bb,n,info)
    if(info==0) then; x=bb(:,1); status=esback_ok
    else; x=0.0_dp; status=esback_singular; end if
  end subroutine solve_linear

  subroutine invert_matrix(a,ainv,status)
    real(dp),intent(in)::a(:,:)
    real(dp),allocatable,intent(out)::ainv(:,:)
    integer,intent(out)::status
    integer::n,info,lwork
    integer,allocatable::ipiv(:)
    real(dp),allocatable::work(:)
    n=size(a,1); allocate(ainv(n,n),ipiv(n)); ainv=a
    call dgetrf(n,n,ainv,n,ipiv,info)
    if(info/=0) then; status=esback_singular; ainv=0.0_dp; return; end if
    lwork=max(1,64*n); allocate(work(lwork))
    call dgetri(n,ainv,n,ipiv,work,lwork,info)
    if(info==0) then; status=esback_ok; else; status=esback_singular; ainv=0.0_dp; end if
  end subroutine invert_matrix

  subroutine least_squares(x,y,b,status,ridge)
    real(dp),intent(in)::x(:,:),y(:)
    real(dp),allocatable,intent(out)::b(:)
    integer,intent(out)::status
    real(dp),intent(in),optional::ridge
    real(dp),allocatable::xtx(:,:),xty(:)
    real(dp)::r
    integer::k,i
    k=size(x,2); allocate(xtx(k,k),xty(k)); xtx=matmul(transpose(x),x); xty=matmul(transpose(x),y)
    r=1.0e-10_dp; if(present(ridge))r=ridge
    do i=1,k; xtx(i,i)=xtx(i,i)+r; end do
    call solve_linear(xtx,xty,b,status)
  end subroutine least_squares

  pure function outer_product(a,b) result(c)
    real(dp),intent(in)::a(:),b(:)
    real(dp)::c(size(a),size(b))
    integer::i
    do i=1,size(a); c(i,:)=a(i)*b; end do
  end function outer_product

  recursive subroutine quicksort(a,lo,hi)
    real(dp),intent(inout)::a(:); integer,intent(in)::lo,hi
    integer::i,j; real(dp)::p,t
    if(lo>=hi)return
    p=a((lo+hi)/2); i=lo; j=hi
    do
      do while(a(i)<p); i=i+1; end do
      do while(a(j)>p); j=j-1; end do
      if(i<=j) then; t=a(i);a(i)=a(j);a(j)=t;i=i+1;j=j-1;end if
      if(i>j)exit
    end do
    if(lo<j)call quicksort(a,lo,j); if(i<hi)call quicksort(a,i,hi)
  end subroutine quicksort

  subroutine sort_real(a)
    real(dp),intent(inout)::a(:)
    if(size(a)>1)call quicksort(a,1,size(a))
  end subroutine sort_real

  function median_real(x) result(m)
    real(dp),intent(in)::x(:); real(dp)::m
    real(dp),allocatable::y(:);integer::n
    n=size(x);allocate(y(n));y=x;call sort_real(y)
    if(mod(n,2)==1)then;m=y((n+1)/2);else;m=0.5_dp*(y(n/2)+y(n/2+1));end if
  end function median_real

  pure real(dp) function interpolate_linear(x,y,x0) result(y0)
    real(dp),intent(in)::x(:),y(:),x0
    integer::i,n
    n=size(x)
    if(x0<=x(1))then;y0=y(1);return;else if(x0>=x(n))then;y0=y(n);return;end if
    do i=1,n-1
      if(x0<=x(i+1))then;y0=y(i)+(y(i+1)-y(i))*(x0-x(i))/(x(i+1)-x(i));return;end if
    end do
    y0=y(n)
  end function interpolate_linear

  subroutine empirical_cdf_values(sample,z,out)
    real(dp),intent(in)::sample(:),z(:)
    real(dp),intent(out)::out(size(z))
    integer::i
    do i=1,size(z); out(i)=real(count(sample<=z(i)),dp)/real(size(sample),dp); end do
  end subroutine empirical_cdf_values
end module esback_math
