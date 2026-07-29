! SPDX-License-Identifier: MIT
module bekks_math
  use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
  use bekks_kinds, only: dp, pi
  implicit none
  private
  public :: normal_cdf, normal_quantile, student_t_cdf, student_t_quantile
  public :: chi_square_cdf, log_gamma_dp, empirical_quantile, sample_kurtosis

contains

  pure elemental real(dp) function normal_cdf(x) result(p)
    real(dp), intent(in) :: x
    p = 0.5_dp * erfc(-x / sqrt(2.0_dp))
  end function normal_cdf

  pure elemental real(dp) function normal_quantile(p) result(x)
    real(dp), intent(in) :: p
    real(dp), parameter :: a(6) = [ &
      -3.969683028665376e+01_dp, 2.209460984245205e+02_dp, &
      -2.759285104469687e+02_dp, 1.383577518672690e+02_dp, &
      -3.066479806614716e+01_dp, 2.506628277459239e+00_dp]
    real(dp), parameter :: b(5) = [ &
      -5.447609879822406e+01_dp, 1.615858368580409e+02_dp, &
      -1.556989798598866e+02_dp, 6.680131188771972e+01_dp, &
      -1.328068155288572e+01_dp]
    real(dp), parameter :: c(6) = [ &
      -7.784894002430293e-03_dp, -3.223964580411365e-01_dp, &
      -2.400758277161838e+00_dp, -2.549732539343734e+00_dp, &
       4.374664141464968e+00_dp,  2.938163982698783e+00_dp]
    real(dp), parameter :: d(4) = [ &
       7.784695709041462e-03_dp, 3.224671290700398e-01_dp, &
       2.445134137142996e+00_dp, 3.754408661907416e+00_dp]
    real(dp), parameter :: plow = 0.02425_dp, phigh = 1.0_dp - plow
    real(dp) :: q, r
    if (p <= 0.0_dp) then
      x = -huge(1.0_dp)
    else if (p >= 1.0_dp) then
      x = huge(1.0_dp)
    else if (p < plow) then
      q = sqrt(-2.0_dp * log(p))
      x = (((((c(1)*q+c(2))*q+c(3))*q+c(4))*q+c(5))*q+c(6)) / &
          ((((d(1)*q+d(2))*q+d(3))*q+d(4))*q+1.0_dp)
    else if (p <= phigh) then
      q = p - 0.5_dp
      r = q*q
      x = (((((a(1)*r+a(2))*r+a(3))*r+a(4))*r+a(5))*r+a(6))*q / &
          (((((b(1)*r+b(2))*r+b(3))*r+b(4))*r+b(5))*r+1.0_dp)
    else
      q = sqrt(-2.0_dp * log(1.0_dp-p))
      x = -(((((c(1)*q+c(2))*q+c(3))*q+c(4))*q+c(5))*q+c(6)) / &
           ((((d(1)*q+d(2))*q+d(3))*q+d(4))*q+1.0_dp)
    end if
  end function normal_quantile

  pure real(dp) function log_gamma_dp(x) result(v)
    real(dp), intent(in) :: x
    v = log_gamma(x)
  end function log_gamma_dp

  pure real(dp) function betacf(a,b,x) result(cf)
    real(dp), intent(in) :: a,b,x
    integer, parameter :: maxit=300
    real(dp), parameter :: eps=3.0e-14_dp, fpmin=1.0e-300_dp
    integer :: m, m2
    real(dp) :: aa,c,d,del,h,qab,qam,qap
    qab=a+b; qap=a+1.0_dp; qam=a-1.0_dp
    c=1.0_dp; d=1.0_dp-qab*x/qap
    if (abs(d)<fpmin) d=fpmin
    d=1.0_dp/d; h=d
    do m=1,maxit
      m2=2*m
      aa=m*(b-m)*x/((qam+m2)*(a+m2))
      d=1.0_dp+aa*d; if(abs(d)<fpmin)d=fpmin
      c=1.0_dp+aa/c; if(abs(c)<fpmin)c=fpmin
      d=1.0_dp/d; h=h*d*c
      aa=-(a+m)*(qab+m)*x/((a+m2)*(qap+m2))
      d=1.0_dp+aa*d; if(abs(d)<fpmin)d=fpmin
      c=1.0_dp+aa/c; if(abs(c)<fpmin)c=fpmin
      d=1.0_dp/d; del=d*c; h=h*del
      if(abs(del-1.0_dp)<=eps) exit
    end do
    cf=h
  end function betacf

  pure real(dp) function regularized_beta(x,a,b) result(v)
    real(dp), intent(in) :: x,a,b
    real(dp) :: bt
    if (x <= 0.0_dp) then
      v=0.0_dp; return
    else if (x >= 1.0_dp) then
      v=1.0_dp; return
    end if
    bt=exp(log_gamma(a+b)-log_gamma(a)-log_gamma(b)+a*log(x)+b*log(1.0_dp-x))
    if (x < (a+1.0_dp)/(a+b+2.0_dp)) then
      v=bt*betacf(a,b,x)/a
    else
      v=1.0_dp-bt*betacf(b,a,1.0_dp-x)/b
    end if
  end function regularized_beta

  pure real(dp) function student_t_cdf(x,nu) result(p)
    real(dp), intent(in) :: x,nu
    real(dp) :: z, ib
    if (nu <= 0.0_dp) then
      p = ieee_value(0.0_dp,ieee_quiet_nan)
      return
    end if
    z=nu/(nu+x*x)
    ib=regularized_beta(z,0.5_dp*nu,0.5_dp)
    if (x >= 0.0_dp) then
      p=1.0_dp-0.5_dp*ib
    else
      p=0.5_dp*ib
    end if
  end function student_t_cdf

  real(dp) function student_t_quantile(p,nu) result(x)
    real(dp), intent(in) :: p,nu
    real(dp) :: lo,hi,mid
    integer :: i
    if (p <= 0.0_dp) then
      x=-huge(1.0_dp); return
    else if (p >= 1.0_dp) then
      x=huge(1.0_dp); return
    end if
    lo=-50.0_dp; hi=50.0_dp
    do i=1,200
      mid=0.5_dp*(lo+hi)
      if(student_t_cdf(mid,nu)<p) then
        lo=mid
      else
        hi=mid
      end if
      if(hi-lo<1.0e-12_dp) exit
    end do
    x=0.5_dp*(lo+hi)
  end function student_t_quantile

  pure real(dp) function gamma_p(a,x) result(p)
    real(dp), intent(in) :: a,x
    integer, parameter :: itmax=300
    real(dp), parameter :: eps=3.0e-14_dp, fpmin=1.0e-300_dp
    integer :: n
    real(dp) :: ap,del,sumg,b,c,d,h,an
    if(x<=0.0_dp) then
      p=0.0_dp; return
    end if
    if(x<a+1.0_dp) then
      ap=a; sumg=1.0_dp/a; del=sumg
      do n=1,itmax
        ap=ap+1.0_dp; del=del*x/ap; sumg=sumg+del
        if(abs(del)<abs(sumg)*eps) exit
      end do
      p=sumg*exp(-x+a*log(x)-log_gamma(a))
    else
      b=x+1.0_dp-a; c=1.0_dp/fpmin; d=1.0_dp/b; h=d
      do n=1,itmax
        an=-real(n,dp)*(real(n,dp)-a)
        b=b+2.0_dp
        d=an*d+b; if(abs(d)<fpmin)d=fpmin
        c=b+an/c; if(abs(c)<fpmin)c=fpmin
        d=1.0_dp/d; del=d*c; h=h*del
        if(abs(del-1.0_dp)<eps) exit
      end do
      p=1.0_dp-exp(-x+a*log(x)-log_gamma(a))*h
    end if
    p=max(0.0_dp,min(1.0_dp,p))
  end function gamma_p

  pure real(dp) function chi_square_cdf(x,df) result(p)
    real(dp), intent(in) :: x
    integer, intent(in) :: df
    if (df <= 0) then
      p=ieee_value(0.0_dp,ieee_quiet_nan)
    else
      p=gamma_p(0.5_dp*real(df,dp),0.5_dp*x)
    end if
  end function chi_square_cdf

  real(dp) function empirical_quantile(x,p) result(q)
    real(dp), intent(in) :: x(:),p
    real(dp), allocatable :: y(:)
    real(dp) :: h, frac
    integer :: n,j
    n=size(x)
    if(n==0) then
      q=ieee_value(0.0_dp,ieee_quiet_nan); return
    end if
    allocate(y(n)); y=x; call sort_in_place(y)
    if(p<=0.0_dp) then
      q=y(1)
    else if(p>=1.0_dp) then
      q=y(n)
    else
      h=1.0_dp+real(n-1,dp)*p
      j=int(floor(h)); frac=h-real(j,dp)
      if(j>=n) then
        q=y(n)
      else
        q=(1.0_dp-frac)*y(j)+frac*y(j+1)
      end if
    end if
  end function empirical_quantile

  subroutine sort_in_place(x)
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
  end subroutine sort_in_place

  real(dp) function sample_kurtosis(x) result(k)
    real(dp), intent(in) :: x(:)
    real(dp) :: m,v
    integer :: n
    n=size(x)
    if(n<4) then
      k=3.0_dp; return
    end if
    m=sum(x)/real(n,dp)
    v=sum((x-m)**2)/real(n,dp)
    if(v<=tiny(1.0_dp)) then
      k=3.0_dp
    else
      k=sum((x-m)**4)/real(n,dp)/(v*v)
    end if
  end function sample_kurtosis

end module bekks_math
