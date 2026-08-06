! SPDX-License-Identifier: GPL-3.0-only
module stochvoltmb_math
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_value, ieee_quiet_nan
  use stochvoltmb_kinds, only : dp, pi, log_two_pi, tiny_dp
  implicit none
  private
  public :: normal_pdf, normal_logpdf, normal_cdf, normal_logcdf, normal_quantile
  public :: student_t_logpdf, student_t_cdf, skew_normal_logpdf, skew_normal_cdf
  public :: owen_t, clamp, finite_real, quantile_type7
  public :: logit_pm1, inv_logit_pm1

contains

  elemental pure logical function finite_real(x) result(ok)
    real(dp), intent(in) :: x
    ok = ieee_is_finite(x)
  end function finite_real

  elemental pure real(dp) function clamp(x, lo, hi) result(y)
    real(dp), intent(in) :: x, lo, hi
    y = min(hi,max(lo,x))
  end function clamp

  pure real(dp) function normal_logpdf(x) result(v)
    real(dp), intent(in) :: x
    v = -0.5_dp*(log_two_pi+x*x)
  end function normal_logpdf

  pure real(dp) function normal_pdf(x) result(v)
    real(dp), intent(in) :: x
    v = exp(normal_logpdf(x))
  end function normal_pdf

  pure real(dp) function normal_cdf(x) result(v)
    real(dp), intent(in) :: x
    v = 0.5_dp*erfc(-x/sqrt(2.0_dp))
  end function normal_cdf

  pure real(dp) function normal_logcdf(x) result(v)
    real(dp), intent(in) :: x
    real(dp) :: z, term, sum, xx
    if (x > -8.0_dp) then
      v = log(max(tiny_dp,normal_cdf(x)))
    else
      xx = -x
      term = 1.0_dp
      sum = 1.0_dp
      term = term*(-1.0_dp)/(xx*xx); sum = sum+term
      term = term*(-3.0_dp)/(xx*xx); sum = sum+term
      term = term*(-5.0_dp)/(xx*xx); sum = sum+term
      z = max(tiny_dp,sum)
      v = normal_logpdf(x)-log(xx)+log(z)
    end if
  end function normal_logcdf

  real(dp) function normal_quantile(p) result(x)
    real(dp), intent(in) :: p
    real(dp), parameter :: a(6) = [ -3.969683028665376d1, 2.209460984245205d2, &
      -2.759285104469687d2, 1.383577518672690d2, -3.066479806614716d1, 2.506628277459239d0 ]
    real(dp), parameter :: b(5) = [ -5.447609879822406d1, 1.615858368580409d2, &
      -1.556989798598866d2, 6.680131188771972d1, -1.328068155288572d1 ]
    real(dp), parameter :: c(6) = [ -7.784894002430293d-3, -3.223964580411365d-1, &
      -2.400758277161838d0, -2.549732539343734d0, 4.374664141464968d0, 2.938163982698783d0 ]
    real(dp), parameter :: d(4) = [ 7.784695709041462d-3, 3.224671290700398d-1, &
      2.445134137142996d0, 3.754408661907416d0 ]
    real(dp), parameter :: plow=0.02425_dp, phigh=1.0_dp-plow
    real(dp) :: q, r, e, u
    if (p <= 0.0_dp) then
      x = -huge(1.0_dp)
      return
    else if (p >= 1.0_dp) then
      x = huge(1.0_dp)
      return
    end if
    if (p < plow) then
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
    e = normal_cdf(x)-p
    u = e/max(tiny_dp,normal_pdf(x))
    x = x-u/(1.0_dp+0.5_dp*x*u)
  end function normal_quantile

  pure real(dp) function student_t_logpdf(x, df) result(v)
    real(dp), intent(in) :: x, df
    v = log_gamma(0.5_dp*(df+1.0_dp))-log_gamma(0.5_dp*df) &
        -0.5_dp*log(df*pi)-0.5_dp*(df+1.0_dp)*log(1.0_dp+x*x/df)
  end function student_t_logpdf


  real(dp) function regularized_beta(x, a, b) result(v)
    real(dp), intent(in) :: x, a, b
    real(dp) :: bt
    if (x <= 0.0_dp) then
      v=0.0_dp
    else if (x >= 1.0_dp) then
      v=1.0_dp
    else
      bt=exp(log_gamma(a+b)-log_gamma(a)-log_gamma(b)+a*log(x)+b*log(1.0_dp-x))
      if (x < (a+1.0_dp)/(a+b+2.0_dp)) then
        v=bt*beta_cf(x,a,b)/a
      else
        v=1.0_dp-bt*beta_cf(1.0_dp-x,b,a)/b
      end if
      v=clamp(v,0.0_dp,1.0_dp)
    end if
  contains
    real(dp) function beta_cf(xx,aa,bb) result(cf)
      real(dp), intent(in) :: xx,aa,bb
      integer, parameter :: maxit=300
      real(dp), parameter :: eps=3.0e-14_dp, fpmin=1.0e-300_dp
      real(dp) :: qab,qap,qam,c,d,h,del,an
      integer :: m,m2
      qab=aa+bb; qap=aa+1.0_dp; qam=aa-1.0_dp
      c=1.0_dp; d=1.0_dp-qab*xx/qap
      if (abs(d)<fpmin) d=fpmin
      d=1.0_dp/d; h=d
      do m=1,maxit
        m2=2*m
        an=real(m,dp)*(bb-real(m,dp))*xx/((qam+real(m2,dp))*(aa+real(m2,dp)))
        d=1.0_dp+an*d; if (abs(d)<fpmin) d=fpmin
        c=1.0_dp+an/c; if (abs(c)<fpmin) c=fpmin
        d=1.0_dp/d; h=h*d*c
        an=-(aa+real(m,dp))*(qab+real(m,dp))*xx/((aa+real(m2,dp))*(qap+real(m2,dp)))
        d=1.0_dp+an*d; if (abs(d)<fpmin) d=fpmin
        c=1.0_dp+an/c; if (abs(c)<fpmin) c=fpmin
        d=1.0_dp/d; del=d*c; h=h*del
        if (abs(del-1.0_dp)<eps) exit
      end do
      cf=h
    end function beta_cf
  end function regularized_beta

  real(dp) function student_t_cdf(x, df) result(v)
    real(dp), intent(in) :: x, df
    real(dp) :: ib
    if (df<=0.0_dp) then
      v=ieee_value(0.0_dp,ieee_quiet_nan)
      return
    end if
    if (abs(x)<=tiny_dp) then
      v=0.5_dp
      return
    end if
    ib=regularized_beta(df/(df+x*x),0.5_dp*df,0.5_dp)
    if (x>0.0_dp) then
      v=1.0_dp-0.5_dp*ib
    else
      v=0.5_dp*ib
    end if
  end function student_t_cdf

  real(dp) function owen_t(h, a) result(value)
    real(dp), intent(in) :: h, a
    real(dp) :: aa, ah
    aa=abs(a); ah=abs(h)
    if (aa<=tiny_dp) then
      value=0.0_dp
    else if (aa<=1.0_dp) then
      value=sign(t_int(ah,aa),a)
    else
      value=sign(0.5_dp*normal_cdf(ah)+normal_cdf(aa*ah)*(0.5_dp-normal_cdf(ah)) &
                 -t_int(aa*ah,1.0_dp/aa),a)
    end if
  contains
    real(dp) function t_int(x,b) result(ans)
      real(dp), intent(in) :: x,b
      integer :: i
      real(dp) :: term,cumulative,alternating_sum,powb,factorial
      cumulative=0.0_dp; alternating_sum=0.0_dp; powb=b; factorial=1.0_dp
      do i=0,60
        if (i>0) factorial=factorial*real(i,dp)
        term=x**(2*i)/(2.0_dp**i*factorial)
        cumulative=cumulative+term
        alternating_sum=alternating_sum+(-1.0_dp)**i* &
          (1.0_dp-exp(-0.5_dp*x*x)*cumulative)*powb/real(2*i+1,dp)
        powb=powb*b*b
        if (abs(term*powb)<1.0e-16_dp) exit
      end do
      ans=(atan(b)-alternating_sum)/(2.0_dp*pi)
    end function t_int
  end function owen_t

  real(dp) function skew_normal_cdf(x, xi, omega, alpha) result(v)
    real(dp), intent(in) :: x,xi,omega,alpha
    real(dp) :: z
    if (omega<=0.0_dp) then
      v=ieee_value(0.0_dp,ieee_quiet_nan)
    else
      z=(x-xi)/omega
      v=clamp(normal_cdf(z)-2.0_dp*owen_t(z,alpha),0.0_dp,1.0_dp)
    end if
  end function skew_normal_cdf

  pure real(dp) function skew_normal_logpdf(x, xi, omega, alpha) result(v)
    real(dp), intent(in) :: x, xi, omega, alpha
    real(dp) :: z
    if (omega <= 0.0_dp) then
      v = -huge(1.0_dp)
      return
    end if
    z = (x-xi)/omega
    v = log(2.0_dp)-log(omega)+normal_logpdf(z)+normal_logcdf(alpha*z)
  end function skew_normal_logpdf

  elemental pure real(dp) function inv_logit_pm1(x) result(y)
    real(dp), intent(in) :: x
    y = tanh(0.5_dp*x)
  end function inv_logit_pm1

  elemental pure real(dp) function logit_pm1(y) result(x)
    real(dp), intent(in) :: y
    real(dp) :: z
    z = clamp(y,-1.0_dp+1.0e-12_dp,1.0_dp-1.0e-12_dp)
    x = log((1.0_dp+z)/(1.0_dp-z))
  end function logit_pm1

  real(dp) function quantile_type7(x, p) result(q)
    real(dp), intent(in) :: x(:), p
    real(dp), allocatable :: z(:)
    real(dp) :: h, g
    integer :: n, j
    n = size(x)
    if (n < 1 .or. p < 0.0_dp .or. p > 1.0_dp) then
      q = ieee_value(0.0_dp,ieee_quiet_nan)
      return
    end if
    allocate(z(n)); z=x
    call sort_real(z)
    if (p <= 0.0_dp) then
      q=z(1)
    else if (p >= 1.0_dp) then
      q=z(n)
    else
      h = 1.0_dp+real(n-1,dp)*p
      j = floor(h)
      g = h-real(j,dp)
      q = (1.0_dp-g)*z(j)+g*z(j+1)
    end if
  contains
    subroutine sort_real(a)
      real(dp), intent(inout) :: a(:)
      integer :: i, k
      real(dp) :: key
      do i=2,size(a)
        key=a(i); k=i-1
        do while (k>=1)
          if (a(k)<=key) exit
          a(k+1)=a(k); k=k-1
        end do
        a(k+1)=key
      end do
    end subroutine sort_real
  end function quantile_type7

end module stochvoltmb_math
