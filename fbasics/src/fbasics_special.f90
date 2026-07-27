! SPDX-License-Identifier: GPL-2.0-or-later
! Numerical translation derived from the GPL-2-or-later fBasics package.
module fbasics_special
  use fbasics_kinds, only: dp, pi, sqrt2, clamp
  implicit none
  private
  public :: normal_pdf, normal_cdf, normal_quantile
  public :: regularized_gamma_p, regularized_gamma_q
  public :: regularized_beta, student_pdf, student_cdf, student_quantile
  public :: chi_square_cdf, chi_square_quantile, f_cdf
  public :: bessel_k_nu, adaptive_simpson
  abstract interface
    function scalar_function(x) result(y)
      import dp
      real(dp), intent(in) :: x
      real(dp) :: y
    end function scalar_function
  end interface
contains
  pure elemental real(dp) function normal_pdf(x) result(f)
    real(dp), intent(in) :: x
    f = exp(-0.5_dp*x*x)/sqrt(2.0_dp*pi)
  end function normal_pdf

  pure elemental real(dp) function normal_cdf(x) result(p)
    real(dp), intent(in) :: x
    p = 0.5_dp*erfc(-x/sqrt2)
  end function normal_cdf

  pure elemental real(dp) function normal_quantile(p) result(x)
    real(dp), intent(in) :: p
    real(dp), parameter :: a(6) = [ -3.969683028665376d1, 2.209460984245205d2, &
      -2.759285104469687d2, 1.383577518672690d2, -3.066479806614716d1, 2.506628277459239d0 ]
    real(dp), parameter :: b(5) = [ -5.447609879822406d1, 1.615858368580409d2, &
      -1.556989798598866d2, 6.680131188771972d1, -1.328068155288572d1 ]
    real(dp), parameter :: c(6) = [ -7.784894002430293d-3, -3.223964580411365d-1, &
      -2.400758277161838d0, -2.549732539343734d0, 4.374664141464968d0, 2.938163982698783d0 ]
    real(dp), parameter :: d(4) = [ 7.784695709041462d-3, 3.224671290700398d-1, &
      2.445134137142996d0, 3.754408661907416d0 ]
    real(dp), parameter :: plow = 0.02425_dp, phigh = 1.0_dp - plow
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
      q = p - 0.5_dp
      r = q*q
      x = (((((a(1)*r+a(2))*r+a(3))*r+a(4))*r+a(5))*r+a(6))*q/ &
          (((((b(1)*r+b(2))*r+b(3))*r+b(4))*r+b(5))*r+1.0_dp)
    else
      q = sqrt(-2.0_dp*log(1.0_dp-p))
      x = -(((((c(1)*q+c(2))*q+c(3))*q+c(4))*q+c(5))*q+c(6))/ &
          ((((d(1)*q+d(2))*q+d(3))*q+d(4))*q+1.0_dp)
    end if
  end function normal_quantile

  real(dp) function regularized_gamma_p(a, x) result(p)
    real(dp), intent(in) :: a, x
    integer, parameter :: itmax = 10000
    real(dp), parameter :: eps = 2.0e-14_dp, fpmin = 1.0e-300_dp
    integer :: n
    real(dp) :: ap, del, sumv, b, c, d, h, an
    if (a <= 0.0_dp .or. x < 0.0_dp) then
      p = 0.0_dp
      return
    end if
    if (x == 0.0_dp) then
      p = 0.0_dp
    else if (x < a + 1.0_dp) then
      ap = a
      sumv = 1.0_dp/a
      del = sumv
      do n = 1, itmax
        ap = ap + 1.0_dp
        del = del*x/ap
        sumv = sumv + del
        if (abs(del) < abs(sumv)*eps) exit
      end do
      p = sumv*exp(-x + a*log(x) - log_gamma(a))
    else
      b = x + 1.0_dp - a
      c = 1.0_dp/fpmin
      d = 1.0_dp/b
      h = d
      do n = 1, itmax
        an = -real(n,dp)*(real(n,dp)-a)
        b = b + 2.0_dp
        d = an*d + b
        if (abs(d) < fpmin) d = fpmin
        c = b + an/c
        if (abs(c) < fpmin) c = fpmin
        d = 1.0_dp/d
        del = d*c
        h = h*del
        if (abs(del-1.0_dp) < eps) exit
      end do
      p = 1.0_dp - exp(-x + a*log(x) - log_gamma(a))*h
    end if
    p = clamp(p, 0.0_dp, 1.0_dp)
  end function regularized_gamma_p

  real(dp) function regularized_gamma_q(a, x) result(q)
    real(dp), intent(in) :: a, x
    q = 1.0_dp - regularized_gamma_p(a, x)
  end function regularized_gamma_q

  real(dp) function beta_cf(a, b, x) result(cf)
    real(dp), intent(in) :: a, b, x
    integer, parameter :: maxit = 10000
    real(dp), parameter :: eps = 3.0e-14_dp, fpmin = 1.0e-300_dp
    integer :: m, m2
    real(dp) :: qab, qap, qam, c, d, h, aa, del
    qab = a+b; qap = a+1.0_dp; qam = a-1.0_dp
    c = 1.0_dp
    d = 1.0_dp - qab*x/qap
    if (abs(d) < fpmin) d = fpmin
    d = 1.0_dp/d
    h = d
    do m = 1, maxit
      m2 = 2*m
      aa = real(m,dp)*(b-real(m,dp))*x/((qam+real(m2,dp))*(a+real(m2,dp)))
      d = 1.0_dp + aa*d; if (abs(d) < fpmin) d = fpmin
      c = 1.0_dp + aa/c; if (abs(c) < fpmin) c = fpmin
      d = 1.0_dp/d; h = h*d*c
      aa = -(a+real(m,dp))*(qab+real(m,dp))*x/((a+real(m2,dp))*(qap+real(m2,dp)))
      d = 1.0_dp + aa*d; if (abs(d) < fpmin) d = fpmin
      c = 1.0_dp + aa/c; if (abs(c) < fpmin) c = fpmin
      d = 1.0_dp/d; del = d*c; h = h*del
      if (abs(del-1.0_dp) < eps) exit
    end do
    cf = h
  end function beta_cf

  real(dp) function regularized_beta(a, b, x) result(v)
    real(dp), intent(in) :: a, b, x
    real(dp) :: bt
    if (x <= 0.0_dp) then
      v = 0.0_dp
      return
    else if (x >= 1.0_dp) then
      v = 1.0_dp
      return
    end if
    bt = exp(log_gamma(a+b)-log_gamma(a)-log_gamma(b)+a*log(x)+b*log(1.0_dp-x))
    if (x < (a+1.0_dp)/(a+b+2.0_dp)) then
      v = bt*beta_cf(a,b,x)/a
    else
      v = 1.0_dp-bt*beta_cf(b,a,1.0_dp-x)/b
    end if
    v = clamp(v, 0.0_dp, 1.0_dp)
  end function regularized_beta

  pure elemental real(dp) function student_pdf(x, df) result(f)
    real(dp), intent(in) :: x, df
    if (df <= 0.0_dp) then
      f = 0.0_dp
    else
      f = exp(log_gamma(0.5_dp*(df+1.0_dp))-log_gamma(0.5_dp*df)) / &
        sqrt(df*pi) * (1.0_dp+x*x/df)**(-0.5_dp*(df+1.0_dp))
    end if
  end function student_pdf

  real(dp) function student_cdf(x, df) result(p)
    real(dp), intent(in) :: x, df
    real(dp) :: z
    if (df <= 0.0_dp) then
      p = 0.0_dp
      return
    end if
    if (x == 0.0_dp) then
      p = 0.5_dp
      return
    end if
    z = df/(df+x*x)
    if (x > 0.0_dp) then
      p = 1.0_dp - 0.5_dp*regularized_beta(0.5_dp*df,0.5_dp,z)
    else
      p = 0.5_dp*regularized_beta(0.5_dp*df,0.5_dp,z)
    end if
  end function student_cdf

  real(dp) function student_quantile(p, df) result(x)
    real(dp), intent(in) :: p, df
    real(dp) :: lo, hi, mid
    integer :: i
    if (p <= 0.0_dp) then; x=-huge(1.0_dp); return; end if
    if (p >= 1.0_dp) then; x= huge(1.0_dp); return; end if
    lo = -1.0_dp; hi = 1.0_dp
    do while (student_cdf(lo,df) > p); lo=2.0_dp*lo; end do
    do while (student_cdf(hi,df) < p); hi=2.0_dp*hi; end do
    do i=1,120
      mid=0.5_dp*(lo+hi)
      if (student_cdf(mid,df) < p) then; lo=mid; else; hi=mid; end if
    end do
    x=0.5_dp*(lo+hi)
  end function student_quantile

  real(dp) function chi_square_cdf(x, df) result(p)
    real(dp), intent(in) :: x, df
    if (x <= 0.0_dp) then; p=0.0_dp; else; p=regularized_gamma_p(0.5_dp*df,0.5_dp*x); end if
  end function chi_square_cdf

  real(dp) function chi_square_quantile(p, df) result(x)
    real(dp), intent(in) :: p, df
    real(dp) :: lo, hi, mid
    integer :: i
    if (p <= 0.0_dp) then; x=0.0_dp; return; end if
    if (p >= 1.0_dp) then; x=huge(1.0_dp); return; end if
    lo=0.0_dp; hi=max(df,1.0_dp)
    do while (chi_square_cdf(hi,df)<p); hi=2.0_dp*hi; end do
    do i=1,120
      mid=0.5_dp*(lo+hi)
      if (chi_square_cdf(mid,df)<p) then; lo=mid; else; hi=mid; end if
    end do
    x=0.5_dp*(lo+hi)
  end function chi_square_quantile

  real(dp) function f_cdf(x, df1, df2) result(p)
    real(dp), intent(in) :: x, df1, df2
    if (x <= 0.0_dp) then
      p=0.0_dp
    else
      p=regularized_beta(0.5_dp*df1,0.5_dp*df2,df1*x/(df1*x+df2))
    end if
  end function f_cdf

  recursive real(dp) function adaptive_simpson(f, a, b, tol, max_depth) result(v)
    procedure(scalar_function) :: f
    real(dp), intent(in) :: a,b,tol
    integer, intent(in), optional :: max_depth
    integer :: md
    real(dp) :: fa, fb, fc, s
    md=20; if (present(max_depth)) md=max_depth
    fa=f(a); fb=f(b); fc=f(0.5_dp*(a+b)); s=(b-a)*(fa+4.0_dp*fc+fb)/6.0_dp
    v=asr(f,a,b,tol,s,fa,fb,fc,md)
  contains
    recursive real(dp) function asr(g,l,r,eps,whole,gl,gr,gc,depth) result(ans)
      procedure(scalar_function) :: g
      real(dp), intent(in) :: l,r,eps,whole,gl,gr,gc
      integer, intent(in) :: depth
      real(dp) :: m,lm,rm,glm,grm,left,right,delta
      m=0.5_dp*(l+r); lm=0.5_dp*(l+m); rm=0.5_dp*(m+r)
      glm=g(lm); grm=g(rm)
      left=(m-l)*(gl+4.0_dp*glm+gc)/6.0_dp
      right=(r-m)*(gc+4.0_dp*grm+gr)/6.0_dp
      delta=left+right-whole
      if (depth<=0 .or. abs(delta)<=15.0_dp*eps) then
        ans=left+right+delta/15.0_dp
      else
        ans=asr(g,l,m,0.5_dp*eps,left,gl,gc,glm,depth-1)+ &
            asr(g,m,r,0.5_dp*eps,right,gc,gr,grm,depth-1)
      end if
    end function asr
  end function adaptive_simpson

  real(dp) function bessel_k_nu(nu, x) result(k)
    real(dp), intent(in) :: nu, x
    real(dp) :: upper, old_value, new_value
    integer :: n
    if (x <= 0.0_dp) then
      k = huge(1.0_dp)
      return
    end if
    if (abs(nu) < 1.0e-13_dp) then
      k = bessel_k0_fast(x)
      return
    else if (abs(abs(nu)-1.0_dp) < 1.0e-13_dp) then
      k = bessel_k1_fast(x)
      return
    else if (abs(abs(nu)-0.5_dp-real(nint(abs(nu)-0.5_dp),dp)) < 1.0e-12_dp) then
      k = bessel_k_half_integer(abs(nu),x)
      return
    end if
    upper = max(8.0_dp, min(30.0_dp, log(1.0e14_dp/x + 1.0_dp) + 3.0_dp))
    n = 64
    old_value = bessel_simpson(nu, x, upper, n)
    do while (n < 32768)
      n = 2*n
      new_value = bessel_simpson(nu, x, upper, n)
      if (abs(new_value-old_value) <= 2.0e-11_dp*max(1.0_dp,abs(new_value))) exit
      old_value = new_value
    end do
    k = new_value
  end function bessel_k_nu

  real(dp) function bessel_i0_fast(x) result(v)
    real(dp),intent(in)::x
    real(dp)::ax,y
    ax=abs(x)
    if(ax<3.75_dp)then
      y=(x/3.75_dp)**2
      v=1.0_dp+y*(3.5156229_dp+y*(3.0899424_dp+y*(1.2067492_dp+y*(0.2659732_dp+y*(0.0360768_dp+y*0.0045813_dp)))))
    else
      y=3.75_dp/ax
      v=exp(ax)/sqrt(ax)*(0.39894228_dp+y*(0.01328592_dp+y*(0.00225319_dp+y*(-0.00157565_dp+y*(0.00916281_dp+y*(-0.02057706_dp+y*(0.02635537_dp+y*(-0.01647633_dp+y*0.00392377_dp))))))))
    end if
  end function bessel_i0_fast

  real(dp) function bessel_i1_fast(x) result(v)
    real(dp),intent(in)::x
    real(dp)::ax,y
    ax=abs(x)
    if(ax<3.75_dp)then
      y=(x/3.75_dp)**2
      v=ax*(0.5_dp+y*(0.87890594_dp+y*(0.51498869_dp+y*(0.15084934_dp+y*(0.02658733_dp+y*(0.00301532_dp+y*0.00032411_dp))))))
    else
      y=3.75_dp/ax
      v=exp(ax)/sqrt(ax)*(0.39894228_dp+y*(-0.03988024_dp+y*(-0.00362018_dp+y*(0.00163801_dp+y*(-0.01031555_dp+y*(0.02282967_dp+y*(-0.02895312_dp+y*(0.01787654_dp-y*0.00420059_dp))))))))
    end if
    if(x<0.0_dp)v=-v
  end function bessel_i1_fast

  real(dp) function bessel_k0_fast(x) result(v)
    real(dp),intent(in)::x
    real(dp)::y
    if(x<=2.0_dp)then
      y=x*x/4.0_dp
      v=-log(x/2.0_dp)*bessel_i0_fast(x)+(-0.57721566_dp+y*(0.42278420_dp+y*(0.23069756_dp+y*(0.03488590_dp+y*(0.00262698_dp+y*(0.00010750_dp+y*0.00000740_dp))))))
    else
      y=2.0_dp/x
      v=exp(-x)/sqrt(x)*(1.25331414_dp+y*(-0.07832358_dp+y*(0.02189568_dp+y*(-0.01062446_dp+y*(0.00587872_dp+y*(-0.00251540_dp+y*0.00053208_dp))))))
    end if
  end function bessel_k0_fast

  real(dp) function bessel_k1_fast(x) result(v)
    real(dp),intent(in)::x
    real(dp)::y
    if(x<=2.0_dp)then
      y=x*x/4.0_dp
      v=log(x/2.0_dp)*bessel_i1_fast(x)+(1.0_dp/x)*(1.0_dp+y*(0.15443144_dp+y*(-0.67278579_dp+y*(-0.18156897_dp+y*(-0.01919402_dp+y*(-0.00110404_dp+y*(-0.00004686_dp)))))))
    else
      y=2.0_dp/x
      v=exp(-x)/sqrt(x)*(1.25331414_dp+y*(0.23498619_dp+y*(-0.03655620_dp+y*(0.01504268_dp+y*(-0.00780353_dp+y*(0.00325614_dp+y*(-0.00068245_dp)))))))
    end if
  end function bessel_k1_fast

  real(dp) function bessel_k_half_integer(nu,x) result(v)
    real(dp),intent(in)::nu,x
    real(dp)::km1,kcur,knext,order
    integer::m,j
    m=nint(nu-0.5_dp)
    km1=sqrt(pi/(2.0_dp*x))*exp(-x)
    if(m==0)then
      v=km1
      return
    end if
    kcur=km1*(1.0_dp+1.0_dp/x)
    if(m==1)then
      v=kcur
      return
    end if
    order=1.5_dp
    do j=2,m
      knext=km1+2.0_dp*order/x*kcur
      km1=kcur
      kcur=knext
      order=order+1.0_dp
    end do
    v=kcur
  end function bessel_k_half_integer

  real(dp) function bessel_simpson(nu, x, upper, n) result(v)
    real(dp), intent(in) :: nu, x, upper
    integer, intent(in) :: n
    real(dp) :: h, t, y
    integer :: i
    h = upper/real(n,dp)
    v = bessel_integrand(0.0_dp,nu,x) + bessel_integrand(upper,nu,x)
    do i=1,n-1
      t = h*real(i,dp)
      y = bessel_integrand(t,nu,x)
      if (mod(i,2)==0) then
        v = v + 2.0_dp*y
      else
        v = v + 4.0_dp*y
      end if
    end do
    v = v*h/3.0_dp
  end function bessel_simpson

  real(dp) function bessel_integrand(t,nu,x) result(y)
    real(dp), intent(in) :: t,nu,x
    real(dp) :: exponent
    exponent = -x*cosh(t) + log(max(cosh(nu*t),tiny(1.0_dp)))
    if (exponent < log(tiny(1.0_dp))) then
      y = 0.0_dp
    else if (exponent > log(huge(1.0_dp))) then
      y = huge(1.0_dp)
    else
      y = exp(exponent)
    end if
  end function bessel_integrand
end module fbasics_special
