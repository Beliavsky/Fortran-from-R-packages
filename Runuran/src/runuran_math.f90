! SPDX-License-Identifier: GPL-2.0-or-later
! Derived from Runuran 0.41 / UNU.RAN by Wolfgang Hoermann and Josef Leydold.
module runuran_math
  use runuran_kinds, only : dp, pi, sqrt2, sqrt2pi, eps_dp, clamp
  implicit none
  private
  public :: normal_pdf, normal_cdf, normal_quantile, reg_gamma_p, reg_gamma_q
  public :: reg_beta, gamma_cdf, gamma_quantile, beta_cdf, beta_quantile
  public :: chisq_cdf, chisq_quantile, student_t_cdf, student_t_quantile
  public :: f_cdf, f_quantile, adaptive_integral, bisect_root
  public :: bessel_k_nu, log_bessel_k_nu, log_gamma_abs_complex
  public :: log1p_safe, expm1_safe

  abstract interface
    function scalar_fun(x) result(y)
      import dp
      real(dp), intent(in) :: x
      real(dp) :: y
    end function scalar_fun
  end interface
contains
  pure real(dp) function log1p_safe(x) result(y)
    real(dp), intent(in) :: x
    if (abs(x) < 1.0e-8_dp) then
      y=x*(1.0_dp-x*(0.5_dp-x*(1.0_dp/3.0_dp-x*0.25_dp)))
    else
      y=log(1.0_dp+x)
    end if
  end function log1p_safe

  pure real(dp) function expm1_safe(x) result(y)
    real(dp), intent(in) :: x
    if (abs(x) < 1.0e-7_dp) then
      y=x*(1.0_dp+x*(0.5_dp+x*(1.0_dp/6.0_dp+x/24.0_dp)))
    else
      y=exp(x)-1.0_dp
    end if
  end function expm1_safe

  pure real(dp) function normal_pdf(x) result(y)
    real(dp), intent(in) :: x
    y=exp(-0.5_dp*x*x)/sqrt2pi
  end function normal_pdf

  pure real(dp) function normal_cdf(x) result(p)
    real(dp), intent(in) :: x
    p=0.5_dp*erfc(-x/sqrt2)
  end function normal_cdf

  pure real(dp) function normal_quantile(p) result(x)
    real(dp), intent(in) :: p
    real(dp), parameter :: a(6)=[-3.969683028665376e1_dp,2.209460984245205e2_dp, &
      -2.759285104469687e2_dp,1.383577518672690e2_dp,-3.066479806614716e1_dp,2.506628277459239_dp]
    real(dp), parameter :: b(5)=[-5.447609879822406e1_dp,1.615858368580409e2_dp, &
      -1.556989798598866e2_dp,6.680131188771972e1_dp,-1.328068155288572e1_dp]
    real(dp), parameter :: c(6)=[-7.784894002430293e-3_dp,-3.223964580411365e-1_dp, &
      -2.400758277161838_dp,-2.549732539343734_dp,4.374664141464968_dp,2.938163982698783_dp]
    real(dp), parameter :: d(4)=[7.784695709041462e-3_dp,3.224671290700398e-1_dp, &
      2.445134137142996_dp,3.754408661907416_dp]
    real(dp) :: q,r,pp
    pp=clamp(p,tiny(1.0_dp),1.0_dp-epsilon(1.0_dp))
    if (pp < 0.02425_dp) then
      q=sqrt(-2.0_dp*log(pp))
      x=(((((c(1)*q+c(2))*q+c(3))*q+c(4))*q+c(5))*q+c(6))/ &
        ((((d(1)*q+d(2))*q+d(3))*q+d(4))*q+1.0_dp)
    else if (pp > 0.97575_dp) then
      q=sqrt(-2.0_dp*log(1.0_dp-pp))
      x=-(((((c(1)*q+c(2))*q+c(3))*q+c(4))*q+c(5))*q+c(6))/ &
        ((((d(1)*q+d(2))*q+d(3))*q+d(4))*q+1.0_dp)
    else
      q=pp-0.5_dp; r=q*q
      x=(((((a(1)*r+a(2))*r+a(3))*r+a(4))*r+a(5))*r+a(6))*q/ &
        (((((b(1)*r+b(2))*r+b(3))*r+b(4))*r+b(5))*r+1.0_dp)
    end if
    ! Halley refinement.
    x=x-(normal_cdf(x)-pp)/(normal_pdf(x)+0.5_dp*x*(normal_cdf(x)-pp))
  end function normal_quantile

  pure real(dp) function reg_gamma_p(a,x) result(p)
    real(dp), intent(in) :: a,x
    integer :: n
    real(dp) :: sum,del,ap,b,c,d,h,an
    if (a <= 0.0_dp .or. x < 0.0_dp) then; p=0.0_dp; return; end if
    if (x <= 0.0_dp) then; p=0.0_dp; return; end if
    if (x < a+1.0_dp) then
      ap=a; del=1.0_dp/a; sum=del
      do n=1,10000
        ap=ap+1.0_dp; del=del*x/ap; sum=sum+del
        if (abs(del) <= abs(sum)*1.0e-15_dp) exit
      end do
      p=sum*exp(-x+a*log(x)-log_gamma(a))
    else
      b=x+1.0_dp-a; c=1.0_dp/tiny(1.0_dp); d=1.0_dp/b; h=d
      do n=1,10000
        an=-real(n,dp)*(real(n,dp)-a)
        b=b+2.0_dp; d=an*d+b
        if (abs(d)<tiny(1.0_dp)) d=tiny(1.0_dp)
        c=b+an/c; if(abs(c)<tiny(1.0_dp)) c=tiny(1.0_dp)
        d=1.0_dp/d; del=d*c; h=h*del
        if (abs(del-1.0_dp) <= 1.0e-15_dp) exit
      end do
      p=1.0_dp-exp(-x+a*log(x)-log_gamma(a))*h
    end if
    p=clamp(p,0.0_dp,1.0_dp)
  end function reg_gamma_p

  pure real(dp) function reg_gamma_q(a,x) result(q)
    real(dp), intent(in) :: a,x
    q=1.0_dp-reg_gamma_p(a,x)
  end function reg_gamma_q

  pure real(dp) function betacf(a,b,x) result(h)
    real(dp), intent(in) :: a,b,x
    integer :: m,m2
    real(dp) :: aa,c,d,del,qab,qam,qap
    qab=a+b; qap=a+1.0_dp; qam=a-1.0_dp
    c=1.0_dp; d=1.0_dp-qab*x/qap
    if(abs(d)<tiny(1.0_dp)) d=tiny(1.0_dp)
    d=1.0_dp/d; h=d
    do m=1,10000
      m2=2*m
      aa=real(m,dp)*(b-real(m,dp))*x/((qam+real(m2,dp))*(a+real(m2,dp)))
      d=1.0_dp+aa*d; if(abs(d)<tiny(1.0_dp))d=tiny(1.0_dp)
      c=1.0_dp+aa/c; if(abs(c)<tiny(1.0_dp))c=tiny(1.0_dp)
      d=1.0_dp/d; h=h*d*c
      aa=-(a+real(m,dp))*(qab+real(m,dp))*x/((a+real(m2,dp))*(qap+real(m2,dp)))
      d=1.0_dp+aa*d; if(abs(d)<tiny(1.0_dp))d=tiny(1.0_dp)
      c=1.0_dp+aa/c; if(abs(c)<tiny(1.0_dp))c=tiny(1.0_dp)
      d=1.0_dp/d; del=d*c; h=h*del
      if(abs(del-1.0_dp)<1.0e-15_dp) exit
    end do
  end function betacf

  pure real(dp) function reg_beta(x,a,b) result(p)
    real(dp), intent(in) :: x,a,b
    real(dp) :: bt
    if (x <= 0.0_dp) then; p=0.0_dp; return; end if
    if (x >= 1.0_dp) then; p=1.0_dp; return; end if
    bt=exp(log_gamma(a+b)-log_gamma(a)-log_gamma(b)+a*log(x)+b*log1p_safe(-x))
    if (x < (a+1.0_dp)/(a+b+2.0_dp)) then
      p=bt*betacf(a,b,x)/a
    else
      p=1.0_dp-bt*betacf(b,a,1.0_dp-x)/b
    end if
    p=clamp(p,0.0_dp,1.0_dp)
  end function reg_beta

  pure real(dp) function gamma_cdf(x,shape,scale) result(p)
    real(dp), intent(in) :: x,shape,scale
    if (x <= 0.0_dp) then; p=0.0_dp; else; p=reg_gamma_p(shape,x/scale); end if
  end function gamma_cdf

  real(dp) function gamma_quantile(p,shape,scale) result(x)
    real(dp), intent(in) :: p,shape,scale
    real(dp) :: lo,hi,mid
    integer :: i
    if(p<=0.0_dp)then;x=0.0_dp;return;end if
    if(p>=1.0_dp)then;x=huge(1.0_dp);return;end if
    lo=0.0_dp; hi=max(scale,shape*scale)
    do while(gamma_cdf(hi,shape,scale)<p .and. hi<huge(1.0_dp)/4.0_dp); hi=2.0_dp*hi; end do
    do i=1,150
      mid=0.5_dp*(lo+hi)
      if(gamma_cdf(mid,shape,scale)<p)then;lo=mid;else;hi=mid;end if
    end do
    x=0.5_dp*(lo+hi)
  end function gamma_quantile

  pure real(dp) function beta_cdf(x,a,b) result(p)
    real(dp), intent(in) :: x,a,b
    p=reg_beta(x,a,b)
  end function beta_cdf

  real(dp) function beta_quantile(p,a,b) result(x)
    real(dp), intent(in)::p,a,b
    real(dp)::lo,hi,mid
    integer::i
    lo=0.0_dp;hi=1.0_dp
    do i=1,120
      mid=0.5_dp*(lo+hi)
      if(reg_beta(mid,a,b)<p)then;lo=mid;else;hi=mid;end if
    end do
    x=0.5_dp*(lo+hi)
  end function beta_quantile

  pure real(dp) function chisq_cdf(x,df) result(p)
    real(dp),intent(in)::x,df
    p=gamma_cdf(x,0.5_dp*df,2.0_dp)
  end function chisq_cdf
  real(dp) function chisq_quantile(p,df) result(x)
    real(dp),intent(in)::p,df
    x=gamma_quantile(p,0.5_dp*df,2.0_dp)
  end function chisq_quantile

  pure real(dp) function student_t_cdf(x,df) result(p)
    real(dp),intent(in)::x,df
    real(dp)::z
    if(abs(x)<=tiny(1.0_dp))then;p=0.5_dp;return;end if
    z=df/(df+x*x)
    if(x>0.0_dp)then;p=1.0_dp-0.5_dp*reg_beta(z,0.5_dp*df,0.5_dp)
    else;p=0.5_dp*reg_beta(z,0.5_dp*df,0.5_dp);end if
  end function student_t_cdf
  real(dp) function student_t_quantile(p,df) result(x)
    real(dp),intent(in)::p,df
    real(dp)::lo,hi,mid
    integer::i
    lo=-1.0_dp;hi=1.0_dp
    do while(student_t_cdf(lo,df)>p);lo=2.0_dp*lo;end do
    do while(student_t_cdf(hi,df)<p);hi=2.0_dp*hi;end do
    do i=1,140
      mid=0.5_dp*(lo+hi)
      if(student_t_cdf(mid,df)<p)then;lo=mid;else;hi=mid;end if
    end do
    x=0.5_dp*(lo+hi)
  end function student_t_quantile

  pure real(dp) function f_cdf(x,df1,df2) result(p)
    real(dp),intent(in)::x,df1,df2
    if(x<=0.0_dp)then;p=0.0_dp;else;p=reg_beta(df1*x/(df1*x+df2),0.5_dp*df1,0.5_dp*df2);end if
  end function f_cdf
  real(dp) function f_quantile(p,df1,df2) result(x)
    real(dp),intent(in)::p,df1,df2
    real(dp)::z
    z=beta_quantile(p,0.5_dp*df1,0.5_dp*df2)
    x=df2*z/(df1*(1.0_dp-z))
  end function f_quantile

  recursive real(dp) function simpson_rec(f,a,b,fa,fb,fm,s,eps,depth) result(v)
    procedure(scalar_fun) :: f
    real(dp),intent(in)::a,b,fa,fb,fm,s,eps
    integer,intent(in)::depth
    real(dp)::m,lm,rm,fl,fr,sl,sr
    m=0.5_dp*(a+b);lm=0.5_dp*(a+m);rm=0.5_dp*(m+b)
    fl=f(lm);fr=f(rm);sl=(m-a)*(fa+4.0_dp*fl+fm)/6.0_dp;sr=(b-m)*(fm+4.0_dp*fr+fb)/6.0_dp
    if(depth<=0 .or. abs(sl+sr-s)<=15.0_dp*eps) then
      v=sl+sr+(sl+sr-s)/15.0_dp
    else
      v=simpson_rec(f,a,m,fa,fm,fl,sl,0.5_dp*eps,depth-1)+ &
        simpson_rec(f,m,b,fm,fb,fr,sr,0.5_dp*eps,depth-1)
    end if
  end function simpson_rec

  real(dp) function adaptive_integral(f,a,b,tol) result(v)
    procedure(scalar_fun) :: f
    real(dp),intent(in)::a,b
    real(dp),intent(in),optional::tol
    real(dp)::fa,fb,fm,s,e
    e=1.0e-10_dp;if(present(tol))e=tol
    fa=f(a);fb=f(b);fm=f(0.5_dp*(a+b));s=(b-a)*(fa+4.0_dp*fm+fb)/6.0_dp
    v=simpson_rec(f,a,b,fa,fb,fm,s,e,25)
  end function adaptive_integral

  real(dp) function bisect_root(f,a,b,tol) result(x)
    procedure(scalar_fun)::f
    real(dp),intent(in)::a,b
    real(dp),intent(in),optional::tol
    real(dp)::lo,hi,mid,fl,fm,e
    integer::i
    e=1.0e-12_dp;if(present(tol))e=tol
    lo=a;hi=b;fl=f(lo)
    do i=1,200
      mid=0.5_dp*(lo+hi);fm=f(mid)
      if(abs(fm)<e .or. abs(hi-lo)<e*(1.0_dp+abs(mid)))exit
      if((fm>=0.0_dp).eqv.(fl>=0.0_dp))then;lo=mid;fl=fm;else;hi=mid;end if
    end do
    x=0.5_dp*(lo+hi)
  end function bisect_root

  pure real(dp) function bessel_integrand(t,x,nu) result(y)
    real(dp),intent(in)::t,x,nu
    real(dp)::z
    z=-x*cosh(t)+abs(nu)*t
    if(z < -740.0_dp)then
      y=0.0_dp
    else
      y=exp(-x*cosh(t))*cosh(nu*t)
    end if
  end function bessel_integrand

  recursive real(dp) function bessel_simpson(x,nu,a,b,fa,fb,fm,s,eps,depth) result(v)
    real(dp),intent(in)::x,nu,a,b,fa,fb,fm,s,eps
    integer,intent(in)::depth
    real(dp)::m,lm,rm,fl,fr,sl,sr
    m=0.5_dp*(a+b)
    lm=0.5_dp*(a+m)
    rm=0.5_dp*(m+b)
    fl=bessel_integrand(lm,x,nu)
    fr=bessel_integrand(rm,x,nu)
    sl=(m-a)*(fa+4.0_dp*fl+fm)/6.0_dp
    sr=(b-m)*(fm+4.0_dp*fr+fb)/6.0_dp
    if(depth<=0 .or. abs(sl+sr-s)<=15.0_dp*eps)then
      v=sl+sr+(sl+sr-s)/15.0_dp
    else
      v=bessel_simpson(x,nu,a,m,fa,fm,fl,sl,0.5_dp*eps,depth-1)+ &
        bessel_simpson(x,nu,m,b,fm,fb,fr,sr,0.5_dp*eps,depth-1)
    end if
  end function bessel_simpson

  real(dp) function bessel_k_nu(x,nu) result(k)
    real(dp),intent(in)::x,nu
    real(dp)::tmax,a,b,fa,fb,fm,s
    if(x<=0.0_dp)then
      k=huge(1.0_dp)
      return
    end if
    tmax=max(12.0_dp,log(max(4.0_dp,80.0_dp/x)))
    a=0.0_dp
    b=tmax
    fa=bessel_integrand(a,x,nu)
    fb=bessel_integrand(b,x,nu)
    fm=bessel_integrand(0.5_dp*(a+b),x,nu)
    s=(b-a)*(fa+4.0_dp*fm+fb)/6.0_dp
    k=bessel_simpson(x,nu,a,b,fa,fb,fm,s,1.0e-10_dp,25)
  end function bessel_k_nu

  real(dp) function log_bessel_k_nu(x,nu) result(lk)
    real(dp),intent(in)::x,nu
    real(dp)::k
    if(x>60.0_dp)then
      lk=0.5_dp*log(pi/(2.0_dp*x))-x+log(1.0_dp+(4.0_dp*nu*nu-1.0_dp)/(8.0_dp*x))
    else
      k=bessel_k_nu(x,nu)
      lk=log(max(k,tiny(1.0_dp)))
    end if
  end function log_bessel_k_nu

  recursive real(dp) function log_gamma_abs_complex(a,y) result(v)
    real(dp),intent(in)::a,y
    complex(dp), parameter :: one=(1.0_dp,0.0_dp)
    real(dp), parameter :: g=7.0_dp
    real(dp), parameter :: c(9)=[0.99999999999980993_dp,676.5203681218851_dp,-1259.1392167224028_dp, &
       771.32342877765313_dp,-176.61502916214059_dp,12.507343278686905_dp,-0.13857109526572012_dp, &
       9.9843695780195716e-6_dp,1.5056327351493116e-7_dp]
    complex(dp)::z,w,t,s
    integer::j
    z=cmplx(a,y,dp)
    if(real(z,dp)<0.5_dp)then
      v=log(pi)-log(abs(sin(pi*z)))-log_gamma_abs_complex(1.0_dp-a,-y)
      return
    end if
    w=z-one;s=cmplx(c(1),0.0_dp,dp)
    do j=2,9;s=s+c(j)/(w+real(j-1,dp));end do
    t=w+g+0.5_dp
    v=real(0.5_dp*log(2.0_dp*pi)+(w+0.5_dp)*log(t)-t+log(s),dp)
  end function log_gamma_abs_complex
end module runuran_math
