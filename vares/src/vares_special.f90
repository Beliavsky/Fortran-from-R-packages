! SPDX-License-Identifier: GPL-2.0-or-later
! Modern Fortran translation of VaRES 1.0.2.
! See NOTICE.md and original/ for attribution and provenance.
module vares_special
  use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
  use vares_kinds, only : dp, pi
  implicit none
  private
  public :: gamma_fn, log_gamma_fn, beta_fn, log_beta_fn, r_sign
  public :: normal_pdf, normal_cdf, normal_quantile
  public :: beta_pdf, beta_cdf, beta_quantile
  public :: gamma_pdf, gamma_cdf, gamma_quantile
  public :: student_t_pdf, student_t_cdf, student_t_quantile
  public :: f_pdf, f_cdf, f_quantile
  public :: lognormal_pdf, lognormal_cdf, lognormal_quantile
  public :: logistic_pdf, logistic_cdf, logistic_quantile
  public :: cauchy_pdf, cauchy_cdf, cauchy_quantile
  public :: uniform_pdf, uniform_cdf, uniform_quantile
  public :: weibull_pdf, weibull_cdf, weibull_quantile
contains
  pure elemental function gamma_fn(x) result(y)
    real(dp), intent(in) :: x
    real(dp) :: y
    y = gamma(x)
  end function gamma_fn

  pure elemental function log_gamma_fn(x) result(y)
    real(dp), intent(in) :: x
    real(dp) :: y
    y = log_gamma(x)
  end function log_gamma_fn

  pure elemental function beta_fn(a, b) result(y)
    real(dp), intent(in) :: a, b
    real(dp) :: y
    y = exp(log_gamma(a) + log_gamma(b) - log_gamma(a + b))
  end function beta_fn

  pure elemental function log_beta_fn(a, b) result(y)
    real(dp), intent(in) :: a, b
    real(dp) :: y
    y = log_gamma(a) + log_gamma(b) - log_gamma(a + b)
  end function log_beta_fn

  pure elemental function r_sign(x) result(y)
    real(dp), intent(in) :: x
    real(dp) :: y
    if (x > 0.0_dp) then
      y = 1.0_dp
    else if (x < 0.0_dp) then
      y = -1.0_dp
    else
      y = 0.0_dp
    end if
  end function r_sign

  pure elemental function normal_pdf(x, mean, sd, log_pdf) result(y)
    real(dp), intent(in) :: x
    real(dp), intent(in), optional :: mean, sd
    logical, intent(in), optional :: log_pdf
    real(dp) :: y, m, s, z
    logical :: lp
    m=0.0_dp
    if (present(mean)) m=mean
    s=1.0_dp
    if (present(sd)) s=sd
    lp=.false.
    if (present(log_pdf)) lp=log_pdf
    z=(x-m)/s
    y=-0.5_dp*z*z-log(s)-0.5_dp*log(2.0_dp*pi)
    if (.not. lp) y=exp(y)
  end function normal_pdf

  pure elemental function normal_cdf(x, mean, sd, log_p, lower_tail) result(y)
    real(dp), intent(in) :: x
    real(dp), intent(in), optional :: mean, sd
    logical, intent(in), optional :: log_p, lower_tail
    real(dp) :: y, m, s
    logical :: lp, lt
    m=0.0_dp
    if (present(mean)) m=mean
    s=1.0_dp
    if (present(sd)) s=sd
    lp=.false.
    if (present(log_p)) lp=log_p
    lt=.true.
    if (present(lower_tail)) lt=lower_tail
    if (lt) then
      y=0.5_dp*erfc(-(x-m)/(s*sqrt(2.0_dp)))
    else
      y=0.5_dp*erfc((x-m)/(s*sqrt(2.0_dp)))
    end if
    if (lp) y=log(y)
  end function normal_cdf

  pure elemental function normal_quantile(p, mean, sd, log_p, lower_tail) result(x)
    real(dp), intent(in) :: p
    real(dp), intent(in), optional :: mean, sd
    logical, intent(in), optional :: log_p, lower_tail
    real(dp) :: x, pp, m, s, q, r
    logical :: lp, lt
    real(dp), parameter :: a1=-3.969683028665376e1_dp,a2=2.209460984245205e2_dp,&
      a3=-2.759285104469687e2_dp,a4=1.383577518672690e2_dp,a5=-3.066479806614716e1_dp,a6=2.506628277459239_dp
    real(dp), parameter :: b1=-5.447609879822406e1_dp,b2=1.615858368580409e2_dp,&
      b3=-1.556989798598866e2_dp,b4=6.680131188771972e1_dp,b5=-1.328068155288572e1_dp
    real(dp), parameter :: c1=-7.784894002430293e-3_dp,c2=-3.223964580411365e-1_dp,&
      c3=-2.400758277161838_dp,c4=-2.549732539343734_dp,c5=4.374664141464968_dp,c6=2.938163982698783_dp
    real(dp), parameter :: d1=7.784695709041462e-3_dp,d2=3.224671290700398e-1_dp,&
      d3=2.445134137142996_dp,d4=3.754408661907416_dp
    real(dp), parameter :: plow=0.02425_dp, phigh=1.0_dp-plow
    m=0.0_dp
    if (present(mean)) m=mean
    s=1.0_dp
    if (present(sd)) s=sd
    lp=.false.
    if (present(log_p)) lp=log_p
    lt=.true.
    if (present(lower_tail)) lt=lower_tail
    pp=p
    if (lp) pp=exp(pp)
    if (.not.lt) pp=1.0_dp-pp
    if (pp<=0.0_dp) then
    x=-huge(1.0_dp)
    return
    else if (pp>=1.0_dp) then
    x=huge(1.0_dp)
    return
    else if (pp<plow) then
      q=sqrt(-2.0_dp*log(pp))
      x=(((((c1*q+c2)*q+c3)*q+c4)*q+c5)*q+c6)/((((d1*q+d2)*q+d3)*q+d4)*q+1.0_dp)
    else if (pp<=phigh) then
      q=pp-0.5_dp
      r=q*q
      x=((((((a1*r+a2)*r+a3)*r+a4)*r+a5)*r+a6)*q)/(((((b1*r+b2)*r+b3)*r+b4)*r+b5)*r+1.0_dp)
    else
      q=sqrt(-2.0_dp*log(1.0_dp-pp))
      x=-(((((c1*q+c2)*q+c3)*q+c4)*q+c5)*q+c6)/((((d1*q+d2)*q+d3)*q+d4)*q+1.0_dp)
    end if
    ! One Halley refinement.
    x=x-(normal_cdf(x)-pp)/normal_pdf(x)/(1.0_dp+0.5_dp*x*(normal_cdf(x)-pp)/normal_pdf(x))
    x=m+s*x
  end function normal_quantile

  pure function beta_cf(a,b,x) result(cf)
    real(dp), intent(in) :: a,b,x
    real(dp) :: cf, qab,qap,qam,c,d,h,aa,del
    integer :: m,m2
    real(dp), parameter :: fpmin=tiny(1.0_dp)/epsilon(1.0_dp)
    qab=a+b
    qap=a+1.0_dp
    qam=a-1.0_dp
    c=1.0_dp
    d=1.0_dp-qab*x/qap
    if(abs(d)<fpmin)d=fpmin
    d=1.0_dp/d
    h=d
    do m=1,300
      m2=2*m
      aa=real(m,dp)*(b-real(m,dp))*x/((qam+real(m2,dp))*(a+real(m2,dp)))
      d=1.0_dp+aa*d
      if(abs(d)<fpmin)d=fpmin
      c=1.0_dp+aa/c
      if(abs(c)<fpmin)c=fpmin
      d=1.0_dp/d
      h=h*d*c
      aa=-(a+real(m,dp))*(qab+real(m,dp))*x/((a+real(m2,dp))*(qap+real(m2,dp)))
      d=1.0_dp+aa*d
      if(abs(d)<fpmin)d=fpmin
      c=1.0_dp+aa/c
      if(abs(c)<fpmin)c=fpmin
      d=1.0_dp/d
      del=d*c
      h=h*del
      if(abs(del-1.0_dp)<8.0_dp*epsilon(1.0_dp)) exit
    end do
    cf=h
  end function beta_cf

  pure elemental function reg_beta(x,a,b) result(y)
    real(dp), intent(in) :: x,a,b
    real(dp) :: y,bt
    if(x<=0.0_dp) then
    y=0.0_dp
    return
    else if(x>=1.0_dp) then
    y=1.0_dp
    return
    end if
    bt=exp(log_gamma(a+b)-log_gamma(a)-log_gamma(b)+a*log(x)+b*log(1.0_dp-x))
    if(x<(a+1.0_dp)/(a+b+2.0_dp)) then
      y=bt*beta_cf(a,b,x)/a
    else
      y=1.0_dp-bt*beta_cf(b,a,1.0_dp-x)/b
    end if
    y=max(0.0_dp,min(1.0_dp,y))
  end function reg_beta

  pure elemental function beta_pdf(x, shape1, shape2, log_pdf) result(y)
    real(dp), intent(in) :: x
    real(dp), intent(in) :: shape1,shape2
    logical, intent(in), optional :: log_pdf
    real(dp) :: y
    logical :: lp
    lp=.false.
    if(present(log_pdf))lp=log_pdf
    if(x<0.0_dp .or. x>1.0_dp) then
      y=merge(-huge(1.0_dp),0.0_dp,lp)
      return
    end if
    y=(shape1-1.0_dp)*log(x)+(shape2-1.0_dp)*log(1.0_dp-x)-log_beta_fn(shape1,shape2)
    if(.not.lp)y=exp(y)
  end function beta_pdf

  pure elemental function beta_cdf(x, shape1, shape2, log_p, lower_tail) result(y)
    real(dp), intent(in) :: x,shape1,shape2
    logical, intent(in), optional :: log_p,lower_tail
    real(dp) :: y
    logical :: lp,lt
    lp=.false.
    if(present(log_p))lp=log_p
    lt=.true.
    if(present(lower_tail))lt=lower_tail
    y=reg_beta(x,shape1,shape2)
    if(.not.lt)y=1.0_dp-y
    if(lp)y=log(y)
  end function beta_cdf

  pure elemental function beta_quantile(p, shape1, shape2, log_p, lower_tail) result(x)
    real(dp), intent(in)::p,shape1,shape2
    logical,intent(in),optional::log_p,lower_tail
    real(dp)::x,pp,lo,hi,mid
    logical::lp,lt
    integer::i
    lp=.false.
    if(present(log_p))lp=log_p
    lt=.true.
    if(present(lower_tail))lt=lower_tail
    pp=p
    if(lp)pp=exp(pp)
    if(.not.lt)pp=1.0_dp-pp
    if(pp<=0.0_dp)then
    x=0.0_dp
    return
    else if(pp>=1.0_dp)then
    x=1.0_dp
    return
    end if
    lo=0.0_dp
    hi=1.0_dp
    do i=1,100
      mid=0.5_dp*(lo+hi)
      if(reg_beta(mid,shape1,shape2)<pp)then
      lo=mid
      else
      hi=mid
      end if
    end do
    x=0.5_dp*(lo+hi)
  end function beta_quantile

  pure elemental function reg_gamma_p(a,x) result(y)
    real(dp),intent(in)::a,x
    real(dp)::y,ap,del,sumv,b,c,d,h,an
    integer::n
    real(dp),parameter::fpmin=tiny(1.0_dp)/epsilon(1.0_dp)
    if(x<=0.0_dp)then
    y=0.0_dp
    return
    end if
    if(x<a+1.0_dp)then
      ap=a
      sumv=1.0_dp/a
      del=sumv
      do n=1,1000
        ap=ap+1.0_dp
        del=del*x/ap
        sumv=sumv+del
        if(abs(del)<abs(sumv)*8.0_dp*epsilon(1.0_dp))exit
      end do
      y=sumv*exp(-x+a*log(x)-log_gamma(a))
    else
      b=x+1.0_dp-a
      c=1.0_dp/fpmin
      d=1.0_dp/b
      h=d
      do n=1,1000
        an=-real(n,dp)*(real(n,dp)-a)
        b=b+2.0_dp
        d=an*d+b
        if(abs(d)<fpmin)d=fpmin
        c=b+an/c
        if(abs(c)<fpmin)c=fpmin
        d=1.0_dp/d
        del=d*c
        h=h*del
        if(abs(del-1.0_dp)<8.0_dp*epsilon(1.0_dp))exit
      end do
      y=1.0_dp-exp(-x+a*log(x)-log_gamma(a))*h
    end if
    y=max(0.0_dp,min(1.0_dp,y))
  end function reg_gamma_p

  pure elemental function gamma_pdf(x, shape, rate, scale, log_pdf) result(y)
    real(dp),intent(in)::x,shape
    real(dp),intent(in),optional::rate,scale
    logical,intent(in),optional::log_pdf
    real(dp)::y,sc
    logical::lp
    sc=1.0_dp
    if(present(scale))sc=scale
    if(present(rate))sc=1.0_dp/rate
    lp=.false.
    if(present(log_pdf))lp=log_pdf
    if(x<0.0_dp)then
    y=merge(-huge(1.0_dp),0.0_dp,lp)
    return
    end if
    y=(shape-1.0_dp)*log(x)-x/sc-log_gamma(shape)-shape*log(sc)
    if(.not.lp)y=exp(y)
  end function gamma_pdf

  pure elemental function gamma_cdf(x, shape, rate, scale, log_p, lower_tail) result(y)
    real(dp),intent(in)::x,shape
    real(dp),intent(in),optional::rate,scale
    logical,intent(in),optional::log_p,lower_tail
    real(dp)::y,sc
    logical::lp,lt
    sc=1.0_dp
    if(present(scale))sc=scale
    if(present(rate))sc=1.0_dp/rate
    lp=.false.
    if(present(log_p))lp=log_p
    lt=.true.
    if(present(lower_tail))lt=lower_tail
    y=reg_gamma_p(shape,x/sc)
    if(.not.lt)y=1.0_dp-y
    if(lp)y=log(y)
  end function gamma_cdf

  pure elemental function gamma_quantile(p, shape, rate, scale, log_p, lower_tail) result(x)
    real(dp),intent(in)::p,shape
    real(dp),intent(in),optional::rate,scale
    logical,intent(in),optional::log_p,lower_tail
    real(dp)::x,pp,sc,lo,hi,mid
    logical::lp,lt
    integer::i
    sc=1.0_dp
    if(present(scale))sc=scale
    if(present(rate))sc=1.0_dp/rate
    lp=.false.
    if(present(log_p))lp=log_p
    lt=.true.
    if(present(lower_tail))lt=lower_tail
    pp=p
    if(lp)pp=exp(pp)
    if(.not.lt)pp=1.0_dp-pp
    if(pp<=0.0_dp)then
    x=0.0_dp
    return
    else if(pp>=1.0_dp)then
    x=huge(1.0_dp)
    return
    end if
    lo=0.0_dp
    hi=max(1.0_dp,shape)
    do while(reg_gamma_p(shape,hi)<pp)
    hi=2.0_dp*hi
    if(hi>huge(1.0_dp)/4.0_dp)exit
    end do
    do i=1,120
      mid=0.5_dp*(lo+hi)
      if(reg_gamma_p(shape,mid)<pp)then
      lo=mid
      else
      hi=mid
      end if
    end do
    x=0.5_dp*(lo+hi)*sc
  end function gamma_quantile

  pure elemental function student_t_pdf(x, df, log_pdf) result(y)
    real(dp),intent(in)::x,df
    logical,intent(in),optional::log_pdf
    real(dp)::y
    logical::lp
    lp=.false.
    if(present(log_pdf))lp=log_pdf
    y=log_gamma(0.5_dp*(df+1.0_dp))-log_gamma(0.5_dp*df)-0.5_dp*log(df*pi)-0.5_dp*(df+1.0_dp)*log(1.0_dp+x*x/df)
    if(.not.lp)y=exp(y)
  end function student_t_pdf

  pure elemental function student_t_cdf(x, df, log_p, lower_tail) result(y)
    real(dp),intent(in)::x,df
    logical,intent(in),optional::log_p,lower_tail
    real(dp)::y,z
    logical::lp,lt
    lp=.false.
    if(present(log_p))lp=log_p
    lt=.true.
    if(present(lower_tail))lt=lower_tail
    z=df/(df+x*x)
    if(x>=0.0_dp)then
    y=1.0_dp-0.5_dp*reg_beta(z,0.5_dp*df,0.5_dp)
    else
    y=0.5_dp*reg_beta(z,0.5_dp*df,0.5_dp)
    end if
    if(.not.lt)y=1.0_dp-y
    if(lp)y=log(y)
  end function student_t_cdf

  pure elemental function student_t_quantile(p, df, log_p, lower_tail) result(x)
    real(dp),intent(in)::p,df
    logical,intent(in),optional::log_p,lower_tail
    real(dp)::x,pp,z
    logical::lp,lt
    lp=.false.
    if(present(log_p))lp=log_p
    lt=.true.
    if(present(lower_tail))lt=lower_tail
    pp=p
    if(lp)pp=exp(pp)
    if(.not.lt)pp=1.0_dp-pp
    if(pp<=0.0_dp)then
    x=-huge(1.0_dp)
    return
    else if(pp>=1.0_dp)then
    x=huge(1.0_dp)
    return
    end if
    if (abs(pp - 0.5_dp) <= epsilon(1.0_dp)) then
      x = 0.0_dp
    else if (pp < 0.5_dp) then
      z = beta_quantile(2.0_dp*pp, 0.5_dp*df, 0.5_dp)
      x = -sqrt(df*(1.0_dp-z)/z)
    else
      z = beta_quantile(2.0_dp*(1.0_dp-pp), 0.5_dp*df, 0.5_dp)
      x = sqrt(df*(1.0_dp-z)/z)
    end if
  end function student_t_quantile

  pure elemental function f_pdf(x, df1, df2, log_pdf) result(y)
    real(dp),intent(in)::x
    real(dp),intent(in),optional::df1,df2
    logical,intent(in),optional::log_pdf
    real(dp)::y,d1,d2
    logical::lp
    d1=1.0_dp
    if(present(df1))d1=df1
    d2=1.0_dp
    if(present(df2))d2=df2
    lp=.false.
    if(present(log_pdf))lp=log_pdf
    y=0.5_dp*d1*log(d1/d2)+(0.5_dp*d1-1.0_dp)*log(x)-0.5_dp*(d1+d2)*log(1.0_dp+d1*x/d2)-log_beta_fn(0.5_dp*d1,0.5_dp*d2)
    if(.not.lp)y=exp(y)
  end function f_pdf

  pure elemental function f_cdf(x, df1, df2, log_p, lower_tail) result(y)
    real(dp),intent(in)::x
    real(dp),intent(in),optional::df1,df2
    logical,intent(in),optional::log_p,lower_tail
    real(dp)::y,d1,d2
    logical::lp,lt
    d1=1.0_dp
    if(present(df1))d1=df1
    d2=1.0_dp
    if(present(df2))d2=df2
    lp=.false.
    if(present(log_p))lp=log_p
    lt=.true.
    if(present(lower_tail))lt=lower_tail
    y=reg_beta(d1*x/(d1*x+d2),0.5_dp*d1,0.5_dp*d2)
    if(.not.lt)y=1.0_dp-y
    if(lp)y=log(y)
  end function f_cdf

  pure elemental function f_quantile(p, df1, df2, log_p, lower_tail) result(x)
    real(dp),intent(in)::p
    real(dp),intent(in),optional::df1,df2
    logical,intent(in),optional::log_p,lower_tail
    real(dp)::x,pp,d1,d2,z
    logical::lp,lt
    d1=1.0_dp
    if(present(df1))d1=df1
    d2=1.0_dp
    if(present(df2))d2=df2
    lp=.false.
    if(present(log_p))lp=log_p
    lt=.true.
    if(present(lower_tail))lt=lower_tail
    pp=p
    if(lp)pp=exp(pp)
    if(.not.lt)pp=1.0_dp-pp
    z=beta_quantile(pp,0.5_dp*d1,0.5_dp*d2)
    x=d2*z/(d1*(1.0_dp-z))
  end function f_quantile

  pure elemental function lognormal_pdf(x, meanlog, sdlog, log_pdf) result(y)
    real(dp),intent(in)::x
    real(dp),intent(in),optional::meanlog,sdlog
    logical,intent(in),optional::log_pdf
    real(dp)::y,m,s
    logical::lp
    m=0.0_dp
    if(present(meanlog))m=meanlog
    s=1.0_dp
    if(present(sdlog))s=sdlog
    lp=.false.
    if(present(log_pdf))lp=log_pdf
    y=normal_pdf(log(x),m,s,.true.)-log(x)
    if(.not.lp)y=exp(y)
  end function lognormal_pdf
  pure elemental function lognormal_cdf(x, meanlog, sdlog, log_p, lower_tail) result(y)
    real(dp),intent(in)::x
    real(dp),intent(in),optional::meanlog,sdlog
    logical,intent(in),optional::log_p,lower_tail
    real(dp)::y,m,s
    m=0.0_dp
    if(present(meanlog))m=meanlog
    s=1.0_dp
    if(present(sdlog))s=sdlog
    y=normal_cdf(log(x),m,s,log_p,lower_tail)
  end function lognormal_cdf
  pure elemental function lognormal_quantile(p, meanlog, sdlog, log_p, lower_tail) result(x)
    real(dp),intent(in)::p
    real(dp),intent(in),optional::meanlog,sdlog
    logical,intent(in),optional::log_p,lower_tail
    real(dp)::x,m,s
    m=0.0_dp
    if(present(meanlog))m=meanlog
    s=1.0_dp
    if(present(sdlog))s=sdlog
    x=exp(normal_quantile(p,m,s,log_p,lower_tail))
  end function lognormal_quantile

  pure elemental function logistic_pdf(x, location, scale, log_pdf) result(y)
    real(dp),intent(in)::x
    real(dp),intent(in),optional::location,scale
    logical,intent(in),optional::log_pdf
    real(dp)::y,m,s,z
    logical::lp
    m=0.0_dp
    if(present(location))m=location
    s=1.0_dp
    if(present(scale))s=scale
    lp=.false.
    if(present(log_pdf))lp=log_pdf
    z=(x-m)/s
    y=-z-log(s)-2.0_dp*log(1.0_dp+exp(-z))
    if(.not.lp)y=exp(y)
  end function logistic_pdf
  pure elemental function logistic_cdf(x, location, scale, log_p, lower_tail) result(y)
    real(dp),intent(in)::x
    real(dp),intent(in),optional::location,scale
    logical,intent(in),optional::log_p,lower_tail
    real(dp)::y,m,s
    logical::lp,lt
    m=0.0_dp
    if(present(location))m=location
    s=1.0_dp
    if(present(scale))s=scale
    lp=.false.
    if(present(log_p))lp=log_p
    lt=.true.
    if(present(lower_tail))lt=lower_tail
    y=1.0_dp/(1.0_dp+exp(-(x-m)/s))
    if(.not.lt)y=1.0_dp-y
    if(lp)y=log(y)
  end function logistic_cdf
  pure elemental function logistic_quantile(p, location, scale, log_p, lower_tail) result(x)
    real(dp),intent(in)::p
    real(dp),intent(in),optional::location,scale
    logical,intent(in),optional::log_p,lower_tail
    real(dp)::x,pp,m,s
    logical::lp,lt
    m=0.0_dp
    if(present(location))m=location
    s=1.0_dp
    if(present(scale))s=scale
    lp=.false.
    if(present(log_p))lp=log_p
    lt=.true.
    if(present(lower_tail))lt=lower_tail
    pp=p
    if(lp)pp=exp(pp)
    if(.not.lt)pp=1.0_dp-pp
    x=m+s*log(pp/(1.0_dp-pp))
  end function logistic_quantile

  pure elemental function cauchy_pdf(x, location, scale, log_pdf) result(y)
    real(dp),intent(in)::x
    real(dp),intent(in),optional::location,scale
    logical,intent(in),optional::log_pdf
    real(dp)::y,m,s,z
    logical::lp
    m=0.0_dp
    if(present(location))m=location
    s=1.0_dp
    if(present(scale))s=scale
    lp=.false.
    if(present(log_pdf))lp=log_pdf
    z=(x-m)/s
    y=-log(pi*s)-log(1.0_dp+z*z)
    if(.not.lp)y=exp(y)
  end function cauchy_pdf
  pure elemental function cauchy_cdf(x, location, scale, log_p, lower_tail) result(y)
    real(dp),intent(in)::x
    real(dp),intent(in),optional::location,scale
    logical,intent(in),optional::log_p,lower_tail
    real(dp)::y,m,s
    logical::lp,lt
    m=0.0_dp
    if(present(location))m=location
    s=1.0_dp
    if(present(scale))s=scale
    lp=.false.
    if(present(log_p))lp=log_p
    lt=.true.
    if(present(lower_tail))lt=lower_tail
    y=0.5_dp+atan((x-m)/s)/pi
    if(.not.lt)y=1.0_dp-y
    if(lp)y=log(y)
  end function cauchy_cdf
  pure elemental function cauchy_quantile(p, location, scale, log_p, lower_tail) result(x)
    real(dp),intent(in)::p
    real(dp),intent(in),optional::location,scale
    logical,intent(in),optional::log_p,lower_tail
    real(dp)::x,pp,m,s
    logical::lp,lt
    m=0.0_dp
    if(present(location))m=location
    s=1.0_dp
    if(present(scale))s=scale
    lp=.false.
    if(present(log_p))lp=log_p
    lt=.true.
    if(present(lower_tail))lt=lower_tail
    pp=p
    if(lp)pp=exp(pp)
    if(.not.lt)pp=1.0_dp-pp
    x=m+s*tan(pi*(pp-0.5_dp))
  end function cauchy_quantile

  pure elemental function uniform_pdf(x, lower, upper, log_pdf) result(y)
    real(dp),intent(in)::x
    real(dp),intent(in),optional::lower,upper
    logical,intent(in),optional::log_pdf
    real(dp)::y,a,b
    logical::lp
    a=0.0_dp
    if(present(lower))a=lower
    b=1.0_dp
    if(present(upper))b=upper
    lp=.false.
    if(present(log_pdf))lp=log_pdf
    if(x<a.or.x>b)then
    y=merge(-huge(1.0_dp),0.0_dp,lp)
    else
    y=-log(b-a)
    if(.not.lp)y=exp(y)
    end if
  end function uniform_pdf
  pure elemental function uniform_cdf(x, lower, upper, log_p, lower_tail) result(y)
    real(dp),intent(in)::x
    real(dp),intent(in),optional::lower,upper
    logical,intent(in),optional::log_p,lower_tail
    real(dp)::y,a,b
    logical::lp,lt
    a=0.0_dp
    if(present(lower))a=lower
    b=1.0_dp
    if(present(upper))b=upper
    lp=.false.
    if(present(log_p))lp=log_p
    lt=.true.
    if(present(lower_tail))lt=lower_tail
    y=max(0.0_dp,min(1.0_dp,(x-a)/(b-a)))
    if(.not.lt)y=1.0_dp-y
    if(lp)y=log(y)
  end function uniform_cdf
  pure elemental function uniform_quantile(p, lower, upper, log_p, lower_tail) result(x)
    real(dp),intent(in)::p
    real(dp),intent(in),optional::lower,upper
    logical,intent(in),optional::log_p,lower_tail
    real(dp)::x,pp,a,b
    logical::lp,lt
    a=0.0_dp
    if(present(lower))a=lower
    b=1.0_dp
    if(present(upper))b=upper
    lp=.false.
    if(present(log_p))lp=log_p
    lt=.true.
    if(present(lower_tail))lt=lower_tail
    pp=p
    if(lp)pp=exp(pp)
    if(.not.lt)pp=1.0_dp-pp
    x=a+(b-a)*pp
  end function uniform_quantile

  pure elemental function weibull_pdf(x, shape, scale, log_pdf) result(y)
    real(dp),intent(in)::x
    real(dp),intent(in),optional::shape,scale
    logical,intent(in),optional::log_pdf
    real(dp)::y,a,s
    logical::lp
    a=1.0_dp
    if(present(shape))a=shape
    s=1.0_dp
    if(present(scale))s=scale
    lp=.false.
    if(present(log_pdf))lp=log_pdf
    if (x < 0.0_dp) then
      y = merge(-huge(1.0_dp), 0.0_dp, lp)
      return
    end if
    y=log(a/s)+(a-1.0_dp)*log(x/s)-(x/s)**a
    if(.not.lp)y=exp(y)
  end function weibull_pdf
  pure elemental function weibull_cdf(x, shape, scale, log_p, lower_tail) result(y)
    real(dp),intent(in)::x
    real(dp),intent(in),optional::shape,scale
    logical,intent(in),optional::log_p,lower_tail
    real(dp)::y,a,s
    logical::lp,lt
    a=1.0_dp
    if(present(shape))a=shape
    s=1.0_dp
    if(present(scale))s=scale
    lp=.false.
    if(present(log_p))lp=log_p
    lt=.true.
    if(present(lower_tail))lt=lower_tail
    if (x < 0.0_dp) then
      y = merge(0.0_dp, 1.0_dp, lt)
    else
      y=1.0_dp-exp(-(x/s)**a)
      if(.not.lt)y=1.0_dp-y
    end if
    if(lp)y=log(y)
  end function weibull_cdf
  pure elemental function weibull_quantile(p, shape, scale, log_p, lower_tail) result(x)
    real(dp),intent(in)::p
    real(dp),intent(in),optional::shape,scale
    logical,intent(in),optional::log_p,lower_tail
    real(dp)::x,pp,a,s
    logical::lp,lt
    a=1.0_dp
    if(present(shape))a=shape
    s=1.0_dp
    if(present(scale))s=scale
    lp=.false.
    if(present(log_p))lp=log_p
    lt=.true.
    if(present(lower_tail))lt=lower_tail
    pp=p
    if(lp)pp=exp(pp)
    if(.not.lt)pp=1.0_dp-pp
    x=s*(-log(1.0_dp-pp))**(1.0_dp/a)
  end function weibull_quantile
end module vares_special
