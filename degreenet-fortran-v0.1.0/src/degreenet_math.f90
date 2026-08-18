! SPDX-License-Identifier: GPL-3.0-or-later
! Based on 'statnet' project software (statnet.org).
module degreenet_math
  use degreenet_kinds, only : dp, pi
  implicit none
  private
  public :: zeta_r, logaddexp, logdiffexp, geom_pmf, geom_sf
  public :: poisson_logpmf, poisson_pmf, poisson_cdf
  public :: nbinom_logpmf, nbinom_pmf, nbinom_cdf, nbinom_sf
  public :: gauss_hermite, adaptive_simpson

  abstract interface
    function scalar_fun(x) result(y)
      import dp
      real(dp), intent(in) :: x
      real(dp) :: y
    end function scalar_fun
  end interface

contains
  real(dp) function zeta_r(s) result(z)
    real(dp), intent(in) :: s
    real(dp), parameter :: b2(8) = [1.0_dp/6.0_dp, -1.0_dp/30.0_dp, &
      1.0_dp/42.0_dp, -1.0_dp/30.0_dp, 5.0_dp/66.0_dp, &
      -691.0_dp/2730.0_dp, 7.0_dp/6.0_dp, -3617.0_dp/510.0_dp]
    integer, parameter :: a = 12
    integer :: m, n, m2
    real(dp) :: a2, p, tail
    if (s <= 1.0_dp) then
      z = huge(1.0_dp)
      return
    end if
    a2 = real(a*a, dp)
    p = s/(2.0_dp*a2)
    tail = 1.0_dp/(s-1.0_dp) + 0.5_dp/real(a,dp) + b2(1)*p
    do m = 1, 7
      m2 = 2*m + 2
      p = p*(s+real(m2-3,dp))*(s+real(m2-2,dp)) / &
          (real(m2-1,dp)*real(m2,dp)*a2)
      tail = tail + p*b2(m+1)
    end do
    z = 1.0_dp + tail/real(a,dp)**(s-1.0_dp)
    do n = 2, a-1
      z = z + real(n,dp)**(-s)
    end do
  end function zeta_r

  pure real(dp) function logaddexp(a,b) result(c)
    real(dp), intent(in) :: a,b
    if (a > b) then
      c = a + log(1.0_dp + exp(b-a))
    else
      c = b + log(1.0_dp + exp(a-b))
    end if
  end function logaddexp

  pure real(dp) function logdiffexp(a,b) result(c)
    real(dp), intent(in) :: a,b
    if (b >= a) then
      c = -huge(1.0_dp)
    else
      c = a + log(1.0_dp-exp(b-a))
    end if
  end function logdiffexp

  pure real(dp) function geom_pmf(k,p) result(v)
    integer, intent(in) :: k
    real(dp), intent(in) :: p
    if (k < 0 .or. p <= 0.0_dp .or. p > 1.0_dp) then
      v = 0.0_dp
    else if (p >= 1.0_dp) then
      v = merge(1.0_dp,0.0_dp,k==0)
    else
      v = p*(1.0_dp-p)**k
    end if
  end function geom_pmf

  pure real(dp) function geom_sf(k,p) result(v)
    integer, intent(in) :: k
    real(dp), intent(in) :: p
    if (k < 0) then
      v = 1.0_dp
    else if (p <= 0.0_dp) then
      v = 1.0_dp
    else if (p >= 1.0_dp) then
      v = 0.0_dp
    else
      v = (1.0_dp-p)**(k+1)
    end if
  end function geom_sf

  pure real(dp) function poisson_logpmf(k,lambda) result(v)
    integer, intent(in) :: k
    real(dp), intent(in) :: lambda
    if (k < 0 .or. lambda < 0.0_dp) then
      v = -huge(1.0_dp)
    else if (lambda <= 0.0_dp) then
      v = merge(0.0_dp,-huge(1.0_dp),k==0)
    else
      v = real(k,dp)*log(lambda)-lambda-log_gamma(real(k+1,dp))
    end if
  end function poisson_logpmf

  pure real(dp) function poisson_pmf(k,lambda) result(v)
    integer, intent(in) :: k
    real(dp), intent(in) :: lambda
    v = exp(poisson_logpmf(k,lambda))
  end function poisson_pmf

  real(dp) function poisson_cdf(k,lambda) result(v)
    integer, intent(in) :: k
    real(dp), intent(in) :: lambda
    integer :: j
    real(dp) :: term
    if (k < 0) then
      v=0.0_dp; return
    end if
    if (lambda <= 0.0_dp) then
      v=1.0_dp; return
    end if
    term=exp(-lambda); v=term
    do j=1,k
      term=term*lambda/real(j,dp); v=v+term
    end do
    v=min(1.0_dp,v)
  end function poisson_cdf

  pure real(dp) function nbinom_logpmf(k,size,p) result(v)
    integer, intent(in) :: k
    real(dp), intent(in) :: size,p
    if (k<0 .or. size<=0.0_dp .or. p<=0.0_dp .or. p>1.0_dp) then
      v=-huge(1.0_dp)
    else if (p>=1.0_dp) then
      v=merge(0.0_dp,-huge(1.0_dp),k==0)
    else
      v=log_gamma(real(k,dp)+size)-log_gamma(size)-log_gamma(real(k+1,dp)) + &
        size*log(p)+real(k,dp)*log(1.0_dp-p)
    end if
  end function nbinom_logpmf

  pure real(dp) function nbinom_pmf(k,size,p) result(v)
    integer,intent(in)::k
    real(dp),intent(in)::size,p
    v=exp(nbinom_logpmf(k,size,p))
  end function nbinom_pmf

  real(dp) function nbinom_cdf(k,size,p) result(v)
    integer,intent(in)::k
    real(dp),intent(in)::size,p
    integer::j
    real(dp)::term
    if(k<0) then; v=0.0_dp; return; end if
    if(size<=0.0_dp.or.p<=0.0_dp) then; v=0.0_dp; return; end if
    if(p>=1.0_dp) then; v=1.0_dp; return; end if
    term=p**size; v=term
    do j=1,k
      term=term*(real(j-1,dp)+size)/real(j,dp)*(1.0_dp-p)
      v=v+term
    end do
    v=min(1.0_dp,max(0.0_dp,v))
  end function nbinom_cdf

  real(dp) function nbinom_sf(k,size,p) result(v)
    integer,intent(in)::k
    real(dp),intent(in)::size,p
    v=max(0.0_dp,1.0_dp-nbinom_cdf(k,size,p))
  end function nbinom_sf

  subroutine gauss_hermite(n,x,w,ok)
    integer,intent(in)::n
    real(dp),intent(out)::x(n),w(n)
    logical,intent(out),optional::ok
    integer::i,j,m,it
    real(dp)::z,z1,p1,p2,p3,pp
    logical::conv
    x=0.0_dp; w=0.0_dp; m=(n+1)/2; conv=.true.; z=0.0_dp
    do i=1,m
      select case(i)
      case(1)
        z=sqrt(real(2*n+1,dp))-2.0_dp*real(2*n+1,dp)**(-1.0_dp/6.0_dp)
      case(2)
        z=z-sqrt(real(n,dp))/z
      case(3,4)
        z=1.9_dp*z-0.9_dp*x(max(1,i-2))
      case default
        z=2.0_dp*z-x(max(1,i-2))
      end select
      do it=1,100
        z1=z; p1=1.0_dp/pi**0.4_dp; p2=0.0_dp
        do j=1,n
          p3=p2; p2=p1
          p1=z*sqrt(2.0_dp/real(j,dp))*p2 - &
             sqrt(real(j-1,dp)/real(j,dp))*p3
        end do
        pp=sqrt(real(2*n,dp))*p2
        z=z1-p1/pp
        if(abs(z-z1)<=1.0e-14_dp) exit
      end do
      if(it>100) conv=.false.
      x(i)=z; x(n+1-i)=-z; w(i)=2.0_dp/(pp*pp); w(n+1-i)=w(i)
    end do
    x=x*sqrt(2.0_dp); w=w/sum(w)
    if(present(ok)) ok=conv
  end subroutine gauss_hermite

  real(dp) function adaptive_simpson(f,a,b,tol,maxdepth) result(ans)
    procedure(scalar_fun)::f
    real(dp),intent(in)::a,b,tol
    integer,intent(in),optional::maxdepth
    integer::md
    real(dp)::fa,fb,fc,s
    md=20; if(present(maxdepth)) md=maxdepth
    fa=f(a); fb=f(b); fc=f(0.5_dp*(a+b))
    s=(b-a)*(fa+4.0_dp*fc+fb)/6.0_dp
    ans=asr(a,b,fa,fb,fc,s,tol,md)
  contains
    recursive real(dp) function asr(aa,bb,faa,fbb,fcc,ss,eps,depth) result(v)
      real(dp),intent(in)::aa,bb,faa,fbb,fcc,ss,eps
      integer,intent(in)::depth
      real(dp)::c,d,e,fd,fe,s1,s2
      c=0.5_dp*(aa+bb); d=0.5_dp*(aa+c); e=0.5_dp*(c+bb)
      fd=f(d); fe=f(e)
      s1=(c-aa)*(faa+4.0_dp*fd+fcc)/6.0_dp
      s2=(bb-c)*(fcc+4.0_dp*fe+fbb)/6.0_dp
      if(depth<=0.or.abs(s1+s2-ss)<=15.0_dp*eps) then
        v=s1+s2+(s1+s2-ss)/15.0_dp
      else
        v=asr(aa,c,faa,fcc,fd,s1,eps/2.0_dp,depth-1)+ &
          asr(c,bb,fcc,fbb,fe,s2,eps/2.0_dp,depth-1)
      end if
    end function asr
  end function adaptive_simpson
end module degreenet_math
