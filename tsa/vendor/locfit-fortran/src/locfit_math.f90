! Derived from locfit src/math.c and src/prob.c, GPL-2-or-later.
module locfit_math
  use locfit_kinds, only : dp
  use locfit_constants, only : sqrt2, sqrpi
  implicit none
  private
  public :: lf_exp, expit, logit, normal_cdf, gamma_p, beta_i, ptail

contains

  pure real(dp) function lf_exp(x) result(y)
    real(dp), intent(in) :: x
    if (x > 700.0_dp) then
      y = 1.014232054735004e304_dp
    else
      y = exp(x)
    end if
  end function lf_exp

  pure real(dp) function expit(x) result(y)
    real(dp), intent(in) :: x
    real(dp) :: u
    if (x < 0.0_dp) then
      u = exp(x)
      y = u/(1.0_dp+u)
    else
      y = 1.0_dp/(1.0_dp+exp(-x))
    end if
  end function expit

  pure real(dp) function logit(x) result(y)
    real(dp), intent(in) :: x
    y = log(x/(1.0_dp-x))
  end function logit

  pure real(dp) function normal_cdf(x, mu, sigma) result(p)
    real(dp), intent(in) :: x
    real(dp), intent(in), optional :: mu, sigma
    real(dp) :: m, s, z
    m = 0.0_dp; if (present(mu)) m = mu
    s = 1.0_dp; if (present(sigma)) s = sigma
    z = (x-m)/s
    if (z >= 0.0_dp) then
      p = 0.5_dp*(1.0_dp+erf(z/sqrt2))
    else
      p = 0.5_dp*erfc(-z/sqrt2)
    end if
  end function normal_cdf

  pure real(dp) function ptail(x) result(z)
    real(dp), intent(in) :: x
    real(dp) :: y, f0
    integer :: j
    y = -1.0_dp/x
    z = y
    j = 0
    do
      f0 = -real(2*j+1,dp)/(x*x)
      if (abs(f0) >= 1.0_dp .or. abs(y) <= 1.0e-10_dp*abs(z)) exit
      y = y*f0
      z = z+y
      j = j+1
    end do
  end function ptail

  pure real(dp) function gamma_p(x, a) result(p)
    ! Regularized lower incomplete gamma P(a,x), matching upstream igamma(x,a).
    real(dp), intent(in) :: x, a
    integer, parameter :: itmax = 10000
    real(dp), parameter :: eps = 10.0_dp*epsilon(1.0_dp), fpmin = tiny(1.0_dp)/eps
    real(dp) :: ap, del, sumv, b, c, d, h, an
    integer :: n
    if (x <= 0.0_dp .or. a <= 0.0_dp) then
      p = 0.0_dp
      return
    end if
    if (x < a+1.0_dp) then
      ap = a
      sumv = 1.0_dp/a
      del = sumv
      do n = 1, itmax
        ap = ap+1.0_dp
        del = del*x/ap
        sumv = sumv+del
        if (abs(del) <= abs(sumv)*eps) exit
      end do
      p = sumv*exp(-x+a*log(x)-log_gamma(a))
    else
      b = x+1.0_dp-a
      c = 1.0_dp/fpmin
      d = 1.0_dp/b
      h = d
      do n = 1, itmax
        an = -real(n,dp)*(real(n,dp)-a)
        b = b+2.0_dp
        d = an*d+b
        if (abs(d) < fpmin) d = fpmin
        c = b+an/c
        if (abs(c) < fpmin) c = fpmin
        d = 1.0_dp/d
        del = d*c
        h = h*del
        if (abs(del-1.0_dp) <= eps) exit
      end do
      p = 1.0_dp-exp(-x+a*log(x)-log_gamma(a))*h
    end if
    p = min(1.0_dp,max(0.0_dp,p))
  end function gamma_p

  pure real(dp) function beta_cf(a,b,x) result(cf)
    real(dp), intent(in) :: a,b,x
    integer, parameter :: maxit=10000
    real(dp), parameter :: eps=10.0_dp*epsilon(1.0_dp), fpmin=tiny(1.0_dp)/eps
    integer :: m, m2
    real(dp) :: aa,c,d,del,h,qab,qam,qap
    qab=a+b; qap=a+1.0_dp; qam=a-1.0_dp
    c=1.0_dp
    d=1.0_dp-qab*x/qap
    if(abs(d)<fpmin)d=fpmin
    d=1.0_dp/d; h=d
    do m=1,maxit
      m2=2*m
      aa=real(m,dp)*(b-real(m,dp))*x/((qam+real(m2,dp))*(a+real(m2,dp)))
      d=1.0_dp+aa*d; if(abs(d)<fpmin)d=fpmin
      c=1.0_dp+aa/c; if(abs(c)<fpmin)c=fpmin
      d=1.0_dp/d; h=h*d*c
      aa=-(a+real(m,dp))*(qab+real(m,dp))*x/((a+real(m2,dp))*(qap+real(m2,dp)))
      d=1.0_dp+aa*d; if(abs(d)<fpmin)d=fpmin
      c=1.0_dp+aa/c; if(abs(c)<fpmin)c=fpmin
      d=1.0_dp/d; del=d*c; h=h*del
      if(abs(del-1.0_dp)<=eps)exit
    end do
    cf=h
  end function beta_cf

  pure real(dp) function beta_i(x,a,b) result(v)
    real(dp), intent(in) :: x,a,b
    real(dp) :: bt
    if(x<=0.0_dp)then
      v=0.0_dp; return
    else if(x>=1.0_dp)then
      v=1.0_dp; return
    end if
    bt=exp(log_gamma(a+b)-log_gamma(a)-log_gamma(b)+a*log(x)+b*log(1.0_dp-x))
    if(x<(a+1.0_dp)/(a+b+2.0_dp))then
      v=bt*beta_cf(a,b,x)/a
    else
      v=1.0_dp-bt*beta_cf(b,a,1.0_dp-x)/b
    end if
    v=min(1.0_dp,max(0.0_dp,v))
  end function beta_i

end module locfit_math
