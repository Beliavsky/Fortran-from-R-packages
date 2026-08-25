module mnormt_special
  use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
  implicit none
  private
  integer, parameter, public :: dp = kind(1.0d0)
  real(dp), parameter, public :: pi_dp = acos(-1.0_dp)
  public :: normal_pdf, normal_cdf, normal_quantile, reg_incomplete_beta
  public :: student_t_pdf, student_t_cdf, gamma_rng, chisq_rng, normal_rng
contains

  pure real(dp) function normal_pdf(x) result(f)
    real(dp), intent(in) :: x
    f = exp(-0.5_dp*x*x) / sqrt(2.0_dp*pi_dp)
  end function normal_pdf

  pure real(dp) function normal_cdf(x) result(p)
    real(dp), intent(in) :: x
    p = 0.5_dp * erfc(-x / sqrt(2.0_dp))
  end function normal_cdf

  pure real(dp) function normal_quantile(p) result(x)
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
       4.374664141464968e+00_dp, 2.938163982698783e+00_dp]
    real(dp), parameter :: d(4) = [ &
       7.784695709041462e-03_dp, 3.224671290700398e-01_dp, &
       2.445134137142996e+00_dp, 3.754408661907416e+00_dp]
    real(dp), parameter :: plow = 0.02425_dp, phigh = 1.0_dp-plow
    real(dp) :: q, r
    if (p <= 0.0_dp) then
      x = -huge(1.0_dp)
    else if (p >= 1.0_dp) then
      x = huge(1.0_dp)
    else if (p < plow) then
      q = sqrt(-2.0_dp*log(p))
      x = (((((c(1)*q+c(2))*q+c(3))*q+c(4))*q+c(5))*q+c(6)) / &
          ((((d(1)*q+d(2))*q+d(3))*q+d(4))*q+1.0_dp)
    else if (p <= phigh) then
      q = p-0.5_dp
      r = q*q
      x = (((((a(1)*r+a(2))*r+a(3))*r+a(4))*r+a(5))*r+a(6))*q / &
          (((((b(1)*r+b(2))*r+b(3))*r+b(4))*r+b(5))*r+1.0_dp)
    else
      q = sqrt(-2.0_dp*log(1.0_dp-p))
      x = -(((((c(1)*q+c(2))*q+c(3))*q+c(4))*q+c(5))*q+c(6)) / &
           ((((d(1)*q+d(2))*q+d(3))*q+d(4))*q+1.0_dp)
    end if
    if (p > 0.0_dp .and. p < 1.0_dp) then
      x = x - (normal_cdf(x)-p) / max(normal_pdf(x), tiny(1.0_dp))
    end if
  end function normal_quantile

  pure real(dp) function betacf(a, b, x) result(cf)
    real(dp), intent(in) :: a, b, x
    integer, parameter :: maxit = 300
    real(dp), parameter :: eps = 3.0e-15_dp, fpmin = 1.0e-300_dp
    integer :: m, m2
    real(dp) :: aa, c, d, del, h, qab, qam, qap
    qab = a+b; qap = a+1.0_dp; qam = a-1.0_dp
    c = 1.0_dp
    d = 1.0_dp - qab*x/qap
    if (abs(d) < fpmin) d = fpmin
    d = 1.0_dp/d; h = d
    do m=1,maxit
      m2 = 2*m
      aa = real(m,dp)*(b-real(m,dp))*x / &
           ((qam+real(m2,dp))*(a+real(m2,dp)))
      d = 1.0_dp + aa*d; if (abs(d)<fpmin) d=fpmin
      c = 1.0_dp + aa/c; if (abs(c)<fpmin) c=fpmin
      d = 1.0_dp/d; h = h*d*c
      aa = -(a+real(m,dp))*(qab+real(m,dp))*x / &
           ((a+real(m2,dp))*(qap+real(m2,dp)))
      d = 1.0_dp + aa*d; if (abs(d)<fpmin) d=fpmin
      c = 1.0_dp + aa/c; if (abs(c)<fpmin) c=fpmin
      d = 1.0_dp/d; del=d*c; h=h*del
      if (abs(del-1.0_dp) <= eps) exit
    end do
    cf = h
  end function betacf

  pure real(dp) function reg_incomplete_beta(x, a, b) result(p)
    real(dp), intent(in) :: x, a, b
    real(dp) :: bt
    if (a <= 0.0_dp .or. b <= 0.0_dp) then
      p = ieee_value(1.0_dp, ieee_quiet_nan); return
    end if
    if (x <= 0.0_dp) then
      p=0.0_dp; return
    else if (x >= 1.0_dp) then
      p=1.0_dp; return
    end if
    bt = exp(log_gamma(a+b)-log_gamma(a)-log_gamma(b)+a*log(x)+b*log(1.0_dp-x))
    if (x < (a+1.0_dp)/(a+b+2.0_dp)) then
      p = bt*betacf(a,b,x)/a
    else
      p = 1.0_dp-bt*betacf(b,a,1.0_dp-x)/b
    end if
    p = min(1.0_dp,max(0.0_dp,p))
  end function reg_incomplete_beta

  pure real(dp) function student_t_pdf(x, df) result(f)
    real(dp), intent(in) :: x, df
    if (df <= 0.0_dp) then
      f = ieee_value(1.0_dp, ieee_quiet_nan)
    else
      f = exp(log_gamma(0.5_dp*(df+1.0_dp))-log_gamma(0.5_dp*df)) / &
          sqrt(df*pi_dp) * (1.0_dp+x*x/df)**(-0.5_dp*(df+1.0_dp))
    end if
  end function student_t_pdf

  pure real(dp) function student_t_cdf(x, df) result(p)
    real(dp), intent(in) :: x, df
    real(dp) :: ib, z
    if (df <= 0.0_dp) then
      p = ieee_value(1.0_dp, ieee_quiet_nan); return
    end if
    if (abs(x) <= tiny(1.0_dp)) then
      p=0.5_dp; return
    end if
    z = df/(df+x*x)
    ib = reg_incomplete_beta(z,0.5_dp*df,0.5_dp)
    if (x > 0.0_dp) then
      p = 1.0_dp-0.5_dp*ib
    else
      p = 0.5_dp*ib
    end if
  end function student_t_cdf

  real(dp) function normal_rng() result(z)
    real(dp) :: u1,u2
    call random_number(u1); call random_number(u2)
    u1=max(u1,tiny(1.0_dp))
    z=sqrt(-2.0_dp*log(u1))*cos(2.0_dp*pi_dp*u2)
  end function normal_rng

  recursive real(dp) function gamma_rng(shape) result(g)
    real(dp), intent(in) :: shape
    real(dp) :: d,c,x,v,u
    if (shape <= 0.0_dp) then
      g = ieee_value(1.0_dp, ieee_quiet_nan); return
    end if
    if (shape < 1.0_dp) then
      call random_number(u)
      g = gamma_rng(shape+1.0_dp)*u**(1.0_dp/shape)
      return
    end if
    d=shape-1.0_dp/3.0_dp; c=1.0_dp/sqrt(9.0_dp*d)
    do
      do
        x=normal_rng(); v=1.0_dp+c*x
        if (v>0.0_dp) exit
      end do
      v=v*v*v; call random_number(u)
      if (u < 1.0_dp-0.0331_dp*x**4) exit
      if (log(max(u,tiny(1.0_dp))) < 0.5_dp*x*x+d*(1.0_dp-v+log(v))) exit
    end do
    g=d*v
  end function gamma_rng

  real(dp) function chisq_rng(df) result(x)
    real(dp), intent(in) :: df
    x = 2.0_dp*gamma_rng(0.5_dp*df)
  end function chisq_rng
end module mnormt_special
