module kdensity_math
  use kdensity_kinds, only : dp
  use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan, ieee_is_finite
  implicit none
  private
  real(dp), parameter, public :: kd_pi = acos(-1.0_dp)
  public :: normal_pdf, normal_cdf, normal_quantile, gamma_pdf, beta_pdf
  public :: mean_value, sample_sd, quantile_type7, adaptive_integral, golden_minimize
  public :: finite_value, nan_value, log_gamma_pdf, log_beta_pdf

  abstract interface
    function scalar_function(x) result(y)
      import dp
      real(dp), intent(in) :: x
      real(dp) :: y
    end function scalar_function
  end interface

contains

  pure elemental function normal_pdf(x) result(y)
    real(dp), intent(in) :: x
    real(dp) :: y
    y = exp(-0.5_dp*x*x) / sqrt(2.0_dp*kd_pi)
  end function normal_pdf

  pure elemental function normal_cdf(x) result(y)
    real(dp), intent(in) :: x
    real(dp) :: y
    y = 0.5_dp * erfc(-x / sqrt(2.0_dp))
  end function normal_cdf

  pure elemental function normal_quantile(p) result(x)
    real(dp), intent(in) :: p
    real(dp) :: x, q, r
    real(dp), parameter :: a(6) = [ &
      -3.969683028665376e1_dp, 2.209460984245205e2_dp, &
      -2.759285104469687e2_dp, 1.383577518672690e2_dp, &
      -3.066479806614716e1_dp, 2.506628277459239_dp ]
    real(dp), parameter :: b(5) = [ &
      -5.447609879822406e1_dp, 1.615858368580409e2_dp, &
      -1.556989798598866e2_dp, 6.680131188771972e1_dp, &
      -1.328068155288572e1_dp ]
    real(dp), parameter :: c(6) = [ &
      -7.784894002430293e-3_dp, -3.223964580411365e-1_dp, &
      -2.400758277161838_dp, -2.549732539343734_dp, &
       4.374664141464968_dp, 2.938163982698783_dp ]
    real(dp), parameter :: d(4) = [ &
       7.784695709041462e-3_dp, 3.224671290700398e-1_dp, &
       2.445134137142996_dp, 3.754408661907416_dp ]
    real(dp), parameter :: plow = 0.02425_dp, phigh = 1.0_dp - plow
    if (p <= 0.0_dp) then
      x = -huge(1.0_dp)
    else if (p >= 1.0_dp) then
      x = huge(1.0_dp)
    else if (p < plow) then
      q = sqrt(-2.0_dp*log(p))
      x = (((((c(1)*q+c(2))*q+c(3))*q+c(4))*q+c(5))*q+c(6)) / &
          ((((d(1)*q+d(2))*q+d(3))*q+d(4))*q+1.0_dp)
    else if (p > phigh) then
      q = sqrt(-2.0_dp*log(1.0_dp-p))
      x = -(((((c(1)*q+c(2))*q+c(3))*q+c(4))*q+c(5))*q+c(6)) / &
           ((((d(1)*q+d(2))*q+d(3))*q+d(4))*q+1.0_dp)
    else
      q = p - 0.5_dp
      r = q*q
      x = (((((a(1)*r+a(2))*r+a(3))*r+a(4))*r+a(5))*r+a(6))*q / &
          (((((b(1)*r+b(2))*r+b(3))*r+b(4))*r+b(5))*r+1.0_dp)
    end if
  end function normal_quantile

  pure elemental function log_gamma_pdf(x, shape, scale) result(y)
    real(dp), intent(in) :: x, shape, scale
    real(dp) :: y
    if (x < 0.0_dp .or. shape <= 0.0_dp .or. scale <= 0.0_dp) then
      y = -huge(1.0_dp)
    else if (x == 0.0_dp) then
      if (shape < 1.0_dp) then
        y = huge(1.0_dp)
      else if (shape == 1.0_dp) then
        y = -log(scale)
      else
        y = -huge(1.0_dp)
      end if
    else
      y = (shape - 1.0_dp)*log(x) - x/scale - log_gamma(shape) - shape*log(scale)
    end if
  end function log_gamma_pdf

  pure elemental function gamma_pdf(x, shape, scale) result(y)
    real(dp), intent(in) :: x, shape, scale
    real(dp) :: y, ly
    ly = log_gamma_pdf(x, shape, scale)
    if (ly <= log(tiny(1.0_dp))) then
      y = 0.0_dp
    else if (ly >= log(huge(1.0_dp))) then
      y = huge(1.0_dp)
    else
      y = exp(ly)
    end if
  end function gamma_pdf

  pure elemental function log_beta_pdf(x, a, b) result(y)
    real(dp), intent(in) :: x, a, b
    real(dp) :: y
    if (a <= 0.0_dp .or. b <= 0.0_dp .or. x < 0.0_dp .or. x > 1.0_dp) then
      y = -huge(1.0_dp)
    else if (x == 0.0_dp) then
      if (a < 1.0_dp) then
        y = huge(1.0_dp)
      else if (a == 1.0_dp) then
        y = -log_gamma(a)-log_gamma(b)+log_gamma(a+b)
      else
        y = -huge(1.0_dp)
      end if
    else if (x == 1.0_dp) then
      if (b < 1.0_dp) then
        y = huge(1.0_dp)
      else if (b == 1.0_dp) then
        y = -log_gamma(a)-log_gamma(b)+log_gamma(a+b)
      else
        y = -huge(1.0_dp)
      end if
    else
      y = (a-1.0_dp)*log(x) + (b-1.0_dp)*log(1.0_dp-x) + &
          log_gamma(a+b)-log_gamma(a)-log_gamma(b)
    end if
  end function log_beta_pdf

  pure elemental function beta_pdf(x, a, b) result(y)
    real(dp), intent(in) :: x, a, b
    real(dp) :: y, ly
    ly = log_beta_pdf(x, a, b)
    if (ly <= log(tiny(1.0_dp))) then
      y = 0.0_dp
    else if (ly >= log(huge(1.0_dp))) then
      y = huge(1.0_dp)
    else
      y = exp(ly)
    end if
  end function beta_pdf

  pure function mean_value(x) result(m)
    real(dp), intent(in) :: x(:)
    real(dp) :: m
    if (size(x) == 0) then
      m = nan_value()
    else
      m = sum(x) / real(size(x), dp)
    end if
  end function mean_value

  pure function sample_sd(x) result(s)
    real(dp), intent(in) :: x(:)
    real(dp) :: s, m
    if (size(x) < 2) then
      s = 0.0_dp
    else
      m = mean_value(x)
      s = sqrt(max(0.0_dp, sum((x-m)**2)/real(size(x)-1,dp)))
    end if
  end function sample_sd

  function quantile_type7(x, p) result(q)
    real(dp), intent(in) :: x(:), p
    real(dp) :: q, h, g
    real(dp), allocatable :: z(:)
    integer :: j
    if (size(x) == 0) then
      q = nan_value(); return
    end if
    z = x
    call sort_real(z)
    if (p <= 0.0_dp) then
      q = z(1); return
    else if (p >= 1.0_dp) then
      q = z(size(z)); return
    end if
    h = 1.0_dp + real(size(z)-1,dp)*p
    j = int(floor(h))
    g = h-real(j,dp)
    if (j >= size(z)) then
      q = z(size(z))
    else
      q = (1.0_dp-g)*z(j)+g*z(j+1)
    end if
  end function quantile_type7

  subroutine sort_real(x)
    real(dp), intent(inout) :: x(:)
    integer :: i, j
    real(dp) :: key
    do i=2,size(x)
      key=x(i); j=i-1
      do while (j>=1)
        if (x(j) <= key) exit
        x(j+1)=x(j); j=j-1
      end do
      x(j+1)=key
    end do
  end subroutine sort_real

  function adaptive_integral(f, lower, upper, tol, status) result(value)
    procedure(scalar_function) :: f
    real(dp), intent(in) :: lower, upper, tol
    integer, intent(out), optional :: status
    real(dp) :: value, a, b
    integer :: istat
    istat = 0
    if (is_finite_bound(lower) .and. is_finite_bound(upper)) then
      a=lower; b=upper
      value = adaptive_simpson(finite_map,1.0e-9_dp,1.0_dp-1.0e-9_dp,tol,20,istat)
    else if (.not. is_finite_bound(lower) .and. .not. is_finite_bound(upper)) then
      value = transformed_integral(f,1,tol,istat)
    else if (is_finite_bound(lower)) then
      value = transformed_integral(f,2,tol,istat,lower)
    else
      value = transformed_integral(f,3,tol,istat,upper)
    end if
    if (present(status)) status=istat
  contains
    pure function is_finite_bound(x) result(ok)
      real(dp), intent(in) :: x
      logical :: ok
      ok = finite_value(x) .and. abs(x) < 0.1_dp*huge(1.0_dp)
    end function is_finite_bound
    function finite_map(t) result(y)
      real(dp), intent(in) :: t
      real(dp) :: y, x, angle, jac
      angle = 0.5_dp*kd_pi*t
      x = a + (b-a)*sin(angle)**2
      jac = (b-a)*0.5_dp*kd_pi*sin(kd_pi*t)
      y = f(x)*jac
      if (.not. finite_value(y)) y=0.0_dp
    end function finite_map
  end function adaptive_integral

  function transformed_integral(f, kind_map, tol, status, endpoint) result(value)
    procedure(scalar_function) :: f
    integer, intent(in) :: kind_map
    real(dp), intent(in) :: tol
    integer, intent(inout) :: status
    real(dp), intent(in), optional :: endpoint
    real(dp) :: value
    value = adaptive_simpson(g,1.0e-9_dp,1.0_dp-1.0e-9_dp,tol,20,status)
  contains
    function g(t) result(y)
      real(dp), intent(in) :: t
      real(dp) :: y, x, jac
      select case(kind_map)
      case(1)
        x=tan(kd_pi*(t-0.5_dp)); jac=kd_pi/(cos(kd_pi*(t-0.5_dp))**2)
      case(2)
        x=endpoint+t/(1.0_dp-t); jac=1.0_dp/(1.0_dp-t)**2
      case default
        x=endpoint-(1.0_dp-t)/t; jac=1.0_dp/t**2
      end select
      y=f(x)*jac
      if (.not. finite_value(y)) y=0.0_dp
    end function g
  end function transformed_integral

  recursive function adaptive_simpson(f,a,b,tol,depth,status) result(v)
    procedure(scalar_function) :: f
    real(dp), intent(in) :: a,b,tol
    integer, intent(in) :: depth
    integer, intent(inout) :: status
    real(dp) :: v, c, fa, fb, fc, whole
    c=0.5_dp*(a+b); fa=f(a); fb=f(b); fc=f(c)
    if (.not. finite_value(fa+fb+fc)) then
      status=1; v=0.0_dp; return
    end if
    whole=(b-a)*(fa+4.0_dp*fc+fb)/6.0_dp
    v=refine(a,b,fa,fb,fc,whole,tol,depth)
  contains
    recursive function refine(x0,x1,f0,f1,fm,s,eps,d) result(r)
      real(dp), intent(in) :: x0,x1,f0,f1,fm,s,eps
      integer, intent(in) :: d
      real(dp) :: r,m,lm,rm,flm,frm,sl,sr
      m=0.5_dp*(x0+x1); lm=0.5_dp*(x0+m); rm=0.5_dp*(m+x1)
      flm=f(lm); frm=f(rm)
      if (.not. finite_value(flm+frm)) then
        status=1; r=s; return
      end if
      sl=(m-x0)*(f0+4.0_dp*flm+fm)/6.0_dp
      sr=(x1-m)*(fm+4.0_dp*frm+f1)/6.0_dp
      if (d<=0 .or. abs(sl+sr-s)<=15.0_dp*eps) then
        r=sl+sr+(sl+sr-s)/15.0_dp
      else
        r=refine(x0,m,f0,fm,flm,sl,0.5_dp*eps,d-1)+ &
          refine(m,x1,fm,f1,frm,sr,0.5_dp*eps,d-1)
      end if
    end function refine
  end function adaptive_simpson

  function golden_minimize(f, lower, upper, tol, status) result(xmin)
    procedure(scalar_function) :: f
    real(dp), intent(in) :: lower, upper, tol
    integer, intent(out), optional :: status
    real(dp) :: xmin, a,b,c,d,fc,fd,phi
    integer :: iter
    phi=(sqrt(5.0_dp)-1.0_dp)/2.0_dp
    a=lower; b=upper; c=b-phi*(b-a); d=a+phi*(b-a)
    fc=f(c); fd=f(d)
    do iter=1,300
      if (abs(b-a) <= tol*(1.0_dp+abs(a)+abs(b))) exit
      if (fc < fd) then
        b=d; d=c; fd=fc; c=b-phi*(b-a); fc=f(c)
      else
        a=c; c=d; fc=fd; d=a+phi*(b-a); fd=f(d)
      end if
    end do
    xmin=0.5_dp*(a+b)
    if (present(status)) status=merge(0,1,iter<300)
  end function golden_minimize

  pure elemental function finite_value(x) result(ok)
    real(dp), intent(in) :: x
    logical :: ok
    ok=ieee_is_finite(x)
  end function finite_value

  pure function nan_value() result(x)
    real(dp) :: x
    x=ieee_value(0.0_dp,ieee_quiet_nan)
  end function nan_value

end module kdensity_math
