! SPDX-License-Identifier: GPL-2.0-only OR GPL-3.0-only
! Modern Fortran translation of computational code from the R package sn.
module sn_univariate
  use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan, ieee_positive_inf
  use sn_kinds, only : dp, pi, sqrt_two_pi, tiny_dp, huge_dp
  use sn_math, only : normal_pdf, normal_logpdf, normal_cdf, normal_logcdf, &
                      normal_quantile, student_t_pdf, student_t_logpdf, &
                      student_t_cdf, adaptive_simpson, clamp_probability, &
                      finite_real
  use sn_mvn, only : bivariate_normal_cdf
  use sn_rng, only : sn_rng_state
  use sn_status, only : sn_ok, sn_invalid_argument, sn_no_convergence
  implicit none
  private

  type, public :: sn_uv_params
    real(dp) :: xi = 0.0_dp
    real(dp) :: omega = 1.0_dp
    real(dp) :: alpha = 0.0_dp
    real(dp) :: tau = 0.0_dp
  end type sn_uv_params

  type, public :: st_uv_params
    real(dp) :: xi = 0.0_dp
    real(dp) :: omega = 1.0_dp
    real(dp) :: alpha = 0.0_dp
    real(dp) :: nu = 10.0_dp
  end type st_uv_params

  type, public :: fournum_result
    real(dp) :: median = 0.0_dp
    real(dp) :: q_deviation = 0.0_dp
    real(dp) :: galton_bowley = 0.0_dp
    real(dp) :: moors = 0.0_dp
    integer :: status = sn_ok
  end type fournum_result

  public :: dsn, psn, qsn, rsn
  public :: dst, pst, qst, rst
  public :: dsc, psc, qsc, rsc
  public :: owen_t, zeta, sn_cumulants, st_cumulants
  public :: delta_from_alpha, alpha_from_delta, b_nu
  public :: dp_to_cp_sn, cp_to_dp_sn, mode_sn, mode_st, mode_sc
  public :: fournum, quantile_type7

contains

  pure real(dp) function delta_from_alpha(alpha) result(delta)
    real(dp), intent(in) :: alpha
    if (.not. finite_real(alpha)) then
      delta = sign(1.0_dp,alpha)
    else
      delta = alpha/sqrt(1.0_dp+alpha*alpha)
    end if
  end function delta_from_alpha

  pure real(dp) function alpha_from_delta(delta) result(alpha)
    real(dp), intent(in) :: delta
    if (abs(delta) >= 1.0_dp) then
      alpha = sign(huge_dp,delta)
    else
      alpha = delta/sqrt(max(tiny_dp,1.0_dp-delta*delta))
    end if
  end function alpha_from_delta

  real(dp) function owen_t(h, a, jmax, cut_point) result(value)
    real(dp), intent(in) :: h, a
    integer, intent(in), optional :: jmax
    real(dp), intent(in), optional :: cut_point
    integer :: jm
    real(dp) :: cp, aa, ah

    jm = 50
    if (present(jmax)) jm = max(4,jmax)
    cp = 8.0_dp
    if (present(cut_point)) cp = cut_point
    aa = abs(a)
    ah = abs(h)
    if (.not. finite_real(aa)) then
      value = sign(0.5_dp*normal_cdf(-ah),a)
    else if (aa <= tiny_dp) then
      value = 0.0_dp
    else if (aa <= 1.0_dp) then
      value = sign(t_int(ah,aa,jm,cp),a)
    else
      value = sign(0.5_dp*normal_cdf(ah)+normal_cdf(aa*ah)*(0.5_dp-normal_cdf(ah)) &
                   -t_int(aa*ah,1.0_dp/aa,jm,cp),a)
    end if
  contains
    real(dp) function t_int(x,b,jm0,cp0) result(ans)
      real(dp), intent(in) :: x,b,cp0
      integer, intent(in) :: jm0
      integer :: i
      real(dp) :: term, cumulative, alternating_sum, powb, factorial
      if (x > cp0) then
        ans = atan(b)*exp(-0.5_dp*x*x*b/max(atan(b),tiny_dp))* &
              (1.0_dp+0.00868_dp*(x*b)**4)/(2.0_dp*pi)
        return
      end if
      cumulative = 0.0_dp
      alternating_sum = 0.0_dp
      powb = b
      factorial = 1.0_dp
      do i=0,jm0
        if (i > 0) factorial = factorial*real(i,dp)
        term = x**(2*i)/(2.0_dp**i*factorial)
        cumulative = cumulative+term
        alternating_sum = alternating_sum + (-1.0_dp)**i* &
          (1.0_dp-exp(-0.5_dp*x*x)*cumulative)*powb/real(2*i+1,dp)
        powb = powb*b*b
      end do
      ans = (atan(b)-alternating_sum)/(2.0_dp*pi)
    end function t_int
  end function owen_t

  real(dp) function dsn(x, xi, omega, alpha, tau, log_pdf) result(value)
    real(dp), intent(in) :: x
    real(dp), intent(in), optional :: xi, omega, alpha, tau
    logical, intent(in), optional :: log_pdf
    real(dp) :: loc, scale, shape, ext, z, lp
    logical :: give_log

    loc = 0.0_dp
    scale = 1.0_dp
    shape = 0.0_dp
    ext = 0.0_dp
    give_log = .false.
    if (present(xi)) loc = xi
    if (present(omega)) scale = omega
    if (present(alpha)) shape = alpha
    if (present(tau)) ext = tau
    if (present(log_pdf)) give_log = log_pdf
    if (scale <= 0.0_dp) then
      value = ieee_value(0.0_dp,ieee_quiet_nan)
      return
    end if
    z = (x-loc)/scale
    if (.not. finite_real(z)) then
      lp = -ieee_value(0.0_dp,ieee_positive_inf)
    else if (finite_real(shape)) then
      lp = normal_logpdf(z)-log(scale)+ &
           normal_logcdf(ext*sqrt(1.0_dp+shape*shape)+shape*z)-normal_logcdf(ext)
    else
      if (sign(1.0_dp,shape)*z+ext > 0.0_dp) then
        lp = normal_logpdf(z)-log(scale)-normal_logcdf(ext)
      else
        lp = -ieee_value(0.0_dp,ieee_positive_inf)
      end if
    end if
    if (give_log) then
      value = lp
    else
      value = exp(lp)
    end if
  end function dsn

  real(dp) function psn(x, xi, omega, alpha, tau) result(value)
    real(dp), intent(in) :: x
    real(dp), intent(in), optional :: xi, omega, alpha, tau
    real(dp) :: loc, scale, shape, ext, z, delta, den

    loc = 0.0_dp
    scale = 1.0_dp
    shape = 0.0_dp
    ext = 0.0_dp
    if (present(xi)) loc = xi
    if (present(omega)) scale = omega
    if (present(alpha)) shape = alpha
    if (present(tau)) ext = tau
    if (scale <= 0.0_dp) then
      value = ieee_value(0.0_dp,ieee_quiet_nan)
      return
    end if
    z = (x-loc)/scale
    if (z <= -38.0_dp) then
      value = 0.0_dp
      return
    else if (z >= 38.0_dp) then
      value = 1.0_dp
      return
    end if
    if (abs(ext) <= tiny_dp) then
      value = clamp_probability(normal_cdf(z)-2.0_dp*owen_t(z,shape))
    else
      delta = delta_from_alpha(shape)
      den = normal_cdf(ext)
      if (den <= tiny_dp) then
        value = ieee_value(0.0_dp,ieee_quiet_nan)
      else
        value = clamp_probability(bivariate_normal_cdf(z,ext,-delta)/den)
      end if
    end if
  end function psn

  real(dp) function qsn(p, xi, omega, alpha, tau, tol, info) result(value)
    real(dp), intent(in) :: p
    real(dp), intent(in), optional :: xi, omega, alpha, tau, tol
    integer, intent(out), optional :: info
    real(dp) :: loc, scale, shape, ext, eps, lo, hi, mid, pmid
    integer :: iter

    loc = 0.0_dp
    scale = 1.0_dp
    shape = 0.0_dp
    ext = 0.0_dp
    eps = 1.0e-10_dp
    if (present(xi)) loc = xi
    if (present(omega)) scale = omega
    if (present(alpha)) shape = alpha
    if (present(tau)) ext = tau
    if (present(tol)) eps = tol
    if (scale <= 0.0_dp .or. p < 0.0_dp .or. p > 1.0_dp) then
      value = ieee_value(0.0_dp,ieee_quiet_nan)
      if (present(info)) info = sn_invalid_argument
      return
    else if (p <= 0.0_dp) then
      value = -ieee_value(0.0_dp,ieee_positive_inf)
      if (present(info)) info = sn_ok
      return
    else if (p >= 1.0_dp) then
      value = ieee_value(0.0_dp,ieee_positive_inf)
      if (present(info)) info = sn_ok
      return
    end if
    lo = -8.0_dp
    hi = 8.0_dp
    do while (psn(lo,alpha=shape,tau=ext) > p)
      lo = 2.0_dp*lo
    end do
    do while (psn(hi,alpha=shape,tau=ext) < p)
      hi = 2.0_dp*hi
    end do
    do iter=1,200
      mid = 0.5_dp*(lo+hi)
      pmid = psn(mid,alpha=shape,tau=ext)
      if (abs(pmid-p) <= eps) then
        lo = mid
        hi = mid
        exit
      else if (pmid < p) then
        lo = mid
      else
        hi = mid
      end if
      if (abs(hi-lo) <= eps*max(1.0_dp,abs(mid))) exit
    end do
    value = loc+scale*0.5_dp*(lo+hi)
    if (present(info)) then
      if (iter > 200) then
        info = sn_no_convergence
      else
        info = sn_ok
      end if
    end if
  end function qsn

  subroutine rsn(rng, x, xi, omega, alpha, tau, info)
    type(sn_rng_state), intent(inout) :: rng
    real(dp), intent(out) :: x(:)
    real(dp), intent(in), optional :: xi, omega, alpha, tau
    integer, intent(out), optional :: info
    real(dp) :: loc, scale, shape, ext, delta, u0, u1
    integer :: i

    loc = 0.0_dp
    scale = 1.0_dp
    shape = 0.0_dp
    ext = 0.0_dp
    if (present(xi)) loc = xi
    if (present(omega)) scale = omega
    if (present(alpha)) shape = alpha
    if (present(tau)) ext = tau
    if (scale <= 0.0_dp) then
      x = ieee_value(0.0_dp,ieee_quiet_nan)
      if (present(info)) info = sn_invalid_argument
      return
    end if
    delta = delta_from_alpha(shape)
    do i=1,size(x)
      if (abs(ext) <= tiny_dp) then
        u0 = abs(rng%normal())
      else
        u0 = rng%truncated_normal_lower(-ext)
      end if
      u1 = rng%normal()
      x(i) = loc+scale*(delta*u0+sqrt(max(0.0_dp,1.0_dp-delta*delta))*u1)
    end do
    if (present(info)) info = sn_ok
  end subroutine rsn

  real(dp) function dst(x, xi, omega, alpha, nu, log_pdf) result(value)
    real(dp), intent(in) :: x
    real(dp), intent(in), optional :: xi, omega, alpha, nu
    logical, intent(in), optional :: log_pdf
    real(dp) :: loc, scale, shape, df, z, lp
    logical :: give_log

    loc = 0.0_dp
    scale = 1.0_dp
    shape = 0.0_dp
    df = 10.0_dp
    give_log = .false.
    if (present(xi)) loc = xi
    if (present(omega)) scale = omega
    if (present(alpha)) shape = alpha
    if (present(nu)) df = nu
    if (present(log_pdf)) give_log = log_pdf
    if (scale <= 0.0_dp .or. df <= 0.0_dp) then
      value = ieee_value(0.0_dp,ieee_quiet_nan)
      return
    end if
    if (df > 1.0e8_dp) then
      value = dsn(x,loc,scale,shape,log_pdf=give_log)
      return
    end if
    if (abs(df-1.0_dp) <= epsilon(1.0_dp)) then
      value = dsc(x,loc,scale,shape,log_pdf=give_log)
      return
    end if
    z = (x-loc)/scale
    lp = log(2.0_dp)+student_t_logpdf(z,df)+ &
         log(max(student_t_cdf(shape*z*sqrt((df+1.0_dp)/(z*z+df)),df+1.0_dp),tiny_dp))-log(scale)
    if (give_log) then
      value = lp
    else
      value = exp(lp)
    end if
  end function dst

  real(dp) function pst(x, xi, omega, alpha, nu, tol) result(value)
    real(dp), intent(in) :: x
    real(dp), intent(in), optional :: xi, omega, alpha, nu, tol
    real(dp) :: loc, scale, shape, df, eps, z, p0, integ, delta

    loc = 0.0_dp
    scale = 1.0_dp
    shape = 0.0_dp
    df = 10.0_dp
    eps = 1.0e-9_dp
    if (present(xi)) loc = xi
    if (present(omega)) scale = omega
    if (present(alpha)) shape = alpha
    if (present(nu)) df = nu
    if (present(tol)) eps = tol
    if (scale <= 0.0_dp .or. df <= 0.0_dp) then
      value = ieee_value(0.0_dp,ieee_quiet_nan)
      return
    end if
    z = (x-loc)/scale
    if (z <= -1.0e6_dp) then
      value = 0.0_dp
      return
    else if (z >= 1.0e6_dp) then
      value = 1.0_dp
      return
    end if
    if (df > 1.0e6_dp) then
      value = psn(z,alpha=shape)
      return
    else if (abs(df-1.0_dp) <= epsilon(1.0_dp)) then
      value = psc(z,alpha=shape)
      return
    else if (abs(shape) <= tiny_dp) then
      value = student_t_cdf(z,df)
      return
    end if
    delta = delta_from_alpha(shape)
    p0 = acos(delta)/pi
    if (abs(z) <= tiny_dp) then
      value = p0
    else if (z > 0.0_dp) then
      integ = adaptive_simpson(std_density,0.0_dp,z,eps,24)
      value = clamp_probability(p0+integ)
    else
      integ = adaptive_simpson(std_density,z,0.0_dp,eps,24)
      value = clamp_probability(p0-integ)
    end if
  contains
    real(dp) function std_density(t) result(y)
      real(dp), intent(in) :: t
      y = dst(t,alpha=shape,nu=df)
    end function std_density
  end function pst

  real(dp) function qst(p, xi, omega, alpha, nu, tol, info) result(value)
    real(dp), intent(in) :: p
    real(dp), intent(in), optional :: xi, omega, alpha, nu, tol
    integer, intent(out), optional :: info
    real(dp) :: loc, scale, shape, df, eps, lo, hi, mid, pmid
    integer :: iter

    loc = 0.0_dp
    scale = 1.0_dp
    shape = 0.0_dp
    df = 10.0_dp
    eps = 1.0e-8_dp
    if (present(xi)) loc = xi
    if (present(omega)) scale = omega
    if (present(alpha)) shape = alpha
    if (present(nu)) df = nu
    if (present(tol)) eps = tol
    if (scale <= 0.0_dp .or. df <= 0.0_dp .or. p < 0.0_dp .or. p > 1.0_dp) then
      value = ieee_value(0.0_dp,ieee_quiet_nan)
      if (present(info)) info = sn_invalid_argument
      return
    else if (p <= 0.0_dp) then
      value = -ieee_value(0.0_dp,ieee_positive_inf)
      if (present(info)) info = sn_ok
      return
    else if (p >= 1.0_dp) then
      value = ieee_value(0.0_dp,ieee_positive_inf)
      if (present(info)) info = sn_ok
      return
    end if
    lo = -8.0_dp
    hi = 8.0_dp
    do while (pst(lo,alpha=shape,nu=df,tol=0.2_dp*eps) > p)
      lo = 2.0_dp*lo
    end do
    do while (pst(hi,alpha=shape,nu=df,tol=0.2_dp*eps) < p)
      hi = 2.0_dp*hi
    end do
    do iter=1,180
      mid = 0.5_dp*(lo+hi)
      pmid = pst(mid,alpha=shape,nu=df,tol=0.2_dp*eps)
      if (abs(pmid-p) <= eps) then
        lo = mid
        hi = mid
        exit
      else if (pmid < p) then
        lo = mid
      else
        hi = mid
      end if
      if (abs(hi-lo) <= eps*max(1.0_dp,abs(mid))) exit
    end do
    value = loc+scale*0.5_dp*(lo+hi)
    if (present(info)) then
      if (iter > 180) then
        info = sn_no_convergence
      else
        info = sn_ok
      end if
    end if
  end function qst

  subroutine rst(rng, x, xi, omega, alpha, nu, info)
    type(sn_rng_state), intent(inout) :: rng
    real(dp), intent(out) :: x(:)
    real(dp), intent(in), optional :: xi, omega, alpha, nu
    integer, intent(out), optional :: info
    real(dp) :: loc, scale, shape, df, v
    integer :: i, ierr

    loc = 0.0_dp
    scale = 1.0_dp
    shape = 0.0_dp
    df = 10.0_dp
    if (present(xi)) loc = xi
    if (present(omega)) scale = omega
    if (present(alpha)) shape = alpha
    if (present(nu)) df = nu
    if (scale <= 0.0_dp .or. df <= 0.0_dp) then
      x = ieee_value(0.0_dp,ieee_quiet_nan)
      if (present(info)) info = sn_invalid_argument
      return
    end if
    call rsn(rng,x,omega=scale,alpha=shape,info=ierr)
    if (ierr /= sn_ok) then
      if (present(info)) info = ierr
      return
    end if
    if (df < 1.0e8_dp) then
      do i=1,size(x)
        v = rng%chi_square(df)/df
        x(i) = loc+x(i)/sqrt(v)
      end do
    else
      x = loc+x
    end if
    if (present(info)) info = sn_ok
  end subroutine rst

  real(dp) function dsc(x, xi, omega, alpha, log_pdf) result(value)
    real(dp), intent(in) :: x
    real(dp), intent(in), optional :: xi, omega, alpha
    logical, intent(in), optional :: log_pdf
    real(dp) :: loc, scale, shape, z, lp
    logical :: give_log

    loc = 0.0_dp
    scale = 1.0_dp
    shape = 0.0_dp
    give_log = .false.
    if (present(xi)) loc = xi
    if (present(omega)) scale = omega
    if (present(alpha)) shape = alpha
    if (present(log_pdf)) give_log = log_pdf
    if (scale <= 0.0_dp) then
      value = ieee_value(0.0_dp,ieee_quiet_nan)
      return
    end if
    z = (x-loc)/scale
    lp = -log(pi*scale)-log(1.0_dp+z*z)+ &
         log(1.0_dp+shape*z/sqrt(1.0_dp+z*z*(1.0_dp+shape*shape)))
    if (give_log) then
      value = lp
    else
      value = exp(lp)
    end if
  end function dsc

  real(dp) function psc(x, xi, omega, alpha) result(value)
    real(dp), intent(in) :: x
    real(dp), intent(in), optional :: xi, omega, alpha
    real(dp) :: loc, scale, shape, z, delta
    loc = 0.0_dp
    scale = 1.0_dp
    shape = 0.0_dp
    if (present(xi)) loc = xi
    if (present(omega)) scale = omega
    if (present(alpha)) shape = alpha
    if (scale <= 0.0_dp) then
      value = ieee_value(0.0_dp,ieee_quiet_nan)
      return
    end if
    z = (x-loc)/scale
    delta = delta_from_alpha(shape)
    value = clamp_probability(atan(z)/pi+acos(delta/sqrt(1.0_dp+z*z))/pi)
  end function psc

  real(dp) function qsc(p, xi, omega, alpha, info) result(value)
    real(dp), intent(in) :: p
    real(dp), intent(in), optional :: xi, omega, alpha
    integer, intent(out), optional :: info
    real(dp) :: loc, scale, shape, u, delta, z
    loc = 0.0_dp
    scale = 1.0_dp
    shape = 0.0_dp
    if (present(xi)) loc = xi
    if (present(omega)) scale = omega
    if (present(alpha)) shape = alpha
    if (scale <= 0.0_dp .or. p < 0.0_dp .or. p > 1.0_dp) then
      value = ieee_value(0.0_dp,ieee_quiet_nan)
      if (present(info)) info = sn_invalid_argument
      return
    else if (p <= 0.0_dp) then
      value = -ieee_value(0.0_dp,ieee_positive_inf)
      if (present(info)) info = sn_ok
      return
    else if (p >= 1.0_dp) then
      value = ieee_value(0.0_dp,ieee_positive_inf)
      if (present(info)) info = sn_ok
      return
    end if
    u = (p-0.5_dp)*pi
    delta = delta_from_alpha(shape)
    z = delta/cos(u)+tan(u)
    value = loc+scale*z
    if (present(info)) info = sn_ok
  end function qsc

  subroutine rsc(rng, x, xi, omega, alpha, info)
    type(sn_rng_state), intent(inout) :: rng
    real(dp), intent(out) :: x(:)
    real(dp), intent(in), optional :: xi, omega, alpha
    integer, intent(out), optional :: info
    real(dp) :: loc, scale, shape
    integer :: i, ierr
    loc = 0.0_dp
    scale = 1.0_dp
    shape = 0.0_dp
    if (present(xi)) loc = xi
    if (present(omega)) scale = omega
    if (present(alpha)) shape = alpha
    call rsn(rng,x,omega=scale,alpha=shape,info=ierr)
    if (ierr /= sn_ok) then
      if (present(info)) info = ierr
      return
    end if
    do i=1,size(x)
      x(i) = loc+x(i)/abs(rng%normal())
    end do
    if (present(info)) info = sn_ok
  end subroutine rsc

  real(dp) function zeta(k, x) result(value)
    integer, intent(in) :: k
    real(dp), intent(in) :: x
    real(dp) :: z1, z2, z3, z4, x2, den

    if (k < 0 .or. k > 5) then
      value = ieee_value(0.0_dp,ieee_quiet_nan)
      return
    end if
    if (k == 0) then
      value = normal_logcdf(x)+log(2.0_dp)
      return
    end if
    if (x > -50.0_dp) then
      z1 = exp(normal_logpdf(x)-normal_logcdf(x))
    else
      x2 = x*x
      den = 1.0_dp-1.0_dp/(x2+2.0_dp)+1.0_dp/((x2+2.0_dp)*(x2+4.0_dp)) &
            -5.0_dp/((x2+2.0_dp)*(x2+4.0_dp)*(x2+6.0_dp)) &
            +9.0_dp/((x2+2.0_dp)*(x2+4.0_dp)*(x2+6.0_dp)*(x2+8.0_dp)) &
            -129.0_dp/((x2+2.0_dp)*(x2+4.0_dp)*(x2+6.0_dp)*(x2+8.0_dp)*(x2+10.0_dp))
      z1 = -x/den
    end if
    if (k == 1) then
      value = z1
      return
    end if
    z2 = -z1*(x+z1)
    if (k == 2) then
      value = z2
      return
    end if
    z3 = -z2*(x+z1)-z1*(1.0_dp+z2)
    if (k == 3) then
      value = z3
      return
    end if
    z4 = -z3*(x+2.0_dp*z1)-2.0_dp*z2*(1.0_dp+z2)
    if (k == 4) then
      value = z4
    else
      value = -z4*(x+2.0_dp*z1)-z3*(3.0_dp+4.0_dp*z2)-2.0_dp*z2*z3
    end if
  end function zeta

  subroutine sn_cumulants(params, n, cumulants, info)
    type(sn_uv_params), intent(in) :: params
    integer, intent(in) :: n
    real(dp), allocatable, intent(out) :: cumulants(:)
    integer, intent(out), optional :: info
    real(dp), allocatable :: kv(:), coeff(:), a(:)
    real(dp) :: delta, factorial
    integer :: nn, halfn, i, j

    if (params%omega <= 0.0_dp .or. n < 1) then
      allocate(cumulants(0))
      if (present(info)) info = sn_invalid_argument
      return
    end if
    nn = n
    allocate(cumulants(nn))
    delta = delta_from_alpha(params%alpha)
    if (abs(params%tau) <= tiny_dp) then
      halfn = (max(nn,2)+1)/2
      allocate(a(2*halfn),coeff(2*halfn),kv(2*halfn))
      a = 0.0_dp
      do i=0,halfn-1
        factorial = exp(log_gamma(real(i+1,dp)))
        a(2*i+1) = (-1.0_dp)**i*sqrt(2.0_dp/pi)/(factorial*2.0_dp**i*real(2*i+1,dp))
      end do
      coeff = a(1)
      do i=2,2*halfn
        coeff(i) = a(i)
        do j=1,i-1
          coeff(i) = coeff(i)-real(j,dp)*coeff(j)*a(i-j)/real(i,dp)
        end do
      end do
      do i=1,2*halfn
        kv(i) = coeff(i)*exp(log_gamma(real(i+1,dp)))
      end do
      do i=1,nn
        cumulants(i) = delta**i*kv(i)
      end do
      cumulants(2) = cumulants(2)+1.0_dp
    else
      do i=1,nn
        if (i <= 5) then
          cumulants(i) = zeta(i,params%tau)*delta**i
        else
          cumulants(i) = ieee_value(0.0_dp,ieee_quiet_nan)
        end if
      end do
      if (nn >= 2) cumulants(2) = cumulants(2)+1.0_dp
    end if
    do i=1,nn
      cumulants(i) = cumulants(i)*params%omega**i
    end do
    cumulants(1) = cumulants(1)+params%xi
    if (present(info)) info = sn_ok
  end subroutine sn_cumulants

  pure real(dp) function b_nu(nu) result(value)
    real(dp), intent(in) :: nu
    if (nu <= 1.0_dp) then
      value = ieee_value(0.0_dp,ieee_quiet_nan)
    else if (nu > 1.0e4_dp) then
      value = sqrt(2.0_dp/pi)*(1.0_dp+0.75_dp/nu+0.78125_dp/(nu*nu))
    else
      value = sqrt(nu/pi)*exp(log_gamma(0.5_dp*(nu-1.0_dp))-log_gamma(0.5_dp*nu))
    end if
  end function b_nu

  subroutine st_cumulants(params, n, cumulants, info)
    type(st_uv_params), intent(in) :: params
    integer, intent(in) :: n
    real(dp), allocatable, intent(out) :: cumulants(:)
    integer, intent(out), optional :: info
    real(dp) :: delta, mu, s2
    integer :: nn, i

    if (params%omega <= 0.0_dp .or. params%nu <= 0.0_dp .or. n < 1) then
      allocate(cumulants(0))
      if (present(info)) info = sn_invalid_argument
      return
    end if
    nn = min(n,4)
    allocate(cumulants(nn))
    cumulants = ieee_value(0.0_dp,ieee_quiet_nan)
    if (params%nu > 1.0e8_dp) then
      sn_block: block
        type(sn_uv_params) :: p
        p = sn_uv_params(params%xi,params%omega,params%alpha,0.0_dp)
        call sn_cumulants(p,nn,cumulants)
      end block sn_block
      if (present(info)) info = sn_ok
      return
    end if
    delta = delta_from_alpha(params%alpha)
    mu = b_nu(params%nu)*delta
    cumulants(1) = mu
    if (nn >= 2 .and. params%nu > 2.0_dp) cumulants(2) = params%nu/(params%nu-2.0_dp)-mu*mu
    if (nn >= 3 .and. params%nu > 3.0_dp) then
      cumulants(3) = mu*((3.0_dp-delta*delta)*params%nu/(params%nu-3.0_dp) &
                     -3.0_dp*params%nu/(params%nu-2.0_dp)+2.0_dp*mu*mu)
    end if
    if (nn >= 4 .and. params%nu > 4.0_dp) then
      s2 = params%nu/(params%nu-2.0_dp)-mu*mu
      cumulants(4) = 3.0_dp*params%nu/(params%nu-2.0_dp)*params%nu/(params%nu-4.0_dp) &
                     -4.0_dp*mu*mu*(3.0_dp-delta*delta)*params%nu/(params%nu-3.0_dp) &
                     +6.0_dp*mu*mu*params%nu/(params%nu-2.0_dp)-3.0_dp*mu**4-3.0_dp*s2*s2
    end if
    do i=1,nn
      cumulants(i) = cumulants(i)*params%omega**i
    end do
    cumulants(1) = cumulants(1)+params%xi
    if (present(info)) info = sn_ok
  end subroutine st_cumulants

  subroutine dp_to_cp_sn(params, mean, sd, skewness, info)
    type(sn_uv_params), intent(in) :: params
    real(dp), intent(out) :: mean, sd, skewness
    integer, intent(out), optional :: info
    real(dp), allocatable :: cum(:)
    if (params%omega <= 0.0_dp) then
      mean = 0.0_dp
      sd = 0.0_dp
      skewness = 0.0_dp
      if (present(info)) info = sn_invalid_argument
      return
    end if
    call sn_cumulants(params,3,cum)
    mean = cum(1)
    sd = sqrt(cum(2))
    skewness = cum(3)/sd**3
    if (present(info)) info = sn_ok
  end subroutine dp_to_cp_sn

  subroutine cp_to_dp_sn(mean, sd, skewness, params, info)
    real(dp), intent(in) :: mean, sd, skewness
    type(sn_uv_params), intent(out) :: params
    integer, intent(out), optional :: info
    real(dp) :: max_skew, r, delta, b, muz, sdz

    max_skew = 0.5_dp*(4.0_dp-pi)*(2.0_dp/(pi-2.0_dp))**1.5_dp
    if (sd <= 0.0_dp .or. abs(skewness) >= max_skew) then
      params = sn_uv_params()
      if (present(info)) info = sn_invalid_argument
      return
    end if
    b = sqrt(2.0_dp/pi)
    r = sign(1.0_dp,skewness)*(2.0_dp*abs(skewness)/(4.0_dp-pi))**(1.0_dp/3.0_dp)
    delta = r/(b*sqrt(1.0_dp+r*r))
    muz = b*delta
    sdz = sqrt(1.0_dp-muz*muz)
    params%omega = sd/sdz
    params%alpha = alpha_from_delta(delta)
    params%xi = mean-params%omega*muz
    params%tau = 0.0_dp
    if (present(info)) info = sn_ok
  end subroutine cp_to_dp_sn

  real(dp) function mode_sn(params) result(mode)
    type(sn_uv_params), intent(in) :: params
    mode = maximize_density(1)
  contains
    real(dp) function objective(x) result(v)
      real(dp), intent(in) :: x
      v = dsn(x,params%xi,params%omega,params%alpha,params%tau,log_pdf=.true.)
    end function objective
    real(dp) function maximize_density(dummy) result(xmax)
      integer, intent(in) :: dummy
      real(dp) :: a,b,c,d,fc,fd,gr
      integer :: iter
      if (dummy < 0) continue
      a = params%xi-10.0_dp*params%omega
      b = params%xi+10.0_dp*params%omega
      gr = 0.5_dp*(sqrt(5.0_dp)-1.0_dp)
      c = b-gr*(b-a)
      d = a+gr*(b-a)
      fc = objective(c)
      fd = objective(d)
      do iter=1,160
        if (fc > fd) then
          b=d; d=c; fd=fc; c=b-gr*(b-a); fc=objective(c)
        else
          a=c; c=d; fc=fd; d=a+gr*(b-a); fd=objective(d)
        end if
      end do
      xmax = 0.5_dp*(a+b)
    end function maximize_density
  end function mode_sn

  real(dp) function mode_st(params) result(mode)
    type(st_uv_params), intent(in) :: params
    real(dp) :: a,b,c,d,fc,fd,gr
    integer :: iter
    a = params%xi-20.0_dp*params%omega
    b = params%xi+20.0_dp*params%omega
    gr = 0.5_dp*(sqrt(5.0_dp)-1.0_dp)
    c = b-gr*(b-a); d = a+gr*(b-a)
    fc = dst(c,params%xi,params%omega,params%alpha,params%nu,log_pdf=.true.)
    fd = dst(d,params%xi,params%omega,params%alpha,params%nu,log_pdf=.true.)
    do iter=1,160
      if (fc > fd) then
        b=d; d=c; fd=fc; c=b-gr*(b-a)
        fc=dst(c,params%xi,params%omega,params%alpha,params%nu,log_pdf=.true.)
      else
        a=c; c=d; fc=fd; d=a+gr*(b-a)
        fd=dst(d,params%xi,params%omega,params%alpha,params%nu,log_pdf=.true.)
      end if
    end do
    mode = 0.5_dp*(a+b)
  end function mode_st

  real(dp) function mode_sc(xi,omega,alpha) result(mode)
    real(dp), intent(in), optional :: xi,omega,alpha
    type(st_uv_params) :: p
    p%xi = 0.0_dp; p%omega=1.0_dp; p%alpha=0.0_dp; p%nu=1.0_dp
    if (present(xi)) p%xi=xi
    if (present(omega)) p%omega=omega
    if (present(alpha)) p%alpha=alpha
    mode = mode_st(p)
  end function mode_sc

  real(dp) function quantile_type7(x, p) result(q)
    real(dp), intent(in) :: x(:), p
    real(dp), allocatable :: y(:)
    real(dp) :: h, frac
    integer :: n, j
    n = size(x)
    if (n == 0 .or. p < 0.0_dp .or. p > 1.0_dp) then
      q = ieee_value(0.0_dp,ieee_quiet_nan)
      return
    end if
    allocate(y(n))
    y = x
    call sort_real(y)
    if (n == 1) then
      q = y(1)
      return
    end if
    h = 1.0_dp+real(n-1,dp)*p
    j = floor(h)
    frac = h-real(j,dp)
    if (j >= n) then
      q = y(n)
    else
      q = (1.0_dp-frac)*y(j)+frac*y(j+1)
    end if
  end function quantile_type7

  function fournum(x) result(out)
    real(dp), intent(in) :: x(:)
    type(fournum_result) :: out
    real(dp) :: q(7), den
    integer :: i
    if (size(x) < 8) then
      out%status = sn_invalid_argument
      return
    end if
    do i=1,7
      q(i) = quantile_type7(x,real(i,dp)/8.0_dp)
    end do
    den = q(6)-q(2)
    if (abs(den) <= tiny_dp) then
      out%status = sn_invalid_argument
      return
    end if
    out%median = q(4)
    out%q_deviation = 0.5_dp*den
    out%galton_bowley = (q(6)-2.0_dp*q(4)+q(2))/den
    out%moors = (q(7)-q(5)+q(3)-q(1))/den
    out%status = sn_ok
  end function fournum

  subroutine sort_real(x)
    real(dp), intent(inout) :: x(:)
    integer :: i,j
    real(dp) :: key
    do i=2,size(x)
      key=x(i)
      j=i-1
      do while (j>=1)
        if (x(j)<=key) exit
        x(j+1)=x(j)
        j=j-1
      end do
      x(j+1)=key
    end do
  end subroutine sort_real

end module sn_univariate
