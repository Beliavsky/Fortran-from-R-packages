! SPDX-License-Identifier: GPL-2.0-or-later
module actuar_continuous
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use actuar_kinds, only : dp, pi, huge_dp
  use actuar_special, only : nan_dp, normal_cdf, regularized_gamma_p, &
    regularized_gamma_q, gamma_quantile, regularized_beta, beta_quantile, &
    log_beta_fn
  use actuar_rng, only : runif, rexp, rgamma_shape, rbeta_ab, rinvgauss_rng
  implicit none
  private

  public :: dpareto, ppareto, qpareto, rpareto, mpareto, levpareto
  public :: dpareto1, ppareto1, qpareto1, rpareto1, mpareto1, levpareto1
  public :: dburr, pburr, qburr, rburr, mburr, levburr
  public :: dgenpareto, pgenpareto, qgenpareto, rgenpareto
  public :: mgenpareto, levgenpareto
  public :: dllogis, pllogis, qllogis, rllogis, mllogis, levllogis
  public :: dinvexp, pinvexp, qinvexp, rinvexp, minvexp, levinvexp
  public :: dinvgamma, pinvgamma, qinvgamma, rinvgamma
  public :: minvgamma, levinvgamma
  public :: dinvweibull, pinvweibull, qinvweibull, rinvweibull
  public :: minvweibull, levinvweibull
  public :: dtrgamma, ptrgamma, qtrgamma, rtrgamma, mtrgamma, levtrgamma
  public :: dgenbeta, pgenbeta, qgenbeta, rgenbeta, mgenbeta, levgenbeta
  public :: dgumbel, pgumbel, qgumbel, rgumbel, mgumbel, mgfgumbel
  public :: dinvgauss, pinvgauss, qinvgauss, rinvgauss
  public :: minvgauss, levinvgauss, mgfinvgauss

contains

  pure function dpareto(x, shape, scale) result(y)
    real(dp), intent(in) :: x, shape, scale
    real(dp) :: y
    if (shape <= 0.0_dp .or. scale <= 0.0_dp) then
      y = nan_dp()
    else if (x < 0.0_dp .or. .not. ieee_is_finite(x)) then
      y = 0.0_dp
    else
      y = shape/scale*(1.0_dp+x/scale)**(-shape-1.0_dp)
    end if
  end function dpareto

  pure function ppareto(x, shape, scale) result(y)
    real(dp), intent(in) :: x, shape, scale
    real(dp) :: y
    if (shape <= 0.0_dp .or. scale <= 0.0_dp) then
      y = nan_dp()
    else if (x <= 0.0_dp) then
      y = 0.0_dp
    else
      y = 1.0_dp-(1.0_dp+x/scale)**(-shape)
    end if
  end function ppareto

  pure function qpareto(p, shape, scale) result(x)
    real(dp), intent(in) :: p, shape, scale
    real(dp) :: x
    if (shape <= 0.0_dp .or. scale <= 0.0_dp .or. p < 0.0_dp .or. p > 1.0_dp) then
      x = nan_dp()
    else if (p == 1.0_dp) then
      x = huge_dp
    else
      x = scale*((1.0_dp-p)**(-1.0_dp/shape)-1.0_dp)
    end if
  end function qpareto

  function rpareto(shape, scale) result(x)
    real(dp), intent(in) :: shape, scale
    real(dp) :: x
    x = qpareto(runif(),shape,scale)
  end function rpareto

  pure function mpareto(order, shape, scale) result(m)
    real(dp), intent(in) :: order, shape, scale
    real(dp) :: m
    if (shape <= 0.0_dp .or. scale <= 0.0_dp) then
      m = nan_dp()
    else if (order <= -1.0_dp .or. order >= shape) then
      m = huge_dp
    else
      m = scale**order*gamma(1.0_dp+order)*gamma(shape-order)/gamma(shape)
    end if
  end function mpareto

  pure function levpareto(limit, shape, scale, order) result(m)
    real(dp), intent(in) :: limit, shape, scale, order
    real(dp) :: m, z, p
    if (shape <= 0.0_dp .or. scale <= 0.0_dp .or. order <= -1.0_dp) then
      m = merge(huge_dp,nan_dp(),order <= -1.0_dp)
    else if (limit <= 0.0_dp) then
      m = 0.0_dp
    else
      z = limit/(limit+scale)
      p = regularized_beta(z,1.0_dp+order,shape-order)
      m = scale**order*gamma(1.0_dp+order)*gamma(shape-order)/gamma(shape)*p + &
          limit**order*(1.0_dp-ppareto(limit,shape,scale))
    end if
  end function levpareto

  pure function dpareto1(x, shape, xmin) result(y)
    real(dp), intent(in) :: x, shape, xmin
    real(dp) :: y
    if (shape <= 0.0_dp .or. xmin <= 0.0_dp) then
      y = nan_dp()
    else if (x < xmin .or. .not. ieee_is_finite(x)) then
      y = 0.0_dp
    else
      y = shape*xmin**shape/x**(shape+1.0_dp)
    end if
  end function dpareto1

  pure function ppareto1(x, shape, xmin) result(y)
    real(dp), intent(in) :: x, shape, xmin
    real(dp) :: y
    if (shape <= 0.0_dp .or. xmin <= 0.0_dp) then
      y = nan_dp()
    else if (x <= xmin) then
      y = 0.0_dp
    else
      y = 1.0_dp-(xmin/x)**shape
    end if
  end function ppareto1

  pure function qpareto1(p, shape, xmin) result(x)
    real(dp), intent(in) :: p, shape, xmin
    real(dp) :: x
    if (shape <= 0.0_dp .or. xmin <= 0.0_dp .or. p < 0.0_dp .or. p > 1.0_dp) then
      x = nan_dp()
    else if (p == 1.0_dp) then
      x = huge_dp
    else
      x = xmin/(1.0_dp-p)**(1.0_dp/shape)
    end if
  end function qpareto1

  function rpareto1(shape, xmin) result(x)
    real(dp), intent(in) :: shape, xmin
    real(dp) :: x
    x = qpareto1(runif(),shape,xmin)
  end function rpareto1

  pure function mpareto1(order, shape, xmin) result(m)
    real(dp), intent(in) :: order, shape, xmin
    real(dp) :: m
    if (shape <= 0.0_dp .or. xmin <= 0.0_dp) then
      m = nan_dp()
    else if (order >= shape) then
      m = huge_dp
    else
      m = shape*xmin**order/(shape-order)
    end if
  end function mpareto1

  pure function levpareto1(limit, shape, xmin, order) result(m)
    real(dp), intent(in) :: limit, shape, xmin, order
    real(dp) :: m, tmp
    if (shape <= 0.0_dp .or. xmin <= 0.0_dp) then
      m = nan_dp()
    else if (limit <= xmin) then
      m = 0.0_dp
    else if (abs(shape-order) < 1.0e-14_dp) then
      m = shape*xmin**shape*log(limit/xmin) + limit**order*(xmin/limit)**shape
    else
      tmp = shape-order
      m = shape*xmin**order/tmp - order*xmin**shape/(tmp*limit**tmp)
    end if
  end function levpareto1

  pure function dburr(x, shape1, shape2, scale) result(y)
    real(dp), intent(in) :: x, shape1, shape2, scale
    real(dp) :: y, v
    if (min(shape1,shape2,scale) <= 0.0_dp) then
      y = nan_dp()
    else if (x < 0.0_dp .or. .not. ieee_is_finite(x)) then
      y = 0.0_dp
    else if (x == 0.0_dp) then
      if (shape2 < 1.0_dp) then
        y = huge_dp
      else if (shape2 > 1.0_dp) then
        y = 0.0_dp
      else
        y = shape1/scale
      end if
    else
      v = (x/scale)**shape2
      y = shape1*shape2*v/(x*(1.0_dp+v)**(shape1+1.0_dp))
    end if
  end function dburr

  pure function pburr(x, shape1, shape2, scale) result(y)
    real(dp), intent(in) :: x, shape1, shape2, scale
    real(dp) :: y
    if (min(shape1,shape2,scale) <= 0.0_dp) then
      y = nan_dp()
    else if (x <= 0.0_dp) then
      y = 0.0_dp
    else
      y = 1.0_dp-(1.0_dp+(x/scale)**shape2)**(-shape1)
    end if
  end function pburr

  pure function qburr(p, shape1, shape2, scale) result(x)
    real(dp), intent(in) :: p, shape1, shape2, scale
    real(dp) :: x
    if (min(shape1,shape2,scale) <= 0.0_dp .or. p < 0.0_dp .or. p > 1.0_dp) then
      x = nan_dp()
    else if (p == 1.0_dp) then
      x = huge_dp
    else
      x = scale*((1.0_dp-p)**(-1.0_dp/shape1)-1.0_dp)**(1.0_dp/shape2)
    end if
  end function qburr

  function rburr(shape1, shape2, scale) result(x)
    real(dp), intent(in) :: shape1, shape2, scale
    real(dp) :: x
    x = qburr(runif(),shape1,shape2,scale)
  end function rburr

  pure function mburr(order, shape1, shape2, scale) result(m)
    real(dp), intent(in) :: order, shape1, shape2, scale
    real(dp) :: m, t
    if (min(shape1,shape2,scale) <= 0.0_dp) then
      m = nan_dp()
    else if (order <= -shape2 .or. order >= shape1*shape2) then
      m = huge_dp
    else
      t = order/shape2
      m = scale**order*gamma(1.0_dp+t)*gamma(shape1-t)/gamma(shape1)
    end if
  end function mburr

  pure function levburr(limit, shape1, shape2, scale, order) result(m)
    real(dp), intent(in) :: limit, shape1, shape2, scale, order
    real(dp) :: m, z, t, full
    if (min(shape1,shape2,scale) <= 0.0_dp .or. order <= -shape2) then
      m = merge(huge_dp,nan_dp(),order <= -shape2)
    else if (limit <= 0.0_dp) then
      m = 0.0_dp
    else
      t = order/shape2
      z = (limit/scale)**shape2/(1.0_dp+(limit/scale)**shape2)
      full = scale**order*gamma(1.0_dp+t)*gamma(shape1-t)/gamma(shape1)
      m = full*regularized_beta(z,1.0_dp+t,shape1-t) + &
          limit**order*(1.0_dp-pburr(limit,shape1,shape2,scale))
    end if
  end function levburr

  pure function dgenpareto(x, shape1, shape2, scale) result(y)
    real(dp), intent(in) :: x, shape1, shape2, scale
    real(dp) :: y, u
    if (min(shape1,shape2,scale) <= 0.0_dp) then
      y = nan_dp()
    else if (x < 0.0_dp .or. .not. ieee_is_finite(x)) then
      y = 0.0_dp
    else if (x == 0.0_dp) then
      if (shape2 < 1.0_dp) then
        y = huge_dp
      else if (shape2 > 1.0_dp) then
        y = 0.0_dp
      else
        y = exp(-log(scale)-log_beta_fn(shape2,shape1))
      end if
    else
      u = x/(x+scale)
      y = u**shape2*(1.0_dp-u)**shape1/(x*exp(log_beta_fn(shape2,shape1)))
    end if
  end function dgenpareto

  pure function pgenpareto(x, shape1, shape2, scale) result(y)
    real(dp), intent(in) :: x, shape1, shape2, scale
    real(dp) :: y
    if (min(shape1,shape2,scale) <= 0.0_dp) then
      y = nan_dp()
    else if (x <= 0.0_dp) then
      y = 0.0_dp
    else
      y = regularized_beta(x/(x+scale),shape2,shape1)
    end if
  end function pgenpareto

  function qgenpareto(p, shape1, shape2, scale) result(x)
    real(dp), intent(in) :: p, shape1, shape2, scale
    real(dp) :: x, u
    if (min(shape1,shape2,scale) <= 0.0_dp .or. p < 0.0_dp .or. p > 1.0_dp) then
      x = nan_dp()
    else if (p == 1.0_dp) then
      x = huge_dp
    else
      u = beta_quantile(p,shape2,shape1)
      x = scale*u/(1.0_dp-u)
    end if
  end function qgenpareto

  function rgenpareto(shape1, shape2, scale) result(x)
    real(dp), intent(in) :: shape1, shape2, scale
    real(dp) :: x, u
    u = rbeta_ab(shape2,shape1)
    x = scale*u/(1.0_dp-u)
  end function rgenpareto

  pure function mgenpareto(order, shape1, shape2, scale) result(m)
    real(dp), intent(in) :: order, shape1, shape2, scale
    real(dp) :: m
    if (min(shape1,shape2,scale) <= 0.0_dp) then
      m = nan_dp()
    else if (order <= -shape2 .or. order >= shape1) then
      m = huge_dp
    else
      m = scale**order*exp(log_beta_fn(shape1-order,shape2+order)- &
          log_beta_fn(shape1,shape2))
    end if
  end function mgenpareto

  pure function levgenpareto(limit, shape1, shape2, scale, order) result(m)
    real(dp), intent(in) :: limit, shape1, shape2, scale, order
    real(dp) :: m, z, full
    if (min(shape1,shape2,scale) <= 0.0_dp .or. order <= -shape2) then
      m = merge(huge_dp,nan_dp(),order <= -shape2)
    else if (limit <= 0.0_dp) then
      m = 0.0_dp
    else
      z = limit/(limit+scale)
      full = scale**order*exp(log_beta_fn(shape1-order,shape2+order)- &
        log_beta_fn(shape1,shape2))
      m = full*regularized_beta(z,shape2+order,shape1-order) + &
          limit**order*(1.0_dp-pgenpareto(limit,shape1,shape2,scale))
    end if
  end function levgenpareto

  pure function dllogis(x, shape, scale) result(y)
    real(dp), intent(in) :: x, shape, scale
    real(dp) :: y, v
    if (shape <= 0.0_dp .or. scale <= 0.0_dp) then
      y = nan_dp()
    else if (x < 0.0_dp .or. .not. ieee_is_finite(x)) then
      y = 0.0_dp
    else if (x == 0.0_dp) then
      if (shape < 1.0_dp) then
        y = huge_dp
      else if (shape > 1.0_dp) then
        y = 0.0_dp
      else
        y = 1.0_dp/scale
      end if
    else
      v = (x/scale)**shape
      y = shape*v/(x*(1.0_dp+v)**2)
    end if
  end function dllogis

  pure function pllogis(x, shape, scale) result(y)
    real(dp), intent(in) :: x, shape, scale
    real(dp) :: y, v
    if (shape <= 0.0_dp .or. scale <= 0.0_dp) then
      y = nan_dp()
    else if (x <= 0.0_dp) then
      y = 0.0_dp
    else
      v = (x/scale)**shape
      y = v/(1.0_dp+v)
    end if
  end function pllogis

  pure function qllogis(p, shape, scale) result(x)
    real(dp), intent(in) :: p, shape, scale
    real(dp) :: x
    if (shape <= 0.0_dp .or. scale <= 0.0_dp .or. p < 0.0_dp .or. p > 1.0_dp) then
      x = nan_dp()
    else if (p == 1.0_dp) then
      x = huge_dp
    else
      x = scale*(p/(1.0_dp-p))**(1.0_dp/shape)
    end if
  end function qllogis

  function rllogis(shape, scale) result(x)
    real(dp), intent(in) :: shape, scale
    real(dp) :: x
    x = qllogis(runif(),shape,scale)
  end function rllogis

  pure function mllogis(order, shape, scale) result(m)
    real(dp), intent(in) :: order, shape, scale
    real(dp) :: m, t
    if (shape <= 0.0_dp .or. scale <= 0.0_dp) then
      m = nan_dp()
    else if (abs(order) >= shape) then
      m = huge_dp
    else
      t = order/shape
      m = scale**order*gamma(1.0_dp+t)*gamma(1.0_dp-t)
    end if
  end function mllogis

  pure function levllogis(limit, shape, scale, order) result(m)
    real(dp), intent(in) :: limit, shape, scale, order
    real(dp) :: m, z, t, full
    if (shape <= 0.0_dp .or. scale <= 0.0_dp .or. order <= -shape) then
      m = merge(huge_dp,nan_dp(),order <= -shape)
    else if (limit <= 0.0_dp) then
      m = 0.0_dp
    else
      z = pllogis(limit,shape,scale)
      t = order/shape
      full = scale**order*gamma(1.0_dp+t)*gamma(1.0_dp-t)
      m = full*regularized_beta(z,1.0_dp+t,1.0_dp-t) + &
          limit**order*(1.0_dp-z)
    end if
  end function levllogis

  pure function dinvexp(x, scale) result(y)
    real(dp), intent(in) :: x, scale
    real(dp) :: y
    if (scale <= 0.0_dp) then
      y = nan_dp()
    else if (x <= 0.0_dp .or. .not. ieee_is_finite(x)) then
      y = 0.0_dp
    else
      y = scale*exp(-scale/x)/(x*x)
    end if
  end function dinvexp

  pure function pinvexp(x, scale) result(y)
    real(dp), intent(in) :: x, scale
    real(dp) :: y
    if (scale <= 0.0_dp) then
      y = nan_dp()
    else if (x <= 0.0_dp) then
      y = 0.0_dp
    else
      y = exp(-scale/x)
    end if
  end function pinvexp

  pure function qinvexp(p, scale) result(x)
    real(dp), intent(in) :: p, scale
    real(dp) :: x
    if (scale <= 0.0_dp .or. p < 0.0_dp .or. p > 1.0_dp) then
      x = nan_dp()
    else if (p == 0.0_dp) then
      x = 0.0_dp
    else if (p == 1.0_dp) then
      x = huge_dp
    else
      x = -scale/log(p)
    end if
  end function qinvexp

  function rinvexp(scale) result(x)
    real(dp), intent(in) :: scale
    real(dp) :: x
    x = scale/rexp()
  end function rinvexp

  pure function minvexp(order, scale) result(m)
    real(dp), intent(in) :: order, scale
    real(dp) :: m
    if (scale <= 0.0_dp) then
      m = nan_dp()
    else if (order >= 1.0_dp) then
      m = huge_dp
    else
      m = scale**order*gamma(1.0_dp-order)
    end if
  end function minvexp

  pure function levinvexp(limit, scale, order) result(m)
    real(dp), intent(in) :: limit, scale, order
    real(dp) :: m, u
    if (scale <= 0.0_dp) then
      m = nan_dp()
    else if (limit <= 0.0_dp) then
      m = 0.0_dp
    else
      u = scale/limit
      m = scale**order*gamma(1.0_dp-order)*regularized_gamma_q(1.0_dp-order,u) + &
          limit**order*(1.0_dp-pinvexp(limit,scale))
    end if
  end function levinvexp

  pure function dinvgamma(x, shape, scale) result(y)
    real(dp), intent(in) :: x, shape, scale
    real(dp) :: y
    if (shape <= 0.0_dp .or. scale <= 0.0_dp) then
      y = nan_dp()
    else if (x <= 0.0_dp .or. .not. ieee_is_finite(x)) then
      y = 0.0_dp
    else
      y = exp(shape*log(scale)-log_gamma(shape)-(shape+1.0_dp)*log(x)-scale/x)
    end if
  end function dinvgamma

  pure function pinvgamma(x, shape, scale) result(y)
    real(dp), intent(in) :: x, shape, scale
    real(dp) :: y
    if (shape <= 0.0_dp .or. scale <= 0.0_dp) then
      y = nan_dp()
    else if (x <= 0.0_dp) then
      y = 0.0_dp
    else
      y = regularized_gamma_q(shape,scale/x)
    end if
  end function pinvgamma

  function qinvgamma(p, shape, scale) result(x)
    real(dp), intent(in) :: p, shape, scale
    real(dp) :: x
    if (shape <= 0.0_dp .or. scale <= 0.0_dp .or. p < 0.0_dp .or. p > 1.0_dp) then
      x = nan_dp()
    else if (p == 0.0_dp) then
      x = 0.0_dp
    else if (p == 1.0_dp) then
      x = huge_dp
    else
      x = scale/gamma_quantile(1.0_dp-p,shape)
    end if
  end function qinvgamma

  function rinvgamma(shape, scale) result(x)
    real(dp), intent(in) :: shape, scale
    real(dp) :: x
    x = scale/rgamma_shape(shape)
  end function rinvgamma

  pure function minvgamma(order, shape, scale) result(m)
    real(dp), intent(in) :: order, shape, scale
    real(dp) :: m
    if (shape <= 0.0_dp .or. scale <= 0.0_dp) then
      m = nan_dp()
    else if (order >= shape) then
      m = huge_dp
    else
      m = scale**order*gamma(shape-order)/gamma(shape)
    end if
  end function minvgamma

  pure function levinvgamma(limit, shape, scale, order) result(m)
    real(dp), intent(in) :: limit, shape, scale, order
    real(dp) :: m, u
    if (shape <= 0.0_dp .or. scale <= 0.0_dp) then
      m = nan_dp()
    else if (limit <= 0.0_dp) then
      m = 0.0_dp
    else if (order >= shape) then
      m = huge_dp
    else
      u = scale/limit
      m = scale**order*gamma(shape-order)/gamma(shape)* &
          regularized_gamma_q(shape-order,u) + &
          limit**order*(1.0_dp-pinvgamma(limit,shape,scale))
    end if
  end function levinvgamma

  pure function dinvweibull(x, shape, scale) result(y)
    real(dp), intent(in) :: x, shape, scale
    real(dp) :: y, u
    if (shape <= 0.0_dp .or. scale <= 0.0_dp) then
      y = nan_dp()
    else if (x <= 0.0_dp .or. .not. ieee_is_finite(x)) then
      y = 0.0_dp
    else
      u = (scale/x)**shape
      y = shape*u*exp(-u)/x
    end if
  end function dinvweibull

  pure function pinvweibull(x, shape, scale) result(y)
    real(dp), intent(in) :: x, shape, scale
    real(dp) :: y
    if (shape <= 0.0_dp .or. scale <= 0.0_dp) then
      y = nan_dp()
    else if (x <= 0.0_dp) then
      y = 0.0_dp
    else
      y = exp(-(scale/x)**shape)
    end if
  end function pinvweibull

  pure function qinvweibull(p, shape, scale) result(x)
    real(dp), intent(in) :: p, shape, scale
    real(dp) :: x
    if (shape <= 0.0_dp .or. scale <= 0.0_dp .or. p < 0.0_dp .or. p > 1.0_dp) then
      x = nan_dp()
    else if (p == 0.0_dp) then
      x = 0.0_dp
    else if (p == 1.0_dp) then
      x = huge_dp
    else
      x = scale*(-log(p))**(-1.0_dp/shape)
    end if
  end function qinvweibull

  function rinvweibull(shape, scale) result(x)
    real(dp), intent(in) :: shape, scale
    real(dp) :: x
    x = scale*rexp()**(-1.0_dp/shape)
  end function rinvweibull

  pure function minvweibull(order, shape, scale) result(m)
    real(dp), intent(in) :: order, shape, scale
    real(dp) :: m
    if (shape <= 0.0_dp .or. scale <= 0.0_dp) then
      m = nan_dp()
    else if (order >= shape) then
      m = huge_dp
    else
      m = scale**order*gamma(1.0_dp-order/shape)
    end if
  end function minvweibull

  pure function levinvweibull(limit, shape, scale, order) result(m)
    real(dp), intent(in) :: limit, shape, scale, order
    real(dp) :: m, u
    if (shape <= 0.0_dp .or. scale <= 0.0_dp) then
      m = nan_dp()
    else if (limit <= 0.0_dp) then
      m = 0.0_dp
    else if (order >= shape) then
      m = huge_dp
    else
      u = (scale/limit)**shape
      m = scale**order*gamma(1.0_dp-order/shape)* &
          regularized_gamma_q(1.0_dp-order/shape,u) + &
          limit**order*(1.0_dp-pinvweibull(limit,shape,scale))
    end if
  end function levinvweibull

  pure function dtrgamma(x, shape1, shape2, scale) result(y)
    real(dp), intent(in) :: x, shape1, shape2, scale
    real(dp) :: y, u
    if (min(shape1,shape2,scale) <= 0.0_dp) then
      y = nan_dp()
    else if (x < 0.0_dp .or. .not. ieee_is_finite(x)) then
      y = 0.0_dp
    else if (x == 0.0_dp) then
      if (shape1*shape2 < 1.0_dp) then
        y = huge_dp
      else if (shape1*shape2 > 1.0_dp) then
        y = 0.0_dp
      else
        y = shape2/(scale*gamma(shape1))
      end if
    else
      u = (x/scale)**shape2
      y = shape2*u**shape1*exp(-u)/(x*gamma(shape1))
    end if
  end function dtrgamma

  pure function ptrgamma(x, shape1, shape2, scale) result(y)
    real(dp), intent(in) :: x, shape1, shape2, scale
    real(dp) :: y
    if (min(shape1,shape2,scale) <= 0.0_dp) then
      y = nan_dp()
    else if (x <= 0.0_dp) then
      y = 0.0_dp
    else
      y = regularized_gamma_p(shape1,(x/scale)**shape2)
    end if
  end function ptrgamma

  function qtrgamma(p, shape1, shape2, scale) result(x)
    real(dp), intent(in) :: p, shape1, shape2, scale
    real(dp) :: x
    if (min(shape1,shape2,scale) <= 0.0_dp .or. p < 0.0_dp .or. p > 1.0_dp) then
      x = nan_dp()
    else
      x = scale*gamma_quantile(p,shape1)**(1.0_dp/shape2)
    end if
  end function qtrgamma

  function rtrgamma(shape1, shape2, scale) result(x)
    real(dp), intent(in) :: shape1, shape2, scale
    real(dp) :: x
    x = scale*rgamma_shape(shape1)**(1.0_dp/shape2)
  end function rtrgamma

  pure function mtrgamma(order, shape1, shape2, scale) result(m)
    real(dp), intent(in) :: order, shape1, shape2, scale
    real(dp) :: m
    if (min(shape1,shape2,scale) <= 0.0_dp) then
      m = nan_dp()
    else if (order <= -shape1*shape2) then
      m = huge_dp
    else
      m = scale**order*gamma(shape1+order/shape2)/gamma(shape1)
    end if
  end function mtrgamma

  pure function levtrgamma(limit, shape1, shape2, scale, order) result(m)
    real(dp), intent(in) :: limit, shape1, shape2, scale, order
    real(dp) :: m, u, a
    if (min(shape1,shape2,scale) <= 0.0_dp) then
      m = nan_dp()
    else if (limit <= 0.0_dp) then
      m = 0.0_dp
    else if (order <= -shape1*shape2) then
      m = huge_dp
    else
      u = (limit/scale)**shape2
      a = shape1+order/shape2
      m = scale**order*gamma(a)/gamma(shape1)*regularized_gamma_p(a,u) + &
          limit**order*(1.0_dp-ptrgamma(limit,shape1,shape2,scale))
    end if
  end function levtrgamma

  pure function dgenbeta(x, shape1, shape2, shape3, scale) result(y)
    real(dp), intent(in) :: x, shape1, shape2, shape3, scale
    real(dp) :: y, u
    if (min(shape1,shape2,shape3,scale) <= 0.0_dp) then
      y = nan_dp()
    else if (x < 0.0_dp .or. x > scale) then
      y = 0.0_dp
    else if (x == 0.0_dp .or. x == scale) then
      y = 0.0_dp
    else
      u = (x/scale)**shape3
      y = shape3*u**shape1*(1.0_dp-u)**(shape2-1.0_dp) / &
          (x*exp(log_beta_fn(shape1,shape2)))
    end if
  end function dgenbeta

  pure function pgenbeta(x, shape1, shape2, shape3, scale) result(y)
    real(dp), intent(in) :: x, shape1, shape2, shape3, scale
    real(dp) :: y
    if (min(shape1,shape2,shape3,scale) <= 0.0_dp) then
      y = nan_dp()
    else if (x <= 0.0_dp) then
      y = 0.0_dp
    else if (x >= scale) then
      y = 1.0_dp
    else
      y = regularized_beta((x/scale)**shape3,shape1,shape2)
    end if
  end function pgenbeta

  function qgenbeta(p, shape1, shape2, shape3, scale) result(x)
    real(dp), intent(in) :: p, shape1, shape2, shape3, scale
    real(dp) :: x
    if (min(shape1,shape2,shape3,scale) <= 0.0_dp .or. p < 0.0_dp .or. p > 1.0_dp) then
      x = nan_dp()
    else
      x = scale*beta_quantile(p,shape1,shape2)**(1.0_dp/shape3)
    end if
  end function qgenbeta

  function rgenbeta(shape1, shape2, shape3, scale) result(x)
    real(dp), intent(in) :: shape1, shape2, shape3, scale
    real(dp) :: x
    x = scale*rbeta_ab(shape1,shape2)**(1.0_dp/shape3)
  end function rgenbeta

  pure function mgenbeta(order, shape1, shape2, shape3, scale) result(m)
    real(dp), intent(in) :: order, shape1, shape2, shape3, scale
    real(dp) :: m
    if (min(shape1,shape2,shape3,scale) <= 0.0_dp) then
      m = nan_dp()
    else if (order <= -shape1*shape3) then
      m = huge_dp
    else
      m = scale**order*exp(log_beta_fn(shape1+order/shape3,shape2)- &
          log_beta_fn(shape1,shape2))
    end if
  end function mgenbeta

  pure function levgenbeta(limit, shape1, shape2, shape3, scale, order) result(m)
    real(dp), intent(in) :: limit, shape1, shape2, shape3, scale, order
    real(dp) :: m, z, full
    if (min(shape1,shape2,shape3,scale) <= 0.0_dp) then
      m = nan_dp()
    else if (limit <= 0.0_dp) then
      m = 0.0_dp
    else if (limit >= scale) then
      m = mgenbeta(order,shape1,shape2,shape3,scale)
    else if (order <= -shape1*shape3) then
      m = huge_dp
    else
      z = (limit/scale)**shape3
      full = scale**order*exp(log_beta_fn(shape1+order/shape3,shape2)- &
          log_beta_fn(shape1,shape2))
      m = full*regularized_beta(z,shape1+order/shape3,shape2) + &
          limit**order*(1.0_dp-pgenbeta(limit,shape1,shape2,shape3,scale))
    end if
  end function levgenbeta

  pure function dgumbel(x, location, scale) result(y)
    real(dp), intent(in) :: x, location, scale
    real(dp) :: y, z
    if (scale <= 0.0_dp) then
      y = nan_dp()
    else
      z = (x-location)/scale
      y = exp(-z-exp(-z))/scale
    end if
  end function dgumbel

  pure function pgumbel(x, location, scale) result(y)
    real(dp), intent(in) :: x, location, scale
    real(dp) :: y
    if (scale <= 0.0_dp) then
      y = nan_dp()
    else
      y = exp(-exp(-(x-location)/scale))
    end if
  end function pgumbel

  pure function qgumbel(p, location, scale) result(x)
    real(dp), intent(in) :: p, location, scale
    real(dp) :: x
    if (scale <= 0.0_dp .or. p < 0.0_dp .or. p > 1.0_dp) then
      x = nan_dp()
    else if (p == 0.0_dp) then
      x = -huge_dp
    else if (p == 1.0_dp) then
      x = huge_dp
    else
      x = location-scale*log(-log(p))
    end if
  end function qgumbel

  function rgumbel(location, scale) result(x)
    real(dp), intent(in) :: location, scale
    real(dp) :: x
    x = qgumbel(runif(),location,scale)
  end function rgumbel

  pure function mgumbel(order, location, scale) result(m)
    integer, intent(in) :: order
    real(dp), intent(in) :: location, scale
    real(dp) :: m, euler, zeta3
    euler = 0.5772156649015329_dp
    zeta3 = 1.2020569031595943_dp
    select case(order)
    case(0); m = 1.0_dp
    case(1); m = location+euler*scale
    case(2); m = (location+euler*scale)**2 + pi*pi*scale*scale/6.0_dp
    case(3)
      m = (location+euler*scale)**3 + 0.5_dp*pi*pi*scale*scale* &
          (location+euler*scale) + 2.0_dp*zeta3*scale**3
    case default; m = nan_dp()
    end select
  end function mgumbel

  pure function mgfgumbel(t, location, scale) result(m)
    real(dp), intent(in) :: t, location, scale
    real(dp) :: m
    if (scale <= 0.0_dp .or. t*scale >= 1.0_dp) then
      m = nan_dp()
    else
      m = exp(location*t)*gamma(1.0_dp-scale*t)
    end if
  end function mgfgumbel

  pure function dinvgauss(x, mu, phi) result(y)
    real(dp), intent(in) :: x, mu, phi
    real(dp) :: y
    if (mu <= 0.0_dp .or. phi <= 0.0_dp) then
      y = nan_dp()
    else if (x <= 0.0_dp .or. .not. ieee_is_finite(x)) then
      y = 0.0_dp
    else
      y = exp(-0.5_dp*(log(2.0_dp*pi*phi)+3.0_dp*log(x)+ &
          ((x-mu)/mu)**2/(phi*x)))
    end if
  end function dinvgauss

  pure function pinvgauss(x, mu, phi) result(y)
    real(dp), intent(in) :: x, mu, phi
    real(dp) :: y, xm, phim, r, a, b
    if (mu <= 0.0_dp .or. phi <= 0.0_dp) then
      y = nan_dp()
    else if (x <= 0.0_dp) then
      y = 0.0_dp
    else
      xm = x/mu
      phim = phi*mu
      r = sqrt(x*phi)
      a = normal_cdf((xm-1.0_dp)/r)
      b = exp(min(700.0_dp,2.0_dp/phim))*normal_cdf(-(xm+1.0_dp)/r)
      y = min(1.0_dp,a+b)
    end if
  end function pinvgauss

  function qinvgauss(p, mu, phi) result(x)
    real(dp), intent(in) :: p, mu, phi
    real(dp) :: x, lo, hi, mid
    integer :: iter
    if (mu <= 0.0_dp .or. phi <= 0.0_dp .or. p < 0.0_dp .or. p > 1.0_dp) then
      x = nan_dp(); return
    end if
    if (p == 0.0_dp) then
      x = 0.0_dp; return
    else if (p == 1.0_dp) then
      x = huge_dp; return
    end if
    lo = 0.0_dp; hi = max(mu,1.0_dp)
    do while (pinvgauss(hi,mu,phi) < p)
      hi = 2.0_dp*hi
      if (hi > huge_dp/4.0_dp) exit
    end do
    do iter = 1, 180
      mid = 0.5_dp*(lo+hi)
      if (pinvgauss(mid,mu,phi) < p) then
        lo = mid
      else
        hi = mid
      end if
      if (abs(hi-lo) < 2.0e-13_dp*max(1.0_dp,mid)) exit
    end do
    x = 0.5_dp*(lo+hi)
  end function qinvgauss

  function rinvgauss(mu, phi) result(x)
    real(dp), intent(in) :: mu, phi
    real(dp) :: x
    x = rinvgauss_rng(mu,phi)
  end function rinvgauss

  pure function minvgauss(order, mu, phi) result(m)
    integer, intent(in) :: order
    real(dp), intent(in) :: mu, phi
    real(dp) :: m, term, s, phir
    integer :: i
    if (mu <= 0.0_dp .or. phi <= 0.0_dp .or. order < 0) then
      m = nan_dp(); return
    end if
    if (order == 0) then
      m = 1.0_dp; return
    end if
    s = 1.0_dp; term = 1.0_dp; phir = phi*mu/2.0_dp
    do i = 1, order-1
      term = term*((real(order+i-1,dp)*real(order-i,dp))/real(i,dp))*phir
      s = s+term
    end do
    m = mu**order*s
  end function minvgauss

  pure function levinvgauss(limit, mu, phi) result(m)
    real(dp), intent(in) :: limit, mu, phi
    real(dp) :: m, xm, phim, r, a, b
    if (mu <= 0.0_dp .or. phi <= 0.0_dp) then
      m = nan_dp()
    else if (limit <= 0.0_dp) then
      m = 0.0_dp
    else
      xm = limit/mu; phim = phi*mu; r = sqrt(limit*phi)
      a = normal_cdf((xm-1.0_dp)/r)
      b = exp(min(700.0_dp,2.0_dp/phim))*normal_cdf(-(xm+1.0_dp)/r)
      m = mu*max(0.0_dp,a-b) + limit*(1.0_dp-pinvgauss(limit,mu,phi))
    end if
  end function levinvgauss

  pure function mgfinvgauss(t, mu, phi) result(m)
    real(dp), intent(in) :: t, mu, phi
    real(dp) :: m
    if (mu <= 0.0_dp .or. phi <= 0.0_dp .or. t > 1.0_dp/(2.0_dp*phi*mu*mu)) then
      m = nan_dp()
    else
      m = exp((1.0_dp-sqrt(1.0_dp-2.0_dp*phi*mu*mu*t))/(phi*mu))
    end if
  end function mgfinvgauss

end module actuar_continuous
