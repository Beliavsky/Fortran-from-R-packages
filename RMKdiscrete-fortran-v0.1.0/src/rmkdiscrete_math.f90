! SPDX-License-Identifier: GPL-2.0-or-later
! Modern Fortran translation of RMKdiscrete 0.1 by Robert M. Kirkpatrick.
! See LICENSE, COPYING, and TRANSLATION_NOTES.md for provenance.
module rmkdiscrete_math
  use, intrinsic :: ieee_arithmetic
  use rmkdiscrete_kinds, only : dp
  implicit none
  private
  public :: qnan, pinf, ninf, real_equal, logsumexp2, kahan_add, log_choose
  public :: dpois, dnbinom_prob, dbinom_prob, rnorm_std, rpois, rgamma_mt, rnbinom_prob, rbinom_simple

contains

  pure real(dp) function qnan() result(x)
    x = ieee_value(0.0_dp, ieee_quiet_nan)
  end function qnan

  pure real(dp) function pinf() result(x)
    x = ieee_value(0.0_dp, ieee_positive_inf)
  end function pinf

  pure real(dp) function ninf() result(x)
    x = ieee_value(0.0_dp, ieee_negative_inf)
  end function ninf

  pure logical function real_equal(a,b) result(eq)
    real(dp), intent(in) :: a,b
    eq=(a<=b .and. a>=b)
  end function real_equal

  pure real(dp) function logsumexp2(a, b) result(v)
    real(dp), intent(in) :: a, b
    real(dp) :: m
    if (real_equal(a,ninf())) then
      v = b
    else if (real_equal(b,ninf())) then
      v = a
    else
      m = max(a, b)
      v = m + log(exp(a-m) + exp(b-m))
    end if
  end function logsumexp2

  pure subroutine kahan_add(sumv, corr, x)
    real(dp), intent(inout) :: sumv, corr
    real(dp), intent(in) :: x
    real(dp) :: y, t
    y = x - corr
    t = sumv + y
    corr = (t - sumv) - y
    sumv = t
  end subroutine kahan_add

  pure real(dp) function log_choose(n, k) result(v)
    integer, intent(in) :: n, k
    if (k < 0 .or. k > n .or. n < 0) then
      v = ninf()
    else
      v = log_gamma(real(n+1,dp)) - log_gamma(real(k+1,dp)) - log_gamma(real(n-k+1,dp))
    end if
  end function log_choose

  pure real(dp) function dpois(x, mu, give_log) result(v)
    integer, intent(in) :: x
    real(dp), intent(in) :: mu
    logical, intent(in), optional :: give_log
    logical :: gl
    real(dp) :: lv
    gl = .false.
    if (present(give_log)) gl = give_log
    if (mu < 0.0_dp .or. x < 0) then
      v = merge(ninf(), 0.0_dp, gl)
      return
    end if
    if (real_equal(mu,0.0_dp)) then
      if (x == 0) then
        v = merge(0.0_dp, 1.0_dp, gl)
      else
        v = merge(ninf(), 0.0_dp, gl)
      end if
      return
    end if
    lv = -mu + real(x,dp)*log(mu) - log_gamma(real(x+1,dp))
    v = merge(lv, exp(lv), gl)
  end function dpois

  pure real(dp) function dnbinom_prob(x, nu, p, give_log) result(v)
    integer, intent(in) :: x
    real(dp), intent(in) :: nu, p
    logical, intent(in), optional :: give_log
    logical :: gl
    real(dp) :: lv
    gl = .false.
    if (present(give_log)) gl = give_log
    if (x < 0 .or. nu < 0.0_dp .or. p <= 0.0_dp .or. p > 1.0_dp) then
      v = merge(ninf(), 0.0_dp, gl)
      return
    end if
    if (real_equal(nu,0.0_dp)) then
      if (x == 0) then
        v = merge(0.0_dp, 1.0_dp, gl)
      else
        v = merge(ninf(), 0.0_dp, gl)
      end if
      return
    end if
    if (real_equal(p,1.0_dp)) then
      if (x == 0) then
        v = merge(0.0_dp, 1.0_dp, gl)
      else
        v = merge(ninf(), 0.0_dp, gl)
      end if
      return
    end if
    lv = log_gamma(real(x,dp)+nu) - log_gamma(nu) - log_gamma(real(x+1,dp)) &
         + nu*log(p) + real(x,dp)*log(1.0_dp-p)
    v = merge(lv, exp(lv), gl)
  end function dnbinom_prob

  pure real(dp) function dbinom_prob(x, n, p, give_log) result(v)
    integer, intent(in) :: x, n
    real(dp), intent(in) :: p
    logical, intent(in), optional :: give_log
    logical :: gl
    real(dp) :: lv
    gl = .false.
    if (present(give_log)) gl = give_log
    if (n < 0 .or. x < 0 .or. x > n .or. p < 0.0_dp .or. p > 1.0_dp) then
      v = merge(ninf(), 0.0_dp, gl)
      return
    end if
    if (real_equal(p,0.0_dp)) then
      if (x == 0) then
        v = merge(0.0_dp, 1.0_dp, gl)
      else
        v = merge(ninf(), 0.0_dp, gl)
      end if
      return
    end if
    if (real_equal(p,1.0_dp)) then
      if (x == n) then
        v = merge(0.0_dp, 1.0_dp, gl)
      else
        v = merge(ninf(), 0.0_dp, gl)
      end if
      return
    end if
    lv = log_choose(n,x) + real(x,dp)*log(p) + real(n-x,dp)*log(1.0_dp-p)
    v = merge(lv, exp(lv), gl)
  end function dbinom_prob

  real(dp) function rnorm_std() result(z)
    real(dp) :: u1, u2
    call random_number(u1)
    call random_number(u2)
    u1 = max(u1, tiny(1.0_dp))
    z = sqrt(-2.0_dp*log(u1))*cos(2.0_dp*acos(-1.0_dp)*u2)
  end function rnorm_std

  integer function rpois(mu) result(k)
    real(dp), intent(in) :: mu
    real(dp) :: l, p, u, sq, alxm, g, y, em, t
    if (mu <= 0.0_dp) then
      k = 0
      return
    end if
    if (mu < 12.0_dp) then
      l = exp(-mu)
      p = 1.0_dp
      k = -1
      do
        k = k + 1
        call random_number(u)
        p = p*u
        if (p <= l) exit
      end do
    else
      sq = sqrt(2.0_dp*mu)
      alxm = log(mu)
      g = mu*alxm - log_gamma(mu+1.0_dp)
      do
        do
          call random_number(u)
          y = tan(acos(-1.0_dp)*u)
          em = sq*y + mu
          if (em >= 0.0_dp) exit
        end do
        em = floor(em)
        t = 0.9_dp*(1.0_dp+y*y)*exp(em*alxm-log_gamma(em+1.0_dp)-g)
        call random_number(u)
        if (u <= t) exit
      end do
      k = int(em)
    end if
  end function rpois

  recursive real(dp) function rgamma_mt(shape, scale) result(x)
    real(dp), intent(in) :: shape, scale
    real(dp) :: d, c, z, v, u
    if (shape <= 0.0_dp .or. scale < 0.0_dp) then
      x = qnan()
      return
    end if
    if (real_equal(scale,0.0_dp)) then
      x = 0.0_dp
      return
    end if
    if (shape < 1.0_dp) then
      call random_number(u)
      x = rgamma_mt(shape+1.0_dp, scale)*u**(1.0_dp/shape)
      return
    end if
    d = shape - 1.0_dp/3.0_dp
    c = 1.0_dp/sqrt(9.0_dp*d)
    do
      do
        z = rnorm_std()
        v = 1.0_dp + c*z
        if (v > 0.0_dp) exit
      end do
      v = v*v*v
      call random_number(u)
      if (u < 1.0_dp - 0.0331_dp*z**4) exit
      if (log(u) < 0.5_dp*z*z + d*(1.0_dp-v+log(v))) exit
    end do
    x = scale*d*v
  end function rgamma_mt

  integer function rnbinom_prob(nu, p) result(x)
    real(dp), intent(in) :: nu, p
    real(dp) :: lam
    if (nu < 0.0_dp .or. p <= 0.0_dp .or. p > 1.0_dp) then
      x = -huge(1)
      return
    end if
    if (real_equal(nu,0.0_dp) .or. real_equal(p,1.0_dp)) then
      x = 0
      return
    end if
    lam = rgamma_mt(nu, (1.0_dp-p)/p)
    x = rpois(lam)
  end function rnbinom_prob

  integer function rbinom_simple(n, p) result(x)
    integer, intent(in) :: n
    real(dp), intent(in) :: p
    integer :: i
    real(dp) :: u
    if (n < 0 .or. p < 0.0_dp .or. p > 1.0_dp) then
      x = -huge(1)
      return
    end if
    x = 0
    do i = 1, n
      call random_number(u)
      if (u <= p) x = x + 1
    end do
  end function rbinom_simple

end module rmkdiscrete_math
