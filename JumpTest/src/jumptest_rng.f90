! SPDX-License-Identifier: MIT
module jumptest_rng
  use jumptest_kinds, only : dp, i8, pi
  implicit none
  private

  integer(i8), parameter :: pm_m = 2147483647_i8
  integer(i8), parameter :: pm_a = 16807_i8
  integer(i8), parameter :: pm_q = 127773_i8
  integer(i8), parameter :: pm_r = 2836_i8

  type, public :: rng_state
    integer(i8) :: state = 104729_i8
    logical :: has_spare = .false.
    real(dp) :: spare = 0.0_dp
  end type rng_state

  public :: rng_seed, rng_uniform, rng_normal, rng_gamma, rng_chisq
  public :: rng_poisson

contains

  subroutine rng_seed(rng, seed)
    type(rng_state), intent(inout) :: rng
    integer(i8), intent(in) :: seed
    integer(i8) :: s

    s = modulo(abs(seed), pm_m - 1_i8) + 1_i8
    rng%state = s
    rng%has_spare = .false.
    rng%spare = 0.0_dp
  end subroutine rng_seed

  function rng_uniform(rng) result(u)
    type(rng_state), intent(inout) :: rng
    real(dp) :: u
    integer(i8) :: hi, lo, test

    hi = rng%state / pm_q
    lo = modulo(rng%state, pm_q)
    test = pm_a*lo - pm_r*hi
    if (test > 0_i8) then
      rng%state = test
    else
      rng%state = test + pm_m
    end if
    u = real(rng%state, dp)/real(pm_m, dp)
  end function rng_uniform

  function rng_normal(rng) result(z)
    type(rng_state), intent(inout) :: rng
    real(dp) :: z
    real(dp) :: u1, u2, radius

    if (rng%has_spare) then
      z = rng%spare
      rng%has_spare = .false.
      return
    end if

    u1 = max(rng_uniform(rng), tiny(1.0_dp))
    u2 = rng_uniform(rng)
    radius = sqrt(-2.0_dp*log(u1))
    z = radius*cos(2.0_dp*pi*u2)
    rng%spare = radius*sin(2.0_dp*pi*u2)
    rng%has_spare = .true.
  end function rng_normal

  recursive function rng_gamma(rng, shape, scale) result(x)
    type(rng_state), intent(inout) :: rng
    real(dp), intent(in) :: shape, scale
    real(dp) :: x
    real(dp) :: d, c, z, v, u

    if (shape <= 0.0_dp .or. scale < 0.0_dp) then
      x = 0.0_dp
      return
    end if
    if (scale <= tiny(1.0_dp)) then
      x = 0.0_dp
      return
    end if

    if (shape < 1.0_dp) then
      u = max(rng_uniform(rng), tiny(1.0_dp))
      x = rng_gamma(rng, shape + 1.0_dp, scale)*u**(1.0_dp/shape)
      return
    end if

    d = shape - 1.0_dp/3.0_dp
    c = 1.0_dp/sqrt(9.0_dp*d)
    do
      z = rng_normal(rng)
      v = 1.0_dp + c*z
      if (v <= 0.0_dp) cycle
      v = v*v*v
      u = rng_uniform(rng)
      if (u < 1.0_dp - 0.0331_dp*z**4) exit
      if (log(max(u, tiny(1.0_dp))) < 0.5_dp*z*z + d*(1.0_dp - v + log(v))) exit
    end do
    x = scale*d*v
  end function rng_gamma

  function rng_chisq(rng, df) result(x)
    type(rng_state), intent(inout) :: rng
    real(dp), intent(in) :: df
    real(dp) :: x

    if (df <= 0.0_dp) then
      x = 0.0_dp
    else
      x = rng_gamma(rng, 0.5_dp*df, 2.0_dp)
    end if
  end function rng_chisq

  function rng_poisson(rng, mean) result(k)
    type(rng_state), intent(inout) :: rng
    real(dp), intent(in) :: mean
    integer :: k
    real(dp) :: product, threshold, y, em, t
    real(dp) :: sq, alxm, g

    if (mean <= 0.0_dp) then
      k = 0
      return
    end if

    if (mean < 12.0_dp) then
      threshold = exp(-mean)
      product = 1.0_dp
      k = -1
      do
        k = k + 1
        product = product*rng_uniform(rng)
        if (product <= threshold) exit
      end do
      return
    end if

    sq = sqrt(2.0_dp*mean)
    alxm = log(mean)
    g = mean*alxm - log_gamma(mean + 1.0_dp)
    do
      y = tan(pi*rng_uniform(rng))
      em = floor(sq*y + mean)
      if (em < 0.0_dp) cycle
      t = 0.9_dp*(1.0_dp + y*y)*exp(em*alxm - log_gamma(em + 1.0_dp) - g)
      if (rng_uniform(rng) <= t) exit
    end do
    k = int(em)
  end function rng_poisson

end module jumptest_rng
