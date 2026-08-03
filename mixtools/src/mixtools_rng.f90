! SPDX-License-Identifier: GPL-2.0-or-later
module mixtools_rng
  use mixtools_kinds, only : dp, pi
  implicit none
  private
  public :: rng_state, rng_seed, random_uniform, random_normal
  public :: random_exponential, random_gamma, random_poisson, random_binomial
  public :: random_weibull, random_dirichlet, random_permutation

  type :: rng_state
    integer(kind=8) :: state = 88172645463393265_8
    logical :: has_spare = .false.
    real(dp) :: spare = 0.0_dp
  end type rng_state
contains
  subroutine rng_seed(rng, seed)
    type(rng_state), intent(inout) :: rng
    integer, intent(in) :: seed
    integer(kind=8) :: z
    z = int(seed,8) + int(z'9E3779B97F4A7C15',8)
    z = ieor(z, shiftr(z,30)) * int(z'BF58476D1CE4E5B9',8)
    z = ieor(z, shiftr(z,27)) * int(z'94D049BB133111EB',8)
    rng%state = ieor(z, shiftr(z,31))
    if (rng%state == 0_8) rng%state = 88172645463393265_8
    rng%has_spare = .false.
  end subroutine rng_seed

  function next_u64(rng) result(x)
    type(rng_state), intent(inout) :: rng
    integer(kind=8) :: x
    x = rng%state
    x = ieor(x, shiftl(x,13))
    x = ieor(x, shiftr(x,7))
    x = ieor(x, shiftl(x,17))
    rng%state = x
  end function next_u64

  function random_uniform(rng) result(u)
    type(rng_state), intent(inout) :: rng
    real(dp) :: u
    integer(kind=8) :: x
    x = next_u64(rng)
    u = real(iand(shiftr(x,11), int(z'001FFFFFFFFFFFFF',8)),dp) / 9007199254740992.0_dp
    u = max(tiny(1.0_dp), min(1.0_dp-epsilon(1.0_dp), u))
  end function random_uniform

  function random_normal(rng) result(z)
    type(rng_state), intent(inout) :: rng
    real(dp) :: z, r, theta
    if (rng%has_spare) then
      z = rng%spare
      rng%has_spare = .false.
      return
    end if
    r = sqrt(-2.0_dp * log(random_uniform(rng)))
    theta = 2.0_dp * pi * random_uniform(rng)
    z = r * cos(theta)
    rng%spare = r * sin(theta)
    rng%has_spare = .true.
  end function random_normal

  function random_exponential(rng, rate) result(x)
    type(rng_state), intent(inout) :: rng
    real(dp), intent(in), optional :: rate
    real(dp) :: x, r
    r = 1.0_dp
    if (present(rate)) r = rate
    x = -log(random_uniform(rng)) / max(r, tiny(1.0_dp))
  end function random_exponential

  recursive function random_gamma(rng, shape, scale) result(x)
    type(rng_state), intent(inout) :: rng
    real(dp), intent(in) :: shape
    real(dp), intent(in), optional :: scale
    real(dp) :: x, sc, d, c, z, u
    sc = 1.0_dp
    if (present(scale)) sc = scale
    if (shape <= 0.0_dp .or. sc <= 0.0_dp) then
      x = 0.0_dp
    else if (shape < 1.0_dp) then
      x = random_gamma(rng, shape + 1.0_dp, sc) * random_uniform(rng)**(1.0_dp/shape)
    else
      d = shape - 1.0_dp/3.0_dp
      c = 1.0_dp / sqrt(9.0_dp*d)
      do
        z = random_normal(rng)
        if (1.0_dp + c*z <= 0.0_dp) cycle
        u = random_uniform(rng)
        if (log(u) < 0.5_dp*z*z + d*(1.0_dp-(1.0_dp+c*z)**3 + log((1.0_dp+c*z)**3))) exit
      end do
      x = sc*d*(1.0_dp+c*z)**3
    end if
  end function random_gamma

  function random_poisson(rng, mean) result(k)
    type(rng_state), intent(inout) :: rng
    real(dp), intent(in) :: mean
    integer :: k
    real(dp) :: p, l, z
    if (mean <= 0.0_dp) then
      k = 0
    else if (mean < 30.0_dp) then
      l = exp(-mean); p = 1.0_dp; k = -1
      do
        k = k + 1; p = p * random_uniform(rng)
        if (p <= l) exit
      end do
    else
      z = mean + sqrt(mean)*random_normal(rng)
      k = max(0, nint(z))
    end if
  end function random_poisson

  function random_binomial(rng, n, prob) result(k)
    type(rng_state), intent(inout) :: rng
    integer, intent(in) :: n
    real(dp), intent(in) :: prob
    integer :: k, i
    k = 0
    do i = 1, max(0,n)
      if (random_uniform(rng) < prob) k = k + 1
    end do
  end function random_binomial

  function random_weibull(rng, shape, scale) result(x)
    type(rng_state), intent(inout) :: rng
    real(dp), intent(in) :: shape, scale
    real(dp) :: x
    x = scale * (-log(random_uniform(rng)))**(1.0_dp/shape)
  end function random_weibull

  subroutine random_dirichlet(rng, alpha, x)
    type(rng_state), intent(inout) :: rng
    real(dp), intent(in) :: alpha(:)
    real(dp), intent(out) :: x(size(alpha))
    integer :: i
    real(dp) :: s
    do i = 1, size(alpha)
      x(i) = random_gamma(rng, alpha(i))
    end do
    s = sum(x)
    if (s > 0.0_dp) x = x/s
  end subroutine random_dirichlet

  subroutine random_permutation(rng, values)
    type(rng_state), intent(inout) :: rng
    integer, intent(inout) :: values(:)
    integer :: i, j, tmp
    do i = size(values), 2, -1
      j = 1 + int(random_uniform(rng)*real(i,dp))
      j = min(i,max(1,j))
      tmp = values(i); values(i) = values(j); values(j) = tmp
    end do
  end subroutine random_permutation
end module mixtools_rng
