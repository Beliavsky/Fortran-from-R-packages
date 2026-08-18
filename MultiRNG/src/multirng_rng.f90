! SPDX-License-Identifier: GPL-2.0-only OR GPL-3.0-only
module multirng_rng
  use multirng_kinds, only : dp
  implicit none
  private

  real(dp), parameter :: pi = acos(-1.0_dp)
  logical, save :: have_spare_normal = .false.
  real(dp), save :: spare_normal = 0.0_dp

  public :: seed_rng
  public :: rng_uniform, rng_normal, rng_gamma, rng_chisq
  public :: rng_poisson, rng_binomial, rng_hypergeometric

contains

  subroutine seed_rng(seed)
    integer, intent(in) :: seed
    integer :: n, i
    integer, allocatable :: put(:)

    call random_seed(size=n)
    allocate(put(n))
    do i = 1, n
      put(i) = modulo(seed + 104729 * i + 7919 * i * i, huge(1) - 1)
      if (put(i) <= 0) put(i) = i
    end do
    call random_seed(put=put)
    have_spare_normal = .false.
  end subroutine seed_rng

  function rng_uniform() result(x)
    real(dp) :: x
    call random_number(x)
    if (x <= 0.0_dp) x = tiny(1.0_dp)
    if (x >= 1.0_dp) x = 1.0_dp - epsilon(1.0_dp)
  end function rng_uniform

  function rng_normal() result(z)
    real(dp) :: z
    real(dp) :: u1, u2, r

    if (have_spare_normal) then
      z = spare_normal
      have_spare_normal = .false.
      return
    end if

    u1 = rng_uniform()
    u2 = rng_uniform()
    r = sqrt(-2.0_dp * log(u1))
    z = r * cos(2.0_dp * pi * u2)
    spare_normal = r * sin(2.0_dp * pi * u2)
    have_spare_normal = .true.
  end function rng_normal

  recursive function rng_gamma(shape, scale) result(x)
    real(dp), intent(in) :: shape, scale
    real(dp) :: x
    real(dp) :: d, c, z, v, u

    if (shape <= 0.0_dp .or. scale <= 0.0_dp) error stop "rng_gamma: parameters must be positive"

    if (shape < 1.0_dp) then
      u = rng_uniform()
      x = rng_gamma(shape + 1.0_dp, scale) * u ** (1.0_dp / shape)
      return
    end if

    d = shape - 1.0_dp / 3.0_dp
    c = 1.0_dp / sqrt(9.0_dp * d)
    do
      do
        z = rng_normal()
        v = 1.0_dp + c * z
        if (v > 0.0_dp) exit
      end do
      v = v * v * v
      u = rng_uniform()
      if (u < 1.0_dp - 0.0331_dp * z ** 4) exit
      if (log(u) < 0.5_dp * z * z + d * (1.0_dp - v + log(v))) exit
    end do
    x = scale * d * v
  end function rng_gamma

  function rng_chisq(df) result(x)
    real(dp), intent(in) :: df
    real(dp) :: x
    if (df <= 0.0_dp) error stop "rng_chisq: df must be positive"
    x = rng_gamma(0.5_dp * df, 2.0_dp)
  end function rng_chisq

  function rng_poisson(lambda) result(k)
    real(dp), intent(in) :: lambda
    integer :: k
    real(dp) :: l, p, u, z, proposal, log_accept
    real(dp) :: c, beta, alpha, k0

    if (lambda < 0.0_dp) error stop "rng_poisson: lambda must be nonnegative"
    if (lambda <= 0.0_dp) then
      k = 0
      return
    end if

    if (lambda < 30.0_dp) then
      l = exp(-lambda)
      p = 1.0_dp
      k = 0
      do
        k = k + 1
        p = p * rng_uniform()
        if (p <= l) exit
      end do
      k = k - 1
      return
    end if

    ! Atkinson's transformed rejection method for moderate/large lambda.
    c = 0.767_dp - 3.36_dp / lambda
    beta = pi / sqrt(3.0_dp * lambda)
    alpha = beta * lambda
    k0 = log(c) - lambda - log(beta)
    do
      u = rng_uniform()
      z = log(u / (1.0_dp - u)) / beta
      proposal = alpha + z
      if (proposal < -0.5_dp) cycle
      k = nint(proposal)
      u = rng_uniform()
      log_accept = alpha - beta * proposal + log(u / (1.0_dp + exp(alpha - beta * proposal)) ** 2)
      if (log_accept <= k0 + real(k, dp) * log(lambda) - log_gamma(real(k + 1, dp))) return
    end do
  end function rng_poisson

  function rng_binomial(n, p) result(k)
    integer, intent(in) :: n
    real(dp), intent(in) :: p
    integer :: k, i

    if (n < 0) error stop "rng_binomial: n must be nonnegative"
    if (p < 0.0_dp .or. p > 1.0_dp) error stop "rng_binomial: p outside [0,1]"
    if (p <= 0.0_dp) then
      k = 0
      return
    end if
    if (p >= 1.0_dp) then
      k = n
      return
    end if

    k = 0
    do i = 1, n
      if (rng_uniform() < p) k = k + 1
    end do
  end function rng_binomial

  function rng_hypergeometric(white, black, draws) result(k)
    integer, intent(in) :: white, black, draws
    integer :: k
    integer :: i, wleft, bleft

    if (white < 0 .or. black < 0 .or. draws < 0) error stop "rng_hypergeometric: negative argument"
    if (draws > white + black) error stop "rng_hypergeometric: too many draws"

    wleft = white
    bleft = black
    k = 0
    do i = 1, draws
      if (wleft <= 0) exit
      if (bleft <= 0) then
        k = k + draws - i + 1
        exit
      end if
      if (rng_uniform() < real(wleft, dp) / real(wleft + bleft, dp)) then
        k = k + 1
        wleft = wleft - 1
      else
        bleft = bleft - 1
      end if
    end do
  end function rng_hypergeometric

end module multirng_rng
