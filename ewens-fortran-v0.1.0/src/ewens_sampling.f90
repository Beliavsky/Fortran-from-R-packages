! SPDX-License-Identifier: MIT
module ewens_sampling
  use ewens_kinds, only : dp, i8
  use copula_random, only : seed_random, random_uniform, random_gamma
  implicit none
  private
  public :: ewens_seed, rewens, gcrp, rgem

contains

  subroutine ewens_seed(seed)
    integer(i8), intent(in) :: seed
    call seed_random(seed)
  end subroutine ewens_seed

  function rewens(n, theta) result(assignment)
    integer, intent(in) :: n
    real(dp), intent(in), optional :: theta
    integer, allocatable :: assignment(:)
    real(dp) :: th

    th = 1.0_dp
    if (present(theta)) th = theta
    assignment = gcrp(n, 0.0_dp, th)
  end function rewens

  function gcrp(n, alpha, theta) result(assignment)
    integer, intent(in) :: n
    real(dp), intent(in), optional :: alpha, theta
    integer, allocatable :: assignment(:)
    integer, allocatable :: table_sizes(:)
    real(dp) :: a, th, total, u, cumulative
    integer :: i, j, chosen, n_tables

    a = 0.0_dp
    th = 1.0_dp
    if (present(alpha)) a = alpha
    if (present(theta)) th = theta

    if (n < 0) error stop "gcrp: n must be non-negative"
    if (a < 0.0_dp .or. a >= 1.0_dp) error stop "gcrp: alpha must be in [0,1)"
    if (th < -a) error stop "gcrp: theta must be >= -alpha"

    allocate(assignment(n))
    if (n == 0) return
    allocate(table_sizes(n))
    table_sizes = 0

    assignment(1) = 1
    table_sizes(1) = 1
    n_tables = 1

    do i = 2, n
      total = th + real(i - 1, dp)
      u = random_uniform() * total
      cumulative = 0.0_dp
      chosen = n_tables + 1

      do j = 1, n_tables
        cumulative = cumulative + real(table_sizes(j), dp) - a
        if (u <= cumulative) then
          chosen = j
          exit
        end if
      end do

      if (chosen == n_tables + 1) then
        n_tables = n_tables + 1
        table_sizes(n_tables) = 1
        assignment(i) = n_tables
      else
        table_sizes(chosen) = table_sizes(chosen) + 1
        assignment(i) = chosen
      end if
    end do
  end function gcrp

  function rgem(alpha, theta, trunc_at) result(shares)
    real(dp), intent(in), optional :: alpha, theta
    integer, intent(in), optional :: trunc_at
    real(dp), allocatable :: shares(:)
    real(dp) :: a, th, v, remaining, b
    integer :: k, m

    a = 0.0_dp
    th = 1.0_dp
    m = 500
    if (present(alpha)) a = alpha
    if (present(theta)) th = theta
    if (present(trunc_at)) m = trunc_at

    if (a < 0.0_dp .or. a >= 1.0_dp) error stop "rgem: alpha must be in [0,1)"
    if (th < -a) error stop "rgem: theta must be >= -alpha"
    if (m < 0) error stop "rgem: trunc_at must be non-negative"

    allocate(shares(m))
    if (m == 0) return

    remaining = 1.0_dp
    do k = 1, m
      b = th + real(k, dp) * a
      if (b <= 0.0_dp) then
        v = 1.0_dp
      else
        v = random_beta(1.0_dp - a, b)
      end if
      shares(k) = remaining * v
      remaining = remaining * (1.0_dp - v)
      if (remaining <= epsilon(1.0_dp)) then
        if (k < m) shares(k + 1:m) = 0.0_dp
        exit
      end if
    end do
  end function rgem

  real(dp) function random_beta(a, b) result(x)
    real(dp), intent(in) :: a, b
    real(dp) :: gx, gy

    if (a <= 0.0_dp .or. b <= 0.0_dp) error stop "random_beta: shapes must be positive"
    gx = random_gamma(a, 1.0_dp)
    gy = random_gamma(b, 1.0_dp)
    x = gx / (gx + gy)
  end function random_beta

end module ewens_sampling
