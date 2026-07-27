! SPDX-License-Identifier: GPL-2.0-or-later
! This file is part of a modern Fortran translation of PerformanceAnalytics.
! Copyright (C) original PerformanceAnalytics authors and translation contributors.
! Distributed under GNU GPL version 2 or, at your option, any later version.
module rng_mod
  use kinds_mod, only: dp
  implicit none
  private
  public :: rng_state, rng_seed, rng_uniform, rng_normal, rng_integer
  public :: cholesky_lower, multivariate_normal_draws

  type :: rng_state
    integer(kind=8) :: state = 88172645463325252_8
    logical :: has_spare = .false.
    real(dp) :: spare = 0.0_dp
  end type rng_state
contains
  subroutine rng_seed(rng, seed)
    type(rng_state), intent(inout) :: rng
    integer(kind=8), intent(in) :: seed
    rng%state = seed
    if (rng%state == 0_8) rng%state = 88172645463325252_8
    rng%has_spare = .false.
    rng%spare = 0.0_dp
  end subroutine rng_seed

  real(dp) function rng_uniform(rng) result(u)
    type(rng_state), intent(inout) :: rng
    integer(kind=8) :: x
    x = rng%state
    x = ieor(x, shiftl(x, 13))
    x = ieor(x, shiftr(x, 7))
    x = ieor(x, shiftl(x, 17))
    rng%state = x
    u = real(iand(x, int(z'001FFFFFFFFFFFFF', kind=8)), dp) / real(int(z'0020000000000000', kind=8), dp)
    if (u <= 0.0_dp) u = epsilon(1.0_dp)
  end function rng_uniform

  real(dp) function rng_normal(rng) result(z)
    type(rng_state), intent(inout) :: rng
    real(dp) :: u1, u2, radius
    real(dp), parameter :: two_pi = 6.283185307179586476925286766559_dp
    if (rng%has_spare) then
      z = rng%spare
      rng%has_spare = .false.
      return
    end if
    u1 = rng_uniform(rng)
    u2 = rng_uniform(rng)
    radius = sqrt(-2.0_dp * log(u1))
    z = radius * cos(two_pi * u2)
    rng%spare = radius * sin(two_pi * u2)
    rng%has_spare = .true.
  end function rng_normal

  integer function rng_integer(rng, lo, hi) result(value)
    type(rng_state), intent(inout) :: rng
    integer, intent(in) :: lo, hi
    if (hi <= lo) then
      value = lo
    else
      value = lo + min(hi-lo, int(rng_uniform(rng) * real(hi-lo+1, dp)))
    end if
  end function rng_integer

  subroutine cholesky_lower(a, l, ok)
    real(dp), intent(in) :: a(:,:)
    real(dp), intent(out) :: l(:,:)
    logical, intent(out) :: ok
    integer :: n, i, j, k
    real(dp) :: s
    n = size(a,1)
    l = 0.0_dp
    ok = size(a,2) == n .and. all(shape(l) == [n,n])
    if (.not. ok) return
    do i = 1, n
      do j = 1, i
        s = a(i,j)
        do k = 1, j-1
          s = s - l(i,k) * l(j,k)
        end do
        if (i == j) then
          if (s <= 100.0_dp * epsilon(1.0_dp)) then
            ok = .false.
            return
          end if
          l(i,j) = sqrt(s)
        else
          l(i,j) = s / l(j,j)
        end if
      end do
    end do
  end subroutine cholesky_lower

  subroutine multivariate_normal_draws(mu, sigma, nsim, seed, draws, ok)
    real(dp), intent(in) :: mu(:), sigma(:,:)
    integer, intent(in) :: nsim
    integer(kind=8), intent(in) :: seed
    real(dp), allocatable, intent(out) :: draws(:,:)
    logical, intent(out) :: ok
    real(dp), allocatable :: l(:,:), adjusted(:,:), z(:)
    type(rng_state) :: rng
    integer :: p, i, j, attempt
    real(dp) :: jitter
    p = size(mu)
    allocate(draws(max(nsim,0),p), l(p,p), adjusted(p,p), z(p))
    draws = 0.0_dp
    if (size(sigma,1) /= p .or. size(sigma,2) /= p .or. nsim < 0) then
      ok = .false.
      return
    end if
    adjusted = 0.5_dp * (sigma + transpose(sigma))
    jitter = max(1.0e-14_dp, 1.0e-12_dp * max(1.0_dp, maxval(abs(adjusted))))
    do attempt = 0, 8
      call cholesky_lower(adjusted, l, ok)
      if (ok) exit
      do j = 1, p
        adjusted(j,j) = adjusted(j,j) + jitter
      end do
      jitter = 10.0_dp * jitter
    end do
    if (.not. ok) return
    call rng_seed(rng, seed)
    do i = 1, nsim
      do j = 1, p
        z(j) = rng_normal(rng)
      end do
      draws(i,:) = mu + matmul(l, z)
    end do
  end subroutine multivariate_normal_draws
end module rng_mod
