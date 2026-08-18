! SPDX-License-Identifier: MIT
module chyper_mle
  use chyper_kinds, only : dp
  implicit none
  private
  public :: mle_s, mle_n, mle_m, loglik_chyper
contains
  real(dp) function loglik_chyper(k, s, n, m) result(ll)
    use chyper_distribution, only : chyper_probabilities
    integer, intent(in) :: k(:), s
    integer, intent(in) :: n(:), m(:)
    real(dp), allocatable :: prob(:)
    integer :: i, stat

    call chyper_probabilities(s, n, m, prob, stat)
    if (stat /= 0) then
      ll = -huge(1.0_dp)
      return
    end if
    ll = 0.0_dp
    do i = 1, size(k)
      if (k(i) < 0 .or. k(i) > ubound(prob, 1)) then
        ll = -huge(1.0_dp)
        return
      end if
      if (prob(k(i)) <= 0.0_dp) then
        ll = -huge(1.0_dp)
        return
      end if
      ll = ll + log(prob(k(i)))
    end do
  end function loglik_chyper

  integer function mle_s(k, n, m) result(best)
    integer, intent(in) :: k(:), n(:), m(:)
    integer :: s, lower
    real(dp) :: cur, nxt
    lower = max(max(0, maxval(m - n)), maxval(k))
    s = lower
    cur = loglik_chyper(k, s, n, m)
    do
      nxt = loglik_chyper(k, s + 1, n, m)
      if (nxt <= cur) exit
      s = s + 1
      cur = nxt
      if (s > 10000000) exit
    end do
    best = s
  end function mle_s

  real(dp) function mle_n(population, k, s, n, m) result(best)
    integer, intent(in) :: population, k(:), s
    integer, intent(in) :: n(:), m(:)
    integer, allocatable :: nt(:)
    integer :: x
    real(dp) :: cur, nxt
    if (all(k == 0)) then
      best = huge(1.0_dp)
      return
    end if
    allocate(nt(size(n)))
    nt = n
    x = max(m(population) - s, 0)
    nt(population) = x
    cur = loglik_chyper(k, s, nt, m)
    do
      nt(population) = x + 1
      nxt = loglik_chyper(k, s, nt, m)
      if (nxt <= cur) exit
      x = x + 1
      cur = nxt
      if (x > 10000000) exit
    end do
    best = real(x, dp)
  end function mle_n

  integer function mle_m(population, k, s, n, m) result(best)
    integer, intent(in) :: population, k(:), s
    integer, intent(in) :: n(:), m(:)
    integer, allocatable :: mt(:)
    integer :: x, upper
    real(dp) :: cur, nxt
    allocate(mt(size(m)))
    mt = m
    x = maxval(k)
    upper = s + n(population)
    mt(population) = x
    cur = loglik_chyper(k, s, n, mt)
    do while (x < upper)
      mt(population) = x + 1
      nxt = loglik_chyper(k, s, n, mt)
      if (nxt <= cur) exit
      x = x + 1
      cur = nxt
    end do
    best = x
  end function mle_m
end module chyper_mle
