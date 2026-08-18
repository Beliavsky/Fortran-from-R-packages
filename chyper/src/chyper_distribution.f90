! SPDX-License-Identifier: MIT
module chyper_distribution
  use chyper_kinds, only : dp
  use chyper_math, only : hypergeom_pmf
  implicit none
  private
  public :: chyper_probabilities, dchyper, pchyper, qchyper, pvalchyper
  public :: dchyper_vec, pchyper_vec, qchyper_vec
contains
  pure logical function valid_model(s, n, m) result(ok)
    integer, intent(in) :: s
    integer, intent(in) :: n(:), m(:)
    integer :: i
    ok = .false.
    if (size(n) < 2 .or. size(m) /= size(n) .or. s < 0) return
    do i = 1, size(n)
      if (n(i) < 0 .or. m(i) < 0) return
      if (m(i) > n(i) + s) return
    end do
    ok = .true.
  end function valid_model

  subroutine chyper_probabilities(s, n, m, prob, status)
    integer, intent(in) :: s
    integer, intent(in) :: n(:), m(:)
    real(dp), allocatable, intent(out) :: prob(:)
    integer, intent(out), optional :: status
    real(dp), allocatable :: old(:), new(:)
    real(dp) :: total
    integer :: i, x, xp, max_old, max_new

    if (.not. valid_model(s, n, m)) then
      allocate(prob(0))
      if (present(status)) status = 1
      return
    end if

    max_old = min(s, m(1))
    allocate(old(0:max_old))
    do x = 0, max_old
      old(x) = hypergeom_pmf(x, s, n(1), m(1))
    end do

    do i = 2, size(n)
      max_new = min(max_old, m(i))
      allocate(new(0:max_new))
      new = 0.0_dp
      do xp = 0, max_old
        if (old(xp) <= 0.0_dp) cycle
        do x = 0, min(xp, m(i))
          new(x) = new(x) + old(xp) * &
              hypergeom_pmf(x, xp, n(i) + s - xp, m(i))
        end do
      end do
      call move_alloc(new, old)
      max_old = max_new
    end do

    total = sum(old)
    if (total > 0.0_dp) old = old / total
    allocate(prob(0:max_old))
    prob = old
    if (present(status)) status = 0
  end subroutine chyper_probabilities

  real(dp) function dchyper(k, s, n, m) result(ans)
    integer, intent(in) :: k, s
    integer, intent(in) :: n(:), m(:)
    real(dp), allocatable :: prob(:)
    integer :: stat
    call chyper_probabilities(s, n, m, prob, stat)
    if (stat /= 0 .or. k < 0 .or. k > ubound(prob, 1)) then
      ans = 0.0_dp
    else
      ans = prob(k)
    end if
  end function dchyper

  subroutine dchyper_vec(k, s, n, m, ans)
    integer, intent(in) :: k(:), s
    integer, intent(in) :: n(:), m(:)
    real(dp), intent(out) :: ans(size(k))
    real(dp), allocatable :: prob(:)
    integer :: i, stat
    call chyper_probabilities(s, n, m, prob, stat)
    ans = 0.0_dp
    if (stat /= 0) return
    do i = 1, size(k)
      if (k(i) >= 0 .and. k(i) <= ubound(prob, 1)) ans(i) = prob(k(i))
    end do
  end subroutine dchyper_vec

  real(dp) function pchyper(k, s, n, m) result(ans)
    integer, intent(in) :: k, s
    integer, intent(in) :: n(:), m(:)
    real(dp), allocatable :: prob(:)
    integer :: stat
    call chyper_probabilities(s, n, m, prob, stat)
    if (stat /= 0) then
      ans = 0.0_dp
    else if (k < 0) then
      ans = 0.0_dp
    else if (k >= ubound(prob, 1)) then
      ans = 1.0_dp
    else
      ans = sum(prob(0:k))
    end if
  end function pchyper

  subroutine pchyper_vec(k, s, n, m, ans)
    integer, intent(in) :: k(:), s
    integer, intent(in) :: n(:), m(:)
    real(dp), intent(out) :: ans(size(k))
    integer :: i
    do i = 1, size(k)
      ans(i) = pchyper(k(i), s, n, m)
    end do
  end subroutine pchyper_vec

  integer function qchyper(p, s, n, m) result(ans)
    real(dp), intent(in) :: p
    integer, intent(in) :: s
    integer, intent(in) :: n(:), m(:)
    real(dp), allocatable :: prob(:)
    real(dp) :: cdf
    integer :: k, stat
    call chyper_probabilities(s, n, m, prob, stat)
    if (stat /= 0) then
      ans = 0
      return
    end if
    if (p <= 0.0_dp) then
      ans = 0
      return
    end if
    if (p >= 1.0_dp) then
      ans = ubound(prob, 1)
      return
    end if
    cdf = 0.0_dp
    do k = 0, ubound(prob, 1)
      cdf = cdf + prob(k)
      if (cdf > p) then
        ans = k
        return
      end if
    end do
    ans = ubound(prob, 1)
  end function qchyper

  subroutine qchyper_vec(p, s, n, m, ans)
    real(dp), intent(in) :: p(:)
    integer, intent(in) :: s
    integer, intent(in) :: n(:), m(:)
    integer, intent(out) :: ans(size(p))
    integer :: i
    do i = 1, size(p)
      ans(i) = qchyper(p(i), s, n, m)
    end do
  end subroutine qchyper_vec

  real(dp) function pvalchyper(k, s, n, m, upper) result(ans)
    integer, intent(in) :: k, s
    integer, intent(in) :: n(:), m(:)
    logical, intent(in), optional :: upper
    logical :: use_upper
    use_upper = .true.
    if (present(upper)) use_upper = upper
    if (use_upper) then
      if (k <= 0) then
        ans = 1.0_dp
      else
        ans = max(0.0_dp, 1.0_dp - pchyper(k - 1, s, n, m))
      end if
    else
      ans = max(0.0_dp, pchyper(k, s, n, m))
    end if
  end function pvalchyper
end module chyper_distribution
