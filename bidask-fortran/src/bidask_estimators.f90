! SPDX-License-Identifier: MIT
! Based on bidask 2.1.5, Copyright (c) 2024 Emanuele Guidotti.
module bidask_estimators
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use bidask_kinds, only: dp
  use bidask_statistics, only: nan_dp, mean_values, sum_values, signed_root, safe_log
  implicit none
  private
  public :: edge, edge_estimate, ar_estimate, cs_estimate, roll_estimate
  public :: ohlc_estimate, estimate_method

  interface edge
    module procedure edge_estimate
  end interface edge

contains

  real(dp) function edge_estimate(open, high, low, close, signed, na_rm) result(spread)
    real(dp), intent(in) :: open(:), high(:), low(:), close(:)
    logical, intent(in), optional :: signed, na_rm
    logical :: keep_sign, remove_na
    integer :: n, i
    real(dp) :: pt, po, pc, nt, e1, e2, v1, v2, vt, s2
    real(dp) :: mr1, mr3, mr5
    real(dp), allocatable :: o(:), h(:), l(:), c(:), m(:)
    real(dp), allocatable :: r1(:), r2(:), r3(:), r4(:), r5(:)
    real(dp), allocatable :: tau(:), po1(:), po2(:), pc1(:), pc2(:)
    real(dp), allocatable :: d1(:), d3(:), d5(:), x1(:), x2(:)

    keep_sign = .false.
    if (present(signed)) keep_sign = signed
    remove_na = .true.
    if (present(na_rm)) remove_na = na_rm
    spread = nan_dp()
    n = size(open)
    if (n < 3) return
    if (size(high) /= n .or. size(low) /= n .or. size(close) /= n) return

    allocate(o(n), h(n), l(n), c(n), m(n))
    do i = 1, n
      o(i) = safe_log(open(i))
      h(i) = safe_log(high(i))
      l(i) = safe_log(low(i))
      c(i) = safe_log(close(i))
      m(i) = (h(i) + l(i)) / 2.0_dp
    end do

    allocate(r1(n - 1), r2(n - 1), r3(n - 1), r4(n - 1), r5(n - 1))
    allocate(tau(n - 1), po1(n - 1), po2(n - 1), pc1(n - 1), pc2(n - 1))
    do i = 1, n - 1
      r1(i) = m(i + 1) - o(i + 1)
      r2(i) = o(i + 1) - m(i)
      r3(i) = m(i + 1) - c(i)
      r4(i) = c(i) - m(i)
      r5(i) = o(i + 1) - c(i)
      if (ieee_is_finite(h(i + 1)) .and. ieee_is_finite(l(i + 1)) .and. &
          ieee_is_finite(c(i))) then
        if (real_different(h(i + 1), l(i + 1)) .or. real_different(l(i + 1), c(i))) then
          tau(i) = 1.0_dp
        else
          tau(i) = 0.0_dp
        end if
      else
        tau(i) = nan_dp()
      end if
      po1(i) = indicator_product(tau(i), real_different(o(i + 1), h(i + 1)))
      po2(i) = indicator_product(tau(i), real_different(o(i + 1), l(i + 1)))
      pc1(i) = indicator_product(tau(i), real_different(c(i), h(i)))
      pc2(i) = indicator_product(tau(i), real_different(c(i), l(i)))
    end do

    pt = mean_values(tau, remove_na)
    po = mean_values(po1, remove_na) + mean_values(po2, remove_na)
    pc = mean_values(pc1, remove_na) + mean_values(pc2, remove_na)
    nt = sum_values(tau, .true.)
    if (.not. ieee_is_finite(pt) .or. .not. ieee_is_finite(po) .or. &
        .not. ieee_is_finite(pc) .or. .not. ieee_is_finite(nt)) return
    if (nt < 2.0_dp .or. real_zero(po) .or. real_zero(pc) .or. real_zero(pt)) return

    mr1 = mean_values(r1, remove_na)
    mr3 = mean_values(r3, remove_na)
    mr5 = mean_values(r5, remove_na)
    if (.not. ieee_is_finite(mr1) .or. .not. ieee_is_finite(mr3) .or. &
        .not. ieee_is_finite(mr5)) return

    allocate(d1(n - 1), d3(n - 1), d5(n - 1), x1(n - 1), x2(n - 1))
    do i = 1, n - 1
      d1(i) = r1(i) - mr1 / pt * tau(i)
      d3(i) = r3(i) - mr3 / pt * tau(i)
      d5(i) = r5(i) - mr5 / pt * tau(i)
      x1(i) = -4.0_dp / po * d1(i) * r2(i) - 4.0_dp / pc * d3(i) * r4(i)
      x2(i) = -4.0_dp / po * d1(i) * r5(i) - 4.0_dp / pc * d5(i) * r4(i)
    end do

    e1 = mean_values(x1, remove_na)
    e2 = mean_values(x2, remove_na)
    v1 = mean_values(x1 * x1, remove_na) - e1 * e1
    v2 = mean_values(x2 * x2, remove_na) - e2 * e2
    if (.not. ieee_is_finite(e1) .or. .not. ieee_is_finite(e2) .or. &
        .not. ieee_is_finite(v1) .or. .not. ieee_is_finite(v2)) return
    vt = v1 + v2
    if (ieee_is_finite(vt) .and. vt > 0.0_dp) then
      s2 = (v2 * e1 + v1 * e2) / vt
    else
      s2 = (e1 + e2) / 2.0_dp
    end if
    spread = signed_root(s2, keep_sign)
  end function edge_estimate

  real(dp) function ar_estimate(high, low, close, two_period, signed, na_rm) result(spread)
    real(dp), intent(in) :: high(:), low(:), close(:)
    logical, intent(in), optional :: two_period, signed, na_rm
    logical :: use_two, keep_sign, remove_na
    integer :: n, i
    real(dp) :: avg
    real(dp), allocatable :: h(:), l(:), c(:), mid(:), s2(:), s(:)

    spread = nan_dp()
    use_two = .false.
    keep_sign = .false.
    remove_na = .false.
    if (present(two_period)) use_two = two_period
    if (present(signed)) keep_sign = signed
    if (present(na_rm)) remove_na = na_rm
    n = size(high)
    if (n < 2 .or. size(low) /= n .or. size(close) /= n) return
    allocate(h(n), l(n), c(n), mid(n), s2(n - 1))
    do i = 1, n
      h(i) = safe_log(high(i))
      l(i) = safe_log(low(i))
      c(i) = safe_log(close(i))
      mid(i) = (h(i) + l(i)) / 2.0_dp
    end do
    do i = 1, n - 1
      s2(i) = 4.0_dp * (c(i) - mid(i)) * (c(i) - mid(i + 1))
    end do
    if (use_two) then
      allocate(s(n - 1))
      do i = 1, n - 1
        if (.not. ieee_is_finite(s2(i))) then
          s(i) = nan_dp()
        else
          s(i) = sqrt(max(s2(i), 0.0_dp))
        end if
      end do
      spread = mean_values(s, remove_na)
    else
      avg = mean_values(s2, remove_na)
      spread = signed_root(avg, .true.)
      if (.not. keep_sign) spread = abs(spread)
    end if
  end function ar_estimate

  real(dp) function cs_estimate(high, low, close, two_period, signed, na_rm) result(spread)
    real(dp), intent(in) :: high(:), low(:), close(:)
    logical, intent(in), optional :: two_period, signed, na_rm
    logical :: use_two, keep_sign, remove_na
    integer :: n, i
    real(dp), parameter :: root2 = 1.414213562373095048801688724210_dp
    real(dp) :: denom
    real(dp), allocatable :: h(:), l(:), c(:), raw(:)
    real(dp) :: gap, ah, al, b, g, a

    spread = nan_dp()
    use_two = .false.
    keep_sign = .false.
    remove_na = .false.
    if (present(two_period)) use_two = two_period
    if (present(signed)) keep_sign = signed
    if (present(na_rm)) remove_na = na_rm
    n = size(high)
    if (n < 2 .or. size(low) /= n .or. size(close) /= n) return
    allocate(h(n), l(n), c(n), raw(n - 1))
    do i = 1, n
      h(i) = safe_log(high(i))
      l(i) = safe_log(low(i))
      c(i) = safe_log(close(i))
    end do
    denom = 3.0_dp - 2.0_dp * root2
    do i = 2, n
      if (.not. all_finite4(h(i), l(i), h(i - 1), l(i - 1)) .or. &
          .not. ieee_is_finite(c(i - 1))) then
        raw(i - 1) = nan_dp()
      else
        gap = max(0.0_dp, c(i - 1) - h(i)) + min(0.0_dp, c(i - 1) - l(i))
        ah = h(i) + gap
        al = l(i) + gap
        b = (h(i) - l(i))**2 + (h(i - 1) - l(i - 1))**2
        g = (max(ah, h(i - 1)) - min(al, l(i - 1)))**2
        a = (sqrt(2.0_dp * b) - sqrt(b)) / denom - sqrt(g / denom)
        raw(i - 1) = 2.0_dp * (exp(a) - 1.0_dp) / (1.0_dp + exp(a))
        if (use_two .and. raw(i - 1) < 0.0_dp) raw(i - 1) = 0.0_dp
      end if
    end do
    spread = mean_values(raw, remove_na)
    if (.not. use_two .and. .not. keep_sign) spread = abs(spread)
  end function cs_estimate

  real(dp) function roll_estimate(close, signed, na_rm) result(spread)
    real(dp), intent(in) :: close(:)
    logical, intent(in), optional :: signed, na_rm
    logical :: keep_sign, remove_na
    integer :: n, i, count_valid
    real(dp) :: m1, m2, m12, s2
    real(dp), allocatable :: c(:), r1(:), r2(:), prod(:)

    spread = nan_dp()
    keep_sign = .false.
    remove_na = .false.
    if (present(signed)) keep_sign = signed
    if (present(na_rm)) remove_na = na_rm
    n = size(close)
    if (n < 3) return
    allocate(c(n), r1(n - 2), r2(n - 2), prod(n - 2))
    do i = 1, n
      c(i) = safe_log(close(i))
    end do
    do i = 3, n
      r1(i - 2) = c(i) - c(i - 1)
      r2(i - 2) = c(i - 1) - c(i - 2)
      prod(i - 2) = r1(i - 2) * r2(i - 2)
    end do
    m1 = mean_values(r1, remove_na)
    m2 = mean_values(r2, remove_na)
    m12 = mean_values(prod, remove_na)
    count_valid = count(ieee_is_finite(r2))
    if (count_valid < 2 .or. .not. ieee_is_finite(m1) .or. &
        .not. ieee_is_finite(m2) .or. .not. ieee_is_finite(m12)) return
    s2 = -4.0_dp * real(count_valid, dp) / real(count_valid - 1, dp) * (m12 - m1 * m2)
    spread = signed_root(s2, .true.)
    if (.not. keep_sign) spread = abs(spread)
  end function roll_estimate

  real(dp) function ohlc_estimate(open, high, low, close, method, signed, na_rm) result(spread)
    real(dp), intent(in) :: open(:), high(:), low(:), close(:)
    character(len=*), intent(in) :: method
    logical, intent(in), optional :: signed, na_rm
    logical :: keep_sign, remove_na
    integer :: n, i, n_parts, start_pos, dot_pos
    character(len=64) :: work
    character(len=16) :: part
    real(dp) :: sum_s2, value
    real(dp), allocatable :: o(:), h(:), l(:), c(:), mid(:), tau(:)
    real(dp) :: pt, nt, po, pc

    spread = nan_dp()
    keep_sign = .false.
    remove_na = .false.
    if (present(signed)) keep_sign = signed
    if (present(na_rm)) remove_na = na_rm
    n = size(open)
    if (n < 2 .or. size(high) /= n .or. size(low) /= n .or. size(close) /= n) return
    allocate(o(n), h(n), l(n), c(n), mid(n), tau(n - 1))
    do i = 1, n
      o(i) = safe_log(open(i))
      h(i) = safe_log(high(i))
      l(i) = safe_log(low(i))
      c(i) = safe_log(close(i))
      mid(i) = (h(i) + l(i)) / 2.0_dp
    end do
    do i = 2, n
      if (ieee_is_finite(h(i)) .and. ieee_is_finite(l(i)) .and. ieee_is_finite(c(i - 1))) then
        if (real_different(h(i), l(i)) .or. real_different(l(i), c(i - 1))) then
          tau(i - 1) = 1.0_dp
        else
          tau(i - 1) = 0.0_dp
        end if
      else
        tau(i - 1) = nan_dp()
      end if
    end do
    pt = mean_values(tau, remove_na)
    nt = sum_values(tau, .true.)
    po = event_probability(tau, o(2:n), h(2:n), remove_na) + &
      event_probability(tau, o(2:n), l(2:n), remove_na)
    pc = event_probability(tau, c(1:n - 1), h(1:n - 1), remove_na) + &
      event_probability(tau, c(1:n - 1), l(1:n - 1), remove_na)
    if (.not. ieee_is_finite(pt) .or. .not. ieee_is_finite(nt)) return
    if (nt < 2.0_dp .or. real_zero(pt)) return

    work = upper_string(trim(method))
    sum_s2 = 0.0_dp
    n_parts = 0
    start_pos = 1
    do
      dot_pos = index(work(start_pos:), '.')
      if (dot_pos == 0) then
        part = adjustl(work(start_pos:len_trim(work)))
      else
        part = adjustl(work(start_pos:start_pos + dot_pos - 2))
      end if
      value = component_s2(trim(part), o, c, mid, tau, pt, nt, po, pc, remove_na)
      if (.not. ieee_is_finite(value)) return
      sum_s2 = sum_s2 + value
      n_parts = n_parts + 1
      if (dot_pos == 0) exit
      start_pos = start_pos + dot_pos
      if (start_pos > len_trim(work)) return
    end do
    if (n_parts == 0) return
    spread = signed_root(sum_s2 / real(n_parts, dp), keep_sign)
  end function ohlc_estimate

  real(dp) function estimate_method(open, high, low, close, method, signed, na_rm) result(spread)
    real(dp), intent(in) :: open(:), high(:), low(:), close(:)
    character(len=*), intent(in) :: method
    logical, intent(in), optional :: signed, na_rm
    logical :: keep_sign, remove_na
    character(len=64) :: name

    keep_sign = .false.
    remove_na = .false.
    if (present(signed)) keep_sign = signed
    if (present(na_rm)) remove_na = na_rm
    name = upper_string(trim(method))
    select case (trim(name))
    case ('EDGE')
      spread = edge_estimate(open, high, low, close, keep_sign, remove_na)
    case ('AR')
      spread = ar_estimate(high, low, close, .false., keep_sign, remove_na)
    case ('AR2')
      spread = ar_estimate(high, low, close, .true., keep_sign, remove_na)
    case ('CS')
      spread = cs_estimate(high, low, close, .false., keep_sign, remove_na)
    case ('CS2')
      spread = cs_estimate(high, low, close, .true., keep_sign, remove_na)
    case ('ROLL')
      spread = roll_estimate(close, keep_sign, remove_na)
    case default
      spread = ohlc_estimate(open, high, low, close, trim(name), keep_sign, remove_na)
    end select
  end function estimate_method

  pure real(dp) function indicator_product(tau, condition) result(value)
    real(dp), intent(in) :: tau
    logical, intent(in) :: condition
    if (.not. ieee_is_finite(tau)) then
      value = nan_dp()
    else if (condition) then
      value = tau
    else
      value = 0.0_dp
    end if
  end function indicator_product

  pure logical function all_finite4(a, b, c, d) result(ok)
    real(dp), intent(in) :: a, b, c, d
    ok = ieee_is_finite(a) .and. ieee_is_finite(b) .and. &
      ieee_is_finite(c) .and. ieee_is_finite(d)
  end function all_finite4

  real(dp) function event_probability(tau, a, b, na_rm) result(probability)
    real(dp), intent(in) :: tau(:), a(:), b(:)
    logical, intent(in) :: na_rm
    integer :: i
    real(dp), allocatable :: x(:)
    allocate(x(size(tau)))
    do i = 1, size(tau)
      x(i) = indicator_product(tau(i), real_different(a(i), b(i)))
    end do
    probability = mean_values(x, na_rm)
  end function event_probability

  real(dp) function component_s2(name, o, c, mid, tau, pt, nt, po, pc, na_rm) result(value)
    character(len=*), intent(in) :: name
    real(dp), intent(in) :: o(:), c(:), mid(:), tau(:)
    real(dp), intent(in) :: pt, nt, po, pc
    logical, intent(in) :: na_rm
    integer :: n
    real(dp), allocatable :: r1(:), r2(:), prod(:), tr2(:)
    real(dp) :: p, mprod, mr1, mtr2

    value = nan_dp()
    n = size(o)
    allocate(r1(n - 1), r2(n - 1), prod(n - 1), tr2(n - 1))
    select case (trim(name))
    case ('OHL')
      r1 = mid(2:n) - o(2:n)
      r2 = o(2:n) - mid(1:n - 1)
      p = po
    case ('OHLC')
      r1 = mid(2:n) - o(2:n)
      r2 = o(2:n) - c(1:n - 1)
      p = po
    case ('CHL')
      r1 = mid(2:n) - c(1:n - 1)
      r2 = c(1:n - 1) - mid(1:n - 1)
      p = pc
    case ('CHLO')
      r1 = o(2:n) - c(1:n - 1)
      r2 = c(1:n - 1) - mid(1:n - 1)
      p = pc
    case default
      return
    end select
    if (.not. ieee_is_finite(p) .or. real_zero(p) .or. nt < 2.0_dp) return
    prod = r1 * r2
    tr2 = tau * r2
    mprod = mean_values(prod, na_rm)
    mr1 = mean_values(r1, na_rm)
    mtr2 = mean_values(tr2, na_rm)
    if (.not. ieee_is_finite(mprod) .or. .not. ieee_is_finite(mr1) .or. &
        .not. ieee_is_finite(mtr2)) return
    value = -8.0_dp / p * (mprod - mr1 * mtr2 / pt)
  end function component_s2

  pure logical function real_different(a, b) result(different)
    real(dp), intent(in) :: a, b
    if (.not. ieee_is_finite(a) .or. .not. ieee_is_finite(b)) then
      different = .false.
    else
      different = a < b .or. a > b
    end if
  end function real_different

  pure logical function real_zero(x) result(is_zero)
    real(dp), intent(in) :: x
    is_zero = ieee_is_finite(x) .and. x <= 0.0_dp .and. x >= 0.0_dp
  end function real_zero

  pure function upper_string(text) result(out)
    character(len=*), intent(in) :: text
    character(len=len(text)) :: out
    integer :: i, code
    do i = 1, len(text)
      code = iachar(text(i:i))
      if (code >= iachar('a') .and. code <= iachar('z')) then
        out(i:i) = achar(code - 32)
      else
        out(i:i) = text(i:i)
      end if
    end do
  end function upper_string

end module bidask_estimators
