! SPDX-License-Identifier: GPL-3.0-or-later
module sgt_special
  use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
  use sgt_kinds, only : dp
  implicit none
  private
  public :: reg_lower_gamma, gamma_quantile, reg_incomplete_beta, beta_quantile
  public :: normal_cdf, invert_matrix
contains
  pure real(dp) function nan_dp() result(x)
    x = ieee_value(0.0_dp, ieee_quiet_nan)
  end function nan_dp

  pure real(dp) function reg_lower_gamma(a, x) result(p)
    real(dp), intent(in) :: a, x
    integer, parameter :: maxit = 10000
    real(dp), parameter :: eps = 5.0e-15_dp
    real(dp), parameter :: fpmin = tiny(1.0_dp) / eps
    integer :: n
    real(dp) :: ap, del, sumv, b, c, d, h, an, gln
    if (a <= 0.0_dp .or. x < 0.0_dp) then
      p = nan_dp()
      return
    end if
    if (x <= 0.0_dp) then
      p = 0.0_dp
      return
    end if
    gln = log_gamma(a)
    if (x < a + 1.0_dp) then
      ap = a
      sumv = 1.0_dp / a
      del = sumv
      do n = 1, maxit
        ap = ap + 1.0_dp
        del = del * x / ap
        sumv = sumv + del
        if (abs(del) <= abs(sumv) * eps) exit
      end do
      p = sumv * exp(-x + a * log(x) - gln)
    else
      b = x + 1.0_dp - a
      c = 1.0_dp / fpmin
      d = 1.0_dp / b
      h = d
      do n = 1, maxit
        an = -real(n, dp) * (real(n, dp) - a)
        b = b + 2.0_dp
        d = an * d + b
        if (abs(d) < fpmin) d = fpmin
        c = b + an / c
        if (abs(c) < fpmin) c = fpmin
        d = 1.0_dp / d
        del = d * c
        h = h * del
        if (abs(del - 1.0_dp) <= eps) exit
      end do
      p = 1.0_dp - exp(-x + a * log(x) - gln) * h
    end if
    p = max(0.0_dp, min(1.0_dp, p))
  end function reg_lower_gamma

  pure real(dp) function betacf(a, b, x) result(cf)
    real(dp), intent(in) :: a, b, x
    integer, parameter :: maxit = 10000
    real(dp), parameter :: eps = 5.0e-15_dp
    real(dp), parameter :: fpmin = tiny(1.0_dp) / eps
    integer :: m, m2
    real(dp) :: aa, c, d, del, h, qab, qam, qap
    qab = a + b
    qap = a + 1.0_dp
    qam = a - 1.0_dp
    c = 1.0_dp
    d = 1.0_dp - qab * x / qap
    if (abs(d) < fpmin) d = fpmin
    d = 1.0_dp / d
    h = d
    do m = 1, maxit
      m2 = 2 * m
      aa = real(m, dp) * (b - real(m, dp)) * x / &
        ((qam + real(m2, dp)) * (a + real(m2, dp)))
      d = 1.0_dp + aa * d
      if (abs(d) < fpmin) d = fpmin
      c = 1.0_dp + aa / c
      if (abs(c) < fpmin) c = fpmin
      d = 1.0_dp / d
      h = h * d * c
      aa = -(a + real(m, dp)) * (qab + real(m, dp)) * x / &
        ((a + real(m2, dp)) * (qap + real(m2, dp)))
      d = 1.0_dp + aa * d
      if (abs(d) < fpmin) d = fpmin
      c = 1.0_dp + aa / c
      if (abs(c) < fpmin) c = fpmin
      d = 1.0_dp / d
      del = d * c
      h = h * del
      if (abs(del - 1.0_dp) <= eps) exit
    end do
    cf = h
  end function betacf

  pure real(dp) function reg_incomplete_beta(x, a, b) result(p)
    real(dp), intent(in) :: x, a, b
    real(dp) :: bt
    if (a <= 0.0_dp .or. b <= 0.0_dp) then
      p = nan_dp()
      return
    end if
    if (x <= 0.0_dp) then
      p = 0.0_dp
      return
    end if
    if (x >= 1.0_dp) then
      p = 1.0_dp
      return
    end if
    bt = exp(log_gamma(a + b) - log_gamma(a) - log_gamma(b) + &
      a * log(x) + b * log(1.0_dp - x))
    if (x < (a + 1.0_dp) / (a + b + 2.0_dp)) then
      p = bt * betacf(a, b, x) / a
    else
      p = 1.0_dp - bt * betacf(b, a, 1.0_dp - x) / b
    end if
    p = max(0.0_dp, min(1.0_dp, p))
  end function reg_incomplete_beta

  pure real(dp) function gamma_quantile(prob, a) result(x)
    real(dp), intent(in) :: prob, a
    integer, parameter :: maxit = 160
    integer :: i
    real(dp) :: lo, hi, mid, pmid
    if (a <= 0.0_dp .or. prob < 0.0_dp .or. prob > 1.0_dp) then
      x = nan_dp()
      return
    end if
    if (prob <= 0.0_dp) then
      x = 0.0_dp
      return
    end if
    if (prob >= 1.0_dp) then
      x = huge(1.0_dp)
      return
    end if
    lo = 0.0_dp
    hi = max(1.0_dp, a)
    do while (reg_lower_gamma(a, hi) < prob .and. hi < huge(hi) / 4.0_dp)
      hi = 2.0_dp * hi
    end do
    do i = 1, maxit
      mid = 0.5_dp * (lo + hi)
      pmid = reg_lower_gamma(a, mid)
      if (pmid < prob) then
        lo = mid
      else
        hi = mid
      end if
      if (abs(hi - lo) <= 2.0e-13_dp * max(1.0_dp, mid)) exit
    end do
    x = 0.5_dp * (lo + hi)
  end function gamma_quantile

  pure real(dp) function beta_quantile(prob, a, b) result(x)
    real(dp), intent(in) :: prob, a, b
    integer, parameter :: maxit = 180
    integer :: i
    real(dp) :: lo, hi, mid, pmid
    if (a <= 0.0_dp .or. b <= 0.0_dp .or. prob < 0.0_dp .or. prob > 1.0_dp) then
      x = nan_dp()
      return
    end if
    if (prob <= 0.0_dp) then
      x = 0.0_dp
      return
    end if
    if (prob >= 1.0_dp) then
      x = 1.0_dp
      return
    end if
    lo = 0.0_dp
    hi = 1.0_dp
    do i = 1, maxit
      mid = 0.5_dp * (lo + hi)
      pmid = reg_incomplete_beta(mid, a, b)
      if (pmid < prob) then
        lo = mid
      else
        hi = mid
      end if
      if (hi - lo <= 5.0e-14_dp) exit
    end do
    x = 0.5_dp * (lo + hi)
  end function beta_quantile

  pure real(dp) function normal_cdf(x) result(p)
    real(dp), intent(in) :: x
    p = 0.5_dp * erfc(-x / sqrt(2.0_dp))
  end function normal_cdf

  subroutine invert_matrix(a, ainv, status)
    real(dp), intent(in) :: a(:,:)
    real(dp), allocatable, intent(out) :: ainv(:,:)
    integer, intent(out) :: status
    real(dp), allocatable :: aug(:,:)
    real(dp) :: piv, factor, temp
    integer :: n, i, j, k, pivot_row
    n = size(a, 1)
    status = 0
    if (size(a, 2) /= n) then
      status = 1
      allocate(ainv(0,0))
      return
    end if
    allocate(aug(n, 2*n), ainv(n,n))
    aug(:, 1:n) = a
    aug(:, n+1:2*n) = 0.0_dp
    do i = 1, n
      aug(i, n+i) = 1.0_dp
    end do
    do i = 1, n
      pivot_row = i
      do k = i + 1, n
        if (abs(aug(k,i)) > abs(aug(pivot_row,i))) pivot_row = k
      end do
      if (abs(aug(pivot_row,i)) <= 100.0_dp * epsilon(1.0_dp)) then
        status = 2
        ainv = 0.0_dp
        return
      end if
      if (pivot_row /= i) then
        do j = 1, 2*n
          temp = aug(i,j)
          aug(i,j) = aug(pivot_row,j)
          aug(pivot_row,j) = temp
        end do
      end if
      piv = aug(i,i)
      aug(i,:) = aug(i,:) / piv
      do k = 1, n
        if (k == i) cycle
        factor = aug(k,i)
        aug(k,:) = aug(k,:) - factor * aug(i,:)
      end do
    end do
    ainv = aug(:, n+1:2*n)
  end subroutine invert_matrix
end module sgt_special
