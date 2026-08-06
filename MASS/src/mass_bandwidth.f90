! SPDX-License-Identifier: GPL-3.0-only
module mass_bandwidth
  use rrcov_kinds, only : dp
  use mass_types, only : mass_success, mass_invalid_argument, mass_no_convergence
  use mass_math, only : pi_dp, sample_variance, type7_quantile
  implicit none
  private
  public :: ucv_bandwidth, bcv_bandwidth, width_sj
contains

  subroutine ucv_bandwidth(x, bandwidth, status, lower, upper, tolerance)
    real(dp), intent(in) :: x(:)
    real(dp), intent(out) :: bandwidth
    integer, intent(out) :: status
    real(dp), intent(in), optional :: lower, upper, tolerance
    real(dp) :: lo, hi, tol, hmax
    integer :: n
    n = size(x)
    if (n < 2) then
      bandwidth = 0.0_dp
      status = mass_invalid_argument
      return
    end if
    hmax = 4.0_dp * 1.144_dp * sqrt(max(sample_variance(x), tiny(1.0_dp))) * &
      real(n, dp)**(-0.2_dp)
    lo = 0.1_dp * hmax
    hi = hmax
    if (present(lower)) lo = lower
    if (present(upper)) hi = upper
    tol = 0.1_dp * lo
    if (present(tolerance)) tol = tolerance
    if (lo <= 0.0_dp .or. hi <= lo) then
      bandwidth = 0.0_dp
      status = mass_invalid_argument
      return
    end if
    call golden_minimize(objective, lo, hi, tol, bandwidth)
    status = mass_success
  contains
    function objective(h) result(value)
      real(dp), intent(in) :: h
      real(dp) :: value
      value = ucv_objective(x, h)
    end function objective
  end subroutine ucv_bandwidth

  subroutine bcv_bandwidth(x, bandwidth, status, lower, upper, tolerance)
    real(dp), intent(in) :: x(:)
    real(dp), intent(out) :: bandwidth
    integer, intent(out) :: status
    real(dp), intent(in), optional :: lower, upper, tolerance
    real(dp) :: lo, hi, tol, hmax
    integer :: n
    n = size(x)
    if (n < 2) then
      bandwidth = 0.0_dp
      status = mass_invalid_argument
      return
    end if
    hmax = 4.0_dp * 1.144_dp * sqrt(max(sample_variance(x), tiny(1.0_dp))) * &
      real(n, dp)**(-0.2_dp)
    lo = 0.1_dp * hmax
    hi = hmax
    if (present(lower)) lo = lower
    if (present(upper)) hi = upper
    tol = 0.1_dp * lo
    if (present(tolerance)) tol = tolerance
    if (lo <= 0.0_dp .or. hi <= lo) then
      bandwidth = 0.0_dp
      status = mass_invalid_argument
      return
    end if
    call golden_minimize(objective, lo, hi, tol, bandwidth)
    status = mass_success
  contains
    function objective(h) result(value)
      real(dp), intent(in) :: h
      real(dp) :: value
      value = bcv_objective(x, h)
    end function objective
  end subroutine bcv_bandwidth

  subroutine width_sj(x, bandwidth, status, method, lower, upper, tolerance)
    real(dp), intent(in) :: x(:)
    real(dp), intent(out) :: bandwidth
    integer, intent(out) :: status
    character(len=*), intent(in), optional :: method
    real(dp), intent(in), optional :: lower, upper, tolerance
    character(len=8) :: selected
    real(dp) :: variance, scale, a, b, c1, td, alpha2, lo, hi, tol
    real(dp) :: flo, fhi, mid, fmid, pilot
    integer :: n, iter

    n = size(x)
    if (n < 2) then
      bandwidth = 0.0_dp
      status = mass_invalid_argument
      return
    end if
    variance = sample_variance(x)
    scale = min(sqrt(max(variance, 0.0_dp)), &
      (type7_quantile(x, 0.75_dp) - type7_quantile(x, 0.25_dp)) / 1.349_dp)
    if (scale <= 0.0_dp) then
      bandwidth = 0.0_dp
      status = mass_invalid_argument
      return
    end if
    a = 1.24_dp * scale * real(n, dp)**(-1.0_dp / 7.0_dp)
    b = 1.23_dp * scale * real(n, dp)**(-1.0_dp / 9.0_dp)
    c1 = 1.0_dp / (2.0_dp * sqrt(pi_dp) * real(n, dp))
    td = -phi6_direct(x, b)
    if (td <= 0.0_dp) then
      bandwidth = 0.0_dp
      status = mass_no_convergence
      return
    end if
    alpha2 = 1.357_dp * (phi4_direct(x, a) / td)**(1.0_dp / 7.0_dp)
    selected = "ste"
    if (present(method)) selected = trim(method)
    if (trim(selected) == "dpi") then
      pilot = (2.394_dp / (real(n, dp) * td))**(1.0_dp / 7.0_dp)
      bandwidth = 4.0_dp * (c1 / phi4_direct(x, pilot))**0.2_dp
      status = mass_success
      return
    end if
    hi = 1.144_dp * sqrt(max(variance, tiny(1.0_dp))) * real(n, dp)**(-0.2_dp)
    lo = 0.1_dp * hi
    if (present(lower)) lo = lower
    if (present(upper)) hi = upper
    tol = 0.1_dp * lo
    if (present(tolerance)) tol = tolerance
    flo = sj_equation(lo, x, alpha2, c1)
    fhi = sj_equation(hi, x, alpha2, c1)
    if (flo * fhi > 0.0_dp) then
      bandwidth = 0.0_dp
      status = mass_no_convergence
      return
    end if
    do iter = 1, 100
      mid = 0.5_dp * (lo + hi)
      fmid = sj_equation(mid, x, alpha2, c1)
      if (abs(hi - lo) <= tol) exit
      if (flo * fmid <= 0.0_dp) then
        hi = mid
        fhi = fmid
      else
        lo = mid
        flo = fmid
      end if
    end do
    bandwidth = 4.0_dp * mid
    status = mass_success
  end subroutine width_sj

  function ucv_objective(x, h) result(value)
    real(dp), intent(in) :: x(:), h
    real(dp) :: value, hh, delta, term, total
    integer :: i, j, n
    n = size(x)
    hh = h / 4.0_dp
    total = 0.0_dp
    do i = 1, n - 1
      do j = i + 1, n
        delta = ((x(i) - x(j)) / hh)**2
        if (delta >= 1000.0_dp) cycle
        term = exp(-delta / 4.0_dp) - sqrt(8.0_dp) * exp(-delta / 2.0_dp)
        total = total + term
      end do
    end do
    value = 1.0_dp / (2.0_dp * real(n, dp) * hh * sqrt(pi_dp)) + &
      total / (real(n, dp)**2 * hh * sqrt(pi_dp))
  end function ucv_objective

  function bcv_objective(x, h) result(value)
    real(dp), intent(in) :: x(:), h
    real(dp) :: value, hh, delta, term, total
    integer :: i, j, n
    n = size(x)
    hh = h / 4.0_dp
    total = 0.0_dp
    do i = 1, n - 1
      do j = i + 1, n
        delta = ((x(i) - x(j)) / hh)**2
        if (delta >= 1000.0_dp) cycle
        term = exp(-delta / 4.0_dp) * &
          (delta * delta - 12.0_dp * delta + 12.0_dp)
        total = total + term
      end do
    end do
    value = 1.0_dp / (2.0_dp * real(n, dp) * hh * sqrt(pi_dp)) + &
      total / (64.0_dp * real(n, dp)**2 * hh * sqrt(pi_dp))
  end function bcv_objective

  function phi4_direct(x, h) result(value)
    real(dp), intent(in) :: x(:), h
    real(dp) :: value, delta, total
    integer :: i, j, n
    n = size(x)
    total = 0.0_dp
    do i = 1, n - 1
      do j = i + 1, n
        delta = ((x(i) - x(j)) / h)**2
        if (delta < 1000.0_dp) total = total + exp(-0.5_dp * delta) * &
          (delta * delta - 6.0_dp * delta + 3.0_dp)
      end do
    end do
    total = 2.0_dp * total + 3.0_dp * real(n, dp)
    value = total / (real(n, dp) * real(n - 1, dp) * h**5 * sqrt(2.0_dp * pi_dp))
  end function phi4_direct

  function phi6_direct(x, h) result(value)
    real(dp), intent(in) :: x(:), h
    real(dp) :: value, delta, total
    integer :: i, j, n
    n = size(x)
    total = 0.0_dp
    do i = 1, n - 1
      do j = i + 1, n
        delta = ((x(i) - x(j)) / h)**2
        if (delta < 1000.0_dp) total = total + exp(-0.5_dp * delta) * &
          (delta**3 - 15.0_dp * delta**2 + 45.0_dp * delta - 15.0_dp)
      end do
    end do
    total = 2.0_dp * total - 15.0_dp * real(n, dp)
    value = total / (real(n, dp) * real(n - 1, dp) * h**7 * sqrt(2.0_dp * pi_dp))
  end function phi6_direct

  function sj_equation(h, x, alpha2, c1) result(value)
    real(dp), intent(in) :: h, x(:), alpha2, c1
    real(dp) :: value, estimate
    estimate = phi4_direct(x, alpha2 * h**(5.0_dp / 7.0_dp))
    value = (c1 / estimate)**0.2_dp - h
  end function sj_equation

  subroutine golden_minimize(fun, a, b, tolerance, minimum)
    interface
      function fun(x) result(value)
        import dp
        real(dp), intent(in) :: x
        real(dp) :: value
      end function fun
    end interface
    real(dp), intent(in) :: a, b, tolerance
    real(dp), intent(out) :: minimum
    real(dp), parameter :: ratio = 0.6180339887498948482_dp
    real(dp) :: left, right, c, d, fc, fd
    integer :: iter
    left = a
    right = b
    c = right - ratio * (right - left)
    d = left + ratio * (right - left)
    fc = fun(c)
    fd = fun(d)
    do iter = 1, 200
      if (right - left <= tolerance) exit
      if (fc < fd) then
        right = d
        d = c
        fd = fc
        c = right - ratio * (right - left)
        fc = fun(c)
      else
        left = c
        c = d
        fc = fd
        d = left + ratio * (right - left)
        fd = fun(d)
      end if
    end do
    minimum = 0.5_dp * (left + right)
  end subroutine golden_minimize

end module mass_bandwidth
