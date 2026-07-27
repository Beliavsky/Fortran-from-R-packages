! SPDX-License-Identifier: GPL-2.0-or-later
! Derived from highfrequency 1.0.2 by Kris Boudt, Jonathan Cornelissen,
! Scott Payseur, Onno Kleen, Emil Sjoerup, and contributors.
module highfrequency_stats
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use highfrequency_kinds, only: dp, pi
  implicit none
  private
  public :: mean_value, sample_variance, median3, normal_cdf, normal_quantile
  public :: normal_pdf, mu_abs_normal, finite_vector, correlation_value
  public :: average_rank, clamp

contains

  pure real(dp) function mean_value(x) result(value)
    real(dp), intent(in) :: x(:)
    if (size(x) == 0) then
      value = 0.0_dp
    else
      value = sum(x) / real(size(x), dp)
    end if
  end function mean_value

  pure real(dp) function sample_variance(x) result(value)
    real(dp), intent(in) :: x(:)
    real(dp) :: m
    if (size(x) < 2) then
      value = 0.0_dp
      return
    end if
    m = mean_value(x)
    value = sum((x - m)**2) / real(size(x) - 1, dp)
  end function sample_variance

  pure real(dp) function correlation_value(x, y) result(value)
    real(dp), intent(in) :: x(:), y(:)
    real(dp) :: mx, my, sx, sy
    integer :: n
    n = min(size(x), size(y))
    if (n < 2) then
      value = 0.0_dp
      return
    end if
    mx = sum(x(:n)) / real(n, dp)
    my = sum(y(:n)) / real(n, dp)
    sx = sum((x(:n) - mx)**2)
    sy = sum((y(:n) - my)**2)
    if (sx <= 0.0_dp .or. sy <= 0.0_dp) then
      value = 0.0_dp
    else
      value = sum((x(:n) - mx) * (y(:n) - my)) / sqrt(sx * sy)
    end if
  end function correlation_value

  pure real(dp) function median3(a, b, c) result(value)
    real(dp), intent(in) :: a, b, c
    value = a + b + c - min(a, min(b, c)) - max(a, max(b, c))
  end function median3

  pure real(dp) function normal_pdf(x) result(value)
    real(dp), intent(in) :: x
    value = exp(-0.5_dp * x*x) / sqrt(2.0_dp*pi)
  end function normal_pdf

  pure real(dp) function normal_cdf(x) result(value)
    real(dp), intent(in) :: x
    value = 0.5_dp * erfc(-x / sqrt(2.0_dp))
  end function normal_cdf

  pure real(dp) function normal_quantile(p) result(x)
    real(dp), intent(in) :: p
    real(dp), parameter :: a(6) = [ &
      -3.969683028665376e+01_dp, 2.209460984245205e+02_dp, &
      -2.759285104469687e+02_dp, 1.383577518672690e+02_dp, &
      -3.066479806614716e+01_dp, 2.506628277459239e+00_dp ]
    real(dp), parameter :: b(5) = [ &
      -5.447609879822406e+01_dp, 1.615858368580409e+02_dp, &
      -1.556989798598866e+02_dp, 6.680131188771972e+01_dp, &
      -1.328068155288572e+01_dp ]
    real(dp), parameter :: c(6) = [ &
      -7.784894002430293e-03_dp, -3.223964580411365e-01_dp, &
      -2.400758277161838e+00_dp, -2.549732539343734e+00_dp, &
       4.374664141464968e+00_dp,  2.938163982698783e+00_dp ]
    real(dp), parameter :: d(4) = [ &
       7.784695709041462e-03_dp, 3.224671290700398e-01_dp, &
       2.445134137142996e+00_dp, 3.754408661907416e+00_dp ]
    real(dp), parameter :: plow = 0.02425_dp, phigh = 1.0_dp - plow
    real(dp) :: q, r, e
    if (p <= 0.0_dp) then
      x = -huge(1.0_dp)
      return
    else if (p >= 1.0_dp) then
      x = huge(1.0_dp)
      return
    end if
    if (p < plow) then
      q = sqrt(-2.0_dp * log(p))
      x = (((((c(1)*q+c(2))*q+c(3))*q+c(4))*q+c(5))*q+c(6)) / &
          ((((d(1)*q+d(2))*q+d(3))*q+d(4))*q+1.0_dp)
    else if (p <= phigh) then
      q = p - 0.5_dp
      r = q*q
      x = (((((a(1)*r+a(2))*r+a(3))*r+a(4))*r+a(5))*r+a(6))*q / &
          (((((b(1)*r+b(2))*r+b(3))*r+b(4))*r+b(5))*r+1.0_dp)
    else
      q = sqrt(-2.0_dp * log(1.0_dp-p))
      x = -(((((c(1)*q+c(2))*q+c(3))*q+c(4))*q+c(5))*q+c(6)) / &
           ((((d(1)*q+d(2))*q+d(3))*q+d(4))*q+1.0_dp)
    end if
    e = normal_cdf(x) - p
    x = x - e / max(normal_pdf(x), tiny(1.0_dp))
  end function normal_quantile

  pure real(dp) function mu_abs_normal(p) result(value)
    real(dp), intent(in) :: p
    value = 2.0_dp**(0.5_dp*p) * gamma(0.5_dp*(p+1.0_dp)) / sqrt(pi)
  end function mu_abs_normal

  pure logical function finite_vector(x) result(ok)
    real(dp), intent(in) :: x(:)
    integer :: i
    ok = .true.
    do i = 1, size(x)
      if (.not. ieee_is_finite(x(i))) then
        ok = .false.
        return
      end if
    end do
  end function finite_vector

  pure real(dp) function clamp(x, lower, upper) result(value)
    real(dp), intent(in) :: x, lower, upper
    value = min(upper, max(lower, x))
  end function clamp

  subroutine average_rank(x, ranks)
    real(dp), intent(in) :: x(:)
    real(dp), intent(out) :: ranks(size(x))
    integer, allocatable :: idx(:)
    integer :: i, j, k, n, tmp
    n = size(x)
    allocate(idx(n))
    idx = [(i, i=1,n)]
    do i = 2, n
      tmp = idx(i)
      j = i - 1
      do while (j >= 1)
        if (x(idx(j)) <= x(tmp)) exit
        idx(j+1) = idx(j)
        j = j - 1
      end do
      idx(j+1) = tmp
    end do
    i = 1
    do while (i <= n)
      j = i
      do while (j < n)
        if (abs(x(idx(j+1))-x(idx(i))) > 0.0_dp) exit
        j = j + 1
      end do
      do k = i, j
        ranks(idx(k)) = 0.5_dp * real(i+j, dp)
      end do
      i = j + 1
    end do
  end subroutine average_rank

end module highfrequency_stats
