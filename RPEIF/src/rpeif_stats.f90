! SPDX-License-Identifier: GPL-3.0-or-later
module rpeif_stats
  use rpeif_kinds, only : dp, pi
  implicit none
  private
  public :: mean_value, sample_sd, population_sd, median_value, mad_value
  public :: quantile_type7, normal_pdf, normal_cdf, normal_quantile
  public :: gaussian_kde_density, lower_partial_moment, upper_partial_moment
  public :: lower_string
contains
  pure function mean_value(x) result(value)
    real(dp), intent(in) :: x(:)
    real(dp) :: value

    if (size(x) == 0) then
      value = 0.0_dp
    else
      value = sum(x) / real(size(x), dp)
    end if
  end function mean_value

  function sample_sd(x) result(value)
    real(dp), intent(in) :: x(:)
    real(dp) :: value, mu
    integer :: n

    n = size(x)
    if (n <= 1) then
      value = 0.0_dp
      return
    end if
    mu = mean_value(x)
    value = sqrt(max(0.0_dp, sum((x - mu) ** 2) / real(n - 1, dp)))
  end function sample_sd

  function population_sd(x) result(value)
    real(dp), intent(in) :: x(:)
    real(dp) :: value, mu
    integer :: n

    n = size(x)
    if (n == 0) then
      value = 0.0_dp
      return
    end if
    mu = mean_value(x)
    value = sqrt(max(0.0_dp, sum((x - mu) ** 2) / real(n, dp)))
  end function population_sd

  function median_value(x) result(value)
    real(dp), intent(in) :: x(:)
    real(dp) :: value

    value = quantile_type7(x, 0.5_dp)
  end function median_value

  function mad_value(x) result(value)
    real(dp), intent(in) :: x(:)
    real(dp) :: value, center
    real(dp), allocatable :: deviations(:)

    if (size(x) == 0) then
      value = 0.0_dp
      return
    end if
    center = median_value(x)
    allocate(deviations(size(x)))
    deviations = abs(x - center)
    value = median_value(deviations) / 0.6744897501960817_dp
  end function mad_value

  function quantile_type7(x, probability) result(value)
    real(dp), intent(in) :: x(:), probability
    real(dp) :: value, h, fraction
    real(dp), allocatable :: work(:)
    integer :: n, lower_index

    n = size(x)
    if (n == 0) then
      value = 0.0_dp
      return
    end if
    allocate(work(n))
    work = x
    call sort_real(work)
    if (probability <= 0.0_dp) then
      value = work(1)
    else if (probability >= 1.0_dp) then
      value = work(n)
    else
      h = 1.0_dp + real(n - 1, dp) * probability
      lower_index = floor(h)
      fraction = h - real(lower_index, dp)
      if (lower_index >= n) then
        value = work(n)
      else
        value = (1.0_dp - fraction) * work(lower_index) + fraction * work(lower_index + 1)
      end if
    end if
  end function quantile_type7

  elemental function normal_pdf(x, mu, sigma) result(value)
    real(dp), intent(in) :: x
    real(dp), intent(in), optional :: mu, sigma
    real(dp) :: value, center, scale, z

    center = 0.0_dp
    if (present(mu)) center = mu
    scale = 1.0_dp
    if (present(sigma)) scale = sigma
    if (scale <= 0.0_dp) then
      value = 0.0_dp
      return
    end if
    z = (x - center) / scale
    value = exp(-0.5_dp * z * z) / (sqrt(2.0_dp * pi) * scale)
  end function normal_pdf

  elemental function normal_cdf(x) result(value)
    real(dp), intent(in) :: x
    real(dp) :: value

    value = 0.5_dp * erfc(-x / sqrt(2.0_dp))
  end function normal_cdf

  function normal_quantile(probability) result(value)
    real(dp), intent(in) :: probability
    real(dp) :: value, q, r, error
    real(dp), parameter :: a1 = -3.969683028665376e+01_dp
    real(dp), parameter :: a2 =  2.209460984245205e+02_dp
    real(dp), parameter :: a3 = -2.759285104469687e+02_dp
    real(dp), parameter :: a4 =  1.383577518672690e+02_dp
    real(dp), parameter :: a5 = -3.066479806614716e+01_dp
    real(dp), parameter :: a6 =  2.506628277459239e+00_dp
    real(dp), parameter :: b1 = -5.447609879822406e+01_dp
    real(dp), parameter :: b2 =  1.615858368580409e+02_dp
    real(dp), parameter :: b3 = -1.556989798598866e+02_dp
    real(dp), parameter :: b4 =  6.680131188771972e+01_dp
    real(dp), parameter :: b5 = -1.328068155288572e+01_dp
    real(dp), parameter :: c1 = -7.784894002430293e-03_dp
    real(dp), parameter :: c2 = -3.223964580411365e-01_dp
    real(dp), parameter :: c3 = -2.400758277161838e+00_dp
    real(dp), parameter :: c4 = -2.549732539343734e+00_dp
    real(dp), parameter :: c5 =  4.374664141464968e+00_dp
    real(dp), parameter :: c6 =  2.938163982698783e+00_dp
    real(dp), parameter :: d1 =  7.784695709041462e-03_dp
    real(dp), parameter :: d2 =  3.224671290700398e-01_dp
    real(dp), parameter :: d3 =  2.445134137142996e+00_dp
    real(dp), parameter :: d4 =  3.754408661907416e+00_dp
    real(dp), parameter :: p_low = 0.02425_dp
    real(dp), parameter :: p_high = 1.0_dp - p_low

    if (probability <= 0.0_dp) then
      value = -huge(1.0_dp)
      return
    else if (probability >= 1.0_dp) then
      value = huge(1.0_dp)
      return
    end if

    if (probability < p_low) then
      q = sqrt(-2.0_dp * log(probability))
      value = (((((c1 * q + c2) * q + c3) * q + c4) * q + c5) * q + c6) / &
        ((((d1 * q + d2) * q + d3) * q + d4) * q + 1.0_dp)
    else if (probability <= p_high) then
      q = probability - 0.5_dp
      r = q * q
      value = (((((a1 * r + a2) * r + a3) * r + a4) * r + a5) * r + a6) * q / &
        (((((b1 * r + b2) * r + b3) * r + b4) * r + b5) * r + 1.0_dp)
    else
      q = sqrt(-2.0_dp * log(1.0_dp - probability))
      value = -(((((c1 * q + c2) * q + c3) * q + c4) * q + c5) * q + c6) / &
        ((((d1 * q + d2) * q + d3) * q + d4) * q + 1.0_dp)
    end if

    error = normal_cdf(value) - probability
    value = value - error / max(normal_pdf(value), tiny(1.0_dp))
  end function normal_quantile

  function gaussian_kde_density(x, point) result(value)
    real(dp), intent(in) :: x(:), point
    real(dp) :: value, bandwidth, hi, lo, iqr
    integer :: n

    n = size(x)
    if (n == 0) then
      value = 0.0_dp
      return
    end if
    hi = sample_sd(x)
    iqr = quantile_type7(x, 0.75_dp) - quantile_type7(x, 0.25_dp)
    lo = min(hi, iqr / 1.34_dp)
    if (lo <= 0.0_dp) lo = hi
    if (lo <= 0.0_dp) lo = abs(x(1))
    if (lo <= 0.0_dp) lo = 1.0_dp
    bandwidth = 0.9_dp * lo * real(n, dp) ** (-0.2_dp)
    value = sum(normal_pdf((point - x) / bandwidth)) / (real(n, dp) * bandwidth)
  end function gaussian_kde_density

  function lower_partial_moment(x, threshold, order) result(value)
    real(dp), intent(in) :: x(:), threshold
    integer, intent(in) :: order
    real(dp) :: value
    integer :: i

    value = 0.0_dp
    if (size(x) == 0) return
    do i = 1, size(x)
      if (x(i) <= threshold) value = value + (threshold - x(i)) ** order
    end do
    value = value / real(size(x), dp)
  end function lower_partial_moment

  function upper_partial_moment(x, threshold, order, source_compatibility) result(value)
    real(dp), intent(in) :: x(:), threshold
    integer, intent(in) :: order
    logical, intent(in), optional :: source_compatibility
    real(dp) :: value
    integer :: i
    logical :: compatible

    compatible = .false.
    if (present(source_compatibility)) compatible = source_compatibility
    value = 0.0_dp
    if (size(x) == 0) return
    do i = 1, size(x)
      if (x(i) >= threshold) then
        if (compatible) then
          value = value + (threshold - x(i)) ** order
        else
          value = value + (x(i) - threshold) ** order
        end if
      end if
    end do
    value = value / real(size(x), dp)
  end function upper_partial_moment

  pure function lower_string(text) result(out)
    character(len=*), intent(in) :: text
    character(len=len(text)) :: out
    integer :: i, code

    out = text
    do i = 1, len(text)
      code = iachar(out(i:i))
      if (code >= iachar('A') .and. code <= iachar('Z')) out(i:i) = achar(code + 32)
    end do
  end function lower_string

  recursive subroutine quicksort(a, left, right)
    real(dp), intent(inout) :: a(:)
    integer, intent(in) :: left, right
    integer :: i, j
    real(dp) :: pivot, tmp

    if (left >= right) return
    i = left
    j = right
    pivot = a((left + right) / 2)
    do
      do while (a(i) < pivot)
        i = i + 1
      end do
      do while (a(j) > pivot)
        j = j - 1
      end do
      if (i <= j) then
        tmp = a(i)
        a(i) = a(j)
        a(j) = tmp
        i = i + 1
        j = j - 1
      end if
      if (i > j) exit
    end do
    if (left < j) call quicksort(a, left, j)
    if (i < right) call quicksort(a, i, right)
  end subroutine quicksort

  subroutine sort_real(a)
    real(dp), intent(inout) :: a(:)
    if (size(a) > 1) call quicksort(a, 1, size(a))
  end subroutine sort_real
end module rpeif_stats
