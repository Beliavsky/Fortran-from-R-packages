! SPDX-License-Identifier: GPL-2.0-or-later
module infoset_stats
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use infoset_kinds, only : dp
  use infoset_status
  implicit none
  private
  public :: sort_real, median_real, quantile_real, normal_cdf, normal_logpdf
  public :: sample_covariance, nearest_positive_definite, left_histogram_risk
  public :: column_means
contains
  subroutine sort_real(x)
    real(dp), intent(inout) :: x(:)
    integer :: i, j
    real(dp) :: key
    do i = 2, size(x)
      key = x(i)
      j = i - 1
      do while (j >= 1)
        if (x(j) <= key) exit
        x(j + 1) = x(j)
        j = j - 1
      end do
      x(j + 1) = key
    end do
  end subroutine sort_real

  function median_real(x) result(value)
    real(dp), intent(in) :: x(:)
    real(dp) :: value
    real(dp), allocatable :: work(:)
    integer :: n
    n = size(x)
    if (n == 0) then
      value = 0.0_dp
      return
    end if
    allocate(work(n))
    work = x
    call sort_real(work)
    if (mod(n, 2) == 1) then
      value = work((n + 1) / 2)
    else
      value = 0.5_dp * (work(n / 2) + work(n / 2 + 1))
    end if
  end function median_real

  function quantile_real(x, probability) result(value)
    real(dp), intent(in) :: x(:)
    real(dp), intent(in) :: probability
    real(dp) :: value, h, fraction
    real(dp), allocatable :: work(:)
    integer :: n, lower
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
      lower = int(floor(h))
      fraction = h - real(lower, dp)
      if (lower >= n) then
        value = work(n)
      else
        value = work(lower) + fraction * (work(lower + 1) - work(lower))
      end if
    end if
  end function quantile_real

  elemental function normal_cdf(x, mean, sd) result(value)
    real(dp), intent(in) :: x, mean, sd
    real(dp) :: value
    if (sd <= 0.0_dp) then
      value = merge(0.0_dp, 1.0_dp, x < mean)
    else
      value = 0.5_dp * erfc(-(x - mean) / (sqrt(2.0_dp) * sd))
    end if
  end function normal_cdf

  elemental function normal_logpdf(x, mean, sd) result(value)
    real(dp), intent(in) :: x, mean, sd
    real(dp) :: value
    if (sd <= 0.0_dp) then
      value = -huge(1.0_dp)
    else
      value = -0.5_dp * log(2.0_dp * acos(-1.0_dp)) - log(sd) &
        - 0.5_dp * ((x - mean) / sd)**2
    end if
  end function normal_logpdf

  function column_means(x) result(means)
    real(dp), intent(in) :: x(:,:)
    real(dp), allocatable :: means(:)
    integer :: n
    n = size(x, 1)
    allocate(means(size(x, 2)))
    if (n > 0) then
      means = sum(x, dim=1) / real(n, dp)
    else
      means = 0.0_dp
    end if
  end function column_means

  subroutine sample_covariance(x, covariance, status)
    real(dp), intent(in) :: x(:,:)
    real(dp), allocatable, intent(out) :: covariance(:,:)
    integer, intent(out) :: status
    real(dp), allocatable :: means(:), centered(:,:)
    integer :: n, p
    n = size(x, 1)
    p = size(x, 2)
    allocate(covariance(p, p))
    covariance = 0.0_dp
    if (n < 2 .or. p < 1 .or. .not. all(ieee_is_finite(x))) then
      status = infoset_invalid_argument
      return
    end if
    means = column_means(x)
    allocate(centered(n, p))
    centered = x - spread(means, 1, n)
    covariance = matmul(transpose(centered), centered) / real(n - 1, dp)
    status = infoset_success
  end subroutine sample_covariance

  subroutine jacobi_eigen(a, eigenvalues, eigenvectors, status)
    real(dp), intent(in) :: a(:,:)
    real(dp), allocatable, intent(out) :: eigenvalues(:), eigenvectors(:,:)
    integer, intent(out) :: status
    real(dp), allocatable :: b(:,:)
    real(dp) :: app, aqq, apq, tau, tangent, cosine, sine
    real(dp) :: bip, biq, maximum
    integer :: n, i, j, p, q, iteration
    n = size(a, 1)
    if (n < 1 .or. size(a, 2) /= n) then
      allocate(eigenvalues(0), eigenvectors(0, 0))
      status = infoset_invalid_argument
      return
    end if
    allocate(b(n, n), eigenvalues(n), eigenvectors(n, n))
    b = 0.5_dp * (a + transpose(a))
    eigenvectors = 0.0_dp
    do i = 1, n
      eigenvectors(i, i) = 1.0_dp
    end do
    if (n == 1) then
      eigenvalues(1) = b(1, 1)
      status = infoset_success
      return
    end if
    do iteration = 1, max(100, 100 * n * n)
      maximum = 0.0_dp
      p = 1
      q = 2
      do i = 1, n - 1
        do j = i + 1, n
          if (abs(b(i, j)) > maximum) then
            maximum = abs(b(i, j))
            p = i
            q = j
          end if
        end do
      end do
      if (maximum <= 1.0e-13_dp * (1.0_dp + maxval(abs(b)))) exit
      app = b(p, p)
      aqq = b(q, q)
      apq = b(p, q)
      tau = (aqq - app) / (2.0_dp * apq)
      if (tau >= 0.0_dp) then
        tangent = 1.0_dp / (tau + sqrt(1.0_dp + tau * tau))
      else
        tangent = -1.0_dp / (-tau + sqrt(1.0_dp + tau * tau))
      end if
      cosine = 1.0_dp / sqrt(1.0_dp + tangent * tangent)
      sine = tangent * cosine
      do i = 1, n
        if (i /= p .and. i /= q) then
          bip = b(i, p)
          biq = b(i, q)
          b(i, p) = cosine * bip - sine * biq
          b(p, i) = b(i, p)
          b(i, q) = sine * bip + cosine * biq
          b(q, i) = b(i, q)
        end if
      end do
      b(p, p) = cosine * cosine * app - 2.0_dp * sine * cosine * apq &
        + sine * sine * aqq
      b(q, q) = sine * sine * app + 2.0_dp * sine * cosine * apq &
        + cosine * cosine * aqq
      b(p, q) = 0.0_dp
      b(q, p) = 0.0_dp
      do i = 1, n
        bip = eigenvectors(i, p)
        biq = eigenvectors(i, q)
        eigenvectors(i, p) = cosine * bip - sine * biq
        eigenvectors(i, q) = sine * bip + cosine * biq
      end do
    end do
    do i = 1, n
      eigenvalues(i) = b(i, i)
    end do
    status = infoset_success
  end subroutine jacobi_eigen

  subroutine nearest_positive_definite(a, result, status, relative_floor)
    real(dp), intent(in) :: a(:,:)
    real(dp), allocatable, intent(out) :: result(:,:)
    integer, intent(out) :: status
    real(dp), intent(in), optional :: relative_floor
    real(dp), allocatable :: eigenvalues(:), eigenvectors(:,:), diagonal(:,:)
    real(dp) :: floor_value, scale, rel
    integer :: n, i
    n = size(a, 1)
    if (n < 1 .or. size(a, 2) /= n .or. .not. all(ieee_is_finite(a))) then
      allocate(result(0, 0))
      status = infoset_invalid_argument
      return
    end if
    rel = 1.0e-8_dp
    if (present(relative_floor)) rel = max(relative_floor, epsilon(1.0_dp))
    call jacobi_eigen(a, eigenvalues, eigenvectors, status)
    if (status /= infoset_success) then
      allocate(result(0, 0))
      return
    end if
    scale = max(1.0_dp, maxval(abs(eigenvalues)))
    floor_value = rel * scale
    eigenvalues = max(eigenvalues, floor_value)
    allocate(diagonal(n, n))
    diagonal = 0.0_dp
    do i = 1, n
      diagonal(i, i) = eigenvalues(i)
    end do
    result = matmul(eigenvectors, matmul(diagonal, transpose(eigenvectors)))
    result = 0.5_dp * (result + transpose(result))
    status = infoset_success
  end subroutine nearest_positive_definite

  subroutine left_histogram_risk(log_returns, log_change_point, risk, status)
    real(dp), intent(in) :: log_returns(:)
    real(dp), intent(in) :: log_change_point
    real(dp), intent(out) :: risk
    integer, intent(out) :: status
    integer, allocatable :: counts(:)
    integer :: n, bins, i, index
    real(dp) :: xmin, xmax, width, lower, upper, fraction_upper
    real(dp) :: density, probability, first_moment
    n = size(log_returns)
    risk = 0.0_dp
    if (n < 2 .or. .not. all(ieee_is_finite(log_returns)) .or. &
        .not. ieee_is_finite(log_change_point)) then
      status = infoset_invalid_argument
      return
    end if
    xmin = minval(log_returns)
    xmax = maxval(log_returns)
    if (log_change_point <= xmin) then
      status = infoset_insufficient_data
      return
    end if
    if (xmax <= xmin + epsilon(1.0_dp)) then
      risk = -xmin
      status = infoset_success
      return
    end if
    bins = max(2, ceiling(log(real(n, dp)) / log(2.0_dp) + 1.0_dp))
    width = (xmax - xmin) / real(bins, dp)
    allocate(counts(bins))
    counts = 0
    do i = 1, n
      index = int((log_returns(i) - xmin) / width) + 1
      index = max(1, min(bins, index))
      counts(index) = counts(index) + 1
    end do
    probability = 0.0_dp
    first_moment = 0.0_dp
    do i = 1, bins
      lower = xmin + real(i - 1, dp) * width
      upper = lower + width
      if (log_change_point <= lower) exit
      fraction_upper = min(upper, log_change_point)
      density = real(counts(i), dp) / (real(n, dp) * width)
      probability = probability + density * (fraction_upper - lower)
      first_moment = first_moment + 0.5_dp * density &
        * (fraction_upper * fraction_upper - lower * lower)
      if (log_change_point < upper) exit
    end do
    if (probability <= sqrt(epsilon(1.0_dp))) then
      status = infoset_insufficient_data
    else
      risk = -first_moment / probability
      status = infoset_success
    end if
  end subroutine left_histogram_risk
end module infoset_stats
