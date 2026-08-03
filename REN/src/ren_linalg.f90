! SPDX-License-Identifier: AGPL-3.0-or-later
! Derived from REN 0.1.0 computational code; see NOTICE.md.
module ren_linalg
  use ren_kinds, only : dp
  use ren_types, only : ren_success, ren_invalid_argument, ren_dimension_error, ren_numerical_error
  use corpcor, only : pseudoinverse, make_positive_definite
  implicit none
  private
  public :: sample_variance, column_variances, covariance_matrix, correlation_matrix
  public :: equality_minimum_variance, long_only_minimum_variance
  public :: least_squares, vector_sd, l1_norm, cumulative_metrics
contains
  pure function sample_variance(x) result(value)
    real(dp), intent(in) :: x(:)
    real(dp) :: value, mean_x
    integer :: n
    n = size(x)
    if (n <= 1) then
      value = 0.0_dp
    else
      mean_x = sum(x) / real(n, dp)
      value = sum((x - mean_x) ** 2) / real(n - 1, dp)
    end if
  end function sample_variance

  pure function vector_sd(x) result(value)
    real(dp), intent(in) :: x(:)
    real(dp) :: value
    value = sqrt(max(0.0_dp, sample_variance(x)))
  end function vector_sd

  pure function column_variances(x) result(value)
    real(dp), intent(in) :: x(:, :)
    real(dp) :: value(size(x, 2))
    integer :: j
    do j = 1, size(x, 2)
      value(j) = sample_variance(x(:, j))
    end do
  end function column_variances

  function covariance_matrix(x, status) result(value)
    real(dp), intent(in) :: x(:, :)
    integer, intent(out), optional :: status
    real(dp), allocatable :: value(:, :)
    real(dp), allocatable :: centered(:, :), mean_x(:)
    integer :: n, p
    n = size(x, 1)
    p = size(x, 2)
    allocate(value(p, p))
    value = 0.0_dp
    if (n <= 1 .or. p < 1) then
      if (present(status)) status = ren_invalid_argument
      return
    end if
    mean_x = sum(x, dim=1) / real(n, dp)
    allocate(centered(n, p))
    centered = x - spread(mean_x, 1, n)
    value = matmul(transpose(centered), centered) / real(n - 1, dp)
    value = 0.5_dp * (value + transpose(value))
    if (present(status)) status = ren_success
  end function covariance_matrix

  function correlation_matrix(x, status) result(value)
    real(dp), intent(in) :: x(:, :)
    integer, intent(out), optional :: status
    real(dp), allocatable :: value(:, :), covariance(:, :), sd(:)
    integer :: i, j, istat, p
    covariance = covariance_matrix(x, istat)
    p = size(covariance, 1)
    allocate(value(p, p), sd(p))
    value = 0.0_dp
    if (istat /= ren_success) then
      if (present(status)) status = istat
      return
    end if
    do i = 1, p
      sd(i) = sqrt(max(0.0_dp, covariance(i, i)))
    end do
    do j = 1, p
      do i = 1, p
        if (sd(i) > 0.0_dp .and. sd(j) > 0.0_dp) value(i, j) = covariance(i, j) / (sd(i) * sd(j))
      end do
      value(j, j) = 1.0_dp
    end do
    value = max(-1.0_dp, min(1.0_dp, value))
    if (present(status)) status = ren_success
  end function correlation_matrix

  subroutine equality_minimum_variance(covariance, weights, status)
    real(dp), intent(in) :: covariance(:, :)
    real(dp), allocatable, intent(out) :: weights(:)
    integer, intent(out) :: status
    real(dp), allocatable :: positive(:, :), inverse(:, :), ones(:), work(:)
    real(dp) :: denominator
    integer :: p, istat
    p = size(covariance, 1)
    allocate(weights(p))
    if (p < 1 .or. size(covariance, 2) /= p) then
      weights = 0.0_dp
      status = ren_dimension_error
      return
    end if
    positive = make_positive_definite(0.5_dp * (covariance + transpose(covariance)), status=istat)
    inverse = pseudoinverse(positive, status=istat)
    allocate(ones(p))
    ones = 1.0_dp
    work = matmul(inverse, ones)
    denominator = sum(work)
    if (abs(denominator) <= 100.0_dp * epsilon(1.0_dp)) then
      weights = 1.0_dp / real(p, dp)
      status = ren_numerical_error
      return
    end if
    weights = work / denominator
    status = ren_success
  end subroutine equality_minimum_variance

  subroutine long_only_minimum_variance(covariance, weights, status, tolerance, max_iterations)
    real(dp), intent(in) :: covariance(:, :)
    real(dp), allocatable, intent(out) :: weights(:)
    integer, intent(out) :: status
    real(dp), intent(in), optional :: tolerance
    integer, intent(in), optional :: max_iterations
    real(dp), allocatable :: q(:, :), current(:), extrapolated(:), next(:), gradient(:)
    real(dp) :: tol, lipschitz, t_old, t_new, error
    integer :: p, iteration, maxit, istat
    p = size(covariance, 1)
    allocate(weights(p))
    if (p < 1 .or. size(covariance, 2) /= p) then
      weights = 0.0_dp
      status = ren_dimension_error
      return
    end if
    q = make_positive_definite(0.5_dp * (covariance + transpose(covariance)), status=istat)
    lipschitz = maxval(sum(abs(q), dim=2))
    if (lipschitz <= tiny(1.0_dp)) then
      weights = 1.0_dp / real(p, dp)
      status = ren_success
      return
    end if
    tol = 1.0e-11_dp
    if (present(tolerance)) tol = max(tolerance, epsilon(1.0_dp))
    maxit = 50000
    if (present(max_iterations)) maxit = max(1, max_iterations)
    allocate(current(p), extrapolated(p), next(p), gradient(p))
    current = 1.0_dp / real(p, dp)
    extrapolated = current
    t_old = 1.0_dp
    do iteration = 1, maxit
      gradient = matmul(q, extrapolated)
      call project_simplex(extrapolated - gradient / lipschitz, next)
      error = sqrt(sum((next - current) ** 2))
      if (error <= tol * max(1.0_dp, sqrt(sum(current ** 2)))) exit
      t_new = 0.5_dp * (1.0_dp + sqrt(1.0_dp + 4.0_dp * t_old * t_old))
      extrapolated = next + ((t_old - 1.0_dp) / t_new) * (next - current)
      current = next
      t_old = t_new
    end do
    weights = next
    where (weights < 10.0_dp * epsilon(1.0_dp)) weights = 0.0_dp
    if (sum(weights) > 0.0_dp) weights = weights / sum(weights)
    if (iteration > maxit) then
      status = ren_numerical_error
    else
      status = ren_success
    end if
  end subroutine long_only_minimum_variance

  subroutine project_simplex(x, projected)
    real(dp), intent(in) :: x(:)
    real(dp), intent(out) :: projected(size(x))
    real(dp), allocatable :: sorted(:)
    real(dp) :: cumulative, theta
    integer :: i, rho, n
    n = size(x)
    sorted = x
    call sort_descending(sorted)
    cumulative = 0.0_dp
    rho = 1
    do i = 1, n
      cumulative = cumulative + sorted(i)
      if (sorted(i) - (cumulative - 1.0_dp) / real(i, dp) > 0.0_dp) rho = i
    end do
    theta = (sum(sorted(1:rho)) - 1.0_dp) / real(rho, dp)
    projected = max(x - theta, 0.0_dp)
  end subroutine project_simplex

  subroutine sort_descending(x)
    real(dp), intent(inout) :: x(:)
    real(dp) :: key
    integer :: i, j
    do i = 2, size(x)
      key = x(i)
      j = i - 1
      do while (j >= 1)
        if (x(j) >= key) exit
        x(j + 1) = x(j)
        j = j - 1
      end do
      x(j + 1) = key
    end do
  end subroutine sort_descending

  subroutine least_squares(x, y, beta, status)
    real(dp), intent(in) :: x(:, :), y(:)
    real(dp), allocatable, intent(out) :: beta(:)
    integer, intent(out) :: status
    real(dp), allocatable :: inverse(:, :)
    integer :: istat
    if (size(x, 1) /= size(y)) then
      allocate(beta(0))
      status = ren_dimension_error
      return
    end if
    if (size(x, 2) == 0) then
      allocate(beta(0))
      status = ren_success
      return
    end if
    inverse = pseudoinverse(x, status=istat)
    beta = matmul(inverse, y)
    if (istat == 0) then
      status = ren_success
    else
      status = ren_numerical_error
    end if
  end subroutine least_squares

  pure function l1_norm(x) result(value)
    real(dp), intent(in) :: x(:)
    real(dp) :: value
    value = sum(abs(x))
  end function l1_norm

  subroutine cumulative_metrics(gross_returns, cumulative_return, wealth, sharpe, volatility, max_drawdown)
    real(dp), intent(in) :: gross_returns(:)
    real(dp), intent(out) :: cumulative_return(size(gross_returns)), wealth(size(gross_returns))
    real(dp), intent(out) :: sharpe, volatility, max_drawdown
    real(dp), allocatable :: net(:)
    real(dp) :: peak, drawdown, sd
    integer :: i, n
    n = size(gross_returns)
    if (n == 0) then
      sharpe = 0.0_dp
      volatility = 0.0_dp
      max_drawdown = 0.0_dp
      return
    end if
    allocate(net(n))
    net = gross_returns - 1.0_dp
    cumulative_return(1) = net(1)
    wealth(1) = gross_returns(1)
    do i = 2, n
      cumulative_return(i) = cumulative_return(i - 1) + net(i)
      wealth(i) = wealth(i - 1) * gross_returns(i)
    end do
    sd = vector_sd(net)
    if (sd > 0.0_dp) then
      sharpe = sqrt(252.0_dp) * (sum(net) / real(n, dp)) / sd
    else
      sharpe = 0.0_dp
    end if
    volatility = sqrt(252.0_dp) * sd * 100.0_dp
    peak = wealth(1)
    max_drawdown = 0.0_dp
    do i = 1, n
      peak = max(peak, wealth(i))
      if (peak > 0.0_dp) then
        drawdown = (peak - wealth(i)) / peak * 100.0_dp
        max_drawdown = max(max_drawdown, drawdown)
      end if
    end do
  end subroutine cumulative_metrics
end module ren_linalg
