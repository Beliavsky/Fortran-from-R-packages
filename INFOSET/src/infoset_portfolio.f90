! SPDX-License-Identifier: GPL-2.0-or-later
module infoset_portfolio
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use infoset_kinds, only : dp
  use infoset_status
  use infoset_types
  use infoset_core, only : create_overlapping_windows
  use infoset_stats, only : column_means, sample_covariance
  use infoset_stats, only : nearest_positive_definite, quantile_real, sort_real
  use quadprog, only : qp_result, qp_success, solve_qp
  implicit none
  private
  public :: ptf_construction, summary_ptf
contains
  function uppercase(text) result(value)
    character(len=*), intent(in) :: text
    character(len=len(text)) :: value
    integer :: i, code
    value = text
    do i = 1, len(text)
      code = iachar(value(i:i))
      if (code >= iachar('a') .and. code <= iachar('z')) then
        value(i:i) = achar(code - iachar('a') + iachar('A'))
      end if
    end do
  end function uppercase

  function log_return_matrix(prices) result(returns)
    real(dp), intent(in) :: prices(:,:)
    real(dp), allocatable :: returns(:,:)
    integer :: n
    n = size(prices, 1)
    allocate(returns(max(0, n - 1), size(prices, 2)))
    if (n >= 2) returns = log(prices(2:n, :) / prices(1:n - 1, :))
  end function log_return_matrix

  subroutine extreme_downside_covariance(returns, covariance, status)
    real(dp), intent(in) :: returns(:,:)
    real(dp), allocatable, intent(out) :: covariance(:,:)
    integer, intent(out) :: status
    real(dp), allocatable :: means(:), downside(:,:)
    real(dp) :: threshold
    integer :: n, p, i, j
    n = size(returns, 1)
    p = size(returns, 2)
    allocate(covariance(p, p))
    covariance = 0.0_dp
    if (n < 2 .or. p < 1) then
      status = infoset_invalid_argument
      return
    end if
    means = column_means(returns)
    allocate(downside(n, p))
    downside = returns - spread(means, 1, n)
    do j = 1, p
      threshold = quantile_real(returns(:, j), 0.05_dp)
      do i = 1, n
        if (returns(i, j) > threshold) downside(i, j) = 0.0_dp
      end do
    end do
    covariance = matmul(transpose(downside), downside) / real(n, dp)
    status = infoset_success
  end subroutine extreme_downside_covariance

  subroutine build_constraints(means, constraints, bounds)
    real(dp), intent(in) :: means(:)
    real(dp), allocatable, intent(out) :: constraints(:,:), bounds(:)
    integer :: n, i
    n = size(means)
    allocate(constraints(n, 2 + 2 * n), bounds(2 + 2 * n))
    constraints = 0.0_dp
    bounds = 0.0_dp
    constraints(:, 1) = 1.0_dp
    constraints(:, 2) = means
    bounds(1) = 1.0_dp
    bounds(2) = sum(means) / real(n, dp)
    do i = 1, n
      constraints(i, 2 + i) = 1.0_dp
      constraints(i, 2 + n + i) = -1.0_dp
      bounds(2 + n + i) = -1.0_dp
    end do
  end subroutine build_constraints

  subroutine ptf_construction(data, window_size, overlap, strategy, result, &
      left_risk, penalty)
    real(dp), intent(in) :: data(:,:)
    integer, intent(in) :: window_size, overlap
    character(len=*), intent(in) :: strategy
    type(portfolio_result), intent(out) :: result
    real(dp), intent(in), optional :: left_risk(:,:)
    real(dp), intent(in), optional :: penalty
    type(window_collection) :: windows
    type(qp_result) :: fit
    real(dp), allocatable :: returns(:,:), tail_returns(:,:), means(:)
    real(dp), allocatable :: covariance(:,:), positive_covariance(:,:)
    real(dp), allocatable :: constraints(:,:), bounds(:), dvec(:), weights(:)
    real(dp) :: lambda_value, portfolio_value
    character(len=:), allocatable :: method
    integer :: number, assets, t, j, start_tail, out_start, out_length, status

    result = portfolio_result()
    method = trim(uppercase(strategy))
    result%strategy = method
    lambda_value = 1.0e-4_dp
    if (present(penalty)) lambda_value = penalty
    windows = create_overlapping_windows(data, window_size, overlap)
    if (windows%status /= infoset_success) then
      allocate(result%weights(0, 0), result%oos_returns(0, 0), result%baseline_value(0))
      result%status = windows%status
      return
    end if
    number = size(windows%values, 3)
    assets = size(data, 2)
    if (number < 2 .or. window_size < 4 .or. overlap >= window_size - 1 &
        .or. lambda_value < 0.0_dp) then
      allocate(result%weights(0, 0), result%oos_returns(0, 0), result%baseline_value(0))
      result%status = infoset_invalid_argument
      return
    end if
    if ((method == 'C_M' .or. method == 'C_EDC') .and. .not. present(left_risk)) then
      allocate(result%weights(0, 0), result%oos_returns(0, 0), result%baseline_value(0))
      result%status = infoset_invalid_argument
      return
    end if
    if (present(left_risk)) then
      if (size(left_risk, 1) /= assets .or. size(left_risk, 2) < number - 1) then
        allocate(result%weights(0, 0), result%oos_returns(0, 0), result%baseline_value(0))
        result%status = infoset_invalid_argument
        return
      end if
    end if
    if (method /= 'M' .and. method /= 'C_M' .and. method /= 'EDC' &
        .and. method /= 'C_EDC') then
      allocate(result%weights(0, 0), result%oos_returns(0, 0), result%baseline_value(0))
      result%status = infoset_invalid_argument
      return
    end if

    out_length = overlap + 1
    allocate(result%weights(assets, number - 1))
    allocate(result%oos_returns(out_length, number - 1))
    allocate(result%baseline_value(number - 1))
    result%weights = 0.0_dp
    result%oos_returns = 0.0_dp
    result%baseline_value = 0.0_dp
    result%status = infoset_success

    do t = 1, number - 1
      returns = log_return_matrix(windows%values(:, :, t))
      means = column_means(returns)
      start_tail = max(1, size(returns, 1) - 2 * (overlap + 1))
      tail_returns = returns(start_tail:size(returns, 1), :)
      if (method == 'EDC' .or. method == 'C_EDC') then
        call extreme_downside_covariance(tail_returns, covariance, status)
      else
        call sample_covariance(tail_returns, covariance, status)
      end if
      if (status /= infoset_success) then
        result%status = status
        return
      end if
      call nearest_positive_definite(covariance, positive_covariance, status)
      if (status /= infoset_success) then
        result%status = status
        return
      end if
      call build_constraints(means, constraints, bounds)
      allocate(dvec(assets))
      dvec = 0.0_dp
      if (method == 'C_M') dvec = -lambda_value * left_risk(:, t)
      if (method == 'C_EDC') dvec = lambda_value * left_risk(:, t)
      fit = solve_qp(positive_covariance, dvec, constraints, bounds, meq=2)
      if (fit%status /= qp_success) then
        result%status = infoset_qp_error
        return
      end if
      weights = fit%solution
      result%weights(:, t) = weights
      portfolio_value = dot_product(weights, windows%values(window_size - 1, :, t))
      result%baseline_value(t) = portfolio_value
      out_start = window_size - 1 - overlap
      do j = 1, out_length
        result%oos_returns(j, t) = (dot_product(weights, &
          windows%values(out_start + j - 1, :, t + 1)) - portfolio_value) &
          / portfolio_value
      end do
      deallocate(dvec)
    end do
  end subroutine ptf_construction

  subroutine summary_ptf(oos_returns, result)
    real(dp), intent(in) :: oos_returns(:,:)
    type(portfolio_summary), intent(out) :: result
    real(dp), allocatable :: values(:)
    integer :: n
    result = portfolio_summary()
    n = size(oos_returns)
    if (n == 0 .or. .not. all(ieee_is_finite(oos_returns))) then
      result%status = infoset_invalid_argument
      return
    end if
    allocate(values(n))
    values = reshape(oos_returns, [n])
    call sort_real(values)
    result%count = n
    result%minimum = values(1)
    result%first_quartile = quantile_real(values, 0.25_dp)
    result%median = quantile_real(values, 0.50_dp)
    result%mean = sum(values) / real(n, dp)
    result%third_quartile = quantile_real(values, 0.75_dp)
    result%maximum = values(n)
    result%status = infoset_success
  end subroutine summary_ptf
end module infoset_portfolio
