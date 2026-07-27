! SPDX-License-Identifier: MIT
! Copyright (c) 2020 RTL Authors
module rtl_portfolio
  use rtl_kinds, only: dp
  use rtl_types, only: frontier_result, lp_result, refinery_result
  use rtl_stats, only: mean_value, sample_sd, kendall_correlation, nearest_psd, seed_random
  implicit none
  private

  public :: efficient_frontier, efficient_frontier_statistics
  public :: simplex_maximize, refinery_lp
  public :: efficientFrontier, refineryLP

  interface efficientFrontier
    module procedure efficient_frontier
  end interface efficientFrontier

  interface refineryLP
    module procedure refinery_lp
  end interface refineryLP

contains

  function efficient_frontier(nsims, prices, expected_returns, seed) result(output)
    integer, intent(in) :: nsims
    real(dp), intent(in) :: prices(:, :)
    real(dp), intent(in), optional :: expected_returns(:)
    integer, intent(in), optional :: seed
    type(frontier_result) :: output
    real(dp), allocatable :: changes(:, :), means(:), sds(:), correlation(:, :), covariance(:, :)
    integer :: n, p, i
    logical :: adjusted

    n = size(prices, 1)
    p = size(prices, 2)
    if (n < 3 .or. p < 1) then
      output%status%ok = .false.
      output%status%message = "efficient_frontier requires at least three rows and one asset"
      return
    end if
    allocate(changes(n - 1, p), means(p), sds(p), correlation(p, p), covariance(p, p))
    changes = prices(2:n, :) - prices(1:n - 1, :)
    do i = 1, p
      means(i) = mean_value(changes(:, i))
      sds(i) = sample_sd(changes(:, i))
    end do
    if (present(expected_returns)) then
      if (size(expected_returns) /= p) then
        output%status%ok = .false.
        output%status%message = "expected_returns length does not match number of assets"
        return
      end if
      means = expected_returns
    end if
    call kendall_correlation(changes, correlation)
    do i = 1, p
      covariance(i, :) = sds(i) * sds * correlation(i, :)
    end do
    call nearest_psd(covariance, adjusted)
    output = efficient_frontier_statistics(nsims, means, covariance, seed)
    if (adjusted .and. output%status%ok) output%status%message = "covariance adjusted to positive semidefinite"
  end function efficient_frontier

  function efficient_frontier_statistics(nsims, expected_returns, covariance, seed) result(output)
    integer, intent(in) :: nsims
    real(dp), intent(in) :: expected_returns(:), covariance(:, :)
    integer, intent(in), optional :: seed
    type(frontier_result) :: output
    real(dp), allocatable :: weights(:)
    real(dp) :: total_weight, variance
    integer :: p, rows, i, j

    p = size(expected_returns)
    if (nsims < 0 .or. p < 1 .or. size(covariance, 1) /= p .or. size(covariance, 2) /= p) then
      output%status%ok = .false.
      output%status%message = "invalid efficient-frontier dimensions"
      return
    end if
    rows = p + nsims
    allocate(output%weights(rows, p), output%expected_return(rows), output%risk(rows), output%sharpe(rows))
    allocate(weights(p))
    output%weights = 0.0_dp
    do i = 1, p
      output%weights(i, i) = 1.0_dp
    end do
    if (present(seed)) call seed_random(seed)
    do i = p + 1, rows
      call random_number(weights)
      total_weight = sum(weights)
      if (total_weight <= 0.0_dp) then
        weights = 1.0_dp / real(p, dp)
      else
        weights = weights / total_weight
      end if
      output%weights(i, :) = weights
    end do
    do i = 1, rows
      output%expected_return(i) = dot_product(output%weights(i, :), expected_returns)
      variance = dot_product(output%weights(i, :), matmul(covariance, output%weights(i, :)))
      output%risk(i) = sqrt(max(0.0_dp, variance))
      if (output%risk(i) > 0.0_dp) then
        output%sharpe(i) = output%expected_return(i) / output%risk(i)
      else
        output%sharpe(i) = 0.0_dp
      end if
    end do
    output%minimum_risk_index = minloc(output%risk, dim=1)
    output%maximum_sharpe_index = maxloc(output%sharpe, dim=1)
    do i = 1, rows
      do j = 1, p
        if (abs(output%weights(i, j)) < 5.0e-16_dp) output%weights(i, j) = 0.0_dp
      end do
    end do
  end function efficient_frontier_statistics

  function simplex_maximize(objective, constraint_matrix, rhs, max_iterations) result(output)
    real(dp), intent(in) :: objective(:), constraint_matrix(:, :), rhs(:)
    integer, intent(in), optional :: max_iterations
    type(lp_result) :: output
    real(dp), allocatable :: tableau(:, :), ratios(:), pivot_row_values(:)
    real(dp) :: pivot_value, best_value
    integer :: m, n, columns, i, j, entering, leaving, iteration, limit, basic_column
    logical :: has_entering

    m = size(rhs)
    n = size(objective)
    if (m < 1 .or. n < 1 .or. size(constraint_matrix, 1) /= m .or. &
        size(constraint_matrix, 2) /= n .or. any(rhs < 0.0_dp)) then
      output%status%ok = .false.
      output%status%message = "simplex requires Ax <= b with b >= 0 and x >= 0"
      return
    end if
    limit = 10000
    if (present(max_iterations)) limit = max_iterations
    columns = n + m + 1
    allocate(tableau(m + 1, columns), ratios(m), pivot_row_values(columns))
    tableau = 0.0_dp
    tableau(1:m, 1:n) = constraint_matrix
    do i = 1, m
      tableau(i, n + i) = 1.0_dp
      tableau(i, columns) = rhs(i)
    end do
    tableau(m + 1, 1:n) = -objective

    do iteration = 1, limit
      entering = 0
      best_value = -1.0e-12_dp
      do j = 1, columns - 1
        if (tableau(m + 1, j) < best_value) then
          best_value = tableau(m + 1, j)
          entering = j
        end if
      end do
      has_entering = entering > 0
      if (.not. has_entering) exit
      ratios = huge(1.0_dp)
      do i = 1, m
        if (tableau(i, entering) > 1.0e-14_dp) ratios(i) = tableau(i, columns) / tableau(i, entering)
      end do
      leaving = minloc(ratios, dim=1)
      if (ratios(leaving) >= huge(1.0_dp) / 2.0_dp) then
        output%status%ok = .false.
        output%status%message = "linear program is unbounded"
        return
      end if
      pivot_value = tableau(leaving, entering)
      tableau(leaving, :) = tableau(leaving, :) / pivot_value
      pivot_row_values = tableau(leaving, :)
      do i = 1, m + 1
        if (i == leaving) cycle
        tableau(i, :) = tableau(i, :) - tableau(i, entering) * pivot_row_values
      end do
    end do
    if (iteration > limit) then
      output%status%ok = .false.
      output%status%message = "simplex iteration limit reached"
      return
    end if
    allocate(output%solution(n))
    output%solution = 0.0_dp
    do j = 1, n
      basic_column = 0
      do i = 1, m
        if (abs(tableau(i, j) - 1.0_dp) <= 1.0e-10_dp) then
          if (basic_column == 0) then
            basic_column = i
          else
            basic_column = -1
          end if
        else if (abs(tableau(i, j)) > 1.0e-10_dp) then
          basic_column = -1
        end if
      end do
      if (basic_column > 0) output%solution(j) = tableau(basic_column, columns)
    end do
    output%objective = tableau(m + 1, columns)
    output%iterations = iteration - 1
  end function simplex_maximize

  function refinery_lp(crude_prices, processing_fees, product_prices, yields, maximum_products) result(output)
    real(dp), intent(in) :: crude_prices(:), processing_fees(:), product_prices(:)
    real(dp), intent(in) :: yields(:, :), maximum_products(:)
    type(refinery_result) :: output
    type(lp_result) :: lp
    integer :: n_crudes

    n_crudes = size(crude_prices)
    if (size(processing_fees) /= n_crudes .or. size(yields, 2) /= n_crudes .or. &
        size(yields, 1) /= size(product_prices) .or. size(maximum_products) /= size(product_prices)) then
      output%status%ok = .false.
      output%status%message = "refinery input dimensions do not match"
      return
    end if
    allocate(output%margin(n_crudes))
    output%margin = matmul(transpose(yields), product_prices) - crude_prices - processing_fees
    lp = simplex_maximize(output%margin, yields, maximum_products)
    if (.not. lp%status%ok) then
      output%status = lp%status
      return
    end if
    output%profit = lp%objective
    output%slate = lp%solution
  end function refinery_lp

end module rtl_portfolio
