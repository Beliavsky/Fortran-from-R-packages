! SPDX-License-Identifier: GPL-3.0-or-later
module frapo_portfolios
  use frapo_kinds, only : dp
  use frapo_types, only : portfolio_result, optimizer_result, frapo_ok, frapo_invalid_input
  use frapo_statistics, only : sample_covariance, covariance_to_correlation, is_symmetric
  use frapo_series, only : return_series, returns_discrete
  use frapo_risk, only : tail_dependence_coefficient, tdc_empirical
  use frapo_optimization, only : solve_simplex_qp, solve_risk_parity, solve_convex_qp
  implicit none
  private

  public :: portfolio_global_minimum_variance
  public :: portfolio_most_diversified
  public :: portfolio_minimum_tail_dependence
  public :: portfolio_equal_risk_contribution
  public :: portfolio_maximum_drawdown
  public :: portfolio_average_drawdown
  public :: portfolio_conditional_drawdown
  public :: portfolio_minimum_conditional_drawdown

  public :: pgmv, pmd, pmtd, perc, pmaxdd, pavedd, pcdar, pmincdar

  interface pgmv
    module procedure portfolio_global_minimum_variance
  end interface
  interface pmd
    module procedure portfolio_most_diversified
  end interface
  interface pmtd
    module procedure portfolio_minimum_tail_dependence
  end interface
  interface perc
    module procedure portfolio_equal_risk_contribution
  end interface
  interface pmaxdd
    module procedure portfolio_maximum_drawdown
  end interface
  interface pavedd
    module procedure portfolio_average_drawdown
  end interface
  interface pcdar
    module procedure portfolio_conditional_drawdown
  end interface
  interface pmincdar
    module procedure portfolio_minimum_conditional_drawdown
  end interface

contains

  function portfolio_global_minimum_variance(returns, percentage) result(result)
    real(dp), intent(in) :: returns(:, :)
    logical, intent(in), optional :: percentage
    type(portfolio_result) :: result
    real(dp), allocatable :: covariance(:, :), weights(:)
    logical :: pct
    integer :: status = frapo_ok, iterations = 0

    pct = .true.
    if (present(percentage)) pct = percentage
    call sample_covariance(returns, covariance, status)
    if (status == frapo_ok) call solve_simplex_qp(covariance, weights, status, iterations)
    result%status = status
    result%iterations = iterations
    result%portfolio_type = 'Global Minimum Variance'
    if (status == frapo_ok) then
      result%weights = weights
      result%objective = dot_product(weights, matmul(covariance, weights))
      if (pct) result%weights = 100.0_dp * result%weights
    else
      allocate(result%weights(size(returns, 2)), source=0.0_dp)
    end if
  end function portfolio_global_minimum_variance

  function portfolio_most_diversified(returns, percentage) result(result)
    real(dp), intent(in) :: returns(:, :)
    logical, intent(in), optional :: percentage
    type(portfolio_result) :: result
    real(dp), allocatable :: covariance(:, :), correlation(:, :), x(:), weights(:), sd(:)
    logical :: pct
    integer :: status, iterations, i, n

    pct = .true.
    if (present(percentage)) pct = percentage
    n = size(returns, 2)
    call sample_covariance(returns, covariance, status)
    if (status == frapo_ok) call covariance_to_correlation(covariance, correlation, status)
    if (status == frapo_ok) call solve_simplex_qp(correlation, x, status, iterations)
    result%status = status
    result%iterations = iterations
    result%portfolio_type = 'Most Diversified'
    allocate(result%weights(n))
    if (status == frapo_ok) then
      allocate(sd(n), weights(n))
      do i = 1, n
        sd(i) = sqrt(covariance(i, i))
      end do
      weights = x / sd
      weights = weights / sum(weights)
      result%weights = weights
      result%objective = dot_product(x, matmul(correlation, x))
      if (pct) result%weights = 100.0_dp * result%weights
    else
      result%weights = 0.0_dp
    end if
  end function portfolio_most_diversified

  function portfolio_minimum_tail_dependence(returns, method, k, percentage) result(result)
    real(dp), intent(in) :: returns(:, :)
    integer, intent(in), optional :: method, k
    logical, intent(in), optional :: percentage
    type(portfolio_result) :: result
    real(dp), allocatable :: covariance(:, :), dependence(:, :), x(:), weights(:), sd(:)
    logical :: pct
    integer :: status, iterations, i, n, meth

    pct = .true.
    if (present(percentage)) pct = percentage
    meth = tdc_empirical
    if (present(method)) meth = method
    n = size(returns, 2)
    dependence = tail_dependence_coefficient(returns, method=meth, k=k, status=status)
    if (status == frapo_ok) call solve_simplex_qp(2.0_dp * dependence, x, status, iterations)
    if (status == frapo_ok) call sample_covariance(returns, covariance, status)
    allocate(result%weights(n))
    result%status = status
    result%iterations = iterations
    result%portfolio_type = 'Minimum Tail Dependent'
    if (status == frapo_ok) then
      allocate(sd(n), weights(n))
      do i = 1, n
        sd(i) = sqrt(covariance(i, i))
      end do
      weights = x / sd
      weights = weights / sum(weights)
      result%weights = weights
      result%objective = dot_product(x, matmul(dependence, x))
      if (pct) result%weights = 100.0_dp * result%weights
    else
      result%weights = 0.0_dp
    end if
  end function portfolio_minimum_tail_dependence

  function portfolio_equal_risk_contribution(covariance, initial, percentage) result(result)
    real(dp), intent(in) :: covariance(:, :)
    real(dp), intent(in), optional :: initial(:)
    logical, intent(in), optional :: percentage
    type(portfolio_result) :: result
    real(dp), allocatable :: weights(:)
    logical :: pct
    integer :: status = frapo_ok, iterations = 0

    pct = .true.
    if (present(percentage)) pct = percentage
    result%portfolio_type = 'Equal Risk Contribution'
    if (.not. is_symmetric(covariance)) then
      result%status = frapo_invalid_input
      allocate(result%weights(size(covariance, 1)), source=0.0_dp)
      return
    end if
    call solve_risk_parity(covariance, weights, initial=initial, status=status, iterations=iterations)
    result%weights = weights
    result%status = status
    result%iterations = iterations
    result%objective = dot_product(weights, matmul(covariance, weights))
    if (pct) result%weights = 100.0_dp * result%weights
  end function portfolio_equal_risk_contribution

  function portfolio_maximum_drawdown(price_data, max_drawdown, soft_budget) result(result)
    real(dp), intent(in) :: price_data(:, :)
    real(dp), intent(in), optional :: max_drawdown
    logical, intent(in), optional :: soft_budget
    type(portfolio_result) :: result
    real(dp), allocatable :: rc(:, :), p(:, :), q(:), aeq(:, :), beq(:), g(:, :), h(:)
    type(optimizer_result) :: opt
    real(dp) :: bound
    logical :: soft
    integer :: n, jn, nv, me, mi, row, i, t

    bound = 0.1_dp
    if (present(max_drawdown)) bound = max_drawdown
    soft = .false.
    if (present(soft_budget)) soft = soft_budget
    if (bound <= 0.0_dp .or. bound >= 1.0_dp) then
      result%status = frapo_invalid_input
      return
    end if
    rc = return_series(price_data, method=returns_discrete, percentage=.false., compound=.true.)
    jn = size(rc, 1); n = size(rc, 2); nv = n + jn
    me = merge(1, 2, soft)
    mi = n + merge(1, 0, soft) + 2 * jn + jn - 1
    allocate(p(nv, nv), q(nv), aeq(me, nv), beq(me), g(mi, nv), h(mi))
    p = 0.0_dp; q = 0.0_dp; aeq = 0.0_dp; beq = 0.0_dp; g = 0.0_dp; h = 0.0_dp
    do i = 1, nv
      p(i, i) = 1.0e-12_dp
    end do
    q(1:n) = -rc(jn, :)
    row = 0
    if (.not. soft) then
      row = row + 1; aeq(row, 1:n) = 1.0_dp; beq(row) = 1.0_dp
    end if
    row = row + 1; aeq(row, n + 1) = 1.0_dp
    row = 0
    do i = 1, n
      row = row + 1; g(row, i) = -1.0_dp
    end do
    if (soft) then
      row = row + 1; g(row, 1:n) = 1.0_dp; h(row) = 1.0_dp
    end if
    do t = 1, jn
      row = row + 1; g(row, 1:n) = -rc(t, :); g(row, n + t) = 1.0_dp; h(row) = bound
    end do
    do t = 1, jn
      row = row + 1; g(row, 1:n) = rc(t, :); g(row, n + t) = -1.0_dp
    end do
    do t = 1, jn - 1
      row = row + 1; g(row, n + t) = 1.0_dp; g(row, n + t + 1) = -1.0_dp
    end do
    call solve_convex_qp(p, q, aeq, beq, g, h, opt, tolerance=1.0e-8_dp, max_iterations=250)
    call fill_drawdown_result(result, opt, rc, n, jn, 'Maximum Drawdown')
  end function portfolio_maximum_drawdown

  function portfolio_average_drawdown(price_data, average_drawdown, soft_budget) result(result)
    real(dp), intent(in) :: price_data(:, :)
    real(dp), intent(in), optional :: average_drawdown
    logical, intent(in), optional :: soft_budget
    type(portfolio_result) :: result
    real(dp), allocatable :: rc(:, :), p(:, :), q(:), aeq(:, :), beq(:), g(:, :), h(:)
    type(optimizer_result) :: opt
    real(dp) :: bound
    logical :: soft
    integer :: n, jn, nv, me, mi, row, i, t

    bound = 0.1_dp
    if (present(average_drawdown)) bound = average_drawdown
    soft = .false.; if (present(soft_budget)) soft = soft_budget
    if (bound <= 0.0_dp .or. bound >= 1.0_dp) then
      result%status = frapo_invalid_input; return
    end if
    rc = return_series(price_data, method=returns_discrete, percentage=.false., compound=.true.)
    jn = size(rc, 1); n = size(rc, 2); nv = n + 2 * jn
    me = jn + merge(0, 1, soft)
    mi = n + merge(1, 0, soft) + 1 + jn + jn - 1
    allocate(p(nv, nv), q(nv), aeq(me, nv), beq(me), g(mi, nv), h(mi))
    p = 0.0_dp; q = 0.0_dp; aeq = 0.0_dp; beq = 0.0_dp; g = 0.0_dp; h = 0.0_dp
    do i = 1, nv; p(i, i) = 1.0e-12_dp; end do
    q(1:n) = -rc(jn, :)
    row = 0
    if (.not. soft) then
      row = row + 1; aeq(row, 1:n) = 1.0_dp; beq(row) = 1.0_dp
    end if
    do t = 1, jn
      row = row + 1
      aeq(row, 1:n) = -rc(t, :); aeq(row, n + t) = 1.0_dp
      aeq(row, n + jn + t) = -1.0_dp
    end do
    row = 0
    do i = 1, n; row = row + 1; g(row, i) = -1.0_dp; end do
    if (soft) then; row = row + 1; g(row, 1:n) = 1.0_dp; h(row) = 1.0_dp; end if
    row = row + 1; g(row, n + jn + 1:n + 2 * jn) = 1.0_dp / real(jn, dp); h(row) = bound
    do t = 1, jn
      row = row + 1; g(row, 1:n) = rc(t, :); g(row, n + t) = -1.0_dp
    end do
    do t = 1, jn - 1
      row = row + 1; g(row, n + t) = 1.0_dp; g(row, n + t + 1) = -1.0_dp
    end do
    call solve_convex_qp(p, q, aeq, beq, g, h, opt, tolerance=1.0e-8_dp, max_iterations=250)
    call fill_drawdown_result(result, opt, rc, n, jn, 'Average Drawdown')
    if (allocated(result%drawdowns)) result%risk_value = sum(result%drawdowns) / real(jn, dp)
  end function portfolio_average_drawdown

  function portfolio_conditional_drawdown(price_data, alpha, bound, soft_budget) result(result)
    real(dp), intent(in) :: price_data(:, :)
    real(dp), intent(in), optional :: alpha, bound
    logical, intent(in), optional :: soft_budget
    type(portfolio_result) :: result
    real(dp) :: a, bnd
    a = 0.95_dp; if (present(alpha)) a = alpha
    bnd = 0.05_dp; if (present(bound)) bnd = bound
    result = conditional_drawdown_core(price_data, a, bnd, soft_budget, .false.)
  end function portfolio_conditional_drawdown

  function portfolio_minimum_conditional_drawdown(price_data, alpha, soft_budget) result(result)
    real(dp), intent(in) :: price_data(:, :)
    real(dp), intent(in), optional :: alpha
    logical, intent(in), optional :: soft_budget
    type(portfolio_result) :: result
    real(dp) :: a
    a = 0.95_dp; if (present(alpha)) a = alpha
    result = conditional_drawdown_core(price_data, a, 0.0_dp, soft_budget, .true.)
  end function portfolio_minimum_conditional_drawdown

  function conditional_drawdown_core(price_data, alpha, bound, soft_budget, minimize_cdar) result(result)
    real(dp), intent(in) :: price_data(:, :), alpha, bound
    logical, intent(in), optional :: soft_budget
    logical, intent(in) :: minimize_cdar
    type(portfolio_result) :: result
    real(dp), allocatable :: rc(:, :), p(:, :), q(:), aeq(:, :), beq(:), g(:, :), h(:)
    type(optimizer_result) :: opt
    logical :: soft
    integer :: n, jn, nv, me, mi, row, i, t, zidx
    real(dp) :: coeff

    soft = .false.; if (present(soft_budget)) soft = soft_budget
    if (alpha <= 0.0_dp .or. alpha >= 1.0_dp .or. (.not. minimize_cdar .and. (bound <= 0.0_dp .or. bound >= 1.0_dp))) then
      result%status = frapo_invalid_input; return
    end if
    rc = return_series(price_data, method=returns_discrete, percentage=.false., compound=.true.)
    jn = size(rc, 1); n = size(rc, 2); nv = n + 2 * jn + 1; zidx = nv
    me = merge(1, 2, soft)
    mi = n + merge(1, 0, soft) + jn + jn + merge(0, 1, minimize_cdar) + jn + jn - 1
    allocate(p(nv, nv), q(nv), aeq(me, nv), beq(me), g(mi, nv), h(mi))
    p = 0.0_dp; q = 0.0_dp; aeq = 0.0_dp; beq = 0.0_dp; g = 0.0_dp; h = 0.0_dp
    do i = 1, nv; p(i, i) = 1.0e-12_dp; end do
    coeff = 1.0_dp / (real(jn, dp) * (1.0_dp - alpha))
    if (minimize_cdar) then
      q(n + jn + 1:n + 2 * jn) = coeff; q(zidx) = 1.0_dp
    else
      q(1:n) = -rc(jn, :)
    end if
    row = 0
    if (.not. soft) then; row = row + 1; aeq(row, 1:n) = 1.0_dp; beq(row) = 1.0_dp; end if
    row = row + 1; aeq(row, n + 1) = 1.0_dp
    row = 0
    do i = 1, n; row = row + 1; g(row, i) = -1.0_dp; end do
    if (soft) then; row = row + 1; g(row, 1:n) = 1.0_dp; h(row) = 1.0_dp; end if
    do t = 1, jn
      row = row + 1
      g(row, 1:n) = -rc(t, :); g(row, n + t) = 1.0_dp
      g(row, n + jn + t) = -1.0_dp; g(row, zidx) = -1.0_dp
    end do
    do t = 1, jn
      row = row + 1; g(row, n + jn + t) = -1.0_dp
    end do
    if (.not. minimize_cdar) then
      row = row + 1; g(row, n + jn + 1:n + 2 * jn) = coeff; g(row, zidx) = 1.0_dp; h(row) = bound
    end if
    do t = 1, jn
      row = row + 1; g(row, 1:n) = rc(t, :); g(row, n + t) = -1.0_dp
    end do
    do t = 1, jn - 1
      row = row + 1; g(row, n + t) = 1.0_dp; g(row, n + t + 1) = -1.0_dp
    end do
    call solve_convex_qp(p, q, aeq, beq, g, h, opt, tolerance=1.0e-8_dp, max_iterations=300)
    call fill_drawdown_result(result, opt, rc, n, jn, &
      merge('Minimum Conditional Drawdown', 'Conditional Drawdown        ', minimize_cdar))
    if (opt%status == frapo_ok) then
      result%threshold = empirical_quantile(result%drawdowns, alpha)
      result%risk_value = cdar_from_drawdowns(result%drawdowns, result%threshold)
      if (minimize_cdar) result%objective = result%risk_value
    end if
  end function conditional_drawdown_core

  subroutine fill_drawdown_result(result, opt, rc, n, jn, label)
    type(portfolio_result), intent(out) :: result
    type(optimizer_result), intent(in) :: opt
    real(dp), intent(in) :: rc(:, :)
    integer, intent(in) :: n, jn
    character(len=*), intent(in) :: label
    real(dp), allocatable :: equity(:), high_water(:)
    integer :: t

    result%status = opt%status
    result%iterations = opt%iterations
    result%portfolio_type = label
    allocate(result%weights(n), result%drawdowns(jn))
    if (opt%status /= frapo_ok) then
      result%weights = 0.0_dp; result%drawdowns = 0.0_dp; return
    end if
    result%weights = opt%x(1:n)
    equity = matmul(rc, result%weights)
    allocate(high_water(jn))
    high_water(1) = max(0.0_dp, equity(1))
    do t = 2, jn
      high_water(t) = max(high_water(t - 1), equity(t))
    end do
    result%drawdowns = high_water - equity
    result%terminal_return = equity(jn)
    result%objective = result%terminal_return
    result%risk_value = maxval(result%drawdowns)
  end subroutine fill_drawdown_result


  pure real(dp) function empirical_quantile(values, probability) result(value)
    real(dp), intent(in) :: values(:), probability
    real(dp) :: sorted(size(values)), temp, position, fraction
    integer :: i, j, lower_index, n

    sorted = values
    n = size(sorted)
    do i = 2, n
      temp = sorted(i)
      j = i - 1
      do while (j >= 1)
        if (sorted(j) <= temp) exit
        sorted(j + 1) = sorted(j)
        j = j - 1
      end do
      sorted(j + 1) = temp
    end do
    if (n == 1) then
      value = sorted(1)
      return
    end if
    position = 1.0_dp + real(n - 1, dp) * min(max(probability, 0.0_dp), 1.0_dp)
    lower_index = min(n - 1, max(1, int(floor(position))))
    fraction = position - real(lower_index, dp)
    value = (1.0_dp - fraction) * sorted(lower_index) + fraction * sorted(lower_index + 1)
  end function empirical_quantile

  pure real(dp) function cdar_from_drawdowns(drawdowns, threshold) result(value)
    real(dp), intent(in) :: drawdowns(:), threshold
    integer :: i, count
    value = 0.0_dp; count = 0
    do i = 1, size(drawdowns)
      if (drawdowns(i) >= threshold) then
        value = value + drawdowns(i); count = count + 1
      end if
    end do
    if (count > 0) value = value / real(count, dp)
  end function cdar_from_drawdowns
end module frapo_portfolios
