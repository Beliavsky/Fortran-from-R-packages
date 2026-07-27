! SPDX-License-Identifier: GPL-3.0-only
! Upstream authors: Jeff Enos, David Kane, and strand contributors.
! Numerical translation of strand 0.2.3.
module strand_optimizer
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use strand_kinds, only : dp
  use strand_simplex, only : solve_bounded_lp
  use strand_types, only : optimizer_config, optimization_result, factor_constraint, &
    category_constraint, lp_result
  implicit none
  private
  public :: optimize_portfolio

contains

  function optimize_portfolio(config, price, expected_dollar_volume, investable, alpha, &
    start_shares, factor_constraints, category_constraints) result(result)
    type(optimizer_config), intent(in) :: config
    real(dp), intent(in) :: price(:), expected_dollar_volume(:), alpha(:, :)
    logical, intent(in) :: investable(:)
    integer, intent(in) :: start_shares(:, :)
    type(factor_constraint), intent(in), optional :: factor_constraints(:)
    type(category_constraint), intent(in), optional :: category_constraints(:)
    type(optimization_result) :: result
    real(dp), allocatable :: objective(:), lower(:), upper(:), amat(:, :), rhs(:), rhs_work(:)
    real(dp), allocatable :: current_nmv(:, :), target_long(:), target_short(:)
    logical, allocatable :: loosenable(:)
    integer, allocatable :: sense(:)
    type(lp_result) :: lp
    real(dp), allocatable :: sequence(:)
    real(dp) :: eps, current_lmv, current_smv, target_lmv, target_smv
    real(dp) :: overall_current_lmv, overall_current_smv, overall_target_lmv, overall_target_smv
    real(dp) :: effective_turnover, cumulative_loosen, coeff
    integer :: n, s, nvar, mmax, m, i, j, k, row, attempt, nf, nc, level
    integer :: block_start
    logical :: dimensions_ok

    eps = 100.0_dp * epsilon(1.0_dp)
    n = size(price)
    s = size(alpha, 2)
    dimensions_ok = n > 0 .and. size(expected_dollar_volume) == n .and. size(investable) == n
    dimensions_ok = dimensions_ok .and. size(alpha, 1) == n
    dimensions_ok = dimensions_ok .and. size(start_shares, 1) == n .and. size(start_shares, 2) == s
    dimensions_ok = dimensions_ok .and. allocated(config%strategies)
    dimensions_ok = dimensions_ok .and. size(config%strategies) == s
    if (.not. dimensions_ok) then
      result%message = 'optimize_portfolio: dimension mismatch'
      return
    end if
    if (any(price <= 0.0_dp) .or. any(expected_dollar_volume < 0.0_dp) .or. &
        any(.not. ieee_is_finite(price)) .or. any(.not. ieee_is_finite(expected_dollar_volume)) .or. &
        any(.not. ieee_is_finite(alpha))) then
      result%message = 'optimize_portfolio: invalid numeric input'
      return
    end if
    do j = 1, s
      if (config%strategies(j)%capital <= 0.0_dp .or. &
          config%strategies(j)%ideal_long_weight < 0.0_dp .or. &
          config%strategies(j)%ideal_short_weight < 0.0_dp) then
        result%message = 'optimize_portfolio: invalid strategy capital or weights'
        return
      end if
    end do

    nf = 0
    if (present(factor_constraints)) nf = size(factor_constraints)
    nc = 0
    if (present(category_constraints)) nc = size(category_constraints)
    nvar = n * s + n
    mmax = 2 * s + 2 * n + 1 + 2 * nf + 2 * n * nc
    allocate(objective(nvar), lower(nvar), upper(nvar), amat(mmax, nvar), rhs(mmax), sense(mmax), loosenable(mmax))
    allocate(current_nmv(n, s), target_long(s), target_short(s))
    allocate(result%max_pos_long(n, s), result%max_pos_short(n, s), result%max_order_gmv(n, s))
    objective = 0.0_dp
    lower = 0.0_dp
    upper = huge(1.0_dp)
    amat = 0.0_dp
    rhs = 0.0_dp
    sense = -1
    loosenable = .false.
    do j = 1, s
      current_nmv(:, j) = real(start_shares(:, j), dp) * price
      call compute_target_weights(config, j, current_nmv(:, j), target_long(j), target_short(j))
      call set_strategy_bounds(config, j, expected_dollar_volume, investable, alpha(:, j), current_nmv(:, j), &
        lower((j - 1) * n + 1:j * n), upper((j - 1) * n + 1:j * n), &
        result%max_pos_long(:, j), result%max_pos_short(:, j), result%max_order_gmv(:, j))
      objective((j - 1) * n + 1:j * n) = -alpha(:, j)
    end do
    lower(n * s + 1:nvar) = 0.0_dp
    upper(n * s + 1:nvar) = huge(1.0_dp)

    m = 0
    overall_current_lmv = 0.0_dp
    overall_current_smv = 0.0_dp
    overall_target_lmv = 0.0_dp
    overall_target_smv = 0.0_dp

    do j = 1, s
      current_lmv = sum(current_nmv(:, j), mask=current_nmv(:, j) > eps)
      current_smv = sum(current_nmv(:, j), mask=current_nmv(:, j) < -eps)
      target_lmv = config%strategies(j)%capital * target_long(j)
      target_smv = -config%strategies(j)%capital * target_short(j)
      overall_current_lmv = overall_current_lmv + current_lmv
      overall_current_smv = overall_current_smv + current_smv
      overall_target_lmv = overall_target_lmv + target_lmv
      overall_target_smv = overall_target_smv + target_smv

      m = m + 1
      block_start = (j - 1) * n
      do i = 1, n
        if (current_nmv(i, j) > eps .or. (abs(current_nmv(i, j)) <= eps .and. alpha(i, j) > 0.0_dp)) then
          amat(m, block_start + i) = 1.0_dp
        end if
      end do
      rhs(m) = target_lmv - current_lmv
      sense(m) = 0

      m = m + 1
      do i = 1, n
        if (current_nmv(i, j) < -eps .or. (abs(current_nmv(i, j)) <= eps .and. alpha(i, j) < 0.0_dp)) then
          amat(m, block_start + i) = 1.0_dp
        end if
      end do
      rhs(m) = target_smv - current_smv
      sense(m) = 0
    end do

    if (present(factor_constraints)) then
      do k = 1, size(factor_constraints)
        j = factor_constraints(k)%strategy
        if (j < 1 .or. j > s .or. .not. allocated(factor_constraints(k)%values)) then
          result%message = 'invalid factor constraint'
          return
        end if
        if (size(factor_constraints(k)%values) /= n) then
          result%message = 'factor constraint dimension mismatch'
          return
        end if
        block_start = (j - 1) * n
        coeff = dot_product(factor_constraints(k)%values, current_nmv(:, j))
        m = m + 1
        amat(m, block_start + 1:block_start + n) = factor_constraints(k)%values
        rhs(m) = factor_constraints(k)%upper_weight * config%strategies(j)%capital - coeff
        sense(m) = -1
        loosenable(m) = rhs(m) < 0.0_dp
        m = m + 1
        amat(m, block_start + 1:block_start + n) = factor_constraints(k)%values
        rhs(m) = factor_constraints(k)%lower_weight * config%strategies(j)%capital - coeff
        sense(m) = 1
        loosenable(m) = rhs(m) > 0.0_dp
      end do
    end if

    if (present(category_constraints)) then
      do k = 1, size(category_constraints)
        j = category_constraints(k)%strategy
        if (j < 1 .or. j > s .or. .not. allocated(category_constraints(k)%category)) then
          result%message = 'invalid category constraint'
          return
        end if
        if (size(category_constraints(k)%category) /= n) then
          result%message = 'category constraint dimension mismatch'
          return
        end if
        block_start = (j - 1) * n
        do i = 1, n
          level = category_constraints(k)%category(i)
          if (.not. first_level(category_constraints(k)%category, i)) cycle
          coeff = sum(current_nmv(:, j), mask=category_constraints(k)%category == level)
          m = m + 1
          do row = 1, n
            if (category_constraints(k)%category(row) == level) amat(m, block_start + row) = 1.0_dp
          end do
          rhs(m) = category_constraints(k)%upper_weight * config%strategies(j)%capital - coeff
          sense(m) = -1
          loosenable(m) = rhs(m) < 0.0_dp
          m = m + 1
          do row = 1, n
            if (category_constraints(k)%category(row) == level) amat(m, block_start + row) = 1.0_dp
          end do
          rhs(m) = category_constraints(k)%lower_weight * config%strategies(j)%capital - coeff
          sense(m) = 1
          loosenable(m) = rhs(m) > 0.0_dp
        end do
      end do
    end if

    do i = 1, n
      m = m + 1
      do j = 1, s
        amat(m, (j - 1) * n + i) = 1.0_dp
      end do
      amat(m, n * s + i) = -1.0_dp
      sense(m) = -1
      rhs(m) = 0.0_dp
      m = m + 1
      do j = 1, s
        amat(m, (j - 1) * n + i) = 1.0_dp
      end do
      amat(m, n * s + i) = 1.0_dp
      sense(m) = 1
      rhs(m) = 0.0_dp
    end do

    if (config%turnover_limit >= 0.0_dp) then
      effective_turnover = max(config%turnover_limit, &
        abs(overall_current_lmv - overall_target_lmv) + abs(overall_current_smv - overall_target_smv))
      m = m + 1
      amat(m, n * s + 1:nvar) = 1.0_dp
      rhs(m) = effective_turnover
      sense(m) = -1
    end if

    if (allocated(config%loosening_sequence)) then
      allocate(sequence(size(config%loosening_sequence)))
      sequence = config%loosening_sequence
    else
      allocate(sequence(4))
      sequence = [0.0_dp, 0.5_dp, 0.5_dp, 1.0_dp]
    end if
    allocate(rhs_work(m), result%loosened_fraction(m))
    rhs_work = rhs(1:m)
    result%loosened_fraction = 0.0_dp
    cumulative_loosen = 0.0_dp
    do attempt = 1, size(sequence)
      if (attempt > 1 .and. sequence(attempt) > 0.0_dp) then
        do row = 1, m
          if (loosenable(row)) then
            rhs_work(row) = rhs_work(row) * (1.0_dp - sequence(attempt))
            result%loosened_fraction(row) = 1.0_dp - &
              (1.0_dp - result%loosened_fraction(row)) * (1.0_dp - sequence(attempt))
          end if
        end do
        cumulative_loosen = 1.0_dp - (1.0_dp - cumulative_loosen) * (1.0_dp - sequence(attempt))
      end if
      lp = solve_bounded_lp(objective, amat(1:m, :), sense(1:m), rhs_work, lower, upper, maxiter=50000)
      if (lp%optimal) exit
    end do

    allocate(result%trade_nmv(n, s), result%order_shares(n, s), result%order_shares_joint(n), &
      result%order_nmv_joint(n), result%target_long_weight(s), result%target_short_weight(s))
    result%trade_nmv = 0.0_dp
    result%order_shares = 0
    result%order_shares_joint = 0
    result%order_nmv_joint = 0.0_dp
    result%target_long_weight = target_long
    result%target_short_weight = target_short
    result%iterations = lp%iterations
    result%message = lp%message
    if (.not. lp%optimal) return

    do j = 1, s
      result%trade_nmv(:, j) = lp%x((j - 1) * n + 1:j * n)
      do i = 1, n
        result%order_shares(i, j) = nint(result%trade_nmv(i, j) / price(i))
        result%trade_nmv(i, j) = real(result%order_shares(i, j), dp) * price(i)
      end do
    end do
    result%order_shares_joint = sum(result%order_shares, dim=2)
    result%order_nmv_joint = real(result%order_shares_joint, dp) * price
    result%objective = -lp%objective
    result%success = .true.
    result%message = 'optimal'
  end function optimize_portfolio

  subroutine compute_target_weights(config, strategy, current_nmv, target_long, target_short)
    type(optimizer_config), intent(in) :: config
    integer, intent(in) :: strategy
    real(dp), intent(in) :: current_nmv(:)
    real(dp), intent(out) :: target_long, target_short
    real(dp) :: current_long, current_short, long_change, short_change, fraction

    if (config%strategies(strategy)%has_target_weights) then
      target_long = config%strategies(strategy)%target_long_weight
      target_short = config%strategies(strategy)%target_short_weight
      return
    end if
    current_long = sum(current_nmv, mask=current_nmv > 0.0_dp) / config%strategies(strategy)%capital
    current_short = abs(sum(current_nmv, mask=current_nmv < 0.0_dp)) / config%strategies(strategy)%capital
    fraction = 1.0_dp
    if (trim(config%target_weight_policy) == 'half-way') fraction = 0.5_dp
    long_change = fraction * (config%strategies(strategy)%ideal_long_weight - current_long)
    short_change = fraction * (config%strategies(strategy)%ideal_short_weight - current_short)
    if (config%max_weight_change >= 0.0_dp) then
      long_change = sign(min(abs(long_change), config%max_weight_change * &
        config%strategies(strategy)%ideal_long_weight), long_change)
      short_change = sign(min(abs(short_change), config%max_weight_change * &
        config%strategies(strategy)%ideal_short_weight), short_change)
    end if
    target_long = current_long + long_change
    target_short = current_short + short_change
  end subroutine compute_target_weights

  subroutine set_strategy_bounds(config, strategy, volume, investable, alpha, current_nmv, &
    lower, upper, max_long, max_short, max_order)
    type(optimizer_config), intent(in) :: config
    integer, intent(in) :: strategy
    real(dp), intent(in) :: volume(:), alpha(:), current_nmv(:)
    logical, intent(in) :: investable(:)
    real(dp), intent(out) :: lower(:), upper(:), max_long(:), max_short(:), max_order(:)
    real(dp) :: ideal_lmv, ideal_smv, eps
    real(dp), allocatable :: pos_upper(:), pos_lower(:)
    integer :: i

    eps = 100.0_dp * epsilon(1.0_dp)
    ideal_lmv = config%strategies(strategy)%capital * config%strategies(strategy)%ideal_long_weight
    ideal_smv = config%strategies(strategy)%capital * config%strategies(strategy)%ideal_short_weight
    max_order = volume * config%strategies(strategy)%trading_limit_pct_adv / 100.0_dp
    max_long = min(volume * config%strategies(strategy)%position_limit_pct_adv / 100.0_dp, &
      ideal_lmv * config%strategies(strategy)%position_limit_pct_lmv / 100.0_dp)
    max_short = -min(volume * config%strategies(strategy)%position_limit_pct_adv / 100.0_dp, &
      ideal_smv * config%strategies(strategy)%position_limit_pct_smv / 100.0_dp)
    allocate(pos_upper(size(volume)), pos_lower(size(volume)))
    pos_upper = max_long
    pos_lower = max_short
    do i = 1, size(volume)
      if (.not. investable(i)) then
        pos_upper(i) = 0.0_dp
        pos_lower(i) = 0.0_dp
      end if
      if (current_nmv(i) < -eps) pos_upper(i) = 0.0_dp
      if (current_nmv(i) > eps) pos_lower(i) = 0.0_dp
      if (abs(current_nmv(i)) <= eps .and. alpha(i) <= 0.0_dp) pos_upper(i) = 0.0_dp
      if (abs(current_nmv(i)) <= eps .and. alpha(i) >= 0.0_dp) pos_lower(i) = 0.0_dp
      pos_upper(i) = min(pos_upper(i), current_nmv(i) + max_order(i))
      pos_lower(i) = max(pos_lower(i), current_nmv(i) - max_order(i))
      if (current_nmv(i) > pos_upper(i)) pos_upper(i) = current_nmv(i)
      if (current_nmv(i) < pos_lower(i)) pos_lower(i) = current_nmv(i)
      lower(i) = pos_lower(i) - current_nmv(i)
      upper(i) = pos_upper(i) - current_nmv(i)
    end do
  end subroutine set_strategy_bounds

  logical function first_level(category, position) result(first)
    integer, intent(in) :: category(:), position
    if (position <= 1) then
      first = .true.
    else
      first = .not. any(category(1:position - 1) == category(position))
    end if
  end function first_level

end module strand_optimizer
