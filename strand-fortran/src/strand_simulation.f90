! SPDX-License-Identifier: GPL-3.0-only
! Upstream authors: Jeff Enos, David Kane, and strand contributors.
! Numerical translation of strand 0.2.3.
module strand_simulation
  use strand_kinds, only : dp
  use strand_optimizer, only : optimize_portfolio
  use strand_types, only : optimizer_config, simulation_config, optimization_result, &
    factor_constraint, category_constraint, day_result, simulation_result
  implicit none
  private
  public :: simulate_day, simulate_portfolio, allocate_market_and_transfer_orders

contains

  function simulate_day(opt_config, sim_config, start_price, end_price, volume_shares, &
    dividend, distribution, adjustment_ratio, expected_dollar_volume, investable, alpha, &
    int_shares, ext_shares, financing_days, factor_constraints, category_constraints, &
    delisting, delisting_return) result(day)
    type(optimizer_config), intent(in) :: opt_config
    type(simulation_config), intent(in) :: sim_config
    real(dp), intent(in) :: start_price(:), end_price(:), volume_shares(:)
    real(dp), intent(in) :: dividend(:), distribution(:), adjustment_ratio(:)
    real(dp), intent(in) :: expected_dollar_volume(:), alpha(:, :)
    logical, intent(in) :: investable(:)
    integer, intent(inout) :: int_shares(:, :), ext_shares(:, :)
    integer, intent(in), optional :: financing_days
    type(factor_constraint), intent(in), optional :: factor_constraints(:)
    type(category_constraint), intent(in), optional :: category_constraints(:)
    logical, intent(in), optional :: delisting(:)
    real(dp), intent(in), optional :: delisting_return(:)
    type(day_result) :: day
    type(optimization_result) :: opt
    integer, allocatable :: start_total(:, :), opt_start(:, :), orders(:, :)
    integer, allocatable :: market_order(:, :), transfer_order(:, :)
    integer, allocatable :: market_fill(:, :), transfer_fill(:, :)
    integer, allocatable :: joint_order(:), buy_total(:), sell_total(:), max_fill(:)
    logical, allocatable :: is_delisting(:)
    real(dp), allocatable :: fill_rate(:), start_nmv(:, :), max_trim(:)
    real(dp), allocatable :: delist_ret(:)
    real(dp) :: trim_gmv, pos_nmv
    integer :: n, s, i, j, days

    n = size(start_price)
    s = size(alpha, 2)
    if (.not. valid_day_dimensions(n, s, end_price, volume_shares, dividend, distribution, &
        adjustment_ratio, expected_dollar_volume, investable, alpha, int_shares, ext_shares)) then
      day%message = 'simulate_day: dimension mismatch'
      return
    end if
    if (any(start_price <= 0.0_dp) .or. any(end_price <= 0.0_dp) .or. &
        any(adjustment_ratio <= 0.0_dp) .or. sim_config%fill_rate_pct_volume <= 0.0_dp) then
      day%message = 'simulate_day: invalid price, ratio, or fill configuration'
      return
    end if

    do j = 1, s
      do i = 1, n
        int_shares(i, j) = nint(real(int_shares(i, j), dp) / adjustment_ratio(i))
        ext_shares(i, j) = nint(real(ext_shares(i, j), dp) / adjustment_ratio(i))
      end do
    end do
    allocate(start_total(n, s), opt_start(n, s), orders(n, s), is_delisting(n), delist_ret(n))
    start_total = int_shares + ext_shares
    opt_start = start_total
    is_delisting = .false.
    delist_ret = 0.0_dp
    if (present(delisting)) then
      if (size(delisting) /= n) then
        day%message = 'simulate_day: delisting dimension mismatch'
        return
      end if
      is_delisting = delisting
      do j = 1, s
        where (is_delisting) opt_start(:, j) = 0
      end do
    end if
    if (present(delisting_return)) then
      if (size(delisting_return) /= n) then
        day%message = 'simulate_day: delisting return dimension mismatch'
        return
      end if
      delist_ret = delisting_return
      if (any(delist_ret < -1.0_dp)) then
        day%message = 'simulate_day: delisting return below -1'
        return
      end if
    end if

    if (present(factor_constraints) .and. present(category_constraints)) then
      opt = optimize_portfolio(opt_config, start_price, expected_dollar_volume, investable, alpha, &
        opt_start, factor_constraints, category_constraints)
    else if (present(factor_constraints)) then
      opt = optimize_portfolio(opt_config, start_price, expected_dollar_volume, investable, alpha, &
        opt_start, factor_constraints=factor_constraints)
    else if (present(category_constraints)) then
      opt = optimize_portfolio(opt_config, start_price, expected_dollar_volume, investable, alpha, &
        opt_start, category_constraints=category_constraints)
    else
      opt = optimize_portfolio(opt_config, start_price, expected_dollar_volume, investable, alpha, opt_start)
    end if
    if (.not. opt%success) then
      day%message = 'simulate_day optimization failed: '//trim(opt%message)
      return
    end if
    orders = opt%order_shares

    if (sim_config%force_trim_factor >= 1.0_dp) then
      allocate(max_trim(n))
      do j = 1, s
        do i = 1, n
          if (is_delisting(i)) cycle
          pos_nmv = real(start_total(i, j), dp) * start_price(i)
          if (start_total(i, j) > 0) then
            max_trim(i) = opt%max_pos_long(i, j) * sim_config%force_trim_factor
          else
            max_trim(i) = opt%max_pos_short(i, j) * sim_config%force_trim_factor
          end if
          if (abs(pos_nmv) > abs(max_trim(i))) then
            trim_gmv = min(abs(pos_nmv) - abs(max_trim(i)), opt%max_order_gmv(i, j))
            if (trim_gmv > 0.0_dp) then
              if (abs(orders(i, j)) < floor(trim_gmv / start_price(i))) then
                orders(i, j) = -sign(1, start_total(i, j)) * int(floor(trim_gmv / start_price(i)))
              end if
            end if
          end if
        end do
      end do
    end if

    if (sim_config%force_exit_non_investable) then
      do j = 1, s
        do i = 1, n
          if (investable(i) .or. is_delisting(i) .or. start_total(i, j) == 0) cycle
          trim_gmv = min(abs(real(start_total(i, j), dp) * start_price(i)), opt%max_order_gmv(i, j))
          orders(i, j) = -sign(1, start_total(i, j)) * int(floor(trim_gmv / start_price(i)))
        end do
      end do
    end if

    do j = 1, s
      do i = 1, n
        if (is_delisting(i)) orders(i, j) = -start_total(i, j)
      end do
    end do

    allocate(joint_order(n), buy_total(n), sell_total(n), max_fill(n), fill_rate(n))
    allocate(market_order(n, s), transfer_order(n, s), market_fill(n, s), transfer_fill(n, s))
    joint_order = sum(orders, dim=2)
    buy_total = 0
    sell_total = 0
    do j = 1, s
      where (orders(:, j) > 0) buy_total = buy_total + orders(:, j)
      where (orders(:, j) < 0) sell_total = sell_total + orders(:, j)
    end do
    max_fill = nint(volume_shares * sim_config%fill_rate_pct_volume / 100.0_dp)
    fill_rate = 1.0_dp
    do i = 1, n
      if (joint_order(i) /= 0 .and. abs(joint_order(i)) > max_fill(i) .and. .not. is_delisting(i)) then
        fill_rate(i) = real(max_fill(i), dp) / real(abs(joint_order(i)), dp)
      end if
    end do
    call allocate_market_and_transfer_orders(orders, joint_order, buy_total, sell_total, &
      market_order, transfer_order)
    do j = 1, s
      do i = 1, n
        market_fill(i, j) = nint(real(market_order(i, j), dp) * fill_rate(i))
        transfer_fill(i, j) = transfer_order(i, j)
      end do
    end do

    allocate(day%start_shares(n, s), day%order_shares(n, s), day%market_order_shares(n, s), &
      day%transfer_order_shares(n, s), day%market_fill_shares(n, s), &
      day%transfer_fill_shares(n, s), day%end_shares(n, s), day%start_nmv(n, s), &
      day%end_nmv(n, s), day%position_pnl(n, s), day%trade_costs(n, s), &
      day%financing_costs(n, s), day%gross_pnl(n, s), day%net_pnl(n, s), day%fill_rate(n))
    day%start_shares = start_total
    day%order_shares = orders
    day%market_order_shares = market_order
    day%transfer_order_shares = transfer_order
    day%market_fill_shares = market_fill
    day%transfer_fill_shares = transfer_fill
    day%fill_rate = fill_rate
    days = 1
    if (present(financing_days)) days = max(1, financing_days)
    allocate(start_nmv(n, s))
    do j = 1, s
      int_shares(:, j) = int_shares(:, j) + transfer_fill(:, j)
      ext_shares(:, j) = ext_shares(:, j) + market_fill(:, j)
      day%end_shares(:, j) = int_shares(:, j) + ext_shares(:, j)
      day%start_nmv(:, j) = real(start_total(:, j), dp) * start_price
      day%end_nmv(:, j) = real(day%end_shares(:, j), dp) * end_price
      day%position_pnl(:, j) = real(start_total(:, j), dp) * &
        (end_price + dividend + distribution - start_price)
      day%trade_costs(:, j) = abs(real(market_fill(:, j), dp) * end_price) * &
        sim_config%transaction_cost_pct / 100.0_dp
      day%financing_costs(:, j) = abs(real(ext_shares(:, j) - market_fill(:, j), dp) * start_price) * &
        sim_config%financing_cost_pct / 100.0_dp / 360.0_dp * real(days, dp)
      day%gross_pnl(:, j) = day%position_pnl(:, j)
      day%net_pnl(:, j) = day%gross_pnl(:, j) - day%trade_costs(:, j) - day%financing_costs(:, j)
      do i = 1, n
        if (is_delisting(i)) then
          day%position_pnl(i, j) = 0.0_dp
          day%trade_costs(i, j) = 0.0_dp
          day%financing_costs(i, j) = 0.0_dp
          day%gross_pnl(i, j) = day%start_nmv(i, j) * delist_ret(i)
          day%net_pnl(i, j) = day%gross_pnl(i, j)
          int_shares(i, j) = 0
          ext_shares(i, j) = 0
          day%end_shares(i, j) = 0
          day%end_nmv(i, j) = 0.0_dp
        end if
      end do
    end do
    day%success = .true.
    day%message = 'ok'
  end function simulate_day

  function simulate_portfolio(opt_config, sim_config, start_price, end_price, volume_shares, &
    dividend, distribution, adjustment_ratio, expected_dollar_volume, investable, alpha, &
    initial_shares, financing_days, factor_constraints, category_constraints, delisting, &
    delisting_return) result(result)
    type(optimizer_config), intent(in) :: opt_config
    type(simulation_config), intent(in) :: sim_config
    real(dp), intent(in) :: start_price(:, :), end_price(:, :), volume_shares(:, :)
    real(dp), intent(in) :: dividend(:, :), distribution(:, :), adjustment_ratio(:, :)
    real(dp), intent(in) :: expected_dollar_volume(:, :), alpha(:, :, :)
    logical, intent(in) :: investable(:, :)
    integer, intent(in) :: initial_shares(:, :)
    integer, intent(in), optional :: financing_days(:)
    type(factor_constraint), intent(in), optional :: factor_constraints(:)
    type(category_constraint), intent(in), optional :: category_constraints(:)
    logical, intent(in), optional :: delisting(:, :)
    real(dp), intent(in), optional :: delisting_return(:, :)
    type(simulation_result) :: result
    integer, allocatable :: int_shares(:, :), ext_shares(:, :)
    logical, allocatable :: day_delisting(:)
    real(dp), allocatable :: day_delisting_return(:)
    integer :: t, n, s, d, fd, j

    t = size(start_price, 1)
    n = size(start_price, 2)
    s = size(initial_shares, 2)
    if (size(initial_shares, 1) /= n .or. size(alpha, 1) /= t .or. size(alpha, 2) /= n .or. &
        size(alpha, 3) /= s .or. any(shape(end_price) /= [t, n]) .or. &
        any(shape(volume_shares) /= [t, n]) .or. any(shape(dividend) /= [t, n]) .or. &
        any(shape(distribution) /= [t, n]) .or. any(shape(adjustment_ratio) /= [t, n]) .or. &
        any(shape(expected_dollar_volume) /= [t, n]) .or. any(shape(investable) /= [t, n])) then
      result%message = 'simulate_portfolio: dimension mismatch'
      return
    end if
    if (present(financing_days)) then
      if (size(financing_days) /= t) then
        result%message = 'simulate_portfolio: financing_days dimension mismatch'
        return
      end if
    end if
    if (present(delisting)) then
      if (any(shape(delisting) /= [t, n])) then
        result%message = 'simulate_portfolio: delisting dimension mismatch'
        return
      end if
    end if
    if (present(delisting_return)) then
      if (any(shape(delisting_return) /= [t, n])) then
        result%message = 'simulate_portfolio: delisting_return dimension mismatch'
        return
      end if
    end if
    allocate(int_shares(n, s), ext_shares(n, s), result%day(t), result%final_shares(n, s))
    allocate(day_delisting(n), day_delisting_return(n))
    allocate(result%gross_pnl(t, s + 1), result%net_pnl(t, s + 1), result%end_gmv(t, s + 1), &
      result%end_nmv(t, s + 1), result%turnover(t, s + 1))
    int_shares = 0
    ext_shares = initial_shares
    result%gross_pnl = 0.0_dp
    result%net_pnl = 0.0_dp
    result%end_gmv = 0.0_dp
    result%end_nmv = 0.0_dp
    result%turnover = 0.0_dp
    do d = 1, t
      fd = 1
      if (present(financing_days)) fd = financing_days(d)
      day_delisting = .false.
      day_delisting_return = 0.0_dp
      if (present(delisting)) day_delisting = delisting(d, :)
      if (present(delisting_return)) day_delisting_return = delisting_return(d, :)
      if (present(factor_constraints) .and. present(category_constraints)) then
        result%day(d) = simulate_day(opt_config, sim_config, start_price(d, :), end_price(d, :), &
          volume_shares(d, :), dividend(d, :), distribution(d, :), adjustment_ratio(d, :), &
          expected_dollar_volume(d, :), investable(d, :), alpha(d, :, :), int_shares, ext_shares, &
          fd, factor_constraints, category_constraints, day_delisting, day_delisting_return)
      else if (present(factor_constraints)) then
        result%day(d) = simulate_day(opt_config, sim_config, start_price(d, :), end_price(d, :), &
          volume_shares(d, :), dividend(d, :), distribution(d, :), adjustment_ratio(d, :), &
          expected_dollar_volume(d, :), investable(d, :), alpha(d, :, :), int_shares, ext_shares, &
          fd, factor_constraints=factor_constraints, delisting=day_delisting, &
          delisting_return=day_delisting_return)
      else if (present(category_constraints)) then
        result%day(d) = simulate_day(opt_config, sim_config, start_price(d, :), end_price(d, :), &
          volume_shares(d, :), dividend(d, :), distribution(d, :), adjustment_ratio(d, :), &
          expected_dollar_volume(d, :), investable(d, :), alpha(d, :, :), int_shares, ext_shares, &
          fd, category_constraints=category_constraints, delisting=day_delisting, &
          delisting_return=day_delisting_return)
      else
        result%day(d) = simulate_day(opt_config, sim_config, start_price(d, :), end_price(d, :), &
          volume_shares(d, :), dividend(d, :), distribution(d, :), adjustment_ratio(d, :), &
          expected_dollar_volume(d, :), investable(d, :), alpha(d, :, :), int_shares, ext_shares, &
          fd, delisting=day_delisting, delisting_return=day_delisting_return)
      end if
      if (.not. result%day(d)%success) then
        result%message = 'day failed: '//trim(result%day(d)%message)
        return
      end if
      do j = 1, s
        result%gross_pnl(d, j) = sum(result%day(d)%gross_pnl(:, j))
        result%net_pnl(d, j) = sum(result%day(d)%net_pnl(:, j))
        result%end_gmv(d, j) = sum(abs(result%day(d)%end_nmv(:, j)))
        result%end_nmv(d, j) = sum(result%day(d)%end_nmv(:, j))
        result%turnover(d, j) = sum(abs(real(result%day(d)%market_fill_shares(:, j), dp) * end_price(d, :)))
      end do
      result%gross_pnl(d, s + 1) = sum(result%gross_pnl(d, 1:s))
      result%net_pnl(d, s + 1) = sum(result%net_pnl(d, 1:s))
      result%end_gmv(d, s + 1) = sum(abs(sum(result%day(d)%end_nmv, dim=2)))
      result%end_nmv(d, s + 1) = sum(result%day(d)%end_nmv)
      result%turnover(d, s + 1) = sum(abs(real(sum(result%day(d)%market_fill_shares, dim=2), dp) * end_price(d, :)))
    end do
    result%final_shares = int_shares + ext_shares
    result%success = .true.
    result%message = 'ok'
  end function simulate_portfolio

  subroutine allocate_market_and_transfer_orders(orders, joint_order, buy_total, sell_total, &
    market_order, transfer_order)
    integer, intent(in) :: orders(:, :), joint_order(:), buy_total(:), sell_total(:)
    integer, intent(out) :: market_order(:, :), transfer_order(:, :)
    integer :: i, j
    real(dp) :: allocation

    market_order = 0
    do j = 1, size(orders, 2)
      do i = 1, size(orders, 1)
        if (orders(i, j) > 0 .and. joint_order(i) > 0 .and. buy_total(i) > 0) then
          allocation = real(orders(i, j), dp) / real(buy_total(i), dp) * real(joint_order(i), dp)
          market_order(i, j) = nint(allocation)
        else if (orders(i, j) < 0 .and. joint_order(i) < 0 .and. sell_total(i) < 0) then
          allocation = abs(real(orders(i, j), dp) / real(sell_total(i), dp)) * real(joint_order(i), dp)
          market_order(i, j) = nint(allocation)
        end if
      end do
    end do
    transfer_order = orders - market_order
  end subroutine allocate_market_and_transfer_orders

  logical function valid_day_dimensions(n, s, end_price, volume, dividend, distribution, ratio, &
    expected_volume, investable, alpha, int_shares, ext_shares) result(ok)
    integer, intent(in) :: n, s
    real(dp), intent(in) :: end_price(:), volume(:), dividend(:), distribution(:), ratio(:), expected_volume(:)
    logical, intent(in) :: investable(:)
    real(dp), intent(in) :: alpha(:, :)
    integer, intent(in) :: int_shares(:, :), ext_shares(:, :)
    ok = size(end_price) == n .and. size(volume) == n .and. size(dividend) == n .and. &
      size(distribution) == n .and. size(ratio) == n .and. size(expected_volume) == n .and. &
      size(investable) == n .and. size(alpha, 1) == n .and. size(alpha, 2) == s .and. &
      size(int_shares, 1) == n .and. size(int_shares, 2) == s .and. &
      size(ext_shares, 1) == n .and. size(ext_shares, 2) == s
  end function valid_day_dimensions

end module strand_simulation
