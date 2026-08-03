! SPDX-License-Identifier: AGPL-3.0-or-later
! Derived from REN 0.1.0 computational code; see NOTICE.md.
module ren_analysis
  use ren_kinds, only : dp
  use ren_types, only : prepared_data_type, analysis_options, analysis_result, cluster_result, &
    ren_success, ren_invalid_argument, ren_dimension_error, ren_method_count, ren_method_names
  use ren_linalg, only : l1_norm, cumulative_metrics
  use ren_portfolio, only : po_cols, po_jm, po_avg, po_gross_exp, po_cov_shrink, buh_clust, &
    po_bhu, po_tzt, po_sw, po_sw_lasso
  implicit none
  private
  public :: prepare_data, perform_analysis, ren_run
contains
  function prepare_data(date, x, start_date, end_date, missing_value) result(prepared)
    integer, intent(in) :: date(:)
    real(dp), intent(in) :: x(:, :)
    integer, intent(in), optional :: start_date, end_date
    real(dp), intent(in), optional :: missing_value
    type(prepared_data_type) :: prepared
    integer, allocatable :: rows(:), columns(:)
    integer :: first_date, last_date, i, j, year0, month0, year_i, month_i, max_month
    real(dp) :: sentinel
    logical, allocatable :: keep_column(:)
    if (size(x, 1) /= size(date) .or. size(x, 2) < 1) then
      prepared%status = ren_dimension_error
      return
    end if
    first_date = minval(date)
    last_date = maxval(date)
    if (present(start_date)) first_date = start_date
    if (present(end_date)) last_date = end_date
    if (first_date > last_date) then
      prepared%status = ren_invalid_argument
      return
    end if
    rows = pack([(i, i=1,size(date))], date >= first_date .and. date <= last_date)
    if (size(rows) == 0) then
      prepared%status = ren_invalid_argument
      return
    end if
    sentinel = -99.99_dp
    if (present(missing_value)) sentinel = missing_value
    allocate(keep_column(size(x, 2)))
    keep_column = .true.
    do j = 1, size(x, 2)
      if (any(abs(x(rows, j) - sentinel) <= 10.0_dp * epsilon(1.0_dp) * max(1.0_dp, abs(sentinel)))) &
        keep_column(j) = .false.
    end do
    columns = pack([(j, j=1,size(x, 2))], keep_column)
    if (size(columns) == 0) then
      prepared%status = ren_invalid_argument
      return
    end if
    prepared%date = date(rows)
    allocate(prepared%x(size(rows), size(columns)))
    do j = 1, size(columns)
      prepared%x(:, j) = x(rows, columns(j))
    end do
    prepared%retained_columns = columns
    allocate(prepared%month(size(rows)))
    year0 = prepared%date(1) / 10000
    month0 = modulo(prepared%date(1) / 100, 100)
    do i = 1, size(rows)
      year_i = prepared%date(i) / 10000
      month_i = modulo(prepared%date(i) / 100, 100)
      prepared%month(i) = (year_i - year0) * 12 + (month_i - month0) + 1
    end do
    max_month = maxval(prepared%month)
    allocate(prepared%count(max_month))
    do i = 1, max_month
      prepared%count(i) = count(prepared%month == i)
    end do
    prepared%status = ren_success
  end function prepare_data

  subroutine perform_analysis(x, month, date, result, options)
    real(dp), intent(in) :: x(:, :)
    integer, intent(in) :: month(:), date(:)
    type(analysis_result), intent(out) :: result
    type(analysis_options), intent(in), optional :: options
    type(analysis_options) :: opt
    type(cluster_result) :: clusters
    real(dp), allocatable :: target(:, :), current(:, :), train(:, :), zero(:), method_weight(:), gross(:), &
      vw_current(:), vw_target(:)
    integer, allocatable :: train_rows(:), current_rows(:), eval_rows(:)
    integer :: n, p, max_month, eval_months, out_count, month_id, eval_index, return_index, &
      method_id, row, istat, b_assets, seed_base
    opt = analysis_options()
    if (present(options)) opt = options
    n = size(x, 1)
    p = size(x, 2)
    if (size(month) /= n .or. size(date) /= n .or. p < 1) then
      result%status = ren_dimension_error
      return
    end if
    max_month = maxval(month)
    if (minval(month) < 1 .or. max_month < 7) then
      result%status = ren_invalid_argument
      return
    end if
    eval_months = max_month - 6
    eval_rows = pack([(row, row=1,n)], month >= 7)
    out_count = size(eval_rows)
    allocate(result%method(ren_method_count))
    result%method = ren_method_names
    allocate(result%weights(eval_months, p, ren_method_count), result%turnover(eval_months, ren_method_count), &
      result%gross_returns(out_count, ren_method_count), result%cumulative_return(out_count, ren_method_count), &
      result%wealth_index(out_count, ren_method_count), result%cumulative_turnover(eval_months, ren_method_count), &
      result%turnover_mean(ren_method_count), result%sharpe_ratio(ren_method_count), &
      result%volatility(ren_method_count), result%max_drawdown(ren_method_count), &
      result%month_date(eval_months), result%return_date(out_count))
    allocate(result%vw_weights(eval_months, p), result%vw_turnover(eval_months), &
      result%vw_gross_returns(out_count), result%vw_cumulative_return(out_count), result%vw_wealth_index(out_count))
    result%weights = 0.0_dp
    result%turnover = 0.0_dp
    result%gross_returns = 1.0_dp
    allocate(target(p, ren_method_count), current(p, ren_method_count), zero(n), gross(p))
    current = 1.0_dp / real(p, dp)
    allocate(vw_current(p), vw_target(p))
    vw_current = 1.0_dp / real(p, dp)
    zero = 0.0_dp
    return_index = 0
    do month_id = 7, max_month
      eval_index = month_id - 6
      train_rows = pack([(row, row=1,n)], month < month_id .and. month >= month_id - 6)
      current_rows = pack([(row, row=1,n)], month == month_id)
      if (size(train_rows) < 3 .or. size(current_rows) == 0) then
        result%status = ren_invalid_argument
        return
      end if
      train = x(train_rows, :)
      deallocate(zero)
      allocate(zero(size(train_rows)))
      zero = 0.0_dp
      clusters = buh_clust(train, opt%variance_tolerance)
      seed_base = opt%random_seed + 10000 * month_id
      call po_cols(zero, train, method_weight, istat, opt%variance_tolerance)
      target(:, 1) = method_weight
      call po_jm(train, method_weight, istat, opt%variance_tolerance)
      target(:, 2) = method_weight
      call po_avg(zero, train, 'LASSO', method_weight, istat, opt%variance_tolerance, seed_base + 3)
      target(:, 3) = method_weight
      call po_avg(zero, train, 'RIDGE', method_weight, istat, opt%variance_tolerance, seed_base + 4)
      target(:, 4) = method_weight
      call po_gross_exp(zero, train, 'NOSHORT', method_weight, istat, &
        opt%variance_tolerance, seed_base + 5)
      target(:, 5) = method_weight
      call po_gross_exp(zero, train, 'EQUAL', method_weight, istat, &
        opt%variance_tolerance, seed_base + 6)
      target(:, 6) = method_weight
      target(:, 7) = 1.0_dp / real(p, dp)
      call po_cov_shrink(zero, train, method_weight, istat, opt%variance_tolerance)
      target(:, 8) = method_weight
      call po_bhu(zero, train, clusters%groups, opt%cluster_repetitions, method_weight, istat, &
        opt%variance_tolerance, seed_base + 9)
        target(:, 9) = method_weight
      call po_tzt(train, 3.0_dp, method_weight, istat, opt%variance_tolerance)
      target(:, 10) = method_weight
      b_assets = max(1, nint(real(p, dp) ** 0.7_dp))
      call po_sw(train, b_assets, opt%stochastic_samples, method_weight, istat, &
        opt%variance_tolerance, seed_base + 11)
      target(:, 11) = method_weight
      call po_sw(train, max(1, size(clusters%groups)), opt%stochastic_samples, method_weight, istat, &
        opt%variance_tolerance, seed_base + 12)
        target(:, 12) = method_weight
      call po_sw_lasso(zero, train, max(1, size(clusters%groups)), opt%stochastic_samples, method_weight, istat, &
        opt%variance_tolerance, seed_base + 13)
        target(:, 13) = method_weight
      result%weights(eval_index, :, :) = target
      do method_id = 1, ren_method_count
        result%turnover(eval_index, method_id) = l1_norm(target(:, method_id) - current(:, method_id))
      end do
      result%month_date(eval_index) = date(current_rows(1))
      vw_target = vw_current
      result%vw_turnover(eval_index) = 0.0_dp
      if (.not. opt%legacy_weight_timing) current = target
      do row = 1, size(current_rows)
        return_index = return_index + 1
        gross = 1.0_dp + x(current_rows(row), :) / 100.0_dp
        result%return_date(return_index) = date(current_rows(row))
        do method_id = 1, ren_method_count
          result%gross_returns(return_index, method_id) = dot_product(gross, current(:, method_id))
          if (abs(result%gross_returns(return_index, method_id)) > tiny(1.0_dp)) then
            current(:, method_id) = gross * current(:, method_id) / result%gross_returns(return_index, method_id)
          end if
        end do
        result%vw_gross_returns(return_index) = dot_product(gross, vw_current)
        if (abs(result%vw_gross_returns(return_index)) > tiny(1.0_dp)) &
          vw_current = gross * vw_current / result%vw_gross_returns(return_index)
      end do
      result%vw_weights(eval_index, :) = vw_current
      result%vw_turnover(eval_index) = l1_norm(vw_current - vw_target)
      if (opt%legacy_weight_timing) current = target
    end do
    do method_id = 1, ren_method_count
      call cumulative_metrics(result%gross_returns(:, method_id), result%cumulative_return(:, method_id), &
        result%wealth_index(:, method_id), result%sharpe_ratio(method_id), result%volatility(method_id), &
        result%max_drawdown(method_id))
      result%turnover_mean(method_id) = sum(result%turnover(:, method_id)) / real(eval_months, dp) * 100.0_dp
      result%cumulative_turnover(1, method_id) = abs(result%turnover(1, method_id))
      do eval_index = 2, eval_months
        result%cumulative_turnover(eval_index, method_id) = result%cumulative_turnover(eval_index - 1, method_id) + &
          abs(result%turnover(eval_index, method_id))
      end do
    end do
    call cumulative_metrics(result%vw_gross_returns, result%vw_cumulative_return, result%vw_wealth_index, &
      result%vw_sharpe_ratio, result%vw_volatility, result%vw_max_drawdown)
    result%vw_to_mean = sum(result%vw_turnover) / real(eval_months, dp) * 100.0_dp
    result%status = ren_success
  end subroutine perform_analysis

  subroutine ren_run(date, x, result, start_date, end_date, missing_value, options)
    integer, intent(in) :: date(:)
    real(dp), intent(in) :: x(:, :)
    type(analysis_result), intent(out) :: result
    integer, intent(in), optional :: start_date, end_date
    real(dp), intent(in), optional :: missing_value
    type(analysis_options), intent(in), optional :: options
    type(prepared_data_type) :: prepared
    prepared = prepare_data(date, x, start_date, end_date, missing_value)
    if (prepared%status /= ren_success) then
      result%status = prepared%status
      return
    end if
    if (present(options)) then
      call perform_analysis(prepared%x, prepared%month, prepared%date, result, options)
    else
      call perform_analysis(prepared%x, prepared%month, prepared%date, result)
    end if
  end subroutine ren_run
end module ren_analysis
