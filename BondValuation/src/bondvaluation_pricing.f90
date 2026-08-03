! SPDX-License-Identifier: GPL-3.0-only
module bondvaluation_pricing
  use ieee_arithmetic, only: ieee_value, ieee_quiet_nan
  use bondvaluation_kinds, only: dp
  use bondvaluation_dates, only: date_type, valid_date, compare_dates, day_diff, &
    find_previous_date, find_next_date
  use bondvaluation_daycount, only: daycount_result, day_count_fraction, &
    year_fraction, dcc_act_act_icma, dcc_bus_252
  use bondvaluation_schedule, only: bond_terms, bond_schedule, build_bond_schedule, &
    previous_coupon_date, next_coupon_date
  implicit none
  private

  type, public :: accrued_interest_result
    real(dp) :: accrued_interest = 0.0_dp
    integer :: days_accrued = 0
    real(dp) :: year_fraction_value = 0.0_dp
    integer :: status = 0
  end type accrued_interest_result

  type, public :: dirty_price_result
    real(dp) :: dirty_price = 0.0_dp
    real(dp) :: clean_price = 0.0_dp
    real(dp) :: accrued_interest = 0.0_dp
    real(dp) :: coupon_payment = 0.0_dp
    type(date_type) :: previous_coupon
    type(date_type) :: settlement_date
    type(date_type) :: next_coupon
    integer :: days_accrued = 0
    integer :: days_in_period = 0
    integer :: status = 0
  end type dirty_price_result

  type, public :: bond_value_result
    real(dp) :: clean_price = 0.0_dp
    real(dp) :: accrued_interest = 0.0_dp
    real(dp) :: dirty_price = 0.0_dp
    real(dp) :: yield_percent = 0.0_dp
    real(dp) :: modified_duration_years = 0.0_dp
    real(dp) :: macaulay_duration_years = 0.0_dp
    real(dp) :: convexity_years = 0.0_dp
    real(dp) :: modified_duration_periods = 0.0_dp
    real(dp) :: macaulay_duration_periods = 0.0_dp
    real(dp) :: convexity_periods = 0.0_dp
    real(dp) :: tau = 0.0_dp
    integer :: iterations = 0
    integer :: status = 0
  end type bond_value_result

  public :: accrued_interest, accr_int, dirty_price, dp_value
  public :: bond_price, bond_val_price, bond_yield, bond_val_yield
  public :: discount_coordinate

contains

  function accrued_interest(start_date, end_date, coupon_rate_percent, dcc, &
                            redemption_value, coupon_frequency, maturity, &
                            next_coupon_year, end_of_month, anniversary_dates) result(out)
    type(date_type), intent(in) :: start_date, end_date
    real(dp), intent(in) :: coupon_rate_percent
    integer, intent(in) :: dcc
    real(dp), intent(in), optional :: redemption_value
    integer, intent(in), optional :: coupon_frequency
    type(date_type), intent(in), optional :: maturity
    integer, intent(in), optional :: next_coupon_year
    logical, intent(in), optional :: end_of_month
    type(date_type), intent(in), optional :: anniversary_dates(:)
    type(accrued_interest_result) :: out
    real(dp) :: rv
    integer :: cpy
    type(date_type) :: mat
    logical :: eom
    integer :: ncy
    type(daycount_result) :: dc

    rv = 100.0_dp
    if (present(redemption_value)) rv = redemption_value
    cpy = 2
    if (present(coupon_frequency)) cpy = coupon_frequency
    mat = end_date
    if (present(maturity)) mat = maturity
    ncy = end_date%year
    if (present(next_coupon_year)) ncy = next_coupon_year
    eom = .false.
    if (present(end_of_month)) eom = end_of_month
    if (present(anniversary_dates)) then
      dc = day_count_fraction(start_date, end_date, dcc, cpy, mat, eom, ncy, anniversary_dates)
    else
      dc = day_count_fraction(start_date, end_date, dcc, cpy, mat, eom, ncy)
    end if
    out%days_accrued = dc%days_accrued
    out%year_fraction_value = dc%fraction
    out%status = dc%status
    if (dcc == dcc_bus_252) then
      out%accrued_interest = rv * ((1.0_dp + coupon_rate_percent / 100.0_dp) ** &
        dc%fraction - 1.0_dp)
    else
      out%accrued_interest = rv * coupon_rate_percent / 100.0_dp * dc%fraction
    end if
  end function accrued_interest

  function accr_int(start_date, end_date, coupon_rate_percent, dcc, &
                    redemption_value, coupon_frequency, maturity, &
                    next_coupon_year, end_of_month, anniversary_dates) result(out)
    type(date_type), intent(in) :: start_date, end_date
    real(dp), intent(in) :: coupon_rate_percent
    integer, intent(in) :: dcc
    real(dp), intent(in), optional :: redemption_value
    integer, intent(in), optional :: coupon_frequency
    type(date_type), intent(in), optional :: maturity
    integer, intent(in), optional :: next_coupon_year
    logical, intent(in), optional :: end_of_month
    type(date_type), intent(in), optional :: anniversary_dates(:)
    type(accrued_interest_result) :: out
    if (present(anniversary_dates)) then
      out = accrued_interest(start_date, end_date, coupon_rate_percent, dcc, &
        redemption_value, coupon_frequency, maturity, next_coupon_year, &
        end_of_month, anniversary_dates)
    else
      out = accrued_interest(start_date, end_date, coupon_rate_percent, dcc, &
        redemption_value, coupon_frequency, maturity, next_coupon_year, end_of_month)
    end if
  end function accr_int

  function dirty_price(clean_price, settlement_date, terms, supplied_schedule) result(out)
    real(dp), intent(in) :: clean_price
    type(date_type), intent(in) :: settlement_date
    type(bond_terms), intent(in) :: terms
    type(bond_schedule), intent(in), optional :: supplied_schedule
    type(dirty_price_result) :: out
    type(bond_schedule) :: schedule
    type(date_type) :: previous_date_value, next_date_value
    type(accrued_interest_result) :: ai
    integer :: inext, coupon_index

    if (present(supplied_schedule)) then
      schedule = supplied_schedule
    else
      call build_bond_schedule(terms, schedule)
    end if
    out%clean_price = clean_price
    out%settlement_date = settlement_date
    if (schedule%status /= 0) then
      out%status = schedule%status
      return
    end if
    if (compare_dates(settlement_date, schedule%real_dates(1)) < 0 .or. &
        compare_dates(settlement_date, terms%maturity_date) >= 0) then
      out%status = 10
      return
    end if
    previous_date_value = previous_coupon_date(settlement_date, schedule)
    next_date_value = next_coupon_date(settlement_date, schedule)
    out%previous_coupon = previous_date_value
    out%next_coupon = next_date_value
    if (schedule%warnings%zero_coupon) then
      out%accrued_interest = 0.0_dp
      out%dirty_price = clean_price
      out%coupon_payment = 0.0_dp
      return
    end if
    ai = accrued_interest(previous_date_value, settlement_date, &
      terms%coupon_rate_percent, terms%day_count_convention, &
      terms%redemption_value, terms%coupon_frequency, terms%maturity_date, &
      next_date_value%year, schedule%end_of_month_used, schedule%anniversary_dates)
    out%accrued_interest = ai%accrued_interest
    out%days_accrued = ai%days_accrued
    out%days_in_period = day_diff(previous_date_value, next_date_value)
    out%dirty_price = clean_price + out%accrued_interest
    inext = find_next_date(settlement_date, schedule%real_dates)
    coupon_index = inext - 1
    if (coupon_index >= 1 .and. coupon_index <= size(schedule%coupon_payments)) &
      out%coupon_payment = schedule%coupon_payments(coupon_index)
  end function dirty_price

  function dp_value(clean_price, settlement_date, terms, supplied_schedule) result(out)
    real(dp), intent(in) :: clean_price
    type(date_type), intent(in) :: settlement_date
    type(bond_terms), intent(in) :: terms
    type(bond_schedule), intent(in), optional :: supplied_schedule
    type(dirty_price_result) :: out
    if (present(supplied_schedule)) then
      out = dirty_price(clean_price, settlement_date, terms, supplied_schedule)
    else
      out = dirty_price(clean_price, settlement_date, terms)
    end if
  end function dp_value

  real(dp) function discount_coordinate(date, dcc, terms, schedule) result(value)
    type(date_type), intent(in) :: date
    integer, intent(in) :: dcc
    type(bond_terms), intent(in) :: terms
    type(bond_schedule), intent(in) :: schedule
    integer :: ip, inext
    real(dp) :: numerator, denominator
    value = 0.0_dp
    ip = find_previous_date(date, schedule%anniversary_dates)
    if (ip == 0) ip = 1
    if (ip >= size(schedule%anniversary_dates)) then
      value = real(size(schedule%anniversary_dates) - 1, dp)
      return
    end if
    inext = ip + 1
    numerator = year_fraction(schedule%anniversary_dates(ip), date, dcc, &
      terms%coupon_frequency, terms%maturity_date, schedule%end_of_month_used, &
      schedule%anniversary_dates(inext)%year, schedule%anniversary_dates)
    denominator = year_fraction(schedule%anniversary_dates(ip), &
      schedule%anniversary_dates(inext), dcc, terms%coupon_frequency, &
      terms%maturity_date, schedule%end_of_month_used, &
      schedule%anniversary_dates(inext)%year, schedule%anniversary_dates)
    if (abs(denominator) <= tiny(1.0_dp)) then
      value = real(ip - 1, dp)
    else
      value = real(ip - 1, dp) + numerator / denominator
    end if
  end function discount_coordinate

  function bond_price(yield_percent, settlement_date, terms, simple_last_period, &
                      calculation_method, supplied_schedule) result(out)
    real(dp), intent(in) :: yield_percent
    type(date_type), intent(in) :: settlement_date
    type(bond_terms), intent(in) :: terms
    logical, intent(in), optional :: simple_last_period
    integer, intent(in), optional :: calculation_method
    type(bond_schedule), intent(in), optional :: supplied_schedule
    type(bond_value_result) :: out
    type(bond_schedule) :: schedule
    type(dirty_price_result) :: dp_result
    logical :: simple_last
    integer :: calc_method, dcc_discount, first_remaining, nremain, i, cpy
    real(dp), allocatable :: cashflows(:), powers(:)
    real(dp) :: a, derivative, second_derivative, settlement_coordinate

    simple_last = .true.
    if (present(simple_last_period)) simple_last = simple_last_period
    calc_method = 1
    if (present(calculation_method)) calc_method = calculation_method
    if (present(supplied_schedule)) then
      schedule = supplied_schedule
    else
      call build_bond_schedule(terms, schedule)
    end if
    if (schedule%status /= 0) then
      out%status = schedule%status
      return
    end if
    if (compare_dates(settlement_date, schedule%real_dates(1)) < 0 .or. &
        compare_dates(settlement_date, terms%maturity_date) >= 0) then
      out%status = 10
      return
    end if
    dp_result = dirty_price(0.0_dp, settlement_date, terms, schedule)
    out%accrued_interest = dp_result%accrued_interest
    out%yield_percent = yield_percent
    cpy = max(terms%coupon_frequency, 1)
    dcc_discount = terms%day_count_convention
    if (calc_method == 0) dcc_discount = dcc_act_act_icma
    settlement_coordinate = discount_coordinate(settlement_date, dcc_discount, terms, schedule)
    out%tau = settlement_coordinate
    first_remaining = find_next_date(settlement_date, schedule%coupon_dates)
    if (first_remaining == 0) then
      out%status = 11
      return
    end if
    nremain = size(schedule%coupon_dates) - first_remaining + 1
    allocate(cashflows(nremain), powers(nremain))
    do i = 1, nremain
      cashflows(i) = schedule%coupon_payments(first_remaining + i - 1)
      if (first_remaining + i - 1 == size(schedule%coupon_dates)) &
        cashflows(i) = cashflows(i) + terms%redemption_value
      powers(i) = discount_coordinate(schedule%coupon_dates(first_remaining + i - 1), &
        dcc_discount, terms, schedule) - settlement_coordinate
    end do

    if (nremain == 1 .and. simple_last) then
      a = 1.0_dp + yield_percent * powers(1) / (100.0_dp * real(cpy, dp))
      if (a <= 0.0_dp) then
        out%status = 12
        return
      end if
      out%dirty_price = cashflows(1) / a
      out%modified_duration_years = cashflows(1) * powers(1) / &
        (real(cpy, dp) * a * a * out%dirty_price)
      out%macaulay_duration_years = out%modified_duration_years * a
      out%convexity_years = cashflows(1) * (powers(1) / real(cpy, dp)) ** 2 / &
        (a ** 3 * out%dirty_price)
      out%modified_duration_periods = out%modified_duration_years * real(cpy, dp)
      out%macaulay_duration_periods = out%macaulay_duration_years * real(cpy, dp)
      out%convexity_periods = out%convexity_years * real(cpy * cpy, dp)
    else
      a = 1.0_dp + yield_percent / (100.0_dp * real(cpy, dp))
      if (a <= 0.0_dp) then
        out%status = 12
        return
      end if
      out%dirty_price = sum(cashflows / a ** powers)
      derivative = sum(-cashflows * powers / a ** (powers + 1.0_dp))
      second_derivative = sum(cashflows * powers * (powers + 1.0_dp) / &
        a ** (powers + 2.0_dp))
      out%modified_duration_periods = -derivative / out%dirty_price
      out%macaulay_duration_periods = out%modified_duration_periods * a
      out%convexity_periods = 0.5_dp * second_derivative / out%dirty_price
      out%modified_duration_years = out%modified_duration_periods / real(cpy, dp)
      out%macaulay_duration_years = out%macaulay_duration_periods / real(cpy, dp)
      out%convexity_years = out%convexity_periods / real(cpy * cpy, dp)
    end if
    out%clean_price = out%dirty_price - out%accrued_interest
  end function bond_price

  function bond_val_price(yield_percent, settlement_date, terms, &
                          simple_last_period, calculation_method, supplied_schedule) result(out)
    real(dp), intent(in) :: yield_percent
    type(date_type), intent(in) :: settlement_date
    type(bond_terms), intent(in) :: terms
    logical, intent(in), optional :: simple_last_period
    integer, intent(in), optional :: calculation_method
    type(bond_schedule), intent(in), optional :: supplied_schedule
    type(bond_value_result) :: out
    if (present(supplied_schedule)) then
      out = bond_price(yield_percent, settlement_date, terms, simple_last_period, &
        calculation_method, supplied_schedule)
    else
      out = bond_price(yield_percent, settlement_date, terms, simple_last_period, calculation_method)
    end if
  end function bond_val_price

  function bond_yield(clean_price, settlement_date, terms, simple_last_period, &
                      precision, calculation_method, supplied_schedule) result(out)
    real(dp), intent(in) :: clean_price
    type(date_type), intent(in) :: settlement_date
    type(bond_terms), intent(in) :: terms
    logical, intent(in), optional :: simple_last_period
    real(dp), intent(in), optional :: precision
    integer, intent(in), optional :: calculation_method
    type(bond_schedule), intent(in), optional :: supplied_schedule
    type(bond_value_result) :: out
    type(bond_schedule) :: schedule
    type(dirty_price_result) :: dp_result
    type(bond_value_result) :: trial
    logical :: simple_last
    real(dp) :: tol, target_dirty, low, high, mid, f_low, f_high, f_mid
    integer :: calc_method, iter, first_remaining

    simple_last = .true.
    if (present(simple_last_period)) simple_last = simple_last_period
    tol = epsilon(1.0_dp) ** 0.75_dp
    if (present(precision)) tol = precision
    calc_method = 1
    if (present(calculation_method)) calc_method = calculation_method
    if (present(supplied_schedule)) then
      schedule = supplied_schedule
    else
      call build_bond_schedule(terms, schedule)
    end if
    if (schedule%status /= 0) then
      out%status = schedule%status
      return
    end if
    dp_result = dirty_price(clean_price, settlement_date, terms, schedule)
    if (dp_result%status /= 0) then
      out%status = dp_result%status
      return
    end if
    target_dirty = dp_result%dirty_price
    first_remaining = find_next_date(settlement_date, schedule%coupon_dates)
    if (first_remaining == size(schedule%coupon_dates) .and. simple_last) then
      trial = bond_price(0.0_dp, settlement_date, terms, simple_last, calc_method, schedule)
      if (trial%status /= 0) then
        out%status = trial%status
        return
      end if
      ! The final cash flow is known from the zero-yield price.
      mid = (trial%dirty_price / target_dirty - 1.0_dp) * &
        100.0_dp * real(max(terms%coupon_frequency, 1), dp) / &
        max(trial%modified_duration_periods, tiny(1.0_dp))
      out = bond_price(mid, settlement_date, terms, simple_last, &
        calc_method, schedule)
      out%yield_percent = mid
      out%iterations = 0
      return
    end if

    low = -99.999999_dp * real(max(terms%coupon_frequency, 1), dp)
    high = 100.0_dp
    trial = bond_price(low, settlement_date, terms, simple_last, calc_method, schedule)
    f_low = trial%dirty_price - target_dirty
    trial = bond_price(high, settlement_date, terms, simple_last, calc_method, schedule)
    f_high = trial%dirty_price - target_dirty
    do while (f_low * f_high > 0.0_dp .and. high < 1.0e7_dp)
      high = high * 2.0_dp
      trial = bond_price(high, settlement_date, terms, simple_last, calc_method, schedule)
      f_high = trial%dirty_price - target_dirty
    end do
    if (f_low * f_high > 0.0_dp) then
      out%status = 20
      return
    end if
    do iter = 1, 300
      mid = 0.5_dp * (low + high)
      trial = bond_price(mid, settlement_date, terms, simple_last, calc_method, schedule)
      f_mid = trial%dirty_price - target_dirty
      if (abs(f_mid) <= tol * max(1.0_dp, target_dirty) .or. &
          abs(high - low) <= tol * max(1.0_dp, abs(mid))) exit
      if (f_low * f_mid <= 0.0_dp) then
        high = mid
        f_high = f_mid
      else
        low = mid
        f_low = f_mid
      end if
    end do
    out = bond_price(mid, settlement_date, terms, simple_last, calc_method, schedule)
    out%clean_price = clean_price
    out%dirty_price = target_dirty
    out%yield_percent = mid
    out%iterations = iter
  end function bond_yield

  function bond_val_yield(clean_price, settlement_date, terms, &
                          simple_last_period, precision, calculation_method, &
                          supplied_schedule) result(out)
    real(dp), intent(in) :: clean_price
    type(date_type), intent(in) :: settlement_date
    type(bond_terms), intent(in) :: terms
    logical, intent(in), optional :: simple_last_period
    real(dp), intent(in), optional :: precision
    integer, intent(in), optional :: calculation_method
    type(bond_schedule), intent(in), optional :: supplied_schedule
    type(bond_value_result) :: out
    if (present(supplied_schedule)) then
      out = bond_yield(clean_price, settlement_date, terms, simple_last_period, &
        precision, calculation_method, supplied_schedule)
    else
      out = bond_yield(clean_price, settlement_date, terms, simple_last_period, &
        precision, calculation_method)
    end if
  end function bond_val_yield

end module bondvaluation_pricing
