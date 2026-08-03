! SPDX-License-Identifier: GPL-3.0-only
module bondvaluation_schedule
  use bondvaluation_kinds, only: dp
  use bondvaluation_dates, only: date_type, valid_date, compare_dates, add_months, &
    day_diff, is_last_day_of_month, unique_sorted_dates, find_previous_date, &
    find_next_date, date_to_serial
  use bondvaluation_daycount, only: day_count_fraction, daycount_result, &
    dcc_bus_252
  implicit none
  private

  integer, parameter, public :: bondvaluation_ok = 0
  integer, parameter, public :: bondvaluation_invalid_date = 1
  integer, parameter, public :: bondvaluation_invalid_frequency = 2
  integer, parameter, public :: bondvaluation_invalid_order = 3
  integer, parameter, public :: bondvaluation_no_schedule = 4

  type, public :: bond_terms
    type(date_type) :: issue_date
    type(date_type) :: maturity_date
    integer :: coupon_frequency = 2
    type(date_type) :: first_coupon_date
    type(date_type) :: last_coupon_date
    type(date_type) :: first_interest_accrual_date
    real(dp) :: redemption_value = 100.0_dp
    real(dp) :: coupon_rate_percent = 0.0_dp
    integer :: day_count_convention = 2
    logical :: end_of_month = .false.
    logical :: infer_end_of_month = .true.
    logical :: regular_cashflows_equal = .false.
  end type bond_terms

  type, public :: schedule_warnings
    logical :: issue_replaced_by_accrual = .false.
    logical :: zero_coupon = .false.
    logical :: first_coupon_dropped = .false.
    logical :: last_coupon_dropped = .false.
    logical :: inconsistent_coupon_dates = .false.
    logical :: end_of_month_overridden = .false.
  end type schedule_warnings

  type, public :: bond_schedule
    type(date_type), allocatable :: real_dates(:)
    real(dp), allocatable :: real_indexes(:)
    type(date_type), allocatable :: coupon_dates(:)
    real(dp), allocatable :: coupon_indexes(:)
    type(date_type), allocatable :: anniversary_dates(:)
    real(dp), allocatable :: anniversary_indexes(:)
    real(dp), allocatable :: coupon_payments(:)
    type(date_type) :: accrual_start
    type(date_type) :: first_coupon_used
    type(date_type) :: last_coupon_used
    type(date_type) :: estimated_first_coupon
    type(date_type) :: estimated_last_coupon
    type(date_type) :: reference_date
    logical :: end_of_month_used = .false.
    character(len=8) :: first_coupon_type = 'NA'
    character(len=8) :: last_coupon_type = 'NA'
    real(dp) :: first_coupon_length = 0.0_dp
    real(dp) :: last_coupon_length = 0.0_dp
    integer :: status = 0
    type(schedule_warnings) :: warnings
  end type bond_schedule

  public :: anniv_dates, build_bond_schedule, period_coordinate
  public :: previous_coupon_date, next_coupon_date

contains

  subroutine anniv_dates(terms, schedule)
    type(bond_terms), intent(in) :: terms
    type(bond_schedule), intent(out) :: schedule
    call build_bond_schedule(terms, schedule)
  end subroutine anniv_dates

  subroutine build_bond_schedule(terms, schedule)
    type(bond_terms), intent(in) :: terms
    type(bond_schedule), intent(out) :: schedule
    type(date_type) :: start_date, anchor, current
    type(date_type), allocatable :: candidates(:), coupons(:), work(:)
    integer :: cpy, months, total_months, nmax, n, i, k
    logical :: first_valid, last_valid, eom
    type(daycount_result) :: dc
    real(dp) :: length_value, regular_payment

    schedule = bond_schedule()
    if (.not. valid_date(terms%issue_date) .or. &
        .not. valid_date(terms%maturity_date)) then
      schedule%status = bondvaluation_invalid_date
      return
    end if
    if (compare_dates(terms%issue_date, terms%maturity_date) >= 0) then
      schedule%status = bondvaluation_invalid_order
      return
    end if
    cpy = terms%coupon_frequency
    if (.not. any(cpy == [0, 1, 2, 3, 4, 6, 12])) then
      schedule%status = bondvaluation_invalid_frequency
      return
    end if
    start_date = terms%issue_date
    if (valid_date(terms%first_interest_accrual_date)) then
      start_date = terms%first_interest_accrual_date
      schedule%warnings%issue_replaced_by_accrual = &
        compare_dates(start_date, terms%issue_date) /= 0
    end if
    if (compare_dates(start_date, terms%maturity_date) >= 0) then
      schedule%status = bondvaluation_invalid_order
      return
    end if
    schedule%accrual_start = start_date

    if (cpy == 0) then
      schedule%warnings%zero_coupon = .true.
      schedule%end_of_month_used = terms%end_of_month
      allocate(schedule%coupon_dates(1), schedule%coupon_payments(1))
      schedule%coupon_dates(1) = terms%maturity_date
      schedule%coupon_payments(1) = 0.0_dp
      allocate(schedule%real_dates(2), schedule%real_indexes(2))
      schedule%real_dates = [start_date, terms%maturity_date]
      schedule%real_indexes = [0.0_dp, 1.0_dp]
      allocate(schedule%coupon_indexes(1))
      schedule%coupon_indexes = [1.0_dp]
      allocate(schedule%anniversary_dates(2), schedule%anniversary_indexes(2))
      schedule%anniversary_dates = schedule%real_dates
      schedule%anniversary_indexes = schedule%real_indexes
      schedule%reference_date = terms%maturity_date
      return
    end if

    months = 12 / cpy
    first_valid = valid_date(terms%first_coupon_date) .and. &
      compare_dates(start_date, terms%first_coupon_date) < 0 .and. &
      compare_dates(terms%first_coupon_date, terms%maturity_date) < 0
    last_valid = valid_date(terms%last_coupon_date) .and. &
      compare_dates(start_date, terms%last_coupon_date) < 0 .and. &
      compare_dates(terms%last_coupon_date, terms%maturity_date) < 0
    if (first_valid .and. last_valid) then
      if (compare_dates(terms%first_coupon_date, terms%last_coupon_date) > 0) then
        first_valid = .false.
        last_valid = .false.
        schedule%warnings%inconsistent_coupon_dates = .true.
      end if
    end if
    if (valid_date(terms%first_coupon_date) .and. .not. first_valid) &
      schedule%warnings%first_coupon_dropped = .true.
    if (valid_date(terms%last_coupon_date) .and. .not. last_valid) &
      schedule%warnings%last_coupon_dropped = .true.

    eom = terms%end_of_month
    if (terms%infer_end_of_month) then
      if (first_valid) then
        eom = is_last_day_of_month(terms%first_coupon_date)
      else if (last_valid) then
        eom = is_last_day_of_month(terms%last_coupon_date)
      else
        eom = is_last_day_of_month(terms%maturity_date)
      end if
      schedule%warnings%end_of_month_overridden = eom .neqv. terms%end_of_month
    end if
    schedule%end_of_month_used = eom

    if (first_valid) then
      anchor = terms%first_coupon_date
    else if (last_valid) then
      anchor = terms%last_coupon_date
    else
      anchor = terms%maturity_date
    end if
    schedule%reference_date = anchor

    total_months = abs((terms%maturity_date%year - start_date%year) * 12 + &
      terms%maturity_date%month - start_date%month)
    nmax = total_months / months + 16
    allocate(candidates(2 * nmax + 1))
    n = 0
    do k = -nmax, nmax
      current = add_months(anchor, k * months, eom, anchor%day)
      if (compare_dates(current, add_months(start_date, -2 * months, eom, anchor%day)) >= 0 .and. &
          compare_dates(current, add_months(terms%maturity_date, 2 * months, eom, anchor%day)) <= 0) then
        n = n + 1
        candidates(n) = current
      end if
    end do
    call unique_sorted_dates(candidates(1:n), schedule%anniversary_dates)
    allocate(schedule%anniversary_indexes(size(schedule%anniversary_dates)))
    do i = 1, size(schedule%anniversary_dates)
      schedule%anniversary_indexes(i) = real(i - 1, dp)
    end do

    allocate(work(nmax + 4))
    n = 0
    if (first_valid) then
      current = terms%first_coupon_date
      do while (compare_dates(current, terms%maturity_date) < 0)
        if (last_valid .and. compare_dates(current, terms%last_coupon_date) > 0) exit
        n = n + 1
        work(n) = current
        if (last_valid .and. compare_dates(current, terms%last_coupon_date) == 0) exit
        current = add_months(terms%first_coupon_date, n * months, eom, anchor%day)
        if (n >= size(work) - 2) exit
      end do
      if (last_valid) then
        if (n == 0 .or. compare_dates(work(n), terms%last_coupon_date) /= 0) then
          n = n + 1
          work(n) = terms%last_coupon_date
        end if
      end if
    else if (last_valid) then
      k = 0
      do
        current = add_months(terms%last_coupon_date, -k * months, eom, anchor%day)
        if (compare_dates(current, start_date) <= 0) exit
        n = n + 1
        work(n) = current
        k = k + 1
        if (n >= size(work) - 2) exit
      end do
    else
      do i = 1, size(schedule%anniversary_dates)
        if (compare_dates(schedule%anniversary_dates(i), start_date) > 0 .and. &
            compare_dates(schedule%anniversary_dates(i), terms%maturity_date) < 0) then
          n = n + 1
          work(n) = schedule%anniversary_dates(i)
        end if
      end do
    end if
    n = n + 1
    work(n) = terms%maturity_date
    call unique_sorted_dates(work(1:n), coupons)
    allocate(schedule%coupon_dates(size(coupons)))
    schedule%coupon_dates = coupons
    schedule%first_coupon_used = schedule%coupon_dates(1)
    if (size(schedule%coupon_dates) > 1) then
      schedule%last_coupon_used = schedule%coupon_dates(size(schedule%coupon_dates) - 1)
    else
      schedule%last_coupon_used = date_type()
    end if
    if (.not. first_valid) schedule%estimated_first_coupon = schedule%first_coupon_used
    if (.not. last_valid .and. size(schedule%coupon_dates) > 1) &
      schedule%estimated_last_coupon = schedule%last_coupon_used

    allocate(schedule%real_dates(size(schedule%coupon_dates) + 1))
    schedule%real_dates(1) = start_date
    schedule%real_dates(2:) = schedule%coupon_dates
    allocate(schedule%real_indexes(size(schedule%real_dates)))
    allocate(schedule%coupon_indexes(size(schedule%coupon_dates)))
    do i = 1, size(schedule%real_dates)
      schedule%real_indexes(i) = period_coordinate(schedule%real_dates(i), &
        schedule%anniversary_dates)
    end do
    schedule%coupon_indexes = schedule%real_indexes(2:)

    schedule%first_coupon_length = schedule%real_indexes(2) - schedule%real_indexes(1)
    schedule%last_coupon_length = schedule%real_indexes(size(schedule%real_indexes)) - &
      schedule%real_indexes(size(schedule%real_indexes) - 1)
    call classify_period(schedule%first_coupon_length, schedule%first_coupon_type)
    call classify_period(schedule%last_coupon_length, schedule%last_coupon_type)

    allocate(schedule%coupon_payments(size(schedule%coupon_dates)))
    regular_payment = terms%redemption_value * terms%coupon_rate_percent / &
      (100.0_dp * real(cpy, dp))
    do i = 1, size(schedule%coupon_dates)
      dc = day_count_fraction(schedule%real_dates(i), schedule%real_dates(i + 1), &
        terms%day_count_convention, cpy, terms%maturity_date, eom, &
        schedule%coupon_dates(i)%year, schedule%anniversary_dates)
      if (terms%day_count_convention == dcc_bus_252) then
        schedule%coupon_payments(i) = terms%redemption_value * &
          ((1.0_dp + terms%coupon_rate_percent / 100.0_dp) ** dc%fraction - 1.0_dp)
      else
        schedule%coupon_payments(i) = terms%redemption_value * &
          terms%coupon_rate_percent / 100.0_dp * dc%fraction
      end if
      length_value = schedule%real_indexes(i + 1) - schedule%real_indexes(i)
      if (terms%regular_cashflows_equal .and. abs(length_value - 1.0_dp) < 1.0e-10_dp .and. &
          terms%day_count_convention /= dcc_bus_252) &
        schedule%coupon_payments(i) = regular_payment
    end do
  end subroutine build_bond_schedule

  subroutine classify_period(length_value, label)
    real(dp), intent(in) :: length_value
    character(len=*), intent(out) :: label
    if (abs(length_value - 1.0_dp) <= 1.0e-8_dp) then
      label = 'regular'
    else if (length_value < 1.0_dp) then
      label = 'short'
    else
      label = 'long'
    end if
  end subroutine classify_period

  real(dp) function period_coordinate(date, anniversary_dates) result(value)
    type(date_type), intent(in) :: date
    type(date_type), intent(in) :: anniversary_dates(:)
    integer :: ip, inext, numerator, denominator
    value = 0.0_dp
    if (size(anniversary_dates) < 2) return
    ip = find_previous_date(date, anniversary_dates)
    if (ip == 0) then
      ip = 1
    else if (ip == size(anniversary_dates)) then
      value = real(ip - 1, dp)
      return
    end if
    inext = ip + 1
    numerator = day_diff(anniversary_dates(ip), date)
    denominator = day_diff(anniversary_dates(ip), anniversary_dates(inext))
    value = real(ip - 1, dp) + real(numerator, dp) / real(denominator, dp)
  end function period_coordinate

  function previous_coupon_date(date, schedule) result(out)
    type(date_type), intent(in) :: date
    type(bond_schedule), intent(in) :: schedule
    type(date_type) :: out
    integer :: i
    i = find_previous_date(date, schedule%real_dates)
    if (i > 0) then
      out = schedule%real_dates(i)
    else
      out = date_type()
    end if
  end function previous_coupon_date

  function next_coupon_date(date, schedule) result(out)
    type(date_type), intent(in) :: date
    type(bond_schedule), intent(in) :: schedule
    type(date_type) :: out
    integer :: i
    i = find_next_date(date, schedule%real_dates)
    if (i > 0) then
      out = schedule%real_dates(i)
    else
      out = date_type()
    end if
  end function next_coupon_date

end module bondvaluation_schedule
