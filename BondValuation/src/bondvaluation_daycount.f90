! SPDX-License-Identifier: GPL-3.0-only
module bondvaluation_daycount
  use bondvaluation_kinds, only: dp
  use bondvaluation_dates, only: date_type, valid_date, leap_year, days_in_year, &
    days_in_month, date_to_serial, serial_to_date, day_diff, compare_dates, &
    add_months, is_last_day_of_month, find_previous_date, find_next_date
  use bondvaluation_brazil_calendar, only: count_nonbusiness_brazil
  implicit none
  private

  integer, parameter, public :: dcc_act_act_isda = 1
  integer, parameter, public :: dcc_act_act_icma = 2
  integer, parameter, public :: dcc_act_act_afb = 3
  integer, parameter, public :: dcc_act_365l = 4
  integer, parameter, public :: dcc_30_360_bond = 5
  integer, parameter, public :: dcc_30e_360 = 6
  integer, parameter, public :: dcc_30e_360_isda = 7
  integer, parameter, public :: dcc_30_360_german = 8
  integer, parameter, public :: dcc_30u_360 = 9
  integer, parameter, public :: dcc_act_365 = 10
  integer, parameter, public :: dcc_act_nl_365 = 11
  integer, parameter, public :: dcc_act_360 = 12
  integer, parameter, public :: dcc_30_365 = 13
  integer, parameter, public :: dcc_act_365_canadian = 14
  integer, parameter, public :: dcc_act_364 = 15
  integer, parameter, public :: dcc_bus_252 = 16

  type, public :: daycount_result
    integer :: days_accrued = 0
    real(dp) :: fraction = 0.0_dp
    integer :: status = 0
  end type daycount_result

  public :: day_count_fraction, year_fraction, day_count_name
  public :: leap_day_inside

contains

  function day_count_name(dcc) result(name)
    integer, intent(in) :: dcc
    character(len=36) :: name
    select case (dcc)
    case (1);  name = 'Actual/Actual (ISDA)'
    case (2);  name = 'Actual/Actual (ICMA)'
    case (3);  name = 'Actual/Actual (AFB)'
    case (4);  name = 'Actual/365L'
    case (5);  name = '30/360 Bond Basis'
    case (6);  name = '30E/360'
    case (7);  name = '30E/360 (ISDA)'
    case (8);  name = '30/360 German'
    case (9);  name = '30U/360 US'
    case (10); name = 'Actual/365 Fixed'
    case (11); name = 'Actual (No Leap)/365'
    case (12); name = 'Actual/360'
    case (13); name = '30/365'
    case (14); name = 'Actual/365 Canadian Bond'
    case (15); name = 'Actual/364'
    case (16); name = 'Business/252 Brazil'
    case default; name = 'Unknown'
    end select
  end function day_count_name

  recursive pure logical function leap_day_inside(first_date, second_date) result(ans)
    type(date_type), intent(in) :: first_date, second_date
    integer :: year, first_serial, second_serial, leap_serial
    ans = .false.
    if (.not. valid_date(first_date) .or. .not. valid_date(second_date)) return
    first_serial = date_to_serial(first_date)
    second_serial = date_to_serial(second_date)
    if (second_serial < first_serial) then
      ans = leap_day_inside(second_date, first_date)
      return
    end if
    do year = first_date%year, second_date%year
      if (leap_year(year)) then
        leap_serial = date_to_serial(date_type(year, 2, 29))
        if (first_serial <= leap_serial .and. leap_serial < second_serial) then
          ans = .true.
          return
        end if
      end if
    end do
  end function leap_day_inside

  function year_fraction(first_date, second_date, dcc, coupon_frequency, &
                         maturity, end_of_month, next_coupon_year, &
                         anniversary_dates) result(value)
    type(date_type), intent(in) :: first_date, second_date
    integer, intent(in) :: dcc
    integer, intent(in), optional :: coupon_frequency
    type(date_type), intent(in), optional :: maturity
    logical, intent(in), optional :: end_of_month
    integer, intent(in), optional :: next_coupon_year
    type(date_type), intent(in), optional :: anniversary_dates(:)
    real(dp) :: value
    type(daycount_result) :: result
    result = day_count_fraction(first_date, second_date, dcc, coupon_frequency, &
      maturity, end_of_month, next_coupon_year, anniversary_dates)
    value = result%fraction
  end function year_fraction

  function day_count_fraction(first_date, second_date, dcc, coupon_frequency, &
                              maturity, end_of_month, next_coupon_year, &
                              anniversary_dates) result(result)
    type(date_type), intent(in) :: first_date, second_date
    integer, intent(in) :: dcc
    integer, intent(in), optional :: coupon_frequency
    type(date_type), intent(in), optional :: maturity
    logical, intent(in), optional :: end_of_month
    integer, intent(in), optional :: next_coupon_year
    type(date_type), intent(in), optional :: anniversary_dates(:)
    type(daycount_result) :: result

    type(date_type) :: a, b, mat
    integer :: sign_value, cpy, next_year
    logical :: eom

    result = daycount_result()
    if (.not. valid_date(first_date) .or. .not. valid_date(second_date)) then
      result%status = 1
      return
    end if
    if (compare_dates(first_date, second_date) == 0) return
    if (compare_dates(first_date, second_date) < 0) then
      a = first_date
      b = second_date
      sign_value = 1
    else
      a = second_date
      b = first_date
      sign_value = -1
    end if
    cpy = 2
    if (present(coupon_frequency)) cpy = coupon_frequency
    if (cpy <= 0) cpy = 1
    eom = .false.
    if (present(end_of_month)) eom = end_of_month
    mat = b
    if (present(maturity)) mat = maturity
    next_year = b%year
    if (present(next_coupon_year)) next_year = next_coupon_year

    select case (dcc)
    case (dcc_act_act_isda)
      call actual_actual_isda(a, b, result)
    case (dcc_act_act_icma)
      if (present(anniversary_dates)) then
        call actual_actual_icma(a, b, cpy, eom, result, anniversary_dates)
      else
        call actual_actual_icma(a, b, cpy, eom, result)
      end if
    case (dcc_act_act_afb)
      call actual_actual_afb(a, b, result)
    case (dcc_act_365l)
      result%days_accrued = day_diff(a, b)
      if (cpy == 1) then
        if (leap_day_inside(a, b)) then
          result%fraction = real(result%days_accrued, dp) / 366.0_dp
        else
          result%fraction = real(result%days_accrued, dp) / 365.0_dp
        end if
      else if (leap_year(next_year)) then
        result%fraction = real(result%days_accrued, dp) / 366.0_dp
      else
        result%fraction = real(result%days_accrued, dp) / 365.0_dp
      end if
    case (dcc_30_360_bond)
      call thirty_360(a, b, dcc, eom, mat, result)
    case (dcc_30e_360)
      call thirty_360(a, b, dcc, eom, mat, result)
    case (dcc_30e_360_isda)
      call thirty_360(a, b, dcc, eom, mat, result)
    case (dcc_30_360_german)
      call thirty_360(a, b, dcc, eom, mat, result)
    case (dcc_30u_360)
      call thirty_360(a, b, dcc, eom, mat, result)
    case (dcc_act_365)
      result%days_accrued = day_diff(a, b)
      result%fraction = real(result%days_accrued, dp) / 365.0_dp
    case (dcc_act_nl_365)
      result%days_accrued = day_diff(a, b)
      ! This preserves the package behavior: at most one leap day is removed.
      if (leap_day_inside(a, b)) result%days_accrued = result%days_accrued - 1
      result%fraction = real(result%days_accrued, dp) / 365.0_dp
    case (dcc_act_360)
      result%days_accrued = day_diff(a, b)
      result%fraction = real(result%days_accrued, dp) / 360.0_dp
    case (dcc_30_365)
      call thirty_360(a, b, dcc, eom, mat, result)
    case (dcc_act_365_canadian)
      if (present(anniversary_dates)) then
        call actual_365_canadian(a, b, cpy, eom, result, anniversary_dates)
      else
        call actual_365_canadian(a, b, cpy, eom, result)
      end if
    case (dcc_act_364)
      result%days_accrued = day_diff(a, b)
      result%fraction = real(result%days_accrued, dp) / 364.0_dp
    case (dcc_bus_252)
      result%days_accrued = day_diff(a, b) - &
        count_nonbusiness_brazil(date_to_serial(a), date_to_serial(b))
      result%fraction = real(result%days_accrued, dp) / 252.0_dp
    case default
      result%status = 2
    end select
    result%days_accrued = sign_value * result%days_accrued
    result%fraction = real(sign_value, dp) * result%fraction
  end function day_count_fraction

  subroutine actual_actual_isda(a, b, result)
    type(date_type), intent(in) :: a, b
    type(daycount_result), intent(inout) :: result
    type(date_type) :: left, right
    integer :: year, n
    real(dp) :: fraction
    n = day_diff(a, b)
    fraction = 0.0_dp
    left = a
    do year = a%year, b%year
      if (year < b%year) then
        right = date_type(year + 1, 1, 1)
      else
        right = b
      end if
      if (compare_dates(right, left) > 0) &
        fraction = fraction + real(day_diff(left, right), dp) / &
                   real(days_in_year(year), dp)
      left = right
    end do
    result%days_accrued = n
    result%fraction = fraction
  end subroutine actual_actual_isda

  subroutine actual_actual_afb(a, b, result)
    type(date_type), intent(in) :: a, b
    type(daycount_result), intent(inout) :: result
    type(date_type) :: cursor, previous_year
    integer :: whole_years, remainder_days, denominator
    cursor = b
    whole_years = 0
    do
      previous_year = date_type(cursor%year - 1, cursor%month, cursor%day)
      if (.not. valid_date(previous_year)) previous_year = date_type(cursor%year - 1, 2, 28)
      if (compare_dates(previous_year, a) < 0) exit
      whole_years = whole_years + 1
      cursor = previous_year
    end do
    remainder_days = day_diff(a, cursor)
    denominator = 365
    if (leap_day_inside(a, cursor)) denominator = 366
    result%days_accrued = day_diff(a, b)
    result%fraction = real(whole_years, dp) + &
      real(remainder_days, dp) / real(denominator, dp)
  end subroutine actual_actual_afb

  subroutine actual_actual_icma(a, b, cpy, eom, result, supplied_dates)
    type(date_type), intent(in) :: a, b
    integer, intent(in) :: cpy
    logical, intent(in) :: eom
    type(daycount_result), intent(inout) :: result
    type(date_type), intent(in), optional :: supplied_dates(:)
    type(date_type), allocatable :: dates(:)
    type(date_type) :: anchor
    integer :: months, i, first_index, last_index, n, k0, k1
    real(dp) :: total

    months = 12 / cpy
    if (present(supplied_dates)) then
      allocate(dates(size(supplied_dates)))
      dates = supplied_dates
    else
      anchor = date_type(a%year, a%month, 15)
      allocate(dates(129))
      n = 0
      do k0 = -64, 64
        n = n + 1
        dates(n) = add_months(anchor, k0 * months, eom, a%day)
      end do
      dates = dates(1:n)
    end if
    first_index = find_previous_date(a, dates)
    if (first_index == 0) first_index = 1
    last_index = find_next_date(b, dates)
    if (last_index == 0) last_index = size(dates)
    total = 0.0_dp
    do i = first_index, last_index - 1
      k0 = max(date_to_serial(a), date_to_serial(dates(i)))
      k1 = min(date_to_serial(b), date_to_serial(dates(i + 1)))
      if (k1 > k0) then
        total = total + real(k1 - k0, dp) / &
          (real(day_diff(dates(i), dates(i + 1)), dp) * real(cpy, dp))
      end if
    end do
    result%days_accrued = day_diff(a, b)
    result%fraction = total
  end subroutine actual_actual_icma

  subroutine actual_365_canadian(a, b, cpy, eom, result, supplied_dates)
    type(date_type), intent(in) :: a, b
    integer, intent(in) :: cpy
    logical, intent(in) :: eom
    type(daycount_result), intent(inout) :: result
    type(date_type), intent(in), optional :: supplied_dates(:)
    type(date_type), allocatable :: dates(:)
    type(date_type) :: anchor
    integer :: months, k, ip, inext, num, den
    real(dp) :: coordinate_a, coordinate_b

    months = 12 / cpy
    if (present(supplied_dates)) then
      allocate(dates(size(supplied_dates)))
      dates = supplied_dates
    else
      anchor = date_type(a%year, a%month, 15)
      allocate(dates(129))
      do k = -64, 64
        dates(k + 65) = add_months(anchor, k * months, eom, a%day)
      end do
    end if
    ip = find_previous_date(a, dates)
    inext = find_next_date(a, dates)
    if (ip == 0 .or. inext == 0) then
      result%status = 3
      return
    end if
    num = day_diff(a, dates(inext))
    den = day_diff(dates(ip), dates(inext))
    if (num == den) then
      coordinate_a = real(ip, dp)
    else if (num < 365 / cpy) then
      coordinate_a = real(inext, dp) - real(cpy * num, dp) / 365.0_dp
    else
      coordinate_a = real(ip, dp) + real(cpy * day_diff(dates(ip), a), dp) / 365.0_dp
    end if
    ip = find_previous_date(b, dates)
    inext = find_next_date(b, dates)
    if (ip == 0) ip = 1
    if (inext == 0) inext = size(dates)
    num = day_diff(dates(ip), b)
    if (num == 0) then
      coordinate_b = real(ip, dp)
    else if (num < 365 / cpy) then
      coordinate_b = real(ip, dp) + real(cpy * num, dp) / 365.0_dp
    else
      coordinate_b = real(inext, dp) - &
        real(cpy * day_diff(b, dates(inext)), dp) / 365.0_dp
    end if
    result%days_accrued = day_diff(a, b)
    result%fraction = (coordinate_b - coordinate_a) / real(cpy, dp)
  end subroutine actual_365_canadian

  subroutine thirty_360(a, b, dcc, eom, maturity, result)
    type(date_type), intent(in) :: a, b, maturity
    integer, intent(in) :: dcc
    logical, intent(in) :: eom
    type(daycount_result), intent(inout) :: result
    integer :: d1, d2
    logical :: first_feb_end, second_feb_end, second_is_maturity

    d1 = a%day
    d2 = b%day
    first_feb_end = a%month == 2 .and. is_last_day_of_month(a)
    second_feb_end = b%month == 2 .and. is_last_day_of_month(b)
    second_is_maturity = compare_dates(b, maturity) == 0
    select case (dcc)
    case (dcc_30_360_bond)
      d1 = min(d1, 30)
      if (b%day == 31 .and. (a%day == 30 .or. a%day == 31)) d2 = 30
    case (dcc_30e_360)
      d1 = min(d1, 30)
      d2 = min(d2, 30)
    case (dcc_30e_360_isda)
      if (is_last_day_of_month(a)) d1 = 30
      if (b%day == 31) then
        d2 = 30
      else if (.not. second_is_maturity .and. second_feb_end) then
        d2 = 30
      end if
    case (dcc_30_360_german)
      if (is_last_day_of_month(a)) d1 = 30
      if (is_last_day_of_month(b)) d2 = 30
    case (dcc_30u_360, dcc_30_365)
      if ((eom .and. first_feb_end) .or. a%day == 31) d1 = 30
      if (eom .and. first_feb_end .and. second_feb_end) then
        d2 = 30
      else if (b%day == 31 .and. (a%day == 30 .or. a%day == 31)) then
        d2 = 30
      end if
    end select
    result%days_accrued = 360 * (b%year - a%year) + &
      30 * (b%month - a%month) + d2 - d1
    if (dcc == dcc_30_365) then
      result%fraction = real(result%days_accrued, dp) / 365.0_dp
    else
      result%fraction = real(result%days_accrued, dp) / 360.0_dp
    end if
  end subroutine thirty_360

end module bondvaluation_daycount
