! SPDX-License-Identifier: GPL-3.0-only AND LicenseRef-ISDA-CDS-Standard-Model
module creditr_dates
  use creditr_kinds, only: dp, creditr_ok, creditr_invalid_input
  implicit none
  private

  type, public :: date_t
    integer :: year = 1970
    integer :: month = 1
    integer :: day = 1
  contains
    procedure :: serial => date_serial_type
  end type date_t

  type, public :: cds_dates_t
    type(date_t) :: trade_date
    type(date_t) :: stepin_date
    type(date_t) :: value_date
    type(date_t) :: start_date
    type(date_t) :: first_coupon_date
    type(date_t) :: penultimate_coupon_date
    type(date_t) :: end_date
    type(date_t) :: backstop_date
    type(date_t) :: base_date
  end type cds_dates_t

  type, public :: conventions_t
    character(len=1) :: bad_day = 'M'
    character(len=12) :: mm_dcc = 'ACT/360'
    character(len=12) :: fixed_dcc = '30/360'
    character(len=12) :: float_dcc = 'ACT/360'
    integer :: fixed_frequency_months = 6
    integer :: float_frequency_months = 3
    character(len=8) :: mm_calendar = 'NONE'
    character(len=8) :: swap_calendar = 'NONE'
  end type conventions_t

  public :: make_date, date_from_serial, date_to_string, parse_date
  public :: add_days, add_months, compare_dates, is_weekend
  public :: following, modified_following, add_business_days
  public :: year_fraction, days_between, add_dates, add_conventions
  public :: is_valid_date

contains

  pure type(date_t) function make_date(year, month, day) result(dt)
    integer, intent(in) :: year, month, day
    dt%year = year
    dt%month = month
    dt%day = day
  end function make_date

  pure logical function is_leap_year(year) result(answer)
    integer, intent(in) :: year
    answer = (mod(year, 4) == 0 .and. mod(year, 100) /= 0) .or. mod(year, 400) == 0
  end function is_leap_year

  pure integer function days_in_month(year, month) result(n)
    integer, intent(in) :: year, month
    integer, parameter :: month_days(12) = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
    if (month < 1 .or. month > 12) then
      n = 0
    else
      n = month_days(month)
      if (month == 2 .and. is_leap_year(year)) n = 29
    end if
  end function days_in_month

  pure logical function is_valid_date(dt) result(answer)
    type(date_t), intent(in) :: dt
    answer = dt%month >= 1 .and. dt%month <= 12
    if (answer) answer = dt%day >= 1 .and. dt%day <= days_in_month(dt%year, dt%month)
  end function is_valid_date

  pure integer function date_to_serial_values(year, month, day) result(serial)
    integer, intent(in) :: year, month, day
    integer :: a, y, m
    a = (14 - month) / 12
    y = year + 4800 - a
    m = month + 12 * a - 3
    serial = day + (153 * m + 2) / 5 + 365 * y + y / 4 - y / 100 + y / 400 - 32045
  end function date_to_serial_values

  pure integer function date_serial_type(self) result(serial)
    class(date_t), intent(in) :: self
    serial = date_to_serial_values(self%year, self%month, self%day)
  end function date_serial_type

  pure type(date_t) function date_from_serial(serial) result(dt)
    integer, intent(in) :: serial
    integer :: a, b, c, d, e, m
    a = serial + 32044
    b = (4 * a + 3) / 146097
    c = a - (146097 * b) / 4
    d = (4 * c + 3) / 1461
    e = c - (1461 * d) / 4
    m = (5 * e + 2) / 153
    dt%day = e - (153 * m + 2) / 5 + 1
    dt%month = m + 3 - 12 * (m / 10)
    dt%year = 100 * b + d - 4800 + m / 10
  end function date_from_serial

  pure type(date_t) function add_days(dt, n) result(out)
    type(date_t), intent(in) :: dt
    integer, intent(in) :: n
    out = date_from_serial(dt%serial() + n)
  end function add_days

  pure type(date_t) function add_months(dt, n) result(out)
    type(date_t), intent(in) :: dt
    integer, intent(in) :: n
    integer :: total, y, m
    total = dt%year * 12 + dt%month - 1 + n
    y = floor_div(total, 12)
    m = total - 12 * y + 1
    out = make_date(y, m, min(dt%day, days_in_month(y, m)))
  end function add_months

  pure integer function floor_div(a, b) result(q)
    integer, intent(in) :: a, b
    q = a / b
    if (a < 0 .and. mod(a, b) /= 0) q = q - 1
  end function floor_div

  pure integer function compare_dates(a, b) result(cmp)
    type(date_t), intent(in) :: a, b
    if (a%serial() < b%serial()) then
      cmp = -1
    else if (a%serial() > b%serial()) then
      cmp = 1
    else
      cmp = 0
    end if
  end function compare_dates

  pure integer function days_between(a, b) result(n)
    type(date_t), intent(in) :: a, b
    n = b%serial() - a%serial()
  end function days_between

  pure logical function is_weekend(dt) result(answer)
    type(date_t), intent(in) :: dt
    integer :: w
    w = modulo(dt%serial(), 7)
    answer = w == 5 .or. w == 6
  end function is_weekend

  pure logical function is_usd_holiday(dt) result(answer)
    type(date_t), intent(in) :: dt
    integer :: s, i
    integer, parameter :: holidays(54) = [ &
      2457024, 2457017, 2456989, 2456973, 2456944, 2456902, 2456843, 2456804, 2456766, &
      2456706, 2456678, 2456659, 2456652, 2456625, 2456608, 2456580, 2456538, 2456478, &
      2456440, 2456381, 2456342, 2456314, 2456294, 2456287, 2456254, 2456244, 2456231, &
      2456209, 2456174, 2456113, 2456076, 2455978, 2455943, 2455929, 2455922, 2455890, &
      2455877, 2455845, 2455810, 2455747, 2455712, 2455674, 2455614, 2455579, 2455555, &
      2455526, 2455512, 2455481, 2455446, 2455383, 2455348, 2455243, 2455215, 2455198]
    s = dt%serial()
    answer = .false.
    do i = 1, size(holidays)
      if (s == holidays(i)) then
        answer = .true.
        return
      end if
    end do
  end function is_usd_holiday

  pure logical function is_jpy_holiday(dt) result(answer)
    type(date_t), intent(in) :: dt
    integer :: s, i
    integer, parameter :: holidays(17) = [2454911,2455096,2455097,2455098,2455278,2455460,2455642, &
      2456007,2456372,2457287,2457288,2457289,2457469,2457833,2458929,2459114,2459115]
    s = dt%serial()
    answer = .false.
    do i = 1, size(holidays)
      if (s == holidays(i)) then
        answer = .true.
        return
      end if
    end do
  end function is_jpy_holiday

  pure logical function is_business_day(dt, currency, include_holidays) result(answer)
    type(date_t), intent(in) :: dt
    character(len=*), intent(in) :: currency
    logical, intent(in) :: include_holidays
    answer = .not. is_weekend(dt)
    if (.not. answer .or. .not. include_holidays) return
    select case (trim(currency))
    case ('USD')
      answer = .not. is_usd_holiday(dt)
    case ('JPY')
      answer = .not. is_jpy_holiday(dt)
    case default
      answer = .true.
    end select
  end function is_business_day

  pure type(date_t) function following(dt, currency, include_holidays) result(out)
    type(date_t), intent(in) :: dt
    character(len=*), intent(in), optional :: currency
    logical, intent(in), optional :: include_holidays
    character(len=3) :: ccy
    logical :: use_holidays
    ccy = '   '
    if (present(currency)) ccy = currency(1:min(3, len_trim(currency)))
    use_holidays = .false.
    if (present(include_holidays)) use_holidays = include_holidays
    out = dt
    do while (.not. is_business_day(out, ccy, use_holidays))
      out = add_days(out, 1)
    end do
  end function following

  pure type(date_t) function modified_following(dt, currency, include_holidays) result(out)
    type(date_t), intent(in) :: dt
    character(len=*), intent(in), optional :: currency
    logical, intent(in), optional :: include_holidays
    type(date_t) :: f
    character(len=3) :: ccy
    logical :: use_holidays
    ccy = '   '
    if (present(currency)) ccy = currency(1:min(3, len_trim(currency)))
    use_holidays = .false.
    if (present(include_holidays)) use_holidays = include_holidays
    f = following(dt, ccy, use_holidays)
    if (f%month == dt%month) then
      out = f
    else
      out = dt
      do
        out = add_days(out, -1)
        if (is_business_day(out, ccy, use_holidays)) exit
      end do
    end if
  end function modified_following

  pure type(date_t) function add_business_days(dt, n, currency, include_holidays) result(out)
    type(date_t), intent(in) :: dt
    integer, intent(in) :: n
    character(len=*), intent(in), optional :: currency
    logical, intent(in), optional :: include_holidays
    integer :: left, direction
    character(len=3) :: ccy
    logical :: use_holidays
    ccy = '   '
    if (present(currency)) ccy = currency(1:min(3, len_trim(currency)))
    use_holidays = .false.
    if (present(include_holidays)) use_holidays = include_holidays
    out = dt
    left = abs(n)
    direction = merge(1, -1, n >= 0)
    do while (left > 0)
      out = add_days(out, direction)
      if (is_business_day(out, ccy, use_holidays)) left = left - 1
    end do
  end function add_business_days

  pure real(kind=dp) function year_fraction(start_date, end_date, convention) result(yf)
    type(date_t), intent(in) :: start_date, end_date
    character(len=*), intent(in) :: convention
    integer :: d1, d2
    select case (trim(adjustl(convention)))
    case ('ACT/360')
      yf = real(days_between(start_date, end_date), dp) / 360.0_dp
    case ('ACT/365', 'ACT/365F')
      yf = real(days_between(start_date, end_date), dp) / 365.0_dp
    case ('30/360')
      d1 = min(start_date%day, 30)
      d2 = end_date%day
      if (d1 == 30) d2 = min(d2, 30)
      yf = real(360 * (end_date%year - start_date%year) + &
        30 * (end_date%month - start_date%month) + d2 - d1, dp) / 360.0_dp
    case default
      yf = real(days_between(start_date, end_date), dp) / 365.0_dp
    end select
  end function year_fraction

  pure subroutine add_conventions(currency, conventions, status)
    character(len=*), intent(in) :: currency
    type(conventions_t), intent(out) :: conventions
    integer, intent(out), optional :: status
    conventions = conventions_t()
    select case (trim(currency))
    case ('USD')
      conventions%fixed_dcc = '30/360'
      conventions%fixed_frequency_months = 6
      conventions%float_frequency_months = 3
    case ('EUR')
      conventions%fixed_dcc = '30/360'
      conventions%fixed_frequency_months = 12
      conventions%float_frequency_months = 6
    case ('JPY')
      conventions%mm_dcc = 'ACT/365'
      conventions%fixed_dcc = 'ACT/365'
      conventions%fixed_frequency_months = 6
      conventions%float_frequency_months = 6
    case default
      if (present(status)) status = creditr_invalid_input
      return
    end select
    if (present(status)) status = creditr_ok
  end subroutine add_conventions

  pure subroutine add_dates(trade_date, currency, dates, tenor_years, maturity, status)
    type(date_t), intent(in) :: trade_date
    character(len=*), intent(in) :: currency
    type(cds_dates_t), intent(out) :: dates
    integer, intent(in), optional :: tenor_years
    type(date_t), intent(in), optional :: maturity
    integer, intent(out), optional :: status
    type(date_t) :: roll, candidate, base
    integer :: quarter_month

    if (.not. is_valid_date(trade_date) .or. (present(tenor_years) .eqv. present(maturity))) then
      if (present(status)) status = creditr_invalid_input
      return
    end if

    dates%trade_date = trade_date
    dates%stepin_date = add_days(trade_date, 1)
    dates%value_date = add_business_days(trade_date, 3, '', .false.)

    base = add_business_days(trade_date, 1, currency, .true.)
    base = add_business_days(base, 1, currency, .true.)
    dates%base_date = base

    quarter_month = 3 * ((trade_date%month - 1) / 3 + 1)
    candidate = make_date(trade_date%year, quarter_month, 20)
    if (candidate%serial() > trade_date%serial()) candidate = add_months(candidate, -3)
    roll = candidate
    dates%start_date = following(roll)
    dates%first_coupon_date = following(add_months(dates%start_date, 3))

    if (present(tenor_years)) then
      dates%end_date = add_months(roll, 12 * tenor_years + 3)
    else
      dates%end_date = maturity
    end if
    dates%penultimate_coupon_date = following(add_months(dates%end_date, -3))
    dates%backstop_date = add_days(trade_date, -60)
    if (present(status)) status = creditr_ok
  end subroutine add_dates

  pure function date_to_string(dt) result(text)
    type(date_t), intent(in) :: dt
    character(len=10) :: text
    write(text, '(i4.4,"-",i2.2,"-",i2.2)') dt%year, dt%month, dt%day
  end function date_to_string

  subroutine parse_date(text, dt, status)
    character(len=*), intent(in) :: text
    type(date_t), intent(out) :: dt
    integer, intent(out), optional :: status
    integer :: ios
    read(text, '(i4,1x,i2,1x,i2)', iostat=ios) dt%year, dt%month, dt%day
    if (ios /= 0 .or. .not. is_valid_date(dt)) then
      if (present(status)) status = creditr_invalid_input
    else
      if (present(status)) status = creditr_ok
    end if
  end subroutine parse_date

end module creditr_dates
