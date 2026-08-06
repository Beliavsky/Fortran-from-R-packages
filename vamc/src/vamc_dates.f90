module vamc_dates
  use vamc_kinds, only: dp, i8
  use vamc_status, only: status_type, vamc_invalid_argument, vamc_date_error
  implicit none
  private

  integer, parameter, public :: dcc_thirty360 = 1
  integer, parameter, public :: dcc_act360 = 2
  integer, parameter, public :: dcc_act365 = 3
  integer, parameter, public :: dcc_actact = 4
  integer, parameter, public :: bdc_actual = 1
  integer, parameter, public :: bdc_following = 2
  integer, parameter, public :: bdc_preceding = 3
  integer, parameter, public :: bdc_modified_preceding = 4
  integer, parameter, public :: bdc_modified_following = 5
  integer, parameter, public :: calendar_general = 1
  integer, parameter, public :: calendar_ny = 2

  type, public :: date_type
    integer :: year = 1970
    integer :: month = 1
    integer :: day = 1
  contains
    procedure :: serial => date_serial_method
    procedure :: to_string => date_to_string
  end type date_type

  interface operator(==)
    module procedure dates_equal
  end interface
  interface operator(/=)
    module procedure dates_not_equal
  end interface
  interface operator(<)
    module procedure date_less
  end interface
  interface operator(<=)
    module procedure date_less_equal
  end interface
  interface operator(>)
    module procedure date_greater
  end interface
  interface operator(>=)
    module procedure date_greater_equal
  end interface
  public :: operator(==), operator(/=), operator(<), operator(<=), operator(>), operator(>=)
  public :: make_date, parse_date, date_from_serial, days_in_month, is_leap_year
  public :: add_days, add_months, add_years, months_between, days_between
  public :: weekday_number, is_weekend, is_holiday_ny, is_business_day, roll_date
  public :: gen_schedule, frac_year, dcc_from_string, bdc_from_string, calendar_from_string
contains
  pure logical function is_leap_year(year)
    integer, intent(in) :: year
    is_leap_year = (mod(year, 4) == 0 .and. mod(year, 100) /= 0) .or. mod(year, 400) == 0
  end function is_leap_year

  pure integer function days_in_month(year, month)
    integer, intent(in) :: year, month
    integer, parameter :: mdays(12) = [31,28,31,30,31,30,31,31,30,31,30,31]
    if (month < 1 .or. month > 12) then
      days_in_month = 0
    else
      days_in_month = mdays(month)
      if (month == 2 .and. is_leap_year(year)) days_in_month = 29
    end if
  end function days_in_month

  pure function make_date(year, month, day) result(date)
    integer, intent(in) :: year, month, day
    type(date_type) :: date
    date%year = year
    date%month = month
    date%day = min(max(day, 1), max(days_in_month(year, month), 1))
  end function make_date

  subroutine parse_date(text, date, status)
    character(len=*), intent(in) :: text
    type(date_type), intent(out) :: date
    type(status_type), intent(inout), optional :: status
    integer :: y, m, d, ios
    if (present(status)) call status%clear()
    read(text, '(i4,1x,i2,1x,i2)', iostat=ios) y, m, d
    if (ios /= 0 .or. m < 1 .or. m > 12 .or. d < 1 .or. d > days_in_month(y, m)) then
      date = date_type()
      if (present(status)) call status%set(vamc_date_error, 'Invalid date; expected YYYY-MM-DD.')
      return
    end if
    date = make_date(y, m, d)
  end subroutine parse_date

  pure integer(i8) function days_from_civil(year, month, day)
    integer, intent(in) :: year, month, day
    integer(i8) :: y, m, era, yoe, doy, doe
    y = int(year, i8)
    m = int(month, i8)
    if (m <= 2_i8) y = y - 1_i8
    if (y >= 0_i8) then
      era = y / 400_i8
    else
      era = (y - 399_i8) / 400_i8
    end if
    yoe = y - era * 400_i8
    if (m > 2_i8) then
      doy = (153_i8 * (m - 3_i8) + 2_i8) / 5_i8 + int(day, i8) - 1_i8
    else
      doy = (153_i8 * (m + 9_i8) + 2_i8) / 5_i8 + int(day, i8) - 1_i8
    end if
    doe = yoe * 365_i8 + yoe / 4_i8 - yoe / 100_i8 + doy
    days_from_civil = era * 146097_i8 + doe - 719468_i8
  end function days_from_civil

  pure integer(i8) function date_serial_method(self)
    class(date_type), intent(in) :: self
    date_serial_method = days_from_civil(self%year, self%month, self%day)
  end function date_serial_method

  pure function date_from_serial(serial) result(date)
    integer(i8), intent(in) :: serial
    type(date_type) :: date
    integer(i8) :: z, era, doe, yoe, y, doy, mp, d, m
    z = serial + 719468_i8
    if (z >= 0_i8) then
      era = z / 146097_i8
    else
      era = (z - 146096_i8) / 146097_i8
    end if
    doe = z - era * 146097_i8
    yoe = (doe - doe / 1460_i8 + doe / 36524_i8 - doe / 146096_i8) / 365_i8
    y = yoe + era * 400_i8
    doy = doe - (365_i8 * yoe + yoe / 4_i8 - yoe / 100_i8)
    mp = (5_i8 * doy + 2_i8) / 153_i8
    d = doy - (153_i8 * mp + 2_i8) / 5_i8 + 1_i8
    if (mp < 10_i8) then
      m = mp + 3_i8
    else
      m = mp - 9_i8
    end if
    if (m <= 2_i8) y = y + 1_i8
    date = make_date(int(y), int(m), int(d))
  end function date_from_serial

  function date_to_string(self) result(text)
    class(date_type), intent(in) :: self
    character(len=10) :: text
    write(text, '(i4.4,"-",i2.2,"-",i2.2)') self%year, self%month, self%day
  end function date_to_string

  pure logical function dates_equal(a, b)
    type(date_type), intent(in) :: a, b
    dates_equal = a%year == b%year .and. a%month == b%month .and. a%day == b%day
  end function dates_equal
  pure logical function dates_not_equal(a, b)
    type(date_type), intent(in) :: a, b
    dates_not_equal = .not. dates_equal(a, b)
  end function dates_not_equal
  pure logical function date_less(a, b)
    type(date_type), intent(in) :: a, b
    date_less = a%serial() < b%serial()
  end function date_less
  pure logical function date_less_equal(a, b)
    type(date_type), intent(in) :: a, b
    date_less_equal = a%serial() <= b%serial()
  end function date_less_equal
  pure logical function date_greater(a, b)
    type(date_type), intent(in) :: a, b
    date_greater = a%serial() > b%serial()
  end function date_greater
  pure logical function date_greater_equal(a, b)
    type(date_type), intent(in) :: a, b
    date_greater_equal = a%serial() >= b%serial()
  end function date_greater_equal

  pure function add_days(date, number) result(out)
    type(date_type), intent(in) :: date
    integer, intent(in) :: number
    type(date_type) :: out
    out = date_from_serial(date%serial() + int(number, i8))
  end function add_days

  pure function add_months(date, number) result(out)
    type(date_type), intent(in) :: date
    integer, intent(in) :: number
    type(date_type) :: out
    integer :: total, y, m
    total = date%year * 12 + date%month - 1 + number
    y = floor(real(total, dp) / 12.0_dp)
    m = total - 12 * y + 1
    out = make_date(y, m, min(date%day, days_in_month(y, m)))
  end function add_months

  pure function add_years(date, number) result(out)
    type(date_type), intent(in) :: date
    integer, intent(in) :: number
    type(date_type) :: out
    out = make_date(date%year + number, date%month, min(date%day, days_in_month(date%year + number, date%month)))
  end function add_years

  pure integer function days_between(a, b)
    type(date_type), intent(in) :: a, b
    days_between = int(b%serial() - a%serial())
  end function days_between

  pure integer function months_between(a, b)
    type(date_type), intent(in) :: a, b
    months_between = 12 * (b%year - a%year) + b%month - a%month
  end function months_between

  pure integer function weekday_number(date)
    type(date_type), intent(in) :: date
    ! Monday=1, ..., Sunday=7. 1970-01-01 was Thursday.
    weekday_number = modulo(int(date%serial()) + 3, 7) + 1
  end function weekday_number

  pure logical function is_weekend(date)
    type(date_type), intent(in) :: date
    integer :: w
    w = weekday_number(date)
    is_weekend = w >= 6
  end function is_weekend

  pure logical function is_holiday_ny(date)
    type(date_type), intent(in) :: date
    integer :: d, m, w
    d = date%day
    m = date%month
    w = weekday_number(date)
    is_holiday_ny = &
      ((d >= 15 .and. d <= 21) .and. w == 1 .and. m == 1) .or. &
      ((d >= 15 .and. d <= 21) .and. w == 1 .and. m == 2) .or. &
      (d >= 25 .and. w == 1 .and. m == 5) .or. &
      ((d == 4 .and. m == 7) .or. (d == 5 .and. w == 1 .and. m == 7) .or. (d == 3 .and. w == 5 .and. m == 7)) .or. &
      (d <= 7 .and. w == 1 .and. m == 9) .or. &
      ((d >= 8 .and. d <= 14) .and. w == 1 .and. m == 10) .or. &
      ((d == 11 .and. m == 11) .or. (d == 12 .and. w == 1 .and. m == 11) .or. (d == 10 .and. w == 5 .and. m == 11)) .or. &
      ((d >= 22 .and. d <= 28) .and. w == 4 .and. m == 11) .or. &
      ((d == 25 .and. m == 12) .or. (d == 26 .and. w == 1 .and. m == 12) .or. (d == 24 .and. w == 5 .and. m == 12)) .or. &
      ((d == 1 .and. m == 1) .or. (d == 2 .and. w == 1 .and. m == 1) .or. (d == 31 .and. w == 5 .and. m == 12))
  end function is_holiday_ny

  logical function is_business_day(date, calendar, holidays, source_compatible_holidays)
    type(date_type), intent(in) :: date
    integer, intent(in) :: calendar
    type(date_type), intent(in), optional :: holidays(:)
    logical, intent(in), optional :: source_compatible_holidays
    logical :: source_mode, listed
    integer :: i
    source_mode = .true.
    if (present(source_compatible_holidays)) source_mode = source_compatible_holidays
    if (present(holidays)) then
      if (size(holidays) > 0) then
        if (date >= min_date(holidays) .and. date <= max_date(holidays)) then
          listed = .false.
          do i = 1, size(holidays)
            if (date == holidays(i)) then
              listed = .true.
              exit
            end if
          end do
          if (source_mode) then
            is_business_day = listed
          else
            is_business_day = .not. listed .and. .not. is_weekend(date)
          end if
          return
        end if
      end if
    end if
    select case (calendar)
    case (calendar_ny)
      is_business_day = .not. is_weekend(date) .and. .not. is_holiday_ny(date)
    case default
      is_business_day = .not. is_weekend(date)
    end select
  end function is_business_day

  pure function min_date(dates) result(out)
    type(date_type), intent(in) :: dates(:)
    type(date_type) :: out
    integer :: i
    out = dates(1)
    do i = 2, size(dates)
      if (dates(i) < out) out = dates(i)
    end do
  end function min_date

  pure function max_date(dates) result(out)
    type(date_type), intent(in) :: dates(:)
    type(date_type) :: out
    integer :: i
    out = dates(1)
    do i = 2, size(dates)
      if (dates(i) > out) out = dates(i)
    end do
  end function max_date

  recursive function roll_date(date, convention, calendar, holidays, source_compatible_holidays) result(out)
    type(date_type), intent(in) :: date
    integer, intent(in) :: convention, calendar
    type(date_type), intent(in), optional :: holidays(:)
    logical, intent(in), optional :: source_compatible_holidays
    type(date_type) :: out, temp
    if (convention == bdc_actual) then
      out = date
      return
    end if
    select case (convention)
    case (bdc_following)
      temp = date
      do while (.not. is_business_day(temp, calendar, holidays, source_compatible_holidays))
        temp = add_days(temp, 1)
      end do
      out = temp
    case (bdc_preceding)
      temp = date
      do while (.not. is_business_day(temp, calendar, holidays, source_compatible_holidays))
        temp = add_days(temp, -1)
      end do
      out = temp
    case (bdc_modified_preceding)
      temp = roll_date(date, bdc_preceding, calendar, holidays, source_compatible_holidays)
      if (temp%month /= date%month) temp = roll_date(date, bdc_following, calendar, holidays, source_compatible_holidays)
      out = temp
    case (bdc_modified_following)
      temp = roll_date(date, bdc_following, calendar, holidays, source_compatible_holidays)
      if (temp%month /= date%month) temp = roll_date(date, bdc_preceding, calendar, holidays, source_compatible_holidays)
      out = temp
    case default
      out = date
    end select
  end function roll_date

  subroutine gen_schedule(settle_date, frequency_months, tenor_years, calendar, convention, schedule, holidays, &
                          source_compatible_holidays, status)
    type(date_type), intent(in) :: settle_date
    integer, intent(in) :: frequency_months, tenor_years, calendar, convention
    type(date_type), allocatable, intent(out) :: schedule(:)
    type(date_type), intent(in), optional :: holidays(:)
    logical, intent(in), optional :: source_compatible_holidays
    type(status_type), intent(inout), optional :: status
    integer :: count, i
    if (present(status)) call status%clear()
    if (frequency_months <= 0 .or. tenor_years < 0 .or. mod(tenor_years * 12, frequency_months) /= 0) then
      allocate(schedule(0))
      if (present(status)) call status%set(vamc_invalid_argument, 'Frequency must divide the tenor in months.')
      return
    end if
    count = tenor_years * 12 / frequency_months
    allocate(schedule(count + 1))
    do i = 0, count
      schedule(i + 1) = roll_date(add_months(settle_date, i * frequency_months), convention, calendar, holidays, &
                                  source_compatible_holidays)
    end do
  end subroutine gen_schedule

  real(dp) function frac_year(date_a, date_b, convention)
    type(date_type), intent(in) :: date_a, date_b
    integer, intent(in) :: convention
    type(date_type) :: a, b, start_y, end_y, lo, hi
    integer :: da, db, leap_days, y, overlap
    if (date_b < date_a) then
      a = date_b
      b = date_a
    else
      a = date_a
      b = date_b
    end if
    select case (convention)
    case (dcc_thirty360)
      da = min(a%day, 30)
      db = b%day
      if (db == 31 .and. da == 30) db = 30
      frac_year = real(b%year - a%year, dp) + real(b%month - a%month, dp) / 12.0_dp + real(db - da, dp) / 360.0_dp
    case (dcc_act360)
      frac_year = real(days_between(a, b), dp) / 360.0_dp
    case (dcc_act365)
      frac_year = real(days_between(a, b), dp) / 365.0_dp
    case (dcc_actact)
      leap_days = 0
      do y = a%year, b%year
        if (.not. is_leap_year(y)) cycle
        start_y = make_date(y, 1, 1)
        end_y = make_date(y + 1, 1, 1)
        if (a > start_y) then
          lo = a
        else
          lo = start_y
        end if
        if (b < end_y) then
          hi = b
        else
          hi = end_y
        end if
        overlap = max(0, days_between(lo, hi))
        leap_days = leap_days + overlap
      end do
      frac_year = real(leap_days, dp) / 366.0_dp + real(days_between(a, b) - leap_days, dp) / 365.0_dp
    case default
      frac_year = 0.0_dp
    end select
  end function frac_year

  integer function dcc_from_string(name)
    character(len=*), intent(in) :: name
    select case (trim(name))
    case ('Thirty360')
      dcc_from_string = dcc_thirty360
    case ('ACT360')
      dcc_from_string = dcc_act360
    case ('ACT365')
      dcc_from_string = dcc_act365
    case ('ACTACT')
      dcc_from_string = dcc_actact
    case default
      dcc_from_string = 0
    end select
  end function dcc_from_string

  integer function bdc_from_string(name)
    character(len=*), intent(in) :: name
    select case (trim(name))
    case ('Actual')
      bdc_from_string = bdc_actual
    case ('Following')
      bdc_from_string = bdc_following
    case ('Preceding')
      bdc_from_string = bdc_preceding
    case ('Modified_Prec')
      bdc_from_string = bdc_modified_preceding
    case ('Modified_Foll')
      bdc_from_string = bdc_modified_following
    case default
      bdc_from_string = 0
    end select
  end function bdc_from_string

  integer function calendar_from_string(name)
    character(len=*), intent(in) :: name
    select case (trim(name))
    case ('General')
      calendar_from_string = calendar_general
    case ('NY')
      calendar_from_string = calendar_ny
    case default
      calendar_from_string = 0
    end select
  end function calendar_from_string
end module vamc_dates
