! SPDX-License-Identifier: GPL-2.0-or-later
module jrvfinance_dates
  use jrvfinance_kinds, only: dp
  use jrvfinance_types, only: date_t, JRV_OK, JRV_INVALID_ARGUMENT
  implicit none
  private
  public :: date, parse_date, date_string, date_ordinal, ordinal_date
  public :: daycount_actual, daycount_30_360, year_fraction, edate
  public :: is_leap_year, last_day_of_month
  public :: operator(==), operator(/=), operator(<), operator(<=), operator(>), operator(>=)

  interface date
    module procedure make_date
  end interface
  interface operator(==)
    module procedure date_eq
  end interface
  interface operator(/=)
    module procedure date_ne
  end interface
  interface operator(<)
    module procedure date_lt
  end interface
  interface operator(<=)
    module procedure date_le
  end interface
  interface operator(>)
    module procedure date_gt
  end interface
  interface operator(>=)
    module procedure date_ge
  end interface
contains
  pure function make_date(year, month, day) result(d)
    integer, intent(in) :: year, month, day
    type(date_t) :: d
    d%year = year
    d%month = month
    d%day = day
  end function make_date

  function parse_date(text, status) result(d)
    character(len=*), intent(in) :: text
    integer, intent(out), optional :: status
    type(date_t) :: d
    integer :: ios
    character(len=10) :: s
    s = text
    read(s(1:4), *, iostat=ios) d%year
    if (ios == 0) read(s(6:7), *, iostat=ios) d%month
    if (ios == 0) read(s(9:10), *, iostat=ios) d%day
    if (ios /= 0 .or. .not. valid_date(d)) then
      d = date_t()
      if (present(status)) status = JRV_INVALID_ARGUMENT
    else
      if (present(status)) status = JRV_OK
    end if
  end function parse_date

  pure function date_string(d) result(text)
    type(date_t), intent(in) :: d
    character(len=10) :: text
    write(text, '(i4.4,"-",i2.2,"-",i2.2)') d%year, d%month, d%day
  end function date_string

  pure logical function is_leap_year(year) result(value)
    integer, intent(in) :: year
    value = modulo(year, 4) == 0 .and. (modulo(year, 100) /= 0 .or. modulo(year, 400) == 0)
  end function is_leap_year

  pure integer function last_day_of_month(month, year) result(last)
    integer, intent(in) :: month, year
    integer, parameter :: mdays(12) = [31,28,31,30,31,30,31,31,30,31,30,31]
    if (month < 1 .or. month > 12) then
      last = 0
    else
      last = mdays(month)
      if (month == 2 .and. is_leap_year(year)) last = 29
    end if
  end function last_day_of_month

  pure logical function valid_date(d) result(value)
    type(date_t), intent(in) :: d
    value = d%month >= 1 .and. d%month <= 12
    if (value) value = d%day >= 1 .and. d%day <= last_day_of_month(d%month, d%year)
  end function valid_date

  pure integer function date_ordinal(d) result(ordinal)
    type(date_t), intent(in) :: d
    integer :: y, era, yoe, doy, doe, m
    y = d%year
    m = d%month
    if (m <= 2) y = y - 1
    era = floor_div(y, 400)
    yoe = y - era*400
    if (m > 2) then
      doy = (153*(m - 3) + 2)/5 + d%day - 1
    else
      doy = (153*(m + 9) + 2)/5 + d%day - 1
    end if
    doe = yoe*365 + yoe/4 - yoe/100 + doy
    ordinal = era*146097 + doe + 719469
  end function date_ordinal

  pure subroutine ordinal_date(ordinal, d)
    integer, intent(in) :: ordinal
    type(date_t), intent(out) :: d
    integer :: z, era, doe, yoe, doy, mp
    z = ordinal - 719469
    era = floor_div(z, 146097)
    doe = z - era*146097
    yoe = (doe - doe/1460 + doe/36524 - doe/146096)/365
    d%year = yoe + era*400
    doy = doe - (365*yoe + yoe/4 - yoe/100)
    mp = (5*doy + 2)/153
    d%day = doy - (153*mp + 2)/5 + 1
    if (mp < 10) then
      d%month = mp + 3
    else
      d%month = mp - 9
    end if
    if (d%month <= 2) d%year = d%year + 1
  end subroutine ordinal_date

  pure integer function floor_div(a, b) result(q)
    integer, intent(in) :: a, b
    q = a/b
    if (modulo(a, b) /= 0 .and. ((a < 0) .neqv. (b < 0))) q = q - 1
  end function floor_div

  pure integer function daycount_actual(d1, d2) result(days)
    type(date_t), intent(in) :: d1, d2
    days = date_ordinal(d2) - date_ordinal(d1)
  end function daycount_actual

  pure integer function daycount_30_360(d1, d2, variant) result(days)
    type(date_t), intent(in) :: d1, d2
    character(len=*), intent(in), optional :: variant
    character(len=8) :: v
    integer :: dd1, dd2, mm1, mm2, yy1, yy2
    v = 'US'
    if (present(variant)) v = upper(trim(variant))
    dd1 = d1%day; dd2 = d2%day
    mm1 = d1%month - 1; mm2 = d2%month - 1
    yy1 = d1%year - 1900; yy2 = d2%year - 1900
    if (trim(v) == 'US' .and. dd2 == 31 .and. dd1 < 30) then
      dd2 = 1
      mm2 = mm2 + 1
    end if
    if (trim(v) == 'IT' .and. mm1 == 2 .and. dd1 > 27) dd1 = 30
    if (trim(v) == 'IT' .and. mm2 == 2 .and. dd2 > 27) dd2 = 30
    days = 360*(yy2 - yy1) + 30*(mm2 - mm1 - 1) + max(0, 30 - dd1) + min(30, dd2)
  end function daycount_30_360

  pure real(dp) function year_fraction(d1, d2, r1, r2, freq, convention) result(value)
    type(date_t), intent(in) :: d1, d2
    type(date_t), intent(in), optional :: r1, r2
    integer, intent(in), optional :: freq
    character(len=*), intent(in), optional :: convention
    type(date_t) :: ref1, ref2
    integer :: use_freq
    character(len=16) :: conv
    use_freq = 2
    if (present(freq)) use_freq = freq
    conv = '30/360'
    if (present(convention)) conv = upper(trim(convention))
    ref1 = d1; ref2 = d2
    if (present(r1)) ref1 = r1
    if (present(r2)) ref2 = r2
    select case (trim(conv))
    case ('ACT/ACT')
      value = real(daycount_actual(d1, d2), dp) / &
        (real(use_freq, dp)*real(daycount_actual(ref1, ref2), dp))
    case ('ACT/360')
      value = real(daycount_actual(d1, d2), dp)/360.0_dp
    case ('30/360E')
      value = real(daycount_30_360(d1, d2, 'EU'), dp)/360.0_dp
    case default
      value = real(daycount_30_360(d1, d2, 'US'), dp)/360.0_dp
    end select
  end function year_fraction

  pure function edate(from, months) result(out)
    type(date_t), intent(in) :: from
    integer, intent(in), optional :: months
    type(date_t) :: out
    integer :: nmonth, zero_month, new_zero, ld, new_ld
    nmonth = 1
    if (present(months)) nmonth = months
    zero_month = from%month - 1
    new_zero = modulo(zero_month + nmonth, 12)
    out%year = from%year + floor_div(zero_month + nmonth - new_zero, 12)
    out%month = new_zero + 1
    ld = last_day_of_month(from%month, from%year)
    new_ld = last_day_of_month(out%month, out%year)
    if (from%day == ld .or. from%day > new_ld) then
      out%day = new_ld
    else
      out%day = from%day
    end if
  end function edate

  pure logical function date_eq(a,b)
    type(date_t), intent(in) :: a,b
    date_eq = date_ordinal(a) == date_ordinal(b)
  end function
  pure logical function date_ne(a,b)
    type(date_t), intent(in) :: a,b
    date_ne = .not. date_eq(a,b)
  end function
  pure logical function date_lt(a,b)
    type(date_t), intent(in) :: a,b
    date_lt = date_ordinal(a) < date_ordinal(b)
  end function
  pure logical function date_le(a,b)
    type(date_t), intent(in) :: a,b
    date_le = date_ordinal(a) <= date_ordinal(b)
  end function
  pure logical function date_gt(a,b)
    type(date_t), intent(in) :: a,b
    date_gt = date_ordinal(a) > date_ordinal(b)
  end function
  pure logical function date_ge(a,b)
    type(date_t), intent(in) :: a,b
    date_ge = date_ordinal(a) >= date_ordinal(b)
  end function

  pure function upper(text) result(out)
    character(len=*), intent(in) :: text
    character(len=len(text)) :: out
    integer :: i, c
    do i=1,len(text)
      c=iachar(text(i:i))
      if(c>=iachar('a').and.c<=iachar('z')) then
        out(i:i)=achar(c-32)
      else
        out(i:i)=text(i:i)
      end if
    end do
  end function upper
end module jrvfinance_dates
