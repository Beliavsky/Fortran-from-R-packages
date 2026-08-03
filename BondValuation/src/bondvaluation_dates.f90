! SPDX-License-Identifier: GPL-3.0-only
module bondvaluation_dates
  use bondvaluation_kinds, only: dp
  implicit none
  private

  type, public :: date_type
    integer :: year = 0
    integer :: month = 0
    integer :: day = 0
  end type date_type

  public :: date_from_ymd, valid_date, leap_year, days_in_month, days_in_year
  public :: date_to_serial, serial_to_date, day_diff, compare_dates
  public :: add_months, last_day_of_month, is_last_day_of_month
  public :: date_to_string, weekday_number, previous_date, next_date
  public :: sort_dates, unique_sorted_dates, find_previous_date, find_next_date

contains

  pure function date_from_ymd(year, month, day) result(date)
    integer, intent(in) :: year, month, day
    type(date_type) :: date
    date%year = year
    date%month = month
    date%day = day
  end function date_from_ymd

  pure logical function leap_year(year) result(ans)
    integer, intent(in) :: year
    ans = modulo(year, 400) == 0 .or. &
          (modulo(year, 4) == 0 .and. modulo(year, 100) /= 0)
  end function leap_year

  pure integer function days_in_month(year, month) result(n)
    integer, intent(in) :: year, month
    integer, parameter :: mdays(12) = [31, 28, 31, 30, 31, 30, &
                                       31, 31, 30, 31, 30, 31]
    if (month < 1 .or. month > 12) then
      n = 0
    else
      n = mdays(month)
      if (month == 2 .and. leap_year(year)) n = 29
    end if
  end function days_in_month

  pure integer function days_in_year(year) result(n)
    integer, intent(in) :: year
    if (leap_year(year)) then
      n = 366
    else
      n = 365
    end if
  end function days_in_year

  pure logical function valid_date(date) result(ans)
    type(date_type), intent(in) :: date
    ans = date%year /= 0 .and. date%month >= 1 .and. date%month <= 12
    if (ans) ans = date%day >= 1 .and. &
                   date%day <= days_in_month(date%year, date%month)
  end function valid_date

  pure integer function date_to_serial(date) result(serial)
    ! Days since 1970-01-01, matching R's Date representation.
    type(date_type), intent(in) :: date
    integer :: y, m, era, yoe, doy, doe
    if (.not. valid_date(date)) then
      serial = -huge(0)
      return
    end if
    y = date%year
    m = date%month
    if (m <= 2) y = y - 1
    if (y >= 0) then
      era = y / 400
    else
      era = (y - 399) / 400
    end if
    yoe = y - era * 400
    if (m > 2) then
      doy = (153 * (m - 3) + 2) / 5 + date%day - 1
    else
      doy = (153 * (m + 9) + 2) / 5 + date%day - 1
    end if
    doe = yoe * 365 + yoe / 4 - yoe / 100 + doy
    serial = era * 146097 + doe - 719468
  end function date_to_serial

  pure function serial_to_date(serial) result(date)
    integer, intent(in) :: serial
    type(date_type) :: date
    integer :: z, era, doe, yoe, y, doy, mp
    z = serial + 719468
    if (z >= 0) then
      era = z / 146097
    else
      era = (z - 146096) / 146097
    end if
    doe = z - era * 146097
    yoe = (doe - doe / 1460 + doe / 36524 - doe / 146096) / 365
    y = yoe + era * 400
    doy = doe - (365 * yoe + yoe / 4 - yoe / 100)
    mp = (5 * doy + 2) / 153
    date%day = doy - (153 * mp + 2) / 5 + 1
    if (mp < 10) then
      date%month = mp + 3
    else
      date%month = mp - 9
    end if
    if (date%month <= 2) y = y + 1
    date%year = y
  end function serial_to_date

  pure integer function day_diff(first_date, second_date) result(n)
    type(date_type), intent(in) :: first_date, second_date
    n = date_to_serial(second_date) - date_to_serial(first_date)
  end function day_diff

  pure integer function compare_dates(a, b) result(cmp)
    type(date_type), intent(in) :: a, b
    integer :: sa, sb
    sa = date_to_serial(a)
    sb = date_to_serial(b)
    if (sa < sb) then
      cmp = -1
    else if (sa > sb) then
      cmp = 1
    else
      cmp = 0
    end if
  end function compare_dates

  pure function last_day_of_month(date) result(out)
    type(date_type), intent(in) :: date
    type(date_type) :: out
    out = date
    if (date%year /= 0 .and. date%month >= 1 .and. date%month <= 12) &
      out%day = days_in_month(date%year, date%month)
  end function last_day_of_month

  pure logical function is_last_day_of_month(date) result(ans)
    type(date_type), intent(in) :: date
    ans = valid_date(date) .and. &
          date%day == days_in_month(date%year, date%month)
  end function is_last_day_of_month

  pure function add_months(date, number_months, end_of_month, reference_day) result(out)
    type(date_type), intent(in) :: date
    integer, intent(in) :: number_months
    logical, intent(in), optional :: end_of_month
    integer, intent(in), optional :: reference_day
    type(date_type) :: out
    integer :: total, y, m, target_day
    logical :: use_eom

    if (.not. valid_date(date)) then
      out = date_type()
      return
    end if
    total = date%year * 12 + date%month - 1 + number_months
    if (total >= 0) then
      y = total / 12
    else
      y = (total - 11) / 12
    end if
    m = total - y * 12 + 1
    use_eom = .false.
    if (present(end_of_month)) use_eom = end_of_month
    if (present(reference_day)) then
      target_day = reference_day
    else
      target_day = date%day
    end if
    out%year = y
    out%month = m
    if (use_eom) then
      out%day = days_in_month(y, m)
    else
      out%day = min(target_day, days_in_month(y, m))
    end if
  end function add_months

  pure function previous_date(date) result(out)
    type(date_type), intent(in) :: date
    type(date_type) :: out
    out = serial_to_date(date_to_serial(date) - 1)
  end function previous_date

  pure function next_date(date) result(out)
    type(date_type), intent(in) :: date
    type(date_type) :: out
    out = serial_to_date(date_to_serial(date) + 1)
  end function next_date

  pure integer function weekday_number(date) result(day_number)
    ! Monday=1, ..., Sunday=7.
    type(date_type), intent(in) :: date
    day_number = modulo(date_to_serial(date) + 3, 7) + 1
  end function weekday_number

  function date_to_string(date) result(text)
    type(date_type), intent(in) :: date
    character(len=10) :: text
    if (valid_date(date)) then
      write(text, '(i4.4,"-",i2.2,"-",i2.2)') date%year, date%month, date%day
    else
      text = 'NA        '
    end if
  end function date_to_string

  subroutine sort_dates(dates)
    type(date_type), intent(inout) :: dates(:)
    integer :: i, j
    type(date_type) :: key
    do i = 2, size(dates)
      key = dates(i)
      j = i - 1
      do while (j >= 1)
        if (compare_dates(dates(j), key) <= 0) exit
        dates(j + 1) = dates(j)
        j = j - 1
      end do
      dates(j + 1) = key
    end do
  end subroutine sort_dates

  subroutine unique_sorted_dates(input, output)
    type(date_type), intent(in) :: input(:)
    type(date_type), allocatable, intent(out) :: output(:)
    type(date_type), allocatable :: work(:)
    integer :: i, n
    if (size(input) == 0) then
      allocate(output(0))
      return
    end if
    work = input
    call sort_dates(work)
    n = 1
    do i = 2, size(work)
      if (compare_dates(work(i), work(n)) /= 0) then
        n = n + 1
        work(n) = work(i)
      end if
    end do
    allocate(output(n))
    output = work(1:n)
  end subroutine unique_sorted_dates

  integer function find_previous_date(date, dates) result(index)
    type(date_type), intent(in) :: date
    type(date_type), intent(in) :: dates(:)
    integer :: i
    index = 0
    do i = 1, size(dates)
      if (compare_dates(dates(i), date) <= 0) then
        index = i
      else
        exit
      end if
    end do
  end function find_previous_date

  integer function find_next_date(date, dates) result(index)
    type(date_type), intent(in) :: date
    type(date_type), intent(in) :: dates(:)
    integer :: i
    index = 0
    do i = 1, size(dates)
      if (compare_dates(dates(i), date) > 0) then
        index = i
        return
      end if
    end do
  end function find_next_date

end module bondvaluation_dates
