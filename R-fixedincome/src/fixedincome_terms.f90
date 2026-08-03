! SPDX-License-Identifier: MIT
module fixedincome_terms
   use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
   use fixedincome_kinds, only : dp
   use fixedincome_types, only : term_t, daycount_t, UNIT_DAY, UNIT_MONTH, UNIT_YEAR, &
      FI_OK, FI_INVALID_ARGUMENT, FI_SIZE_MISMATCH
   implicit none
   private

   public :: term, term_units, parse_term, daycount, dib
   public :: todays, tomonths, toyears, convert_term
   public :: unit_from_string, unit_name, shift, term_difference
   public :: gregorian_to_ordinal, ordinal_to_gregorian, offset_date

   interface term
      module procedure term_scalar
      module procedure term_vector
   end interface term

   interface shift
      module procedure shift_real
      module procedure shift_term
   end interface shift

contains

   pure integer function unit_from_string(name) result(unit_code)
      character(len=*), intent(in) :: name
      character(len=:), allocatable :: s
      s = lower(trim(adjustl(name)))
      if (len(s) > 1 .and. s(len(s):len(s)) == 's') s = s(:len(s)-1)
      select case (s)
      case ('day')
         unit_code = UNIT_DAY
      case ('month')
         unit_code = UNIT_MONTH
      case ('year')
         unit_code = UNIT_YEAR
      case default
         unit_code = 0
      end select
   end function unit_from_string

   pure function unit_name(unit_code, plural) result(name)
      integer, intent(in) :: unit_code
      logical, intent(in), optional :: plural
      character(len=8) :: name
      logical :: use_plural
      use_plural = .false.
      if (present(plural)) use_plural = plural
      select case (unit_code)
      case (UNIT_DAY)
         name = 'day'
      case (UNIT_MONTH)
         name = 'month'
      case (UNIT_YEAR)
         name = 'year'
      case default
         name = 'unknown'
      end select
      if (use_plural .and. trim(name) /= 'unknown') name = trim(name)//'s'
   end function unit_name

   function term_scalar(value, units, status) result(t)
      real(dp), intent(in) :: value
      character(len=*), intent(in) :: units
      integer, intent(out), optional :: status
      type(term_t) :: t
      integer :: code
      code = unit_from_string(units)
      allocate(t%value(1), t%unit(1))
      t%value = value
      t%unit = code
      if (present(status)) then
         if (code == 0) then
            status = FI_INVALID_ARGUMENT
         else
            status = FI_OK
         end if
      end if
   end function term_scalar

   function term_vector(value, units, status) result(t)
      real(dp), intent(in) :: value(:)
      character(len=*), intent(in) :: units
      integer, intent(out), optional :: status
      type(term_t) :: t
      integer :: code
      code = unit_from_string(units)
      allocate(t%value(size(value)), t%unit(size(value)))
      t%value = value
      t%unit = code
      if (present(status)) then
         if (code == 0) then
            status = FI_INVALID_ARGUMENT
         else
            status = FI_OK
         end if
      end if
   end function term_vector

   function term_units(value, units, status) result(t)
      real(dp), intent(in) :: value(:)
      character(len=*), intent(in) :: units(:)
      integer, intent(out), optional :: status
      type(term_t) :: t
      integer :: i
      allocate(t%value(size(value)), t%unit(size(value)))
      t%value = value
      if (size(value) /= size(units)) then
         t%unit = 0
         if (present(status)) status = FI_SIZE_MISMATCH
         return
      end if
      do i = 1, size(value)
         t%unit(i) = unit_from_string(units(i))
      end do
      if (present(status)) then
         if (any(t%unit == 0)) then
            status = FI_INVALID_ARGUMENT
         else
            status = FI_OK
         end if
      end if
   end function term_units

   function parse_term(text, status) result(t)
      character(len=*), intent(in) :: text
      integer, intent(out), optional :: status
      type(term_t) :: t
      character(len=64) :: units
      real(dp) :: value
      integer :: ios
      read(text, *, iostat=ios) value, units
      if (ios /= 0) then
         allocate(t%value(1), t%unit(1))
         t%value = ieee_value(0.0_dp, ieee_quiet_nan)
         t%unit = 0
         if (present(status)) status = FI_INVALID_ARGUMENT
         return
      end if
      t = term(value, trim(units), status)
   end function parse_term

   function daycount(specification, status) result(dc)
      character(len=*), intent(in) :: specification
      integer, intent(out), optional :: status
      type(daycount_t) :: dc
      integer :: slash, ios, base
      slash = index(trim(specification), '/')
      if (slash <= 1 .or. slash >= len_trim(specification)) then
         if (present(status)) status = FI_INVALID_ARGUMENT
         dc%specification = trim(specification)
         dc%days_in_base = 0
         return
      end if
      read(specification(slash+1:), *, iostat=ios) base
      if (ios /= 0 .or. base <= 0) then
         if (present(status)) status = FI_INVALID_ARGUMENT
         dc%specification = trim(specification)
         dc%days_in_base = 0
         return
      end if
      dc%specification = trim(specification)
      dc%days_in_base = base
      if (present(status)) status = FI_OK
   end function daycount

   pure integer function dib(dc) result(base)
      type(daycount_t), intent(in) :: dc
      base = dc%days_in_base
   end function dib

   function convert_term(dc, input, target_unit, status) result(output)
      type(daycount_t), intent(in) :: dc
      type(term_t), intent(in) :: input
      integer, intent(in) :: target_unit
      integer, intent(out), optional :: status
      type(term_t) :: output
      integer :: i
      real(dp) :: days_per_month, days_value
      allocate(output%value(input%size()), output%unit(input%size()))
      if (dc%days_in_base <= 0 .or. target_unit < UNIT_DAY .or. target_unit > UNIT_YEAR) then
         output%value = ieee_value(0.0_dp, ieee_quiet_nan)
         output%unit = 0
         if (present(status)) status = FI_INVALID_ARGUMENT
         return
      end if
      days_per_month = real(dc%days_in_base, dp) / 12.0_dp
      do i = 1, input%size()
         select case (input%unit(i))
         case (UNIT_DAY)
            days_value = input%value(i)
         case (UNIT_MONTH)
            days_value = input%value(i) * days_per_month
         case (UNIT_YEAR)
            days_value = input%value(i) * real(dc%days_in_base, dp)
         case default
            output%value(i) = ieee_value(0.0_dp, ieee_quiet_nan)
            output%unit(i) = 0
            cycle
         end select
         select case (target_unit)
         case (UNIT_DAY)
            output%value(i) = days_value
         case (UNIT_MONTH)
            output%value(i) = days_value / days_per_month
         case (UNIT_YEAR)
            output%value(i) = days_value / real(dc%days_in_base, dp)
         end select
         output%unit(i) = target_unit
      end do
      if (present(status)) then
         if (any(output%unit == 0)) then
            status = FI_INVALID_ARGUMENT
         else
            status = FI_OK
         end if
      end if
   end function convert_term

   function todays(dc, input, status) result(output)
      type(daycount_t), intent(in) :: dc
      type(term_t), intent(in) :: input
      integer, intent(out), optional :: status
      type(term_t) :: output
      output = convert_term(dc, input, UNIT_DAY, status)
   end function todays

   function tomonths(dc, input, status) result(output)
      type(daycount_t), intent(in) :: dc
      type(term_t), intent(in) :: input
      integer, intent(out), optional :: status
      type(term_t) :: output
      output = convert_term(dc, input, UNIT_MONTH, status)
   end function tomonths

   function toyears(dc, input, status) result(output)
      type(daycount_t), intent(in) :: dc
      type(term_t), intent(in) :: input
      integer, intent(out), optional :: status
      type(term_t) :: output
      output = convert_term(dc, input, UNIT_YEAR, status)
   end function toyears

   function shift_real(x, k, fill, status) result(y)
      real(dp), intent(in) :: x(:)
      integer, intent(in), optional :: k
      real(dp), intent(in), optional :: fill
      integer, intent(out), optional :: status
      real(dp), allocatable :: y(:)
      integer :: nshift, n
      real(dp) :: fill_value
      n = size(x)
      nshift = 1
      if (present(k)) nshift = k
      fill_value = ieee_value(0.0_dp, ieee_quiet_nan)
      if (present(fill)) fill_value = fill
      allocate(y(n))
      if (abs(nshift) > n) then
         y = fill_value
         if (present(status)) status = FI_INVALID_ARGUMENT
         return
      end if
      if (nshift == 0) then
         y = x
      else if (nshift > 0) then
         y(:nshift) = fill_value
         if (nshift < n) y(nshift+1:) = x(:n-nshift)
      else
         nshift = abs(nshift)
         if (nshift < n) y(:n-nshift) = x(nshift+1:)
         y(n-nshift+1:) = fill_value
      end if
      if (present(status)) status = FI_OK
   end function shift_real

   function shift_term(x, k, fill, status) result(y)
      type(term_t), intent(in) :: x
      integer, intent(in), optional :: k
      real(dp), intent(in), optional :: fill
      integer, intent(out), optional :: status
      type(term_t) :: y
      integer :: nshift, n, fill_unit
      real(dp) :: fill_value
      n = x%size()
      nshift = 1
      if (present(k)) nshift = k
      fill_value = ieee_value(0.0_dp, ieee_quiet_nan)
      if (present(fill)) fill_value = fill
      allocate(y%value(n), y%unit(n))
      if (n > 0) then
         fill_unit = x%unit(1)
      else
         fill_unit = UNIT_DAY
      end if
      if (abs(nshift) > n) then
         y%value = fill_value
         y%unit = fill_unit
         if (present(status)) status = FI_INVALID_ARGUMENT
         return
      end if
      if (nshift == 0) then
         y = x
      else if (nshift > 0) then
         y%value(:nshift) = fill_value
         y%unit(:nshift) = fill_unit
         if (nshift < n) then
            y%value(nshift+1:) = x%value(:n-nshift)
            y%unit(nshift+1:) = x%unit(:n-nshift)
         end if
      else
         nshift = abs(nshift)
         if (nshift < n) then
            y%value(:n-nshift) = x%value(nshift+1:)
            y%unit(:n-nshift) = x%unit(nshift+1:)
         end if
         y%value(n-nshift+1:) = fill_value
         y%unit(n-nshift+1:) = fill_unit
      end if
      if (present(status)) status = FI_OK
   end function shift_term

   function term_difference(x, lag, fill, status) result(y)
      type(term_t), intent(in) :: x
      integer, intent(in), optional :: lag
      logical, intent(in), optional :: fill
      integer, intent(out), optional :: status
      type(term_t) :: y
      integer :: use_lag, n
      logical :: use_fill
      real(dp) :: nan_value
      use_lag = 1
      if (present(lag)) use_lag = lag
      use_fill = .false.
      if (present(fill)) use_fill = fill
      n = x%size()
      if (use_lag < 1 .or. use_lag >= n) then
         allocate(y%value(0), y%unit(0))
         if (present(status)) status = FI_INVALID_ARGUMENT
         return
      end if
      if (any(x%unit /= x%unit(1))) then
         allocate(y%value(0), y%unit(0))
         if (present(status)) status = FI_INVALID_ARGUMENT
         return
      end if
      if (use_fill) then
         allocate(y%value(n), y%unit(n))
         nan_value = ieee_value(0.0_dp, ieee_quiet_nan)
         y%value(:use_lag) = nan_value
         y%value(use_lag+1:) = x%value(use_lag+1:) - x%value(:n-use_lag)
         y%unit = x%unit(1)
      else
         allocate(y%value(n-use_lag), y%unit(n-use_lag))
         y%value = x%value(use_lag+1:) - x%value(:n-use_lag)
         y%unit = x%unit(1)
      end if
      if (present(status)) status = FI_OK
   end function term_difference

   pure integer function gregorian_to_ordinal(year, month, day) result(ordinal)
      integer, intent(in) :: year, month, day
      integer :: y, era, yoe, doy, doe, m
      y = year
      m = month
      if (m <= 2) y = y - 1
      era = floor_div(y, 400)
      yoe = y - era * 400
      if (m > 2) then
         doy = (153 * (m - 3) + 2) / 5 + day - 1
      else
         doy = (153 * (m + 9) + 2) / 5 + day - 1
      end if
      doe = yoe * 365 + yoe / 4 - yoe / 100 + doy
      ordinal = era * 146097 + doe + 719469
   end function gregorian_to_ordinal

   pure subroutine ordinal_to_gregorian(ordinal, year, month, day)
      integer, intent(in) :: ordinal
      integer, intent(out) :: year, month, day
      integer :: z, era, doe, yoe, doy, mp
      z = ordinal - 719469
      era = floor_div(z, 146097)
      doe = z - era * 146097
      yoe = (doe - doe / 1460 + doe / 36524 - doe / 146096) / 365
      year = yoe + era * 400
      doy = doe - (365 * yoe + yoe / 4 - yoe / 100)
      mp = (5 * doy + 2) / 153
      day = doy - (153 * mp + 2) / 5 + 1
      if (mp < 10) then
         month = mp + 3
      else
         month = mp - 9
      end if
      if (month <= 2) year = year + 1
   end subroutine ordinal_to_gregorian

   function offset_date(reference_date, days, calendar, status) result(date_out)
      integer, intent(in) :: reference_date
      integer, intent(in) :: days
      character(len=*), intent(in), optional :: calendar
      integer, intent(out), optional :: status
      integer :: date_out, remaining, direction
      character(len=64) :: cal
      cal = 'actual'
      if (present(calendar)) cal = lower(trim(calendar))
      if (cal == 'actual') then
         date_out = reference_date + days
         if (present(status)) status = FI_OK
         return
      end if
      if (cal /= 'weekdays' .and. cal /= 'business') then
         date_out = reference_date
         if (present(status)) status = FI_INVALID_ARGUMENT
         return
      end if
      date_out = reference_date
      remaining = abs(days)
      direction = merge(1, -1, days >= 0)
      do while (remaining > 0)
         date_out = date_out + direction
         if (weekday(date_out) <= 5) remaining = remaining - 1
      end do
      if (present(status)) status = FI_OK
   end function offset_date

   pure integer function weekday(ordinal) result(wd)
      integer, intent(in) :: ordinal
      wd = modulo(ordinal, 7) + 1
   end function weekday

   pure integer function floor_div(a, b) result(q)
      integer, intent(in) :: a, b
      q = a / b
      if (mod(a, b) /= 0 .and. ((a < 0) .neqv. (b < 0))) q = q - 1
   end function floor_div

   pure function lower(text) result(out)
      character(len=*), intent(in) :: text
      character(len=len(text)) :: out
      integer :: i, c
      do i = 1, len(text)
         c = iachar(text(i:i))
         if (c >= iachar('A') .and. c <= iachar('Z')) then
            out(i:i) = achar(c + 32)
         else
            out(i:i) = text(i:i)
         end if
      end do
   end function lower

end module fixedincome_terms
