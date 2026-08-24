! Date support adapted from Beliavsky/DataFrame under the MIT License.
module comparison_date
   implicit none
   private

   public :: date, valid, date_from_iso, date_from_basic
   public :: day_of_week, day_of_year, quarter, is_month_end
   public :: operator(+), operator(-), operator(==), operator(/=)
   public :: operator(<), operator(<=), operator(>), operator(>=)

   type :: date
      integer :: year = 0
      integer :: month = 0
      integer :: day = 0
   contains
      procedure :: to_string
   end type date

   interface operator(+)
      module procedure add_days_right
      module procedure add_days_left
   end interface

   interface operator(-)
      module procedure subtract_days
      module procedure difference_days
   end interface

   interface operator(==)
      module procedure equal_dates
   end interface

   interface operator(/=)
      module procedure unequal_dates
   end interface

   interface operator(<)
      module procedure earlier_date
   end interface

   interface operator(<=)
      module procedure earlier_or_equal_date
   end interface

   interface operator(>)
      module procedure later_date
   end interface

   interface operator(>=)
      module procedure later_or_equal_date
   end interface

contains

   pure function to_string(this) result(text)
      class(date), intent(in) :: this
      character(len=10) :: text

      write(text, '(i4.4,"-",i2.2,"-",i2.2)') this%year, this%month, this%day
   end function to_string

   pure elemental logical function valid(value)
      type(date), intent(in) :: value

      valid = value%month >= 1 .and. value%month <= 12
      if (valid) valid = value%day >= 1 .and. &
         value%day <= days_in_month(value%year, value%month)
   end function valid

   pure elemental integer function days_in_month(year, month)
      integer, intent(in) :: year, month

      select case (month)
      case (1, 3, 5, 7, 8, 10, 12)
         days_in_month = 31
      case (4, 6, 9, 11)
         days_in_month = 30
      case (2)
         if ((mod(year, 4) == 0 .and. mod(year, 100) /= 0) .or. &
             mod(year, 400) == 0) then
            days_in_month = 29
         else
            days_in_month = 28
         end if
      case default
         days_in_month = 0
      end select
   end function days_in_month

   pure function date_from_iso(text) result(value)
      character(len=*), intent(in) :: text
      type(date) :: value
      character(len=:), allocatable :: trimmed
      integer :: io

      value = date()
      trimmed = trim(adjustl(text))
      if (len(trimmed) /= 10) return
      if (trimmed(5:5) /= '-' .or. trimmed(8:8) /= '-') return
      read(trimmed(1:4), '(i4)', iostat=io) value%year
      if (io /= 0) then
         value = date()
         return
      end if
      read(trimmed(6:7), '(i2)', iostat=io) value%month
      if (io /= 0) then
         value = date()
         return
      end if
      read(trimmed(9:10), '(i2)', iostat=io) value%day
      if (io /= 0 .or. .not. valid(value)) value = date()
   end function date_from_iso

   pure function date_from_basic(text) result(value)
      character(len=*), intent(in) :: text
      type(date) :: value
      character(len=:), allocatable :: trimmed
      integer :: io

      value = date()
      trimmed = trim(adjustl(text))
      if (len(trimmed) /= 8) return
      read(trimmed(1:4), '(i4)', iostat=io) value%year
      if (io /= 0) then
         value = date()
         return
      end if
      read(trimmed(5:6), '(i2)', iostat=io) value%month
      if (io /= 0) then
         value = date()
         return
      end if
      read(trimmed(7:8), '(i2)', iostat=io) value%day
      if (io /= 0 .or. .not. valid(value)) value = date()
   end function date_from_basic

   pure elemental integer function day_of_week(value)
      ! ISO convention: Monday=1, ..., Sunday=7.
      type(date), intent(in) :: value

      if (valid(value)) then
         day_of_week = modulo(day_number(value) + 3, 7) + 1
      else
         day_of_week = 0
      end if
   end function day_of_week

   pure elemental integer function day_of_year(value)
      type(date), intent(in) :: value

      if (valid(value)) then
         day_of_year = day_number(value) - day_number(date(value%year, 1, 1)) + 1
      else
         day_of_year = 0
      end if
   end function day_of_year

   pure elemental integer function quarter(value)
      type(date), intent(in) :: value

      if (valid(value)) then
         quarter = (value%month - 1)/3 + 1
      else
         quarter = 0
      end if
   end function quarter

   pure elemental logical function is_month_end(value)
      type(date), intent(in) :: value

      is_month_end = valid(value) .and. &
         value%day == days_in_month(value%year, value%month)
   end function is_month_end

   pure elemental type(date) function add_days_right(value, number_days)
      type(date), intent(in) :: value
      integer, intent(in) :: number_days

      if (valid(value)) then
         add_days_right = from_day_number(day_number(value) + number_days)
      else
         add_days_right = date()
      end if
   end function add_days_right

   pure elemental type(date) function add_days_left(number_days, value)
      integer, intent(in) :: number_days
      type(date), intent(in) :: value

      add_days_left = add_days_right(value, number_days)
   end function add_days_left

   pure elemental type(date) function subtract_days(value, number_days)
      type(date), intent(in) :: value
      integer, intent(in) :: number_days

      subtract_days = add_days_right(value, -number_days)
   end function subtract_days

   pure elemental integer function difference_days(left, right)
      type(date), intent(in) :: left, right

      if (valid(left) .and. valid(right)) then
         difference_days = day_number(left) - day_number(right)
      else
         difference_days = 0
      end if
   end function difference_days

   pure elemental logical function equal_dates(left, right)
      type(date), intent(in) :: left, right

      equal_dates = left%year == right%year .and. &
         left%month == right%month .and. left%day == right%day
   end function equal_dates

   pure elemental logical function unequal_dates(left, right)
      type(date), intent(in) :: left, right

      unequal_dates = .not. equal_dates(left, right)
   end function unequal_dates

   pure elemental logical function earlier_date(left, right)
      type(date), intent(in) :: left, right

      earlier_date = left%year < right%year .or. &
         (left%year == right%year .and. &
          (left%month < right%month .or. &
           (left%month == right%month .and. left%day < right%day)))
   end function earlier_date

   pure elemental logical function earlier_or_equal_date(left, right)
      type(date), intent(in) :: left, right

      earlier_or_equal_date = earlier_date(left, right) .or. equal_dates(left, right)
   end function earlier_or_equal_date

   pure elemental logical function later_date(left, right)
      type(date), intent(in) :: left, right

      later_date = .not. earlier_or_equal_date(left, right)
   end function later_date

   pure elemental logical function later_or_equal_date(left, right)
      type(date), intent(in) :: left, right

      later_or_equal_date = .not. earlier_date(left, right)
   end function later_or_equal_date

   pure elemental integer function day_number(value)
      type(date), intent(in) :: value
      integer :: year, month, era, year_of_era, day_of_year_zero, day_of_era

      year = value%year
      month = value%month
      if (month <= 2) year = year - 1
      era = floor_div(year, 400)
      year_of_era = year - era*400
      if (month > 2) then
         day_of_year_zero = (153*(month - 3) + 2)/5 + value%day - 1
      else
         day_of_year_zero = (153*(month + 9) + 2)/5 + value%day - 1
      end if
      day_of_era = year_of_era*365 + year_of_era/4 - year_of_era/100 + &
         day_of_year_zero
      day_number = era*146097 + day_of_era - 719468
   end function day_number

   pure elemental type(date) function from_day_number(number)
      integer, intent(in) :: number
      integer :: adjusted, era, day_of_era, year_of_era, year
      integer :: day_of_year_zero, month_prime, month, day

      adjusted = number + 719468
      era = floor_div(adjusted, 146097)
      day_of_era = adjusted - era*146097
      year_of_era = (day_of_era - day_of_era/1460 + day_of_era/36524 - &
         day_of_era/146096)/365
      year = year_of_era + era*400
      day_of_year_zero = day_of_era - &
         (365*year_of_era + year_of_era/4 - year_of_era/100)
      month_prime = (5*day_of_year_zero + 2)/153
      day = day_of_year_zero - (153*month_prime + 2)/5 + 1
      if (month_prime < 10) then
         month = month_prime + 3
      else
         month = month_prime - 9
      end if
      if (month <= 2) year = year + 1
      from_day_number = date(year, month, day)
   end function from_day_number

   pure elemental integer function floor_div(numerator, denominator)
      integer, intent(in) :: numerator, denominator

      floor_div = numerator/denominator
      if (mod(numerator, denominator) < 0) floor_div = floor_div - 1
   end function floor_div

end module comparison_date
