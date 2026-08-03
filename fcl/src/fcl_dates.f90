! MIT License. Copyright (c) 2024 fcl authors.
module fcl_dates
   use fcl_kinds, only : dp
   implicit none
   private
   public :: date_type, make_date, add_months, year_frac, days_between

   type :: date_type
      integer :: year = 1970
      integer :: month = 1
      integer :: day = 1
   contains
      procedure :: serial => date_serial
   end type date_type

contains

   pure function make_date(year, month, day) result(date)
      integer, intent(in) :: year, month, day
      type(date_type) :: date
      date%year = year
      date%month = month
      date%day = day
   end function make_date

   pure integer function days_in_month(year, month) result(n)
      integer, intent(in) :: year, month
      integer, parameter :: mdays(12) = [31,28,31,30,31,30,31,31,30,31,30,31]
      n = mdays(month)
      if (month == 2 .and. is_leap_year(year)) n = 29
   end function days_in_month

   pure logical function is_leap_year(year)
      integer, intent(in) :: year
      is_leap_year = mod(year, 400) == 0 .or. &
         (mod(year, 4) == 0 .and. mod(year, 100) /= 0)
   end function is_leap_year

   pure function add_months(ref_date, months) result(date)
      type(date_type), intent(in) :: ref_date
      integer, intent(in) :: months
      type(date_type) :: date
      integer :: total, year, month, day
      total = ref_date%year * 12 + ref_date%month + months
      year = (total - 1) / 12
      month = modulo(total - 1, 12) + 1
      day = min(ref_date%day, days_in_month(year, month))
      date = make_date(year, month, day)
   end function add_months

   pure real(dp) function year_frac(d1, d0) result(frac)
      type(date_type), intent(in) :: d1, d0
      frac = real(d1%year - d0%year, dp) + &
         real(d1%month - d0%month, dp) / 12.0_dp + &
         real(d1%day - d0%day, dp) / 365.0_dp
   end function year_frac

   pure integer function date_serial(self) result(jdn)
      class(date_type), intent(in) :: self
      integer :: a, y, m
      a = (14 - self%month) / 12
      y = self%year + 4800 - a
      m = self%month + 12 * a - 3
      jdn = self%day + (153 * m + 2) / 5 + 365 * y + y / 4 - y / 100 + y / 400 - 32045
   end function date_serial

   pure integer function days_between(d1, d0) result(days)
      type(date_type), intent(in) :: d1, d0
      days = d1%serial() - d0%serial()
   end function days_between

end module fcl_dates
