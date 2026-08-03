! SPDX-License-Identifier: GPL-2.0-or-later
module fints_dates
   use fints_kinds, only : dp
   use fints_status, only : fints_ok, fints_invalid_input
   use fints_types, only : yearmon_result
   implicit none
   private
   public :: as_yearmon2

   interface as_yearmon2
      module procedure as_yearmon2_real
      module procedure as_yearmon2_integer
   end interface as_yearmon2

contains

   subroutine as_yearmon2_real(x, result)
      real(dp), intent(in) :: x(:)
      type(yearmon_result), intent(out) :: result
      integer :: i, year_value, month_value

      result = yearmon_result()
      allocate(result%year(size(x)), result%month(size(x)))
      do i = 1, size(x)
         if (abs(100.0_dp * x(i) - anint(100.0_dp * x(i))) > &
            100.0_dp * epsilon(1.0_dp) * max(1.0_dp, abs(x(i)))) then
            result%status = fints_invalid_input
            return
         end if
         year_value = floor(x(i))
         month_value = nint(100.0_dp * (x(i) - real(year_value, dp)))
         if (month_value < 1 .or. month_value > 12) then
            result%status = fints_invalid_input
            return
         end if
         result%year(i) = year_value
         result%month(i) = month_value
      end do
      call count_duplicates(result)
      result%status = fints_ok
   end subroutine as_yearmon2_real

   subroutine as_yearmon2_integer(yyyymm, result)
      integer, intent(in) :: yyyymm(:)
      type(yearmon_result), intent(out) :: result
      integer :: i

      result = yearmon_result()
      allocate(result%year(size(yyyymm)), result%month(size(yyyymm)))
      do i = 1, size(yyyymm)
         result%year(i) = yyyymm(i) / 100
         result%month(i) = modulo(yyyymm(i), 100)
         if (result%month(i) < 1 .or. result%month(i) > 12) then
            result%status = fints_invalid_input
            return
         end if
      end do
      call count_duplicates(result)
      result%status = fints_ok
   end subroutine as_yearmon2_integer

   subroutine count_duplicates(result)
      type(yearmon_result), intent(inout) :: result
      integer :: i, j

      result%duplicate_count = 0
      do i = 1, size(result%year)
         do j = 1, i - 1
            if (result%year(i) == result%year(j) .and. &
               result%month(i) == result%month(j)) then
               result%duplicate_count = result%duplicate_count + 1
               exit
            end if
         end do
      end do
      result%converted = result%duplicate_count == 0
   end subroutine count_duplicates

end module fints_dates
