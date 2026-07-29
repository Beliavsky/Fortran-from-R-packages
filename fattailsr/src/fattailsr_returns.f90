! FatTailsR modern Fortran translation
! Copyright (C) 2014-2026 Patrice Kiener
! Licensed under GPL-2.0-only. See COPYING.
module fattailsr_returns
   use fattailsr_kinds, only : dp
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
   implicit none
   private
   public :: elevate, replace_nonfinite, price_returns

contains

   elemental pure function elevate(x, e) result(y)
      real(dp), intent(in) :: x, e
      real(dp) :: y
      y = (x + sqrt(x*x + e*e))/2.0_dp
   end function elevate

   pure subroutine replace_nonfinite(x, y)
      real(dp), intent(in) :: x(:)
      real(dp), intent(out) :: y(size(x))
      integer :: i, first
      first = 0
      do i = 1, size(x)
         if (ieee_is_finite(x(i))) then
            first = i
            exit
         end if
      end do
      if (first == 0) then
         y = 0.0_dp
         return
      end if
      y(1:first) = x(first)
      do i = first + 1, size(x)
         if (ieee_is_finite(x(i))) then
            y(i) = x(i)
         else
            y(i) = y(i-1)
         end if
      end do
   end subroutine replace_nonfinite

   pure subroutine price_returns(prices, returns, log_returns, multiplier, e, fill_missing, status)
      real(dp), intent(in) :: prices(:)
      real(dp), intent(out) :: returns(size(prices))
      logical, intent(in), optional :: log_returns, fill_missing
      real(dp), intent(in), optional :: multiplier, e
      integer, intent(out), optional :: status
      real(dp) :: x(size(prices)), scale, focal
      logical :: use_log, fill
      integer :: i

      use_log = .true.
      fill = .true.
      scale = 100.0_dp
      focal = -1.0_dp
      if (present(log_returns)) use_log = log_returns
      if (present(fill_missing)) fill = fill_missing
      if (present(multiplier)) scale = multiplier
      if (present(e)) focal = e
      if (present(status)) status = 0
      if (size(prices) == 0) return

      if (fill) then
         call replace_nonfinite(prices, x)
      else
         x = prices
      end if
      if (present(e)) x = elevate(x, focal)
      returns(1) = 0.0_dp
      do i = 2, size(x)
         if (.not. ieee_is_finite(x(i)) .or. .not. ieee_is_finite(x(i-1))) then
            returns(i) = 0.0_dp
            if (present(status)) status = 1
         else if (use_log) then
            if (x(i) <= 0.0_dp .or. x(i-1) <= 0.0_dp) then
               returns(i) = 0.0_dp
               if (present(status)) status = 2
            else
               returns(i) = scale*(log(x(i)) - log(x(i-1)))
            end if
         else
            if (abs(x(i-1)) <= tiny(1.0_dp)) then
               returns(i) = 0.0_dp
               if (present(status)) status = 3
            else
               returns(i) = scale*(x(i)-x(i-1))/x(i-1)
            end if
         end if
      end do
   end subroutine price_returns

end module fattailsr_returns
