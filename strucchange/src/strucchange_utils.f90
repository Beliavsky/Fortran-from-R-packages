! SPDX-License-Identifier: GPL-2.0-only OR GPL-3.0-only
! Derived from the R package strucchange 1.6-0. See NOTICE.md and UPSTREAM.md.
module strucchange_utils
   use r_kinds, only : dp
   implicit none
   private
   public :: cumulative_sum
   public :: linear_interp
   public :: sample_standard_deviation
contains
   pure function cumulative_sum(x) result(y)
      real(dp), intent(in) :: x(:)
      real(dp) :: y(size(x))
      integer :: i

      if (size(x) == 0) return
      y(1) = x(1)
      do i = 2, size(x)
         y(i) = y(i - 1) + x(i)
      end do
   end function cumulative_sum

   pure real(dp) function sample_standard_deviation(x, df) result(value)
      real(dp), intent(in) :: x(:)
      integer, intent(in), optional :: df
      real(dp) :: center
      integer :: denominator

      if (size(x) == 0) then
         value = 0.0_dp
         return
      end if
      center = sum(x) / real(size(x), dp)
      denominator = size(x) - 1
      if (present(df)) denominator = df
      if (denominator <= 0) then
         value = 0.0_dp
      else
         value = sqrt(sum((x - center) ** 2) / real(denominator, dp))
      end if
   end function sample_standard_deviation

   pure real(dp) function linear_interp(x_grid, y_grid, x, clamp) result(y)
      real(dp), intent(in) :: x_grid(:), y_grid(:), x
      logical, intent(in), optional :: clamp
      logical :: use_clamp
      integer :: i, n

      n = size(x_grid)
      use_clamp = .true.
      if (present(clamp)) use_clamp = clamp
      if (n == 0 .or. size(y_grid) /= n) then
         y = 0.0_dp
         return
      end if
      if (n == 1) then
         y = y_grid(1)
         return
      end if
      if (x <= x_grid(1)) then
         if (use_clamp) then
            y = y_grid(1)
         else
            y = y_grid(1) + (x - x_grid(1)) * (y_grid(2) - y_grid(1)) / &
               (x_grid(2) - x_grid(1))
         end if
         return
      end if
      if (x >= x_grid(n)) then
         if (use_clamp) then
            y = y_grid(n)
         else
            y = y_grid(n - 1) + (x - x_grid(n - 1)) * &
               (y_grid(n) - y_grid(n - 1)) / (x_grid(n) - x_grid(n - 1))
         end if
         return
      end if
      do i = 1, n - 1
         if (x <= x_grid(i + 1)) then
            y = y_grid(i) + (x - x_grid(i)) * (y_grid(i + 1) - y_grid(i)) / &
               (x_grid(i + 1) - x_grid(i))
            return
         end if
      end do
      y = y_grid(n)
   end function linear_interp
end module strucchange_utils
