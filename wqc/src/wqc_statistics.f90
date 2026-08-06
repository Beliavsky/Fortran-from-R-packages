! SPDX-License-Identifier: GPL-3.0-only
module wqc_statistics
   use wqc_kinds, only : dp
   implicit none
   private

   public :: mean_value
   public :: sample_sd
   public :: quantile_type7
   public :: quantile_type7_sorted
   public :: sort_real

contains

   pure function mean_value(x) result(value)
      real(dp), intent(in) :: x(:)
      real(dp) :: value

      if (size(x) == 0) then
         value = 0.0_dp
      else
         value = sum(x) / real(size(x), dp)
      end if
   end function mean_value


   pure function sample_sd(x) result(value)
      real(dp), intent(in) :: x(:)
      real(dp) :: value

      real(dp) :: xbar
      integer :: n

      n = size(x)
      if (n < 2) then
         value = 0.0_dp
         return
      end if
      xbar = mean_value(x)
      value = sqrt(sum((x - xbar)**2) / real(n - 1, dp))
   end function sample_sd


   function quantile_type7(x, probability) result(value)
      real(dp), intent(in) :: x(:)
      real(dp), intent(in) :: probability
      real(dp) :: value

      real(dp), allocatable :: work(:)

      allocate(work(size(x)))
      work = x
      call sort_real(work)
      value = quantile_type7_sorted(work, probability)
   end function quantile_type7


   pure function quantile_type7_sorted(sorted_x, probability) result(value)
      real(dp), intent(in) :: sorted_x(:)
      real(dp), intent(in) :: probability
      real(dp) :: value

      real(dp) :: fraction, h
      integer :: j, n

      n = size(sorted_x)
      if (n == 0) then
         value = 0.0_dp
      else if (probability <= 0.0_dp) then
         value = sorted_x(1)
      else if (probability >= 1.0_dp) then
         value = sorted_x(n)
      else
         h = 1.0_dp + real(n - 1, dp) * probability
         j = int(floor(h))
         fraction = h - real(j, dp)
         if (j >= n) then
            value = sorted_x(n)
         else
            value = (1.0_dp - fraction) * sorted_x(j) + fraction * sorted_x(j + 1)
         end if
      end if
   end function quantile_type7_sorted


   recursive subroutine sort_real(x)
      real(dp), intent(inout) :: x(:)

      real(dp) :: pivot, temp
      integer :: i, j, n

      n = size(x)
      if (n <= 1) return
      pivot = x((n + 1) / 2)
      i = 1
      j = n
      do
         do while (x(i) < pivot)
            i = i + 1
         end do
         do while (x(j) > pivot)
            j = j - 1
         end do
         if (i <= j) then
            temp = x(i)
            x(i) = x(j)
            x(j) = temp
            i = i + 1
            j = j - 1
         end if
         if (i > j) exit
      end do
      if (j > 1) call sort_real(x(:j))
      if (i < n) call sort_real(x(i:))
   end subroutine sort_real

end module wqc_statistics
