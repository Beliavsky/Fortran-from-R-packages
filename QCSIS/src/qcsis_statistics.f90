module qcsis_statistics
   use qcsis_kinds, only : dp
   implicit none
   private

   public :: mean_value
   public :: sample_sd
   public :: quantile_type7
   public :: quantile_type7_sorted
   public :: descending_order
   public :: sort_real

contains

   pure function mean_value(x) result(xbar)
      real(dp), intent(in) :: x(:)
      real(dp) :: xbar

      if (size(x) == 0) then
         xbar = 0.0_dp
      else
         xbar = sum(x) / real(size(x), dp)
      end if
   end function mean_value


   pure function sample_sd(x) result(s)
      real(dp), intent(in) :: x(:)
      real(dp) :: s
      real(dp) :: xbar
      integer :: n

      n = size(x)
      if (n < 2) then
         s = 0.0_dp
         return
      end if

      xbar = mean_value(x)
      s = sqrt(sum((x - xbar)**2) / real(n - 1, dp))
   end function sample_sd


   function quantile_type7(x, probability) result(q)
      real(dp), intent(in) :: x(:)
      real(dp), intent(in) :: probability
      real(dp) :: q

      real(dp), allocatable :: work(:)

      allocate(work(size(x)))
      work = x
      call sort_real(work)
      q = quantile_type7_sorted(work, probability)
   end function quantile_type7


   pure function quantile_type7_sorted(sorted_x, probability) result(q)
      real(dp), intent(in) :: sorted_x(:)
      real(dp), intent(in) :: probability
      real(dp) :: q

      real(dp) :: gamma, h
      integer :: j, n

      n = size(sorted_x)
      if (n == 0) then
         q = 0.0_dp
      else if (probability <= 0.0_dp) then
         q = sorted_x(1)
      else if (probability >= 1.0_dp) then
         q = sorted_x(n)
      else
         h = 1.0_dp + real(n - 1, dp) * probability
         j = int(floor(h))
         gamma = h - real(j, dp)
         if (j >= n) then
            q = sorted_x(n)
         else
            q = (1.0_dp - gamma) * sorted_x(j) + gamma * sorted_x(j + 1)
         end if
      end if
   end function quantile_type7_sorted


   function descending_order(values) result(order)
      real(dp), intent(in) :: values(:)
      integer, allocatable :: order(:)

      integer, allocatable :: work(:)
      integer :: i, n

      n = size(values)
      allocate(order(n), work(n))
      order = [(i, i = 1, n)]
      if (n > 1) call merge_sort_indices(values, order, work, 1, n)
   end function descending_order


   recursive subroutine merge_sort_indices(values, order, work, left, right)
      real(dp), intent(in) :: values(:)
      integer, intent(inout) :: order(:)
      integer, intent(inout) :: work(:)
      integer, intent(in) :: left, right

      integer :: i, j, k, middle

      if (left >= right) return

      middle = left + (right - left) / 2
      call merge_sort_indices(values, order, work, left, middle)
      call merge_sort_indices(values, order, work, middle + 1, right)

      i = left
      j = middle + 1
      do k = left, right
         if (i > middle) then
            work(k) = order(j)
            j = j + 1
         else if (j > right) then
            work(k) = order(i)
            i = i + 1
         else if (values(order(i)) >= values(order(j))) then
            ! Prefer the left item on ties, preserving original index order.
            work(k) = order(i)
            i = i + 1
         else
            work(k) = order(j)
            j = j + 1
         end if
      end do
      order(left:right) = work(left:right)
   end subroutine merge_sort_indices


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

end module qcsis_statistics
