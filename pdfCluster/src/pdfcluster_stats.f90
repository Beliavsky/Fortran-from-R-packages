module pdfcluster_stats
   use pdfcluster_kinds, only : dp
   implicit none
   private

   public :: sample_sd_columns, quantile_type7, sort_real, top_two_indices

contains

   pure subroutine sample_sd_columns(x, sd)
      real(dp), intent(in) :: x(:, :) !! Data matrix with observations in rows and variables in columns.
      real(dp), intent(out) :: sd(:) !! Sample standard deviations, one value for each column of x.
      real(dp) :: mu
      integer :: j
      integer :: n

      n = size(x, 1)
      if (size(sd) /= size(x, 2)) return
      if (n <= 1) then
         sd = 0.0_dp
         return
      end if
      do j = 1, size(x, 2)
         mu = sum(x(:, j))/real(n, dp)
         sd(j) = sqrt(sum((x(:, j) - mu)**2)/real(n - 1, dp))
      end do
   end subroutine sample_sd_columns

   pure subroutine sort_real(x, sorted)
      real(dp), intent(in) :: x(:) !! Values to sort in nondecreasing order.
      real(dp), intent(out) :: sorted(:) !! Sorted copy of x with the same extent.
      real(dp) :: key
      integer :: i
      integer :: j

      if (size(sorted) /= size(x)) return
      sorted = x
      do i = 2, size(sorted)
         key = sorted(i)
         j = i - 1
         do while (j >= 1)
            if (sorted(j) <= key) exit
            sorted(j + 1) = sorted(j)
            j = j - 1
         end do
         sorted(j + 1) = key
      end do
   end subroutine sort_real

   pure real(dp) function quantile_type7(x, prob) result(q)
      real(dp), intent(in) :: x(:) !! Sample whose R type-7 quantile is requested.
      real(dp), intent(in) :: prob !! Probability in the closed interval [0,1].
      real(dp), allocatable :: y(:)
      real(dp) :: h
      real(dp) :: gamma
      integer :: j
      integer :: n

      n = size(x)
      if (n == 0) then
         q = 0.0_dp
         return
      end if
      allocate(y(n))
      call sort_real(x, y)
      if (prob <= 0.0_dp) then
         q = y(1)
         return
      end if
      if (prob >= 1.0_dp) then
         q = y(n)
         return
      end if
      h = 1.0_dp + real(n - 1, dp)*prob
      j = int(floor(h))
      gamma = h - real(j, dp)
      if (j >= n) then
         q = y(n)
      else
         q = (1.0_dp - gamma)*y(j) + gamma*y(j + 1)
      end if
   end function quantile_type7

   pure subroutine top_two_indices(x, first, second)
      real(dp), intent(in) :: x(:) !! Score vector for which the two largest entries are requested.
      integer, intent(out) :: first !! One-based index of the largest score, with first-index tie breaking.
      integer, intent(out) :: second !! One-based index of the second-largest score.
      integer :: i

      first = 0
      second = 0
      if (size(x) == 0) return
      first = 1
      do i = 2, size(x)
         if (x(i) > x(first)) first = i
      end do
      if (size(x) == 1) then
         second = first
         return
      end if
      if (first == 1) then
         second = 2
      else
         second = 1
      end if
      do i = 1, size(x)
         if (i == first) cycle
         if (x(i) > x(second)) second = i
      end do
   end subroutine top_two_indices

end module pdfcluster_stats
