module psych_interpolation
   use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
   use r_kinds, only : dp
   implicit none
   private
   public :: interpolated_quantile, interpolated_quartiles, interpolated_values
contains
   pure real(dp) function interpolated_quantile(x, q, width) result(value)
      real(dp), intent(in) :: x(:) !! Sample values; IEEE NaNs are omitted.
      real(dp), intent(in) :: q !! Requested quantile strictly between zero and one.
      real(dp), intent(in), optional :: width !! Category width for tied-value interpolation; defaults to one.
      real(dp), allocatable :: work(:)
      real(dp) :: w, n2, midpoint
      integer :: i, j, n, low_index, high_index, below, above, at_mid
      n = count(.not. ieee_is_nan(x))
      if (n == 0) then
         value = 0.0_dp
         return
      end if
      allocate(work(n))
      j = 0
      do i = 1, size(x)
         if (.not. ieee_is_nan(x(i))) then
            j = j + 1
            work(j) = x(i)
         end if
      end do
      call sort_real(work)
      w = 1.0_dp
      if (present(width)) w = width
      n2 = real(n + 1, dp) * q
      low_index = max(1, min(n, int(floor(n2))))
      high_index = max(1, min(n, int(ceiling(n2))))
      midpoint = 0.5_dp * (work(low_index) + work(high_index))
      below = count(work < midpoint)
      above = count(work > midpoint)
      at_mid = n - below - above
      if (at_mid > 1) then
         value = midpoint - 0.5_dp * w + w * (real(n, dp) * q - real(below, dp)) / real(at_mid, dp)
      else
         value = midpoint
      end if
   end function interpolated_quantile

   pure subroutine interpolated_quartiles(x, quartiles, width)
      real(dp), intent(in) :: x(:) !! Sample values; IEEE NaNs are omitted.
      real(dp), intent(out) :: quartiles(3) !! Interpolated first quartile, median, and third quartile.
      real(dp), intent(in), optional :: width !! Category width for tied-value interpolation; defaults to one.
      quartiles(1) = interpolated_quantile(x, 0.25_dp, width)
      quartiles(2) = interpolated_quantile(x, 0.50_dp, width)
      quartiles(3) = interpolated_quantile(x, 0.75_dp, width)
   end subroutine interpolated_quartiles

   pure subroutine interpolated_values(x, values, width)
      real(dp), intent(in) :: x(:) !! Sample values to replace by within-category evenly spaced ranks.
      real(dp), allocatable, intent(out) :: values(:) !! Sorted interpolated within-category positions.
      real(dp), intent(in), optional :: width !! Category width for interpolation; defaults to one.
      real(dp), allocatable :: work(:)
      real(dp) :: w
      integer :: i, j, k, first, last, count_same
      w = 1.0_dp
      if (present(width)) w = width
      allocate(work(size(x)), values(size(x)))
      work = x
      call sort_real(work)
      first = 1
      k = 0
      do while (first <= size(work))
         last = first
         do while (last < size(work))
            if (abs(work(last + 1) - work(first)) > epsilon(1.0_dp) * &
               max(1.0_dp, abs(work(last + 1)), abs(work(first)))) exit
            last = last + 1
         end do
         count_same = last - first + 1
         do j = 1, count_same
            k = k + 1
            values(k) = real(first, dp) - 1.5_dp * w + w * real(j, dp) / real(count_same + 1, dp)
         end do
         first = last + 1
      end do
      do i = k + 1, size(values)
         values(i) = 0.0_dp
      end do
   end subroutine interpolated_values

   pure subroutine sort_real(x)
      real(dp), intent(inout) :: x(:) !! Real vector sorted in ascending order in place.
      real(dp) :: key
      integer :: i, j
      do i = 2, size(x)
         key = x(i)
         j = i - 1
         do while (j >= 1)
            if (x(j) <= key) exit
            x(j + 1) = x(j)
            j = j - 1
         end do
         x(j + 1) = key
      end do
   end subroutine sort_real
end module psych_interpolation
