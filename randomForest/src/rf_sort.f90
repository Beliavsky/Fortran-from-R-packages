! SPDX-License-Identifier: GPL-2.0-or-later
module rf_sort
   use r_kinds, only : dp
   implicit none
   private
   public :: sort_indices_by_values, sort_categories_by_values

contains

   subroutine sort_indices_by_values(values, indices)
      real(dp), intent(in) :: values(:)
      integer, intent(inout) :: indices(:)

      if (size(indices) > 1) call quicksort_indices(values, indices, 1, size(indices))
   end subroutine sort_indices_by_values

   recursive subroutine quicksort_indices(values, indices, lo, hi)
      real(dp), intent(in) :: values(:)
      integer, intent(inout) :: indices(:)
      integer, intent(in) :: lo, hi
      integer :: i, j, tmp, pivot_index
      real(dp) :: pivot

      if (lo >= hi) return
      pivot_index = indices((lo + hi) / 2)
      pivot = values(pivot_index)
      i = lo
      j = hi
      do
         do while (values(indices(i)) < pivot)
            i = i + 1
         end do
         do while (values(indices(j)) > pivot)
            j = j - 1
         end do
         if (i > j) exit
         tmp = indices(i)
         indices(i) = indices(j)
         indices(j) = tmp
         i = i + 1
         j = j - 1
         if (i > hi .or. j < lo) exit
      end do
      if (lo < j) call quicksort_indices(values, indices, lo, j)
      if (i < hi) call quicksort_indices(values, indices, i, hi)
   end subroutine quicksort_indices

   subroutine sort_categories_by_values(values, categories, n)
      real(dp), intent(in) :: values(:)
      integer, intent(out) :: categories(:)
      integer, intent(in) :: n
      integer :: i, j, key
      real(dp) :: key_value

      categories(1:n) = [(i, i = 1, n)]
      do i = 2, n
         key = categories(i)
         key_value = values(key)
         j = i - 1
         do while (j >= 1)
            if (values(categories(j)) <= key_value) exit
            categories(j + 1) = categories(j)
            j = j - 1
         end do
         categories(j + 1) = key
      end do
   end subroutine sort_categories_by_values

end module rf_sort
