! SPDX-License-Identifier: MIT
! SPDX-FileComment: Shared key utilities for the Fortran tidyr translation.
module tidyr_support
   use vctrs, only: vctrs_character, vctrs_data_frame, vctrs_integer, vctrs_logical, vctrs_real
   implicit none
   private

   public :: column_indices, row_keys_equal

contains

   pure function column_indices(data, names) result(indices)
      !! Resolves exact, unique column names to one-based positions.
      type(vctrs_data_frame), intent(in) :: data !! Data frame to search.
      character(len=*), intent(in) :: names(:)   !! Requested exact names.
      integer, allocatable :: indices(:)
      integer :: i, j

      allocate (indices(size(names)))
      do j = 1, size(names)
         indices(j) = data%column_index(names(j))
         if (indices(j) == 0) error stop "tidyr: column not found: " // trim(names(j))
         do i = 1, j - 1
            if (indices(i) == indices(j)) error stop "tidyr: a column was selected more than once"
         end do
      end do
   end function column_indices

   pure logical function row_keys_equal(data, columns, left, right) result(equal)
      !! Compares two rows across selected columns, treating paired missing values as equal.
      type(vctrs_data_frame), intent(in) :: data !! Data frame containing the rows.
      integer, intent(in) :: columns(:)          !! Key-column positions.
      integer, intent(in) :: left                !! First row position.
      integer, intent(in) :: right               !! Second row position.
      integer :: column, j

      equal = .true.
      do j = 1, size(columns)
         column = columns(j)
         if (data%columns(column)%missing(left) .neqv. data%columns(column)%missing(right)) then
            equal = .false.
            return
         end if
         if (data%columns(column)%missing(left)) cycle
         select case (data%columns(column)%type_code)
         case (vctrs_integer)
            equal = data%columns(column)%integer_values(left) == &
               data%columns(column)%integer_values(right)
         case (vctrs_real)
            equal = data%columns(column)%real_values(left) <= &
               data%columns(column)%real_values(right) .and. &
               data%columns(column)%real_values(right) <= data%columns(column)%real_values(left)
         case (vctrs_logical)
            equal = data%columns(column)%logical_values(left) .eqv. &
               data%columns(column)%logical_values(right)
         case (vctrs_character)
            equal = data%columns(column)%character_values(left) == &
               data%columns(column)%character_values(right)
         case default
            error stop "row_keys_equal: unsupported vector type"
         end select
         if (.not. equal) return
      end do
   end function row_keys_equal

end module tidyr_support
