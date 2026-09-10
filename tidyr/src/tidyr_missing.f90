! SPDX-License-Identifier: MIT
! SPDX-FileComment: Missing-value operations for the Fortran tidyr translation.
module tidyr_missing
   !! Implements row removal, scalar replacement, and directional filling.
   use tibble, only: tibble_column, tibble_type, filter_rows, new_tibble
   use vctrs, only: vec_copy_element, vctrs_character, vctrs_integer, vctrs_logical, vctrs_real
   implicit none
   private

   public :: drop_na, fill, replace_na

contains

   function drop_na(data, names) result(output)
      !! Removes rows missing in any selected column.
      type(tibble_type), intent(in) :: data       !! Input table.
      character(len=*), intent(in), optional :: names(:) !! Columns to inspect; all by default.
      type(tibble_type) :: output
      logical, allocatable :: keep(:)
      integer :: i, j

      allocate (keep(data%nrow()), source=.true.)
      if (present(names)) then
         do i = 1, size(names)
            j = data%column_index(names(i))
            if (j == 0) error stop 'drop_na: column not found'
            keep = keep .and. .not. data%columns(j)%missing
         end do
      else
         do j = 1, data%ncol()
            keep = keep .and. .not. data%columns(j)%missing
         end do
      end if
      output = filter_rows(data, keep)
   end function drop_na

   function replace_na(data, replacements) result(output)
      !! Replaces missing values using named size-one vectors.
      type(tibble_type), intent(in) :: data          !! Input table.
      type(tibble_column), intent(in) :: replacements(:) !! Named scalar replacements.
      type(tibble_type) :: output
      integer :: i, j, row

      output = data
      do i = 1, size(replacements)
         j = output%column_index(replacements(i)%name)
         if (j == 0) error stop 'replace_na: column not found'
         if (replacements(i)%size() /= 1) error stop 'replace_na: replacements must have size one'
         if (replacements(i)%type_code /= output%columns(j)%type_code) then
            error stop 'replace_na: replacement type differs from column type'
         end if
         if (replacements(i)%missing(1)) error stop 'replace_na: replacement cannot be missing'
         do row = 1, output%nrow()
            if (output%columns(j)%missing(row)) then
               call vec_copy_element(replacements(i), 1, output%columns(j), row)
            end if
         end do
      end do
   end function replace_na

   function fill(data, names, direction) result(output)
      !! Fills missing values from neighboring observations in selected columns.
      type(tibble_type), intent(in) :: data    !! Input table.
      character(len=*), intent(in) :: names(:) !! Columns to fill.
      character(len=*), intent(in), optional :: direction !! `down`, `up`, `downup`, or `updown`.
      type(tibble_type) :: output
      character(len=:), allocatable :: mode
      integer :: i, j

      mode = 'down'
      if (present(direction)) mode = trim(direction)
      output = data
      do i = 1, size(names)
         j = output%column_index(names(i))
         if (j == 0) error stop 'fill: column not found'
         select case (mode)
         case ('down')
            call fill_down(output%columns(j))
         case ('up')
            call fill_up(output%columns(j))
         case ('downup')
            call fill_down(output%columns(j))
            call fill_up(output%columns(j))
         case ('updown')
            call fill_up(output%columns(j))
            call fill_down(output%columns(j))
         case default
            error stop 'fill: unsupported direction'
         end select
      end do
   end function fill

   pure elemental subroutine fill_down(vector)
      !! Carries the most recent nonmissing value downward.
      type(tibble_column), intent(inout) :: vector !! Column to modify.
      integer :: previous, row
      previous = 0
      do row = 1, vector%size()
         if (.not. vector%missing(row)) then
            previous = row
         else if (previous > 0) then
            call vec_copy_element(vector, previous, vector, row)
         end if
      end do
   end subroutine fill_down

   pure elemental subroutine fill_up(vector)
      !! Carries the next nonmissing value upward.
      type(tibble_column), intent(inout) :: vector !! Column to modify.
      integer :: following, row
      following = 0
      do row = vector%size(), 1, -1
         if (.not. vector%missing(row)) then
            following = row
         else if (following > 0) then
            call vec_copy_element(vector, following, vector, row)
         end if
      end do
   end subroutine fill_up

end module tidyr_missing
