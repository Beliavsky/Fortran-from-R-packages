! SPDX-License-Identifier: MIT
module dplyr_sets
   use dplyr_joins, only: anti_join, semi_join
   use dplyr_verbs, only: bind_rows, distinct
   use tibble, only: tibble_type
   implicit none
   private

   public :: intersect, setdiff, setequal, symdiff, union, union_all

contains

   function union(x, y) result(out)
      !! Returns distinct rows occurring in either schema-identical table.
      type(tibble_type), intent(in) :: x !! First input table.
      type(tibble_type), intent(in) :: y !! Second input table.
      type(tibble_type) :: out
      out = distinct(bind_rows([x, y]))
   end function union

   function union_all(x, y) result(out)
      !! Concatenates all rows of two schema-identical tables.
      type(tibble_type), intent(in) :: x !! First input table.
      type(tibble_type), intent(in) :: y !! Second input table.
      type(tibble_type) :: out
      out = bind_rows([x, y])
   end function union_all

   function intersect(x, y) result(out)
      !! Returns distinct rows occurring in both schema-identical tables.
      type(tibble_type), intent(in) :: x !! First input table.
      type(tibble_type), intent(in) :: y !! Second input table.
      type(tibble_type) :: out
      character(len=:), allocatable :: names(:)

      names = all_column_names(x, y)
      out = distinct(semi_join(x, y, names))
   end function intersect

   function setdiff(x, y) result(out)
      !! Returns distinct rows of `x` not occurring in `y`.
      type(tibble_type), intent(in) :: x !! First input table.
      type(tibble_type), intent(in) :: y !! Second input table.
      type(tibble_type) :: out
      character(len=:), allocatable :: names(:)

      names = all_column_names(x, y)
      out = distinct(anti_join(x, y, names))
   end function setdiff

   logical function setequal(x, y) result(equal)
      !! Reports whether two tables contain the same distinct rows.
      type(tibble_type), intent(in) :: x !! First input table.
      type(tibble_type), intent(in) :: y !! Second input table.
      type(tibble_type) :: left_only, right_only

      left_only = setdiff(x, y)
      right_only = setdiff(y, x)
      equal = left_only%nrow() == 0 .and. right_only%nrow() == 0
   end function setequal

   function symdiff(x, y) result(out)
      !! Returns distinct rows occurring in exactly one input table.
      type(tibble_type), intent(in) :: x !! First input table.
      type(tibble_type), intent(in) :: y !! Second input table.
      type(tibble_type) :: out
      out = union_all(setdiff(x, y), setdiff(y, x))
   end function symdiff

   function all_column_names(x, y) result(names)
      !! Returns common ordered names after validating identical schemas.
      type(tibble_type), intent(in) :: x !! First table to validate.
      type(tibble_type), intent(in) :: y !! Second table to validate.
      character(len=:), allocatable :: names(:)
      integer :: j, width

      if (x%ncol() /= y%ncol()) error stop "set operation: column counts differ"
      width = 1
      do j = 1, x%ncol()
         if (x%columns(j)%name /= y%columns(j)%name) then
            error stop "set operation: column names or order differ"
         end if
         if (x%columns(j)%type_code /= y%columns(j)%type_code) then
            error stop "set operation: column types differ"
         end if
         width = max(width, len(x%columns(j)%name))
      end do
      allocate (character(len=width) :: names(x%ncol()))
      do j = 1, x%ncol()
         names(j) = x%columns(j)%name
      end do
   end function all_column_names

end module dplyr_sets
