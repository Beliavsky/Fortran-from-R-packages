! SPDX-License-Identifier: MIT
module dplyr_joins
   use tibble, only: filter_rows, new_tibble, tibble_column, tibble_type
   use vctrs, only: dp, vctr_type, vec_cast, vec_copy_element, vec_init, vec_ptype2, &
      vctrs_character, vctrs_integer, vctrs_logical, vctrs_real
   implicit none
   private

   public :: anti_join, cross_join, full_join, inner_join, left_join, right_join, semi_join

contains

   function inner_join(x, y, by) result(out)
      !! Returns all matching equality-join row pairs.
      type(tibble_type), intent(in) :: x !! Left table.
      type(tibble_type), intent(in) :: y !! Right table.
      character(len=*), intent(in) :: by(:) !! Same-name equality keys.
      type(tibble_type) :: out
      out = mutating_join(x, y, by, "inner")
   end function inner_join

   function left_join(x, y, by) result(out)
      !! Returns every left row and all matching right rows.
      type(tibble_type), intent(in) :: x !! Left table.
      type(tibble_type), intent(in) :: y !! Right table.
      character(len=*), intent(in) :: by(:) !! Same-name equality keys.
      type(tibble_type) :: out
      out = mutating_join(x, y, by, "left")
   end function left_join

   function right_join(x, y, by) result(out)
      !! Returns all matched pairs followed by unmatched right rows.
      type(tibble_type), intent(in) :: x !! Left table.
      type(tibble_type), intent(in) :: y !! Right table.
      character(len=*), intent(in) :: by(:) !! Same-name equality keys.
      type(tibble_type) :: out
      out = mutating_join(x, y, by, "right")
   end function right_join

   function full_join(x, y, by) result(out)
      !! Returns matched pairs and unmatched rows from both inputs.
      type(tibble_type), intent(in) :: x !! Left table.
      type(tibble_type), intent(in) :: y !! Right table.
      character(len=*), intent(in) :: by(:) !! Same-name equality keys.
      type(tibble_type) :: out
      out = mutating_join(x, y, by, "full")
   end function full_join

   function mutating_join(x, y, by, kind) result(out)
      !! Implements same-name equality joins with stable nested-loop matching.
      type(tibble_type), intent(in) :: x !! Left table.
      type(tibble_type), intent(in) :: y !! Right table.
      character(len=*), intent(in) :: by(:) !! Same-name equality keys.
      character(len=*), intent(in) :: kind  !! inner, left, right, or full.
      type(tibble_type) :: out
      type(tibble_column), allocatable :: columns(:)
      integer, allocatable :: x_keys(:), x_rows(:), y_keys(:), y_rows(:)
      logical, allocatable :: matched_y(:)
      integer :: capacity, i, j, number_of_pairs
      logical :: matched

      call resolve_join_keys(x, y, by, x_keys, y_keys)
      if (x%nrow() > 0 .and. y%nrow() > 0) then
         if (x%nrow() > (huge(capacity) - x%nrow() - y%nrow()) / y%nrow()) then
            error stop "join: candidate-pair count overflows default integer"
         end if
      end if
      capacity = x%nrow() * y%nrow() + x%nrow() + y%nrow()
      allocate (x_rows(capacity), y_rows(capacity))
      allocate (matched_y(y%nrow()), source=.false.)
      number_of_pairs = 0
      do i = 1, x%nrow()
         matched = .false.
         do j = 1, y%nrow()
            if (join_rows_equal(x, y, x_keys, y_keys, i, j)) then
               matched = .true.
               matched_y(j) = .true.
               number_of_pairs = number_of_pairs + 1
               x_rows(number_of_pairs) = i
               y_rows(number_of_pairs) = j
            end if
         end do
         if (.not. matched .and. (kind == "left" .or. kind == "full")) then
            number_of_pairs = number_of_pairs + 1
            x_rows(number_of_pairs) = i
            y_rows(number_of_pairs) = 0
         end if
      end do
      if (kind == "right" .or. kind == "full") then
         do j = 1, y%nrow()
            if (.not. matched_y(j)) then
               number_of_pairs = number_of_pairs + 1
               x_rows(number_of_pairs) = 0
               y_rows(number_of_pairs) = j
            end if
         end do
      end if
      call materialize_join(x, y, x_keys, y_keys, x_rows(:number_of_pairs), &
         y_rows(:number_of_pairs), columns)
      out = new_tibble(columns, nrow=number_of_pairs)
   end function mutating_join

   subroutine materialize_join(x, y, x_keys, y_keys, x_rows, y_rows, columns)
      !! Builds output columns from row-pair indices and coalesced common-type keys.
      type(tibble_type), intent(in) :: x !! Left table.
      type(tibble_type), intent(in) :: y !! Right table.
      integer, intent(in) :: x_keys(:) !! Left key positions.
      integer, intent(in) :: y_keys(:) !! Right key positions.
      integer, intent(in) :: x_rows(:) !! Left row positions; zero means absent.
      integer, intent(in) :: y_rows(:) !! Right row positions; zero means absent.
      type(tibble_column), allocatable, intent(out) :: columns(:) !! Materialized columns.
      type(vctr_type) :: left_source, right_source
      logical, allocatable :: y_is_key(:)
      integer :: common_type, i, j, key_position, output_column

      allocate (y_is_key(y%ncol()), source=.false.)
      y_is_key(y_keys) = .true.
      allocate (columns(x%ncol() + count(.not. y_is_key)))
      do j = 1, x%ncol()
         key_position = find_position(x_keys, j)
         if (key_position > 0) then
            common_type = vec_ptype2(x%columns(j), y%columns(y_keys(key_position)))
            if (common_type == 0) error stop "join: key columns have incompatible types"
            left_source = vec_cast(x%columns(j), common_type)
            right_source = vec_cast(y%columns(y_keys(key_position)), common_type)
            columns(j) = vec_init(left_source, size(x_rows), missing=.true.)
            do i = 1, size(x_rows)
               if (x_rows(i) > 0) then
                  call vec_copy_element(left_source, x_rows(i), columns(j), i)
               else
                  call vec_copy_element(right_source, y_rows(i), columns(j), i)
               end if
            end do
         else
            columns(j) = vec_init(x%columns(j), size(x_rows), missing=.true.)
            do i = 1, size(x_rows)
               if (x_rows(i) > 0) call vec_copy_element(x%columns(j), x_rows(i), columns(j), i)
            end do
         end if
      end do
      output_column = x%ncol()
      do j = 1, y%ncol()
         if (y_is_key(j)) cycle
         output_column = output_column + 1
         columns(output_column) = vec_init(y%columns(j), size(y_rows), missing=.true.)
         if (x%has_name(columns(output_column)%name)) then
            columns(output_column)%name = columns(output_column)%name // ".y"
         end if
         do i = 1, output_column - 1
            if (columns(i)%name == columns(output_column)%name) then
               error stop "join: suffix did not produce a unique output name"
            end if
         end do
         do i = 1, size(y_rows)
            if (y_rows(i) > 0) then
               call vec_copy_element(y%columns(j), y_rows(i), columns(output_column), i)
            end if
         end do
      end do
   end subroutine materialize_join

   function semi_join(x, y, by) result(out)
      !! Keeps left rows having at least one matching right row.
      type(tibble_type), intent(in) :: x !! Left table.
      type(tibble_type), intent(in) :: y !! Right table.
      character(len=*), intent(in) :: by(:) !! Same-name equality keys.
      type(tibble_type) :: out
      logical, allocatable :: keep(:)

      keep = matching_left_rows(x, y, by)
      out = filter_rows(x, keep)
   end function semi_join

   function anti_join(x, y, by) result(out)
      !! Keeps left rows having no matching right row.
      type(tibble_type), intent(in) :: x !! Left table.
      type(tibble_type), intent(in) :: y !! Right table.
      character(len=*), intent(in) :: by(:) !! Same-name equality keys.
      type(tibble_type) :: out
      logical, allocatable :: keep(:)

      keep = .not. matching_left_rows(x, y, by)
      out = filter_rows(x, keep)
   end function anti_join

   function matching_left_rows(x, y, by) result(matched)
      !! Identifies left rows with at least one equal-key right row.
      type(tibble_type), intent(in) :: x !! Left table.
      type(tibble_type), intent(in) :: y !! Right table.
      character(len=*), intent(in) :: by(:) !! Same-name equality keys.
      logical, allocatable :: matched(:)
      integer, allocatable :: x_keys(:), y_keys(:)
      integer :: i, j

      call resolve_join_keys(x, y, by, x_keys, y_keys)
      allocate (matched(x%nrow()), source=.false.)
      do i = 1, x%nrow()
         do j = 1, y%nrow()
            if (join_rows_equal(x, y, x_keys, y_keys, i, j)) then
               matched(i) = .true.
               exit
            end if
         end do
      end do
   end function matching_left_rows

   function cross_join(x, y) result(out)
      !! Returns the Cartesian product of two tables.
      type(tibble_type), intent(in) :: x !! Left table.
      type(tibble_type), intent(in) :: y !! Right table.
      type(tibble_type) :: out
      type(tibble_column), allocatable :: columns(:)
      integer, allocatable :: empty_keys(:), x_rows(:), y_rows(:)
      integer :: i, j, k, number_of_pairs

      if (x%nrow() > 0) then
         if (y%nrow() > huge(number_of_pairs) / x%nrow()) then
            error stop "cross_join: row count overflows default integer"
         end if
      end if
      number_of_pairs = x%nrow() * y%nrow()
      allocate (x_rows(number_of_pairs), y_rows(number_of_pairs), empty_keys(0))
      k = 0
      do i = 1, x%nrow()
         do j = 1, y%nrow()
            k = k + 1
            x_rows(k) = i
            y_rows(k) = j
         end do
      end do
      call materialize_join(x, y, empty_keys, empty_keys, x_rows, y_rows, columns)
      out = new_tibble(columns, nrow=number_of_pairs)
   end function cross_join

   subroutine resolve_join_keys(x, y, by, x_keys, y_keys)
      !! Resolves and validates same-name equality keys in two tables.
      type(tibble_type), intent(in) :: x !! Left table to search.
      type(tibble_type), intent(in) :: y !! Right table to search.
      character(len=*), intent(in) :: by(:) !! Same-name key list.
      integer, allocatable, intent(out) :: x_keys(:) !! Resolved left key positions.
      integer, allocatable, intent(out) :: y_keys(:) !! Resolved right key positions.
      integer :: common_type, i, j

      if (size(by) == 0) error stop "join: at least one key is required"
      allocate (x_keys(size(by)), y_keys(size(by)))
      do j = 1, size(by)
         x_keys(j) = x%column_index(by(j))
         y_keys(j) = y%column_index(by(j))
         if (x_keys(j) == 0 .or. y_keys(j) == 0) error stop "join: key column not found"
         common_type = vec_ptype2(x%columns(x_keys(j)), y%columns(y_keys(j)))
         if (common_type == 0) error stop "join: key columns have incompatible types"
         do i = 1, j - 1
            if (x_keys(i) == x_keys(j)) error stop "join: duplicate key name"
         end do
      end do
   end subroutine resolve_join_keys

   pure logical function join_rows_equal(x, y, x_keys, y_keys, x_row, y_row) result(equal)
      !! Compares two rows across compatible join keys, with paired missing values equal.
      type(tibble_type), intent(in) :: x !! Left table containing keys.
      type(tibble_type), intent(in) :: y !! Right table containing keys.
      integer, intent(in) :: x_keys(:) !! Left key positions.
      integer, intent(in) :: y_keys(:) !! Corresponding right key positions.
      integer, intent(in) :: x_row !! Left row position.
      integer, intent(in) :: y_row !! Right row position.
      integer :: j

      equal = .true.
      do j = 1, size(x_keys)
         if (.not. join_keys_equal( &
             x%columns(x_keys(j)), y%columns(y_keys(j)), x_row, y_row)) then
            equal = .false.
            return
         end if
      end do
   end function join_rows_equal

   pure elemental logical function join_keys_equal(left, right, left_row, right_row) result(equal)
      !! Compares compatible scalar key values using vctrs numeric promotion.
      type(vctr_type), intent(in) :: left  !! Left key vector.
      type(vctr_type), intent(in) :: right !! Right key vector.
      integer, intent(in) :: left_row  !! Left element position.
      integer, intent(in) :: right_row !! Right element position.
      integer :: common_type

      if (left%missing(left_row) .neqv. right%missing(right_row)) then
         equal = .false.
         return
      end if
      if (left%missing(left_row)) then
         equal = .true.
         return
      end if
      common_type = vec_ptype2(left, right)
      select case (common_type)
      case (vctrs_integer)
         equal = integer_key(left, left_row) == integer_key(right, right_row)
      case (vctrs_real)
         equal = real_key(left, left_row) <= real_key(right, right_row) .and. &
            real_key(right, right_row) <= real_key(left, left_row)
      case (vctrs_character)
         equal = left%character_values(left_row) == right%character_values(right_row)
      case default
         error stop "join: unsupported key type"
      end select
   end function join_keys_equal

   pure elemental integer function integer_key(vector, row) result(value)
      !! Converts an integer-compatible key element to integer.
      type(vctr_type), intent(in) :: vector !! Key vector.
      integer, intent(in) :: row            !! Element position.

      if (vector%type_code == vctrs_integer) then
         value = vector%integer_values(row)
      else if (vector%type_code == vctrs_logical) then
         value = merge(1, 0, vector%logical_values(row))
      else
         error stop "join: key is not integer-compatible"
      end if
   end function integer_key

   pure elemental real(dp) function real_key(vector, row) result(value)
      !! Converts a numeric key element to double precision.
      type(vctr_type), intent(in) :: vector !! Key vector.
      integer, intent(in) :: row            !! Element position.

      select case (vector%type_code)
      case (vctrs_real)
         value = vector%real_values(row)
      case (vctrs_integer)
         value = real(vector%integer_values(row), dp)
      case (vctrs_logical)
         value = merge(1.0_dp, 0.0_dp, vector%logical_values(row))
      case default
         error stop "join: key is not numeric"
      end select
   end function real_key

   pure integer function find_position(values, value) result(position)
      !! Finds an integer in an array or returns zero.
      integer, intent(in) :: values(:) !! Values to search.
      integer, intent(in) :: value     !! Value to find.
      integer :: i

      position = 0
      do i = 1, size(values)
         if (values(i) == value) then
            position = i
            return
         end if
      end do
   end function find_position

end module dplyr_joins
