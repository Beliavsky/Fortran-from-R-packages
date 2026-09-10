! SPDX-License-Identifier: MIT
module dplyr_groups
   use tibble, only: new_tibble, slice_rows, tibble_column, tibble_type
   use vctrs, only: dp, new_vctr, vctr_type, vec_init, &
      vctrs_character, vctrs_integer, vctrs_logical, vctrs_real
   implicit none
   private

   type, public :: grouped_df
      !! Stores a table and stable first-occurrence grouping metadata.
      type(tibble_type) :: data
      integer, allocatable :: keys(:)
      integer, allocatable :: group_id(:)
      integer, allocatable :: first_rows(:)
      integer, allocatable :: sizes(:)
   end type grouped_df

   type, public :: summary_spec
      !! Describes one built-in grouped aggregation.
      character(len=:), allocatable :: name
      character(len=:), allocatable :: column
      character(len=:), allocatable :: operation
      logical :: na_rm = .false.
   end type summary_spec

   public :: count_rows, group_by, group_indices, group_keys, group_size, n_groups
   public :: summarise, summary, ungroup

contains

   function group_by(data, names) result(grouped)
      !! Groups rows by exact columns in stable first-occurrence order.
      type(tibble_type), intent(in) :: data !! Input table.
      character(len=*), intent(in) :: names(:) !! Exact grouping columns.
      type(grouped_df) :: grouped
      integer, allocatable :: first_buffer(:), size_buffer(:)
      integer :: group, i, j, number_of_groups

      if (size(names) == 0) error stop "group_by: at least one grouping column is required"
      grouped%data = data
      allocate (grouped%keys(size(names)), grouped%group_id(data%nrow()))
      do j = 1, size(names)
         grouped%keys(j) = data%column_index(names(j))
         if (grouped%keys(j) == 0) error stop "group_by: column not found: " // trim(names(j))
         if (j > 1) then
            if (any(grouped%keys(:j - 1) == grouped%keys(j))) then
               error stop "group_by: duplicate grouping column"
            end if
         end if
      end do
      allocate (first_buffer(max(1, data%nrow())), size_buffer(max(1, data%nrow())), source=0)
      number_of_groups = 0
      do i = 1, data%nrow()
         group = 0
         do j = 1, number_of_groups
            if (rows_equal(data, grouped%keys, i, first_buffer(j))) then
               group = j
               exit
            end if
         end do
         if (group == 0) then
            number_of_groups = number_of_groups + 1
            group = number_of_groups
            first_buffer(group) = i
         end if
         grouped%group_id(i) = group
         size_buffer(group) = size_buffer(group) + 1
      end do
      allocate (grouped%first_rows(number_of_groups), grouped%sizes(number_of_groups))
      if (number_of_groups > 0) then
         grouped%first_rows = first_buffer(:number_of_groups)
         grouped%sizes = size_buffer(:number_of_groups)
      end if
   end function group_by

   pure elemental function ungroup(grouped) result(data)
      !! Returns the underlying ungrouped table.
      type(grouped_df), intent(in) :: grouped !! Grouped table.
      type(tibble_type) :: data
      data = grouped%data
   end function ungroup

   pure elemental integer function n_groups(grouped) result(number)
      !! Returns the number of distinct groups.
      type(grouped_df), intent(in) :: grouped !! Grouped table.
      number = size(grouped%first_rows)
   end function n_groups

   pure function group_indices(grouped) result(indices)
      !! Returns the one-based group identifier for each input row.
      type(grouped_df), intent(in) :: grouped !! Grouped table.
      integer, allocatable :: indices(:)
      indices = grouped%group_id
   end function group_indices

   pure function group_size(grouped) result(sizes)
      !! Returns the row count of each group.
      type(grouped_df), intent(in) :: grouped !! Grouped table.
      integer, allocatable :: sizes(:)
      sizes = grouped%sizes
   end function group_size

   function group_keys(grouped) result(keys)
      !! Returns one row of grouping-key values per group.
      type(grouped_df), intent(in) :: grouped !! Grouped table.
      type(tibble_type) :: keys
      type(tibble_column), allocatable :: columns(:)
      integer :: j

      allocate (columns(size(grouped%keys)))
      do j = 1, size(grouped%keys)
         columns(j) = slice_vector(grouped%data%columns(grouped%keys(j)), grouped%first_rows)
      end do
      keys = new_tibble(columns, nrow=n_groups(grouped))
   end function group_keys

   pure elemental function summary(name, column, operation, na_rm) result(specification)
      !! Constructs a built-in grouped summary specification.
      character(len=*), intent(in) :: name      !! Output column name.
      character(len=*), intent(in) :: column    !! Input column; empty for operation `n`.
      character(len=*), intent(in) :: operation !! n, sum, mean, min, or max.
      logical, intent(in), optional :: na_rm    !! Whether missing values are removed.
      type(summary_spec) :: specification

      if (len_trim(name) == 0) error stop "summary: output name cannot be empty"
      specification%name = trim(name)
      specification%column = trim(column)
      specification%operation = trim(operation)
      if (present(na_rm)) specification%na_rm = na_rm
   end function summary

   function summarise(grouped, specifications) result(out)
      !! Reduces each group using built-in summary specifications.
      type(grouped_df), intent(in) :: grouped !! Grouped table.
      type(summary_spec), intent(in) :: specifications(:) !! Requested aggregations.
      type(tibble_type) :: out, keys
      type(tibble_column), allocatable :: columns(:)
      integer :: j

      keys = group_keys(grouped)
      allocate (columns(keys%ncol() + size(specifications)))
      if (keys%ncol() > 0) columns(:keys%ncol()) = keys%columns
      do j = 1, size(specifications)
         columns(keys%ncol() + j) = evaluate_summary(grouped, specifications(j))
      end do
      out = new_tibble(columns, nrow=n_groups(grouped))
   end function summarise

   function count_rows(data, names, name) result(out)
      !! Counts rows for each distinct combination of grouping columns.
      type(tibble_type), intent(in) :: data !! Input table.
      character(len=*), intent(in) :: names(:) !! Exact grouping columns.
      character(len=*), intent(in), optional :: name !! Count-column name; default `n`.
      type(tibble_type) :: out
      type(grouped_df) :: grouped
      character(len=:), allocatable :: count_name

      count_name = "n"
      if (present(name)) count_name = trim(name)
      grouped = group_by(data, names)
      out = summarise(grouped, [summary(count_name, "", "n")])
   end function count_rows

   function evaluate_summary(grouped, specification) result(vector)
      !! Evaluates one specification for every group.
      type(grouped_df), intent(in) :: grouped !! Grouped table.
      type(summary_spec), intent(in) :: specification !! Aggregation request.
      type(vctr_type) :: vector, source, prototype
      integer :: column, group

      if (specification%operation == "n") then
         vector = new_vctr(specification%name, grouped%sizes)
         return
      end if
      column = grouped%data%column_index(specification%column)
      if (column == 0) error stop "summarise: column not found: " // specification%column
      source = grouped%data%columns(column)
      if (source%type_code /= vctrs_integer .and. source%type_code /= vctrs_real) then
         error stop "summarise: built-in reductions require an integer or real column"
      end if
      if (specification%operation == "mean") then
         prototype = new_vctr(specification%name, 0.0_dp, 1)
      else
         prototype = source
         prototype%name = specification%name
      end if
      vector = vec_init(prototype, n_groups(grouped), missing=.false.)
      do group = 1, n_groups(grouped)
         call reduce_group(source, grouped%group_id, group, specification, vector)
      end do
   end function evaluate_summary

   subroutine reduce_group(source, group_id, group, specification, result)
      !! Reduces one group's nonmissing numeric elements into a result position.
      type(vctr_type), intent(in) :: source !! Numeric source vector.
      integer, intent(in) :: group_id(:)    !! Group identifier per row.
      integer, intent(in) :: group          !! Group being reduced.
      type(summary_spec), intent(in) :: specification !! Aggregation request.
      type(vctr_type), intent(inout) :: result !! Result vector to populate.
      integer :: count_values, i, integer_value
      real(dp) :: real_value, total
      logical :: has_missing

      count_values = 0
      has_missing = .false.
      total = 0.0_dp
      integer_value = 0
      real_value = 0.0_dp
      do i = 1, source%size()
         if (group_id(i) /= group) cycle
         if (source%missing(i)) then
            has_missing = .true.
            cycle
         end if
         count_values = count_values + 1
         if (source%type_code == vctrs_integer) then
            call update_integer(source%integer_values(i), count_values, specification%operation, &
               integer_value, total)
         else
            call update_real(source%real_values(i), count_values, specification%operation, &
               real_value, total)
         end if
      end do
      if (has_missing .and. .not. specification%na_rm) then
         result%missing(group) = .true.
      else if (count_values == 0 .and. specification%operation /= "sum") then
         result%missing(group) = .true.
      else if (specification%operation == "mean") then
         result%real_values(group) = total / real(count_values, dp)
      else if (source%type_code == vctrs_integer) then
         if (specification%operation == "sum") integer_value = int(total)
         result%integer_values(group) = integer_value
      else
         if (specification%operation == "sum") real_value = total
         result%real_values(group) = real_value
      end if
   end subroutine reduce_group

   pure elemental subroutine update_integer(value, count_values, operation, aggregate, total)
      !! Incorporates one integer into a built-in reduction.
      integer, intent(in) :: value        !! Value to incorporate.
      integer, intent(in) :: count_values !! Count including this value.
      character(len=*), intent(in) :: operation !! Reduction name.
      integer, intent(inout) :: aggregate !! Minimum or maximum accumulator.
      real(dp), intent(inout) :: total    !! Sum accumulator.

      select case (operation)
      case ("sum", "mean")
         total = total + real(value, dp)
      case ("min")
         if (count_values == 1) aggregate = value
         if (count_values > 1) aggregate = min(aggregate, value)
      case ("max")
         if (count_values == 1) aggregate = value
         if (count_values > 1) aggregate = max(aggregate, value)
      case default
         error stop "summarise: operation must be n, sum, mean, min, or max"
      end select
   end subroutine update_integer

   pure elemental subroutine update_real(value, count_values, operation, aggregate, total)
      !! Incorporates one real value into a built-in reduction.
      real(dp), intent(in) :: value       !! Value to incorporate.
      integer, intent(in) :: count_values !! Count including this value.
      character(len=*), intent(in) :: operation !! Reduction name.
      real(dp), intent(inout) :: aggregate !! Minimum or maximum accumulator.
      real(dp), intent(inout) :: total     !! Sum accumulator.

      select case (operation)
      case ("sum", "mean")
         total = total + value
      case ("min")
         if (count_values == 1) aggregate = value
         if (count_values > 1) aggregate = min(aggregate, value)
      case ("max")
         if (count_values == 1) aggregate = value
         if (count_values > 1) aggregate = max(aggregate, value)
      case default
         error stop "summarise: operation must be n, sum, mean, min, or max"
      end select
   end subroutine update_real

   pure function slice_vector(vector, rows) result(out)
      !! Slices a vector without introducing a module dependency cycle.
      use vctrs, only: vec_slice
      type(vctr_type), intent(in) :: vector !! Source vector.
      integer, intent(in) :: rows(:)        !! One-based positions.
      type(vctr_type) :: out
      out = vec_slice(vector, rows)
   end function slice_vector

   pure logical function rows_equal(data, keys, left, right) result(equal)
      !! Compares two rows across grouping keys, treating paired missing values as equal.
      type(tibble_type), intent(in) :: data !! Table containing keys.
      integer, intent(in) :: keys(:)        !! Key-column positions.
      integer, intent(in) :: left  !! First row position.
      integer, intent(in) :: right !! Second row position.
      integer :: j

      equal = .true.
      do j = 1, size(keys)
         if (.not. elements_equal(data%columns(keys(j)), left, right)) then
            equal = .false.
            return
         end if
      end do
   end function rows_equal

   pure elemental logical function elements_equal(vector, left, right) result(equal)
      !! Compares two elements for grouping equality.
      type(vctr_type), intent(in) :: vector !! Vector containing values.
      integer, intent(in) :: left  !! First element position.
      integer, intent(in) :: right !! Second element position.

      if (vector%missing(left) .neqv. vector%missing(right)) then
         equal = .false.
      else if (vector%missing(left)) then
         equal = .true.
      else
         select case (vector%type_code)
         case (vctrs_integer)
            equal = vector%integer_values(left) == vector%integer_values(right)
         case (vctrs_real)
            equal = vector%real_values(left) <= vector%real_values(right) .and. &
               vector%real_values(right) <= vector%real_values(left)
         case (vctrs_logical)
            equal = vector%logical_values(left) .eqv. vector%logical_values(right)
         case (vctrs_character)
            equal = vector%character_values(left) == vector%character_values(right)
         case default
            error stop "group_by: unsupported key type"
         end select
      end if
   end function elements_equal

end module dplyr_groups
