! SPDX-License-Identifier: MIT
! SPDX-FileComment: Table verbs for the Fortran dplyr translation.
module dplyr_verbs
   use dplyr_expression, only: evaluate_expression
   use tibble, only: add_column, filter_rows, new_tibble, replace_column, select_columns, &
      slice_rows, tibble_column, tibble_type
   use vctrs, only: vctr_type, vec_cbind, vec_recycle, vec_rbind, vctrs_character, &
      vctrs_integer, vctrs_logical, vctrs_real
   implicit none
   private

   public :: arrange, bind_cols, bind_rows, distinct, filter, mutate, pull, relocate
   public :: rename, select, slice, slice_head, slice_tail

   interface filter
      module procedure filter_vector, filter_logical, filter_expression
   end interface filter

   interface mutate
      module procedure mutate_one, mutate_many
   end interface mutate

   interface select
      module procedure select_names, select_indices
   end interface select

   interface rename
      module procedure rename_names, rename_indices
   end interface rename

   interface relocate
      module procedure relocate_names, relocate_indices
   end interface relocate

contains

   function filter_vector(data, condition) result(out)
      !! Keeps rows selected by a logical vector, treating missing conditions as false.
      type(tibble_type), intent(in) :: data !! Input table.
      type(vctr_type), intent(in) :: condition !! Logical condition vector.
      type(tibble_type) :: out
      logical, allocatable :: keep(:)
      type(vctr_type) :: recycled

      if (condition%type_code /= vctrs_logical) error stop "filter: condition must be logical"
      recycled = vec_recycle(condition, data%nrow())
      keep = recycled%logical_values .and. .not. recycled%missing
      out = filter_rows(data, keep)
   end function filter_vector

   function filter_logical(data, condition) result(out)
      !! Keeps rows selected by a conformable intrinsic logical array.
      type(tibble_type), intent(in) :: data !! Input table.
      logical, intent(in) :: condition(:)   !! Row-selection mask `(nrow)`.
      type(tibble_type) :: out

      if (size(condition) /= data%nrow()) error stop "filter: condition size differs from row count"
      out = filter_rows(data, condition)
   end function filter_logical

   function filter_expression(data, expression) result(out)
      !! Parses a concise expression and keeps rows where it is true and nonmissing.
      type(tibble_type), intent(in) :: data !! Input table.
      character(len=*), intent(in) :: expression !! Filter expression.
      type(tibble_type) :: out
      type(vctr_type) :: condition

      condition = evaluate_expression(data, expression)
      out = filter_vector(data, condition)
   end function filter_expression

   function select_names(data, names) result(out)
      !! Selects exact columns in requested order.
      type(tibble_type), intent(in) :: data !! Input table.
      character(len=*), intent(in) :: names(:) !! Exact output column names.
      type(tibble_type) :: out

      out = select_columns(data, names)
   end function select_names

   function select_indices(data, indices) result(out)
      !! Selects columns by one-based positions returned by tidyselect.
      type(tibble_type), intent(in) :: data !! Input table.
      integer, intent(in) :: indices(:)     !! Selected positions.
      type(tibble_type) :: out
      if (any(indices < 1) .or. any(indices > data%ncol())) error stop "select: index out of range"
      if (has_duplicates(indices)) error stop "select: duplicate index"
      out = new_tibble(data%columns(indices), nrow=data%nrow())
   end function select_indices

   function rename_names(data, old_names, new_names) result(out)
      !! Renames exact columns and validates uniqueness of the resulting schema.
      type(tibble_type), intent(in) :: data !! Input table.
      character(len=*), intent(in) :: old_names(:) !! Existing exact names.
      character(len=*), intent(in) :: new_names(:) !! Corresponding new names.
      type(tibble_type) :: out
      type(tibble_column), allocatable :: columns(:)
      integer :: index, j

      if (size(old_names) /= size(new_names)) error stop "rename: name-array sizes differ"
      columns = data%columns
      do j = 1, size(old_names)
         index = data%column_index(old_names(j))
         if (index == 0) error stop "rename: column not found: " // trim(old_names(j))
         if (len_trim(new_names(j)) == 0) error stop "rename: new name cannot be empty"
         columns(index)%name = trim(new_names(j))
      end do
      out = new_tibble(columns, nrow=data%nrow())
   end function rename_names

   function rename_indices(data, indices, new_names) result(out)
      !! Renames columns by one-based positions returned by tidyselect.
      type(tibble_type), intent(in) :: data !! Input table.
      integer, intent(in) :: indices(:)     !! Positions to rename.
      character(len=*), intent(in) :: new_names(:) !! Replacement names.
      type(tibble_type) :: out
      type(tibble_column), allocatable :: columns(:)
      integer :: j
      if (size(indices) /= size(new_names)) error stop "rename: size mismatch"
      if (any(indices < 1) .or. any(indices > data%ncol())) error stop "rename: index out of range"
      columns = data%columns
      do j = 1, size(indices)
         columns(indices(j))%name = trim(new_names(j))
      end do
      out = new_tibble(columns, nrow=data%nrow())
   end function rename_indices

   function relocate_names(data, names, before) result(out)
      !! Moves selected columns before a retained column, or to the front by default.
      type(tibble_type), intent(in) :: data !! Input table.
      character(len=*), intent(in) :: names(:) !! Exact columns to move.
      character(len=*), intent(in), optional :: before !! Retained insertion anchor.
      type(tibble_type) :: out
      integer, allocatable :: order(:), selected(:)
      logical, allocatable :: move(:)
      integer :: anchor, i, j, k

      selected = resolve_columns(data, names)
      allocate (move(data%ncol()), source=.false.)
      move(selected) = .true.
      anchor = 1
      if (present(before)) then
         anchor = data%column_index(before)
         if (anchor == 0) error stop "relocate: anchor column not found"
         if (move(anchor)) error stop "relocate: anchor cannot be among moved columns"
      end if
      allocate (order(data%ncol()))
      k = 0
      do j = 1, data%ncol()
         if (j == anchor) then
            do i = 1, size(selected)
               k = k + 1
               order(k) = selected(i)
            end do
         end if
         if (.not. move(j)) then
            k = k + 1
            order(k) = j
         end if
      end do
      out = new_tibble(data%columns(order), nrow=data%nrow())
   end function relocate_names

   function relocate_indices(data, selected, before) result(out)
      !! Relocates columns using one-based selected and anchor positions.
      type(tibble_type), intent(in) :: data !! Input table.
      integer, intent(in) :: selected(:)    !! Positions to move.
      integer, intent(in), optional :: before !! Retained anchor; first by default.
      type(tibble_type) :: out
      integer, allocatable :: order(:)
      logical, allocatable :: move(:)
      integer :: anchor, i, j, k
      if (any(selected < 1) .or. any(selected > data%ncol())) error stop "relocate: index out of range"
      allocate (move(data%ncol()), source=.false.)
      move(selected) = .true.
      anchor = 1
      if (present(before)) anchor = before
      if (anchor < 1 .or. anchor > data%ncol()) error stop "relocate: anchor out of range"
      if (move(anchor)) error stop "relocate: anchor cannot be moved"
      allocate (order(data%ncol()))
      k = 0
      do j = 1, data%ncol()
         if (j == anchor) then
            do i = 1, size(selected)
               k = k + 1
               order(k) = selected(i)
            end do
         end if
         if (.not. move(j)) then
            k = k + 1
            order(k) = j
         end if
      end do
      out = new_tibble(data%columns(order), nrow=data%nrow())
   end function relocate_indices

   function mutate_one(data, new_column) result(out)
      !! Adds a named column or replaces an existing column of conformable size.
      type(tibble_type), intent(in) :: data !! Input table.
      type(vctr_type), intent(in) :: new_column !! Named replacement or addition.
      type(tibble_type) :: out

      if (new_column%size() /= data%nrow()) error stop "mutate: column size differs from row count"
      if (data%has_name(new_column%name)) then
         out = replace_column(data, new_column)
      else
         out = add_column(data, new_column)
      end if
   end function mutate_one

   function mutate_many(data, new_columns) result(out)
      !! Applies named column additions or replacements from left to right.
      type(tibble_type), intent(in) :: data !! Input table.
      type(vctr_type), intent(in) :: new_columns(:) !! Named columns to install.
      type(tibble_type) :: out
      integer :: j

      out = data
      do j = 1, size(new_columns)
         out = mutate_one(out, new_columns(j))
      end do
   end function mutate_many

   function slice(data, rows) result(out)
      !! Selects rows by one-based positions, preserving repetition and order.
      type(tibble_type), intent(in) :: data !! Input table.
      integer, intent(in) :: rows(:)        !! One-based row positions.
      type(tibble_type) :: out

      out = slice_rows(data, rows)
   end function slice

   function slice_head(data, n) result(out)
      !! Keeps up to the first `n` rows.
      type(tibble_type), intent(in) :: data !! Input table.
      integer, intent(in) :: n             !! Nonnegative maximum row count.
      type(tibble_type) :: out
      integer :: count_rows, i
      integer, allocatable :: rows(:)

      if (n < 0) error stop "slice_head: n must be nonnegative"
      count_rows = min(n, data%nrow())
      allocate (rows(count_rows))
      do i = 1, count_rows
         rows(i) = i
      end do
      out = slice_rows(data, rows)
   end function slice_head

   function slice_tail(data, n) result(out)
      !! Keeps up to the last `n` rows.
      type(tibble_type), intent(in) :: data !! Input table.
      integer, intent(in) :: n             !! Nonnegative maximum row count.
      type(tibble_type) :: out
      integer :: count_rows, first_row, i
      integer, allocatable :: rows(:)

      if (n < 0) error stop "slice_tail: n must be nonnegative"
      count_rows = min(n, data%nrow())
      first_row = data%nrow() - count_rows + 1
      allocate (rows(count_rows))
      do i = 1, count_rows
         rows(i) = first_row + i - 1
      end do
      out = slice_rows(data, rows)
   end function slice_tail

   function arrange(data, names, descending) result(out)
      !! Stably orders rows by selected columns, with missing values always last.
      type(tibble_type), intent(in) :: data !! Input table.
      character(len=*), intent(in) :: names(:) !! Exact sort-key names.
      logical, intent(in), optional :: descending(:) !! Per-key descending flags.
      type(tibble_type) :: out
      integer, allocatable :: keys(:), order(:)
      logical, allocatable :: reverse(:)
      integer :: i, j, saved

      keys = resolve_columns(data, names)
      allocate (reverse(size(keys)), source=.false.)
      if (present(descending)) then
         if (size(descending) /= size(keys)) error stop "arrange: descending size differs from keys"
         reverse = descending
      end if
      allocate (order(data%nrow()))
      do i = 1, data%nrow()
         order(i) = i
      end do
      do i = 2, data%nrow()
         saved = order(i)
         j = i - 1
         do while (j >= 1)
            if (.not. row_precedes(data, saved, order(j), keys, reverse)) exit
            order(j + 1) = order(j)
            j = j - 1
         end do
         order(j + 1) = saved
      end do
      out = slice_rows(data, order)
   end function arrange

   function distinct(data, names, keep_all) result(out)
      !! Retains the first row for each distinct selected-key combination.
      type(tibble_type), intent(in) :: data !! Input table.
      character(len=*), intent(in), optional :: names(:) !! Keys; default all columns.
      logical, intent(in), optional :: keep_all !! Whether to keep non-key columns.
      type(tibble_type) :: out, source
      integer, allocatable :: keys(:), rows(:), buffer(:)
      logical :: retain_all, duplicate
      integer :: i, j, number_kept

      if (present(names)) then
         keys = resolve_columns(data, names)
      else
         allocate (keys(data%ncol()))
         do i = 1, data%ncol()
            keys(i) = i
         end do
      end if
      retain_all = .false.
      if (present(keep_all)) retain_all = keep_all
      allocate (buffer(data%nrow()))
      number_kept = 0
      do i = 1, data%nrow()
         duplicate = .false.
         do j = 1, number_kept
            if (rows_equal(data, keys, i, buffer(j))) then
               duplicate = .true.
               exit
            end if
         end do
         if (.not. duplicate) then
            number_kept = number_kept + 1
            buffer(number_kept) = i
         end if
      end do
      allocate (rows(number_kept))
      if (number_kept > 0) rows = buffer(:number_kept)
      if (retain_all .or. .not. present(names)) then
         source = data
      else
         source = new_tibble(data%columns(keys), nrow=data%nrow())
      end if
      out = slice_rows(source, rows)
   end function distinct

   function bind_rows(frames) result(out)
      !! Row-binds frames with identical ordered names and compatible column types.
      type(tibble_type), intent(in) :: frames(:) !! Frames to combine.
      type(tibble_type) :: out
      out = vec_rbind(frames)
   end function bind_rows

   function bind_cols(frames) result(out)
      !! Column-binds frames with unique names and size-one recycling.
      type(tibble_type), intent(in) :: frames(:) !! Frames to combine.
      type(tibble_type) :: out
      out = vec_cbind(frames)
   end function bind_cols

   pure elemental function pull(data, name) result(vector)
      !! Returns one exact named column.
      type(tibble_type), intent(in) :: data !! Input table.
      character(len=*), intent(in) :: name  !! Exact column name.
      type(vctr_type) :: vector
      integer :: index

      index = data%column_index(name)
      if (index == 0) error stop "pull: column not found: " // name
      vector = data%columns(index)
   end function pull

   pure function resolve_columns(data, names) result(indices)
      !! Resolves unique exact column names to positions.
      type(tibble_type), intent(in) :: data !! Table to search.
      character(len=*), intent(in) :: names(:) !! Exact names.
      integer, allocatable :: indices(:)
      integer :: i, j

      allocate (indices(size(names)))
      do j = 1, size(names)
         indices(j) = data%column_index(names(j))
         if (indices(j) == 0) error stop "column not found: " // trim(names(j))
         do i = 1, j - 1
            if (indices(i) == indices(j)) error stop "column selected more than once"
         end do
      end do
   end function resolve_columns

   pure logical function has_duplicates(indices) result(duplicate)
      !! Reports whether an integer selection repeats a position.
      integer, intent(in) :: indices(:) !! Positions to inspect.
      integer :: i, j
      duplicate = .false.
      do i = 1, size(indices)
         do j = 1, i - 1
            if (indices(i) == indices(j)) then
               duplicate = .true.
               return
            end if
         end do
      end do
   end function has_duplicates

   pure logical function rows_equal(data, keys, left, right) result(equal)
      !! Compares two rows across key columns, treating paired missing values as equal.
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
      !! Compares two positions of one vector under grouping equality.
      type(vctr_type), intent(in) :: vector !! Vector containing values.
      integer, intent(in) :: left  !! First element position.
      integer, intent(in) :: right !! Second element position.

      if (vector%missing(left) .neqv. vector%missing(right)) then
         equal = .false.
         return
      end if
      if (vector%missing(left)) then
         equal = .true.
         return
      end if
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
         error stop "row comparison: unsupported type"
      end select
   end function elements_equal

   pure logical function row_precedes(data, left, right, keys, descending) result(precedes)
      !! Applies lexicographic row ordering across keys.
      type(tibble_type), intent(in) :: data !! Table containing keys.
      integer, intent(in) :: left  !! Candidate preceding row.
      integer, intent(in) :: right !! Candidate following row.
      integer, intent(in) :: keys(:)        !! Key-column positions.
      logical, intent(in) :: descending(:)  !! Per-key descending flags.
      integer :: comparison, j

      precedes = .false.
      do j = 1, size(keys)
         comparison = compare_elements(data%columns(keys(j)), left, right)
         if (comparison == 0) cycle
         if (comparison == 2) return
         if (comparison == -2) then
            precedes = .true.
            return
         end if
         precedes = merge(comparison > 0, comparison < 0, descending(j))
         return
      end do
   end function row_precedes

   pure elemental integer function compare_elements(vector, left, right) result(comparison)
      !! Returns -1/0/1 for values and -2/2 for missing-before/value-before-missing.
      type(vctr_type), intent(in) :: vector !! Vector containing values.
      integer, intent(in) :: left  !! First element position.
      integer, intent(in) :: right !! Second element position.

      if (vector%missing(left)) then
         comparison = merge(0, 2, vector%missing(right))
         return
      else if (vector%missing(right)) then
         comparison = -2
         return
      end if
      if (elements_equal(vector, left, right)) then
         comparison = 0
         return
      end if
      select case (vector%type_code)
      case (vctrs_integer)
         comparison = merge(-1, 1, vector%integer_values(left) < vector%integer_values(right))
      case (vctrs_real)
         comparison = merge(-1, 1, vector%real_values(left) < vector%real_values(right))
      case (vctrs_logical)
         comparison = merge(-1, 1, .not. vector%logical_values(left))
      case (vctrs_character)
         comparison = merge(-1, 1, vector%character_values(left) < vector%character_values(right))
      case default
         error stop "arrange: unsupported type"
      end select
   end function compare_elements

end module dplyr_verbs
