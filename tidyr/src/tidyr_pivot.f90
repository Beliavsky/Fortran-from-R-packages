! SPDX-License-Identifier: MIT
! SPDX-FileComment: Pivot operations for the Fortran tidyr translation.
module tidyr_pivot
   !! Reshapes rectangular tables between long and wide representations.
   use tibble, only: make_column, new_tibble, tibble_character, tibble_column, tibble_type
   use tidyr_support, only: column_indices, row_keys_equal
   use vctrs, only: dp, vec_c, vec_copy_element, vec_slice
   implicit none
   private

   public :: gather, pivot_longer, pivot_wider, spread

   interface pivot_longer
      module procedure pivot_longer_names, pivot_longer_indices
   end interface pivot_longer

contains

   function pivot_longer_names(data, columns, names_to, values_to, cols_vary) result(output)
      !! Stacks selected columns while repeating all unselected identifier columns.
      type(tibble_type), intent(in) :: data    !! Input wide table.
      character(len=*), intent(in) :: columns(:) !! Columns to stack.
      character(len=*), intent(in), optional :: names_to  !! Output name column; default `name`.
      character(len=*), intent(in), optional :: values_to !! Output value column; default `value`.
      character(len=*), intent(in), optional :: cols_vary !! `fastest` or `slowest` row ordering.
      type(tibble_type) :: output
      type(tibble_column), allocatable :: out_columns(:), pieces(:)
      character(len=:), allocatable :: name_column, value_column
      character(len=:), allocatable :: labels(:)
      integer, allocatable :: selected(:), identifiers(:), rows(:)
      logical, allocatable :: is_selected(:)
      integer :: i, j, k, max_name, total
      character(len=:), allocatable :: ordering

      if (size(columns) == 0) error stop 'pivot_longer: no columns selected'
      name_column = 'name'
      value_column = 'value'
      if (present(names_to)) name_column = trim(names_to)
      if (present(values_to)) value_column = trim(values_to)
      ordering = 'fastest'
      if (present(cols_vary)) ordering = trim(cols_vary)
      if (ordering /= 'fastest' .and. ordering /= 'slowest') then
         error stop 'pivot_longer: cols_vary must be fastest or slowest'
      end if
      selected = column_indices(data, columns)
      allocate (is_selected(data%ncol()), source=.false.)
      do j = 1, size(columns)
         is_selected(selected(j)) = .true.
      end do
      identifiers = pack([(j, j = 1, data%ncol())], .not. is_selected)
      total = data%nrow() * size(selected)
      allocate (rows(total))
      max_name = max(1, maxval([(len_trim(data%columns(selected(j))%name), j = 1, size(selected))]))
      allocate (character(len=max_name) :: labels(total))
      allocate (pieces(total), out_columns(size(identifiers) + 2))
      k = 0
      if (ordering == 'fastest') then
         do i = 1, data%nrow()
            do j = 1, size(selected)
               k = k + 1
               rows(k) = i
               labels(k) = data%columns(selected(j))%name
               pieces(k) = vec_slice(data%columns(selected(j)), [i])
            end do
         end do
      else
         do j = 1, size(selected)
            do i = 1, data%nrow()
               k = k + 1
               rows(k) = i
               labels(k) = data%columns(selected(j))%name
               pieces(k) = vec_slice(data%columns(selected(j)), [i])
            end do
         end do
      end if
      do j = 1, size(identifiers)
         out_columns(j) = vec_slice(data%columns(identifiers(j)), rows)
      end do
      out_columns(size(identifiers) + 1) = make_column(name_column, labels)
      out_columns(size(identifiers) + 2) = vec_c(pieces)
      out_columns(size(identifiers) + 2)%name = value_column
      output = new_tibble(out_columns, nrow=total)
   end function pivot_longer_names

   function pivot_longer_indices(data, indices, names_to, values_to, cols_vary) result(output)
      !! Stacks positions returned by tidyselect.
      type(tibble_type), intent(in) :: data !! Input wide table.
      integer, intent(in) :: indices(:)     !! One-based columns to stack.
      character(len=*), intent(in), optional :: names_to !! Output name column.
      character(len=*), intent(in), optional :: values_to !! Output value column.
      character(len=*), intent(in), optional :: cols_vary !! Row ordering.
      type(tibble_type) :: output
      character(len=:), allocatable :: selected_names(:)
      integer :: j, width
      if (any(indices < 1) .or. any(indices > data%ncol())) error stop 'pivot_longer: index out of range'
      width = 1
      do j = 1, size(indices)
         width = max(width, len(data%columns(indices(j))%name))
      end do
      allocate (character(len=width) :: selected_names(size(indices)))
      do j = 1, size(indices)
         selected_names(j) = data%columns(indices(j))%name
      end do
      output = pivot_longer_names(data, selected_names, names_to, values_to, cols_vary)
   end function pivot_longer_indices

   function gather(data, key, value, columns) result(output)
      !! Legacy alias for `pivot_longer`.
      type(tibble_type), intent(in) :: data !! Input wide table.
      character(len=*), intent(in) :: key   !! Output key-column name.
      character(len=*), intent(in) :: value !! Output value-column name.
      character(len=*), intent(in) :: columns(:) !! Columns to stack.
      type(tibble_type) :: output
      output = pivot_longer_names(data, columns, names_to=key, values_to=value)
   end function gather

   function pivot_wider(data, name_column, value_column, identifiers) result(output)
      !! Widens a unique long table using character names and one value column.
      type(tibble_type), intent(in) :: data      !! Input long table.
      character(len=*), intent(in) :: identifiers(:) !! Identifier columns.
      character(len=*), intent(in) :: name_column !! Character column supplying output names.
      character(len=*), intent(in) :: value_column !! Column supplying cell values.
      type(tibble_type) :: output
      type(tibble_column), allocatable :: out_columns(:)
      integer, allocatable :: id_indices(:), representatives(:), name_rows(:)
      logical, allocatable :: is_new_id(:), is_new_name(:)
      integer :: i, j, k, id_count, name_count, name_index, value_index, target

      name_index = data%column_index(name_column)
      value_index = data%column_index(value_column)
      if (name_index == 0 .or. value_index == 0) error stop 'pivot_wider: source column not found'
      if (data%columns(name_index)%type_code /= tibble_character) then
         error stop 'pivot_wider: names_from must be character'
      end if
      id_indices = column_indices(data, identifiers)
      allocate (is_new_id(data%nrow()), source=.true.)
      allocate (is_new_name(data%nrow()), source=.true.)
      do i = 1, data%nrow()
         do j = 1, i - 1
            if (row_keys_equal(data, id_indices, i, j)) then
               is_new_id(i) = .false.
               exit
            end if
         end do
         do j = 1, i - 1
            if (same_element(data%columns(name_index), i, data%columns(name_index), j)) then
               is_new_name(i) = .false.
               exit
            end if
         end do
      end do
      representatives = pack([(i, i = 1, data%nrow())], is_new_id)
      name_rows = pack([(i, i = 1, data%nrow())], is_new_name)
      id_count = size(representatives)
      name_count = size(name_rows)
      allocate (out_columns(size(identifiers) + name_count))
      do j = 1, size(identifiers)
         out_columns(j) = vec_slice(data%columns(id_indices(j)), representatives)
      end do
      do j = 1, name_count
         call initialize_missing_vector(data%columns(value_index), &
            trim(data%columns(name_index)%character_values(name_rows(j))), id_count, &
            out_columns(size(identifiers) + j))
      end do
      do i = 1, data%nrow()
         target = find_key_row(data, id_indices, i, representatives)
         k = find_name_row(data%columns(name_index), i, name_rows)
         if (.not. out_columns(size(identifiers) + k)%missing(target)) then
            error stop 'pivot_wider: duplicate identifier/name combination'
         end if
         call vec_copy_element(data%columns(value_index), i, out_columns(size(identifiers) + k), target)
      end do
      call move_alloc(out_columns, output%columns)
      output%number_of_rows = id_count
   end function pivot_wider

   function spread(data, key, value, id_cols) result(output)
      !! Legacy alias for `pivot_wider`.
      type(tibble_type), intent(in) :: data !! Input long table.
      character(len=*), intent(in) :: key   !! Column supplying output names.
      character(len=*), intent(in) :: value !! Column supplying values.
      character(len=*), intent(in) :: id_cols(:) !! Identifier columns.
      type(tibble_type) :: output
      output = pivot_wider(data, key, value, id_cols)
   end function spread

   pure elemental logical function same_element(left, left_row, right, right_row) result(equal)
      !! Compares two equal-type vector elements, treating two missing values as equal.
      type(tibble_column), intent(in) :: left  !! First column.
      integer, intent(in) :: left_row         !! First row.
      type(tibble_column), intent(in) :: right !! Second column.
      integer, intent(in) :: right_row        !! Second row.
      if (left%type_code /= right%type_code) then
         equal = .false.
      else if (left%missing(left_row) .or. right%missing(right_row)) then
         equal = left%missing(left_row) .and. right%missing(right_row)
      else
         select case (left%type_code)
         case (1)
            equal = left%integer_values(left_row) == right%integer_values(right_row)
         case (2)
            equal = left%real_values(left_row) <= right%real_values(right_row) .and. &
                    right%real_values(right_row) <= left%real_values(left_row)
         case (3)
            equal = left%logical_values(left_row) .eqv. right%logical_values(right_row)
         case (4)
            equal = left%character_values(left_row) == right%character_values(right_row)
         case default
            equal = .false.
         end select
      end if
   end function same_element

   pure integer function find_key_row(data, indices, row, representatives) result(location)
      !! Finds the output row for one input identifier key.
      type(tibble_type), intent(in) :: data !! Table containing keys.
      integer, intent(in) :: indices(:)     !! Identifier column indices.
      integer, intent(in) :: row            !! Input row.
      integer, intent(in) :: representatives(:) !! Representative input rows.
      integer :: j
      location = 0
      do j = 1, size(representatives)
         if (row_keys_equal(data, indices, row, representatives(j))) then
            location = j
            return
         end if
      end do
   end function find_key_row

   pure integer function find_name_row(vector, row, representatives) result(location)
      !! Finds the output column for one input name.
      type(tibble_column), intent(in) :: vector !! Name column.
      integer, intent(in) :: row                !! Input row.
      integer, intent(in) :: representatives(:) !! Representative name rows.
      integer :: j
      location = 0
      do j = 1, size(representatives)
         if (same_element(vector, row, vector, representatives(j))) then
            location = j
            return
         end if
      end do
   end function find_name_row

   pure subroutine initialize_missing_vector(prototype, name, n, vector)
      !! Initializes missing output storage without a derived-type function temporary.
      type(tibble_column), intent(in) :: prototype !! Vector supplying the storage type.
      character(len=*), intent(in) :: name         !! Output vector name.
      integer, intent(in) :: n                     !! Output vector size.
      type(tibble_column), intent(out) :: vector   !! Initialized missing vector.
      vector%name = name
      vector%type_code = prototype%type_code
      allocate (vector%missing(n), source=.true.)
      select case (prototype%type_code)
      case (1)
         allocate (vector%integer_values(n), source=0)
      case (2)
         allocate (vector%real_values(n), source=0.0_dp)
      case (3)
         allocate (vector%logical_values(n), source=.false.)
      case (4)
         allocate (character(len=len(prototype%character_values)) :: vector%character_values(n))
         vector%character_values = ''
      case default
         error stop 'pivot_wider: unsupported value type'
      end select
   end subroutine initialize_missing_vector

end module tidyr_pivot
