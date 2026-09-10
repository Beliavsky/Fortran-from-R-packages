! SPDX-License-Identifier: MIT
module tibble_operations
   use vctrs, only: vec_insert, vec_is, vec_slice
   use tibble_types, only: dp, make_column, tibble_character, tibble_column, tibble_integer, &
      tibble_logical, tibble_real, tibble_type
   implicit none
   private

   public :: add_column, add_rows, as_tibble, deframe_character, deframe_integer
   public :: deframe_logical, deframe_real, drop_columns, enframe, filter_rows
   public :: get_character, get_integer, get_logical, get_missing, get_real, glimpse
   public :: new_tibble
   public :: print_tibble, repair_names, replace_column, select_columns, slice_columns
   public :: slice_rows, validate_tibble

   interface as_tibble
      module procedure as_tibble_integer_matrix, as_tibble_real_matrix
      module procedure as_tibble_logical_matrix, as_tibble_character_matrix
   end interface as_tibble

   interface enframe
      module procedure enframe_integer, enframe_real, enframe_logical, enframe_character
   end interface enframe

contains

   function new_tibble(columns, nrow, name_repair) result(table)
      !! Constructs a rectangular heterogeneous table and checks its invariants.
      type(tibble_column), intent(in) :: columns(:) !! Typed columns to store.
      integer, intent(in), optional :: nrow         !! Row count for a zero-column table.
      character(len=*), intent(in), optional :: &
         name_repair !! `check_unique`, `unique`, or `minimal`.
      type(tibble_type) :: table
      character(len=:), allocatable :: message, repair
      integer :: j

      table%columns = columns
      if (size(columns) > 0) then
         table%number_of_rows = columns(1)%size()
         if (present(nrow)) then
            if (nrow /= table%number_of_rows) then
               error stop "new_tibble: nrow disagrees with column size"
            end if
         end if
      else
         table%number_of_rows = 0
         if (present(nrow)) table%number_of_rows = nrow
         if (table%number_of_rows < 0) error stop "new_tibble: nrow must be nonnegative"
      end if

      repair = "check_unique"
      if (present(name_repair)) repair = trim(name_repair)
      select case (repair)
      case ("check_unique")
         if (.not. validate_tibble(table, message)) error stop "new_tibble: " // message
      case ("unique")
         call repair_names(table)
         if (.not. validate_tibble(table, message)) error stop "new_tibble: " // message
      case ("minimal")
         do j = 1, table%ncol()
            if (.not. allocated(table%columns(j)%name)) table%columns(j)%name = ""
         end do
         if (.not. validate_tibble(table, message, check_names=.false.)) then
            error stop "new_tibble: " // message
         end if
      case default
         error stop "new_tibble: unsupported name-repair mode"
      end select
   end function new_tibble

   logical function validate_tibble(table, message, check_names) result(valid)
      !! Checks column allocation, types, names, masks, and common row sizes.
      type(tibble_type), intent(in) :: table        !! Table to validate.
      character(len=:), allocatable, intent(out), optional :: message !! Failure explanation.
      logical, intent(in), optional :: check_names !! Whether names must be nonempty and unique.
      logical :: strict_names
      integer :: i, j

      valid = .false.
      strict_names = .true.
      if (present(check_names)) strict_names = check_names
      if (table%number_of_rows < 0) then
         call fail("negative row count")
         return
      end if
      do j = 1, table%ncol()
         if (.not. allocated(table%columns(j)%name)) then
            call fail("unallocated column name")
            return
         end if
         if (strict_names .and. len(table%columns(j)%name) == 0) then
            call fail("empty column name")
            return
         end if
         if (strict_names) then
            do i = 1, j - 1
               if (table%columns(i)%name == table%columns(j)%name) then
                  call fail("duplicate column name: " // table%columns(j)%name)
                  return
               end if
            end do
         end if
         if (.not. allocated(table%columns(j)%missing)) then
            call fail("unallocated missing-value mask: " // table%columns(j)%name)
            return
         end if
         if (size(table%columns(j)%missing) /= table%number_of_rows) then
            call fail("missing-value mask has wrong size: " // table%columns(j)%name)
            return
         end if
         if (.not. vec_is(table%columns(j))) then
            call fail("invalid storage for column: " // table%columns(j)%name)
            return
         end if
         if (table%columns(j)%size() /= table%number_of_rows) then
            call fail("column size differs from table row count: " // table%columns(j)%name)
            return
         end if
      end do
      valid = .true.
      if (present(message)) message = ""

   contains

      subroutine fail(text)
         !! Stores a validation failure message when requested.
         character(len=*), intent(in) :: text !! Explanation to return.

         if (present(message)) message = text
      end subroutine fail

   end function validate_tibble

   pure elemental subroutine repair_names(table)
      !! Makes empty and duplicate names unique using stable `...j` suffixes.
      type(tibble_type), intent(inout) :: table !! Table whose names are repaired.
      character(len=:), allocatable :: candidate, original
      integer :: i, j, suffix_number
      logical :: duplicate

      do j = 1, table%ncol()
         if (allocated(table%columns(j)%name)) then
            original = table%columns(j)%name
         else
            original = "..."
         end if
         if (len(original) == 0) original = "..."
         candidate = original
         suffix_number = j
         do
            duplicate = .false.
            if (len(candidate) == 0) duplicate = .true.
            if (.not. duplicate .and. j > 1) then
               do i = 1, j - 1
                  if (table%columns(i)%name == candidate) then
                     duplicate = .true.
                     exit
                  end if
               end do
            end if
            if (.not. duplicate) exit
            candidate = original // "..." // positive_integer_string(suffix_number)
            suffix_number = suffix_number + 1
         end do
         table%columns(j)%name = candidate
      end do
   end subroutine repair_names

   pure function positive_integer_string(value) result(text)
      !! Converts a positive integer to decimal text without performing I/O.
      integer, intent(in) :: value !! Positive integer to format.
      character(len=:), allocatable :: text
      integer :: digits, i, remaining

      if (value < 1) error stop "positive_integer_string: value must be positive"
      digits = 1
      remaining = value
      do while (remaining >= 10)
         digits = digits + 1
         remaining = remaining / 10
      end do
      allocate (character(len=digits) :: text)
      remaining = value
      do i = digits, 1, -1
         text(i:i) = achar(iachar('0') + mod(remaining, 10))
         remaining = remaining / 10
      end do
   end function positive_integer_string

   function as_tibble_integer_matrix(names, values) result(table)
      !! Converts an integer matrix, interpreted as rows by columns, to a tibble.
      character(len=*), intent(in) :: names(:) !! Column names.
      integer, intent(in) :: values(:,:)        !! Matrix values `(nrow,ncol)`.
      type(tibble_type) :: table
      type(tibble_column), allocatable :: columns(:)
      integer :: j

      if (size(names) /= size(values, 2)) then
         error stop "as_tibble: name count differs from matrix columns"
      end if
      allocate (columns(size(names)))
      do j = 1, size(names)
         columns(j) = make_column(names(j), values(:, j))
      end do
      table = new_tibble(columns, nrow=size(values, 1))
   end function as_tibble_integer_matrix

   function as_tibble_real_matrix(names, values) result(table)
      !! Converts a real matrix, interpreted as rows by columns, to a tibble.
      character(len=*), intent(in) :: names(:) !! Column names.
      real(dp), intent(in) :: values(:,:)       !! Matrix values `(nrow,ncol)`.
      type(tibble_type) :: table
      type(tibble_column), allocatable :: columns(:)
      integer :: j

      if (size(names) /= size(values, 2)) then
         error stop "as_tibble: name count differs from matrix columns"
      end if
      allocate (columns(size(names)))
      do j = 1, size(names)
         columns(j) = make_column(names(j), values(:, j))
      end do
      table = new_tibble(columns, nrow=size(values, 1))
   end function as_tibble_real_matrix

   function as_tibble_logical_matrix(names, values) result(table)
      !! Converts a logical matrix, interpreted as rows by columns, to a tibble.
      character(len=*), intent(in) :: names(:) !! Column names.
      logical, intent(in) :: values(:,:)        !! Matrix values `(nrow,ncol)`.
      type(tibble_type) :: table
      type(tibble_column), allocatable :: columns(:)
      integer :: j

      if (size(names) /= size(values, 2)) then
         error stop "as_tibble: name count differs from matrix columns"
      end if
      allocate (columns(size(names)))
      do j = 1, size(names)
         columns(j) = make_column(names(j), values(:, j))
      end do
      table = new_tibble(columns, nrow=size(values, 1))
   end function as_tibble_logical_matrix

   function as_tibble_character_matrix(names, values) result(table)
      !! Converts a character matrix, interpreted as rows by columns, to a tibble.
      character(len=*), intent(in) :: names(:)  !! Column names.
      character(len=*), intent(in) :: values(:,:) !! Matrix values `(nrow,ncol)`.
      type(tibble_type) :: table
      type(tibble_column), allocatable :: columns(:)
      integer :: j

      if (size(names) /= size(values, 2)) then
         error stop "as_tibble: name count differs from matrix columns"
      end if
      allocate (columns(size(names)))
      do j = 1, size(names)
         columns(j) = make_column(names(j), values(:, j))
      end do
      table = new_tibble(columns, nrow=size(values, 1))
   end function as_tibble_character_matrix

   pure function get_integer(table, name) result(values)
      !! Extracts an integer column by exact name.
      type(tibble_type), intent(in) :: table      !! Table containing the column.
      character(len=*), intent(in) :: name        !! Exact column name.
      integer, allocatable :: values(:)
      integer :: index

      index = require_column(table, name, tibble_integer)
      values = table%columns(index)%integer_values
   end function get_integer

   pure function get_real(table, name) result(values)
      !! Extracts a real column by exact name.
      type(tibble_type), intent(in) :: table      !! Table containing the column.
      character(len=*), intent(in) :: name        !! Exact column name.
      real(dp), allocatable :: values(:)
      integer :: index

      index = require_column(table, name, tibble_real)
      values = table%columns(index)%real_values
   end function get_real

   pure function get_logical(table, name) result(values)
      !! Extracts a logical column by exact name.
      type(tibble_type), intent(in) :: table      !! Table containing the column.
      character(len=*), intent(in) :: name        !! Exact column name.
      logical, allocatable :: values(:)
      integer :: index

      index = require_column(table, name, tibble_logical)
      values = table%columns(index)%logical_values
   end function get_logical

   pure function get_character(table, name) result(values)
      !! Extracts a character column by exact name.
      type(tibble_type), intent(in) :: table      !! Table containing the column.
      character(len=*), intent(in) :: name        !! Exact column name.
      character(len=:), allocatable :: values(:)
      integer :: index

      index = require_column(table, name, tibble_character)
      values = table%columns(index)%character_values
   end function get_character

   pure function get_missing(table, name) result(missing)
      !! Extracts the missing-value mask for a column of any supported type.
      type(tibble_type), intent(in) :: table !! Table containing the column.
      character(len=*), intent(in) :: name   !! Exact column name.
      logical, allocatable :: missing(:)
      integer :: index

      index = table%column_index(name)
      if (index == 0) error stop "column not found: " // name
      missing = table%columns(index)%missing
   end function get_missing

   pure elemental integer function require_column(table, name, type_code) result(index)
      !! Locates a column and enforces its requested static type.
      type(tibble_type), intent(in) :: table !! Table to search.
      character(len=*), intent(in) :: name   !! Exact column name.
      integer, intent(in) :: type_code       !! Required column type tag.

      index = table%column_index(name)
      if (index == 0) error stop "column not found: " // name
      if (table%columns(index)%type_code /= type_code) then
         error stop "column has an incompatible type: " // name
      end if
   end function require_column

   function slice_rows(table, indices) result(out)
      !! Selects rows by one-based integer indices, preserving order and repetition.
      type(tibble_type), intent(in) :: table !! Source table.
      integer, intent(in) :: indices(:)      !! One-based row indices.
      type(tibble_type) :: out
      type(tibble_column), allocatable :: columns(:)
      integer :: j

      if (any(indices < 1) .or. any(indices > table%nrow())) then
         error stop "slice_rows: row index out of bounds"
      end if
      allocate (columns(table%ncol()))
      do j = 1, table%ncol()
         columns(j) = vec_slice(table%columns(j), indices)
      end do
      out = new_tibble(columns, nrow=size(indices), name_repair="minimal")
   end function slice_rows

   function filter_rows(table, keep) result(out)
      !! Selects rows for which a conformable logical mask is true.
      type(tibble_type), intent(in) :: table !! Source table.
      logical, intent(in) :: keep(:)         !! Row-selection mask `(nrow)`.
      type(tibble_type) :: out
      integer, allocatable :: indices(:)
      integer :: i, k

      if (size(keep) /= table%nrow()) error stop "filter_rows: mask size differs from row count"
      allocate (indices(count(keep)))
      k = 0
      do i = 1, size(keep)
         if (keep(i)) then
            k = k + 1
            indices(k) = i
         end if
      end do
      out = slice_rows(table, indices)
   end function filter_rows

   function slice_columns(table, indices) result(out)
      !! Selects columns by one-based indices, preserving order and repetition.
      type(tibble_type), intent(in) :: table !! Source table.
      integer, intent(in) :: indices(:)      !! One-based column indices.
      type(tibble_type) :: out

      if (any(indices < 1) .or. any(indices > table%ncol())) then
         error stop "slice_columns: column index out of bounds"
      end if
      out = new_tibble(table%columns(indices), nrow=table%nrow(), name_repair="minimal")
   end function slice_columns

   function select_columns(table, names) result(out)
      !! Selects columns by exact names without partial matching.
      type(tibble_type), intent(in) :: table !! Source table.
      character(len=*), intent(in) :: names(:) !! Names in requested output order.
      type(tibble_type) :: out
      integer, allocatable :: indices(:)
      integer :: j

      allocate (indices(size(names)))
      do j = 1, size(names)
         indices(j) = table%column_index(names(j))
         if (indices(j) == 0) error stop "select_columns: column not found: " // names(j)
      end do
      out = slice_columns(table, indices)
   end function select_columns

   function drop_columns(table, names) result(out)
      !! Removes named columns while preserving all other columns and their order.
      type(tibble_type), intent(in) :: table !! Source table.
      character(len=*), intent(in) :: names(:) !! Exact names to remove.
      type(tibble_type) :: out
      integer, allocatable :: indices(:)
      logical, allocatable :: keep(:)
      integer :: i, j, k

      allocate (keep(table%ncol()), source=.true.)
      do i = 1, size(names)
         j = table%column_index(names(i))
         if (j == 0) error stop "drop_columns: column not found: " // names(i)
         keep(j) = .false.
      end do
      allocate (indices(count(keep)))
      k = 0
      do j = 1, table%ncol()
         if (keep(j)) then
            k = k + 1
            indices(k) = j
         end if
      end do
      out = slice_columns(table, indices)
   end function drop_columns

   function add_column(table, column, before) result(out)
      !! Inserts a new column, rejecting duplicate names and incompatible sizes.
      type(tibble_type), intent(in) :: table   !! Source table.
      type(tibble_column), intent(in) :: column !! Column to insert.
      integer, intent(in), optional :: before  !! One-based insertion position; default is append.
      type(tibble_type) :: out
      type(tibble_column), allocatable :: columns(:)
      integer :: position

      if (table%has_name(column%name)) then
         error stop "add_column: column already exists: " // column%name
      end if
      if (column%size() /= table%nrow()) error stop "add_column: column size differs from row count"
      position = table%ncol() + 1
      if (present(before)) position = before
      if (position < 1 .or. position > table%ncol() + 1) then
         error stop "add_column: position out of bounds"
      end if
      allocate (columns(table%ncol() + 1))
      if (position > 1) columns(:position - 1) = table%columns(:position - 1)
      columns(position) = column
      if (position <= table%ncol()) columns(position + 1:) = table%columns(position:)
      out = new_tibble(columns, nrow=table%nrow())
   end function add_column

   pure elemental function replace_column(table, column) result(out)
      !! Replaces an existing column by exact name while allowing its type to change.
      type(tibble_type), intent(in) :: table    !! Source table.
      type(tibble_column), intent(in) :: column !! Replacement column.
      type(tibble_type) :: out
      integer :: index

      index = table%column_index(column%name)
      if (index == 0) error stop "replace_column: column not found: " // column%name
      if (column%size() /= table%nrow()) then
         error stop "replace_column: column size differs from row count"
      end if
      out = table
      out%columns(index) = column
   end function replace_column

   function add_rows(table, rows, before) result(out)
      !! Inserts rows from a schema-compatible tibble.
      type(tibble_type), intent(in) :: table !! Destination table.
      type(tibble_type), intent(in) :: rows  !! Rows with identical names and types.
      integer, intent(in), optional :: before !! One-based insertion row; default is append.
      type(tibble_type) :: out
      type(tibble_column), allocatable :: columns(:)
      integer :: j, position

      if (rows%ncol() /= table%ncol()) error stop "add_rows: column counts differ"
      do j = 1, table%ncol()
         if (rows%columns(j)%name /= table%columns(j)%name) then
            error stop "add_rows: column names or order differ"
         end if
      end do
      position = table%nrow() + 1
      if (present(before)) position = before
      if (position < 1 .or. position > table%nrow() + 1) then
         error stop "add_rows: position out of bounds"
      end if
      allocate (columns(table%ncol()))
      do j = 1, table%ncol()
         columns(j) = vec_insert(table%columns(j), rows%columns(j), position)
         columns(j)%name = table%columns(j)%name
      end do
      out = new_tibble(columns, nrow=table%nrow() + rows%nrow())
   end function add_rows

   function enframe_integer(values, names, name, value) result(table)
      !! Converts an integer vector and optional labels to a one- or two-column tibble.
      integer, intent(in) :: values(:)             !! Values to frame.
      character(len=*), intent(in), optional :: names(:) !! Optional element names.
      character(len=*), intent(in), optional :: name  !! Output label-column name.
      character(len=*), intent(in), optional :: value !! Output value-column name.
      type(tibble_type) :: table

      table = enframe_column(make_column(output_value_name(value), values), names, name)
   end function enframe_integer

   function enframe_real(values, names, name, value) result(table)
      !! Converts a real vector and optional labels to a one- or two-column tibble.
      real(dp), intent(in) :: values(:)            !! Values to frame.
      character(len=*), intent(in), optional :: names(:) !! Optional element names.
      character(len=*), intent(in), optional :: name  !! Output label-column name.
      character(len=*), intent(in), optional :: value !! Output value-column name.
      type(tibble_type) :: table

      table = enframe_column(make_column(output_value_name(value), values), names, name)
   end function enframe_real

   function enframe_logical(values, names, name, value) result(table)
      !! Converts a logical vector and optional labels to a one- or two-column tibble.
      logical, intent(in) :: values(:)             !! Values to frame.
      character(len=*), intent(in), optional :: names(:) !! Optional element names.
      character(len=*), intent(in), optional :: name  !! Output label-column name.
      character(len=*), intent(in), optional :: value !! Output value-column name.
      type(tibble_type) :: table

      table = enframe_column(make_column(output_value_name(value), values), names, name)
   end function enframe_logical

   function enframe_character(values, names, name, value) result(table)
      !! Converts a character vector and optional labels to a one- or two-column tibble.
      character(len=*), intent(in) :: values(:)    !! Values to frame.
      character(len=*), intent(in), optional :: names(:) !! Optional element names.
      character(len=*), intent(in), optional :: name  !! Output label-column name.
      character(len=*), intent(in), optional :: value !! Output value-column name.
      type(tibble_type) :: table

      table = enframe_column(make_column(output_value_name(value), values), names, name)
   end function enframe_character

   function enframe_column(value_column, names, name) result(table)
      !! Implements common typed-vector framing behavior.
      type(tibble_column), intent(in) :: value_column !! Value column.
      character(len=*), intent(in), optional :: names(:) !! Optional labels.
      character(len=*), intent(in), optional :: name     !! Label-column name.
      type(tibble_type) :: table
      type(tibble_column) :: columns(2)
      character(len=:), allocatable :: label_name
      integer, allocatable :: sequence(:)
      integer :: i

      label_name = "name"
      if (present(name)) label_name = name
      if (present(names)) then
         if (size(names) /= value_column%size()) then
            error stop "enframe: names and values have different sizes"
         end if
         columns(1) = make_column(label_name, names)
      else
         allocate (sequence(value_column%size()))
         do i = 1, size(sequence)
            sequence(i) = i
         end do
         columns(1) = make_column(label_name, sequence)
      end if
      columns(2) = value_column
      table = new_tibble(columns)
   end function enframe_column

   pure function output_value_name(value) result(name)
      !! Selects the requested value-column name or the tibble default.
      character(len=*), intent(in), optional :: value !! Optional name.
      character(len=:), allocatable :: name

      name = "value"
      if (present(value)) name = value
   end function output_value_name

   function deframe_integer(table, names) result(values)
      !! Converts a one- or two-column tibble whose value column is integer to a vector.
      type(tibble_type), intent(in) :: table !! One- or two-column table.
      character(len=:), allocatable, intent(out), optional :: names(:) !! Converted labels.
      integer, allocatable :: values(:)

      call check_deframe(table, tibble_integer, names)
      values = table%columns(table%ncol())%integer_values
   end function deframe_integer

   function deframe_real(table, names) result(values)
      !! Converts a one- or two-column tibble whose value column is real to a vector.
      type(tibble_type), intent(in) :: table !! One- or two-column table.
      character(len=:), allocatable, intent(out), optional :: names(:) !! Converted labels.
      real(dp), allocatable :: values(:)

      call check_deframe(table, tibble_real, names)
      values = table%columns(table%ncol())%real_values
   end function deframe_real

   function deframe_logical(table, names) result(values)
      !! Converts a one- or two-column tibble whose value column is logical to a vector.
      type(tibble_type), intent(in) :: table !! One- or two-column table.
      character(len=:), allocatable, intent(out), optional :: names(:) !! Converted labels.
      logical, allocatable :: values(:)

      call check_deframe(table, tibble_logical, names)
      values = table%columns(table%ncol())%logical_values
   end function deframe_logical

   function deframe_character(table, names) result(values)
      !! Converts a one- or two-column tibble whose value column is character to a vector.
      type(tibble_type), intent(in) :: table !! One- or two-column table.
      character(len=:), allocatable, intent(out), optional :: names(:) !! Converted labels.
      character(len=:), allocatable :: values(:)

      call check_deframe(table, tibble_character, names)
      values = table%columns(table%ncol())%character_values
   end function deframe_character

   subroutine check_deframe(table, value_type, names)
      !! Validates deframe structure and returns character labels when present.
      type(tibble_type), intent(in) :: table !! One- or two-column table.
      integer, intent(in) :: value_type      !! Required type of final column.
      character(len=:), allocatable, intent(out), optional :: names(:) !! Converted labels.
      integer :: i

      if (table%ncol() < 1 .or. table%ncol() > 2) then
         error stop "deframe: table must have one or two columns"
      end if
      if (table%columns(table%ncol())%type_code /= value_type) then
         error stop "deframe: unexpected value-column type"
      end if
      if (present(names)) then
         if (table%ncol() == 1) then
            allocate (character(len=0) :: names(0))
         else
            select case (table%columns(1)%type_code)
            case (tibble_character)
               names = table%columns(1)%character_values
            case default
               allocate (character(len=64) :: names(table%nrow()))
               do i = 1, table%nrow()
                  names(i) = cell_text(table%columns(1), i)
               end do
            end select
         end if
      end if
   end subroutine check_deframe

   subroutine print_tibble(table, n)
      !! Prints a compact type-aware preview resembling the R tibble display.
      type(tibble_type), intent(in) :: table !! Table to print.
      integer, intent(in), optional :: n     !! Maximum rows to display; default ten.
      integer :: i, j, rows_to_print

      rows_to_print = min(10, table%nrow())
      if (present(n)) rows_to_print = min(table%nrow(), max(0, n))
      write (*, '(a,i0,a,i0)') "# A tibble: ", table%nrow(), " x ", table%ncol()
      if (table%ncol() == 0) return
      write (*, '(*(a,1x))') (table%columns(j)%name, j=1, table%ncol())
      write (*, '(*(a,1x))') ("<" // table%columns(j)%type_name() // ">", j=1, table%ncol())
      do i = 1, rows_to_print
         write (*, '(*(a,1x))') (cell_text(table%columns(j), i), j=1, table%ncol())
      end do
      if (rows_to_print < table%nrow()) then
         write (*, '(a,i0,a)') &
            "# ... with ", table%nrow() - rows_to_print, " more rows"
      end if
   end subroutine print_tibble

   subroutine glimpse(table, n)
      !! Prints table dimensions and a transposed compact preview of each column.
      type(tibble_type), intent(in) :: table !! Table to inspect.
      integer, intent(in), optional :: n     !! Maximum values per column; default five.
      integer :: i, j, values_to_print

      values_to_print = 5
      if (present(n)) values_to_print = max(0, n)
      write (*, '(a,i0,a,i0,a)') "Rows: ", table%nrow(), "; Columns: ", table%ncol(), ""
      do j = 1, table%ncol()
         write (*, '(a,a,a)', advance='no') &
            "$ ", table%columns(j)%name, " <" // table%columns(j)%type_name() // "> "
         do i = 1, min(values_to_print, table%nrow())
            if (i > 1) write (*, '(a)', advance='no') ", "
            write (*, '(a)', advance='no') cell_text(table%columns(j), i)
         end do
         if (values_to_print < table%nrow()) write (*, '(a)', advance='no') ", ..."
         write (*, *)
      end do
   end subroutine glimpse

   function cell_text(column, row) result(text)
      !! Formats one table cell, respecting its missing-value mask and type.
      type(tibble_column), intent(in) :: column !! Column containing the cell.
      integer, intent(in) :: row                !! One-based row index.
      character(len=:), allocatable :: text
      character(len=128) :: buffer

      if (column%missing(row)) then
         text = "NA"
         return
      end if
      select case (column%type_code)
      case (tibble_integer)
         write (buffer, '(i0)') column%integer_values(row)
         text = trim(buffer)
      case (tibble_real)
         write (buffer, '(g0.8)') column%real_values(row)
         text = trim(buffer)
      case (tibble_logical)
         if (column%logical_values(row)) then
            text = "TRUE"
         else
            text = "FALSE"
         end if
      case (tibble_character)
         text = trim(column%character_values(row))
      case default
         text = "<?>"
      end select
   end function cell_text

end module tibble_operations
