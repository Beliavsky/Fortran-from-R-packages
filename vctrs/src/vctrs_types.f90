! SPDX-License-Identifier: MIT
module vctrs_types
   use, intrinsic :: iso_fortran_env, only: real64
   implicit none
   private

   integer, parameter, public :: dp = real64, vctrs_integer = 1, vctrs_real = 2, &
      vctrs_logical = 3, vctrs_character = 4

   type, public :: vctr_type
      !! Stores a homogeneous typed vector and an explicit missing-value mask.
      character(len=:), allocatable :: name
      integer :: type_code = 0
      integer, allocatable :: integer_values(:)
      real(dp), allocatable :: real_values(:)
      logical, allocatable :: logical_values(:)
      character(len=:), allocatable :: character_values(:)
      logical, allocatable :: missing(:)
   contains
      procedure :: size => vctr_size
      procedure :: type_name => vctr_type_name
   end type vctr_type

   type, public :: vctrs_data_frame
      !! Stores a rectangular collection of named vectors.
      type(vctr_type), allocatable :: columns(:)
      integer :: number_of_rows = 0
   contains
      procedure :: nrow => data_frame_number_of_rows
      procedure :: ncol => data_frame_number_of_columns
      procedure :: column_index => data_frame_column_index
      procedure :: has_name => data_frame_has_name
   end type vctrs_data_frame

   public :: new_vctr

   interface new_vctr
      module procedure new_integer_vctr, new_integer_scalar_vctr
      module procedure new_real_vctr, new_real_scalar_vctr
      module procedure new_logical_vctr, new_logical_scalar_vctr
      module procedure new_character_vctr, new_character_scalar_vctr
   end interface new_vctr

contains

   pure function new_integer_vctr(name, values, missing) result(vector)
      !! Constructs an integer vector without coercing its values.
      character(len=*), intent(in) :: name        !! Optional semantic name.
      integer, intent(in) :: values(:)            !! Vector values.
      logical, intent(in), optional :: missing(:) !! Missing-value mask.
      type(vctr_type) :: vector

      vector%name = trim(name)
      vector%type_code = vctrs_integer
      vector%integer_values = values
      call set_missing_mask(vector, size(values), missing)
   end function new_integer_vctr

   pure elemental function new_integer_scalar_vctr(name, value, n, missing) result(vector)
      !! Constructs an integer vector by recycling a scalar to a requested size.
      character(len=*), intent(in) :: name !! Optional semantic name.
      integer, intent(in) :: value         !! Scalar value to recycle.
      integer, intent(in) :: n             !! Nonnegative output size.
      logical, intent(in), optional :: missing !! Whether every output value is missing.
      type(vctr_type) :: vector

      if (n < 0) error stop "new_vctr: vector size must be nonnegative"
      vector%name = trim(name)
      vector%type_code = vctrs_integer
      allocate (vector%integer_values(n), source=value)
      allocate (vector%missing(n), source=.false.)
      if (present(missing)) vector%missing = missing
   end function new_integer_scalar_vctr

   pure function new_real_vctr(name, values, missing) result(vector)
      !! Constructs a double-precision real vector without coercing its values.
      character(len=*), intent(in) :: name        !! Optional semantic name.
      real(dp), intent(in) :: values(:)           !! Vector values.
      logical, intent(in), optional :: missing(:) !! Missing-value mask.
      type(vctr_type) :: vector

      vector%name = trim(name)
      vector%type_code = vctrs_real
      vector%real_values = values
      call set_missing_mask(vector, size(values), missing)
   end function new_real_vctr

   pure elemental function new_real_scalar_vctr(name, value, n, missing) result(vector)
      !! Constructs a real vector by recycling a scalar to a requested size.
      character(len=*), intent(in) :: name !! Optional semantic name.
      real(dp), intent(in) :: value        !! Scalar value to recycle.
      integer, intent(in) :: n             !! Nonnegative output size.
      logical, intent(in), optional :: missing !! Whether every output value is missing.
      type(vctr_type) :: vector

      if (n < 0) error stop "new_vctr: vector size must be nonnegative"
      vector%name = trim(name)
      vector%type_code = vctrs_real
      allocate (vector%real_values(n), source=value)
      allocate (vector%missing(n), source=.false.)
      if (present(missing)) vector%missing = missing
   end function new_real_scalar_vctr

   pure function new_logical_vctr(name, values, missing) result(vector)
      !! Constructs a logical vector without coercing its values.
      character(len=*), intent(in) :: name        !! Optional semantic name.
      logical, intent(in) :: values(:)            !! Vector values.
      logical, intent(in), optional :: missing(:) !! Missing-value mask.
      type(vctr_type) :: vector

      vector%name = trim(name)
      vector%type_code = vctrs_logical
      vector%logical_values = values
      call set_missing_mask(vector, size(values), missing)
   end function new_logical_vctr

   pure elemental function new_logical_scalar_vctr(name, value, n, missing) result(vector)
      !! Constructs a logical vector by recycling a scalar to a requested size.
      character(len=*), intent(in) :: name !! Optional semantic name.
      logical, intent(in) :: value         !! Scalar value to recycle.
      integer, intent(in) :: n             !! Nonnegative output size.
      logical, intent(in), optional :: missing !! Whether every output value is missing.
      type(vctr_type) :: vector

      if (n < 0) error stop "new_vctr: vector size must be nonnegative"
      vector%name = trim(name)
      vector%type_code = vctrs_logical
      allocate (vector%logical_values(n), source=value)
      allocate (vector%missing(n), source=.false.)
      if (present(missing)) vector%missing = missing
   end function new_logical_scalar_vctr

   pure function new_character_vctr(name, values, missing) result(vector)
      !! Constructs a character vector without coercing its values.
      character(len=*), intent(in) :: name        !! Optional semantic name.
      character(len=*), intent(in) :: values(:)   !! Vector values.
      logical, intent(in), optional :: missing(:) !! Missing-value mask.
      type(vctr_type) :: vector

      vector%name = trim(name)
      vector%type_code = vctrs_character
      vector%character_values = values
      call set_missing_mask(vector, size(values), missing)
   end function new_character_vctr

   pure elemental function new_character_scalar_vctr(name, value, n, missing) result(vector)
      !! Constructs a character vector by recycling a scalar to a requested size.
      character(len=*), intent(in) :: name  !! Optional semantic name.
      character(len=*), intent(in) :: value !! Scalar value to recycle.
      integer, intent(in) :: n              !! Nonnegative output size.
      logical, intent(in), optional :: missing !! Whether every output value is missing.
      type(vctr_type) :: vector

      if (n < 0) error stop "new_vctr: vector size must be nonnegative"
      vector%name = trim(name)
      vector%type_code = vctrs_character
      allocate (character(len=len(value)) :: vector%character_values(n))
      vector%character_values = value
      allocate (vector%missing(n), source=.false.)
      if (present(missing)) vector%missing = missing
   end function new_character_scalar_vctr

   pure subroutine set_missing_mask(vector, n, missing)
      !! Validates and stores a missing-value mask.
      type(vctr_type), intent(inout) :: vector !! Vector receiving the mask.
      integer, intent(in) :: n                 !! Required mask size.
      logical, intent(in), optional :: missing(:) !! Caller-supplied mask.

      allocate (vector%missing(n), source=.false.)
      if (present(missing)) then
         if (size(missing) /= n) error stop "new_vctr: missing mask has the wrong size"
         vector%missing = missing
      end if
   end subroutine set_missing_mask

   pure elemental integer function vctr_size(self) result(n)
      !! Returns the number of elements in a vector.
      class(vctr_type), intent(in) :: self !! Vector to inspect.

      n = 0
      select case (self%type_code)
      case (vctrs_integer)
         if (allocated(self%integer_values)) n = size(self%integer_values)
      case (vctrs_real)
         if (allocated(self%real_values)) n = size(self%real_values)
      case (vctrs_logical)
         if (allocated(self%logical_values)) n = size(self%logical_values)
      case (vctrs_character)
         if (allocated(self%character_values)) n = size(self%character_values)
      end select
   end function vctr_size

   pure elemental function vctr_type_name(self) result(name)
      !! Returns a compact type abbreviation.
      class(vctr_type), intent(in) :: self !! Vector to inspect.
      character(len=3) :: name

      select case (self%type_code)
      case (vctrs_integer)
         name = "int"
      case (vctrs_real)
         name = "dbl"
      case (vctrs_logical)
         name = "lgl"
      case (vctrs_character)
         name = "chr"
      case default
         name = "???"
      end select
   end function vctr_type_name

   pure elemental integer function data_frame_number_of_rows(self) result(n)
      !! Returns the data-frame row count, including for zero-column frames.
      class(vctrs_data_frame), intent(in) :: self !! Frame to inspect.

      n = self%number_of_rows
   end function data_frame_number_of_rows

   pure elemental integer function data_frame_number_of_columns(self) result(n)
      !! Returns the data-frame column count.
      class(vctrs_data_frame), intent(in) :: self !! Frame to inspect.

      n = 0
      if (allocated(self%columns)) n = size(self%columns)
   end function data_frame_number_of_columns

   pure elemental integer function data_frame_column_index(self, name) result(index)
      !! Finds an exact column name without partial matching.
      class(vctrs_data_frame), intent(in) :: self !! Frame to search.
      character(len=*), intent(in) :: name        !! Exact column name.
      integer :: j

      index = 0
      do j = 1, self%ncol()
         if (self%columns(j)%name == name) then
            index = j
            return
         end if
      end do
   end function data_frame_column_index

   pure elemental logical function data_frame_has_name(self, name) result(found)
      !! Reports whether an exact column name exists.
      class(vctrs_data_frame), intent(in) :: self !! Frame to search.
      character(len=*), intent(in) :: name        !! Exact column name.

      found = self%column_index(name) > 0
   end function data_frame_has_name

end module vctrs_types
