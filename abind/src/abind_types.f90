module abind_types
   use iso_fortran_env, only : real64
   implicit none
   private

   integer, parameter, public :: dp = real64

   type, public :: string_vector
      character(len=:), allocatable :: values(:)
   end type string_vector

   type, public :: int_vector
      integer, allocatable :: values(:)
   end type int_vector

   type, public :: array_value
      real(dp), allocatable :: data(:)
      integer, allocatable :: dims(:)
      type(string_vector), allocatable :: dimnames(:)
      character(len=:), allocatable :: dimname_names(:)
   end type array_value

   public :: init_array
   public :: set_dimnames
   public :: set_dimname_names
   public :: clear_dimnames
   public :: array_rank
   public :: array_size
   public :: has_dimnames
   public :: linear_index
   public :: unravel_index
   public :: make_int_vector

contains

   pure integer function array_rank(x) result(rank_value)
      type(array_value), intent(in) :: x !! Array whose number of dimensions is requested.

      if (allocated(x%dims)) then
         rank_value = size(x%dims)
      else
         rank_value = 0
      end if
   end function array_rank

   pure integer function array_size(x) result(n)
      type(array_value), intent(in) :: x !! Array whose stored element count is requested.

      if (allocated(x%data)) then
         n = size(x%data)
      else
         n = 0
      end if
   end function array_size

   subroutine init_array(x, data, dims, status)
      type(array_value), intent(out) :: x !! Array object initialized by this routine.
      real(dp), intent(in) :: data(:) !! Flat column-major element values to store.
      integer, intent(in) :: dims(:) !! Dimension lengths; their product must equal SIZE(data).
      integer, intent(out), optional :: status !! Zero on success; nonzero when dimensions are invalid.
      integer :: expected

      if (present(status)) status = 0
      if (any(dims < 0)) then
         if (present(status)) status = 1
         return
      end if
      expected = product_or_one(dims)
      if (expected /= size(data)) then
         if (present(status)) status = 2
         return
      end if
      x%data = data
      x%dims = dims
      allocate(x%dimnames(size(dims)))
      allocate(character(len=1) :: x%dimname_names(size(dims)))
      x%dimname_names = ''
   end subroutine init_array

   subroutine clear_dimnames(x)
      type(array_value), intent(inout) :: x !! Array whose dimension labels are cleared while retaining shape/data.
      integer :: i

      if (allocated(x%dimnames)) then
         do i = 1, size(x%dimnames)
            if (allocated(x%dimnames(i)%values)) deallocate(x%dimnames(i)%values)
         end do
      end if
      if (allocated(x%dimname_names)) x%dimname_names = ''
   end subroutine clear_dimnames

   subroutine set_dimnames(x, dim_number, names, status)
      type(array_value), intent(inout) :: x !! Array receiving labels for one dimension.
      integer, intent(in) :: dim_number !! One-based dimension number whose element labels are set.
      character(len=*), intent(in) :: names(:) !! Labels; length must equal the selected dimension length.
      integer, intent(out), optional :: status !! Zero on success; nonzero for an invalid dimension or label count.
      integer :: width

      if (present(status)) status = 0
      if (.not. allocated(x%dims)) then
         if (present(status)) status = 1
         return
      end if
      if (dim_number < 1 .or. dim_number > size(x%dims)) then
         if (present(status)) status = 2
         return
      end if
      if (size(names) /= x%dims(dim_number)) then
         if (present(status)) status = 3
         return
      end if
      if (.not. allocated(x%dimnames)) allocate(x%dimnames(size(x%dims)))
      width = max(1, len(names))
      if (allocated(x%dimnames(dim_number)%values)) deallocate(x%dimnames(dim_number)%values)
      allocate(character(len=width) :: x%dimnames(dim_number)%values(size(names)))
      x%dimnames(dim_number)%values = names
   end subroutine set_dimnames

   subroutine set_dimname_names(x, names, status)
      type(array_value), intent(inout) :: x !! Array receiving names for its dimensions.
      character(len=*), intent(in) :: names(:) !! Dimension names, one per dimension.
      integer, intent(out), optional :: status !! Zero on success; nonzero when the count does not match array rank.
      integer :: width

      if (present(status)) status = 0
      if (.not. allocated(x%dims)) then
         if (present(status)) status = 1
         return
      end if
      if (size(names) /= size(x%dims)) then
         if (present(status)) status = 2
         return
      end if
      width = max(1, len(names))
      if (allocated(x%dimname_names)) deallocate(x%dimname_names)
      allocate(character(len=width) :: x%dimname_names(size(names)))
      x%dimname_names = names
   end subroutine set_dimname_names

   pure logical function has_dimnames(x, dim_number) result(has_names)
      type(array_value), intent(in) :: x !! Array queried for labels on one dimension.
      integer, intent(in) :: dim_number !! One-based dimension number to query.

      has_names = .false.
      if (.not. allocated(x%dimnames)) return
      if (dim_number < 1 .or. dim_number > size(x%dimnames)) return
      if (.not. allocated(x%dimnames(dim_number)%values)) return
      has_names = size(x%dimnames(dim_number)%values) == x%dims(dim_number)
   end function has_dimnames

   pure integer function linear_index(coords, dims) result(idx)
      integer, intent(in) :: coords(:) !! One-based multidimensional coordinates.
      integer, intent(in) :: dims(:) !! Dimension lengths defining column-major strides.
      integer :: i
      integer :: stride

      if (size(coords) == 0) then
         idx = 1
         return
      end if
      idx = coords(1)
      stride = dims(1)
      do i = 2, size(dims)
         idx = idx + (coords(i) - 1) * stride
         stride = stride * dims(i)
      end do
   end function linear_index

   pure subroutine unravel_index(index_value, dims, coords)
      integer, intent(in) :: index_value !! One-based flat column-major element index.
      integer, intent(in) :: dims(:) !! Dimension lengths defining the array shape.
      integer, intent(out) :: coords(:) !! One-based coordinates corresponding to INDEX_VALUE.
      integer :: i
      integer :: q

      if (size(dims) == 0) return
      q = index_value - 1
      do i = 1, size(dims)
         if (dims(i) > 0) then
            coords(i) = mod(q, dims(i)) + 1
            q = q / dims(i)
         else
            coords(i) = 1
         end if
      end do
   end subroutine unravel_index

   pure function make_int_vector(values) result(v)
      integer, intent(in) :: values(:) !! Integer values copied into a selector vector.
      type(int_vector) :: v

      allocate(v%values(size(values)))
      v%values = values
   end function make_int_vector

   pure integer function product_or_one(dims) result(value)
      integer, intent(in) :: dims(:) !! Dimension vector whose product is needed; rank-zero gives one.

      if (size(dims) == 0) then
         value = 1
      else
         value = product(dims)
      end if
   end function product_or_one

end module abind_types
