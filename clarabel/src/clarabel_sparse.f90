! SPDX-License-Identifier: Apache-2.0
module clarabel_sparse
   use, intrinsic :: iso_c_binding, only : c_size_t
   use clarabel_kinds, only : dp
   implicit none
   private

   type, public :: csc_matrix
      integer(c_size_t) :: nrows = 0_c_size_t
      integer(c_size_t) :: ncols = 0_c_size_t
      integer(c_size_t), allocatable :: colptr(:)
      integer(c_size_t), allocatable :: rowind(:)
      real(dp), allocatable :: values(:)
   contains
      procedure :: nnz => csc_nnz
      procedure :: validate => csc_validate
      procedure :: to_dense => csc_to_dense
   end type csc_matrix

   public :: csc_from_dense, csc_from_symmetric_upper, csc_from_arrays, csc_from_triplets, csc_empty

contains

   function csc_empty(nrows, ncols) result(a)
      integer, intent(in) :: nrows, ncols
      type(csc_matrix) :: a
      a%nrows = int(nrows, c_size_t)
      a%ncols = int(ncols, c_size_t)
      allocate(a%colptr(ncols + 1), a%rowind(0), a%values(0))
      a%colptr = 0_c_size_t
   end function csc_empty

   function csc_from_dense(x, tolerance) result(a)
      real(dp), intent(in) :: x(:, :)
      real(dp), intent(in), optional :: tolerance
      type(csc_matrix) :: a
      real(dp) :: tol
      integer :: i, j, k, nz

      tol = 0.0_dp
      if (present(tolerance)) tol = max(0.0_dp, tolerance)
      nz = count(abs(x) > tol)
      a%nrows = int(size(x, 1), c_size_t)
      a%ncols = int(size(x, 2), c_size_t)
      allocate(a%colptr(size(x, 2) + 1), a%rowind(nz), a%values(nz))
      k = 0
      a%colptr(1) = 0_c_size_t
      do j = 1, size(x, 2)
         do i = 1, size(x, 1)
            if (abs(x(i, j)) > tol) then
               k = k + 1
               a%rowind(k) = int(i - 1, c_size_t)
               a%values(k) = x(i, j)
            end if
         end do
         a%colptr(j + 1) = int(k, c_size_t)
      end do
   end function csc_from_dense

   function csc_from_symmetric_upper(x, tolerance) result(a)
      real(dp), intent(in) :: x(:, :)
      real(dp), intent(in), optional :: tolerance
      type(csc_matrix) :: a
      real(dp) :: tol, scale
      integer :: i, j, k, nz, n

      n = size(x, 1)
      if (size(x, 2) /= n) error stop "csc_from_symmetric_upper: matrix must be square"
      scale = max(1.0_dp, maxval(abs(x)))
      if (maxval(abs(x - transpose(x))) > 100.0_dp * epsilon(1.0_dp) * scale) &
         error stop "csc_from_symmetric_upper: matrix is not symmetric"
      tol = 0.0_dp
      if (present(tolerance)) tol = max(0.0_dp, tolerance)
      nz = 0
      do j = 1, n
         do i = 1, j
            if (abs(x(i, j)) > tol) nz = nz + 1
         end do
      end do
      a%nrows = int(n, c_size_t)
      a%ncols = int(n, c_size_t)
      allocate(a%colptr(n + 1), a%rowind(nz), a%values(nz))
      k = 0
      a%colptr(1) = 0_c_size_t
      do j = 1, n
         do i = 1, j
            if (abs(x(i, j)) > tol) then
               k = k + 1
               a%rowind(k) = int(i - 1, c_size_t)
               a%values(k) = x(i, j)
            end if
         end do
         a%colptr(j + 1) = int(k, c_size_t)
      end do
   end function csc_from_symmetric_upper


   function csc_from_arrays(nrows, ncols, colptr, rowind, values, index_base) result(a)
      integer, intent(in) :: nrows, ncols
      integer, intent(in) :: colptr(:), rowind(:)
      real(dp), intent(in) :: values(:)
      integer, intent(in), optional :: index_base
      type(csc_matrix) :: a
      integer :: base
      logical :: ok
      character(len=:), allocatable :: message

      base = 0
      if (present(index_base)) base = index_base
      if (base /= 0 .and. base /= 1) error stop "csc_from_arrays: index_base must be zero or one"
      if (nrows < 0 .or. ncols < 0) error stop "csc_from_arrays: dimensions cannot be negative"
      if (size(colptr) /= ncols + 1) error stop "csc_from_arrays: colptr has the wrong length"
      if (size(rowind) /= size(values)) error stop "csc_from_arrays: rowind and values lengths differ"
      a%nrows = int(nrows, c_size_t)
      a%ncols = int(ncols, c_size_t)
      allocate(a%colptr(size(colptr)), a%rowind(size(rowind)), a%values(size(values)))
      a%colptr = int(colptr - base, c_size_t)
      a%rowind = int(rowind - base, c_size_t)
      a%values = values
      call a%validate(ok, message)
      if (.not. ok) error stop "csc_from_arrays: " // message
   end function csc_from_arrays

   function csc_from_triplets(nrows, ncols, rows, cols, values, one_based, tolerance) result(a)
      integer, intent(in) :: nrows, ncols
      integer, intent(in) :: rows(:), cols(:)
      real(dp), intent(in) :: values(:)
      logical, intent(in), optional :: one_based
      real(dp), intent(in), optional :: tolerance
      type(csc_matrix) :: a
      integer, allocatable :: count_col(:), next(:), row_work(:), col_work(:)
      real(dp), allocatable :: value_work(:), value_unique(:)
      integer(c_size_t), allocatable :: row_unique(:)
      integer :: base, i, j, k, nkeep, pos, first, last, rtmp, nunique
      real(dp) :: tol, vtmp
      logical :: use_one_based

      if (nrows < 0 .or. ncols < 0) error stop "csc_from_triplets: dimensions cannot be negative"
      if (size(rows) /= size(cols) .or. size(rows) /= size(values)) &
         error stop "csc_from_triplets: triplet arrays have different lengths"
      use_one_based = .true.
      if (present(one_based)) use_one_based = one_based
      base = merge(1, 0, use_one_based)
      tol = 0.0_dp
      if (present(tolerance)) tol = max(0.0_dp, tolerance)

      nkeep = count(abs(values) > tol)
      allocate(count_col(max(0, ncols)), source=0)
      allocate(row_work(nkeep), col_work(nkeep), value_work(nkeep))
      k = 0
      do i = 1, size(values)
         if (abs(values(i)) <= tol) cycle
         if (rows(i) - base < 0 .or. rows(i) - base >= nrows) &
            error stop "csc_from_triplets: row index is outside the matrix"
         if (cols(i) - base < 0 .or. cols(i) - base >= ncols) &
            error stop "csc_from_triplets: column index is outside the matrix"
         k = k + 1
         row_work(k) = rows(i) - base
         col_work(k) = cols(i) - base
         value_work(k) = values(i)
         count_col(col_work(k) + 1) = count_col(col_work(k) + 1) + 1
      end do

      a%nrows = int(nrows, c_size_t)
      a%ncols = int(ncols, c_size_t)
      allocate(a%colptr(ncols + 1))
      a%colptr(1) = 0_c_size_t
      do j = 1, ncols
         a%colptr(j + 1) = a%colptr(j) + int(count_col(j), c_size_t)
      end do
      allocate(next(max(0, ncols)))
      do j = 1, ncols
         next(j) = int(a%colptr(j)) + 1
      end do
      allocate(row_unique(nkeep), value_unique(nkeep))
      do i = 1, nkeep
         j = col_work(i) + 1
         pos = next(j)
         row_unique(pos) = int(row_work(i), c_size_t)
         value_unique(pos) = value_work(i)
         next(j) = pos + 1
      end do

      ! Sort rows within each column.  Columns are usually short, so insertion
      ! sort avoids another global ordering allocation.
      do j = 1, ncols
         first = int(a%colptr(j)) + 1
         last = int(a%colptr(j + 1))
         do i = first + 1, last
            rtmp = int(row_unique(i))
            vtmp = value_unique(i)
            k = i - 1
            do while (k >= first)
               if (int(row_unique(k)) <= rtmp) exit
               row_unique(k + 1) = row_unique(k)
               value_unique(k + 1) = value_unique(k)
               k = k - 1
            end do
            row_unique(k + 1) = int(rtmp, c_size_t)
            value_unique(k + 1) = vtmp
         end do
      end do

      ! Aggregate duplicate entries and rebuild column pointers.
      allocate(a%rowind(nkeep), a%values(nkeep))
      nunique = 0
      do j = 1, ncols
         first = int(a%colptr(j)) + 1
         last = int(a%colptr(j + 1))
         a%colptr(j) = int(nunique, c_size_t)
         i = first
         do while (i <= last)
            rtmp = int(row_unique(i))
            vtmp = value_unique(i)
            i = i + 1
            do while (i <= last)
               if (int(row_unique(i)) /= rtmp) exit
               vtmp = vtmp + value_unique(i)
               i = i + 1
            end do
            if (abs(vtmp) > tol) then
               nunique = nunique + 1
               a%rowind(nunique) = int(rtmp, c_size_t)
               a%values(nunique) = vtmp
            end if
         end do
      end do
      a%colptr(ncols + 1) = int(nunique, c_size_t)
      if (nunique < nkeep) then
         a%rowind = a%rowind(:nunique)
         a%values = a%values(:nunique)
      end if
   end function csc_from_triplets

   pure integer function csc_nnz(self) result(n)
      class(csc_matrix), intent(in) :: self
      if (allocated(self%values)) then
         n = size(self%values)
      else
         n = 0
      end if
   end function csc_nnz

   subroutine csc_validate(self, ok, message)
      class(csc_matrix), intent(in) :: self
      logical, intent(out) :: ok
      character(len=:), allocatable, intent(out) :: message
      integer :: j, first, last

      ok = .false.
      message = ""
      if (.not. allocated(self%colptr) .or. .not. allocated(self%rowind) .or. &
          .not. allocated(self%values)) then
         message = "CSC arrays are not allocated"
         return
      end if
      if (size(self%colptr) /= int(self%ncols) + 1) then
         message = "colptr has the wrong length"
         return
      end if
      if (size(self%rowind) /= size(self%values)) then
         message = "rowind and values lengths differ"
         return
      end if
      if (self%colptr(1) /= 0_c_size_t .or. self%colptr(size(self%colptr)) /= int(size(self%values), c_size_t)) then
         message = "colptr endpoints are invalid"
         return
      end if
      if (any(self%colptr(2:) < self%colptr(:size(self%colptr)-1))) then
         message = "colptr is not nondecreasing"
         return
      end if
      if (any(self%rowind < 0_c_size_t) .or. any(self%rowind >= self%nrows)) then
         message = "row index is outside the matrix"
         return
      end if
      do j = 1, int(self%ncols)
         first = int(self%colptr(j)) + 1
         last = int(self%colptr(j + 1))
         if (last >= first + 1) then
            if (any(self%rowind(first + 1:last) <= self%rowind(first:last - 1))) then
               message = "row indices must be strictly increasing within each column"
               return
            end if
         end if
      end do
      ok = .true.
   end subroutine csc_validate

   function csc_to_dense(self, symmetric_upper) result(x)
      class(csc_matrix), intent(in) :: self
      logical, intent(in), optional :: symmetric_upper
      real(dp), allocatable :: x(:, :)
      logical :: sym
      integer :: i, j, k

      sym = .false.
      if (present(symmetric_upper)) sym = symmetric_upper
      allocate(x(int(self%nrows), int(self%ncols)), source=0.0_dp)
      do j = 1, int(self%ncols)
         do k = int(self%colptr(j)) + 1, int(self%colptr(j + 1))
            i = int(self%rowind(k)) + 1
            x(i, j) = self%values(k)
            if (sym .and. i /= j) x(j, i) = self%values(k)
         end do
      end do
   end function csc_to_dense

end module clarabel_sparse
