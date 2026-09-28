! SPDX-License-Identifier: GPL-3.0-or-later
! Sparse matrix support translated from sanic's dgCMatrix-facing operations.
module sanic_sparse
   use sanic_kinds, only: dp
   implicit none
   private

   type, public :: sparse_matrix
      integer :: nrow = 0
      integer :: ncol = 0
      integer, allocatable :: row(:)
      integer, allocatable :: col(:)
      real(dp), allocatable :: value(:)
   end type sparse_matrix

   public :: sparsify
   public :: sparse_from_triplets
   public :: sparse_to_dense
   public :: sparse_matvec
   public :: sparse_transpose_matvec
   public :: sparse_matmat
   public :: sparse_is_symmetric
   public :: canonical_sparse_entries

contains

   pure function sparsify(x, tol) result(a)
      !! Convert a dense matrix to package-local coordinate sparse storage.
      real(dp), intent(in) :: x(:, :) !! Dense matrix to convert; its full shape is retained.
      real(dp), intent(in), optional :: tol !! Entries with absolute value <= tol are omitted; default is zero.
      type(sparse_matrix) :: a
      real(dp) :: threshold
      integer :: i, j, k, nnz

      threshold = 0.0_dp
      if (present(tol)) threshold = max(0.0_dp, tol)

      a%nrow = size(x, 1)
      a%ncol = size(x, 2)
      nnz = count(abs(x) > threshold)
      allocate(a%row(nnz), a%col(nnz), a%value(nnz))

      k = 0
      do j = 1, a%ncol
         do i = 1, a%nrow
            if (abs(x(i, j)) > threshold) then
               k = k + 1
               a%row(k) = i
               a%col(k) = j
               a%value(k) = x(i, j)
            end if
         end do
      end do
   end function sparsify

   pure function sparse_from_triplets(nrow, ncol, row, col, value) result(a)
      !! Construct coordinate sparse storage from one-based triplets.
      integer, intent(in) :: nrow !! Number of matrix rows; must be nonnegative.
      integer, intent(in) :: ncol !! Number of matrix columns; must be nonnegative.
      integer, intent(in) :: row(:) !! One-based row index for each stored entry.
      integer, intent(in) :: col(:) !! One-based column index for each stored entry.
      real(dp), intent(in) :: value(:) !! Stored values corresponding to row and col.
      type(sparse_matrix) :: a
      integer :: k, nvalid

      a%nrow = max(0, nrow)
      a%ncol = max(0, ncol)
      if (size(row) /= size(col) .or. size(row) /= size(value)) then
         allocate(a%row(0), a%col(0), a%value(0))
         return
      end if

      nvalid = 0
      do k = 1, size(value)
         if (row(k) >= 1 .and. row(k) <= a%nrow .and. &
             col(k) >= 1 .and. col(k) <= a%ncol .and. abs(value(k)) > 0.0_dp) then
            nvalid = nvalid + 1
         end if
      end do

      allocate(a%row(nvalid), a%col(nvalid), a%value(nvalid))
      nvalid = 0
      do k = 1, size(value)
         if (row(k) >= 1 .and. row(k) <= a%nrow .and. &
             col(k) >= 1 .and. col(k) <= a%ncol .and. abs(value(k)) > 0.0_dp) then
            nvalid = nvalid + 1
            a%row(nvalid) = row(k)
            a%col(nvalid) = col(k)
            a%value(nvalid) = value(k)
         end if
      end do
   end function sparse_from_triplets

   pure function sparse_to_dense(a) result(x)
      !! Expand coordinate sparse storage to a dense matrix, summing duplicates.
      type(sparse_matrix), intent(in) :: a !! Sparse matrix to expand; duplicate coordinates are summed.
      real(dp), allocatable :: x(:, :)
      integer :: k

      allocate(x(a%nrow, a%ncol))
      x = 0.0_dp
      if (.not. allocated(a%value)) return
      do k = 1, size(a%value)
         if (a%row(k) >= 1 .and. a%row(k) <= a%nrow .and. &
             a%col(k) >= 1 .and. a%col(k) <= a%ncol) then
            x(a%row(k), a%col(k)) = x(a%row(k), a%col(k)) + a%value(k)
         end if
      end do
   end function sparse_to_dense

   pure subroutine sparse_matvec(a, x, y, info)
      !! Compute y = A*x directly from sparse coordinate storage.
      type(sparse_matrix), intent(in) :: a !! Sparse coefficient matrix with one-based stored coordinates.
      real(dp), intent(in) :: x(:) !! Input vector of length a%ncol.
      real(dp), allocatable, intent(out) :: y(:) !! Product vector of length a%nrow.
      integer, intent(out) :: info !! Status: zero on success, one for incompatible dimensions.
      integer :: k

      info = 0
      if (size(x) /= a%ncol) then
         allocate(y(0))
         info = 1
         return
      end if

      allocate(y(a%nrow))
      y = 0.0_dp
      if (.not. allocated(a%value)) return
      do k = 1, size(a%value)
         y(a%row(k)) = y(a%row(k)) + a%value(k)*x(a%col(k))
      end do
   end subroutine sparse_matvec

   pure subroutine sparse_transpose_matvec(a, x, y, info)
      !! Compute y = transpose(A)*x directly from sparse coordinate storage.
      type(sparse_matrix), intent(in) :: a !! Sparse coefficient matrix with one-based stored coordinates.
      real(dp), intent(in) :: x(:) !! Input vector of length a%nrow.
      real(dp), allocatable, intent(out) :: y(:) !! Product vector of length a%ncol.
      integer, intent(out) :: info !! Status: zero on success, one for incompatible dimensions.
      integer :: k

      info = 0
      if (size(x) /= a%nrow) then
         allocate(y(0))
         info = 1
         return
      end if

      allocate(y(a%ncol))
      y = 0.0_dp
      if (.not. allocated(a%value)) return
      do k = 1, size(a%value)
         y(a%col(k)) = y(a%col(k)) + a%value(k)*x(a%row(k))
      end do
   end subroutine sparse_transpose_matvec

   pure subroutine sparse_matmat(a, x, y, info)
      !! Compute Y = A*X directly from sparse coordinate storage.
      type(sparse_matrix), intent(in) :: a !! Sparse coefficient matrix with one-based stored coordinates.
      real(dp), intent(in) :: x(:, :) !! Dense right factor with a%ncol rows.
      real(dp), allocatable, intent(out) :: y(:, :) !! Dense product with a%nrow rows.
      integer, intent(out) :: info !! Status: zero on success, one for incompatible dimensions.
      integer :: j, status
      real(dp), allocatable :: column(:)

      info = 0
      if (size(x, 1) /= a%ncol) then
         allocate(y(0, 0))
         info = 1
         return
      end if

      allocate(y(a%nrow, size(x, 2)))
      do j = 1, size(x, 2)
         call sparse_matvec(a, x(:, j), column, status)
         if (status /= 0) then
            info = status
            return
         end if
         y(:, j) = column
      end do
   end subroutine sparse_matmat

   pure logical function sparse_is_symmetric(a, tol) result(ok)
      !! Test sparse symmetry from canonical duplicate-summed coordinates without dense storage or repeated full scans.
      type(sparse_matrix), intent(in) :: a !! Sparse matrix to test; duplicate coordinates are summed before comparison.
      real(dp), intent(in), optional :: tol !! Relative Frobenius tolerance; zero or absent uses machine epsilon.
      integer, allocatable :: col_index(:), row_index(:)
      real(dp), allocatable :: values(:)
      real(dp) :: difference_squared, matrix_squared, threshold, transposed_value
      integer :: k, position

      if (a%nrow /= a%ncol) then
         ok = .false.
         return
      end if
      threshold = epsilon(1.0_dp)
      if (present(tol)) then
         if (tol > 0.0_dp) threshold = tol
      end if
      call canonical_sparse_entries(a, row_index, col_index, values)
      if (size(values) == 0) then
         ok = .true.
         return
      end if
      matrix_squared = sum(values*values)
      if (matrix_squared <= tiny(1.0_dp)) then
         ok = .true.
         return
      end if
      difference_squared = 0.0_dp
      do k = 1, size(values)
         position = coordinate_find(row_index, col_index, col_index(k), row_index(k))
         if (position > 0) then
            transposed_value = values(position)
            difference_squared = difference_squared + (values(k) - transposed_value)**2
         else
            difference_squared = difference_squared + values(k)*values(k)
            if (row_index(k) /= col_index(k)) difference_squared = difference_squared + values(k)*values(k)
         end if
      end do
      ok = sqrt(difference_squared) <= threshold*sqrt(matrix_squared)
   end function sparse_is_symmetric

   pure subroutine canonical_sparse_entries(a, row_index, col_index, values)
      !! Return valid coordinate entries sorted by row/column with duplicate coordinates summed.
      type(sparse_matrix), intent(in) :: a !! Coordinate sparse matrix whose stored entries are canonicalized.
      integer, allocatable, intent(out) :: row_index(:) !! Sorted one-based row coordinates after duplicate summation.
      integer, allocatable, intent(out) :: col_index(:) !! Sorted one-based column coordinates after duplicate summation.
      real(dp), allocatable, intent(out) :: values(:) !! Summed values aligned with row_index and col_index.
      integer, allocatable :: rows_work(:), cols_work(:), rows_final(:), cols_final(:)
      real(dp), allocatable :: values_work(:), values_final(:)
      integer :: count_value, k, kept, nstored

      nstored = 0
      if (allocated(a%value)) nstored = size(a%value)
      allocate(rows_work(nstored), cols_work(nstored), values_work(nstored))
      count_value = 0
      do k = 1, nstored
         if (a%row(k) < 1 .or. a%row(k) > a%nrow) cycle
         if (a%col(k) < 1 .or. a%col(k) > a%ncol) cycle
         count_value = count_value + 1
         rows_work(count_value) = a%row(k)
         cols_work(count_value) = a%col(k)
         values_work(count_value) = a%value(k)
      end do
      if (count_value == 0) then
         allocate(row_index(0), col_index(0), values(0))
         return
      end if
      call sort_coordinate_triplets(rows_work, cols_work, values_work, 1, count_value)

      allocate(rows_final(count_value), cols_final(count_value), values_final(count_value))
      kept = 0
      do k = 1, count_value
         if (kept > 0) then
            if (rows_work(k) == rows_final(kept) .and. cols_work(k) == cols_final(kept)) then
               values_final(kept) = values_final(kept) + values_work(k)
               cycle
            end if
         end if
         kept = kept + 1
         rows_final(kept) = rows_work(k)
         cols_final(kept) = cols_work(k)
         values_final(kept) = values_work(k)
      end do
      row_index = rows_final(1:kept)
      col_index = cols_final(1:kept)
      values = values_final(1:kept)
   end subroutine canonical_sparse_entries

   recursive pure subroutine sort_coordinate_triplets(rows, cols, values, left, right)
      !! Quicksort coordinate triplets lexicographically by row then column.
      integer, intent(inout) :: rows(:) !! Row coordinates reordered in place.
      integer, intent(inout) :: cols(:) !! Column coordinates reordered in place with rows.
      real(dp), intent(inout) :: values(:) !! Numerical entries reordered in place with coordinates.
      integer, intent(in) :: left !! First array position in the partition.
      integer, intent(in) :: right !! Last array position in the partition.
      integer :: i, j, pivot_col, pivot_row, temp_integer
      real(dp) :: temp_value

      if (left >= right) return
      pivot_row = rows((left + right)/2)
      pivot_col = cols((left + right)/2)
      i = left
      j = right
      do
         do while (coordinate_less(rows(i), cols(i), pivot_row, pivot_col))
            i = i + 1
         end do
         do while (coordinate_less(pivot_row, pivot_col, rows(j), cols(j)))
            j = j - 1
         end do
         if (i > j) exit
         temp_integer = rows(i)
         rows(i) = rows(j)
         rows(j) = temp_integer
         temp_integer = cols(i)
         cols(i) = cols(j)
         cols(j) = temp_integer
         temp_value = values(i)
         values(i) = values(j)
         values(j) = temp_value
         i = i + 1
         j = j - 1
         if (i > j) exit
      end do
      if (left < j) call sort_coordinate_triplets(rows, cols, values, left, j)
      if (i < right) call sort_coordinate_triplets(rows, cols, values, i, right)
   end subroutine sort_coordinate_triplets

   pure logical function coordinate_less(row_a, col_a, row_b, col_b) result(less)
      !! Compare two one-based coordinates lexicographically.
      integer, intent(in) :: row_a !! Row coordinate of the left comparison operand.
      integer, intent(in) :: col_a !! Column coordinate of the left comparison operand.
      integer, intent(in) :: row_b !! Row coordinate of the right comparison operand.
      integer, intent(in) :: col_b !! Column coordinate of the right comparison operand.

      less = row_a < row_b .or. (row_a == row_b .and. col_a < col_b)
   end function coordinate_less

   pure integer function coordinate_find(rows, cols, target_row, target_col) result(position)
      !! Binary-search one coordinate in canonical row-major coordinate arrays.
      integer, intent(in) :: rows(:) !! Sorted row coordinates.
      integer, intent(in) :: cols(:) !! Sorted column coordinates aligned with rows.
      integer, intent(in) :: target_row !! One-based row coordinate sought.
      integer, intent(in) :: target_col !! One-based column coordinate sought.
      integer :: left, middle, right

      position = 0
      left = 1
      right = size(rows)
      do while (left <= right)
         middle = (left + right)/2
         if (rows(middle) == target_row .and. cols(middle) == target_col) then
            position = middle
            return
         end if
         if (coordinate_less(rows(middle), cols(middle), target_row, target_col)) then
            left = middle + 1
         else
            right = middle - 1
         end if
      end do
   end function coordinate_find

end module sanic_sparse
