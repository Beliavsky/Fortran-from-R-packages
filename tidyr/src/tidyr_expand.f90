! SPDX-License-Identifier: MIT
! SPDX-FileComment: Expansion operations for the Fortran tidyr translation.
module tidyr_expand
   !! Implements Cartesian grids, sorted unique crossings, sequences, and row counts.
   use tibble, only: drop_columns, make_column, new_tibble, tibble_column, tibble_type
   use vctrs, only: vec_order, vec_slice, vec_unique
   implicit none
   private

   public :: crossing, expand_grid, full_seq, uncount

   interface uncount
      module procedure uncount_weights, uncount_column
   end interface uncount

   interface full_seq
      module procedure full_seq_integer, full_seq_real
   end interface full_seq

contains

   function expand_grid(columns) result(output)
      !! Forms the Cartesian product of named input vectors.
      type(tibble_column), intent(in) :: columns(:) !! Vectors defining dimensions.
      type(tibble_type) :: output
      type(tibble_column), allocatable :: out_columns(:)
      integer, allocatable :: indices(:)
      integer :: block, i, j, total

      if (size(columns) == 0) then
         output = new_tibble([tibble_column ::], nrow=1)
         return
      end if
      total = 1
      do j = 1, size(columns)
         if (columns(j)%size() == 0) then
            total = 0
            exit
         end if
         if (total > huge(total) / columns(j)%size()) error stop 'expand_grid: result is too large'
         total = total * columns(j)%size()
      end do
      allocate (out_columns(size(columns)), indices(total))
      if (total == 0) then
         do j = 1, size(columns)
            out_columns(j) = vec_slice(columns(j), indices)
         end do
         output = new_tibble(out_columns, nrow=0)
         return
      end if
      block = total
      do j = 1, size(columns)
         block = block / columns(j)%size()
         do i = 1, total
            indices(i) = modulo((i - 1) / block, columns(j)%size()) + 1
         end do
         out_columns(j) = vec_slice(columns(j), indices)
      end do
      output = new_tibble(out_columns, nrow=total)
   end function expand_grid

   function crossing(columns) result(output)
      !! Forms a Cartesian product after sorting and deduplicating each vector.
      type(tibble_column), intent(in) :: columns(:) !! Vectors defining dimensions.
      type(tibble_type) :: output
      type(tibble_column), allocatable :: unique_columns(:)
      integer, allocatable :: order(:)
      integer :: j
      allocate (unique_columns(size(columns)))
      do j = 1, size(columns)
         unique_columns(j) = vec_unique(columns(j))
         order = vec_order(unique_columns(j))
         unique_columns(j) = vec_slice(unique_columns(j), order)
      end do
      output = expand_grid(unique_columns)
   end function crossing

   pure function full_seq_integer(values, period) result(sequence)
      !! Completes an integer sequence at a positive period.
      integer, intent(in) :: values(:) !! Existing sequence values.
      integer, intent(in) :: period    !! Positive step.
      integer, allocatable :: sequence(:)
      integer :: first, j, last, n
      if (size(values) == 0) then
         allocate (sequence(0))
         return
      end if
      if (period <= 0) error stop 'full_seq: period must be positive'
      first = minval(values)
      last = maxval(values)
      if (modulo(last - first, period) /= 0) error stop 'full_seq: range is not aligned to period'
      n = (last - first) / period + 1
      allocate (sequence(n))
      sequence = [(first + (j - 1) * period, j = 1, n)]
   end function full_seq_integer

   pure function full_seq_real(values, period, tolerance) result(sequence)
      !! Completes a real sequence while checking alignment within tolerance.
      real(kind(1.0d0)), intent(in) :: values(:) !! Existing sequence values.
      real(kind(1.0d0)), intent(in) :: period    !! Positive step.
      real(kind(1.0d0)), intent(in), optional :: tolerance !! Alignment tolerance.
      real(kind(1.0d0)), allocatable :: sequence(:)
      real(kind(1.0d0)) :: first, last, tol
      integer :: j, n
      if (size(values) == 0) then
         allocate (sequence(0))
         return
      end if
      if (period <= 0.0d0) error stop 'full_seq: period must be positive'
      tol = 1.0d-6
      if (present(tolerance)) tol = tolerance
      first = minval(values)
      last = maxval(values)
      n = nint((last - first) / period) + 1
      if (abs(first + real(n - 1, kind(1.0d0)) * period - last) > tol) then
         error stop 'full_seq: range is not aligned to period'
      end if
      allocate (sequence(n))
      sequence = [(first + real(j - 1, kind(1.0d0)) * period, j = 1, n)]
   end function full_seq_real

   function uncount_weights(data, weights, id) result(output)
      !! Replicates rows according to nonnegative integer weights.
      type(tibble_type), intent(in) :: data !! Input table.
      integer, intent(in) :: weights(:)     !! Repetition count per row `(nrow)`.
      character(len=*), intent(in), optional :: id !! Optional within-row sequence column.
      type(tibble_type) :: output
      type(tibble_column), allocatable :: columns(:)
      integer, allocatable :: ids(:), rows(:)
      integer :: i, j, k, total
      if (size(weights) /= data%nrow()) error stop 'uncount: weight count differs from row count'
      if (any(weights < 0)) error stop 'uncount: weights must be nonnegative'
      total = sum(weights)
      allocate (rows(total), ids(total))
      k = 0
      do i = 1, data%nrow()
         do j = 1, weights(i)
            k = k + 1
            rows(k) = i
            ids(k) = j
         end do
      end do
      allocate (columns(data%ncol()))
      do j = 1, data%ncol()
         columns(j) = vec_slice(data%columns(j), rows)
      end do
      output = new_tibble(columns, nrow=total)
      if (present(id)) then
         output = new_tibble([output%columns, make_column(id, ids)], nrow=total)
      end if
   end function uncount_weights

   function uncount_column(data, weights, id) result(output)
      !! Replicates rows using a named nonmissing integer weight column and removes it.
      type(tibble_type), intent(in) :: data !! Input table.
      character(len=*), intent(in) :: weights !! Integer weight-column name.
      character(len=*), intent(in), optional :: id !! Optional within-row sequence column.
      type(tibble_type) :: output, without_weights
      integer :: index
      index = data%column_index(weights)
      if (index == 0) error stop 'uncount: weight column not found'
      if (data%columns(index)%type_code /= 1) error stop 'uncount: weight column must be integer'
      if (any(data%columns(index)%missing)) error stop 'uncount: weights cannot be missing'
      without_weights = drop_columns(data, [character(len=len(weights)) :: weights])
      output = uncount_weights(without_weights, data%columns(index)%integer_values, id)
   end function uncount_column

end module tidyr_expand
