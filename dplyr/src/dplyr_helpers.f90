! SPDX-License-Identifier: MIT
module dplyr_helpers
   use dplyr_vector_ops
   use vctrs, only: dp, vctr_type, vec_c, vec_cast, vec_copy_element, vec_identify_runs, &
      vec_init, vec_ptype_common, vec_recycle_common, vec_slice, vec_unique_count
   implicit none
   private

   public :: between, coalesce, consecutive_id, first, lag, last, lead, n_distinct, na_if, near, nth

contains

   pure elemental function between(vector, left, right) result(out)
      !! Tests inclusively whether numeric values lie between scalar bounds.
      type(vctr_type), intent(in) :: vector !! Numeric values.
      real(dp), intent(in) :: left          !! Inclusive lower bound.
      real(dp), intent(in) :: right         !! Inclusive upper bound.
      type(vctr_type) :: out

      out = (vector >= left) .and. (vector <= right)
   end function between

   pure elemental function near(left, right, tolerance) result(out)
      !! Compares numeric vectors using an absolute tolerance.
      type(vctr_type), intent(in) :: left  !! Left numeric operand.
      type(vctr_type), intent(in) :: right !! Right numeric operand.
      real(dp), intent(in), optional :: tolerance !! Absolute tolerance.
      type(vctr_type) :: out
      real(dp) :: tol

      tol = sqrt(epsilon(1.0_dp))
      if (present(tolerance)) tol = tolerance
      if (tol < 0.0_dp) error stop "near: tolerance must be nonnegative"
      out = vector_abs(left - right) <= tol
   end function near

   pure function coalesce(vectors) result(out)
      !! Selects the first nonmissing value at each recycled position.
      type(vctr_type), intent(in) :: vectors(:) !! Compatible input vectors.
      type(vctr_type) :: out
      type(vctr_type), allocatable :: values(:)
      integer :: i, j, type_code

      if (size(vectors) == 0) error stop "coalesce: at least one vector is required"
      type_code = vec_ptype_common(vectors)
      if (type_code == 0) error stop "coalesce: vectors have incompatible types"
      values = vec_recycle_common(vectors)
      do j = 1, size(values)
         values(j) = vec_cast(values(j), type_code)
      end do
      out = vec_init(values(1), values(1)%size(), missing=.true.)
      do i = 1, out%size()
         do j = 1, size(values)
            if (.not. values(j)%missing(i)) then
               call vec_copy_element(values(j), i, out, i)
               exit
            end if
         end do
      end do
      out%name = vectors(1)%name
   end function coalesce

   pure elemental function na_if(vector, value) result(out)
      !! Replaces values equal to a scalar or recycled vector with missing values.
      type(vctr_type), intent(in) :: vector !! Vector to modify.
      type(vctr_type), intent(in) :: value  !! Comparison value, usually size one.
      type(vctr_type) :: out
      type(vctr_type) :: equal

      equal = vector == value
      if (equal%size() /= vector%size()) error stop "na_if: value has incompatible size"
      out = vector
      where (.not. equal%missing .and. equal%logical_values) out%missing = .true.
   end function na_if

   pure elemental function lag(vector, n) result(out)
      !! Shifts values toward higher indices and pads the front with missing values.
      type(vctr_type), intent(in) :: vector !! Vector to shift.
      integer, intent(in), optional :: n    !! Nonnegative lag; default one.
      type(vctr_type) :: out
      integer :: amount, i

      amount = 1
      if (present(n)) amount = n
      if (amount < 0) error stop "lag: n must be nonnegative"
      out = vec_init(vector, vector%size(), missing=.true.)
      do i = amount + 1, vector%size()
         call vec_copy_element(vector, i - amount, out, i)
      end do
   end function lag

   pure elemental function lead(vector, n) result(out)
      !! Shifts values toward lower indices and pads the end with missing values.
      type(vctr_type), intent(in) :: vector !! Vector to shift.
      integer, intent(in), optional :: n    !! Nonnegative lead; default one.
      type(vctr_type) :: out
      integer :: amount, i

      amount = 1
      if (present(n)) amount = n
      if (amount < 0) error stop "lead: n must be nonnegative"
      out = vec_init(vector, vector%size(), missing=.true.)
      do i = 1, max(0, vector%size() - amount)
         call vec_copy_element(vector, i + amount, out, i)
      end do
   end function lead

   pure elemental function nth(vector, position) result(out)
      !! Returns one position as a size-one vector; negative positions count from the end.
      type(vctr_type), intent(in) :: vector !! Vector to index.
      integer, intent(in) :: position       !! Nonzero position.
      type(vctr_type) :: out
      integer :: index

      if (position == 0) error stop "nth: position cannot be zero"
      index = position
      if (position < 0) index = vector%size() + position + 1
      if (index < 1 .or. index > vector%size()) then
         out = vec_init(vector, 1, missing=.true.)
      else
         out = vec_slice(vector, [index])
      end if
   end function nth

   pure elemental function first(vector) result(out)
      !! Returns the first element as a size-one vector, or missing for an empty input.
      type(vctr_type), intent(in) :: vector !! Vector to inspect.
      type(vctr_type) :: out
      out = nth(vector, 1)
   end function first

   pure elemental function last(vector) result(out)
      !! Returns the last element as a size-one vector, or missing for an empty input.
      type(vctr_type), intent(in) :: vector !! Vector to inspect.
      type(vctr_type) :: out
      out = nth(vector, -1)
   end function last

   pure elemental integer function n_distinct(vector) result(number)
      !! Counts distinct values, including one group for missing values.
      type(vctr_type), intent(in) :: vector !! Vector to inspect.
      number = vec_unique_count(vector)
   end function n_distinct

   pure function consecutive_id(vectors) result(ids)
      !! Assigns a new identifier whenever any selected value changes.
      type(vctr_type), intent(in) :: vectors(:) !! Equal-size vectors defining runs.
      integer, allocatable :: ids(:)
      integer, allocatable :: runs(:,:)
      integer :: i, j

      if (size(vectors) == 0) then
         allocate (ids(0))
         return
      end if
      allocate (ids(vectors(1)%size()), source=1)
      allocate (runs(size(ids), size(vectors)))
      do j = 1, size(vectors)
         if (vectors(j)%size() /= size(ids)) error stop "consecutive_id: vector sizes differ"
         runs(:, j) = vec_identify_runs(vectors(j))
      end do
      do i = 2, size(ids)
         do j = 1, size(vectors)
            if (runs(i, j) /= runs(i - 1, j)) then
               ids(i:) = ids(i:) + 1
               exit
            end if
         end do
      end do
   end function consecutive_id

end module dplyr_helpers
