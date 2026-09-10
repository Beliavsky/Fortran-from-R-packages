! SPDX-License-Identifier: MIT
module vctrs_core
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use vctrs_types, only: dp, new_vctr, vctr_type, vctrs_character, vctrs_integer, &
      vctrs_logical, vctrs_real
   implicit none
   private

   public :: vec_any_missing, vec_cast, vec_c, vec_compare, vec_copy_element, vec_detect_missing
   public :: vec_duplicate_any, vec_duplicate_detect, vec_equal, vec_group_id
   public :: vec_identify_runs, vec_if_else, vec_in, vec_match, vec_order
   public :: vec_init, vec_insert
   public :: vec_ptype2, vec_ptype_common, vec_recycle, vec_recycle_common
   public :: vec_rep, vec_rep_each, vec_size, vec_size_common, vec_slice
   public :: vec_sort, vec_unique, vec_unique_count, vec_unique_loc, vec_is

contains

   pure elemental integer function vec_size(vector) result(n)
      !! Returns a vector's size.
      type(vctr_type), intent(in) :: vector !! Vector to inspect.

      n = vector%size()
   end function vec_size

   pure elemental logical function vec_is(vector) result(valid)
      !! Checks a vector's type tag, active storage, and missing-value mask.
      type(vctr_type), intent(in) :: vector !! Vector to validate.
      integer :: allocated_count

      allocated_count = count([allocated(vector%integer_values), &
         allocated(vector%real_values), allocated(vector%logical_values), &
         allocated(vector%character_values)])
      valid = allocated_count == 1
      select case (vector%type_code)
      case (vctrs_integer)
         valid = valid .and. allocated(vector%integer_values)
      case (vctrs_real)
         valid = valid .and. allocated(vector%real_values)
      case (vctrs_logical)
         valid = valid .and. allocated(vector%logical_values)
      case (vctrs_character)
         valid = valid .and. allocated(vector%character_values)
      case default
         valid = .false.
      end select
      valid = valid .and. allocated(vector%missing)
      if (allocated(vector%missing)) valid = valid .and. size(vector%missing) == vector%size()
      valid = valid .and. allocated(vector%name)
   end function vec_is

   pure function vec_detect_missing(vector) result(missing)
      !! Returns the explicit elementwise missing-value mask.
      type(vctr_type), intent(in) :: vector !! Vector to inspect.
      logical, allocatable :: missing(:)

      if (.not. vec_is(vector)) error stop "vec_detect_missing: invalid vector"
      missing = vector%missing
   end function vec_detect_missing

   pure elemental logical function vec_any_missing(vector) result(any_missing)
      !! Reports whether a vector contains one or more missing values.
      type(vctr_type), intent(in) :: vector !! Vector to inspect.

      if (.not. vec_is(vector)) error stop "vec_any_missing: invalid vector"
      any_missing = any(vector%missing)
   end function vec_any_missing

   pure elemental integer function vec_ptype2(x, y) result(type_code)
      !! Determines the common type of two vectors, or zero when incompatible.
      type(vctr_type), intent(in) :: x !! First vector.
      type(vctr_type), intent(in) :: y !! Second vector.

      if (.not. vec_is(x) .or. .not. vec_is(y)) error stop "vec_ptype2: invalid vector"
      if (x%type_code == y%type_code) then
         type_code = x%type_code
      else if (is_numeric_type(x%type_code) .and. is_numeric_type(y%type_code)) then
         if (x%type_code == vctrs_real .or. y%type_code == vctrs_real) then
            type_code = vctrs_real
         else
            type_code = vctrs_integer
         end if
      else
         type_code = 0
      end if
   end function vec_ptype2

   pure elemental logical function is_numeric_type(type_code) result(is_numeric)
      !! Reports whether a type participates in logical/integer/real promotion.
      integer, intent(in) :: type_code !! Type tag to inspect.

      is_numeric = any(type_code == [vctrs_logical, vctrs_integer, vctrs_real])
   end function is_numeric_type

   pure integer function vec_ptype_common(vectors) result(type_code)
      !! Determines the common type of a vector collection.
      type(vctr_type), intent(in) :: vectors(:) !! Vectors to combine.
      integer :: j

      type_code = 0
      if (size(vectors) == 0) return
      if (.not. vec_is(vectors(1))) error stop "vec_ptype_common: invalid vector"
      type_code = vectors(1)%type_code
      do j = 2, size(vectors)
         if (.not. vec_is(vectors(j))) error stop "vec_ptype_common: invalid vector"
         type_code = common_type_codes(type_code, vectors(j)%type_code)
         if (type_code == 0) return
      end do
   end function vec_ptype_common

   pure elemental integer function common_type_codes(x_type, y_type) result(type_code)
      !! Determines a common type directly from two valid type tags.
      integer, intent(in) :: x_type !! First type tag.
      integer, intent(in) :: y_type !! Second type tag.

      if (x_type == y_type) then
         type_code = x_type
      else if (is_numeric_type(x_type) .and. is_numeric_type(y_type)) then
         if (x_type == vctrs_real .or. y_type == vctrs_real) then
            type_code = vctrs_real
         else
            type_code = vctrs_integer
         end if
      else
         type_code = 0
      end if
   end function common_type_codes

   pure elemental function vec_cast(vector, type_code) result(out)
      !! Casts to a requested type when all nonmissing values are representable.
      type(vctr_type), intent(in) :: vector !! Vector to cast.
      integer, intent(in) :: type_code      !! Requested output type tag.
      type(vctr_type) :: out
      integer, allocatable :: integer_values(:)
      integer :: i

      if (.not. vec_is(vector)) error stop "vec_cast: invalid vector"
      if (vector%type_code == type_code) then
         out = vector
         return
      end if
      select case (type_code)
      case (vctrs_real)
         select case (vector%type_code)
         case (vctrs_integer)
            out = new_vctr(vector%name, real(vector%integer_values, dp), vector%missing)
         case (vctrs_logical)
            out = new_vctr(vector%name, merge(1.0_dp, 0.0_dp, vector%logical_values), &
               vector%missing)
         case default
            error stop "vec_cast: incompatible cast to real"
         end select
      case (vctrs_integer)
         select case (vector%type_code)
         case (vctrs_logical)
            out = new_vctr(vector%name, merge(1, 0, vector%logical_values), vector%missing)
         case (vctrs_real)
            allocate (integer_values(vector%size()), source=0)
            do i = 1, vector%size()
               if (vector%missing(i)) cycle
               if (.not. ieee_is_finite(vector%real_values(i))) then
                  error stop "vec_cast: non-finite real cannot be cast to integer"
               end if
               if (abs(vector%real_values(i) - anint(vector%real_values(i))) > 0.0_dp) then
                  error stop "vec_cast: lossy real-to-integer cast"
               end if
               if (vector%real_values(i) < real(-huge(0), dp) .or. &
                   vector%real_values(i) > real(huge(0), dp)) then
                  error stop "vec_cast: real is outside the integer range"
               end if
               integer_values(i) = int(vector%real_values(i))
            end do
            out = new_vctr(vector%name, integer_values, vector%missing)
         case default
            error stop "vec_cast: incompatible cast to integer"
         end select
      case (vctrs_logical)
         select case (vector%type_code)
         case (vctrs_integer)
            do i = 1, vector%size()
               if (.not. vector%missing(i) .and. &
                   vector%integer_values(i) /= 0 .and. vector%integer_values(i) /= 1) then
                  error stop "vec_cast: lossy integer-to-logical cast"
               end if
            end do
            out = new_vctr(vector%name, vector%integer_values /= 0, vector%missing)
         case (vctrs_real)
            do i = 1, vector%size()
               if (.not. vector%missing(i) .and. abs(vector%real_values(i)) > 0.0_dp .and. &
                   abs(vector%real_values(i) - 1.0_dp) > 0.0_dp) then
                  error stop "vec_cast: lossy real-to-logical cast"
               end if
            end do
            out = new_vctr(vector%name, abs(vector%real_values) > 0.0_dp, vector%missing)
         case default
            error stop "vec_cast: incompatible cast to logical"
         end select
      case (vctrs_character)
         error stop "vec_cast: only character vectors cast to character"
      case default
         error stop "vec_cast: invalid requested type"
      end select
   end function vec_cast

   pure integer function vec_size_common(vectors) result(common_size)
      !! Computes the common size under size-one recycling rules.
      type(vctr_type), intent(in) :: vectors(:) !! Vectors to reconcile.
      integer :: j, current_size

      common_size = 1
      if (size(vectors) == 0) then
         common_size = 0
         return
      end if
      do j = 1, size(vectors)
         if (.not. vec_is(vectors(j))) error stop "vec_size_common: invalid vector"
         current_size = vectors(j)%size()
         if (current_size == 1) cycle
         if (common_size == 1) then
            common_size = current_size
         else if (current_size /= common_size) then
            error stop "vec_size_common: incompatible vector sizes"
         end if
      end do
   end function vec_size_common

   pure elemental function vec_recycle(vector, size_out) result(out)
      !! Recycles a size-one vector or preserves a vector of the requested size.
      type(vctr_type), intent(in) :: vector !! Vector to recycle.
      integer, intent(in) :: size_out       !! Required nonnegative size.
      type(vctr_type) :: out

      if (.not. vec_is(vector)) error stop "vec_recycle: invalid vector"
      if (size_out < 0) error stop "vec_recycle: output size must be nonnegative"
      if (vector%size() == size_out) then
         out = vector
      else if (vector%size() == 1) then
         select case (vector%type_code)
         case (vctrs_integer)
            out = new_vctr(vector%name, vector%integer_values(1), size_out, vector%missing(1))
         case (vctrs_real)
            out = new_vctr(vector%name, vector%real_values(1), size_out, vector%missing(1))
         case (vctrs_logical)
            out = new_vctr(vector%name, vector%logical_values(1), size_out, vector%missing(1))
         case (vctrs_character)
            out = new_vctr(vector%name, vector%character_values(1), size_out, vector%missing(1))
         end select
      else
         error stop "vec_recycle: only size-one vectors can be recycled"
      end if
   end function vec_recycle

   pure function vec_recycle_common(vectors) result(out)
      !! Recycles all inputs to their common size.
      type(vctr_type), intent(in) :: vectors(:) !! Vectors to recycle.
      type(vctr_type), allocatable :: out(:)
      integer :: common_size, j

      common_size = vec_size_common(vectors)
      allocate (out(size(vectors)))
      do j = 1, size(vectors)
         out(j) = vec_recycle(vectors(j), common_size)
      end do
   end function vec_recycle_common

   pure function vec_slice(vector, indices) result(out)
      !! Selects elements by one-based indices, preserving order and repetition.
      type(vctr_type), intent(in) :: vector !! Vector to slice.
      integer, intent(in) :: indices(:)     !! One-based indices.
      type(vctr_type) :: out

      if (.not. vec_is(vector)) error stop "vec_slice: invalid vector"
      if (any(indices < 1) .or. any(indices > vector%size())) then
         error stop "vec_slice: index out of bounds"
      end if
      out%name = vector%name
      out%type_code = vector%type_code
      out%missing = vector%missing(indices)
      select case (vector%type_code)
      case (vctrs_integer)
         out%integer_values = vector%integer_values(indices)
      case (vctrs_real)
         out%real_values = vector%real_values(indices)
      case (vctrs_logical)
         out%logical_values = vector%logical_values(indices)
      case (vctrs_character)
         out%character_values = vector%character_values(indices)
      end select
   end function vec_slice

   pure function vec_c(vectors) result(out)
      !! Concatenates vectors after casting them to a common type.
      type(vctr_type), intent(in) :: vectors(:) !! Vectors to concatenate.
      type(vctr_type) :: out
      type(vctr_type) :: converted
      character(len=:), allocatable :: character_values(:)
      logical, allocatable :: missing(:)
      integer, allocatable :: integer_values(:)
      real(dp), allocatable :: real_values(:)
      logical, allocatable :: logical_values(:)
      integer :: first, j, last, max_length, total, type_code

      if (size(vectors) == 0) error stop "vec_c: at least one vector is required"
      type_code = vec_ptype_common(vectors)
      if (type_code == 0) error stop "vec_c: vectors have no common type"
      total = sum(vec_size(vectors))
      allocate (missing(total))
      select case (type_code)
      case (vctrs_integer)
         allocate (integer_values(total))
      case (vctrs_real)
         allocate (real_values(total))
      case (vctrs_logical)
         allocate (logical_values(total))
      case (vctrs_character)
         max_length = 0
         do j = 1, size(vectors)
            max_length = max(max_length, len(vectors(j)%character_values))
         end do
         allocate (character(len=max_length) :: character_values(total))
      end select
      first = 1
      do j = 1, size(vectors)
         converted = vec_cast(vectors(j), type_code)
         last = first + converted%size() - 1
         missing(first:last) = converted%missing
         select case (type_code)
         case (vctrs_integer)
            integer_values(first:last) = converted%integer_values
         case (vctrs_real)
            real_values(first:last) = converted%real_values
         case (vctrs_logical)
            logical_values(first:last) = converted%logical_values
         case (vctrs_character)
            character_values(first:last) = converted%character_values
         end select
         first = last + 1
      end do
      select case (type_code)
      case (vctrs_integer)
         out = new_vctr(vectors(1)%name, integer_values, missing)
      case (vctrs_real)
         out = new_vctr(vectors(1)%name, real_values, missing)
      case (vctrs_logical)
         out = new_vctr(vectors(1)%name, logical_values, missing)
      case (vctrs_character)
         out = new_vctr(vectors(1)%name, character_values, missing)
      end select
   end function vec_c

   pure elemental function vec_insert(vector, values, before) result(out)
      !! Inserts compatible values before a one-based vector position.
      type(vctr_type), intent(in) :: vector !! Existing vector.
      type(vctr_type), intent(in) :: values !! Values to insert.
      integer, intent(in) :: before         !! Position in `1:size+1`.
      type(vctr_type) :: out
      type(vctr_type) :: inputs(2), converted_values
      character(len=:), allocatable :: character_values(:)
      logical, allocatable :: missing(:)
      integer, allocatable :: integer_values(:)
      real(dp), allocatable :: real_values(:)
      logical, allocatable :: logical_values(:)
      integer :: common_type, first_count, total

      if (before < 1 .or. before > vector%size() + 1) then
         error stop "vec_insert: position out of bounds"
      end if
      inputs = [vector, values]
      common_type = vec_ptype_common(inputs)
      if (common_type == 0) error stop "vec_insert: vectors have no common type"
      out = vec_cast(vector, common_type)
      converted_values = vec_cast(values, common_type)
      first_count = before - 1
      total = out%size() + converted_values%size()
      allocate (missing(total))
      missing(:first_count) = out%missing(:first_count)
      missing(first_count + 1:first_count + converted_values%size()) = &
         converted_values%missing
      missing(first_count + converted_values%size() + 1:) = out%missing(before:)
      select case (common_type)
      case (vctrs_integer)
         allocate (integer_values(total))
         integer_values(:first_count) = out%integer_values(:first_count)
         integer_values(first_count + 1:first_count + converted_values%size()) = &
            converted_values%integer_values
         integer_values(first_count + converted_values%size() + 1:) = &
            out%integer_values(before:)
         out = new_vctr(vector%name, integer_values, missing)
      case (vctrs_real)
         allocate (real_values(total))
         real_values(:first_count) = out%real_values(:first_count)
         real_values(first_count + 1:first_count + converted_values%size()) = &
            converted_values%real_values
         real_values(first_count + converted_values%size() + 1:) = out%real_values(before:)
         out = new_vctr(vector%name, real_values, missing)
      case (vctrs_logical)
         allocate (logical_values(total))
         logical_values(:first_count) = out%logical_values(:first_count)
         logical_values(first_count + 1:first_count + converted_values%size()) = &
            converted_values%logical_values
         logical_values(first_count + converted_values%size() + 1:) = &
            out%logical_values(before:)
         out = new_vctr(vector%name, logical_values, missing)
      case (vctrs_character)
         allocate (character(len=max(len(out%character_values), &
            len(converted_values%character_values))) :: character_values(total))
         character_values(:first_count) = out%character_values(:first_count)
         character_values(first_count + 1:first_count + converted_values%size()) = &
            converted_values%character_values
         character_values(first_count + converted_values%size() + 1:) = &
            out%character_values(before:)
         out = new_vctr(vector%name, character_values, missing)
      end select
   end function vec_insert

   pure elemental function vec_equal(x, y, na_equal) result(out)
      !! Compares vectors elementwise after common-type casting and recycling.
      type(vctr_type), intent(in) :: x         !! First vector.
      type(vctr_type), intent(in) :: y         !! Second vector.
      logical, intent(in), optional :: na_equal !! Whether two missing values compare equal.
      type(vctr_type) :: out
      type(vctr_type) :: inputs(2), recycled(2), left, right
      logical, allocatable :: equal(:), missing(:)
      logical :: compare_missing
      integer :: i, type_code

      inputs = [x, y]
      type_code = vec_ptype_common(inputs)
      if (type_code == 0) error stop "vec_equal: vectors have no common type"
      recycled = vec_recycle_common(inputs)
      left = vec_cast(recycled(1), type_code)
      right = vec_cast(recycled(2), type_code)
      allocate (equal(left%size()), source=.false.)
      allocate (missing(left%size()), source=.false.)
      compare_missing = .false.
      if (present(na_equal)) compare_missing = na_equal
      do i = 1, left%size()
         if (left%missing(i) .or. right%missing(i)) then
            if (compare_missing) equal(i) = left%missing(i) .and. right%missing(i)
            if (.not. compare_missing) missing(i) = .true.
         else
            equal(i) = elements_equal(left, i, right, i)
         end if
      end do
      out = new_vctr("", equal, missing)
   end function vec_equal

   pure elemental function vec_compare(x, y) result(out)
      !! Returns elementwise -1, 0, or 1 after common casting and recycling.
      type(vctr_type), intent(in) :: x !! First vector.
      type(vctr_type), intent(in) :: y !! Second vector.
      type(vctr_type) :: out
      type(vctr_type) :: inputs(2), recycled(2), left, right
      integer, allocatable :: comparison(:)
      logical, allocatable :: missing(:)
      integer :: i, type_code

      inputs = [x, y]
      type_code = vec_ptype_common(inputs)
      if (type_code == 0) error stop "vec_compare: vectors have no common type"
      recycled = vec_recycle_common(inputs)
      left = vec_cast(recycled(1), type_code)
      right = vec_cast(recycled(2), type_code)
      allocate (comparison(left%size()), source=0)
      missing = left%missing .or. right%missing
      do i = 1, left%size()
         if (missing(i)) cycle
         if (elements_less(left, i, right, i)) comparison(i) = -1
         if (elements_less(right, i, left, i)) comparison(i) = 1
      end do
      out = new_vctr("", comparison, missing)
   end function vec_compare

   pure elemental logical function elements_equal(x, ix, y, iy) result(equal)
      !! Compares two nonmissing elements of equal type.
      type(vctr_type), intent(in) :: x !! First vector.
      integer, intent(in) :: ix        !! First element index.
      type(vctr_type), intent(in) :: y !! Second vector.
      integer, intent(in) :: iy        !! Second element index.

      select case (x%type_code)
      case (vctrs_integer)
         equal = x%integer_values(ix) == y%integer_values(iy)
      case (vctrs_real)
         equal = real_values_equal(x%real_values(ix), y%real_values(iy))
      case (vctrs_logical)
         equal = x%logical_values(ix) .eqv. y%logical_values(iy)
      case (vctrs_character)
         equal = x%character_values(ix) == y%character_values(iy)
      case default
         equal = .false.
      end select
   end function elements_equal

   pure elemental logical function real_values_equal(x, y) result(equal)
      !! Compares real values exactly without compiler-ambiguous equality syntax.
      real(dp), intent(in) :: x !! First value.
      real(dp), intent(in) :: y !! Second value.

      equal = x <= y .and. x >= y
   end function real_values_equal

   pure elemental logical function elements_less(x, ix, y, iy) result(less)
      !! Orders two nonmissing elements of equal type.
      type(vctr_type), intent(in) :: x !! First vector.
      integer, intent(in) :: ix        !! First element index.
      type(vctr_type), intent(in) :: y !! Second vector.
      integer, intent(in) :: iy        !! Second element index.

      select case (x%type_code)
      case (vctrs_integer)
         less = x%integer_values(ix) < y%integer_values(iy)
      case (vctrs_real)
         less = x%real_values(ix) < y%real_values(iy)
      case (vctrs_logical)
         less = .not. x%logical_values(ix) .and. y%logical_values(iy)
      case (vctrs_character)
         less = x%character_values(ix) < y%character_values(iy)
      case default
         less = .false.
      end select
   end function elements_less

   pure elemental function vec_match(needles, haystack, na_equal) result(out)
      !! Finds the first matching haystack location for each needle.
      type(vctr_type), intent(in) :: needles  !! Values to locate.
      type(vctr_type), intent(in) :: haystack !! Values to search.
      logical, intent(in), optional :: na_equal !! Whether missing values match.
      type(vctr_type) :: out
      type(vctr_type) :: inputs(2), left, right
      integer, allocatable :: locations(:)
      logical, allocatable :: missing(:)
      logical :: match_missing
      integer :: i, j, type_code

      inputs = [needles, haystack]
      type_code = vec_ptype_common(inputs)
      if (type_code == 0) error stop "vec_match: vectors have no common type"
      left = vec_cast(needles, type_code)
      right = vec_cast(haystack, type_code)
      allocate (locations(left%size()), source=0)
      allocate (missing(left%size()), source=.true.)
      match_missing = .true.
      if (present(na_equal)) match_missing = na_equal
      do i = 1, left%size()
         do j = 1, right%size()
            if (left%missing(i) .or. right%missing(j)) then
               if (.not. match_missing .or. .not. left%missing(i) .or. &
                   .not. right%missing(j)) cycle
            else if (.not. elements_equal(left, i, right, j)) then
               cycle
            end if
            locations(i) = j
            missing(i) = .false.
            exit
         end do
      end do
      out = new_vctr("", locations, missing)
   end function vec_match

   pure elemental function vec_in(needles, haystack, na_equal) result(out)
      !! Reports whether each needle occurs in the haystack.
      type(vctr_type), intent(in) :: needles  !! Values to locate.
      type(vctr_type), intent(in) :: haystack !! Values to search.
      logical, intent(in), optional :: na_equal !! Whether missing values match.
      type(vctr_type) :: out
      type(vctr_type) :: locations

      locations = vec_match(needles, haystack, na_equal)
      out = new_vctr("", .not. locations%missing)
   end function vec_in

   pure function vec_unique_loc(vector) result(locations)
      !! Returns positions of first occurrences, treating missing values as equal.
      type(vctr_type), intent(in) :: vector !! Vector to inspect.
      integer, allocatable :: locations(:)
      logical, allocatable :: keep(:)
      integer :: i, j

      if (.not. vec_is(vector)) error stop "vec_unique_loc: invalid vector"
      allocate (keep(vector%size()), source=.true.)
      do i = 2, vector%size()
         do j = 1, i - 1
            if (same_or_both_missing(vector, i, vector, j)) then
               keep(i) = .false.
               exit
            end if
         end do
      end do
      allocate (locations(count(keep)))
      j = 0
      do i = 1, vector%size()
         if (keep(i)) then
            j = j + 1
            locations(j) = i
         end if
      end do
   end function vec_unique_loc

   pure elemental function vec_unique(vector) result(out)
      !! Returns unique values in first-occurrence order.
      type(vctr_type), intent(in) :: vector !! Vector to deduplicate.
      type(vctr_type) :: out

      out = vec_slice(vector, vec_unique_loc(vector))
   end function vec_unique

   pure elemental integer function vec_unique_count(vector) result(n)
      !! Counts distinct values, treating all missing values as one value.
      type(vctr_type), intent(in) :: vector !! Vector to inspect.

      n = size(vec_unique_loc(vector))
   end function vec_unique_count

   pure function vec_duplicate_detect(vector) result(duplicate)
      !! Marks every member of a duplicated value group.
      type(vctr_type), intent(in) :: vector !! Vector to inspect.
      logical, allocatable :: duplicate(:)
      integer :: i, j

      if (.not. vec_is(vector)) error stop "vec_duplicate_detect: invalid vector"
      allocate (duplicate(vector%size()), source=.false.)
      do i = 1, vector%size()
         do j = i + 1, vector%size()
            if (same_or_both_missing(vector, i, vector, j)) then
               duplicate(i) = .true.
               duplicate(j) = .true.
            end if
         end do
      end do
   end function vec_duplicate_detect

   pure elemental logical function vec_duplicate_any(vector) result(any_duplicate)
      !! Reports whether any value, including missing, occurs more than once.
      type(vctr_type), intent(in) :: vector !! Vector to inspect.

      any_duplicate = any(vec_duplicate_detect(vector))
   end function vec_duplicate_any

   pure elemental logical function same_or_both_missing(x, ix, y, iy) result(same)
      !! Compares elements with two missing values considered equal.
      type(vctr_type), intent(in) :: x !! First vector.
      integer, intent(in) :: ix        !! First element index.
      type(vctr_type), intent(in) :: y !! Second vector.
      integer, intent(in) :: iy        !! Second element index.

      if (x%missing(ix) .or. y%missing(iy)) then
         same = x%missing(ix) .and. y%missing(iy)
      else
         same = elements_equal(x, ix, y, iy)
      end if
   end function same_or_both_missing

   pure function vec_order(vector, descending) result(order)
      !! Returns a stable order with missing values last.
      type(vctr_type), intent(in) :: vector !! Vector to order.
      logical, intent(in), optional :: descending !! Reverse nonmissing order.
      integer, allocatable :: order(:)
      logical :: reverse
      integer :: i, j, key

      if (.not. vec_is(vector)) error stop "vec_order: invalid vector"
      reverse = .false.
      if (present(descending)) reverse = descending
      allocate (order(vector%size()))
      do i = 1, size(order)
         order(i) = i
      end do
      do i = 2, size(order)
         key = order(i)
         j = i - 1
         do while (j >= 1)
            if (.not. ordered_before(vector, key, order(j), reverse)) exit
            order(j + 1) = order(j)
            j = j - 1
         end do
         order(j + 1) = key
      end do
   end function vec_order

   pure elemental logical function ordered_before(vector, left, right, descending) result(before)
      !! Compares two locations for stable ordering with missing values last.
      type(vctr_type), intent(in) :: vector !! Vector being ordered.
      integer, intent(in) :: left           !! Candidate location.
      integer, intent(in) :: right          !! Existing location.
      logical, intent(in) :: descending     !! Whether nonmissing values descend.

      if (vector%missing(left)) then
         before = .false.
      else if (vector%missing(right)) then
         before = .true.
      else if (descending) then
         before = elements_less(vector, right, vector, left)
      else
         before = elements_less(vector, left, vector, right)
      end if
   end function ordered_before

   pure elemental function vec_sort(vector, descending) result(out)
      !! Returns a stable sorted copy with missing values last.
      type(vctr_type), intent(in) :: vector !! Vector to sort.
      logical, intent(in), optional :: descending !! Reverse nonmissing order.
      type(vctr_type) :: out

      out = vec_slice(vector, vec_order(vector, descending))
   end function vec_sort

   pure function vec_group_id(vector) result(groups)
      !! Assigns dense group identifiers in first-occurrence order.
      type(vctr_type), intent(in) :: vector !! Vector to group.
      integer, allocatable :: groups(:)
      integer :: i, j, next_group

      if (.not. vec_is(vector)) error stop "vec_group_id: invalid vector"
      allocate (groups(vector%size()), source=0)
      next_group = 0
      do i = 1, vector%size()
         do j = 1, i - 1
            if (same_or_both_missing(vector, i, vector, j)) then
               groups(i) = groups(j)
               exit
            end if
         end do
         if (groups(i) == 0) then
            next_group = next_group + 1
            groups(i) = next_group
         end if
      end do
   end function vec_group_id

   pure function vec_identify_runs(vector) result(runs)
      !! Assigns identifiers to consecutive runs of equal values.
      type(vctr_type), intent(in) :: vector !! Vector to inspect.
      integer, allocatable :: runs(:)
      integer :: i

      if (.not. vec_is(vector)) error stop "vec_identify_runs: invalid vector"
      allocate (runs(vector%size()))
      if (vector%size() == 0) return
      runs(1) = 1
      do i = 2, vector%size()
         runs(i) = runs(i - 1)
         if (.not. same_or_both_missing(vector, i - 1, vector, i)) runs(i) = runs(i) + 1
      end do
   end function vec_identify_runs

   pure elemental function vec_rep(vector, times) result(out)
      !! Repeats an entire vector a nonnegative number of times.
      type(vctr_type), intent(in) :: vector !! Vector to repeat.
      integer, intent(in) :: times          !! Repetition count.
      type(vctr_type) :: out
      integer, allocatable :: indices(:)
      integer :: i, j, k

      if (times < 0) error stop "vec_rep: times must be nonnegative"
      allocate (indices(vector%size() * times))
      k = 0
      do j = 1, times
         do i = 1, vector%size()
            k = k + 1
            indices(k) = i
         end do
      end do
      out = vec_slice(vector, indices)
   end function vec_rep

   pure elemental function vec_rep_each(vector, times) result(out)
      !! Repeats each vector element a nonnegative number of times.
      type(vctr_type), intent(in) :: vector !! Vector to repeat.
      integer, intent(in) :: times          !! Per-element repetition count.
      type(vctr_type) :: out
      integer, allocatable :: indices(:)
      integer :: i, j, k

      if (times < 0) error stop "vec_rep_each: times must be nonnegative"
      allocate (indices(vector%size() * times))
      k = 0
      do i = 1, vector%size()
         do j = 1, times
            k = k + 1
            indices(k) = i
         end do
      end do
      out = vec_slice(vector, indices)
   end function vec_rep_each

   pure elemental function vec_if_else(condition, true_value, false_value) result(out)
      !! Selects type-stable values under a logical condition with recycling.
      type(vctr_type), intent(in) :: condition   !! Logical condition vector.
      type(vctr_type), intent(in) :: true_value  !! Values selected for true elements.
      type(vctr_type), intent(in) :: false_value !! Values selected for false elements.
      type(vctr_type) :: out
      type(vctr_type) :: values(2), recycled_values(2), selected
      type(vctr_type) :: recycled_condition
      integer :: common_size, i, type_code

      if (condition%type_code /= vctrs_logical) error stop "vec_if_else: condition must be logical"
      values = [true_value, false_value]
      type_code = vec_ptype_common(values)
      if (type_code == 0) error stop "vec_if_else: values have no common type"
      common_size = vec_size_common([condition, true_value, false_value])
      recycled_condition = vec_recycle(condition, common_size)
      recycled_values(1) = vec_cast(vec_recycle(true_value, common_size), type_code)
      recycled_values(2) = vec_cast(vec_recycle(false_value, common_size), type_code)
      selected = recycled_values(1)
      selected%name = ""
      do i = 1, common_size
         if (recycled_condition%missing(i)) then
            selected%missing(i) = .true.
         else if (.not. recycled_condition%logical_values(i)) then
            call vec_copy_element(recycled_values(2), i, selected, i)
         end if
      end do
      out = selected
   end function vec_if_else

   pure elemental subroutine vec_copy_element(source, source_index, target, target_index)
      !! Copies one value and its missingness between equal-type vectors.
      type(vctr_type), intent(in) :: source !! Source vector.
      integer, intent(in) :: source_index   !! Source location.
      type(vctr_type), intent(inout) :: target !! Destination vector.
      integer, intent(in) :: target_index      !! Destination location.

      if (source%type_code /= target%type_code) then
         error stop "vec_copy_element: source and target types differ"
      end if
      target%missing(target_index) = source%missing(source_index)
      select case (source%type_code)
      case (vctrs_integer)
         target%integer_values(target_index) = source%integer_values(source_index)
      case (vctrs_real)
         target%real_values(target_index) = source%real_values(source_index)
      case (vctrs_logical)
         target%logical_values(target_index) = source%logical_values(source_index)
      case (vctrs_character)
         target%character_values(target_index) = source%character_values(source_index)
      end select
   end subroutine vec_copy_element

   pure elemental function vec_init(prototype, n, missing) result(out)
      !! Allocates a same-type vector of requested size with uniform missingness.
      type(vctr_type), intent(in) :: prototype !! Vector supplying type and name.
      integer, intent(in) :: n                 !! Nonnegative output size.
      logical, intent(in), optional :: missing !! Initial mask value; default false.
      type(vctr_type) :: out
      logical :: initial_missing

      if (n < 0) error stop "vec_init: size must be nonnegative"
      initial_missing = .false.
      if (present(missing)) initial_missing = missing
      out%name = prototype%name
      out%type_code = prototype%type_code
      allocate (out%missing(n), source=initial_missing)
      select case (prototype%type_code)
      case (vctrs_integer)
         allocate (out%integer_values(n), source=0)
      case (vctrs_real)
         allocate (out%real_values(n), source=0.0_dp)
      case (vctrs_logical)
         allocate (out%logical_values(n), source=.false.)
      case (vctrs_character)
         if (.not. allocated(prototype%character_values)) then
            error stop "vec_init: character prototype has no storage"
         end if
         allocate (character(len=len(prototype%character_values)) :: out%character_values(n))
         out%character_values = ""
      case default
         error stop "vec_init: unsupported vector type"
      end select
   end function vec_init

end module vctrs_core
