module car_data
   use car_kinds, only: dp
   use car_utils, only: median_valid
   use, intrinsic :: ieee_arithmetic, only: ieee_is_nan
   implicit none
   private

   public :: compare_coefs_numeric
   public :: match_coefficient_names
   public :: recode_numeric
   public :: some_indices
   public :: strings2factors_codes
   public :: tapply_numeric
   public :: which_names

contains

   pure subroutine which_names(object_names, query_names, indices, matched, status)
      character(len=*), intent(in) :: object_names(:) !! Names available on the source object, searched by exact trimmed equality.
      character(len=*), intent(in) :: query_names(:) !! Requested names returned in query order when an exact match exists.
      integer, allocatable, intent(out) :: indices(:) !! One-based first-match positions for matched query names; unmatched queries are omitted.
      logical, allocatable, intent(out) :: matched(:) !! Per-query flag indicating whether an exact source-name match was found.
      integer, intent(out) :: status !! Zero on success; one when either name vector is empty.
      integer :: i
      integer :: j
      integer :: k
      integer :: n_match

      if (size(object_names) == 0 .or. size(query_names) == 0) then
         allocate(indices(0), matched(size(query_names)))
         matched = .false.
         status = 1
         return
      end if

      allocate(matched(size(query_names)))
      matched = .false.
      n_match = 0
      do j = 1, size(query_names)
         do i = 1, size(object_names)
            if (trim(object_names(i)) == trim(query_names(j))) then
               matched(j) = .true.
               n_match = n_match + 1
               exit
            end if
         end do
      end do

      allocate(indices(n_match))
      k = 0
      do j = 1, size(query_names)
         if (.not. matched(j)) cycle
         do i = 1, size(object_names)
            if (trim(object_names(i)) == trim(query_names(j))) then
               k = k + 1
               indices(k) = i
               exit
            end if
         end do
      end do
      status = 0
   end subroutine which_names

   subroutine some_indices(n_items, n_select, indices, status, random_values)
      integer, intent(in) :: n_items !! Number of available elements or rows; must be nonnegative.
      integer, intent(in) :: n_select !! Requested sample size; values larger than n_items select all available items.
      integer, allocatable, intent(out) :: indices(:) !! Sorted one-based sample without replacement, matching the ordering returned by R some().
      integer, intent(out) :: status !! Zero on success; nonzero for invalid sizes or supplied random values outside [0,1).
      real(dp), intent(in), optional :: random_values(:) !! Optional deterministic uniforms in [0,1), one per selected item; omitted uses random_number.
      integer, allocatable :: pool(:)
      integer :: i
      integer :: j
      integer :: k
      integer :: key
      integer :: m
      real(dp) :: u

      if (n_items < 0 .or. n_select < 0) then
         allocate(indices(0))
         status = 1
         return
      end if
      m = min(n_items, n_select)
      if (present(random_values)) then
         if (size(random_values) < m) then
            allocate(indices(0))
            status = 2
            return
         end if
         if (m > 0) then
            if (any(random_values(1:m) < 0.0_dp) .or. any(random_values(1:m) >= 1.0_dp)) then
               allocate(indices(0))
               status = 3
               return
            end if
         end if
      end if

      allocate(indices(m), pool(n_items))
      do i = 1, n_items
         pool(i) = i
      end do
      do i = 1, m
         if (present(random_values)) then
            u = random_values(i)
         else
            call random_number(u)
         end if
         j = i + int(u * real(n_items - i + 1, dp))
         j = min(j, n_items)
         k = pool(i)
         pool(i) = pool(j)
         pool(j) = k
         indices(i) = pool(i)
      end do

      do i = 2, m
         key = indices(i)
         j = i - 1
         do while (j >= 1)
            if (indices(j) <= key) exit
            indices(j + 1) = indices(j)
            j = j - 1
         end do
         indices(j + 1) = key
      end do
      status = 0
   end subroutine some_indices

   pure subroutine strings2factors_codes(object, codes, levels, n_levels, converted, status, which, exclude_unique)
      character(len=*), intent(in) :: object(:, :) !! Character data with observations in rows and variables in columns.
      integer, allocatable, intent(out) :: codes(:, :) !! One-based factor codes for converted columns; zero marks columns not converted.
      character(len=:), allocatable, intent(out) :: levels(:, :) !! Sorted factor levels by row and column; unused entries are blank.
      integer, allocatable, intent(out) :: n_levels(:) !! Number of factor levels retained for each input column.
      logical, allocatable, intent(out) :: converted(:) !! True for columns converted to factor codes under the requested rules.
      integer, intent(out) :: status !! Zero on success; nonzero for an invalid selected-column index.
      integer, intent(in), optional :: which(:) !! Optional one-based columns to convert; omitted selects every character column.
      logical, intent(in), optional :: exclude_unique !! If true, all-unique columns are left unconverted; default is true as in R.
      character(len=len(object)) :: key
      character(len=len(object)), allocatable :: work(:)
      logical :: exclude_all_unique
      logical, allocatable :: selected(:)
      integer :: i
      integer :: j
      integer :: k
      integer :: n_unique
      integer :: nrow
      integer :: ncol

      nrow = size(object, 1)
      ncol = size(object, 2)
      allocate(codes(nrow, ncol), n_levels(ncol), converted(ncol), selected(ncol))
      allocate(character(len=len(object)) :: levels(max(1, nrow), ncol))
      codes = 0
      n_levels = 0
      converted = .false.
      selected = .true.
      levels = ''
      exclude_all_unique = .true.
      if (present(exclude_unique)) exclude_all_unique = exclude_unique

      if (present(which)) then
         selected = .false.
         do i = 1, size(which)
            if (which(i) < 1 .or. which(i) > ncol) then
               status = 1
               return
            end if
            selected(which(i)) = .true.
         end do
      end if

      allocate(work(nrow))
      do j = 1, ncol
         if (.not. selected(j)) cycle
         n_unique = 0
         do i = 1, nrow
            if (n_unique == 0) then
               n_unique = 1
               work(1) = object(i, j)
            else if (.not. any([(trim(work(k)) == trim(object(i, j)), k = 1, n_unique)])) then
               n_unique = n_unique + 1
               work(n_unique) = object(i, j)
            end if
         end do
         if (exclude_all_unique .and. n_unique == nrow) cycle

         do i = 2, n_unique
            key = work(i)
            k = i - 1
            do while (k >= 1)
               if (trim(work(k)) <= trim(key)) exit
               work(k + 1) = work(k)
               k = k - 1
            end do
            work(k + 1) = key
         end do
         converted(j) = .true.
         n_levels(j) = n_unique
         levels(1:n_unique, j) = work(1:n_unique)
         do i = 1, nrow
            do k = 1, n_unique
               if (trim(object(i, j)) == trim(work(k))) then
                  codes(i, j) = k
                  exit
               end if
            end do
         end do
      end do
      status = 0
   end subroutine strings2factors_codes

   pure subroutine compare_coefs_numeric(coefficients, covariance, std_errors, z_values, p_values, status)
      real(dp), intent(in) :: coefficients(:, :) !! Coefficient estimates with coefficients in rows and models in columns; NaNs mark absent estimates.
      real(dp), intent(in) :: covariance(:, :, :) !! Per-model coefficient covariance matrices with shape p by p by n_models.
      real(dp), allocatable, intent(out) :: std_errors(:, :) !! Square-root diagonal standard errors with the same shape as coefficients.
      real(dp), allocatable, intent(out) :: z_values(:, :) !! Estimate divided by standard error for each available coefficient.
      real(dp), allocatable, intent(out) :: p_values(:, :) !! Two-sided standard-normal p-values corresponding to z_values.
      integer, intent(out) :: status !! Zero on success; nonzero for incompatible dimensions or a negative covariance diagonal.
      real(dp) :: nan_value
      integer :: i
      integer :: j
      integer :: p
      integer :: n_models

      p = size(coefficients, 1)
      n_models = size(coefficients, 2)
      if (size(covariance, 1) /= p .or. size(covariance, 2) /= p .or. size(covariance, 3) /= n_models) then
         allocate(std_errors(0, 0), z_values(0, 0), p_values(0, 0))
         status = 1
         return
      end if
      do j = 1, n_models
         do i = 1, p
            if (.not. ieee_is_nan(coefficients(i, j)) .and. covariance(i, i, j) < 0.0_dp) then
               allocate(std_errors(0, 0), z_values(0, 0), p_values(0, 0))
               status = 2
               return
            end if
         end do
      end do

      allocate(std_errors(p, n_models), z_values(p, n_models), p_values(p, n_models))
      nan_value = ieee_nan()
      do j = 1, n_models
         do i = 1, p
            if (ieee_is_nan(coefficients(i, j))) then
               std_errors(i, j) = nan_value
               z_values(i, j) = nan_value
               p_values(i, j) = nan_value
            else
               std_errors(i, j) = sqrt(max(0.0_dp, covariance(i, i, j)))
               if (std_errors(i, j) <= tiny(1.0_dp)) then
                  z_values(i, j) = sign(huge(1.0_dp), coefficients(i, j))
                  p_values(i, j) = 0.0_dp
               else
                  z_values(i, j) = coefficients(i, j) / std_errors(i, j)
                  p_values(i, j) = erfc(abs(z_values(i, j)) / sqrt(2.0_dp))
               end if
            end if
         end do
      end do
      status = 0
   end subroutine compare_coefs_numeric

   pure subroutine match_coefficient_names(names, pattern, indices, status, case_sensitive)
      character(len=*), intent(in) :: names(:) !! Coefficient names searched for the requested literal pattern.
      character(len=*), intent(in) :: pattern !! Literal substring to match; regular-expression syntax is not interpreted.
      integer, allocatable, intent(out) :: indices(:) !! One-based indices of names containing pattern, in original order.
      integer, intent(out) :: status !! Zero on success; one when the pattern is empty.
      logical, intent(in), optional :: case_sensitive !! True for case-sensitive matching; default is true.
      character(len=:), allocatable :: needle
      character(len=:), allocatable :: candidate
      integer :: count_matches
      integer :: i
      integer :: j
      logical :: sensitive

      sensitive = .true.
      if (present(case_sensitive)) sensitive = case_sensitive
      needle = trim(pattern)
      if (.not. sensitive) needle = lowercase_text(needle)
      if (len_trim(needle) == 0) then
         allocate(indices(0))
         status = 1
         return
      end if

      count_matches = 0
      do i = 1, size(names)
         candidate = trim(names(i))
         if (.not. sensitive) candidate = lowercase_text(candidate)
         if (index(candidate, needle) > 0) count_matches = count_matches + 1
      end do
      allocate(indices(count_matches))
      j = 0
      do i = 1, size(names)
         candidate = trim(names(i))
         if (.not. sensitive) candidate = lowercase_text(candidate)
         if (index(candidate, needle) > 0) then
            j = j + 1
            indices(j) = i
         end if
      end do
      status = 0
   end subroutine match_coefficient_names

   pure function lowercase_text(text) result(lower)
      character(len=*), intent(in) :: text !! Text converted to ASCII lowercase for optional case-insensitive matching.
      character(len=len(text)) :: lower
      integer :: code
      integer :: i

      lower = text
      do i = 1, len(text)
         code = iachar(text(i:i))
         if (code >= iachar('A') .and. code <= iachar('Z')) lower(i:i) = achar(code + 32)
      end do
   end function lowercase_text

   subroutine recode_numeric(values, lower, upper, target, result, status, otherwise)
      real(dp), intent(in) :: values(:) !! Numeric values to recode; NaNs are preserved unless an otherwise value is supplied.
      real(dp), intent(in) :: lower(:) !! Inclusive lower bounds for ordered recode intervals.
      real(dp), intent(in) :: upper(:) !! Inclusive upper bounds for ordered recode intervals.
      real(dp), intent(in) :: target(:) !! Replacement value corresponding to each lower/upper interval.
      real(dp), allocatable, intent(out) :: result(:) !! Recoded numeric values, with later rules overriding earlier rules.
      integer, intent(out) :: status !! Zero on success; nonzero when rule arrays have incompatible lengths or bounds.
      real(dp), intent(in), optional :: otherwise !! Optional value assigned to all observations before interval rules are applied.
      integer :: i, j

      status = 0
      if (size(lower) /= size(upper) .or. size(lower) /= size(target) .or. any(lower > upper)) then
         status = 1
         allocate(result(0))
         return
      end if
      allocate(result(size(values)))
      if (present(otherwise)) then
         result = otherwise
      else
         result = values
      end if
      do j = 1, size(lower)
         do i = 1, size(values)
            if (.not. ieee_is_nan(values(i))) then
               if (values(i) >= lower(j) .and. values(i) <= upper(j)) result(i) = target(j)
            end if
         end do
      end do
   end subroutine recode_numeric

   subroutine tapply_numeric(values, group, operation, result, counts, status, na_rm)
      real(dp), intent(in) :: values(:) !! Numeric values to split by positive integer group labels.
      integer, intent(in) :: group(:) !! Positive integer group labels, one per value.
      character(len=*), intent(in) :: operation !! Reduction name: sum, mean, or median.
      real(dp), allocatable, intent(out) :: result(:) !! Reduction for each group number from one through max(group).
      integer, allocatable, intent(out) :: counts(:) !! Number of observations used in each group.
      integer, intent(out) :: status !! Zero on success; nonzero for invalid group labels or operation.
      logical, intent(in), optional :: na_rm !! If true, omit NaNs; default false, matching tapply-style propagation.
      logical :: remove_na
      logical, allocatable :: mask(:)
      integer :: g, ng, n_valid, median_status
      real(dp) :: median_value
      character(len=8) :: op

      status = 0
      if (size(group) /= size(values) .or. size(values) == 0 .or. any(group <= 0)) then
         status = 1
         allocate(result(0), counts(0))
         return
      end if
      op = lowercase(adjustl(operation))
      if (trim(op) /= "sum" .and. trim(op) /= "mean" .and. trim(op) /= "median") then
         status = 2
         allocate(result(0), counts(0))
         return
      end if
      remove_na = .false.
      if (present(na_rm)) remove_na = na_rm
      ng = maxval(group)
      allocate(result(ng), counts(ng), mask(size(values)))
      result = 0.0_dp
      counts = 0
      do g = 1, ng
         mask = group == g
         if (remove_na) mask = mask .and. (.not. ieee_is_nan(values))
         counts(g) = count(mask)
         if (counts(g) == 0) cycle
         if (.not. remove_na .and. any(ieee_is_nan(pack(values, mask)))) then
            result(g) = ieee_nan()
            cycle
         end if
         select case (trim(op))
         case ("sum")
            result(g) = sum(pack(values, mask))
         case ("mean")
            result(g) = sum(pack(values, mask)) / real(counts(g), dp)
         case ("median")
            call median_valid(pack(values, mask), median_value, n_valid, median_status)
            if (median_status /= 0) then
               result(g) = ieee_nan()
            else
               result(g) = median_value
            end if
         end select
      end do
   end subroutine tapply_numeric

   pure function ieee_nan() result(value)
      use, intrinsic :: ieee_arithmetic, only: ieee_quiet_nan, ieee_value
      real(dp) :: value
      value = ieee_value(0.0_dp, ieee_quiet_nan)
   end function ieee_nan

   pure function lowercase(text) result(lower)
      character(len=*), intent(in) :: text !! Option text converted to lowercase.
      character(len=len(text)) :: lower
      integer :: i, code

      lower = text
      do i = 1, len(text)
         code = iachar(text(i:i))
         if (code >= iachar('A') .and. code <= iachar('Z')) lower(i:i) = achar(code + 32)
      end do
   end function lowercase

end module car_data
