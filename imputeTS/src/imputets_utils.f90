module imputets_utils
   use, intrinsic :: ieee_arithmetic, only : ieee_is_nan, ieee_quiet_nan, ieee_value
   use imputets_kinds, only : dp
   implicit none
   private

   public :: apply_maxgap
   public :: count_observed
   public :: lower_string
   public :: missing_value
   public :: observed_values
   public :: require_observed
   public :: sorted_copy

contains

   pure elemental logical function is_missing(x) result(missing)
      real(dp), intent(in) :: x !! Value tested for an IEEE NaN used as the missing-value marker.
      missing = ieee_is_nan(x)
   end function is_missing

   pure function missing_value() result(value)
      real(dp) :: value
      value = ieee_value(0.0_dp, ieee_quiet_nan)
   end function missing_value

   pure integer function count_observed(x) result(nobs)
      real(dp), intent(in) :: x(:) !! Series whose finite-or-infinite non-NaN observations are counted.
      integer :: i

      nobs = 0
      do i = 1, size(x)
         if (.not. is_missing(x(i))) nobs = nobs + 1
      end do
   end function count_observed

   pure subroutine require_observed(x, minimum, routine_name)
      real(dp), intent(in) :: x(:) !! Series to validate for a minimum number of non-NaN observations.
      integer, intent(in) :: minimum !! Minimum required count of non-NaN observations; must be nonnegative.
      character(len=*), intent(in) :: routine_name !! Procedure name used in any validation diagnostic.

      if (count_observed(x) < minimum) then
         error stop trim(routine_name) // ': too few non-NaN observations'
      end if
   end subroutine require_observed

   pure function observed_values(x) result(values)
      real(dp), intent(in) :: x(:) !! Series from which NaN entries are removed while preserving order.
      real(dp), allocatable :: values(:)
      integer :: i
      integer :: j
      integer :: nobs

      nobs = count_observed(x)
      allocate(values(nobs))
      j = 0
      do i = 1, size(x)
         if (.not. is_missing(x(i))) then
            j = j + 1
            values(j) = x(i)
         end if
      end do
   end function observed_values

   pure function sorted_copy(x) result(y)
      real(dp), intent(in) :: x(:) !! Values to sort in ascending numerical order.
      real(dp), allocatable :: y(:)
      integer :: i
      integer :: j
      real(dp) :: key

      y = x
      do i = 2, size(y)
         key = y(i)
         j = i - 1
         do while (j >= 1)
            if (y(j) <= key) exit
            y(j + 1) = y(j)
            j = j - 1
         end do
         y(j + 1) = key
      end do
   end function sorted_copy

   pure function lower_string(text) result(lower)
      character(len=*), intent(in) :: text !! Text converted to lower-case ASCII without otherwise changing characters.
      character(len=len(text)) :: lower
      integer :: i
      integer :: code

      lower = text
      do i = 1, len(text)
         code = iachar(text(i:i))
         if (code >= iachar('A') .and. code <= iachar('Z')) lower(i:i) = achar(code + 32)
      end do
   end function lower_string

   pure subroutine apply_maxgap(original, imputed, maxgap)
      real(dp), intent(in) :: original(:) !! Original series whose NaN runs determine which imputations may be retained.
      real(dp), intent(inout) :: imputed(:) !! Imputed series; long original gaps are restored to NaN in place.
      integer, intent(in), optional :: maxgap !! Largest NaN-run length to impute; absent or negative means unlimited.
      integer :: i
      integer :: j
      integer :: limit
      integer :: run_length
      real(dp) :: nan

      if (.not. present(maxgap)) return
      limit = maxgap
      if (limit < 0) return
      if (size(imputed) /= size(original)) error stop 'apply_maxgap: shape mismatch'

      nan = missing_value()
      i = 1
      do while (i <= size(original))
         if (.not. is_missing(original(i))) then
            i = i + 1
            cycle
         end if
         j = i
         do while (j <= size(original))
            if (.not. is_missing(original(j))) exit
            j = j + 1
         end do
         run_length = j - i
         if (run_length > limit) imputed(i:j - 1) = nan
         i = j
      end do
   end subroutine apply_maxgap

end module imputets_utils
