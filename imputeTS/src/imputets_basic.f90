module imputets_basic
   use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
   use, intrinsic :: iso_fortran_env, only : int64
   use imputets_kinds, only : dp
   use imputets_utils, only : apply_maxgap, count_observed, lower_string, observed_values, sorted_copy
   implicit none
   private

   public :: na_locf
   public :: na_ma
   public :: na_mean
   public :: na_random
   public :: na_remove
   public :: na_replace

contains

   pure function na_remove(x) result(out)
      real(dp), intent(in) :: x(:) !! Numeric series from which all NaN missing values are removed.
      real(dp), allocatable :: out(:)

      if (count_observed(x) == 0) error stop 'na_remove: input contains only NaN values'
      out = observed_values(x)
   end function na_remove

   pure function na_replace(x, fill, maxgap) result(out)
      real(dp), intent(in) :: x(:) !! Numeric series whose NaN entries are candidates for replacement.
      real(dp), intent(in), optional :: fill !! Replacement value; defaults to zero and may itself be NaN.
      integer, intent(in), optional :: maxgap !! Largest NaN-run length to replace; absent or negative means unlimited.
      real(dp), allocatable :: out(:)
      real(dp) :: replacement
      integer :: i

      replacement = 0.0_dp
      if (present(fill)) replacement = fill
      out = x
      do i = 1, size(out)
         if (ieee_is_nan(out(i))) out(i) = replacement
      end do
      call apply_maxgap(x, out, maxgap)
   end function na_replace

   pure function na_mean(x, option, maxgap) result(out)
      real(dp), intent(in) :: x(:) !! Numeric series whose NaN entries are replaced by a global location statistic.
      character(len=*), intent(in), optional :: option !! Mean type: mean, median, mode, geometric, or harmonic.
      integer, intent(in), optional :: maxgap !! Largest NaN-run length to impute; absent or negative means unlimited.
      real(dp), allocatable :: out(:)
      real(dp), allocatable :: values(:)
      real(dp), allocatable :: sorted(:)
      character(len=:), allocatable :: method
      real(dp) :: replacement
      integer :: best_count
      integer :: count_here
      integer :: i
      integer :: j
      integer :: n

      if (count_observed(x) < 1) error stop 'na_mean: at least one non-NaN observation is required'
      method = 'mean'
      if (present(option)) method = trim(lower_string(option))
      values = observed_values(x)
      n = size(values)

      select case (method)
      case ('mean')
         replacement = sum(values) / real(n, dp)
      case ('median')
         allocate(sorted(n))
         sorted = sorted_copy(values)
         if (mod(n, 2) == 0) then
            replacement = 0.5_dp * (sorted(n / 2) + sorted(n / 2 + 1))
         else
            replacement = sorted((n + 1) / 2)
         end if
      case ('mode')
         allocate(sorted(n))
         sorted = sorted_copy(values)
         replacement = sorted(1)
         best_count = 0
         i = 1
         do while (i <= n)
            j = i + 1
            do while (j <= n)
               if (transfer(sorted(j), 0_int64) /= transfer(sorted(i), 0_int64)) exit
               j = j + 1
            end do
            count_here = j - i
            if (count_here > best_count) then
               best_count = count_here
               replacement = sorted(i)
            end if
            i = j
         end do
      case ('geometric')
         if (any(values <= 0.0_dp)) error stop 'na_mean: geometric mean requires positive observations'
         replacement = exp(sum(log(values)) / real(n, dp))
      case ('harmonic')
         if (any(values <= 0.0_dp)) error stop 'na_mean: harmonic mean requires positive observations'
         replacement = real(n, dp) / sum(1.0_dp / values)
      case default
         error stop 'na_mean: option must be mean, median, mode, geometric, or harmonic'
      end select

      out = x
      do i = 1, size(out)
         if (ieee_is_nan(out(i))) out(i) = replacement
      end do
      call apply_maxgap(x, out, maxgap)
   end function na_mean

   pure function na_locf(x, option, na_remaining, maxgap) result(out)
      real(dp), intent(in) :: x(:) !! Numeric series whose NaN entries are filled from neighboring observations.
      character(len=*), intent(in), optional :: option !! Direction: locf for forward carry or nocb for backward carry.
      character(len=*), intent(in), optional :: na_remaining !! Edge policy: rev, keep, rm, or mean; defaults to rev.
      integer, intent(in), optional :: maxgap !! Largest NaN-run length to impute; absent or negative means unlimited.
      real(dp), allocatable :: out(:)
      real(dp), allocatable :: compact(:)
      character(len=:), allocatable :: direction
      character(len=:), allocatable :: remaining
      integer :: i

      if (count_observed(x) < 1) error stop 'na_locf: at least one non-NaN observation is required'
      direction = 'locf'
      if (present(option)) direction = trim(lower_string(option))
      remaining = 'rev'
      if (present(na_remaining)) remaining = trim(lower_string(na_remaining))
      if (direction /= 'locf' .and. direction /= 'nocb') error stop 'na_locf: option must be locf or nocb'

      out = x
      if (direction == 'locf') then
         do i = 2, size(out)
            if (ieee_is_nan(out(i)) .and. .not. ieee_is_nan(out(i - 1))) out(i) = out(i - 1)
         end do
      else
         do i = size(out) - 1, 1, -1
            if (ieee_is_nan(out(i)) .and. .not. ieee_is_nan(out(i + 1))) out(i) = out(i + 1)
         end do
      end if

      if (any_nan(out)) then
         select case (remaining)
         case ('keep')
         case ('rev')
            if (direction == 'locf') then
               do i = size(out) - 1, 1, -1
                  if (ieee_is_nan(out(i)) .and. .not. ieee_is_nan(out(i + 1))) out(i) = out(i + 1)
               end do
            else
               do i = 2, size(out)
                  if (ieee_is_nan(out(i)) .and. .not. ieee_is_nan(out(i - 1))) out(i) = out(i - 1)
               end do
            end if
         case ('mean')
            out = na_mean(out)
         case ('rm')
            compact = observed_values(out)
            out = compact
            return
         case default
            error stop 'na_locf: na_remaining must be keep, rev, rm, or mean'
         end select
      end if
      call apply_maxgap(x, out, maxgap)
   end function na_locf

   pure function na_ma(x, k, weighting, maxgap) result(out)
      real(dp), intent(in) :: x(:) !! Numeric series whose NaN entries are replaced from a centered moving window.
      integer, intent(in), optional :: k !! Initial half-window size in observations; defaults to four and must be positive.
      character(len=*), intent(in), optional :: weighting !! Weighting rule: simple, linear, or exponential; default exponential.
      integer, intent(in), optional :: maxgap !! Largest NaN-run length to impute; absent or negative means unlimited.
      real(dp), allocatable :: out(:)
      character(len=:), allocatable :: method
      integer :: half_window
      integer :: i
      integer :: j
      integer :: left
      integer :: nobs
      integer :: right
      real(dp) :: distance
      real(dp) :: value_sum
      real(dp) :: weight
      real(dp) :: weight_sum

      if (count_observed(x) < 2) error stop 'na_ma: at least two non-NaN observations are required'
      half_window = 4
      if (present(k)) half_window = k
      if (half_window < 1) error stop 'na_ma: k must be positive'
      method = 'exponential'
      if (present(weighting)) method = trim(lower_string(weighting))
      if (method /= 'simple' .and. method /= 'linear' .and. method /= 'exponential') then
         error stop 'na_ma: weighting must be simple, linear, or exponential'
      end if

      out = x
      do i = 1, size(x)
         if (.not. ieee_is_nan(x(i))) cycle
         left = max(1, i - half_window)
         right = min(size(x), i + half_window)
         nobs = count_observed(x(left:right))
         do while (nobs < 2)
            left = max(1, left - 1)
            right = min(size(x), right + 1)
            nobs = count_observed(x(left:right))
         end do

         value_sum = 0.0_dp
         weight_sum = 0.0_dp
         do j = left, right
            if (ieee_is_nan(x(j))) cycle
            distance = real(abs(j - i), dp)
            select case (method)
            case ('simple')
               weight = 1.0_dp
            case ('linear')
               weight = 1.0_dp / (distance + 1.0_dp)
            case ('exponential')
               weight = 1.0_dp / 2.0_dp ** abs(j - i)
            end select
            value_sum = value_sum + weight * x(j)
            weight_sum = weight_sum + weight
         end do
         out(i) = value_sum / weight_sum
      end do
      call apply_maxgap(x, out, maxgap)
   end function na_ma

   function na_random(x, lower_bound, upper_bound, maxgap) result(out)
      real(dp), intent(in) :: x(:) !! Numeric series whose NaN entries are replaced by uniform random draws.
      real(dp), intent(in), optional :: lower_bound !! Lower uniform bound; defaults to the minimum observed value.
      real(dp), intent(in), optional :: upper_bound !! Upper uniform bound; defaults to the maximum observed value.
      integer, intent(in), optional :: maxgap !! Largest NaN-run length to impute; absent or negative means unlimited.
      real(dp), allocatable :: out(:)
      real(dp), allocatable :: values(:)
      real(dp) :: lower
      real(dp) :: upper
      real(dp) :: u
      integer :: i

      if (.not. present(lower_bound) .or. .not. present(upper_bound)) then
         if (count_observed(x) < 2) error stop 'na_random: at least two observations are required for default bounds'
      end if
      values = observed_values(x)
      if (present(lower_bound)) then
         lower = lower_bound
      else
         lower = minval(values)
      end if
      if (present(upper_bound)) then
         upper = upper_bound
      else
         upper = maxval(values)
      end if
      if (lower >= upper) error stop 'na_random: lower_bound must be smaller than upper_bound'

      out = x
      do i = 1, size(out)
         if (ieee_is_nan(out(i))) then
            call random_number(u)
            out(i) = lower + (upper - lower) * u
         end if
      end do
      call apply_maxgap(x, out, maxgap)
   end function na_random

   pure logical function any_nan(x) result(found)
      real(dp), intent(in) :: x(:) !! Series tested for at least one NaN missing value.
      integer :: i

      found = .false.
      do i = 1, size(x)
         if (ieee_is_nan(x(i))) then
            found = .true.
            return
         end if
      end do
   end function any_nan

end module imputets_basic
