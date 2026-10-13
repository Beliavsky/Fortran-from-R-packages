! SPDX-License-Identifier: Artistic-2.0
! Modern Fortran translation of matrixStats 1.5.0; see NOTICE.md.
module matrixstats_common
   use r_kinds, only: dp
   use r_ordering, only: r_sort_values_in_place
   use, intrinsic :: ieee_arithmetic, only: ieee_is_nan, ieee_is_finite, ieee_value, &
      ieee_quiet_nan, ieee_positive_inf, ieee_negative_inf
   implicit none
   private
   integer, parameter, public :: ms_na_index = -huge(0)
   public :: dp, ms_nan, ms_inf, ms_remove, ms_indices, ms_select, ms_vector
   public :: ms_sum, ms_mean, ms_var, ms_product, ms_extreme, ms_log_sum_exp
   public :: ms_count, ms_all, ms_any, ms_difference, ms_cumulative
contains
   pure real(dp) function ms_nan() result(y)
      y = ieee_value(0.0_dp, ieee_quiet_nan)
   end function

   pure real(dp) function ms_inf(negative) result(y)
      logical, intent(in), optional :: negative !! Whether to return negative infinity; default false.
      y = ieee_value(0.0_dp, ieee_positive_inf)
      if (present(negative)) then
         if (negative) y = ieee_value(0.0_dp, ieee_negative_inf)
      end if
   end function

   pure logical function ms_remove(na_rm) result(y)
      logical, intent(in), optional :: na_rm !! Remove NaN observations when true; absent means false.
      y = .false.
      if (present(na_rm)) y = na_rm
   end function

   pure function ms_indices(n, idxs) result(y)
      integer, intent(in) :: n !! Length of the dimension being selected; nonnegative.
      integer, intent(in), optional :: idxs(:) !! R indices: positive select, negative exclude, zero drops, ms_na_index missing.
      integer, allocatable :: y(:)
      logical, allocatable :: keep(:)
      integer :: i, k
      if (.not. present(idxs)) then
         y = [(i, i = 1, n)]
      else if (any(idxs < 0 .and. idxs /= ms_na_index)) then
         if (any(idxs > 0) .or. any(idxs == ms_na_index)) then
            y = [ms_na_index]
            return
         end if
         allocate(keep(n))
         keep = .true.
         do i = 1, size(idxs)
            k = -idxs(i)
            if (k >= 1 .and. k <= n) keep(k) = .false.
         end do
         y = pack([(i, i = 1, n)], keep)
      else
         y = pack(idxs, idxs /= 0)
         where (y < 1 .or. y > n) y = ms_na_index
      end if
   end function

   pure function ms_vector(x, idxs, na_rm) result(y)
      real(dp), intent(in) :: x(:) !! Source numeric vector; NaN represents R NA or NaN.
      integer, intent(in), optional :: idxs(:) !! Optional R-style vector selection; see ms_indices.
      logical, intent(in), optional :: na_rm !! Drop missing selected observations when true.
      real(dp), allocatable :: y(:)
      integer, allocatable :: ix(:)
      integer :: i
      ix = ms_indices(size(x), idxs)
      allocate(y(size(ix)))
      y = ms_nan()
      do i = 1, size(ix)
         if (ix(i) /= ms_na_index) y(i) = x(ix(i))
      end do
      if (ms_remove(na_rm)) y = pack(y, .not. ieee_is_nan(y))
   end function

   pure function ms_select(x, rows, cols) result(y)
      real(dp), intent(in) :: x(:, :) !! Source matrix with observations in either dimension.
      integer, intent(in), optional :: rows(:) !! Optional R-style selection of rows; see ms_indices.
      integer, intent(in), optional :: cols(:) !! Optional R-style selection of columns; see ms_indices.
      real(dp), allocatable :: y(:, :)
      integer, allocatable :: ir(:), ic(:)
      integer :: i, j
      ir = ms_indices(size(x, 1), rows)
      ic = ms_indices(size(x, 2), cols)
      allocate(y(size(ir), size(ic)))
      y = ms_nan()
      do j = 1, size(ic)
         if (ic(j) == ms_na_index) cycle
         do i = 1, size(ir)
            if (ir(i) /= ms_na_index) y(i, j) = x(ir(i), ic(j))
         end do
      end do
   end function

   pure recursive real(dp) function ms_sum(x, na_rm) result(y)
      real(dp), intent(in) :: x(:) !! Numeric observations, possibly containing NaN or infinities.
      logical, intent(in), optional :: na_rm !! Ignore NaNs when true; default false.
      real(dp) :: correction, next, value, largest
      integer :: i
      y = 0.0_dp
      correction = 0.0_dp
      do i = 1, size(x)
         if (ieee_is_nan(x(i))) then
            if (ms_remove(na_rm)) cycle
            y = ms_nan()
            return
         end if
         value = x(i)
         if (.not. ieee_is_finite(value) .or. .not. ieee_is_finite(y)) then
            y = y + value
            correction = 0.0_dp
         else
            next = y + value
            if (.not. ieee_is_finite(next)) then
               if (all(ieee_is_finite(x))) then
                  largest = maxval(abs(x))
                  y = ms_sum(x / largest) * largest
                  return
               end if
               y = next
               correction = 0.0_dp
               cycle
            end if
            if (abs(y) >= abs(value)) then
               correction = correction + ((y - next) + value)
            else
               correction = correction + ((value - next) + y)
            end if
            y = next
         end if
      end do
      y = y + correction
   end function

   pure real(dp) function ms_mean(x, na_rm, refine) result(y)
      real(dp), intent(in) :: x(:) !! Numeric sample, possibly containing NaN or infinities.
      logical, intent(in), optional :: na_rm !! Ignore missing observations when true.
      logical, intent(in), optional :: refine !! Apply residual correction when true; default true.
      real(dp), allocatable :: a(:)
      logical :: extra
      real(dp) :: correction
      a = ms_vector(x, na_rm = na_rm)
      y = ms_nan()
      if (size(a) == 0) return
      y = ms_sum(a) / real(size(a), dp)
      if (.not. ieee_is_finite(y) .and. all(ieee_is_finite(a))) y = ms_sum(a / real(size(a), dp))
      extra = .true.
      if (present(refine)) extra = refine
      if (extra .and. ieee_is_finite(y)) then
         correction = ms_sum(a - y) / real(size(a), dp)
         if (ieee_is_finite(correction)) y = y + correction
      end if
   end function

   pure real(dp) function ms_var(x, na_rm, refine, center) result(y)
      real(dp), intent(in) :: x(:) !! Numeric sample; at least two retained observations are needed.
      logical, intent(in), optional :: na_rm !! Remove NaN observations when true.
      logical, intent(in), optional :: refine !! Refine internally estimated mean; default true.
      real(dp), intent(in), optional :: center !! Prespecified sample mean; not independently validated.
      real(dp), allocatable :: a(:)
      real(dp) :: mu
      a = ms_vector(x, na_rm = na_rm)
      y = ms_nan()
      if (size(a) < 2) return
      mu = ms_mean(a, refine = refine)
      if (present(center)) mu = center
      if (.not. ieee_is_finite(mu)) return
      y = ms_sum((a - mu)**2) / real(size(a) - 1, dp)
   end function

   pure real(dp) function ms_product(x, na_rm, method) result(y)
      real(dp), intent(in) :: x(:) !! Numeric factors; empty product equals one.
      logical, intent(in), optional :: na_rm !! Ignore NaN factors when true.
      character(*), intent(in), optional :: method !! 'direct' (default) or sign-preserving 'expSumLog'.
      real(dp), allocatable :: a(:)
      real(dp) :: z, mantissa
      logical :: logarithmic
      integer :: i, power
      a = ms_vector(x, na_rm = na_rm)
      logarithmic = .false.
      if (present(method)) logarithmic = method == 'expSumLog'
      y = 1.0_dp
      if (any(ieee_is_nan(a))) then
         y = ms_nan()
      else if (logarithmic) then
         z = ms_sum(log(abs(a)))
         y = exp(z)
         if (mod(count(a < 0.0_dp), 2) == 1) y = -y
      else if (any(a == 0.0_dp)) then
         y = 0.0_dp
         if (any(.not. ieee_is_finite(a))) y = ms_nan()
      else if (any(.not. ieee_is_finite(a))) then
         y = ms_inf()
         if (mod(count(a < 0.0_dp), 2) == 1) y = -y
      else
         mantissa = 1.0_dp
         power = 0
         do i = 1, size(a)
            mantissa = mantissa * fraction(a(i))
            power = power + exponent(a(i)) + exponent(mantissa)
            mantissa = fraction(mantissa)
         end do
         y = scale(mantissa, power)
      end if
   end function

   pure real(dp) function ms_extreme(x, maximum, na_rm) result(y)
      real(dp), intent(in) :: x(:) !! Sample for a minimum or maximum.
      logical, intent(in) :: maximum !! True requests maximum; false requests minimum.
      logical, intent(in), optional :: na_rm !! Ignore NaNs when true.
      real(dp), allocatable :: a(:)
      a = ms_vector(x, na_rm = na_rm)
      y = ms_inf(maximum)
      if (any(ieee_is_nan(a))) then
         y = ms_nan()
      else if (size(a) > 0) then
         if (maximum) then
            y = maxval(a)
         else
            y = minval(a)
         end if
      end if
   end function

   pure real(dp) function ms_log_sum_exp(x, na_rm) result(y)
      real(dp), intent(in) :: x(:) !! Logarithms of nonnegative terms, with negative infinity allowed.
      logical, intent(in), optional :: na_rm !! Remove NaN log terms when true.
      real(dp), allocatable :: a(:)
      real(dp) :: m
      a = ms_vector(x, na_rm = na_rm)
      y = ms_inf(.true.)
      if (any(ieee_is_nan(a))) then
         y = ms_nan()
         return
      end if
      if (size(a) == 0) return
      m = maxval(a)
      y = m
      if (.not. ieee_is_finite(m)) return
      y = m + log(ms_sum(exp(a - m)))
   end function

   pure real(dp) function ms_count(x, value, na_rm) result(y)
      real(dp), intent(in) :: x(:) !! Values to compare with the target.
      real(dp), intent(in), optional :: value !! Exact target; default one, and NaN counts all missing values.
      logical, intent(in), optional :: na_rm !! Ignore missing comparisons when true.
      real(dp) :: v
      v = 1.0_dp
      if (present(value)) v = value
      if (ieee_is_nan(v)) then
         y = real(count(ieee_is_nan(x)), dp)
      else if (any(ieee_is_nan(x)) .and. .not. ms_remove(na_rm)) then
         y = ms_nan()
      else
         y = real(count(x == v), dp)
      end if
   end function

   pure real(dp) function ms_all(x, value, na_rm) result(y)
      real(dp), intent(in) :: x(:) !! Numeric observations; output is 0, 1, or NaN for unknown.
      real(dp), intent(in), optional :: value !! Explicit numeric target; omitted tests nonzero, NaN tests missingness.
      logical, intent(in), optional :: na_rm !! Ignore NaN observations when true unless target is NaN.
      real(dp) :: v
      v = 1.0_dp
      if (present(value)) v = value
      y = 1.0_dp
      if (.not. present(value)) then
         if (any(x == 0.0_dp)) then
            y = 0.0_dp
         else if (any(ieee_is_nan(x)) .and. .not. ms_remove(na_rm)) then
            y = ms_nan()
         end if
      else if (ieee_is_nan(v)) then
         if (.not. all(ieee_is_nan(x))) y = 0.0_dp
      else if (any(.not. ieee_is_nan(x) .and. x /= v)) then
         y = 0.0_dp
      else if (any(ieee_is_nan(x)) .and. .not. ms_remove(na_rm)) then
         y = ms_nan()
      end if
   end function

   pure real(dp) function ms_any(x, value, na_rm) result(y)
      real(dp), intent(in) :: x(:) !! Numeric observations; output is 0, 1, or NaN for unknown.
      real(dp), intent(in), optional :: value !! Explicit numeric target; omitted tests nonzero, NaN tests missingness.
      logical, intent(in), optional :: na_rm !! Ignore NaN observations when true unless target is NaN.
      real(dp) :: v
      v = 1.0_dp
      if (present(value)) v = value
      y = 0.0_dp
      if (.not. present(value)) then
         if (any(.not. ieee_is_nan(x) .and. x /= 0.0_dp)) then
            y = 1.0_dp
         else if (any(ieee_is_nan(x)) .and. .not. ms_remove(na_rm)) then
            y = ms_nan()
         end if
      else if (ieee_is_nan(v)) then
         if (any(ieee_is_nan(x))) y = 1.0_dp
      else if (any(.not. ieee_is_nan(x) .and. x == v)) then
         y = 1.0_dp
      else if (any(ieee_is_nan(x)) .and. .not. ms_remove(na_rm)) then
         y = ms_nan()
      end if
   end function

   pure function ms_difference(x, lag, differences) result(y)
      real(dp), intent(in) :: x(:) !! Ordered observations to difference.
      integer, intent(in), optional :: lag !! Positive index separation; default one.
      integer, intent(in), optional :: differences !! Nonnegative number of repeated differences; default one.
      real(dp), allocatable :: y(:)
      integer :: l, d, i, n
      l = 1
      d = 1
      if (present(lag)) l = lag
      if (present(differences)) d = differences
      if (l < 1 .or. d < 0) then
         allocate(y(0))
         return
      end if
      y = x
      do i = 1, d
         n = max(0, size(y) - l)
         if (n == 0) then
            y = [real(dp) ::]
            exit
         end if
         y = y(l + 1:l + n) - y(1:n)
      end do
   end function

   pure function ms_cumulative(x, operation) result(y)
      real(dp), intent(in) :: x(:) !! Ordered observations; missing values propagate from first occurrence.
      character(*), intent(in) :: operation !! One of 'sum', 'prod', 'min', or 'max'.
      real(dp) :: y(size(x))
      integer :: i
      if (size(x) == 0) return
      y(1) = x(1)
      do i = 2, size(x)
         if (ieee_is_nan(y(i - 1)) .or. ieee_is_nan(x(i))) then
            y(i) = ms_nan()
            cycle
         end if
         select case (operation)
         case ('sum')
            y(i) = y(i - 1) + x(i)
         case ('prod')
            y(i) = y(i - 1) * x(i)
         case ('min')
            y(i) = min(y(i - 1), x(i))
         case ('max')
            y(i) = max(y(i - 1), x(i))
         case default
            y(i) = ms_nan()
         end select
      end do
   end function
end module matrixstats_common
