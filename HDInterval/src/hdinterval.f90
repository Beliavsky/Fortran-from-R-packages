! SPDX-License-Identifier: GPL-3.0-only
! Numerical translation of HDInterval 0.2.4 by Mike Meredith and John Kruschke.
! See NOTICE.md for upstream attribution and translation provenance.
module hdinterval
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_value, ieee_quiet_nan, ieee_positive_inf
   use r_kinds, only: dp
   use r_sorting, only: r_sort
   implicit none
   private

   integer, parameter, public :: hdi_ok = 0
   integer, parameter, public :: hdi_invalid_mass = 1
   integer, parameter, public :: hdi_invalid_data = 2
   integer, parameter, public :: hdi_invalid_callback = 3
   integer, parameter, public :: hdi_no_interval = 4
   integer, parameter, public :: hdi_no_convergence = 5
   integer, parameter, public :: hdi_invalid_options = 6
   integer, parameter, public :: hdi_unassociated_callback = 7
   public :: dp, scalar_function_interface
   public :: hdi_samples, hdi_columns, hdi_density, hdi_function
   public :: inverse_cdf, inverse_cdf_scalar, inverse_cdf_vector

   abstract interface
      pure function scalar_function_interface(x) result(value)
         import dp
         real(dp), intent(in) :: x !! Probability for a quantile callback, or coordinate for a CDF callback.
         real(dp) :: value
      end function scalar_function_interface
   end interface

   interface inverse_cdf
      !! Numerically invert a continuous CDF at one probability or a vector of probabilities.
      module procedure inverse_cdf_scalar, inverse_cdf_vector
   end interface inverse_cdf

contains

   pure elemental logical function valid_mass(mass) result(valid)
      real(dp), intent(in) :: mass !! Requested probability mass, strictly between zero and one.
      valid = .false.
      if (.not. ieee_is_finite(mass)) return
      valid = mass > 0.0_dp .and. mass < 1.0_dp
   end function valid_mass

   pure elemental function nonnegative_difference(lower, upper) result(width)
      !! Subtract ordered endpoints without raising floating-point overflow.
      real(dp), intent(in) :: lower !! Non-NaN lower coordinate; may be negative infinity.
      real(dp), intent(in) :: upper !! Non-NaN upper coordinate at least as large as lower.
      real(dp) :: width

      width = ieee_value(0.0_dp, ieee_positive_inf)
      if (.not. ieee_is_finite(lower) .or. .not. ieee_is_finite(upper)) return
      if (lower < 0.0_dp) then
         if (upper > huge(1.0_dp) + lower) return
      end if
      width = upper - lower
   end function nonnegative_difference

   pure subroutine hdi_samples(samples, interval, cred_mass, status)
      !! Shortest ordered-sample interval; removes NaNs, retains infinities, and resolves ties first.
      real(dp), intent(in) :: samples(:) !! Draws in arbitrary order; NaNs are omitted and infinities retained.
      real(dp), intent(out) :: interval(2) !! Lower and upper endpoints; NaNs if no interval can be obtained.
      real(dp), intent(in), optional :: cred_mass !! Probability mass in (0,1); defaults to 0.95.
      integer, intent(out), optional :: status !! Zero for success, invalid mass, or no interval for an empty sample.
      real(dp), allocatable :: ordered(:)
      real(dp) :: mass, width, best_width
      integer :: n, span, i, best

      interval = ieee_value(0.0_dp, ieee_quiet_nan)
      if (present(status)) status = hdi_invalid_mass
      mass = 0.95_dp
      if (present(cred_mass)) mass = cred_mass
      if (.not. valid_mass(mass)) return
      call r_sort(samples, ordered, na_rm = .true.)
      n = size(ordered)
      if (present(status)) status = hdi_no_interval
      if (n == 0) return
      span = min(n - 1, floor(real(n, dp) * mass))
      best = 0
      best_width = huge(1.0_dp)
      do i = 1, n - span
         ! R ignores NaN widths from Inf - Inf without discarding the infinite draws.
         if (.not. ieee_is_finite(ordered(i)) .and. ordered(i) == ordered(i + span)) cycle
         width = nonnegative_difference(ordered(i), ordered(i + span))
         if (best == 0) then
            best = i
            best_width = width
         else if (width < best_width) then
            best = i
            best_width = width
         end if
      end do
      if (best == 0) then
         interval = [ordered(1), ordered(n)]
      else
         interval = [ordered(best), ordered(best + span)]
      end if
      if (present(status)) status = hdi_ok
   end subroutine hdi_samples

   pure subroutine hdi_columns(samples, intervals, cred_mass, statuses)
      !! Applies the sample HDI independently to matrix columns, as hdi.matrix does.
      real(dp), intent(in) :: samples(:, :) !! Matrix with draws in rows and variables in columns.
      real(dp), intent(out) :: intervals(2, size(samples, 2)) !! Lower and upper endpoints for each variable.
      real(dp), intent(in), optional :: cred_mass !! Shared probability mass in (0,1); defaults to 0.95.
      integer, intent(out), optional :: statuses(size(samples, 2)) !! Per-column status codes, including empty columns.
      integer :: j, code

      do j = 1, size(samples, 2)
         call hdi_samples(samples(:, j), intervals(:, j), cred_mass, code)
         if (present(statuses)) statuses(j) = code
      end do
   end subroutine hdi_columns

   pure subroutine hdi_density(x, y, intervals, height, cred_mass, allow_split, status, is_hdi, achieved_mass)
      !! Threshold an equally spaced density grid, optionally retaining disjoint intervals.
      real(dp), intent(in) :: x(:) !! At least two finite, increasing, equally spaced density-grid coordinates.
      real(dp), intent(in) :: y(:) !! Finite nonnegative density ordinates matching x; need not integrate to one.
      real(dp), allocatable, intent(out) :: intervals(:, :) !! Shape (2,k) endpoints, or (2,0) on failure.
      real(dp), intent(out) :: height !! Threshold in the original density units; NaN on invalid input.
      real(dp), intent(in), optional :: cred_mass !! Desired discrete grid mass in (0,1); defaults to 0.95.
      logical, intent(in), optional :: allow_split !! Retain separate modal intervals when true; defaults to false.
      integer, intent(out), optional :: status !! Zero for success, invalid mass/data, or no fallback interval.
      logical, intent(out), optional :: is_hdi !! False if a contiguous credible-interval fallback is used or on error.
      real(dp), intent(out), optional :: achieved_mass !! Selected grid mass; fallback uses the cumulative mass difference.
      real(dp), allocatable :: weights(:), sorted(:), cumulative(:)
      logical, allocatable :: selected(:)
      real(dp) :: mass, ymax, total, partial, threshold, dx, spacing_tolerance, target
      real(dp) :: lower_sum, upper_sum, selected_total
      logical :: split
      integer :: n, i, j, groups, k, best_width, ties, left, right, middle, width

      allocate(intervals(2, 0))
      height = ieee_value(0.0_dp, ieee_quiet_nan)
      if (present(achieved_mass)) achieved_mass = height
      if (present(is_hdi)) is_hdi = .false.
      if (present(status)) status = hdi_invalid_mass
      mass = 0.95_dp
      if (present(cred_mass)) mass = cred_mass
      if (.not. valid_mass(mass)) return
      if (present(status)) status = hdi_invalid_data
      n = size(x)
      if (n < 2 .or. size(y) /= n) return
      if (any(.not. ieee_is_finite(x)) .or. any(.not. ieee_is_finite(y))) return
      if (any(y < 0.0_dp)) return
      ymax = maxval(y)
      if (ymax <= 0.0_dp) return
      if (x(2) <= x(1)) return
      dx = nonnegative_difference(x(1), x(2))
      if (.not. ieee_is_finite(dx)) return
      if (dx <= 0.0_dp) return
      spacing_tolerance = 128.0_dp * epsilon(1.0_dp) * max(maxval(abs(x)), dx)
      do i = 2, n
         if (x(i) <= x(i - 1)) return
         if (.not. ieee_is_finite(nonnegative_difference(x(i - 1), x(i)))) return
         if (abs((x(i) - x(i - 1)) - dx) > spacing_tolerance) return
      end do
      ! Scaling prevents overflow and makes the calculation invariant to density normalization.
      weights = y / ymax
      total = sum(weights)
      call r_sort(weights, sorted)
      partial = 0.0_dp
      threshold = sorted(1)
      do i = n, 1, -1
         partial = partial + sorted(i)
         threshold = sorted(i)
         if (partial >= total * mass) exit
      end do
      height = threshold * ymax
      selected = weights >= threshold
      groups = 0
      do i = 1, n
         if (.not. selected(i)) cycle
         if (i == 1) then
            groups = groups + 1
         else if (.not. selected(i - 1)) then
            groups = groups + 1
         end if
      end do
      split = .false.
      if (present(allow_split)) split = allow_split
      if (groups > 1 .and. .not. split) then
         allocate(cumulative(n))
         partial = 0.0_dp
         do i = 1, n
            partial = partial + weights(i)
            cumulative(i) = partial / total
         end do
         best_width = n
         ties = 0
         lower_sum = 0.0_dp
         upper_sum = 0.0_dp
         selected_total = 0.0_dp
         do i = 1, n
            if (cumulative(i) >= 1.0_dp - mass) exit
            target = cumulative(i) + mass
            left = i + 1
            right = n
            do while (left < right)
               middle = left + (right - left) / 2
               if (cumulative(middle) > target) then
                  right = middle
               else
                  left = middle + 1
               end if
            end do
            j = left
            if (j > n) cycle
            if (cumulative(j) <= target) cycle
            width = j - i
            if (width < best_width) then
               best_width = width
               ties = 0
               lower_sum = 0.0_dp
               upper_sum = 0.0_dp
               selected_total = 0.0_dp
            end if
            if (width == best_width) then
               ties = ties + 1
               ! Online means avoid overflow from summing large coordinates.
               lower_sum = lower_sum * (real(ties - 1, dp) / real(ties, dp)) + x(i) / real(ties, dp)
               upper_sum = upper_sum * (real(ties - 1, dp) / real(ties, dp)) + x(j) / real(ties, dp)
               selected_total = selected_total + (cumulative(j) - cumulative(i))
            end if
         end do
         if (ties == 0) then
            if (present(status)) status = hdi_no_interval
            return
         end if
         deallocate(intervals)
         allocate(intervals(2, 1))
         intervals(:, 1) = [lower_sum, upper_sum]
         if (present(achieved_mass)) achieved_mass = selected_total / real(ties, dp)
      else
         deallocate(intervals)
         allocate(intervals(2, groups))
         k = 0
         i = 1
         do while (i <= n)
            if (.not. selected(i)) then
               i = i + 1
               cycle
            end if
            k = k + 1
            intervals(1, k) = x(i)
            j = i
            do while (j < n)
               if (.not. selected(j + 1)) exit
               j = j + 1
            end do
            intervals(2, k) = x(j)
            i = j + 1
         end do
         if (present(is_hdi)) is_hdi = .true.
         if (present(achieved_mass)) achieved_mass = sum(weights, mask = selected) / total
      end if
      if (present(status)) status = hdi_ok
   end subroutine hdi_density

   pure subroutine hdi_function(icdf, interval, cred_mass, tol, status, lower_tail, max_iterations)
      !! Minimize inverse-CDF interval width by bounded golden-section search.
      procedure(scalar_function_interface), pointer, intent(in) :: icdf !! Associated pure, nondecreasing quantile callback.
      real(dp), intent(out) :: interval(2) !! Finite lower and upper HDI endpoints, or NaNs on failure.
      real(dp), intent(in), optional :: cred_mass !! Probability mass in (0,1); defaults to 0.95.
      real(dp), intent(in), optional :: tol !! Positive finite absolute tolerance on tail probability; defaults to 1e-8.
      integer, intent(out), optional :: status !! Zero on success, otherwise one of the public hdi status constants.
      real(dp), intent(out), optional :: lower_tail !! Optimized lower-tail probability; NaN on failure.
      integer, intent(in), optional :: max_iterations !! Positive golden-section iteration budget; defaults to 200.
      real(dp), parameter :: ratio = 0.6180339887498948482_dp
      real(dp) :: mass, tolerance, median, a, b, c, d, fc, fd, location
      integer :: limit, iteration

      interval = ieee_value(0.0_dp, ieee_quiet_nan)
      if (present(lower_tail)) lower_tail = interval(1)
      if (present(status)) status = hdi_invalid_mass
      mass = 0.95_dp
      if (present(cred_mass)) mass = cred_mass
      if (.not. valid_mass(mass)) return
      if (present(status)) status = hdi_invalid_options
      tolerance = 1e-8_dp
      if (present(tol)) tolerance = tol
      limit = 200
      if (present(max_iterations)) limit = max_iterations
      if (.not. ieee_is_finite(tolerance)) return
      if (tolerance <= 0.0_dp .or. limit <= 0) return
      if (present(status)) status = hdi_unassociated_callback
      if (.not. associated(icdf)) return
      if (present(status)) status = hdi_invalid_callback
      median = icdf(0.5_dp)
      if (.not. ieee_is_finite(median)) return
      a = 0.0_dp
      b = 1.0_dp - mass
      c = b - ratio * (b - a)
      d = a + ratio * (b - a)
      fc = interval_width(c)
      fd = interval_width(d)
      do iteration = 1, limit
         if (.not. ieee_is_finite(fc) .or. .not. ieee_is_finite(fd)) return
         if (b - a <= max(tolerance, 4.0_dp * epsilon(1.0_dp))) exit
         if (fc <= fd) then
            b = d
            d = c
            fd = fc
            c = b - ratio * (b - a)
            fc = interval_width(c)
         else
            a = c
            c = d
            fc = fd
            d = a + ratio * (b - a)
            fd = interval_width(d)
         end if
      end do
      if (iteration > limit) then
         if (present(status)) status = hdi_no_convergence
         return
      end if
      location = 0.5_dp * (a + b)
      interval = [icdf(location), icdf(location + mass)]
      if (any(.not. ieee_is_finite(interval)) .or. interval(1) > interval(2)) then
         interval = ieee_value(0.0_dp, ieee_quiet_nan)
         return
      end if
      if (present(lower_tail)) lower_tail = location
      if (present(status)) status = hdi_ok

   contains

      pure function interval_width(probability) result(width)
         real(dp), intent(in) :: probability !! Lower-tail probability probed by the bounded golden-section search.
         real(dp) :: width, lo, hi

         lo = icdf(probability)
         hi = icdf(probability + mass)
         width = ieee_value(0.0_dp, ieee_quiet_nan)
         if (.not. ieee_is_finite(lo) .or. .not. ieee_is_finite(hi)) return
         if (lo > hi) return
         width = nonnegative_difference(lo, hi)
      end function interval_width

   end subroutine hdi_function

   pure subroutine inverse_cdf_scalar(p, cdf, quantile, tol, lower, upper, status, max_iterations, max_extensions)
      !! Invert a CDF after applying HDInterval's [1e-6, 1-1e-6] probability clipping.
      real(dp), intent(in) :: p !! Probability in [0,1]; endpoints are clipped, and NaNs are rejected.
      procedure(scalar_function_interface), pointer, intent(in) :: cdf !! Associated pure continuous nondecreasing CDF callback.
      real(dp), intent(out) :: quantile !! Finite root estimate, or NaN if invalid or not converged.
      real(dp), intent(in), optional :: tol !! Positive finite absolute x tolerance; defaults to epsilon(dp)**0.25.
      real(dp), intent(in), optional :: lower !! Finite initial lower bracket; defaults to -10.
      real(dp), intent(in), optional :: upper !! Finite initial upper bracket greater than lower; defaults to 10.
      integer, intent(out), optional :: status !! Zero on success, otherwise a public hdi status code.
      integer, intent(in), optional :: max_iterations !! Positive bisection iteration budget; defaults to 1000.
      integer, intent(in), optional :: max_extensions !! Nonnegative symmetric bracket-extension budget; defaults to 100.
      real(dp) :: probability, tolerance, a, b, fa, fb, mid, fm, step, half_width
      integer :: limit, extensions, iteration

      quantile = ieee_value(0.0_dp, ieee_quiet_nan)
      if (present(status)) status = hdi_invalid_data
      if (.not. ieee_is_finite(p)) return
      if (p < 0.0_dp .or. p > 1.0_dp) return
      probability = max(1e-6_dp, min(1.0_dp - 1e-6_dp, p))
      tolerance = epsilon(1.0_dp) ** 0.25_dp
      if (present(tol)) tolerance = tol
      a = -10.0_dp
      if (present(lower)) a = lower
      b = 10.0_dp
      if (present(upper)) b = upper
      limit = 1000
      if (present(max_iterations)) limit = max_iterations
      extensions = 100
      if (present(max_extensions)) extensions = max_extensions
      if (present(status)) status = hdi_invalid_options
      if (.not. ieee_is_finite(tolerance) .or. .not. ieee_is_finite(a) .or. .not. ieee_is_finite(b)) return
      if (tolerance <= 0.0_dp .or. a >= b .or. limit <= 0 .or. extensions < 0) return
      if (present(status)) status = hdi_unassociated_callback
      if (.not. associated(cdf)) return
      if (present(status)) status = hdi_invalid_callback
      fa = checked_cdf(a)
      fb = checked_cdf(b)
      if (.not. ieee_is_finite(fa) .or. .not. ieee_is_finite(fb)) return
      if (fa > fb) return
      fa = fa - probability
      fb = fb - probability
      step = max(1.0_dp, abs(0.5_dp * b - 0.5_dp * a))
      do iteration = 0, extensions
         if (fa == 0.0_dp) then
            quantile = a
            if (present(status)) status = hdi_ok
            return
         else if (fb == 0.0_dp) then
            quantile = b
            if (present(status)) status = hdi_ok
            return
         end if
         if (fa < 0.0_dp .and. fb > 0.0_dp) exit
         if (iteration == extensions) then
            if (present(status)) status = hdi_no_convergence
            return
         end if
         a = max(-huge(1.0_dp) + step, a) - step
         b = min(huge(1.0_dp) - step, b) + step
         fa = checked_cdf(a) - probability
         fb = checked_cdf(b) - probability
         if (.not. ieee_is_finite(fa) .or. .not. ieee_is_finite(fb)) return
         if (fa > fb) return
         step = min(0.5_dp * huge(1.0_dp), step) * 2.0_dp
      end do
      do iteration = 1, limit
         mid = 0.5_dp * a + 0.5_dp * b
         half_width = 0.5_dp * b - 0.5_dp * a
         fm = checked_cdf(mid) - probability
         if (.not. ieee_is_finite(fm)) return
         if (fm < fa .or. fm > fb) return
         if (fm == 0.0_dp .or. half_width <= tolerance .or. mid == a .or. mid == b) then
            quantile = mid
            if (present(status)) status = hdi_ok
            return
         end if
         if (fm < 0.0_dp) then
            a = mid
            fa = fm
         else
            b = mid
            fb = fm
         end if
      end do
      if (present(status)) status = hdi_no_convergence

   contains

      pure function checked_cdf(coordinate) result(value)
         real(dp), intent(in) :: coordinate !! Real coordinate at which a probability-valued CDF is required.
         real(dp) :: value

         value = ieee_value(0.0_dp, ieee_quiet_nan)
         if (.not. ieee_is_finite(coordinate)) return
         value = cdf(coordinate)
         if (.not. ieee_is_finite(value)) return
         if (value < 0.0_dp .or. value > 1.0_dp) value = ieee_value(0.0_dp, ieee_quiet_nan)
      end function checked_cdf

   end subroutine inverse_cdf_scalar

   pure subroutine inverse_cdf_vector(p, cdf, quantiles, tol, lower, upper, status, max_iterations, &
                                     max_extensions, element_status)
      !! Invert each probability independently, retaining successful elements on partial failure.
      real(dp), intent(in) :: p(:) !! Probabilities in [0,1]; endpoints receive upstream clipping.
      procedure(scalar_function_interface), pointer, intent(in) :: cdf !! Associated pure continuous nondecreasing CDF callback.
      real(dp), intent(out) :: quantiles(size(p)) !! Root estimates in p order; failed elements are NaNs.
      real(dp), intent(in), optional :: tol !! Positive finite absolute root tolerance; defaults to epsilon(dp)**0.25.
      real(dp), intent(in), optional :: lower !! Finite initial lower bracket for every root; defaults to -10.
      real(dp), intent(in), optional :: upper !! Finite initial upper bracket for every root; defaults to 10.
      integer, intent(out), optional :: status !! First nonzero element status, or zero including an empty vector.
      integer, intent(in), optional :: max_iterations !! Positive per-root bisection budget; defaults to 1000.
      integer, intent(in), optional :: max_extensions !! Nonnegative per-root symmetric expansion budget; defaults to 100.
      integer, intent(out), optional :: element_status(size(p)) !! Independent status code for each supplied probability.
      integer :: i, code, first_error

      first_error = hdi_ok
      do i = 1, size(p)
         call inverse_cdf_scalar(p(i), cdf, quantiles(i), tol, lower, upper, code, max_iterations, max_extensions)
         if (present(element_status)) element_status(i) = code
         if (first_error == hdi_ok .and. code /= hdi_ok) first_error = code
      end do
      if (present(status)) status = first_error
   end subroutine inverse_cdf_vector

end module hdinterval
