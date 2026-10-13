! SPDX-License-Identifier: GPL-2.0-or-later
! Corrected AS 217, translated from diptest 0.77-2 src/dip.c and R/*.R.
! P. M. Hartigan (1985); Dario Ringach; Martin Maechler, ETH Zurich.
! Retains Maechler's 1994 termination fix and Yong Lu's 2003 correction.
! See NOTICE.md and upstream/src/dip.c for the complete historical notice.
module diptest_api
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan, ieee_value, ieee_quiet_nan
   use r_kinds, only : dp
   use r_sorting, only : r_sort
   use diptest_table, only : table_sizes, table_probabilities, table_quantiles
   implicit none
   private
   public :: dp, dip, dip_compute, dip_test, dip_test_simulate, dip_p_value
   public :: dip_result, dip_test_result
   public :: dip_ok, dip_empty_sample, dip_unsorted, dip_nonfinite, dip_invalid_simulation, dip_invalid_argument

   integer, parameter :: dip_ok = 0
   integer, parameter :: dip_empty_sample = 1
   integer, parameter :: dip_unsorted = 2
   integer, parameter :: dip_nonfinite = 3
   integer, parameter :: dip_invalid_simulation = 4
   integer, parameter :: dip_invalid_argument = 5

   type :: dip_result
      real(dp) :: statistic = 0.0_dp
      real(dp) :: lower = 0.0_dp
      real(dp) :: upper = 0.0_dp
      integer :: n = 0
      integer :: status = dip_ok
      integer :: low = 0
      integer :: high = 0
      real(dp), allocatable :: x(:)
      integer, allocatable :: gcm(:), lcm(:), mn(:), mj(:), full_gcm(:), full_lcm(:)
   end type dip_result

   type :: dip_test_result
      real(dp) :: statistic = 0.0_dp
      real(dp) :: p_value = 0.0_dp
      integer :: n = 0
      integer :: status = dip_ok
      integer :: replicates = 0
      logical :: simulated = .false.
      logical :: asymptotic = .false.
   end type dip_test_result

contains

   pure function dip(x, min_is_zero) result(statistic)
      real(dp), intent(in) :: x(:) !! Finite observations in arbitrary order; an empty or nonfinite sample returns NaN.
      logical, intent(in), optional :: min_is_zero !! If true allow a zero minimum; default false uses Hartigan's 1/(2*n) floor.
      real(dp) :: statistic
      type(dip_result) :: fit
      call dip_compute(x, fit, min_is_zero)
      statistic = fit%statistic
   end function dip

   pure subroutine dip_compute(x, result, min_is_zero, sorted)
      real(dp), intent(in) :: x(:) !! Finite observations; copied and sorted unless sorted is explicitly true.
      type(dip_result), intent(out) :: result !! Statistic, sorted data, 1-based modal endpoints, both hulls and predecessor arrays.
      logical, intent(in), optional :: min_is_zero !! Allow zero dip when true; default false retains the original 1/(2*n) floor.
      logical, intent(in), optional :: sorted !! If true require ascending input and return dip_unsorted on a decreasing pair.
      real(dp), allocatable :: work(:)
      logical :: zero, ordered
      integer :: n

      result%statistic = ieee_value(0.0_dp, ieee_quiet_nan)
      result%lower = result%statistic
      result%upper = result%statistic
      n = size(x)
      result%n = n
      zero = .false.
      if (present(min_is_zero)) zero = min_is_zero
      ordered = .false.
      if (present(sorted)) ordered = sorted
      if (n == 0) then
         result%status = dip_empty_sample
         return
      end if
      if (any(.not. ieee_is_finite(x))) then
         result%status = dip_nonfinite
         return
      end if
      if (ordered) then
         if (any(x(2:n) < x(1:n - 1))) then
            result%status = dip_unsorted
            return
         end if
         result%x = x
      else
         call r_sort(x, result%x)
      end if
      allocate(result%mn(n), result%mj(n))
      result%low = 1
      result%high = n
      if (n == 1 .or. result%x(n) == result%x(1)) then
         result%statistic = 0.5_dp / real(n, dp)
         if (zero) result%statistic = 0.0_dp
         result%mn = 1
         result%mj = n
         if (n == 1) then
            result%gcm = [1]
            result%lcm = [1]
         else
            result%gcm = [n, 1]
            result%lcm = [1, n]
         end if
         result%full_gcm = result%gcm
         result%full_lcm = result%lcm
      else
         work = result%x
         call dip_kernel(work, zero, result)
      end if
      result%lower = result%x(result%low)
      result%upper = result%x(result%high)
   end subroutine dip_compute

   pure subroutine dip_kernel(x, zero, result)
      real(dp), intent(in) :: x(:) !! Ascending, nonconstant finite working data; local scaled differences preserve extreme-range clusters.
      logical, intent(in) :: zero !! Whether the statistic may start at zero rather than the classical scaled value one.
      type(dip_result), intent(inout) :: result !! Allocated predecessor arrays and final statistic, modal interval and hull nodes.
      integer :: n, j, k, mnj, mnmnj, mjk, mjmjk, low, high, ng, nl, i, ig, ih, ix, iv
      integer :: gcmi, lcmi, previous, jb, je, jj
      integer :: gcm(size(x)), lcm(size(x))
      real(dp) :: current, distance, delta, dip_lower, dip_upper, largest, candidate

      n = size(x)
      low = 1
      high = n
      current = 1.0_dp
      if (zero) current = 0.0_dp
      result%mn(1) = 1
      do j = 2, n
         result%mn(j) = j - 1
         do
            mnj = result%mn(j)
            mnmnj = result%mn(mnj)
            if (mnj == 1) exit
            if (slope_less(x(j), x(mnj), mnj - mnmnj, x(mnj), x(mnmnj), j - mnj)) exit
            result%mn(j) = mnmnj
         end do
      end do
      result%mj(n) = n
      do k = n - 1, 1, -1
         result%mj(k) = k + 1
         do
            mjk = result%mj(k)
            mjmjk = result%mj(mjk)
            if (mjk == n) exit
            if (slope_less(x(mjk), x(k), mjmjk - mjk, x(mjmjk), x(mjk), mjk - k)) exit
            result%mj(k) = mjmjk
         end do
      end do
      call collect_hulls(result%mn, result%mj, 1, n, gcm, lcm, ng, nl)
      result%full_gcm = gcm(:ng)
      result%full_lcm = lcm(:nl)

      do
         call collect_hulls(result%mn, result%mj, low, high, gcm, lcm, ng, nl)
         ig = ng
         ix = ig - 1
         ih = nl
         iv = 2
         distance = 0.0_dp
         if (ng /= 2 .or. nl /= 2) then
            do
               gcmi = gcm(ix)
               lcmi = lcm(iv)
               if (gcmi > lcmi) then
                  previous = gcm(ix + 1)
                  delta = real(lcmi - previous + 1, dp) - &
                     interval_fraction(x(lcmi), x(previous), x(gcmi)) * real(gcmi - previous, dp)
                  iv = iv + 1
                  if (delta >= distance) then
                     distance = delta
                     ig = ix + 1
                     ih = iv - 1
                  end if
               else
                  previous = lcm(iv - 1)
                  ! Yong Lu's corrected parentheses: subtraction is outside the quotient.
                  delta = interval_fraction(x(gcmi), x(previous), x(lcmi)) * real(lcmi - previous, dp) &
                     - real(gcmi - previous - 1, dp)
                  ix = ix - 1
                  if (delta >= distance) then
                     distance = delta
                     ig = ix + 1
                     ih = iv
                  end if
               end if
               ix = max(1, ix)
               iv = min(nl, iv)
               if (gcm(ix) == lcm(iv)) exit
            end do
         else
            distance = 1.0_dp
            if (zero) distance = 0.0_dp
         end if
         if (distance < current) exit

         dip_lower = 0.0_dp
         do j = ig, ng - 1
            largest = 1.0_dp
            jb = gcm(j + 1)
            je = gcm(j)
            if (je - jb > 1 .and. x(je) /= x(jb)) then
               do jj = jb, je
                  candidate = real(jj - jb + 1, dp) - &
                     interval_fraction(x(jj), x(jb), x(je)) * real(je - jb, dp)
                  largest = max(largest, candidate)
               end do
            end if
            dip_lower = max(dip_lower, largest)
         end do
         dip_upper = 0.0_dp
         do j = ih, nl - 1
            largest = 1.0_dp
            jb = lcm(j)
            je = lcm(j + 1)
            if (je - jb > 1 .and. x(je) /= x(jb)) then
               do jj = jb, je
                  candidate = interval_fraction(x(jj), x(jb), x(je)) * real(je - jb, dp) &
                     - real(jj - jb - 1, dp)
                  largest = max(largest, candidate)
               end do
            end if
            dip_upper = max(dip_upper, largest)
         end do
         current = max(current, dip_lower, dip_upper)
         ! Maechler's 1994 termination fix is essential for unimodal inputs.
         if (low == gcm(ig) .and. high == lcm(ih)) exit
         low = gcm(ig)
         high = lcm(ih)
      end do
      result%statistic = current / (2.0_dp * real(n, dp))
      result%low = low
      result%high = high
      result%gcm = gcm(:ng)
      result%lcm = lcm(:nl)
   end subroutine dip_kernel

   pure elemental subroutine positive_difference(a, b, part, power)
      real(dp), intent(in) :: a !! Finite upper coordinate, at least b.
      real(dp), intent(in) :: b !! Finite lower coordinate; either sign and subnormal values are supported.
      real(dp), intent(out) :: part !! Nonnegative normalized binary fraction of a-b, or zero when a equals b.
      integer, intent(out) :: power !! Binary exponent such that a-b equals part times 2**power without forming that product.
      real(dp) :: difference
      integer :: adjustment
      power = exponent(max(abs(a), abs(b)))
      difference = scale(a, -power) - scale(b, -power)
      adjustment = exponent(difference)
      part = fraction(difference)
      power = power + adjustment
   end subroutine positive_difference

   pure elemental function interval_fraction(y, a, b) result(ratio)
      real(dp), intent(in) :: y !! Coordinate inside the finite closed interval [a,b].
      real(dp), intent(in) :: a !! Lower interval endpoint, strictly less than b.
      real(dp), intent(in) :: b !! Upper interval endpoint; subtraction from a need not fit in dp.
      real(dp) :: ratio, numerator, denominator
      integer :: numerator_power, denominator_power
      call positive_difference(y, a, numerator, numerator_power)
      call positive_difference(b, a, denominator, denominator_power)
      ratio = scale(numerator / denominator, numerator_power - denominator_power)
   end function interval_fraction

   pure elemental function slope_less(a, b, left_weight, c, d, right_weight) result(less)
      real(dp), intent(in) :: a !! Upper endpoint of the first finite coordinate difference.
      real(dp), intent(in) :: b !! Lower endpoint of the first coordinate difference.
      integer, intent(in) :: left_weight !! Positive index separation multiplying the first difference.
      real(dp), intent(in) :: c !! Upper endpoint of the second finite coordinate difference.
      real(dp), intent(in) :: d !! Lower endpoint of the second coordinate difference.
      integer, intent(in) :: right_weight !! Positive index separation multiplying the second difference.
      logical :: less
      real(dp) :: left, right
      integer :: left_power, right_power, adjustment
      call positive_difference(a, b, left, left_power)
      call positive_difference(c, d, right, right_power)
      if (left == 0.0_dp .or. right == 0.0_dp) then
         less = left < right
         return
      end if
      left = left * real(left_weight, dp)
      right = right * real(right_weight, dp)
      adjustment = exponent(left)
      left_power = left_power + adjustment
      left = fraction(left)
      adjustment = exponent(right)
      right_power = right_power + adjustment
      right = fraction(right)
      if (left_power == right_power) then
         less = left < right
      else
         less = left_power < right_power
      end if
   end function slope_less

   pure subroutine collect_hulls(mn, mj, low, high, gcm, lcm, ng, nl)
      integer, intent(in) :: mn(:) !! Convex-minorant predecessor indices over the full sorted sample.
      integer, intent(in) :: mj(:) !! Concave-majorant successor indices over the full sorted sample.
      integer, intent(in) :: low !! Inclusive lower 1-based endpoint of the current modal interval.
      integer, intent(in) :: high !! Inclusive upper 1-based endpoint of the current modal interval.
      integer, intent(out) :: gcm(:) !! Workspace of at least sample size, filled from high to low.
      integer, intent(out) :: lcm(:) !! Workspace of at least sample size, filled from low to high.
      integer, intent(out) :: ng !! Number of meaningful convex-minorant nodes in gcm.
      integer, intent(out) :: nl !! Number of meaningful concave-majorant nodes in lcm.
      ng = 1
      gcm(1) = high
      do while (gcm(ng) > low)
         gcm(ng + 1) = mn(gcm(ng))
         ng = ng + 1
      end do
      nl = 1
      lcm(1) = low
      do while (lcm(nl) < high)
         lcm(nl + 1) = mj(lcm(nl))
         nl = nl + 1
      end do
   end subroutine collect_hulls

   pure elemental subroutine dip_p_value(statistic, n, p_value, status, asymptotic)
      real(dp), intent(in) :: statistic !! Observed dip in [0, 0.5]; this wider interval includes the n=1 classical floor.
      integer, intent(in) :: n !! Positive number of observations after removal of missing values.
      real(dp), intent(out) :: p_value !! Table-interpolated upper-tail probability, or NaN on invalid arguments.
      integer, intent(out), optional :: status !! dip_ok on success, otherwise dip_invalid_argument.
      logical, intent(out), optional :: asymptotic !! True when n exceeds 72000 and the final sqrt(n)-scaled table row is used.
      real(dp) :: weight, scaled_dip, probabilities(size(table_probabilities)), knots(size(table_probabilities))
      real(dp) :: ordinate
      integer :: lower, upper, j, k, m

      p_value = ieee_value(0.0_dp, ieee_quiet_nan)
      if (present(status)) status = dip_invalid_argument
      if (present(asymptotic)) asymptotic = .false.
      if (n <= 0 .or. .not. ieee_is_finite(statistic)) return
      if (statistic < 0.0_dp .or. statistic > 0.5_dp) return
      if (present(status)) status = dip_ok
      if (n <= 3) then
         p_value = 1.0_dp
         return
      end if
      lower = size(table_sizes)
      upper = lower
      weight = 0.0_dp
      if (n < table_sizes(lower)) then
         lower = 1
         do while (table_sizes(lower + 1) <= n)
            lower = lower + 1
         end do
         upper = lower + 1
         weight = real(n - table_sizes(lower), dp) / real(table_sizes(upper) - table_sizes(lower), dp)
      end if
      if (present(asymptotic)) asymptotic = n > table_sizes(size(table_sizes))
      knots = sqrt(real(table_sizes(lower), dp)) * table_quantiles(lower, :)
      knots = knots + weight * (sqrt(real(table_sizes(upper), dp)) * table_quantiles(upper, :) - knots)
      probabilities = table_probabilities
      ! stats::approx defaults to ties=mean: repeated quantile knots must be collapsed.
      m = 0
      j = 1
      do while (j <= size(knots))
         k = j
         do while (k < size(knots))
            if (knots(k + 1) /= knots(j)) exit
            k = k + 1
         end do
         m = m + 1
         knots(m) = knots(j)
         probabilities(m) = sum(table_probabilities(j:k)) / real(k - j + 1, dp)
         j = k + 1
      end do
      scaled_dip = sqrt(real(n, dp)) * statistic
      if (scaled_dip <= knots(1)) then
         ordinate = probabilities(1)
      else if (scaled_dip >= knots(m)) then
         ordinate = probabilities(m)
      else
         j = 1
         do while (knots(j + 1) < scaled_dip)
            j = j + 1
         end do
         weight = (scaled_dip - knots(j)) / (knots(j + 1) - knots(j))
         ordinate = probabilities(j) + weight * (probabilities(j + 1) - probabilities(j))
      end if
      p_value = max(0.0_dp, min(1.0_dp, 1.0_dp - ordinate))
   end subroutine dip_p_value

   pure subroutine dip_test(x, result, uniform_samples)
      real(dp), intent(in) :: x(:) !! Observations in arbitrary order; NaNs are omitted as in R complete.cases, infinities rejected.
      type(dip_test_result), intent(out) :: result !! Observed statistic, p-value, retained sample size, status and calibration metadata.
      real(dp), intent(in), optional :: uniform_samples(:, :) !! Null draws: n cleaned observations by B positive replicates, in [0,1].
      real(dp), allocatable :: clean(:)
      type(dip_result) :: fit
      real(dp) :: null_dip
      integer :: b, exceedances

      clean = pack(x, .not. ieee_is_nan(x))
      call dip_compute(clean, fit)
      result%n = fit%n
      result%statistic = fit%statistic
      result%status = fit%status
      result%p_value = ieee_value(0.0_dp, ieee_quiet_nan)
      if (fit%status /= dip_ok) return
      if (fit%n <= 3) then
         result%p_value = 1.0_dp
         return
      end if
      if (present(uniform_samples)) then
         result%simulated = .true.
         result%replicates = size(uniform_samples, 2)
         if (size(uniform_samples, 1) /= fit%n .or. result%replicates < 1) then
            result%status = dip_invalid_simulation
            return
         end if
         if (any(.not. ieee_is_finite(uniform_samples))) then
            result%status = dip_invalid_simulation
            return
         end if
         if (any(uniform_samples < 0.0_dp) .or. any(uniform_samples > 1.0_dp)) then
            result%status = dip_invalid_simulation
            return
         end if
         exceedances = 0
         do b = 1, result%replicates
            null_dip = dip(uniform_samples(:, b))
            if (fit%statistic <= null_dip) exceedances = exceedances + 1
         end do
         result%p_value = real(exceedances, dp) / real(result%replicates, dp)
      else
         call dip_p_value(fit%statistic, fit%n, result%p_value, result%status, result%asymptotic)
      end if
   end subroutine dip_test

   subroutine dip_test_simulate(x, result, replicates)
      real(dp), intent(in) :: x(:) !! Observations with NaNs removed; infinities and an empty retained sample are errors.
      type(dip_test_result), intent(out) :: result !! Simulation-based test using the compiler's intrinsic random-number generator.
      integer, intent(in), optional :: replicates !! Number of Monte Carlo replicates, positive, default 2000; ignored for n<=3.
      real(dp), allocatable :: sample(:)
      integer :: count_replicates, b, exceedances
      count_replicates = 2000
      if (present(replicates)) count_replicates = replicates
      call dip_test(x, result)
      if (result%status /= dip_ok .or. result%n <= 3) return
      if (count_replicates < 1) then
         result%status = dip_invalid_simulation
         result%p_value = ieee_value(0.0_dp, ieee_quiet_nan)
         return
      end if
      result%simulated = .true.
      result%asymptotic = .false.
      result%replicates = count_replicates
      allocate(sample(result%n))
      exceedances = 0
      do b = 1, count_replicates
         call random_number(sample)
         if (result%statistic <= dip(sample)) exceedances = exceedances + 1
      end do
      result%p_value = real(exceedances, dp) / real(count_replicates, dp)
   end subroutine dip_test_simulate

end module diptest_api
