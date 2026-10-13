! SPDX-License-Identifier: GPL-3.0-or-later
! Modern Fortran translation of infotheo 1.2.0.1 by Patrick E. Meyer.
module infotheo
   use, intrinsic :: ieee_arithmetic, only : ieee_is_nan, ieee_is_finite, ieee_quiet_nan, ieee_value
   use r_kinds, only : dp
   use r_special, only : r_digamma
   use r_ordering, only : r_order
   use r_sorting, only : r_sort
   implicit none
   private

   integer, parameter, public :: na_integer = -huge(0) - 1
   integer, parameter, public :: infotheo_ok = 0
   integer, parameter, public :: infotheo_invalid = 1
   integer, parameter, public :: max_interaction_variables = 20
   public :: dp, entropy, condentropy, mutinformation, mutinformation_matrix
   public :: condinformation, multiinformation, interinformation, natstobits, discretize

   interface entropy
      module procedure entropy_matrix, entropy_vector
   end interface entropy
contains

   pure function entropy_vector(x, method) result(h)
      integer, intent(in) :: x(:) !! Discrete observations; na_integer marks a missing observation.
      character(len=*), intent(in), optional :: method !! Estimator: emp (default), mm, sg, sg_r, or shrink.
      real(dp) :: h
      h = entropy_matrix(reshape(x, [size(x), 1]), method)
   end function entropy_vector

   pure function entropy_matrix(x, method) result(h)
      integer, intent(in) :: x(:, :) !! Observations by variables; rows with any na_integer are omitted.
      character(len=*), intent(in), optional :: method !! Estimator: emp (default), mm, sg, sg_r, or shrink.
      real(dp) :: h
      integer, allocatable :: rows(:), permutation(:), frequencies(:)
      integer :: n, k, i, j
      real(dp) :: samples, states, beta, lambda, concentration, squared_probability
      character(len=:), allocatable :: estimator

      h = ieee_value(0.0_dp, ieee_quiet_nan)
      estimator = 'emp'
      if (present(method)) estimator = trim(method)
      select case (estimator)
      case ('emp', 'mm', 'sg', 'sg_r', 'shrink')
      case default
         return
      end select
      n = count(all(x /= na_integer, dim = 2))
      if (n == 0) return
      if (size(x, 2) == 0) then
         h = 0.0_dp
         return
      end if
      allocate(rows(n), frequencies(n))
      j = 0
      do i = 1, size(x, 1)
         if (any(x(i, :) == na_integer)) cycle
         j = j + 1
         rows(j) = i
      end do
      ! Stable sorts, least significant variable first, give lexicographic row order.
      ! Every default integer is exactly representable by core's binary64 dp.
      do j = size(x, 2), 1, -1
         call r_order(real(x(rows, j), dp), permutation)
         rows = rows(permutation)
      end do
      frequencies = 0
      k = 1
      frequencies(k) = 1
      do i = 2, n
         if (any(x(rows(i), :) /= x(rows(i - 1), :))) k = k + 1
         frequencies(k) = frequencies(k) + 1
      end do
      samples = real(n, dp)
      states = real(k, dp)
      if (k == 1) then
         h = 0.0_dp
         return
      end if
      select case (estimator)
      case ('emp', 'mm')
         h = -sum((real(frequencies(:k), dp) / samples) * log(real(frequencies(:k), dp) / samples))
         if (estimator == 'mm') h = h + (states - 1.0_dp) / (2.0_dp * samples)
      case ('sg', 'sg_r')
         beta = 1.0_dp / states
         if (estimator == 'sg_r') beta = 0.0_dp
         h = dirichlet_entropy(real(frequencies(:k), dp) + beta)
      case ('shrink')
         ! Algebraically identical to upstream, without integer-square overflow.
         squared_probability = sum((real(frequencies(:k), dp) / samples) ** 2)
         concentration = states * squared_probability - 1.0_dp
         if (concentration <= 0.0_dp) then
            h = log(states)
            return
         end if
         lambda = states * (1.0_dp - squared_probability) / ((samples - 1.0_dp) * concentration)
         lambda = max(0.0_dp, min(1.0_dp, lambda))
         if (lambda >= 1.0_dp) then
            h = log(states)
         else
            beta = (lambda / (1.0_dp - lambda)) * samples / states
            h = dirichlet_entropy(real(frequencies(:k), dp) + beta)
         end if
      end select
   end function entropy_matrix

   pure function dirichlet_entropy(alpha) result(h)
      real(dp), intent(in) :: alpha(:) !! Positive posterior Dirichlet parameters of observed states.
      real(dp) :: h, total
      total = sum(alpha)
      h = sum((alpha / total) * (r_digamma(total + 1.0_dp) - r_digamma(alpha + 1.0_dp)))
   end function dirichlet_entropy

   pure function condentropy(x, y, method) result(h)
      integer, intent(in) :: x(:, :) !! Target observations by variables, using na_integer for missing cells.
      integer, intent(in), optional :: y(:, :) !! Conditioning block with the same row count; absent means H(X).
      character(len=*), intent(in), optional :: method !! Entropy estimator; default emp, with sg_r for upstream SG compatibility.
      real(dp) :: h
      if (.not. present(y)) then
         h = entropy(x, method)
         return
      end if
      h = ieee_value(0.0_dp, ieee_quiet_nan)
      if (size(x, 1) /= size(y, 1)) return
      h = entropy(join_columns(x, y), method) - entropy(y, method)
   end function condentropy

   pure function mutinformation(x, y, method) result(h)
      integer, intent(in) :: x(:, :) !! First discrete variable block, observations in rows; missing cells are na_integer.
      integer, intent(in) :: y(:, :) !! Second block with the same row count as X.
      character(len=*), intent(in), optional :: method !! Entropy estimator; default emp, with sg_r for upstream SG compatibility.
      real(dp) :: h
      h = ieee_value(0.0_dp, ieee_quiet_nan)
      if (size(x, 1) /= size(y, 1)) return
      h = entropy(x, method) + entropy(y, method) - entropy(join_columns(x, y), method)
      ! The scalar R operation truncates negative estimator values; its matrix operation does not.
      if (.not. ieee_is_nan(h)) h = max(0.0_dp, h)
   end function mutinformation

   pure function mutinformation_matrix(x, method) result(h)
      integer, intent(in) :: x(:, :) !! Discrete observations by variables, with na_integer for missing cells.
      character(len=*), intent(in), optional :: method !! Entropy estimator; default emp, with sg_r for upstream SG compatibility.
      real(dp), allocatable :: h(:, :)
      integer :: i, j, p
      p = size(x, 2)
      allocate(h(p, p))
      do i = 1, p
         h(i, i) = entropy(x(:, i), method)
      end do
      do i = 1, p
         do j = 1, i - 1
            h(i, j) = h(i, i) + h(j, j) - entropy(x(:, [i, j]), method)
            h(j, i) = h(i, j)
         end do
      end do
   end function mutinformation_matrix

   pure function condinformation(x, y, s, method) result(h)
      integer, intent(in) :: x(:, :) !! First variable block; rows are observations and na_integer marks missing cells.
      integer, intent(in) :: y(:, :) !! Second variable block with the same number of observations.
      integer, intent(in), optional :: s(:, :) !! Complete conditioning block; absent means scalar mutual information.
      character(len=*), intent(in), optional :: method !! Entropy estimator; default emp, with sg_r for upstream SG compatibility.
      real(dp) :: h
      if (.not. present(s)) then
         h = mutinformation(x, y, method)
         return
      end if
      h = ieee_value(0.0_dp, ieee_quiet_nan)
      if (size(x, 1) /= size(y, 1)) return
      if (size(x, 1) /= size(s, 1)) return
      h = entropy(join_columns(s, x), method) + entropy(join_columns(s, y), method) &
         - entropy(s, method) - entropy(join_columns(join_columns(s, x), y), method)
   end function condinformation

   pure function multiinformation(x, method) result(h)
      integer, intent(in) :: x(:, :) !! Discrete observations by variables; missing cells are na_integer.
      character(len=*), intent(in), optional :: method !! Entropy estimator; default emp, with sg_r for upstream SG compatibility.
      real(dp) :: h
      integer :: j
      h = -entropy(x, method)
      do j = 1, size(x, 2)
         h = h + entropy(x(:, j), method)
      end do
   end function multiinformation

   pure function interinformation(x, method) result(h)
      integer, intent(in) :: x(:, :) !! Observations by at most max_interaction_variables discrete variables.
      character(len=*), intent(in), optional :: method !! Entropy estimator; default emp, with sg_r for upstream SG compatibility.
      real(dp) :: h
      integer :: mask, p, j, k, columns(size(x, 2))
      p = size(x, 2)
      h = ieee_value(0.0_dp, ieee_quiet_nan)
      if (p > max_interaction_variables) return
      if (p == 0) then
         h = entropy(x, method)
         return
      end if
      h = 0.0_dp
      do mask = 1, 2 ** p - 1
         k = 0
         do j = 1, p
            if (.not. btest(mask, j - 1)) cycle
            k = k + 1
            columns(k) = j
         end do
         if (mod(k, 2) == 1) then
            h = h + entropy(x(:, columns(:k)), method)
         else
            h = h - entropy(x(:, columns(:k)), method)
         end if
      end do
   end function interinformation

   pure function join_columns(x, y) result(xy)
      integer, intent(in) :: x(:, :) !! Left matrix with an already-validated matching row count.
      integer, intent(in) :: y(:, :) !! Right matrix whose columns follow the columns of X.
      integer :: xy(size(x, 1), size(x, 2) + size(y, 2))
      xy(:, :size(x, 2)) = x
      xy(:, size(x, 2) + 1:) = y
   end function join_columns

   pure elemental function natstobits(h) result(bits)
      real(dp), intent(in) :: h !! Information in nats; IEEE NaNs and infinities propagate normally.
      real(dp) :: bits
      bits = h / log(2.0_dp)
   end function natstobits

   pure subroutine discretize(x, bins, disc, nbins, status)
      real(dp), intent(in) :: x(:, :) !! Continuous observations by variables; NaNs represent missing cells; infinities invalid.
      integer, allocatable, intent(out) :: bins(:, :) !! Same shape as X, labels 1:nbins and na_integer for missing/invalid cells.
      character(len=*), intent(in), optional :: disc !! equalfreq (default), equalwidth, or globalequalwidth.
      integer, intent(in), optional :: nbins !! Positive bin count, default floor(nrows ** (1/3)); global method requires at least two.
      integer, intent(out), optional :: status !! infotheo_ok on success; infotheo_invalid for invalid method, count, or infinity.
      character(len=:), allocatable :: strategy
      real(dp), allocatable :: sorted(:), cuts(:)
      real(dp) :: lo, hi, value, fraction
      integer :: b, i, j, n, frequency, remainder, position, k

      allocate(bins(size(x, 1), size(x, 2)))
      bins = na_integer
      if (present(status)) status = infotheo_invalid
      strategy = 'equalfreq'
      if (present(disc)) strategy = trim(disc)
      select case (strategy)
      case ('equalfreq', 'equalwidth', 'globalequalwidth')
      case default
         return
      end select
      b = int(real(size(x, 1), dp) ** (1.0_dp / 3.0_dp))
      if (present(nbins)) b = nbins
      if (b < 1) return
      if (strategy == 'globalequalwidth' .and. b < 2) return
      do j = 1, size(x, 2)
         do i = 1, size(x, 1)
            if (ieee_is_nan(x(i, j))) cycle
            if (.not. ieee_is_finite(x(i, j))) return
         end do
      end do
      if (strategy == 'globalequalwidth') then
         if (all(ieee_is_nan(x))) then
            if (present(status)) status = infotheo_ok
            return
         end if
         lo = minval(x, mask = .not. ieee_is_nan(x))
         hi = maxval(x, mask = .not. ieee_is_nan(x))
         ! R cut() uses right-closed intervals and expands the outside bounds.
         ! Outside-bound expansion cannot change finite, in-range observations.
         if (hi == lo) then
            where (.not. ieee_is_nan(x)) bins = b / 2 + mod(b, 2)
         else
            do j = 1, size(x, 2)
               do i = 1, size(x, 1)
                  if (ieee_is_nan(x(i, j))) cycle
                  fraction = range_fraction(x(i, j), lo, hi)
                  bins(i, j) = max(1, min(b, ceiling(fraction * real(b, dp))))
               end do
            end do
         end if
      else
         ! Validate all column counts before changing the all-missing failure result.
         if (strategy == 'equalfreq') then
            do j = 1, size(x, 2)
               n = count(.not. ieee_is_nan(x(:, j)))
               if (n > 0 .and. b > n) return
            end do
         end if
         do j = 1, size(x, 2)
            n = count(.not. ieee_is_nan(x(:, j)))
            if (n == 0) cycle
            if (strategy == 'equalfreq') then
               call r_sort(x(:, j), sorted, na_rm = .true.)
               allocate(cuts(max(0, b - 1)))
               frequency = n / b
               remainder = mod(n, b)
               position = frequency
               do k = 1, b - 1
                  if (remainder > 0) then
                     cuts(k) = sorted(position + 1)
                     remainder = remainder - 1
                  else
                     cuts(k) = sorted(position)
                  end if
                  position = position + frequency
               end do
               do i = 1, size(x, 1)
                  if (ieee_is_nan(x(i, j))) cycle
                  bins(i, j) = first_upper_cut(x(i, j), cuts)
               end do
               deallocate(cuts)
            else
               lo = minval(x(:, j), mask = .not. ieee_is_nan(x(:, j)))
               hi = maxval(x(:, j), mask = .not. ieee_is_nan(x(:, j)))
               do i = 1, size(x, 1)
                  value = x(i, j)
                  if (ieee_is_nan(value)) cycle
                  if (hi == lo) then
                     bins(i, j) = 1
                  else
                     fraction = range_fraction(value, lo, hi)
                     if (fraction >= 1.0_dp) then
                        bins(i, j) = b
                     else
                        bins(i, j) = min(b, 1 + int(fraction * real(b, dp)))
                     end if
                  end if
               end do
            end if
         end do
      end if
      if (present(status)) status = infotheo_ok
   end subroutine discretize

   pure function range_fraction(value, lo, hi) result(fraction)
      real(dp), intent(in) :: value !! Finite observation between lo and hi, inclusive.
      real(dp), intent(in) :: lo !! Finite lower range endpoint, strictly smaller than hi.
      real(dp), intent(in) :: hi !! Finite upper range endpoint, strictly larger than lo.
      real(dp) :: fraction, scale
      if (lo < 0.0_dp .and. hi > 0.0_dp) then
         scale = max(abs(lo), abs(hi))
         fraction = (value / scale - lo / scale) / (hi / scale - lo / scale)
      else
         fraction = (value - lo) / (hi - lo)
      end if
      fraction = max(0.0_dp, min(1.0_dp, fraction))
   end function range_fraction

   pure function first_upper_cut(value, cuts) result(bin)
      real(dp), intent(in) :: value !! Nonmissing finite observation to assign to a right-closed bin.
      real(dp), intent(in) :: cuts(:) !! Nondecreasing internal cutpoints; equal values stay in the earliest matching bin.
      integer :: bin, left, right, middle
      left = 1
      right = size(cuts)
      do while (left <= right)
         middle = left + (right - left) / 2
         if (value <= cuts(middle)) then
            right = middle - 1
         else
            left = middle + 1
         end if
      end do
      bin = left
   end function first_upper_cut
end module infotheo
